use clap::Parser;
use dirs::home_dir;
use genesis::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::Deserialize;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// A simple program to initialize variables from a TOML file
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the .toml configuration file (positional argument)
    config: String,
}

/// Struct to represent the "output" section in the TOML file
#[derive(Debug, Deserialize)]
struct InputConfig {
    mol_path: String,
    par_path: String,
    pos_path: Option<String>,
    vel_path: Option<String>,
    rst_path: Option<String>,
}

/// Struct to represent the "output" section in the TOML file
#[derive(Debug, Deserialize)]
struct OutputConfig {
    csv_path: String,
    dcd_path: String,
    log_path: String,
    rst_path: String,
    xyz_path: String,
    csv_freq: u32,
    dcd_freq: u32,
    rst_freq: u32,
}

/// Struct to represent the "dynamics" section in the TOML file
#[derive(Debug, Deserialize)]
struct DynamicsConfig {
    time_step: f32,
    num_steps: u32,
    temperature: f64,
    remove_com_v_freq: u32,
}

/// Struct to represent the "boundary" section in the TOML file
#[derive(Debug, Deserialize)]
struct BoundaryConfig {
    length_x: f32,
    length_y: f32,
    length_z: f32,
}

/// Struct to represent the "rng" section in the TOML file
#[derive(Debug, Deserialize)]
struct RngConfig {
    seed: u64,
}

/// Main config struct combining all sections
#[derive(Debug, Deserialize)]
struct Config {
    input: InputConfig,
    output: OutputConfig,
    dynamics: DynamicsConfig,
    boundary: BoundaryConfig,
    rng: RngConfig,
}

/// Expand `~` to the home directory
fn expand_tilde<P: AsRef<Path>>(path: P) -> PathBuf {
    let path_str = path.as_ref().to_string_lossy().to_string();
    if path_str.starts_with("~") {
        if let Some(home) = home_dir() {
            return PathBuf::from(path_str.replacen("~", &home.to_string_lossy(), 1));
        }
    }
    PathBuf::from(path_str)
}

/// Resolve the path to handle relative vs. absolute paths
fn resolve_path<P: AsRef<Path>>(path: P, base_dir: &Path) -> PathBuf {
    let path = expand_tilde(Path::new(path.as_ref()));
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base_dir.join(path)
    }
}

fn main() {
    // Initialize section timers
    let mut setup_timer = Timer::default();
    let mut output_timer = Timer::default();
    let mut dynamics_timer = Timer::default();

    ////////////////////////////////////////
    setup_timer.start();
    //**************************************

    // Parse command-line arguments
    let args = Args::parse();

    // Read the TOML file from the given path
    let config_path = Path::new(&args.config);
    let config_content = match fs::read_to_string(config_path) {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Failed to read the config file {:?}: {}", config_path, e);
            std::process::exit(1);
        }
    };

    // Parse the TOML file into the `Config` struct
    let config: Config = match toml::from_str(&config_content) {
        Ok(config) => config,
        Err(e) => {
            eprintln!("Failed to parse the config file: {}", e);
            std::process::exit(1);
        }
    };

    // Get the directory where the config file is located
    let config_dir = config_path.parent().unwrap_or_else(|| Path::new(""));

    // Resolve input file paths
    let inp_mol_path = resolve_path(&config.input.mol_path, config_dir);
    let inp_par_path = resolve_path(&config.input.par_path, config_dir);
    let inp_pos_path = &config.input.pos_path.map(|p| resolve_path(p, config_dir));
    let inp_vel_path = &config.input.vel_path.map(|p| resolve_path(p, config_dir));
    let inp_rst_path = &config.input.rst_path.map(|p| resolve_path(p, config_dir));

    // Resolve output file paths
    let out_csv_path = resolve_path(&config.output.csv_path, config_dir);
    let out_dcd_path = resolve_path(&config.output.dcd_path, config_dir);
    let out_log_path = resolve_path(&config.output.log_path, config_dir);
    let out_rst_path = resolve_path(&config.output.rst_path, config_dir);
    let out_xyz_path = resolve_path(&config.output.xyz_path, config_dir);

    // Initialize reporters
    let mut csv_reporter = match CSVReporter::with_path(out_csv_path) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create CSV file: {}", e);
            std::process::exit(1);
        }
    };

    let mut dcd_reporter = match DCDReporter::with_path(out_dcd_path) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create DCD file: {}", e);
            std::process::exit(1);
        }
    };

    let mut log_reporter = match LOGReporter::with_path(out_log_path) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create LOG file: {}", e);
            std::process::exit(1);
        }
    };

    let mut rst_reporter = match RSTReporter::with_path(out_rst_path) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create RST file: {}", e);
            std::process::exit(1);
        }
    };

    let mut xyz_reporter = match XYZReporter::with_path(out_xyz_path) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create XYZ file: {}", e);
            std::process::exit(1);
        }
    };

    // Initialize observers
    let mut et_obs = TotalEnergyObserver::new();
    let mut ue_obs = PotentialEnergyObserver::new();
    let mut ke_obs = KineticEnergyObserver::new();
    let mut te_obs = TemperatureObserver::new();
    let mut vi_obs = VirialObserver::new();
    let mut vo_obs = VolumeObserver::new();
    let mut pr_obs = PressureObserver::new();
    let mut df_obs = DegreesOfFreedomObserver::new();

    let mut rng: StdRng = StdRng::seed_from_u64(config.rng.seed);
    let e: f32 = 1.003;
    let s: f32 = 0.340;

    // Initialize the system
    let mut system = match inp_rst_path {
        Some(rst_path) => {
            // Deserialize the system from the RST file
            let file = match File::open(rst_path) {
                Ok(file) => file,
                Err(e) => {
                    eprintln!("Failed to open RST file: {}", e);
                    std::process::exit(1);
                }
            };
            let reader = BufReader::new(file);
            match serde_json::from_reader(reader) {
                Ok(system) => system,
                Err(e) => {
                    eprintln!("Failed to deserialize RST file: {}", e);
                    std::process::exit(1);
                }
            }
        }
        None => {
            // Initialize system builder
            let builder = System::builder();

            // Set the boundary
            let builder = builder.with_cuboid_boundary(
                config.boundary.length_x,
                config.boundary.length_y,
                config.boundary.length_z,
            );

            // Set the masses, charges
            let mol_parser = match MolParser::with_path(inp_mol_path) {
                Ok(parser) => parser,
                Err(e) => {
                    eprintln!("Failed to open MOL file: {}", e);
                    std::process::exit(1);
                }
            };
            let parser_result = match mol_parser.parse() {
                Ok(parser_result) => parser_result,
                Err(e) => {
                    eprintln!("Failed to parse POS file: {}", e);
                    std::process::exit(1);
                }
            };
            let builder = builder.with_masses(parser_result.masses);
            let builder = builder.with_charges(parser_result.charges);
            let builder = builder.with_classes(parser_result.classes);
            let builder = builder.with_names(parser_result.names);

            // Set the positions
            let builder = match inp_pos_path {
                Some(pos_path) => {
                    let pos_parser = match PosParser::with_path(pos_path) {
                        Ok(parser) => parser,
                        Err(e) => {
                            eprintln!("Failed to open POS file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    let parser_result = match pos_parser.parse() {
                        Ok(parser_result) => parser_result,
                        Err(e) => {
                            eprintln!("Failed to parse POS file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    builder.with_positions(parser_result.positions)
                }
                None => builder,
            };

            // Set the velocities
            let builder = match inp_vel_path {
                Some(vel_path) => {
                    let vel_parser = match VelParser::with_path(vel_path) {
                        Ok(parser) => parser,
                        Err(e) => {
                            eprintln!("Failed to open VEL file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    let velocities = match vel_parser.parse() {
                        Ok(velocities) => velocities,
                        Err(e) => {
                            eprintln!("Failed to parse VEL file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    builder.with_velocities(velocities)
                }
                None => builder.with_random_velocities(config.dynamics.temperature, &mut rng),
            };

            // Build the system
            builder.build()
        }
    };

    // Remove the center of mass velocity
    system.remove_v_com();

    // Observe degrees of freedom
    df_obs.observe(&system);

    // Check output frequencies
    if config.dynamics.num_steps % config.output.csv_freq != 0 {
        eprintln!(
            "Warning: num_steps is not a multiple of csv_freq: {} % {} != 0",
            config.dynamics.num_steps, config.output.csv_freq
        );
    }

    if config.dynamics.num_steps % config.output.dcd_freq != 0 {
        eprintln!(
            "Warning: num_steps is not a multiple of dcd_freq: {} % {} != 0",
            config.dynamics.num_steps, config.output.dcd_freq
        );
    }

    if config.dynamics.num_steps % config.output.rst_freq != 0 {
        eprintln!(
            "Warning: num_steps is not a multiple of rst_freq: {} % {} != 0",
            config.dynamics.num_steps, config.output.rst_freq
        );
    }

    // Setup integrator
    let dt: f32 = config.dynamics.time_step;
    let dt_half: f32 = dt * 0.5;

    //**************************************
    setup_timer.stop();
    ////////////////////////////////////////

    ////////////////////////////////////////
    output_timer.start();
    //**************************************

    dcd_reporter
        .write_header(
            system.n as u32,
            config.dynamics.num_steps / config.output.dcd_freq + 1,
        )
        .expect("Failed to write DCD header");
    dcd_reporter
        .write_report(&system)
        .expect("Failed to write DCD frame");
    csv_reporter
        .write_header()
        .expect("Failed to write CSV header");

    ke_obs.observe(&system);
    ue_obs.observe(&system, e, s);
    et_obs.observe(&ke_obs, &ue_obs);
    te_obs.observe(&ke_obs, &df_obs);
    vi_obs.observe(&system, e, s);
    vo_obs.observe(&system);
    pr_obs.observe(&vo_obs, &vi_obs, &te_obs, &df_obs);
    csv_reporter
        .write_report(
            0, &et_obs, &ue_obs, &ke_obs, &te_obs, &vi_obs, &vo_obs, &pr_obs,
        )
        .expect("Failed to write CSV frame");

    //**************************************
    output_timer.stop();
    ////////////////////////////////////////

    for step in 1..=config.dynamics.num_steps {
        ////////////////////////////////////////
        dynamics_timer.start();
        //**************************************

        for i in 0..system.n {
            system.r[i] += system.v[i] * dt_half;
        }

        for f in system.f.iter_mut() {
            f.fill(0.0);
        }

        compute_force(&mut system, e, s);

        for i in 0..system.n {
            system.v[i] += system.f[i] * dt / system.m[i];
            system.r[i] += system.v[i] * dt_half;
        }

        if step % config.dynamics.remove_com_v_freq == 0 {
            system.remove_v_com();
        }

        //**************************************
        dynamics_timer.stop();
        ////////////////////////////////////////

        ////////////////////////////////////////
        output_timer.start();
        //**************************************

        if step % config.output.dcd_freq == 0 {
            match dcd_reporter.write_report(&system) {
                Ok(_) => (),
                Err(e) => {
                    eprintln!("Failed to write DCD frame: {}", e);
                    std::process::exit(1);
                }
            }
        }

        if step % config.output.csv_freq == 0 {
            ke_obs.observe(&system);
            ue_obs.observe(&system, e, s);
            et_obs.observe(&ke_obs, &ue_obs);
            te_obs.observe(&ke_obs, &df_obs);
            vi_obs.observe(&system, e, s);
            vo_obs.observe(&system);
            pr_obs.observe(&vo_obs, &vi_obs, &te_obs, &df_obs);
            match csv_reporter.write_report(
                step, &et_obs, &ue_obs, &ke_obs, &te_obs, &vi_obs, &vo_obs, &pr_obs,
            ) {
                Ok(_) => (),
                Err(e) => {
                    eprintln!("Failed to write CSV frame: {}", e);
                    std::process::exit(1);
                }
            }
        }

        if step % config.output.rst_freq == 0 {
            match rst_reporter.write_report(&system) {
                Ok(_) => (),
                Err(e) => {
                    eprintln!("Failed to write RST frame: {}", e);
                    std::process::exit(1);
                }
            }
        }

        //**************************************
        output_timer.stop();
        ////////////////////////////////////////
    }

    ////////////////////////////////////////
    output_timer.start();
    //**************************************

    match xyz_reporter.write_report(&system) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Failed to write XYZ frame: {}", e);
            std::process::exit(1);
        }
    }

    //**************************************
    output_timer.stop();
    ////////////////////////////////////////

    // Report the results
    log_reporter
        .write_report(
            setup_timer.elapsed(),
            output_timer.elapsed(),
            dynamics_timer.elapsed(),
        )
        .expect("Failed to write LOG report");
}

fn compute_force(system: &mut System, e: f32, s: f32) {
    let s2: f32 = s * s;

    for i in 0..system.n {
        let ri = system.r[i];
        for j in (i + 1)..system.n {
            let mut dr = system.r[j] - ri;
            let offset = dr.component_div(&system.b).map(|x| x.round());
            dr -= system.b.component_mul(&offset);

            let r2: f32 = 1.0 / dr.norm_squared();
            let c2: f32 = s2 * r2;
            let c4: f32 = c2 * c2;
            let c6: f32 = c4 * c2;

            let force: f32 = 48.0 * e * c6 * (c6 - 0.5) * r2;

            system.f[i] -= force * dr;
            system.f[j] += force * dr;
        }
    }
}
