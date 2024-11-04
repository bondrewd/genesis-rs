use clap::Parser;
use dirs::home_dir;
use genesis::observer::{
    DegreesOfFreedomObserver, KineticEnergyObserver, PotentialEnergyObserver, PressureObserver,
    TemperatureObserver, TotalEnergyObserver, VirialObserver, VolumeObserver,
};
use genesis::reporter::csv::CSVReporter;
use genesis::reporter::dcd::DCDReporter;
use genesis::reporter::log::LOGReporter;
use genesis::reporter::xyz::XYZReporter;
use genesis::system::System;
use genesis::timer::Timer;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

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

/// A simple program to initialize variables from a TOML file
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the .toml configuration file (positional argument)
    config: String,
}

/// Struct to represent the "output" section in the TOML file
#[derive(Debug, Deserialize)]
struct OutputConfig {
    out_csv_path: String,
    out_dcd_path: String,
    out_log_path: String,
    out_xyz_path: String,
    out_csv_freq: u32,
    out_dcd_freq: u32,
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

    // Resolve the output file paths
    let out_csv_path = resolve_path(&config.output.out_csv_path, config_dir);
    let out_dcd_path = resolve_path(&config.output.out_dcd_path, config_dir);
    let out_log_path = resolve_path(&config.output.out_log_path, config_dir);
    let out_xyz_path = resolve_path(&config.output.out_xyz_path, config_dir);

    // Initialize writers
    let mut csv_reporter = CSVReporter::with_path(out_csv_path).expect("Failed to create CSV file");
    let mut dcd_reporter = DCDReporter::with_path(out_dcd_path).expect("Failed to create DCD file");
    let mut log_reporter = LOGReporter::with_path(out_log_path).expect("Failed to create LOG file");
    let mut xyz_reporter = XYZReporter::with_path(out_xyz_path).expect("Failed to create XYZ file");

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

    let n: usize = 100;
    let e: f32 = 1.003;
    let s: f32 = 0.340;
    let mut system = System::builder()
        .n(n)
        .with_cuboid_boundary(
            config.boundary.length_x,
            config.boundary.length_y,
            config.boundary.length_z,
        )
        .m(vec![39.948; n])
        .with_random_charges(&mut rng)
        .with_random_positions(&mut rng)
        .with_random_velocities(config.dynamics.temperature, &mut rng)
        .build();

    // Remove the center of mass velocity
    system.remove_v_com();

    // Observe degrees of freedom
    df_obs.observe(&system);

    assert!(
        config.output.out_csv_freq % config.dynamics.remove_com_v_freq == 0,
        "out_ene_freq is not a multiple of rem_com_freq"
    );
    assert!(
        config.output.out_dcd_freq % config.dynamics.remove_com_v_freq == 0,
        "out_dcd_freq is not a multiple of rem_com_freq"
    );

    let dt: f32 = config.dynamics.time_step;
    let dt_half: f32 = dt * 0.5;
    let n_steps: u32 = config.dynamics.num_steps;

    //**************************************
    setup_timer.stop();
    ////////////////////////////////////////

    ////////////////////////////////////////
    output_timer.start();
    //**************************************

    dcd_reporter
        .write_header(n as u32, n_steps / config.output.out_dcd_freq + 1)
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

    for step in 1..=n_steps {
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

        if step % config.output.out_dcd_freq == 0 {
            dcd_reporter
                .write_report(&system)
                .expect("Failed to write DCD frame");
        }

        if step % config.output.out_csv_freq == 0 {
            ke_obs.observe(&system);
            ue_obs.observe(&system, e, s);
            et_obs.observe(&ke_obs, &ue_obs);
            te_obs.observe(&ke_obs, &df_obs);
            vi_obs.observe(&system, e, s);
            vo_obs.observe(&system);
            pr_obs.observe(&vo_obs, &vi_obs, &te_obs, &df_obs);
            csv_reporter
                .write_report(
                    step, &et_obs, &ue_obs, &ke_obs, &te_obs, &vi_obs, &vo_obs, &pr_obs,
                )
                .expect("Failed to write CSV frame");
        }

        //**************************************
        output_timer.stop();
        ////////////////////////////////////////
    }

    ////////////////////////////////////////
    output_timer.start();
    //**************************************

    xyz_reporter
        .write_report(&system)
        .expect("Failed to write XYZ frame");

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
