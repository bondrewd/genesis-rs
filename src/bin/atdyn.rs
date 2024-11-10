use clap::Parser;
use genesis::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// A simple program to initialize variables from a TOML file
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the .toml configuration file (positional argument)
    config: String,
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
    let config = Config::try_from(Path::new(&args.config)).unwrap();

    // Initialize reporters
    let mut csv_reporter = match CSVReporter::with_path(config.output.csv_path.unwrap()) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create CSV file: {}", e);
            std::process::exit(1);
        }
    };

    let mut dcd_reporter = match DCDReporter::with_path(config.output.dcd_path.unwrap()) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create DCD file: {}", e);
            std::process::exit(1);
        }
    };

    let mut log_reporter = match LOGReporter::with_path(config.output.log_path.unwrap()) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create LOG file: {}", e);
            std::process::exit(1);
        }
    };

    let mut rst_reporter = match RSTReporter::with_path(config.output.rst_path.unwrap()) {
        Ok(reporter) => reporter,
        Err(e) => {
            eprintln!("Failed to create RST file: {}", e);
            std::process::exit(1);
        }
    };

    let mut xyz_reporter = match XYZReporter::with_path(config.output.xyz_path.unwrap()) {
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

    // Initialize the system
    let mut system = match config.input.rst_path {
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
            let builder = match config.input.mol_path {
                Some(mol_path) => {
                    let mol_parser = match MolParser::with_path(mol_path) {
                        Ok(parser) => parser,
                        Err(e) => {
                            eprintln!("Failed to open MOL file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    let parser_result = match mol_parser.parse() {
                        Ok(parser_result) => parser_result,
                        Err(e) => {
                            eprintln!("Failed to parse MOL file: {}", e);
                            std::process::exit(1);
                        }
                    };
                    builder
                        .with_masses(parser_result.masses)
                        .with_charges(parser_result.charges)
                        .with_classes(parser_result.classes)
                        .with_names(parser_result.names)
                }
                None => {
                    eprintln!("Error: MOL file is required when no RST file is provided");
                    std::process::exit(1);
                }
            };

            // Set the positions
            let builder = match config.input.pos_path {
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
                None => {
                    eprintln!("Error: POS file is required when no RST file is provided");
                    std::process::exit(1);
                }
            };

            // Set the velocities
            let builder = match config.input.vel_path {
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
            let mut system = builder.build();

            // Remove the center of mass velocity
            system.remove_v_com();

            system
        }
    };

    // Observe degrees of freedom
    df_obs.observe(&system);

    // Initialize force field
    let par_parser = match ParParser::with_path(config.input.par_path.unwrap()) {
        Ok(parser) => parser,
        Err(e) => {
            eprintln!("Failed to open PAR file: {}", e);
            std::process::exit(1);
        }
    };

    let par = match par_parser.parse() {
        Ok(ff) => ff,
        Err(e) => {
            eprintln!("Failed to parse PAR file: {}", e);
            std::process::exit(1);
        }
    };

    let ff = ForceField::with_parameters(par);

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
    ue_obs.observe(&system, &ff);
    et_obs.observe(&ke_obs, &ue_obs);
    te_obs.observe(&ke_obs, &df_obs);
    vi_obs.observe(&system, &ff);
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

        ff.update_force(&mut system);

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
            ue_obs.observe(&system, &ff);
            et_obs.observe(&ke_obs, &ue_obs);
            te_obs.observe(&ke_obs, &df_obs);
            vi_obs.observe(&system, &ff);
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
