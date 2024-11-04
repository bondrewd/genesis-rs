use clap::Parser;
use dirs::home_dir;
use genesis_rs::reporter::dcd::{DCDHeader, DCDReporter};
use genesis_rs::reporter::xyz::XYZReporter;
use genesis_rs::system::System;
use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

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

fn compute_potential_energy(system: &System, e: f32, s: f32) -> f64 {
    let mut potential_energy: f64 = 0.0;
    let s2: f32 = s * s;

    for i in 0..system.n {
        let ri = system.r[i];
        for rj in system.r.iter().skip(i + 1) {
            let mut dr = rj - ri;
            let offset = dr.component_div(&system.b).map(|x| x.round());
            dr -= system.b.component_mul(&offset);

            let r2: f32 = 1.0 / dr.norm_squared();
            let c2: f32 = s2 * r2;
            let c4: f32 = c2 * c2;
            let c6: f32 = c4 * c2;

            let energy: f32 = 4.0 * e * c6 * (c6 - 1.0);
            potential_energy += energy as f64;
        }
    }

    potential_energy // kJ/mol
}

fn compute_virial(system: &System, e: f32, s: f32) -> f64 {
    let mut total_virial: f64 = 0.0;
    let s2: f32 = s * s;

    for i in 0..system.n {
        let ri = system.r[i];
        for rj in system.r.iter().skip(i + 1) {
            let mut dr = rj - ri;
            let offset = dr.component_div(&system.b).map(|x| x.round());
            dr -= system.b.component_mul(&offset);

            let r2: f32 = 1.0 / dr.norm_squared();
            let c2: f32 = s2 * r2;
            let c4: f32 = c2 * c2;
            let c6: f32 = c4 * c2;

            let force: f32 = 48.0 * e * c6 * (c6 - 0.5) * r2;

            let virial: f32 = force / r2;
            total_virial += virial as f64;
        }
    }

    total_virial
}

fn compute_kinetic_energy(system: &System) -> f64 {
    let mut kinetic_energy: f64 = 0.0;

    for i in 0..system.n {
        let energy: f32 = 0.5_f32 * system.m[i] * system.v[i].dot(&system.v[i]);
        kinetic_energy += energy as f64;
    }

    kinetic_energy // kJ/mol
}

fn compute_temperature(kinetic_energy: f64, dof: u32) -> f64 {
    let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
    let temperature: f64 = 2.0 * kinetic_energy / (dof as f64 * boltzmann);
    temperature
}

fn compute_volume(system: &System) -> f64 {
    (system.b.x * system.b.y * system.b.z) as f64
}

fn compute_pressure(virial: f64, temperature: f64, volume: f64, dof: u32) -> f64 {
    let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
    let pressure: f64 = (virial + 3.0 * dof as f64 * boltzmann * temperature) / (3.0 * volume);
    pressure
}

fn remove_v_com(system: &mut System) {
    let mut v_com: Vector3<f32> = Vector3::zeros();
    let mut m_tot: f32 = 0.0;

    for i in 0..system.n {
        v_com += system.m[i] * system.v[i];
        m_tot += system.m[i];
    }

    v_com /= m_tot;

    for v in system.v.iter_mut() {
        *v -= v_com;
    }
}

fn display_output_header() {
    println!(
        "{:>10} {:>18} {:>18} {:>18} {:>18} {:>18} {:>18} {:>18}",
        "Step", "Total", "Potential", "Kinetic", "Temperature", "Virial", "Volume", "Pressure",
    );
}

#[allow(clippy::too_many_arguments)]
fn display_output(
    step: u32,
    total_energy: f64,
    potential_energy: f64,
    kinetic_energy: f64,
    temperature: f64,
    virial: f64,
    volume: f64,
    pressure: f64,
) {
    println!(
        "{:>10} {:>18.6} {:>18.6} {:>18.6} {:>18.6} {:>18.6} {:>18.6} {:>18.6}",
        step, total_energy, potential_energy, kinetic_energy, temperature, virial, volume, pressure,
    );
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
    out_xyz_path: String,
    out_dcd_path: String,
}

/// Struct to represent the "dynamics" section in the TOML file
#[derive(Debug, Deserialize)]
struct DynamicsConfig {
    time_step: f32,
    num_steps: u32,
    temperature: f64,
}

/// Struct to represent the "boundary" section in the TOML file
#[derive(Debug, Deserialize)]
struct BoundaryConfig {
    length_x: f32,
    length_y: f32,
    length_z: f32,
}

/// Main config struct combining all sections
#[derive(Debug, Deserialize)]
struct Config {
    output: OutputConfig,
    dynamics: DynamicsConfig,
    boundary: BoundaryConfig,
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

fn report_profile(
    setup_time: Duration,
    output_time: Duration,
    dynamics_time: Duration,
    total_time: Duration,
) {
    // Convert durations to seconds for easier formatting
    let setup_secs = setup_time.as_secs_f64();
    let output_secs = output_time.as_secs_f64();
    let dynamics_secs = dynamics_time.as_secs_f64();
    let total_secs = total_time.as_secs_f64();

    // Calculate percentages
    let setup_percent = (setup_secs / total_secs) * 100.0;
    let output_percent = (output_secs / total_secs) * 100.0;
    let dynamics_percent = (dynamics_secs / total_secs) * 100.0;

    // Print the profile report
    println!("Profile Report:");
    println!("setup    = {:>5.1}% {:.3}s ", setup_percent, setup_secs);
    println!("output   = {:>5.1}% {:.3}s ", output_percent, output_secs);
    println!(
        "dynamics = {:>5.1}% {:.3}s",
        dynamics_percent, dynamics_secs
    );
    println!("total    = 100.0% {:.3}s", total_secs);
}

fn main() {
    // Variables to store the duration of each section
    let mut setup_time = Duration::new(0, 0);
    let mut output_time = Duration::new(0, 0);
    let mut dynamics_time = Duration::new(0, 0);

    ////////////////////////////////////////
    let start = Instant::now();
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
    let out_xyz_path = resolve_path(&config.output.out_xyz_path, config_dir);
    let out_dcd_path = resolve_path(&config.output.out_dcd_path, config_dir);

    // Initialize writers
    let mut xyz_reporter = XYZReporter::with_path(out_xyz_path).expect("Failed to create XYZ file");
    let mut dcd_reporter = DCDReporter::with_path(out_dcd_path).expect("Failed to create DCD file");

    let mut rng: StdRng = StdRng::seed_from_u64(0);

    let n: usize = 100;
    let e: f32 = 1.003;
    let s: f32 = 0.340;
    let mut system = System::builder()
        .n(n)
        .b(Vector3::new(
            config.boundary.length_x,
            config.boundary.length_y,
            config.boundary.length_z,
        ))
        .m(vec![39.948; n])
        .with_random_charges(&mut rng)
        .with_random_positions(&mut rng)
        .with_random_velocities(config.dynamics.temperature, &mut rng)
        .build();
    let dof: u32 = 3 * system.n as u32 - 3;

    remove_v_com(&mut system);

    let out_ene_freq: u32 = 1000;
    let out_dcd_freq: u32 = 10;
    let rem_com_freq: u32 = 10;
    assert!(
        out_ene_freq % rem_com_freq == 0,
        "out_ene_freq is not a multiple of rem_com_freq"
    );
    assert!(
        out_dcd_freq % rem_com_freq == 0,
        "out_dcd_freq is not a multiple of rem_com_freq"
    );

    let dt: f32 = config.dynamics.time_step;
    let dt_half: f32 = dt * 0.5;
    let n_steps: u32 = config.dynamics.num_steps;

    //**************************************
    setup_time += start.elapsed();
    ////////////////////////////////////////

    ////////////////////////////////////////
    let start = Instant::now();
    //**************************************

    display_output_header();
    let mut potential_energy: f64 = compute_potential_energy(&system, e, s);
    let mut kinetic_energy: f64 = compute_kinetic_energy(&system);
    let mut total_energy: f64 = potential_energy + kinetic_energy;
    let mut temperature: f64 = compute_temperature(kinetic_energy, dof);
    let mut virial: f64 = compute_virial(&system, e, s);
    let mut volume: f64 = compute_volume(&system);
    let mut pressure: f64 = compute_pressure(virial, temperature, volume, dof);
    display_output(
        0,
        total_energy,
        potential_energy,
        kinetic_energy,
        temperature,
        virial,
        volume,
        pressure,
    );

    dcd_reporter
        .write_header(DCDHeader::new(n as u32, n_steps / out_dcd_freq + 1))
        .expect("Failed to write DCD header");
    dcd_reporter
        .write_dcd_frame(&system)
        .expect("Failed to write DCD frame");

    //**************************************
    output_time += start.elapsed();
    ////////////////////////////////////////

    for step in 1..=n_steps {
        ////////////////////////////////////////
        let start = Instant::now();
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

        if step % rem_com_freq == 0 {
            remove_v_com(&mut system);
        }

        //**************************************
        dynamics_time += start.elapsed();
        ////////////////////////////////////////

        ////////////////////////////////////////
        let start = Instant::now();
        //**************************************

        if step % out_dcd_freq == 0 {
            dcd_reporter
                .write_dcd_frame(&system)
                .expect("Failed to write DCD frame");
        }

        if step % out_ene_freq == 0 {
            potential_energy = compute_potential_energy(&system, e, s);
            kinetic_energy = compute_kinetic_energy(&system);
            total_energy = potential_energy + kinetic_energy;
            temperature = compute_temperature(kinetic_energy, dof);
            virial = compute_virial(&system, e, s);
            volume = compute_volume(&system);
            pressure = compute_pressure(virial, temperature, volume, dof);
            display_output(
                step,
                total_energy,
                potential_energy,
                kinetic_energy,
                temperature,
                virial,
                volume,
                pressure,
            );
        }

        //**************************************
        output_time += start.elapsed();
        ////////////////////////////////////////
    }

    ////////////////////////////////////////
    let start = Instant::now();
    //**************************************

    xyz_reporter
        .write_frame(&system)
        .expect("Failed to write XYZ frame");

    //**************************************
    output_time += start.elapsed();
    ////////////////////////////////////////

    // Calculate total time
    let total_time = setup_time + output_time + dynamics_time;

    // Report the results
    report_profile(setup_time, output_time, dynamics_time, total_time);
}
