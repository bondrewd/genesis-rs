use byteorder::{LittleEndian, WriteBytesExt};
use clap::Parser;
use dirs::home_dir;
use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use serde::Deserialize;
use std::fs;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

#[derive(Debug)]
struct System {
    n: usize,
    b: Vector3<f32>,
    m: Vec<f32>,
    q: Vec<f32>,
    r: Vec<Vector3<f32>>,
    v: Vec<Vector3<f32>>,
    f: Vec<Vector3<f32>>,
}

impl System {
    pub fn new(
        n: usize,
        b: Vector3<f32>,
        m: Vec<f32>,
        q: Vec<f32>,
        r: Vec<Vector3<f32>>,
        v: Vec<Vector3<f32>>,
        f: Vec<Vector3<f32>>,
    ) -> Self {
        System {
            n,
            b,
            m,
            q,
            r,
            v,
            f,
        }
    }
}

#[derive(Debug, Default)]
struct SystemBuilder {
    n: Option<usize>,
    b: Option<Vector3<f32>>,
    m: Option<Vec<f32>>,
    q: Option<Vec<f32>>,
    r: Option<Vec<Vector3<f32>>>,
    v: Option<Vec<Vector3<f32>>>,
    f: Option<Vec<Vector3<f32>>>,
}

impl SystemBuilder {
    pub fn new() -> Self {
        SystemBuilder::default()
    }

    pub fn n(mut self, n: usize) -> Self {
        self.n = Some(n);
        self
    }

    pub fn b(mut self, b: Vector3<f32>) -> Self {
        self.b = Some(b);
        self
    }

    pub fn m(mut self, m: Vec<f32>) -> Self {
        self.m = Some(m);
        self
    }

    pub fn q(mut self, q: Vec<f32>) -> Self {
        self.q = Some(q);
        self
    }

    pub fn r(mut self, r: Vec<Vector3<f32>>) -> Self {
        self.r = Some(r);
        self
    }

    pub fn v(mut self, v: Vec<Vector3<f32>>) -> Self {
        self.v = Some(v);
        self
    }

    pub fn f(mut self, f: Vec<Vector3<f32>>) -> Self {
        self.f = Some(f);
        self
    }

    pub fn build(self) -> System {
        let n = self.n.expect("n is required");
        let b = self.b.expect("b is required");
        let m = self.m.expect("m is required");
        let q = self.q.expect("q is required");
        let r = self.r.expect("r is required");
        let v = self.v.expect("v is required");
        let f = self.f.expect("f is required");

        assert!(n == m.len());
        assert!(n == q.len());
        assert!(n == r.len());
        assert!(n == v.len());
        assert!(n == f.len());
        assert!(b.x > 0.0);
        assert!(b.y > 0.0);
        assert!(b.z > 0.0);

        System::new(n, b, m, q, r, v, f)
    }
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

fn initialize_velocities(system: &mut System, temperature: f64, rng: &mut StdRng) {
    let boltzmann: f64 = 8.314462618; // kJ/(mol*K)

    for i in 0..system.n {
        let sigma: f64 = (boltzmann * temperature / system.m[i] as f64).sqrt();
        let normal_dist = Normal::new(0.0, sigma).unwrap();
        system.v[i][0] = normal_dist.sample(rng) as f32;
        system.v[i][1] = normal_dist.sample(rng) as f32;
        system.v[i][2] = normal_dist.sample(rng) as f32;
    }
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

fn write_xyz_frame(file: &mut File, system: &System) -> io::Result<()> {
    // Write the number of atoms (coordinates) at the top of the file
    writeln!(file, "{}", system.n)?;

    // Write a comment line (can be empty or hold some metadata)
    writeln!(file, "XYZ file generated by genesis-rs")?;

    // Write each coordinate with an atom type (e.g., "C" for carbon in this example)
    for r in system.r.iter() {
        writeln!(file, "Ar {:.6} {:.6} {:.6}", r[0], r[1], r[2])?;
    }

    Ok(())
}

#[derive(Debug)]
struct DcdHeader {
    num_atoms: u32,
    num_frames: u32,
}

fn write_dcd_header(file: &mut File, header: &DcdHeader) -> io::Result<()> {
    // Block 1
    // Block size start
    file.write_u32::<LittleEndian>(84)?;
    // 01 - 04: "CORD" magic number
    file.write_all(b"CORD")?;
    // 05 - 08: Number of frames
    file.write_u32::<LittleEndian>(header.num_frames)?;
    // 09 - 12: Unused (first step)
    file.write_u32::<LittleEndian>(0)?;
    // 13 - 16: Unused (output period)
    file.write_u32::<LittleEndian>(0)?;
    // 17 - 20: Unused (number of time steps)
    file.write_u32::<LittleEndian>(0)?;
    // 21 - 40: Unused
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    // 41 - 44: Unused (time step)
    file.write_f32::<LittleEndian>(0.0)?;
    // 45 - 48: Unit cell flag
    file.write_u32::<LittleEndian>(1)?;
    // 49 - 80: Unused
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u32::<LittleEndian>(0)?;
    // 81 - 84: Version
    file.write_u32::<LittleEndian>(24)?;
    // Block size end
    file.write_u32::<LittleEndian>(84)?;

    // Block 2
    // Block size start
    file.write_u32::<LittleEndian>(84)?;
    // 01 - 04: Number of lines
    file.write_u32::<LittleEndian>(1)?;
    // 05 - 84: Comment
    file.write_all(b"                ")?;
    file.write_all(b"                ")?;
    file.write_all(b"                ")?;
    file.write_all(b"                ")?;
    file.write_all(b"                ")?;
    // Block size end
    file.write_u32::<LittleEndian>(84)?;

    // Block 3
    // Block size start
    file.write_u32::<LittleEndian>(4)?;
    // 01 - 04: Number of particles
    file.write_u32::<LittleEndian>(header.num_atoms)?;
    // Block size end
    file.write_u32::<LittleEndian>(4)?;

    Ok(())
}

fn write_dcd_frame(file: &mut File, system: &System) -> io::Result<()> {
    // Block size start
    file.write_u32::<LittleEndian>(48)?;
    // Write X coordinates
    file.write_f64::<LittleEndian>(system.b.x as f64)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(system.b.y as f64)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(system.b.z as f64)?;
    // Block size end
    file.write_u32::<LittleEndian>(48)?;

    // For each frame, DCD stores X, Y, and Z coordinates in separate chunks
    let block_size = system.n as u32 * 4;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for r in system.r.iter() {
        file.write_f32::<LittleEndian>(r[0])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for r in system.r.iter() {
        file.write_f32::<LittleEndian>(r[1])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for r in system.r.iter() {
        file.write_f32::<LittleEndian>(r[2])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    Ok(())
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

    // Initialize the output files
    let mut out_dcd_file: File = File::create(out_dcd_path).expect("Failed to create DCD file");
    let mut out_xyz_file: File = File::create(out_xyz_path).expect("Failed to create XYZ file");

    let n: usize = 100;
    let e: f32 = 1.003;
    let s: f32 = 0.340;
    let mut system = SystemBuilder::new()
        .n(n)
        .b(Vector3::new(
            config.boundary.length_x,
            config.boundary.length_y,
            config.boundary.length_z,
        ))
        .m(vec![39.948; n])
        .q(vec![0.0; n])
        .r(vec![Vector3::zeros(); n])
        .v(vec![Vector3::zeros(); n])
        .f(vec![Vector3::zeros(); n])
        .build();
    let dof: u32 = 3 * system.n as u32 - 3;

    let mut rng: StdRng = StdRng::seed_from_u64(0);

    for r in system.r.iter_mut() {
        r[0] = rng.gen_range(0.0..system.b.x);
        r[1] = rng.gen_range(0.0..system.b.y);
        r[2] = rng.gen_range(0.0..system.b.z);
    }

    for q in system.q.iter_mut() {
        *q = rng.gen_range(-1.0..1.0);
    }

    let target_temperature: f64 = config.dynamics.temperature;
    initialize_velocities(&mut system, target_temperature, &mut rng);
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

    let header = DcdHeader {
        num_atoms: n as u32,
        num_frames: n_steps / out_dcd_freq + 1,
    };

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

    write_dcd_header(&mut out_dcd_file, &header).expect("Failed to write DCD header");
    write_dcd_frame(&mut out_dcd_file, &system).expect("Failed to write DCD frame");

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
            write_dcd_frame(&mut out_dcd_file, &system).expect("Failed to write DCD frame");
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

    write_xyz_frame(&mut out_xyz_file, &system).expect("Failed to write XYZ file");

    //**************************************
    output_time += start.elapsed();
    ////////////////////////////////////////

    // Calculate total time
    let total_time = setup_time + output_time + dynamics_time;

    // Report the results
    report_profile(setup_time, output_time, dynamics_time, total_time);
}
