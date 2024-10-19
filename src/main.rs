use byteorder::{LittleEndian, WriteBytesExt};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::fs::File;
use std::io::{self, Write};

fn compute_force(r: &[[f32; 4]], f: &mut [[f32; 3]], e: f32, s: f32, b: [f32; 3]) {
    let s2: f32 = s * s;

    for i in 0..r.len() {
        let xi: f32 = r[i][0];
        let yi: f32 = r[i][1];
        let zi: f32 = r[i][2];
        for j in (i + 1)..r.len() {
            let mut dx: f32 = r[j][0] - xi;
            let mut dy: f32 = r[j][1] - yi;
            let mut dz: f32 = r[j][2] - zi;
            dx -= b[0] * (dx / b[0]).round();
            dy -= b[1] * (dy / b[1]).round();
            dz -= b[2] * (dz / b[2]).round();

            let r2: f32 = 1.0 / (dx * dx + dy * dy + dz * dz);
            let c2: f32 = s2 * r2;
            let c4: f32 = c2 * c2;
            let c6: f32 = c4 * c2;

            let force: f32 = 48.0 * e * c6 * (c6 - 0.5) * r2;

            f[i][0] -= force * dx;
            f[i][1] -= force * dy;
            f[i][2] -= force * dz;

            f[j][0] += force * dx;
            f[j][1] += force * dy;
            f[j][2] += force * dz;
        }
    }
}

fn compute_potential_energy(r: &[[f32; 4]], e: f32, s: f32, b: [f32; 3]) -> f64 {
    let mut potential_energy: f64 = 0.0;
    let s2: f32 = s * s;

    for i in 0..r.len() {
        let xi: f32 = r[i][0];
        let yi: f32 = r[i][1];
        let zi: f32 = r[i][2];
        for rj in r.iter().skip(i + 1) {
            let mut dx: f32 = rj[0] - xi;
            let mut dy: f32 = rj[1] - yi;
            let mut dz: f32 = rj[2] - zi;
            dx -= b[0] * (dx / b[0]).round();
            dy -= b[1] * (dy / b[1]).round();
            dz -= b[2] * (dz / b[2]).round();

            let r2: f32 = 1.0 / (dx * dx + dy * dy + dz * dz);
            let c2: f32 = s2 * r2;
            let c4: f32 = c2 * c2;
            let c6: f32 = c4 * c2;

            let energy: f32 = 4.0 * e * c6 * (c6 - 1.0);
            potential_energy += energy as f64;
        }
    }

    potential_energy // kJ/mol
}

fn compute_virial(r: &[[f32; 4]], e: f32, s: f32, b: [f32; 3]) -> f64 {
    let mut total_virial: f64 = 0.0;
    let s2: f32 = s * s;

    for i in 0..r.len() {
        let xi: f32 = r[i][0];
        let yi: f32 = r[i][1];
        let zi: f32 = r[i][2];
        for rj in r.iter().skip(i + 1) {
            let mut dx: f32 = rj[0] - xi;
            let mut dy: f32 = rj[1] - yi;
            let mut dz: f32 = rj[2] - zi;
            dx -= b[0] * (dx / b[0]).round();
            dy -= b[1] * (dy / b[1]).round();
            dz -= b[2] * (dz / b[2]).round();

            let r2: f32 = 1.0 / (dx * dx + dy * dy + dz * dz);
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

fn compute_kinetic_energy(v: &[[f32; 3]], m: &[f32]) -> f64 {
    let mut kinetic_energy: f64 = 0.0;

    for i in 0..v.len() {
        let vx2: f32 = v[i][0] * v[i][0];
        let vy2: f32 = v[i][1] * v[i][1];
        let vz2: f32 = v[i][2] * v[i][2];

        let energy: f32 = 0.5_f32 * m[i] * (vx2 + vy2 + vz2);
        kinetic_energy += energy as f64;
    }

    kinetic_energy // kJ/mol
}

fn compute_temperature(kinetic_energy: f64, dof: u32) -> f64 {
    let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
    let temperature: f64 = 2.0 * kinetic_energy / (dof as f64 * boltzmann);
    temperature
}

fn compute_volume(b: [f32; 3]) -> f64 {
    (b[0] * b[1] * b[2]) as f64
}

fn compute_pressure(virial: f64, temperature: f64, volume: f64, dof: u32) -> f64 {
    let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
    let pressure: f64 = (virial + 3.0 * dof as f64 * boltzmann * temperature) / (3.0 * volume);
    pressure
}

fn remove_v_com(v: &mut [[f32; 3]], m: &[f32]) {
    let mut v_com: [f32; 3] = [0.0; 3];
    let mut m_tot: f32 = 0.0;

    for i in 0..v.len() {
        v_com[0] += m[i] * v[i][0];
        v_com[1] += m[i] * v[i][1];
        v_com[2] += m[i] * v[i][2];
        m_tot += m[i];
    }

    v_com[0] /= m_tot;
    v_com[1] /= m_tot;
    v_com[2] /= m_tot;

    for vi in v.iter_mut() {
        vi[0] -= v_com[0];
        vi[1] -= v_com[1];
        vi[2] -= v_com[2];
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

fn write_xyz_frame(file: &mut File, coordinates: &[[f32; 4]]) -> io::Result<()> {
    // Write the number of atoms (coordinates) at the top of the file
    writeln!(file, "{}", coordinates.len())?;

    // Write a comment line (can be empty or hold some metadata)
    writeln!(file, "XYZ file generated by genesis-rs")?;

    // Write each coordinate with an atom type (e.g., "C" for carbon in this example)
    for coord in coordinates {
        writeln!(file, "Ar {:.6} {:.6} {:.6}", coord[0], coord[1], coord[2])?;
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

fn write_dcd_frame(file: &mut File, r: &[[f32; 4]], b: [f32; 3]) -> io::Result<()> {
    // Block size start
    file.write_u32::<LittleEndian>(48)?;
    // Write X coordinates
    file.write_f64::<LittleEndian>(b[0] as f64)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(b[1] as f64)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(0.0)?;
    file.write_f64::<LittleEndian>(b[2] as f64)?;
    // Block size end
    file.write_u32::<LittleEndian>(48)?;

    // For each frame, DCD stores X, Y, and Z coordinates in separate chunks
    let block_size = r.len() as u32 * 4;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for coord in r {
        file.write_f32::<LittleEndian>(coord[0])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for coord in r {
        file.write_f32::<LittleEndian>(coord[1])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    // Write X coordinates
    // Block size start
    file.write_u32::<LittleEndian>(block_size)?;
    // Write X coordinates
    for coord in r {
        file.write_f32::<LittleEndian>(coord[2])?;
    }
    // Block size end
    file.write_u32::<LittleEndian>(block_size)?;

    Ok(())
}

fn main() {
    let mut out_dcd: File = File::create("out.dcd").expect("Failed to create DCD file");
    let mut out_xyz: File = File::create("out.xyz").expect("Failed to create XYZ file");

    let n: usize = 100;
    let e: f32 = 1.003;
    let s: f32 = 0.340;
    let b: [f32; 3] = [20.0; 3];
    let b_half: [f32; 3] = [b[0] * 0.5, b[1] * 0.5, b[2] * 0.5];
    let dof: u32 = 3 * n as u32 - 3;
    let m: Vec<f32> = vec![39.948; n];
    let mut r: Vec<[f32; 4]> = vec![[0.0; 4]; n];
    let mut v: Vec<[f32; 3]> = vec![[0.0; 3]; n];
    let mut f: Vec<[f32; 3]> = vec![[0.0; 3]; n];

    let mut rng: StdRng = StdRng::seed_from_u64(0);

    for i in 0..n {
        r[i][0] = rng.gen_range(-b_half[0]..b_half[0]);
        r[i][1] = rng.gen_range(-b_half[1]..b_half[1]);
        r[i][2] = rng.gen_range(-b_half[2]..b_half[2]);
        r[i][3] = rng.gen_range(-1.0..1.0);

        v[i][0] = rng.gen_range(-15.0..15.0);
        v[i][1] = rng.gen_range(-15.0..15.0);
        v[i][2] = rng.gen_range(-15.0..15.0);
    }

    let out_ene_freq: u32 = 1000;
    let out_dcd_freq: u32 = 10;
    let rem_com_freq: u32 = 10;
    assert!(
        out_ene_freq % rem_com_freq == 0,
        "out_ene_freq is not a multiple of rem_com_freq"
    );

    remove_v_com(&mut v, &m);

    display_output_header();
    let mut potential_energy: f64 = compute_potential_energy(&r, e, s, b);
    let mut kinetic_energy: f64 = compute_kinetic_energy(&v, &m);
    let mut total_energy: f64 = potential_energy + kinetic_energy;
    let mut temperature: f64 = compute_temperature(kinetic_energy, dof);
    let mut virial: f64 = compute_virial(&r, e, s, b);
    let mut volume: f64 = compute_volume(b);
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

    let dt: f32 = 0.001;
    let dt_half: f32 = dt * 0.5;
    let n_steps: u32 = 10000;

    let header = DcdHeader {
        num_atoms: n as u32,
        num_frames: n_steps / out_dcd_freq + 1,
    };

    write_dcd_header(&mut out_dcd, &header).expect("Failed to write DCD header");
    write_dcd_frame(&mut out_dcd, &r, b).expect("Failed to write DCD frame");

    for step in 1..=n_steps {
        for i in 0..n {
            r[i][0] += v[i][0] * dt_half;
            r[i][1] += v[i][1] * dt_half;
            r[i][2] += v[i][2] * dt_half;
        }

        for fi in f.iter_mut() {
            fi[0] = 0.0;
            fi[1] = 0.0;
            fi[2] = 0.0;
        }

        compute_force(&r, &mut f, e, s, b);

        for i in 0..n {
            v[i][0] += f[i][0] * dt / m[i];
            v[i][1] += f[i][1] * dt / m[i];
            v[i][2] += f[i][2] * dt / m[i];

            r[i][0] += v[i][0] * dt_half;
            r[i][1] += v[i][1] * dt_half;
            r[i][2] += v[i][2] * dt_half;
        }

        if step % rem_com_freq == 0 {
            remove_v_com(&mut v, &m);
        }

        if step % out_dcd_freq == 0 {
            write_dcd_frame(&mut out_dcd, &r, b).expect("Failed to write DCD frame");
        }

        if step % out_ene_freq == 0 {
            potential_energy = compute_potential_energy(&r, e, s, b);
            kinetic_energy = compute_kinetic_energy(&v, &m);
            total_energy = potential_energy + kinetic_energy;
            temperature = compute_temperature(kinetic_energy, dof);
            virial = compute_virial(&r, e, s, b);
            volume = compute_volume(b);
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
    }

    write_xyz_frame(&mut out_xyz, &r).expect("Failed to write XYZ file");
}
