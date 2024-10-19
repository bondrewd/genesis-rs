use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn compute_force(r: &Vec<[f32; 4]>, f: &mut Vec<[f32; 3]>, e: f32, s: f32) {
    let s2: f32 = s * s;

    for i in 0..r.len() {
        let xi = r[i][0];
        let yi = r[i][1];
        let zi = r[i][2];
        for j in (i + 1)..r.len() {
            let dx = r[j][0] - xi;
            let dy = r[j][1] - yi;
            let dz = r[j][2] - zi;
            let r2 = 1.0_f32 / (dx * dx + dy * dy + dz * dz);
            let c2 = s2 * r2;
            let c4 = c2 * c2;
            let c6 = c4 * c2;

            let force = 48.0_f32 * e * c6 * (c6 - 0.5_f32) * r2;

            f[i][0] -= force * dx;
            f[i][1] -= force * dy;
            f[i][2] -= force * dz;

            f[j][0] += force * dx;
            f[j][1] += force * dy;
            f[j][2] += force * dz;
        }
    }
}

fn compute_potential_energy(r: &Vec<[f32; 4]>, e: f32, s: f32) -> f64 {
    let mut potential_energy: f64 = 0.0;
    let s2: f32 = s * s;

    for i in 0..r.len() {
        let xi = r[i][0];
        let yi = r[i][1];
        let zi = r[i][2];
        for j in (i + 1)..r.len() {
            let dx = r[j][0] - xi;
            let dy = r[j][1] - yi;
            let dz = r[j][2] - zi;
            let r2 = 1.0_f32 / (dx * dx + dy * dy + dz * dz);
            let c2 = s2 * r2;
            let c4 = c2 * c2;
            let c6 = c4 * c2;

            let energy: f32 = 4.0_f32 * e * c6 * (c6 - 1.0_f32);
            potential_energy += energy as f64;
        }
    }

    potential_energy
}

fn compute_kinetic_energy(v: &Vec<[f32; 3]>, m: &Vec<f32>) -> f64 {
    let mut kinetic_energy: f64 = 0.0;

    for i in 0..v.len() {
        let vx2 = v[i][0] * v[i][0];
        let vy2 = v[i][1] * v[i][1];
        let vz2 = v[i][2] * v[i][2];

        let energy: f32 = 0.5_f32 * m[i] * (vx2 + vy2 + vz2);
        kinetic_energy += energy as f64;
    }

    kinetic_energy
}

fn display_output_header() {
    println!(
        "{:>8} {:>14} {:>14} {:>14}",
        "Step", "Total", "Potential", "Kinetic"
    );
}

fn display_output(step: u32, total_energy: f64, potential_energy: f64, kinetic_energy: f64) {
    println!(
        "{:>8} {:>14.6} {:>14.6} {:>14.6}",
        step, total_energy, potential_energy, kinetic_energy,
    );
}

fn main() {
    let n: usize = 100;
    let e: f32 = 1.0;
    let s: f32 = 1.0;
    let m: Vec<f32> = vec![1.0_f32; n];
    let mut r: Vec<[f32; 4]> = vec![[0.0_f32; 4]; n];
    let mut v: Vec<[f32; 3]> = vec![[0.0_f32; 3]; n];
    let mut f: Vec<[f32; 3]> = vec![[0.0_f32; 3]; n];

    let mut rng = StdRng::seed_from_u64(0);

    for i in 0..n {
        r[i][0] = rng.gen_range(-5.0..5.0);
        r[i][1] = rng.gen_range(-5.0..5.0);
        r[i][2] = rng.gen_range(-5.0..5.0);
        r[i][3] = rng.gen_range(-1.0..1.0);

        v[i][0] = rng.gen_range(-0.5..0.5);
        v[i][1] = rng.gen_range(-0.5..0.5);
        v[i][2] = rng.gen_range(-0.5..0.5);
    }

    display_output_header();
    let mut potential_energy: f64 = compute_potential_energy(&r, e, s);
    let mut kinetic_energy: f64 = compute_kinetic_energy(&v, &m);
    let mut total_energy: f64 = potential_energy + kinetic_energy;
    display_output(0, total_energy, potential_energy, kinetic_energy);

    let dt: f32 = 0.00001;
    let dt_half: f32 = dt * 0.5_f32;
    let n_steps: u32 = 1000;
    let ene_freq: u32 = 100;
    for step in 1..=n_steps {
        for i in 0..n {
            r[i][0] += v[i][0] * dt_half;
            r[i][1] += v[i][1] * dt_half;
            r[i][2] += v[i][2] * dt_half;
        }

        for i in 0..n {
            f[i][0] = 0.0_f32;
            f[i][1] = 0.0_f32;
            f[i][2] = 0.0_f32;
        }

        compute_force(&r, &mut f, e, s);

        for i in 0..n {
            v[i][0] += f[i][0] * dt / m[i];
            v[i][1] += f[i][1] * dt / m[i];
            v[i][2] += f[i][2] * dt / m[i];

            r[i][0] += v[i][0] * dt_half;
            r[i][1] += v[i][1] * dt_half;
            r[i][2] += v[i][2] * dt_half;
        }

        if step % ene_freq == 0 {
            potential_energy = compute_potential_energy(&r, e, s);
            kinetic_energy = compute_kinetic_energy(&v, &m);
            total_energy = potential_energy + kinetic_energy;
            display_output(step, total_energy, potential_energy, kinetic_energy);
        }
    }
}
