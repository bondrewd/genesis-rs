use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn compute_force(r: &Vec<[f32; 4]>, f: &mut Vec<[f32; 3]>, e: f32, s: f32) {
    let s2 = s * s;

    for i in 0..r.len() {
        f[i] = [0.0, 0.0, 0.0];
        for j in 0..r.len() {
            if i != j {
                let dr = [r[j][0] - r[i][0], r[j][1] - r[i][1], r[j][2] - r[i][2]];
                let r2 = 1.0_f32 / (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
                let c2 = s2 * r2;
                let c4 = c2 * c2;
                let c8 = c4 * c4;

                let force = 48.0 * e * c8 * c4 * (c2 - 0.5) * r2;

                f[i][0] -= force * dr[0];
                f[i][1] -= force * dr[1];
                f[i][2] -= force * dr[2];

                f[j][0] += force * dr[0];
                f[j][1] += force * dr[1];
                f[j][2] += force * dr[2];
            }
        }
    }
}

fn main() {
    let n = 100;
    let e = 1.0_f32;
    let s = 1.0_f32;
    let m: Vec<f32> = vec![1.0_f32; n];
    let mut r: Vec<[f32; 4]> = vec![[0.0_f32; 4]; n];
    let mut v: Vec<[f32; 3]> = vec![[0.0_f32; 3]; n];
    let mut f: Vec<[f32; 3]> = vec![[0.0_f32; 3]; n];

    let mut rng = StdRng::seed_from_u64(0);
    for i in 0..n {
        r[i][0] = rng.gen_range(0.0..10.0);
        r[i][1] = rng.gen_range(0.0..10.0);
        r[i][2] = rng.gen_range(0.0..10.0);
        r[i][3] = rng.gen_range(-1.0..1.0);

        v[i][0] = rng.gen_range(-1.0..1.0);
        v[i][1] = rng.gen_range(-1.0..1.0);
        v[i][2] = rng.gen_range(-1.0..1.0);
    }

    let dt = 0.01_f32;
    let dt_half = dt * 0.5;
    let n_steps = 1000;
    for _ in 0..n_steps {
        for i in 0..n {
            r[i][0] += v[i][0] * dt_half;
            r[i][1] += v[i][1] * dt_half;
            r[i][2] += v[i][2] * dt_half;
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
    }

    println!("r = {:?}", r);
    println!("v = {:?}", v);
    println!("f = {:?}", f);
}
