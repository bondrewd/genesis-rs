use std::vec;

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Deserialize, Serialize)]
pub struct System {
    pub n: usize,
    pub b: Vector3<f32>,
    pub m: Vec<f32>,
    pub q: Vec<f32>,
    pub c: Vec<usize>,
    pub r: Vec<Vector3<f32>>,
    pub v: Vec<Vector3<f32>>,
    pub f: Vec<Vector3<f32>>,
    pub names: Vec<String>,
}

impl System {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        n: usize,
        b: Vector3<f32>,
        m: Vec<f32>,
        q: Vec<f32>,
        c: Vec<usize>,
        r: Vec<Vector3<f32>>,
        v: Vec<Vector3<f32>>,
        f: Vec<Vector3<f32>>,
        names: Vec<String>,
    ) -> Self {
        System {
            n,
            b,
            m,
            q,
            c,
            r,
            v,
            f,
            names,
        }
    }

    pub fn builder() -> SystemBuilder {
        SystemBuilder::new()
    }

    pub fn remove_v_com(&mut self) {
        let mut v_com: Vector3<f32> = Vector3::zeros();
        let mut m_tot: f32 = 0.0;

        for i in 0..self.n {
            v_com += self.m[i] * self.v[i];
            m_tot += self.m[i];
        }

        v_com /= m_tot;

        for v in self.v.iter_mut() {
            *v -= v_com;
        }
    }
}

#[derive(Debug, Default)]
pub struct SystemBuilder {
    n: Option<usize>,
    b: Option<Vector3<f32>>,
    m: Option<Vec<f32>>,
    q: Option<Vec<f32>>,
    c: Option<Vec<usize>>,
    r: Option<Vec<Vector3<f32>>>,
    v: Option<Vec<Vector3<f32>>>,
    f: Option<Vec<Vector3<f32>>>,
    names: Option<Vec<String>>,
}

impl SystemBuilder {
    pub fn new() -> Self {
        SystemBuilder::default()
    }

    pub fn with_masses(mut self, masses: Vec<f32>) -> Self {
        self.n = self.n.or(Some(masses.len()));
        self.m = Some(masses);
        self
    }

    pub fn with_charges(mut self, charges: Vec<f32>) -> Self {
        self.n = self.n.or(Some(charges.len()));
        self.q = Some(charges);
        self
    }

    pub fn with_classes(mut self, classes: Vec<usize>) -> Self {
        self.n = self.n.or(Some(classes.len()));
        self.c = Some(classes);
        self
    }

    pub fn with_cuboid_boundary(mut self, bx: f32, by: f32, bz: f32) -> Self {
        self.b = Some(Vector3::new(bx, by, bz));
        self
    }

    pub fn with_positions(mut self, positions: Vec<Vector3<f32>>) -> Self {
        self.n = self.n.or(Some(positions.len()));
        self.r = Some(positions);
        self
    }

    pub fn with_velocities(mut self, velocities: Vec<Vector3<f32>>) -> Self {
        self.n = self.n.or(Some(velocities.len()));
        self.v = Some(velocities);
        self
    }

    pub fn with_random_velocities(mut self, temperature: f64, rng: &mut StdRng) -> Self {
        let boltzmann: f64 = 8.314462618; // kJ/(mol*K)

        let m = self.m.as_ref().expect("m is required");
        let mut v = vec![Vector3::zeros(); self.n.expect("n is required")];
        for i in 0..self.n.expect("n is required") {
            let sigma: f64 = (boltzmann * temperature / m[i] as f64).sqrt();
            let normal_dist = Normal::new(0.0, sigma).unwrap();
            v[i][0] = normal_dist.sample(rng) as f32;
            v[i][1] = normal_dist.sample(rng) as f32;
            v[i][2] = normal_dist.sample(rng) as f32;
        }

        self.v = Some(v);
        self
    }

    pub fn with_forces(mut self, forces: Vec<Vector3<f32>>) -> Self {
        self.n = self.n.or(Some(forces.len()));
        self.f = Some(forces);
        self
    }

    pub fn with_names(mut self, names: Vec<String>) -> Self {
        self.n = self.n.or(Some(names.len()));
        self.names = Some(names);
        self
    }

    pub fn build(self) -> System {
        let n = self.n.expect("n is required");
        let b = self.b.expect("b is required");
        let m = self.m.expect("m is required");
        let q = self.q.expect("q is required");
        let c = self.c.expect("c is required");
        let r = self.r.expect("r is required");
        let v = self.v.expect("v is required");
        let f = self.f.unwrap_or(vec![Vector3::zeros(); n]);
        let names = self.names.unwrap_or(vec!["".to_string(); n]);

        assert!(n == m.len());
        assert!(n == q.len());
        assert!(n == r.len());
        assert!(n == v.len());
        assert!(n == f.len());
        assert!(n == names.len());
        assert!(b.x > 0.0);
        assert!(b.y > 0.0);
        assert!(b.z > 0.0);

        System::new(n, b, m, q, c, r, v, f, names)
    }
}
