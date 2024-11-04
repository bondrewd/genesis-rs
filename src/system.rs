use std::vec;

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::Rng;
use rand_distr::{Distribution, Normal};

#[derive(Debug)]
pub struct System {
    pub n: usize,
    pub b: Vector3<f32>,
    pub m: Vec<f32>,
    pub q: Vec<f32>,
    pub r: Vec<Vector3<f32>>,
    pub v: Vec<Vector3<f32>>,
    pub f: Vec<Vector3<f32>>,
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

    pub fn with_cuboid_boundary(mut self, bx: f32, by: f32, bz: f32) -> Self {
        self.b = Some(Vector3::new(bx, by, bz));
        self
    }

    pub fn with_random_positions(mut self, rng: &mut StdRng) -> Self {
        let n = self.n.expect("n is required");
        let b = self.b.expect("b is required");
        let r: Vec<Vector3<f32>> = (0..n)
            .map(|_| {
                Vector3::new(
                    rng.gen_range(0.0..b.x),
                    rng.gen_range(0.0..b.y),
                    rng.gen_range(0.0..b.z),
                )
            })
            .collect();

        self.r = Some(r);
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

    pub fn with_random_charges(mut self, rng: &mut StdRng) -> Self {
        let n = self.n.expect("n is required");
        let mut q = vec![0.0; n];
        for qi in q.iter_mut() {
            *qi = rng.gen_range(-1.0..1.0);
        }

        self.q = Some(q);
        self
    }

    pub fn build(self) -> System {
        let n = self.n.expect("n is required");
        let b = self.b.expect("b is required");
        let m = self.m.expect("m is required");
        let q = self.q.expect("q is required");
        let r = self.r.expect("r is required");
        let v = self.v.expect("v is required");
        let f = self.f.unwrap_or(vec![Vector3::zeros(); n]);

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
