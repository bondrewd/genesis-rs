use crate::parser::par::ParParserResult;
use crate::system::System;
use nalgebra::DMatrix;

// Prelude for ff module
pub mod prelude {
    pub use crate::ff::{ForceField, LennardJones, LennardJonesItem};
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct LennardJonesItem {
    pub epsilon: f32,
    pub sigma: f32,
}

impl LennardJonesItem {
    pub fn new(epsilon: f32, sigma: f32) -> Self {
        Self { epsilon, sigma }
    }
}

#[derive(Debug)]
pub struct LennardJones {
    pub parameters: DMatrix<LennardJonesItem>,
}

impl LennardJones {
    pub fn new(parameters: DMatrix<LennardJonesItem>) -> Self {
        Self { parameters }
    }

    pub fn update_force(&self, system: &mut System) {
        for i in 0..system.n {
            let ri = system.r[i];
            let ci = system.c[i];
            for j in (i + 1)..system.n {
                let cj = system.c[j];
                let lj = self.parameters[(ci, cj)];
                let e = lj.epsilon;
                let s = lj.sigma;
                let s2 = s * s;

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

    pub fn compute_energy(&self, system: &System) -> f64 {
        let mut potential_energy: f64 = 0.0;

        for i in 0..system.n {
            let ri = system.r[i];
            let ci = system.c[i];
            for j in (i + 1)..system.n {
                let cj = system.c[j];
                let lj = self.parameters[(ci, cj)];
                let e = lj.epsilon;
                let s = lj.sigma;
                let s2 = s * s;

                let mut dr = system.r[j] - ri;
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

        potential_energy
    }

    pub fn compute_virial(&self, system: &System) -> f64 {
        let mut total_virial: f64 = 0.0;

        for i in 0..system.n {
            let ri = system.r[i];
            let ci = system.c[i];
            for j in (i + 1)..system.n {
                let cj = system.c[j];
                let lj = self.parameters[(ci, cj)];
                let e = lj.epsilon;
                let s = lj.sigma;
                let s2 = s * s;

                let mut dr = system.r[j] - ri;
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
}

#[derive(Debug, Default)]
pub struct ForceField {
    pub lennard_jones: Option<LennardJones>,
}

impl ForceField {
    pub fn new(lennard_jones: Option<LennardJones>) -> Self {
        Self { lennard_jones }
    }

    pub fn with_parameters(parameters: ParParserResult) -> Self {
        let mut ff = ForceField::default();
        if let Some(lj) = parameters.lennard_jones {
            ff.lennard_jones = Some(LennardJones::new(lj));
        }

        ff
    }

    pub fn update_force(&self, system: &mut System) {
        if let Some(lj) = &self.lennard_jones {
            lj.update_force(system);
        }
    }

    pub fn compute_energy(&self, system: &System) -> f64 {
        if let Some(lj) = &self.lennard_jones {
            lj.compute_energy(system)
        } else {
            0.0
        }
    }

    pub fn compute_virial(&self, system: &System) -> f64 {
        if let Some(lj) = &self.lennard_jones {
            lj.compute_virial(system)
        } else {
            0.0
        }
    }
}
