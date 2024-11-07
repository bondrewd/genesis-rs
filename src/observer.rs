use crate::ff::ForceField;
use crate::system::System;

#[derive(Debug, Default)]
pub struct PotentialEnergyObserver {
    observation: Option<f64>,
}

impl PotentialEnergyObserver {
    pub fn new() -> Self {
        PotentialEnergyObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, system: &System, ff: &ForceField) {
        let mut potential_energy: f64 = 0.0;

        for i in 0..system.n {
            let ri = system.r[i];
            let ci = system.c[i];
            for j in (i + 1)..system.n {
                let cj = system.c[j];
                let lj = ff.lj[&(ci, cj)];
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

        self.observation = Some(potential_energy); // kJ/mol
    }
}

#[derive(Debug, Default)]
pub struct TotalEnergyObserver {
    observation: Option<f64>,
}

impl TotalEnergyObserver {
    pub fn new() -> Self {
        TotalEnergyObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, ke_obs: &KineticEnergyObserver, ue_obs: &PotentialEnergyObserver) {
        let ke = ke_obs.last_observation();
        let ue = ue_obs.last_observation();
        if let (Some(ke), Some(ue)) = (ke, ue) {
            let total_energy = ke + ue;
            self.observation = Some(total_energy);
        } else {
            self.observation = None;
        }
    }
}

#[derive(Debug, Default)]
pub struct VirialObserver {
    observation: Option<f64>,
}

impl VirialObserver {
    pub fn new() -> Self {
        VirialObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, system: &System, ff: &ForceField) {
        let mut total_virial: f64 = 0.0;

        for i in 0..system.n {
            let ri = system.r[i];
            let ci = system.c[i];
            for j in (i + 1)..system.n {
                let cj = system.c[j];
                let lj = ff.lj[&(ci, cj)];
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

        self.observation = Some(total_virial); // kJ/mol
    }
}

#[derive(Debug, Default)]
pub struct KineticEnergyObserver {
    observation: Option<f64>,
}

impl KineticEnergyObserver {
    pub fn new() -> Self {
        KineticEnergyObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, system: &System) {
        let mut kinetic_energy: f64 = 0.0;

        for i in 0..system.n {
            let energy: f32 = 0.5_f32 * system.m[i] * system.v[i].dot(&system.v[i]);
            kinetic_energy += energy as f64;
        }

        self.observation = Some(kinetic_energy); // kJ/mol
    }
}

#[derive(Debug, Default)]
pub struct DegreesOfFreedomObserver {
    observation: Option<u32>,
}

impl DegreesOfFreedomObserver {
    pub fn new() -> Self {
        DegreesOfFreedomObserver::default()
    }

    pub fn last_observation(&self) -> Option<u32> {
        self.observation
    }

    pub fn observe(&mut self, system: &System) {
        let dof: u32 = 3 * system.n as u32 - 3;
        self.observation = Some(dof);
    }
}

#[derive(Debug, Default)]
pub struct TemperatureObserver {
    observation: Option<f64>,
}

impl TemperatureObserver {
    pub fn new() -> Self {
        TemperatureObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, ke_obs: &KineticEnergyObserver, df_obs: &DegreesOfFreedomObserver) {
        let ke = ke_obs.last_observation();
        let df = df_obs.last_observation();
        if let (Some(ke), Some(df)) = (ke, df) {
            let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
            let temperature = 2.0 * ke / (df as f64 * boltzmann);
            self.observation = Some(temperature);
        } else {
            self.observation = None;
        }
    }
}

#[derive(Debug, Default)]
pub struct VolumeObserver {
    observation: Option<f64>,
}

impl VolumeObserver {
    pub fn new() -> Self {
        VolumeObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(&mut self, system: &System) {
        let volume: f64 = system.b.x as f64 * system.b.y as f64 * system.b.z as f64;
        self.observation = Some(volume);
    }
}

#[derive(Debug, Default)]
pub struct PressureObserver {
    observation: Option<f64>,
}

impl PressureObserver {
    pub fn new() -> Self {
        PressureObserver::default()
    }

    pub fn last_observation(&self) -> Option<f64> {
        self.observation
    }

    pub fn observe(
        &mut self,
        vo_obs: &VolumeObserver,
        vi_obs: &VirialObserver,
        te_obs: &TemperatureObserver,
        df_obs: &DegreesOfFreedomObserver,
    ) {
        let vo = vo_obs.last_observation();
        let vi = vi_obs.last_observation();
        let te = te_obs.last_observation();
        let df = df_obs.last_observation();
        if let (Some(vo), Some(vi), Some(te), Some(df)) = (vo, vi, te, df) {
            let boltzmann: f64 = 8.314462618; // kJ/(mol*K)
            let pressure: f64 = (vi + 3.0 * df as f64 * boltzmann * te) / (3.0 * vo);
            self.observation = Some(pressure);
        } else {
            self.observation = None;
        }
    }
}
