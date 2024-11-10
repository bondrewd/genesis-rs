use crate::ff::ForceField;
use crate::system::System;

// Prelude for observer module
pub mod prelude {
    pub use crate::observer::{
        DegreesOfFreedomObserver, GeneralObserver, KineticEnergyObserver, Observer,
        PotentialEnergyObserver, VirialObserver, VolumeObserver,
    };
}

pub type GeneralObserver = Observer<
    KineticEnergyObserver,
    PotentialEnergyObserver,
    VirialObserver,
    DegreesOfFreedomObserver,
    VolumeObserver,
>;

pub trait Observable {
    type Observation;
    fn observe(&mut self, system: &System, ff: &ForceField) -> Self::Observation;
    fn reset(&mut self);
}

#[derive(Debug, Default)]
pub struct Observer<K, U, R, D, V>
where
    K: Observable + Default,
    U: Observable + Default,
    R: Observable + Default,
    D: Observable + Default,
    V: Observable + Default,
{
    pub kinetic_energy: K,
    pub potential_energy: U,
    pub virial: R,
    pub degrees_of_freedom: D,
    pub volume: V,
}

impl<K, U, R, D, V> Observer<K, U, R, D, V>
where
    K: Observable + Default,
    U: Observable + Default,
    R: Observable + Default,
    D: Observable + Default,
    V: Observable + Default,
    K::Observation: Into<f64>,
    U::Observation: Into<f64>,
    R::Observation: Into<f64>,
    D::Observation: Into<f64>,
    V::Observation: Into<f64>,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self) {
        self.kinetic_energy.reset();
        self.potential_energy.reset();
        self.virial.reset();
        self.degrees_of_freedom.reset();
        self.volume.reset();
    }

    pub fn kinetic_energy(&mut self, system: &System, ff: &ForceField) -> K::Observation {
        self.kinetic_energy.observe(system, ff)
    }

    pub fn potential_energy(&mut self, system: &System, ff: &ForceField) -> U::Observation {
        self.potential_energy.observe(system, ff)
    }

    pub fn virial(&mut self, system: &System, ff: &ForceField) -> R::Observation {
        self.virial.observe(system, ff)
    }

    pub fn degrees_of_freedom(&mut self, system: &System, ff: &ForceField) -> D::Observation {
        self.degrees_of_freedom.observe(system, ff)
    }

    pub fn volume(&mut self, system: &System, ff: &ForceField) -> V::Observation {
        self.volume.observe(system, ff)
    }

    pub fn total_energy(&mut self, system: &System, ff: &ForceField) -> f64 {
        let ke: f64 = self.kinetic_energy(system, ff).into();
        let pe: f64 = self.potential_energy(system, ff).into();
        ke + pe
    }

    pub fn temperature(&mut self, system: &System, ff: &ForceField) -> f64 {
        let ke: f64 = self.kinetic_energy(system, ff).into();
        let df: f64 = self.degrees_of_freedom(system, ff).into();
        let boltzmann: f64 = 8.314462618;
        2.0 * ke / (df as f64 * boltzmann)
    }

    pub fn pressure(&mut self, system: &System, ff: &ForceField) -> f64 {
        let vo: f64 = self.volume(system, ff).into();
        let vi: f64 = self.virial(system, ff).into();
        let te: f64 = self.temperature(system, ff);
        let df: f64 = self.degrees_of_freedom(system, ff).into();
        let boltzmann: f64 = 8.314462618;
        (vi + 3.0 * df * boltzmann * te) / (3.0 * vo)
    }
}

#[derive(Debug, Default)]
pub struct KineticEnergyObserver {
    observation: Option<f64>,
}

impl Observable for KineticEnergyObserver {
    type Observation = f64;

    fn observe(&mut self, system: &System, _ff: &ForceField) -> Self::Observation {
        match self.observation {
            Some(kinetic_energy) => kinetic_energy,
            None => {
                let mut kinetic_energy: f64 = 0.0;

                for i in 0..system.n {
                    let energy: f32 = 0.5_f32 * system.m[i] * system.v[i].dot(&system.v[i]);
                    kinetic_energy += energy as f64;
                }

                self.observation = Some(kinetic_energy);
                kinetic_energy
            }
        }
    }

    fn reset(&mut self) {
        self.observation = None;
    }
}

#[derive(Debug, Default)]
pub struct PotentialEnergyObserver {
    observation: Option<f64>,
}

impl Observable for PotentialEnergyObserver {
    type Observation = f64;

    fn observe(&mut self, system: &System, ff: &ForceField) -> Self::Observation {
        match self.observation {
            Some(potential_energy) => potential_energy,
            None => {
                let potential_energy = ff.compute_energy(system);
                self.observation = Some(potential_energy);
                potential_energy
            }
        }
    }

    fn reset(&mut self) {
        self.observation = None;
    }
}

#[derive(Debug, Default)]
pub struct VirialObserver {
    observation: Option<f64>,
}

impl Observable for VirialObserver {
    type Observation = f64;

    fn observe(&mut self, system: &System, ff: &ForceField) -> Self::Observation {
        match self.observation {
            Some(virial) => virial,
            None => {
                let virial = ff.compute_virial(system);
                self.observation = Some(virial);
                virial
            }
        }
    }

    fn reset(&mut self) {
        self.observation = None;
    }
}

#[derive(Debug, Default)]
pub struct DegreesOfFreedomObserver {
    observation: Option<u32>,
}

impl Observable for DegreesOfFreedomObserver {
    type Observation = u32;

    fn observe(&mut self, system: &System, _ff: &ForceField) -> Self::Observation {
        match self.observation {
            Some(dof) => dof,
            None => {
                let dof: u32 = 3 * system.n as u32 - 3;
                self.observation = Some(dof);
                dof
            }
        }
    }

    fn reset(&mut self) {
        self.observation = None;
    }
}

#[derive(Debug, Default)]
pub struct VolumeObserver {
    observation: Option<f64>,
}

impl Observable for VolumeObserver {
    type Observation = f64;

    fn observe(&mut self, system: &System, _ff: &ForceField) -> Self::Observation {
        match self.observation {
            Some(volume) => volume,
            None => {
                let volume = system.b.x as f64 * system.b.y as f64 * system.b.z as f64;
                self.observation = Some(volume);
                volume
            }
        }
    }

    fn reset(&mut self) {
        self.observation = None;
    }
}
