use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct LennardJonesItem {
    pub epsilon: f32,
    pub sigma: f32,
}

impl LennardJonesItem {
    pub fn new(epsilon: f32, sigma: f32) -> Self {
        Self { epsilon, sigma }
    }
}

pub type LennardJonesParameters = HashMap<(usize, usize), LennardJonesItem>;

#[derive(Debug)]
pub struct ForceField {
    pub lj: LennardJonesParameters,
}

impl ForceField {
    pub fn new(lj: LennardJonesParameters) -> Self {
        Self { lj }
    }
}
