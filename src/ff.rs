use nalgebra::DMatrix;

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
pub struct ForceField {
    pub lj: DMatrix<LennardJonesItem>,
}

impl ForceField {
    pub fn new(lj: DMatrix<LennardJonesItem>) -> Self {
        Self { lj }
    }
}
