use dirs::home_dir;
use serde::Deserialize;
use std::path::{Path, PathBuf};

#[derive(Debug, Default, Deserialize)]
pub struct InputConfig {
    pub par_path: Option<PathBuf>,
    pub mol_path: Option<PathBuf>,
    pub pos_path: Option<PathBuf>,
    pub vel_path: Option<PathBuf>,
    pub rst_path: Option<PathBuf>,
}

#[derive(Debug, Default, Deserialize)]
pub struct OutputConfig {
    pub csv_path: Option<PathBuf>,
    pub dcd_path: Option<PathBuf>,
    pub log_path: Option<PathBuf>,
    pub rst_path: Option<PathBuf>,
    pub xyz_path: Option<PathBuf>,
    pub csv_freq: u32,
    pub dcd_freq: u32,
    pub rst_freq: u32,
}

#[derive(Debug, Default, Deserialize)]
pub struct DynamicsConfig {
    pub time_step: f32,
    pub num_steps: u32,
    pub temperature: f64,
    pub remove_com_v_freq: u32,
}

#[derive(Debug, Default, Deserialize)]
pub struct BoundaryConfig {
    pub length_x: f32,
    pub length_y: f32,
    pub length_z: f32,
}

#[derive(Debug, Default, Deserialize)]
pub struct RngConfig {
    pub seed: u64,
}

#[derive(Debug, Default, Deserialize)]
pub struct Config {
    pub input: InputConfig,
    pub output: OutputConfig,
    pub dynamics: DynamicsConfig,
    pub boundary: BoundaryConfig,
    pub rng: RngConfig,
}

impl TryFrom<&Path> for Config {
    type Error = Box<dyn std::error::Error>;

    fn try_from(path: &Path) -> Result<Self, Self::Error> {
        let parent = path.parent().map(|p| p.to_path_buf()).unwrap();
        let config = std::fs::read_to_string(path)?;
        let mut config: Config = toml::from_str(&config)?;

        config.input.par_path = config.input.par_path.map(|p| resolve_path(&p, &parent));
        config.input.mol_path = config.input.mol_path.map(|p| resolve_path(&p, &parent));
        config.input.pos_path = config.input.pos_path.map(|p| resolve_path(&p, &parent));
        config.input.vel_path = config.input.vel_path.map(|p| resolve_path(&p, &parent));
        config.input.rst_path = config.input.rst_path.map(|p| resolve_path(&p, &parent));

        config.output.csv_path = config.output.csv_path.map(|p| resolve_path(&p, &parent));
        config.output.dcd_path = config.output.dcd_path.map(|p| resolve_path(&p, &parent));
        config.output.log_path = config.output.log_path.map(|p| resolve_path(&p, &parent));
        config.output.rst_path = config.output.rst_path.map(|p| resolve_path(&p, &parent));
        config.output.xyz_path = config.output.xyz_path.map(|p| resolve_path(&p, &parent));

        Ok(config)
    }
}

pub fn resolve_path<P: AsRef<Path>>(path: P, parent: P) -> PathBuf {
    let path = PathBuf::from(path.as_ref());
    let path = path
        .strip_prefix("~")
        .map(|p| home_dir().map(|home| home.join(p)))
        .ok()
        .flatten()
        .unwrap_or(path);
    if path.is_absolute() {
        path
    } else {
        parent.as_ref().join(path)
    }
}
