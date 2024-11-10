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

        // Resolve paths
        let resolver = |p: PathBuf| -> PathBuf {
            if p.is_absolute() {
                p
            } else {
                parent.join(p)
            }
        };

        config.input.par_path = config.input.par_path.map(resolver);
        config.input.mol_path = config.input.mol_path.map(resolver);
        config.input.pos_path = config.input.pos_path.map(resolver);
        config.input.vel_path = config.input.vel_path.map(resolver);
        config.input.rst_path = config.input.rst_path.map(resolver);

        config.output.csv_path = config.output.csv_path.map(resolver);
        config.output.dcd_path = config.output.dcd_path.map(resolver);
        config.output.log_path = config.output.log_path.map(resolver);
        config.output.rst_path = config.output.rst_path.map(resolver);
        config.output.xyz_path = config.output.xyz_path.map(resolver);

        // Check output frequencies
        if config.dynamics.num_steps % config.output.csv_freq != 0 {
            return Err("num_steps in [dynamics] is not a multiple of csv_freq in [output]".into());
        }

        if config.dynamics.num_steps % config.output.dcd_freq != 0 {
            return Err("num_steps in [dynamics] is not a multiple of dcd_freq in [output]".into());
        }

        if config.dynamics.num_steps % config.output.rst_freq != 0 {
            return Err("num_steps in [dynamics] is not a multiple of rst_freq in [output]".into());
        }

        Ok(config)
    }
}
