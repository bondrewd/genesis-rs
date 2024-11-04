use crate::observer::{
    KineticEnergyObserver, PotentialEnergyObserver, PressureObserver, TemperatureObserver,
    TotalEnergyObserver, VirialObserver, VolumeObserver,
};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(Debug)]
pub struct CSVReporter {
    file: File,
}

impl CSVReporter {
    pub fn new(file: File) -> Self {
        CSVReporter { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(CSVReporter::new(file))
    }

    pub fn write_header(&mut self) -> io::Result<()> {
        // Init buffer
        let mut buffer = Vec::new();

        // Write the header
        buffer.extend_from_slice(
            format!(
                "{},{},{},{},{},{},{},{}\n",
                "Step",
                "TotalEnergy",
                "PotentialEnergy",
                "KineticEnergy",
                "Temperature",
                "Virial",
                "Volume",
                "Pressure",
            )
            .as_bytes(),
        );

        // Flush buffer
        self.file.write_all(&buffer)?;

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn write_report(
        &mut self,
        step: u32,
        et_obs: &TotalEnergyObserver,
        ue_obs: &PotentialEnergyObserver,
        ke_obs: &KineticEnergyObserver,
        te_obs: &TemperatureObserver,
        vi_obs: &VirialObserver,
        vo_obs: &VolumeObserver,
        pr_obs: &PressureObserver,
    ) -> io::Result<()> {
        // Init buffer
        let mut buffer = Vec::new();

        // Write the observations
        buffer.extend_from_slice(
            format!(
                "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                step,
                et_obs.last_observation().unwrap(),
                ue_obs.last_observation().unwrap(),
                ke_obs.last_observation().unwrap(),
                te_obs.last_observation().unwrap(),
                vi_obs.last_observation().unwrap(),
                vo_obs.last_observation().unwrap(),
                pr_obs.last_observation().unwrap(),
            )
            .as_bytes(),
        );

        // Flush buffer
        self.file.write_all(&buffer)?;

        Ok(())
    }
}
