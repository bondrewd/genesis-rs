use crate::ff::ForceField;
use crate::observer::GeneralObserver;
use crate::system::System;
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

    pub fn write_report(
        &mut self,
        step: u32,
        observer: &mut GeneralObserver,
        system: &System,
        ff: &ForceField,
    ) -> io::Result<()> {
        // Init buffer
        let mut buffer = Vec::new();

        // Write the observations
        buffer.extend_from_slice(
            format!(
                "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                step,
                observer.total_energy(system, ff),
                observer.potential_energy(system, ff),
                observer.kinetic_energy(system, ff),
                observer.temperature(system, ff),
                observer.virial(system, ff),
                observer.volume(system, ff),
                observer.pressure(system, ff),
            )
            .as_bytes(),
        );

        // Flush buffer
        self.file.write_all(&buffer)?;

        Ok(())
    }
}
