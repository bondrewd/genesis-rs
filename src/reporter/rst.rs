use crate::system::System;
use serde_json;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(Debug)]
pub struct RSTReporter {
    file: File,
}

impl RSTReporter {
    pub fn new(file: File) -> Self {
        RSTReporter { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(RSTReporter::new(file))
    }

    pub fn write_report(&mut self, system: &System) -> io::Result<()> {
        // Serialize the system to JSON and write it to the file
        let serialized_system = serde_json::to_string(system)?;
        self.file.write_all(serialized_system.as_bytes())?;
        Ok(())
    }
}
