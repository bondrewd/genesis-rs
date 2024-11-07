use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct MolParserResult {
    pub names: Vec<String>,
    pub masses: Vec<f32>,
    pub charges: Vec<f32>,
    pub classes: Vec<usize>,
}

impl MolParserResult {
    pub fn new(
        names: Vec<String>,
        masses: Vec<f32>,
        charges: Vec<f32>,
        classes: Vec<usize>,
    ) -> Self {
        Self {
            names,
            masses,
            charges,
            classes,
        }
    }
}

#[derive(Debug)]
pub struct MolParser {
    pub file: File,
}

impl MolParser {
    pub fn new(file: File) -> Self {
        Self { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        Ok(Self::new(file))
    }

    pub fn parse(&self) -> Result<MolParserResult, io::Error> {
        // Initialize buffer
        let mut name_buffer = Vec::new();
        let mut mass_buffer = Vec::new();
        let mut charge_buffer = Vec::new();
        let mut class_buffer = Vec::new();

        // Read file line by line
        for line in BufReader::new(&self.file).lines() {
            let line = line?;

            // Skip comments
            if line.starts_with("#") {
                continue;
            }

            // Split line into tokens
            let mut tokens = line.split_whitespace();

            // Parse the first token: name
            name_buffer.push(
                tokens
                    .next()
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token x"))?
                    .to_owned(),
            );

            // Parse the second token: mass
            mass_buffer.push(
                tokens
                    .next()
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token y"))?
                    .parse::<f32>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            );

            // Parse the third token: charge
            charge_buffer.push(
                tokens
                    .next()
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token z"))?
                    .parse::<f32>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            );

            // Parse the fourth token: class
            class_buffer.push(
                tokens
                    .next()
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token z"))?
                    .parse::<usize>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
            );
        }

        Ok(MolParserResult::new(
            name_buffer,
            mass_buffer,
            charge_buffer,
            class_buffer,
        ))
    }
}
