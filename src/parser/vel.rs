use nalgebra::Vector3;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct VelParser {
    pub file: File,
}

impl VelParser {
    pub fn new(file: File) -> Self {
        Self { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        Ok(Self::new(file))
    }

    pub fn parse(&self) -> Result<Vec<Vector3<f32>>, io::Error> {
        // Initialize buffer
        let mut buffer = Vec::new();

        // Read file line by line
        for line in BufReader::new(&self.file).lines() {
            let line = line?;

            // Skip comments
            if line.starts_with("#") {
                continue;
            }

            // Split line into tokens
            let mut tokens = line.split_whitespace();

            // Parse the first token: vx
            let x: f32 = tokens
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token x"))?
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            // Parse the second token: vy
            let y: f32 = tokens
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token y"))?
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            // Parse the third token: vz
            let z: f32 = tokens
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing token z"))?
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            // Push the parsed values to the buffer
            buffer.push(Vector3::new(x, y, z));
        }

        Ok(buffer)
    }
}
