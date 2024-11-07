use crate::ff::{LennardJonesItem, LennardJonesParameters};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct ParParserResult {
    pub lj: LennardJonesParameters,
}

impl ParParserResult {
    pub fn new(lj: LennardJonesParameters) -> Self {
        Self { lj }
    }
}

#[derive(Debug)]
enum ParLineKind {
    LennardJones,
}

impl ParLineKind {
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "LennardJones" | "lj" => Some(Self::LennardJones),
            _ => None,
        }
    }
}

#[derive(Debug)]
pub struct ParParser {
    pub file: File,
}

impl ParParser {
    pub fn new(file: File) -> Self {
        Self { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        Ok(Self::new(file))
    }

    pub fn parse(&self) -> Result<ParParserResult, io::Error> {
        // Initialize buffer
        let mut lj_parameters = HashMap::new();

        // Read file line by line
        for line in BufReader::new(&self.file).lines() {
            let line = line?;

            // Skip comments
            if line.starts_with("#") {
                continue;
            }

            // Split line into tokens
            let mut tokens = line.split_whitespace();

            // Skip empty lines
            if tokens.clone().count() == 0 {
                continue;
            }

            // Parse line based on the first token
            match ParLineKind::from_str(tokens.next().unwrap_or("")) {
                Some(ParLineKind::LennardJones) => {
                    let c1 = tokens
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "Missing token class1")
                        })?
                        .parse::<usize>()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    let c2 = tokens
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "Missing token class2")
                        })?
                        .parse::<usize>()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    let e = tokens
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "Missing token epsilon")
                        })?
                        .parse::<f32>()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    let s = tokens
                        .next()
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "Missing token sigma")
                        })?
                        .parse::<f32>()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    lj_parameters
                        .entry((c1, c2))
                        .or_insert(LennardJonesItem::new(e, s));
                    lj_parameters
                        .entry((c2, c1))
                        .or_insert(LennardJonesItem::new(e, s));
                }
                None => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Unknown line kind",
                    ));
                }
            }
        }

        // Add missing Lennard-Jones parameters by using combination rules:
        // sigma_ij = (sigma_ii + sigma_jj) / 2
        // epsilon_ij = sqrt(epsilon_ii * epsilon_jj)
        for (c1, c2) in lj_parameters
            .keys()
            .cloned()
            .collect::<Vec<(usize, usize)>>()
        {
            if !lj_parameters.contains_key(&(c1, c2)) {
                let lj1 = &lj_parameters[&(c1, c1)];
                let lj2 = &lj_parameters[&(c2, c2)];
                let e = (lj1.epsilon * lj2.epsilon).sqrt();
                let s = (lj1.sigma + lj2.sigma) / 2.0;
                lj_parameters.insert((c1, c2), LennardJonesItem::new(e, s));
                lj_parameters.insert((c2, c1), LennardJonesItem::new(e, s));
            }
        }

        Ok(ParParserResult::new(lj_parameters))
    }
}
