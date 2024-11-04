use crate::system::System;
use byteorder::{LittleEndian, WriteBytesExt};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(Debug)]
pub struct DCDReporter {
    file: File,
    header_flag: bool,
}

impl DCDReporter {
    pub fn new(file: File, header_flag: bool) -> Self {
        DCDReporter { file, header_flag }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(DCDReporter::new(file, false))
    }

    pub fn write_header(&mut self, num_atoms: u32, num_frames: u32) -> io::Result<()> {
        if !self.header_flag {
            // Init buffer
            let mut buffer = Vec::new();

            // Block 1
            // Block size start
            buffer.write_u32::<LittleEndian>(84)?;
            // 01 - 04: "CORD" magic number
            buffer.write_all(b"CORD")?;
            // 05 - 08: Number of frames
            buffer.write_u32::<LittleEndian>(num_frames)?;
            // 09 - 12: Unused (first step)
            buffer.write_u32::<LittleEndian>(0)?;
            // 13 - 16: Unused (output period)
            buffer.write_u32::<LittleEndian>(0)?;
            // 17 - 20: Unused (number of time steps)
            buffer.write_u32::<LittleEndian>(0)?;
            // 21 - 40: Unused
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            // 41 - 44: Unused (time step)
            buffer.write_f32::<LittleEndian>(0.0)?;
            // 45 - 48: Unit cell flag
            buffer.write_u32::<LittleEndian>(1)?;
            // 49 - 80: Unused
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            buffer.write_u32::<LittleEndian>(0)?;
            // 81 - 84: Version
            buffer.write_u32::<LittleEndian>(24)?;
            // Block size end
            buffer.write_u32::<LittleEndian>(84)?;

            // Block 2
            // Block size start
            buffer.write_u32::<LittleEndian>(84)?;
            // 01 - 04: Number of lines
            buffer.write_u32::<LittleEndian>(1)?;
            // 05 - 84: Comment
            buffer.write_all(b"                ")?;
            buffer.write_all(b"                ")?;
            buffer.write_all(b"                ")?;
            buffer.write_all(b"                ")?;
            buffer.write_all(b"                ")?;
            // Block size end
            buffer.write_u32::<LittleEndian>(84)?;

            // Block 3
            // Block size start
            buffer.write_u32::<LittleEndian>(4)?;
            // 01 - 04: Number of particles
            buffer.write_u32::<LittleEndian>(num_atoms)?;
            // Block size end
            buffer.write_u32::<LittleEndian>(4)?;

            self.file.write_all(&buffer)?;
            self.header_flag = true;
        }
        Ok(())
    }

    pub fn write_report(&mut self, system: &System) -> io::Result<()> {
        // Check if the header has been written
        if !self.header_flag {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "DCD header must be written before writing frames",
            ));
        }

        // Init buffer
        let mut buffer = Vec::new();

        // Block size start
        buffer.write_u32::<LittleEndian>(48)?;
        // Write X coordinates
        buffer.write_f64::<LittleEndian>(system.b.x as f64)?;
        buffer.write_f64::<LittleEndian>(0.0)?;
        buffer.write_f64::<LittleEndian>(system.b.y as f64)?;
        buffer.write_f64::<LittleEndian>(0.0)?;
        buffer.write_f64::<LittleEndian>(0.0)?;
        buffer.write_f64::<LittleEndian>(system.b.z as f64)?;
        // Block size end
        buffer.write_u32::<LittleEndian>(48)?;

        // For each frame, DCD stores X, Y, and Z coordinates in separate chunks
        let block_size = system.n as u32 * 4;

        // Write X coordinates
        // Block size start
        buffer.write_u32::<LittleEndian>(block_size)?;
        // Write X coordinates
        for r in system.r.iter() {
            buffer.write_f32::<LittleEndian>(r[0])?;
        }
        // Block size end
        buffer.write_u32::<LittleEndian>(block_size)?;

        // Write X coordinates
        // Block size start
        buffer.write_u32::<LittleEndian>(block_size)?;
        // Write X coordinates
        for r in system.r.iter() {
            buffer.write_f32::<LittleEndian>(r[1])?;
        }
        // Block size end
        buffer.write_u32::<LittleEndian>(block_size)?;

        // Write X coordinates
        // Block size start
        buffer.write_u32::<LittleEndian>(block_size)?;
        // Write X coordinates
        for r in system.r.iter() {
            buffer.write_f32::<LittleEndian>(r[2])?;
        }
        // Block size end
        buffer.write_u32::<LittleEndian>(block_size)?;

        // Flush buffer
        self.file.write_all(&buffer)?;

        Ok(())
    }
}
