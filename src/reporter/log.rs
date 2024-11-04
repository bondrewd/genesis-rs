use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(Debug)]
pub struct LOGReporter {
    file: File,
}

impl LOGReporter {
    pub fn new(file: File) -> Self {
        LOGReporter { file }
    }

    pub fn with_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(LOGReporter::new(file))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn write_report(
        &mut self,
        setup_time: std::time::Duration,
        output_time: std::time::Duration,
        dynamics_time: std::time::Duration,
        total_time: std::time::Duration,
    ) -> io::Result<()> {
        // Init buffer
        let mut buffer = Vec::new();

        // Convert durations to seconds for easier formatting
        let setup_secs = setup_time.as_secs_f64();
        let output_secs = output_time.as_secs_f64();
        let dynamics_secs = dynamics_time.as_secs_f64();
        let total_secs = total_time.as_secs_f64();

        // Calculate percentages
        let setup_percent = (setup_secs / total_secs) * 100.0;
        let output_percent = (output_secs / total_secs) * 100.0;
        let dynamics_percent = (dynamics_secs / total_secs) * 100.0;

        // Print the profile report
        buffer.extend_from_slice("Profiler:\n".as_bytes());
        buffer.extend_from_slice(
            format!("setup    = {:>5.1}% {:.3}s\n", setup_percent, setup_secs).as_bytes(),
        );
        buffer.extend_from_slice(
            format!("output   = {:>5.1}% {:.3}s\n", output_percent, output_secs).as_bytes(),
        );
        buffer.extend_from_slice(
            format!(
                "dynamics = {:>5.1}% {:.3}s\n",
                dynamics_percent, dynamics_secs
            )
            .as_bytes(),
        );
        buffer.extend_from_slice(format!("total    = 100.0% {:.3}s\n", total_secs).as_bytes());

        // Flush buffer
        self.file.write_all(&buffer)?;

        Ok(())
    }
}
