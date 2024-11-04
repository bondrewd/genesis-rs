use std::time::{Duration, Instant};

#[derive(Debug)]
pub enum TimerState {
    Running,
    Stopped,
}

#[derive(Debug)]
pub struct Timer {
    time: Duration,
    instant: Instant,
    state: TimerState,
}

impl Default for Timer {
    fn default() -> Self {
        Self {
            time: Duration::from_secs(0),
            instant: Instant::now(),
            state: TimerState::Stopped,
        }
    }
}

impl Timer {
    pub fn new(time: Duration, instant: Instant, state: TimerState) -> Self {
        Timer {
            time,
            instant,
            state,
        }
    }

    pub fn start(&mut self) {
        match self.state {
            TimerState::Running => (),
            TimerState::Stopped => {
                self.instant = Instant::now();
                self.state = TimerState::Running;
            }
        }
    }

    pub fn stop(&mut self) {
        match self.state {
            TimerState::Running => {
                self.time += self.instant.elapsed();
                self.state = TimerState::Stopped;
            }
            TimerState::Stopped => (),
        }
    }

    pub fn reset(&mut self) {
        self.time = Duration::from_secs(0);
        self.state = TimerState::Stopped;
    }

    pub fn elapsed(&self) -> Duration {
        match self.state {
            TimerState::Running => self.time + self.instant.elapsed(),
            TimerState::Stopped => self.time,
        }
    }
}
