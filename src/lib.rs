pub mod alexander_table;
pub mod arm_type;
pub mod batch;
pub mod config;
pub mod error;
pub mod geometry;
pub mod hull;
pub mod io;
pub mod kmt;
pub mod knotsize;
pub mod knottype;
pub mod point;
pub mod polynomial;
#[cfg(feature = "python")]
pub mod python_module;

// Re-export primary public API
pub use alexander_table::AlexanderTable;
pub use arm_type::get_31knot_arm_type;
pub use batch::{
    process_frame, process_frames_parallel, process_frames_streaming, FrameResult, InputFormat,
};
pub use config::KnotConfig;
pub use error::{KnotError, Result};
pub use io::{read_data_lammps, read_data_xyz, read_data_xyz_frames};
pub use knotsize::{find_knot_core, KnotCoreResult};
pub use knottype::get_knottype;
pub use point::Point3;
