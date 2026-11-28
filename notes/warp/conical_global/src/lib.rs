//! conical_global: the model of notes/warp/conical_global.md on multishoot, multisolver and multidouble.
pub mod connect;
pub mod hp;
pub mod json;
#[cfg(feature = "locate_lib")]
pub mod loc;
pub mod model;
pub mod offgrid;
pub mod planted;
pub mod poly;
pub mod scan;
pub mod tables;
pub mod tables_hp;
pub mod walk;
