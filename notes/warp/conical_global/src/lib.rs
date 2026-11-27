//! conical_global: the model of notes/warp/conical_global.md on multishoot, multisolver and multidouble.
#![cfg_attr(feature = "implicit_lib", allow(dead_code, unused_imports))]
pub mod connect;
pub mod hp;
#[cfg(feature = "implicit_lib")]
pub mod hp_lib;
pub mod json;
pub mod model;
pub mod offgrid;
pub mod planted;
pub mod poly;
#[cfg(feature = "implicit_lib")]
pub mod poly_lib;
pub mod scan;
pub mod tables;
pub mod tables_hp;
pub mod walk;
