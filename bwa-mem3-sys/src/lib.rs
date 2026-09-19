#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(clippy::all)]

#[cfg(feature = "regenerate-bindings")]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
#[cfg(not(feature = "regenerate-bindings"))]
include!("bindings.rs");

pub mod build_info;

#[cfg(test)]
#[path = "../build_support/compiler_floor.rs"]
mod compiler_floor;
