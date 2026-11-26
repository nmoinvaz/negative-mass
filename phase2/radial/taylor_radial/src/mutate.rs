//! Deliberate errors for the controls: `mutate=NAME` on the command line turns one on. Each check of
//! taylor_polar.py has a control that must fail under its mutation.
use std::sync::OnceLock;

static M: OnceLock<String> = OnceLock::new();

pub fn set(s: &str) {
    let _ = M.set(s.to_string());
}

pub fn on(name: &str) -> bool {
    M.get().map(|s| s == name).unwrap_or(false)
}
