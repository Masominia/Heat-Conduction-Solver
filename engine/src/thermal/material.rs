//! Thermal properties.

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Material {
    /// Thermal conductivity, W / (m K).
    pub conductivity: f64,
    /// Volumetric heat capacity, J / (m^3 K).
    pub capacity: f64,
}

pub const SILICON: Material = Material { conductivity: 148.0, capacity: 1.63e6 };
pub const SILICA: Material = Material { conductivity: 1.38, capacity: 1.61e6 };
pub const NITRIDE: Material = Material { conductivity: 20.0, capacity: 2.2e6 };
pub const TITANIUM_NITRIDE: Material = Material { conductivity: 28.0, capacity: 3.2e6 };
pub const ALUMINIUM: Material = Material { conductivity: 237.0, capacity: 2.42e6 };
pub const AIR: Material = Material { conductivity: 0.026, capacity: 1.2e3 };

/// si, sio2, sin, tin, al, air (and their long names).
pub fn by_name(text: &str) -> Option<Material> {
    match text {
        "si" | "silicon" => Some(SILICON),
        "sio2" | "silica" | "oxide" => Some(SILICA),
        "sin" | "nitride" => Some(NITRIDE),
        "tin" => Some(TITANIUM_NITRIDE),
        "al" | "aluminium" => Some(ALUMINIUM),
        "air" => Some(AIR),
        _ => None,
    }
}

impl Material {
    /// Thermal diffusivity in um^2 per microsecond (numerically equal to m^2/s times 1e6).
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / self.capacity * 1e6
    }
}
