/// Unified configuration for all knot analysis operations.
///
/// Consolidates all hyperparameters and mode flags that were previously
/// scattered across `KnotOptions`, `KnotSizeOptions`, and hardcoded constants.
#[derive(Clone, Debug)]
pub struct KnotConfig {
    // ── Mode flags ──
    /// Open chain (`false`) or ring/closed chain (`true`).
    pub is_ring: bool,
    /// Apply KMT simplification before computing polynomial (faster, same result).
    pub faster: bool,
    /// Print debug information to stderr.
    pub debug: bool,

    // ── Hull hyperparameters ──
    /// Threshold for considering a point "on" a hull face.
    /// Points closer than this to a face use center-away extension instead of
    /// normal-based projection. (C++ `kHullPlaneEpsilon`)
    pub hull_plane_epsilon: f64,
    /// Scale factor for endpoint extension through the convex hull.
    /// Larger values push new endpoints farther from the chain. (C++ `kExtendFactor`)
    pub extend_factor: f64,

    // ── Knot core search hyperparameters ──
    /// Number of rotation offsets to try when searching for knot core in ring mode.
    /// Each rotation shifts by `n / num_rotations` points. (C++ hardcoded as `4`)
    pub num_rotations: u32,
}

impl Default for KnotConfig {
    fn default() -> Self {
        KnotConfig {
            is_ring: false,
            faster: false,
            debug: false,
            hull_plane_epsilon: 5e-3,
            extend_factor: 100.0,
            num_rotations: 4,
        }
    }
}
