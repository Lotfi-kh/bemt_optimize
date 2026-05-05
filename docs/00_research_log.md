# Research Log

## 2025-05 — BEMT module integration (PX4 SITL, x500 Gazebo)

### Scope
Develop and validate a Blade Element Momentum Theory (BEMT) aerodynamic state-extraction module
for the PX4 x500 quadrotor SITL environment.  The module runs onboard as a PX4 task, subscribes
to uORB topics, and prints per-rotor inflow and nondimensional coefficient summaries for offline
analysis.

### Physical model
- Propeller: 10-inch (0.254 m diameter), 2-bladed.  Geometry from Holybro x500 v2.
- Blade section table: 39 radial stations, chord interpolated at 0.7R for Reynolds number.
- Dynamic viscosity: Sutherland's law derived from `vehicle_air_data.baro_temp_celcius`.
- Density: `vehicle_air_data.rho`; ISA sea-level fallback (1.225 kg/m³) when absent.
- Induced velocity: single-shot momentum-theory correction (not iterated).

### Frame convention
- PX4 body frame: FRD (x-forward, y-right, z-down).
- `vehicle_attitude.q` = ned_from_body quaternion.  `v_body = DCM(q)^T × v_ned`.
- Rotor z-offset: −0.06 m in FRD (z-down).  Source: x500_base Gazebo SDF, rotor links at
  z = +0.06 m in SDF FLU (z-up).

### RPM source
- ESC-reported RPM (`esc_status.esc[i].esc_rpm`) is the sole RPM source.
- No RPM fabrication from command signals in the offline exporter.
- If `esc_status` is absent, or if RPM is missing / non-positive, the exporter leaves
  RPM-dependent outputs (`J`, `J_n`, `J_p`, `Re_07`) as `NaN`.

### Axial convention
- `v_normal` is a signed axial inflow quantity in the rotor-axis convention used by the baseline model.
- `J_n` is also signed because it is defined directly from `v_normal / (nD)`.
- `alpha_disk = atan2(v_inplane, abs(v_normal))`, so the denominator uses axial magnitude.

### Key findings (SITL baseline)
- ESC RPM available in SITL only after `bemt_optimize start` is run post-arming.
- Gazebo esc_rpm topic reports revolutions per minute (already converted from rad/s
  in `gz_bridge` after commit `b282194`).
- Quaternion copied atomically: if any component is non-finite, identity fallback is used.
- Stack size set to 4000 bytes to accommodate BEMT computation.

### Offline tool
`tools/ulog_to_csv.py` reads a PX4 `.ulg` log and produces:
- `output/x500_baseline_timeseries.csv` — flat, single-rate CSV.
- `output/x500_baseline_validation_report.md` — data-quality report.
- The flat CSV is synchronised onto the `vehicle_local_position` timeline using
  nearest-timestamp `merge_asof` joins.
- Current tolerances: 20 ms for fast topics and 100 ms for slow topics (`wind`, `battery_status`).

Tests: `tests/test_ulog_to_csv.py` covers numerical helpers plus a higher-level
`build_timeseries()` merge/provenance check, with no real `.ulg` or `pyulog` required.
