# Assumptions and Limitations

## Physical modelling assumptions

### Propeller geometry
- Blade section table (39 stations, chord and twist) matches the Holybro x500 v2 propeller.
  If a different propeller is installed, constants in `bemt_model.hpp` must be updated.
- Hub radius: 0.02032 m (used only as integration lower bound; no hub correction applied).

### Air properties
- Dynamic viscosity: Sutherland's law, air constants, from `vehicle_air_data.baro_temp_celcius`.
  Engineering approximation; accuracy is sufficient for Reynolds number estimation but not for
  precision viscous flow calculations.
- Air density: `vehicle_air_data.rho` preferred; ISA sea-level 1.225 kg/m³ fallback.

### Wind
- Horizontal wind from `wind` topic when available.
- **Vertical wind assumed zero at all times.** PX4 does not log vertical wind.
  This introduces inflow error in conditions with significant vertical gusts.

### Induced velocity
- Single-shot momentum-theory correction: one pass only.  No fixed-point iteration.
- May underpredict converged induced velocity at high disk loading (e.g., hover).
- Prandtl tip-loss factor is applied but the induced velocity is not iterated to convergence.

### RPM source
- `esc_status.esc_rpm` is the only accepted RPM source.
- No RPM estimation from motor commands.  If `esc_status` is absent, J and Re_07 are NaN/zero.
- ESC RPM in SITL reflects simulated motor dynamics.  Real-hardware RPM may differ.

### Rotor positions
- Four rotors at arm positions derived from `ARM_LENGTH_M = 0.175 m` and 45° arm angles.
- z-offset: −0.06 m in FRD body frame (rotors are 0.06 m above base_link in the SDF).
- No flexing, vibration, or mounting tolerance accounted for.

## Software assumptions

### Quaternion validity
- Input rejected if any component is non-finite or |‖q‖ − 1| > 0.01.
- Identity fallback in `bemt_optimizer.cpp` if any component is non-finite.

### Coordinate frames
- `vehicle_attitude.q` must be ned_from_body (standard PX4 convention).
  If the quaternion convention changes, `rotate_ned_to_body()` in `bemt_model.cpp`
  and `tools/ulog_to_csv.py` must both be updated.

### Timestamp synchronisation (offline tool)
- Topics merged by nearest timestamp (`merge_asof`) with a 20 ms tolerance for fast topics
  and 100 ms for slow topics (wind, battery).  Gaps larger than the tolerance produce NaN.

## Out of scope

- Hardware-in-the-loop (HIL) or real-flight validation.
- Multi-rotor aerodynamic interaction (rotor wake interference).
- Blade flapping or lead-lag dynamics.
- Aerodynamic validation against wind-tunnel data.
- Thrust or torque coefficient (c_t, c_q) validation — these are printed as provisional
  and must not be cited as validated results.
