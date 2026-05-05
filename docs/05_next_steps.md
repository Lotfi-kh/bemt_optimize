# Next Steps

## Immediate (SITL validation)

- [ ] Collect a representative SITL flight log covering hover, forward flight, and descent.
- [ ] Run `tools/ulog_to_csv.py` on the log and review the validation report for:
  - Quaternion norm quality (expect all samples within tolerance).
  - ESC RPM availability (expect `rpm_source = esc_status`).
  - Missing RPM handling (expect no command-based RPM fabrication; RPM-dependent outputs stay `NaN` where RPM is unavailable).
  - Signed axial convention (expect `J_n` sign to follow `v_normal`, while `alpha_disk` uses `abs(v_normal)`).
  - J range in forward flight (expect J ≈ 0.05–0.5 for typical airspeeds).
  - Re_07 range (expect ~50 000–120 000 at typical RPMs).
- [ ] Plot J vs. airspeed and Re_07 vs. RPM from the CSV to confirm expected trends.

## Model validation

- [ ] Iterate induced-velocity correction (`calculate_induced_axial_velocity`) to convergence
  and compare against single-shot values to quantify the approximation error.
- [ ] Validate c_t against momentum theory (T = 2 ρ A vi²) using SITL thrust data if available.
- [ ] Document `c_t`, `c_q`, `c_p` as provisional until a reference dataset is obtained.

## Hardware transition

- [ ] Verify ESC RPM reporting on the target real-hardware platform.
- [ ] Confirm propeller geometry constants match the installed propeller (not just the SITL model).
- [ ] Measure or obtain real-flight air-data for density and temperature validation.
- [ ] Account for rotor-wake interaction between adjacent rotors (currently neglected).

## Software improvements

- [ ] Add an iterated induced-velocity solver as an option (controlled by a compile-time flag
  or parameter), retaining the single-shot version as the default.
- [ ] Publish BEMT outputs to a custom uORB topic for logging and downstream use.
- [ ] Add a parameter (`BEMT_ENABLE`) to gate the module at runtime.
- [ ] Write integration tests that run the PX4 SITL for a fixed duration and verify
  BEMT output ranges against expected values.

## Offline tool

- [ ] Add a `--plot` flag to `tools/ulog_to_csv.py` to generate J vs. time and
  Re_07 vs. time figures using matplotlib.
- [ ] Support multi-instance ESC status logs (multiple `esc_status` instances).
- [ ] Add a `--validate` mode that exits non-zero if the validation report contains warnings.
