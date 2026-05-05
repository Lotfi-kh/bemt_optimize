# Validation Results

**Scope:** PX4/Gazebo x500 SITL only.  These results do **not** validate aerodynamic truth,
hardware equivalence, or final BEMT accuracy.  All figures are provisional.

## SITL smoke-test observations (hover)

| Variable | Observed range (hover) | Notes |
|---|---|---|
| ESC RPM (motor 0) | 5 800 – 6 200 RPM | Varies with commanded throttle |
| `v_inf` (motor 0) | 0.00 – 0.15 m/s | Near-zero in hover; small rigid-body contribution |
| `J` (motor 0) | 0.000 – 0.006 | Expected ≈ 0 in hover |
| `Re_07` (motor 0) | ~80 000 | Dominated by blade-tip speed at hover RPM |
| `alpha_disk` (motor 0) | 0.00 – 0.05 rad | Near-zero in hover |

All values computed from ESC-reported RPM (`esc_status.esc_rpm`).  No fabrication.

## Unit-test coverage

`tests/test_ulog_to_csv.py` — 74 tests, all passing.  Key checks:

| Test class | Verified |
|---|---|
| `TestSutherlandViscosity` | Sutherland formula; reference value; warm/cold monotonicity |
| `TestRotateNedToBody` | Identity passthrough; 90° yaw; vectorised shape; magnitude preservation |
| `TestEulerAngles` | Identity → zero angles; 90° roll; 90° yaw |
| `TestHover` | J = 0; Re_07 matches closed-form hover formula |
| `TestEdgewiseFlight` | alpha_disk = π/2; J = V/(nD) thesis formula; J_p = J |
| `TestDescentFlight` | v_normal = V_desc; alpha_disk = 0; J = |J_n| |
| `TestMissingRpm` (parametrised) | NaN/0/negative RPM → J, Re_07 = NaN; v_inf finite |
| `TestQuaternionInputs` | All-motor symmetry at hover; finite output for zero velocity |
| `TestTimestampMonotonicity` | Sorted/shuffled/duplicate timestamp detection |
| `TestSchema` | SI suffixes; motor column count; advance-ratio column count |
| `TestProvenance` | Provenance flags; no c_t/c_q/c_p columns |

## Known limitations

- SITL ESC RPM reflects simulated motor dynamics; not validated against real hardware.
- Single-shot induced-velocity correction underestimates at high disk loading (hover).
  See comment in `bemt_model.cpp` `calculate()`.
- Quaternion norm validation rejects inputs with |‖q‖ − 1| > 0.01.
- No wind measurement available in typical SITL logs; tool assumes zero wind when absent.
