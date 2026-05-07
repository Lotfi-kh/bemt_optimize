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

All RPM-dependent values are computed from ESC-reported RPM (`esc_status.esc_rpm`) when available.
No RPM is fabricated from actuator commands. If RPM is missing or invalid, the exporter leaves
`J`, `J_n`, `J_p`, and `Re_07` as `NaN`.

`v_normal` and `J_n` use a signed axial convention. `alpha_disk` is computed as
`atan2(v_inplane, abs(v_normal))`, so the denominator uses axial magnitude rather than signed axial velocity.

## Unit-test coverage

`tests/test_ulog_to_csv.py` covers both numerical helpers and higher-level exporter assembly. Key checks:

| Test class | Verified |
|---|---|
| `TestSutherlandViscosity` | Sutherland formula; reference value; warm/cold monotonicity |
| `TestRotateNedToBody` | Identity passthrough; 90° yaw; vectorised shape; magnitude preservation |
| `TestEulerAngles` | Identity → zero angles; 90° roll; 90° yaw |
| `TestHover` | J = 0; Re_07 matches closed-form hover formula |
| `TestEdgewiseFlight` | alpha_disk = π/2; J = V/(nD) thesis formula; J_p = J |
| `TestDescentFlight` | v_normal = V_desc; alpha_disk = 0; J = |J_n| |
| `TestClimbFlight` | Signed axial convention: climb gives negative `v_normal` and negative `J_n` |
| `TestMissingRpm` (parametrised) | NaN/0/negative RPM → `J`, `J_n`, `J_p`, `Re_07` = `NaN`; `v_inf` and `alpha_disk` remain finite |
| `TestQuaternionInputs` | All-motor symmetry at hover; finite output for zero velocity |
| `TestBuildTimeseries` | `vehicle_local_position` timeline sync, merged output columns, provenance flags, and missing-topic defaults together |
| `TestBuildTimeseriesWithBemtTopic` | Preferred path: `bemt_rotor_state` present; rotor_state_source correct; kinematic/corrected columns present; legacy columns populated from kinematic signed axial |
| `TestBuildTimeseriesFallbackReconstruction` | Fallback path: topic absent; rotor_state_source = reconstruction; corrected columns absent |
| `TestBemtTopicPreference` | v_inf and re_07 sourced from BEMT topic, not offline reconstruction |
| `TestMergeBemtHelpers` | Unit tests for _merge_bemt_legacy and _merge_bemt_new_cols |
| `TestTimestampMonotonicity` | Sorted/shuffled/duplicate timestamp detection |
| `TestSchema` | SI suffixes; motor column count; advance-ratio column count; kinematic/corrected schema |
| `TestProvenance` | Provenance flags; no c_t/c_q/c_p columns; rotor_state_source; has_bemt_topic |

## Known limitations

- SITL ESC RPM reflects simulated motor dynamics; not validated against real hardware.
- Single-shot induced-velocity correction underestimates at high disk loading (hover).
  See comment in `bemt_model.cpp` `calculate()`.
- Quaternion norm validation rejects inputs with |‖q‖ − 1| > 0.01.
- No wind measurement available in typical SITL logs; tool assumes zero wind when absent.
- Flat CSV synchronisation uses nearest-timestamp joins onto the `vehicle_local_position` timeline;
  this is a baseline engineering assumption, not a validated estimator reconstruction.
- `bemt_rotor_state` is absent from logs collected before this topic was introduced; the
  Python exporter falls back to reconstruction for those logs automatically.
- Corrected columns (`corrected_*`) are only available when the `bemt_rotor_state` topic
  is present; they are absent from fallback-path CSVs.
