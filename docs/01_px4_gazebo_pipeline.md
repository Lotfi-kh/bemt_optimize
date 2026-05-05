# PX4 / Gazebo Simulation Pipeline

## Architecture overview

```
Gazebo (gz-sim)
  └─ x500 model  ──►  gz_bridge (PX4 module)
                            │
                    uORB topics published
                            │
              ┌─────────────┼─────────────────┐
              │             │                 │
        esc_status   vehicle_attitude   vehicle_local_position
        (RPM, V, I)   (quaternion)      (NED position/velocity)
              │             │                 │
              └─────────────┼─────────────────┘
                            ▼
                    bemt_optimize (PX4 module)
                            │
                    bemt::calculate(input, output)
                            │
                    PX4_INFO log lines (1 Hz)
```

## Topic details

| Topic | Field | Unit | Notes |
|---|---|---|---|
| `esc_status` | `esc[i].esc_rpm` | RPM | Converted from rad/s in gz_bridge (commit b282194) |
| `vehicle_attitude` | `q[0..3]` | — | ned_from_body quaternion [w, x, y, z] |
| `vehicle_local_position` | `vx, vy, vz` | m/s | NED frame |
| `vehicle_air_data` | `rho` | kg/m³ | ISA fallback if absent |
| `vehicle_air_data` | `baro_temp_celcius` | °C | Used for Sutherland viscosity |
| `wind` | `windspeed_north/east` | m/s | Set to 0 if topic absent |
| `battery_status` | `voltage_v` | V | Used for motor voltage fallback |

## Coordinate frames

### NED earth frame
North = +x, East = +y, Down = +z.

### PX4 FRD body frame
Forward = +x, Right = +y, Down = +z.  Matches NED at zero roll/pitch/yaw.

### SDF FLU frame (Gazebo model files)
Forward = +x, Left = +y, Up = +z.  Rotor z-offset in SDF: +0.06 m → −0.06 m in FRD.

### Quaternion convention
`vehicle_attitude.q` is the rotation from FRD body to NED (ned_from_body).
To rotate a NED vector to body frame: `v_body = DCM(q)^T × v_ned`.

## Running BEMT in SITL

```bash
# Terminal 1 — start PX4 SITL with Gazebo
make px4_sitl gz_x500

# Terminal 2 — once simulation is running and armed
listener esc_status          # verify ESC RPM is non-zero
bemt_optimize start
bemt_optimize status
```

Expected output (1 Hz):
```
[bemt] src=esc_status valid=4/4 raw_esc_rpm_m0=5832 rpm_m0=5832.0 ...
[bemt] vx=0.00 vy=0.00 vz=0.00 wind_n=0.00 wind_e=0.00 [m/s NED]
[bemt] kin V0=0.00 J=0.000 ...
```

## Offline post-processing

```bash
# Export a flight log to CSV + validation report
python tools/ulog_to_csv.py --ulog /path/to/log.ulg --output-dir output/

# Run tests (no .ulg or pyulog required)
python -m pytest tests/ -v
```

Baseline exporter conventions:
- The flat CSV is synchronised onto the `vehicle_local_position` timeline.
- Topic joins use `merge_asof(..., direction="nearest")` as an engineering baseline assumption.
- Current tolerances are 20 ms for fast topics and 100 ms for slow topics (`wind`, `battery_status`).
- `v_normal` and `J_n` are signed axial quantities.
- `alpha_disk` uses `abs(v_normal)` in the denominator.
- ESC RPM is the only accepted RPM source; no RPM is inferred from actuator commands.
- If RPM is missing or invalid, `J`, `J_n`, `J_p`, and `Re_07` remain `NaN`.
