#!/usr/bin/env python3
"""
tools/ulog_to_csv.py  –  PX4 SITL ULog → flat CSV + markdown validation report.

Engineering classification: offline post-processing tool.
Scope: baseline x500 Gazebo SITL only.
No aerodynamic validation claims are made.

Usage:
    python tools/ulog_to_csv.py --ulog /path/to/log.ulg [--output-dir output/]

Outputs:
    <output-dir>/x500_baseline_timeseries.csv
    <output-dir>/x500_baseline_validation_report.md
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from pyulog import ULog as _ULog
except ImportError:  # allow tests to import without pyulog installed
    _ULog = None  # type: ignore[assignment,misc]

# ---------------------------------------------------------------------------
# Physical constants — ported from bemt_model.hpp (src/modules/bemt_optimizer).
# These are configuration values for the default x500 Gazebo SITL model.
# ---------------------------------------------------------------------------

MOTOR_COUNT: int = 4
BLADE_COUNT: int = 2
PROP_DIAMETER_M: float = 0.254       # 10-inch propeller
PROP_RADIUS_M: float = 0.127
HUB_RADIUS_M: float = 0.02032
ARM_LENGTH_M: float = 0.175

# Rotor z-offset in PX4 FRD body frame (z-down).
# Source: x500_base Gazebo SDF; all four rotor links are at z = +0.06 m in SDF FLU (z-up),
# mapping to −0.06 m in PX4 FRD (z-down). See bemt_model.hpp kRotorPosZBodyM for derivation.
ROTOR_POS_Z_BODY_M: float = -0.06

# Rotor positions in PX4 FRD body frame [x, y, z] (metres).
ROTOR_POSITIONS: np.ndarray = np.array(
    [
        [ ARM_LENGTH_M * 0.7071,  ARM_LENGTH_M * 0.7071, ROTOR_POS_Z_BODY_M],  # motor 0
        [-ARM_LENGTH_M * 0.7071, -ARM_LENGTH_M * 0.7071, ROTOR_POS_Z_BODY_M],  # motor 1
        [ ARM_LENGTH_M * 0.7071, -ARM_LENGTH_M * 0.7071, ROTOR_POS_Z_BODY_M],  # motor 2
        [-ARM_LENGTH_M * 0.7071,  ARM_LENGTH_M * 0.7071, ROTOR_POS_Z_BODY_M],  # motor 3
    ],
    dtype=float,
)

# Rotor spin axis in PX4 FRD body frame. Rotors thrust upward (−z direction in FRD).
ROTOR_AXIS: np.ndarray = np.array([0.0, 0.0, -1.0])

# Blade-section table — 39 sections, ported from bemt_model.hpp.
SECTION_RADIUS_M: np.ndarray = np.array([
    0.02032, 0.02184, 0.02337, 0.02489, 0.02642, 0.02794, 0.02946, 0.03112,
    0.03401, 0.03703, 0.04006, 0.04308, 0.04610, 0.04912, 0.05215, 0.05517,
    0.05819, 0.06121, 0.06424, 0.06726, 0.07028, 0.07330, 0.07633, 0.07935,
    0.08237, 0.08539, 0.08842, 0.09144, 0.09446, 0.09749, 0.10051, 0.10353,
    0.10655, 0.10958, 0.11260, 0.11562, 0.11864, 0.12165, 0.12448,
], dtype=float)

SECTION_CHORD_M: np.ndarray = np.array([
    0.01867, 0.01978, 0.02081, 0.02178, 0.02268, 0.02352, 0.02428, 0.02504,
    0.02617, 0.02710, 0.02777, 0.02818, 0.02834, 0.02825, 0.02801, 0.02771,
    0.02736, 0.02694, 0.02648, 0.02596, 0.02540, 0.02478, 0.02412, 0.02341,
    0.02266, 0.02188, 0.02105, 0.02018, 0.01928, 0.01835, 0.01739, 0.01639,
    0.01537, 0.01432, 0.01325, 0.01216, 0.01104, 0.00973, 0.00713,
], dtype=float)

# Chord at 0.7R interpolated from the blade-section table — not hardcoded.
_RADIUS_07_M: float = 0.7 * PROP_RADIUS_M  # 0.0889 m
CHORD_07_M: float = float(np.interp(_RADIUS_07_M, SECTION_RADIUS_M, SECTION_CHORD_M))

# Sutherland's law constants for air.
_SUTH_T_REF_K: float = 273.15
_SUTH_MU_REF: float = 1.716e-5   # Pa·s
_SUTH_C: float = 110.4            # K

# ISA sea-level fallbacks (used only when vehicle_air_data is absent).
_ISA_RHO: float = 1.225    # kg/m³
_ISA_TEMP_C: float = 15.0  # °C

# Nearest-timestamp merge tolerance.
_SYNC_TOL_US: int = 20_000   # 20 ms, for fast topics
_SLOW_TOL_US: int = 100_000  # 100 ms, for wind/battery


# ---------------------------------------------------------------------------
# Pure numerical helpers — no I/O; importable for unit testing without pyulog.
# ---------------------------------------------------------------------------

def sutherland_viscosity(temp_c: float) -> float:
    """Dynamic viscosity [Pa·s] from temperature [°C] via Sutherland's law for air."""
    T_k = temp_c + 273.15
    if T_k <= 0.0:
        return _SUTH_MU_REF
    ratio = T_k / _SUTH_T_REF_K
    return _SUTH_MU_REF * (ratio ** 1.5) * (_SUTH_T_REF_K + _SUTH_C) / (T_k + _SUTH_C)


def rotate_ned_to_body(q_wxyz: np.ndarray, v_ned: np.ndarray) -> np.ndarray:
    """Rotate vector(s) from NED to PX4 FRD body frame.

    Per PX4 convention, ``vehicle_attitude.q`` represents the rotation from the
    FRD body frame to the NED earth frame (ned_from_body). Therefore:

        v_body = DCM(q)^T * v_ned

    This matches ``rotate_ned_to_body()`` in ``bemt_model.cpp``.

    Args:
        q_wxyz: quaternion(s) [w, x, y, z], shape (4,) or (N, 4).
        v_ned:  vector(s) [x, y, z],         shape (3,) or (N, 3).

    Returns:
        v_body of the same leading shape as v_ned.
    """
    qw = q_wxyz[..., 0]; qx = q_wxyz[..., 1]
    qy = q_wxyz[..., 2]; qz = q_wxyz[..., 3]

    r11 = 1.0 - 2.0 * (qy * qy + qz * qz)
    r12 = 2.0 * (qx * qy - qw * qz)
    r13 = 2.0 * (qx * qz + qw * qy)
    r21 = 2.0 * (qx * qy + qw * qz)
    r22 = 1.0 - 2.0 * (qx * qx + qz * qz)
    r23 = 2.0 * (qy * qz - qw * qx)
    r31 = 2.0 * (qx * qz - qw * qy)
    r32 = 2.0 * (qy * qz + qw * qx)
    r33 = 1.0 - 2.0 * (qx * qx + qy * qy)

    # DCM^T * v_ned  →  v_body
    vx = r11 * v_ned[..., 0] + r21 * v_ned[..., 1] + r31 * v_ned[..., 2]
    vy = r12 * v_ned[..., 0] + r22 * v_ned[..., 1] + r32 * v_ned[..., 2]
    vz = r13 * v_ned[..., 0] + r23 * v_ned[..., 1] + r33 * v_ned[..., 2]
    return np.stack([vx, vy, vz], axis=-1)


def quat_to_euler_rad(
    q_wxyz: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert [w, x, y, z] quaternion(s) to ZYX Euler angles (roll, pitch, yaw) in radians."""
    qw = q_wxyz[..., 0]; qx = q_wxyz[..., 1]
    qy = q_wxyz[..., 2]; qz = q_wxyz[..., 3]
    roll = np.arctan2(2.0 * (qw * qx + qy * qz), 1.0 - 2.0 * (qx * qx + qy * qy))
    pitch = np.arcsin(np.clip(2.0 * (qw * qy - qz * qx), -1.0, 1.0))
    yaw = np.arctan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy * qy + qz * qz))
    return roll, pitch, yaw


def compute_rotor_inflow_state(
    motor_idx: int,
    air_vel_ned: np.ndarray,      # (N, 3)  NED air velocity = vehicle_vel - wind
    omega_body_rad_s: np.ndarray, # (N, 3)  body angular rates [roll, pitch, yaw]
    q_wxyz: np.ndarray,           # (N, 4)  attitude quaternion [w, x, y, z]
    rpm: np.ndarray,              # (N,)    motor RPM; NaN where unavailable
    rho_kg_m3: np.ndarray,        # (N,)    air density
    mu_pa_s: np.ndarray,          # (N,)    dynamic viscosity
) -> Dict[str, np.ndarray]:
    """Compute nondimensional rotor state for one motor (vectorised over N samples).

    Thesis baseline equations:
        J        = V_inf  / (n D)
        J_n      = V_normal / (n D)           (signed: inherits sign from axial component)
        J_p      = V_inplane / (n D)
        Re_07    = rho * V_rel_07 * c_07 / mu
        V_rel_07 = sqrt((omega * 0.7R)^2 + V_inf^2)
        alpha_disk = atan2(V_inplane, |V_normal|)

    If rpm is NaN or <= 0, J, J_n, J_p, and Re_07 are set to NaN.
    V_inf, V_normal, V_inplane, and alpha_disk are always computed (no RPM required).

    Returns:
        Dict with keys: v_inf, v_normal, v_inplane, alpha_disk, j, j_n, j_p, re_07.
        All arrays have shape (N,).
    """
    r_pos = ROTOR_POSITIONS[motor_idx]  # (3,)

    # Rotate air velocity from NED to body frame.
    air_vel_body = rotate_ned_to_body(q_wxyz, air_vel_ned)  # (N, 3)

    # Add rigid-body velocity at rotor position: omega × r_pos.
    rigid_body_vel = np.cross(omega_body_rad_s, r_pos)      # (N, 3)
    local_air_vel_body = air_vel_body + rigid_body_vel       # (N, 3)

    # Inflow = −local_air_vel_body (air flowing into rotor disk).
    inflow = -local_air_vel_body  # (N, 3)

    # Decompose into axial and in-plane components.
    # dot(inflow, [0, 0, -1]) = -inflow_z
    signed_axial = np.dot(inflow, ROTOR_AXIS)               # (N,)
    axial_vec = signed_axial[:, None] * ROTOR_AXIS[None, :] # (N, 3)
    inplane_vec = inflow - axial_vec                        # (N, 3)
    v_inplane = np.linalg.norm(inplane_vec, axis=1)         # (N,)
    v_normal = signed_axial                                 # signed; NaN for alpha uses |.|
    v_inf = np.sqrt(signed_axial ** 2 + v_inplane ** 2)
    alpha_disk = np.arctan2(v_inplane, np.abs(signed_axial))

    out: Dict[str, np.ndarray] = {
        "v_inf": v_inf,
        "v_normal": v_normal,
        "v_inplane": v_inplane,
        "alpha_disk": alpha_disk,
    }

    # RPM-dependent quantities: NaN where RPM is unavailable or zero.
    valid_rpm = np.isfinite(rpm) & (rpm > 0.0)
    n_rev_s = np.where(valid_rpm, rpm / 60.0, np.nan)
    denom = n_rev_s * PROP_DIAMETER_M

    out["j"]   = np.where(valid_rpm, v_inf    / denom, np.nan)
    out["j_n"] = np.where(valid_rpm, v_normal / denom, np.nan)
    out["j_p"] = np.where(valid_rpm, v_inplane / denom, np.nan)

    omega_rad_s = np.where(valid_rpm, rpm * (2.0 * math.pi / 60.0), np.nan)
    v_tan_07 = omega_rad_s * _RADIUS_07_M
    v_rel_07 = np.sqrt(v_tan_07 ** 2 + v_inf ** 2)
    out["re_07"] = np.where(valid_rpm, rho_kg_m3 * v_rel_07 * CHORD_07_M / mu_pa_s, np.nan)

    return out


# ---------------------------------------------------------------------------
# ULog loading and topic extraction.
# ---------------------------------------------------------------------------

def _load_ulog(path: Path) -> "_ULog":
    if _ULog is None:
        raise ImportError(
            "pyulog is required to load .ulg files. "
            "Install it with: pip install pyulog"
        )
    return _ULog(str(path))


def _topic_df(ulog: "_ULog", name: str, instance: int = 0) -> Optional[pd.DataFrame]:
    """Return a DataFrame for one topic, or None if the topic is absent."""
    try:
        dataset = ulog.get_dataset(name, instance)
    except Exception:
        return None
    if dataset is None:
        return None
    df = pd.DataFrame(dataset.data)
    if "timestamp" not in df.columns:
        return None
    return df.sort_values("timestamp").reset_index(drop=True)


def _asof(
    base: pd.DataFrame,
    other: pd.DataFrame,
    cols: list,
    renames: Optional[dict] = None,
    tolerance_us: int = _SYNC_TOL_US,
) -> pd.DataFrame:
    """Left-join selected columns from other onto base by nearest timestamp."""
    keep = ["timestamp"] + [c for c in cols if c in other.columns]
    right = other[keep].copy()
    if renames:
        right = right.rename(columns=renames)
    return pd.merge_asof(
        base.sort_values("timestamp"),
        right.sort_values("timestamp"),
        on="timestamp",
        direction="nearest",
        tolerance=tolerance_us,
    )


# ---------------------------------------------------------------------------
# Main pipeline.
# ---------------------------------------------------------------------------

def build_timeseries(ulog_path: Path) -> Tuple[pd.DataFrame, dict]:
    """Load a ULog and return (timeseries_df, metadata_dict).

    All required topics are extracted and synchronised by nearest timestamp
    onto the vehicle_local_position timeline.  Missing optional topics are
    documented in metadata and replaced by safe defaults (NaN or zero) rather
    than fabricated values.
    """
    ulog = _load_ulog(ulog_path)

    meta: dict = {
        "input_file": str(ulog_path),
        "topics_found": [],
        "topics_missing": [],
    }

    def _try(name: str, instance: int = 0) -> Optional[pd.DataFrame]:
        df = _topic_df(ulog, name, instance)
        (meta["topics_found"] if df is not None else meta["topics_missing"]).append(name)
        return df

    lpos    = _try("vehicle_local_position")
    att     = _try("vehicle_attitude")
    ang_vel = _try("vehicle_angular_velocity")
    air_dat = _try("vehicle_air_data")
    wind_df = _try("wind")
    batt_df = _try("battery_status")
    motors  = _try("actuator_motors")
    esc     = _try("esc_status")

    if lpos is None or att is None:
        raise RuntimeError(
            "vehicle_local_position and vehicle_attitude are required but missing from the log."
        )

    # --- Base: vehicle_local_position ---
    base = lpos[["timestamp", "x", "y", "z", "vx", "vy", "vz"]].copy()
    base = base.rename(columns={
        "x": "position_n_m", "y": "position_e_m", "z": "position_d_m",
        "vx": "velocity_n_m_s", "vy": "velocity_e_m_s", "vz": "velocity_d_m_s",
    })

    # --- Attitude ---
    base = _asof(base, att,
                 ["q[0]", "q[1]", "q[2]", "q[3]"],
                 {"q[0]": "q_w", "q[1]": "q_x", "q[2]": "q_y", "q[3]": "q_z"})

    # --- Angular velocity ---
    if ang_vel is not None:
        base = _asof(base, ang_vel,
                     ["xyz[0]", "xyz[1]", "xyz[2]"],
                     {"xyz[0]": "roll_rate_rad_s",
                      "xyz[1]": "pitch_rate_rad_s",
                      "xyz[2]": "yaw_rate_rad_s"})
    else:
        base["roll_rate_rad_s"] = 0.0
        base["pitch_rate_rad_s"] = 0.0
        base["yaw_rate_rad_s"] = 0.0

    # --- Air data ---
    has_air_data = air_dat is not None
    if has_air_data:
        base = _asof(base, air_dat,
                     ["rho", "baro_temp_celcius"],
                     {"rho": "air_density_kg_m3", "baro_temp_celcius": "air_temperature_c"})
    else:
        base["air_density_kg_m3"] = _ISA_RHO
        base["air_temperature_c"] = _ISA_TEMP_C

    # --- Wind ---
    has_wind = wind_df is not None
    if has_wind:
        base = _asof(base, wind_df,
                     ["windspeed_north", "windspeed_east"],
                     {"windspeed_north": "wind_n_m_s", "windspeed_east": "wind_e_m_s"},
                     tolerance_us=_SLOW_TOL_US)
    else:
        base["wind_n_m_s"] = 0.0
        base["wind_e_m_s"] = 0.0
    # Vertical wind is not logged by PX4; assumed zero.
    base["wind_d_m_s"] = 0.0

    # --- Battery ---
    has_battery = batt_df is not None
    if has_battery:
        base = _asof(base, batt_df,
                     ["voltage_v", "current_a", "remaining"],
                     {"voltage_v": "battery_voltage_v",
                      "current_a": "battery_current_a",
                      "remaining": "battery_remaining_frac"},
                     tolerance_us=_SLOW_TOL_US)
    else:
        base["battery_voltage_v"] = float("nan")
        base["battery_current_a"] = float("nan")
        base["battery_remaining_frac"] = float("nan")

    # --- Motor commands ---
    if motors is not None:
        motor_cmd_cols = {f"control[{i}]": f"motor_cmd_{i}" for i in range(MOTOR_COUNT)}
        present = [c for c in motor_cmd_cols if c in motors.columns]
        if present:
            base = _asof(base, motors, present, motor_cmd_cols)
    for i in range(MOTOR_COUNT):
        if f"motor_cmd_{i}" not in base.columns:
            base[f"motor_cmd_{i}"] = float("nan")

    # --- ESC status (preferred RPM source) ---
    has_esc = esc is not None
    rpm_source = "esc_status" if has_esc else "unavailable"

    if has_esc:
        # Collect all available ESC fields in one merge pass.
        esc_renames: dict = {}
        esc_cols: list = []
        for i in range(MOTOR_COUNT):
            for field, dest in [
                (f"esc[{i}].esc_rpm", f"motor_rpm_{i}"),
                (f"esc[{i}].esc_voltage", f"motor_voltage_v_{i}"),
                (f"esc[{i}].esc_current", f"motor_current_a_{i}"),
            ]:
                if field in esc.columns:
                    esc_cols.append(field)
                    esc_renames[field] = dest
        if esc_cols:
            base = _asof(base, esc, esc_cols, esc_renames)

    for i in range(MOTOR_COUNT):
        for col in (f"motor_rpm_{i}", f"motor_voltage_v_{i}", f"motor_current_a_{i}"):
            if col not in base.columns:
                base[col] = float("nan")

    # --- Derived: dynamic viscosity from temperature ---
    base["dynamic_viscosity_pa_s"] = base["air_temperature_c"].apply(sutherland_viscosity)

    # --- Derived: Euler angles ---
    q_arr = base[["q_w", "q_x", "q_y", "q_z"]].to_numpy(dtype=float)
    roll, pitch, yaw = quat_to_euler_rad(q_arr)
    base["roll_rad"] = roll
    base["pitch_rad"] = pitch
    base["yaw_rad"] = yaw

    # --- Per-rotor nondimensional state ---
    air_vel_ned = (
        base[["velocity_n_m_s", "velocity_e_m_s", "velocity_d_m_s"]].to_numpy(dtype=float)
        - base[["wind_n_m_s", "wind_e_m_s", "wind_d_m_s"]].to_numpy(dtype=float)
    )
    omega_body = base[["roll_rate_rad_s", "pitch_rate_rad_s", "yaw_rate_rad_s"]].to_numpy(dtype=float)
    rho_arr = base["air_density_kg_m3"].to_numpy(dtype=float)
    mu_arr = base["dynamic_viscosity_pa_s"].to_numpy(dtype=float)

    for i in range(MOTOR_COUNT):
        rpm_arr = base[f"motor_rpm_{i}"].to_numpy(dtype=float)
        state = compute_rotor_inflow_state(
            i, air_vel_ned, omega_body, q_arr, rpm_arr, rho_arr, mu_arr
        )
        base[f"v_inf_{i}_m_s"]    = state["v_inf"]
        base[f"v_normal_{i}_m_s"] = state["v_normal"]
        base[f"v_inplane_{i}_m_s"] = state["v_inplane"]
        base[f"alpha_disk_{i}_rad"] = state["alpha_disk"]
        base[f"j_{i}"]    = state["j"]
        base[f"j_n_{i}"]  = state["j_n"]
        base[f"j_p_{i}"]  = state["j_p"]
        base[f"re_07_{i}"] = state["re_07"]

    # --- Provenance flags ---
    base["has_wind"] = has_wind
    base["has_air_data"] = has_air_data
    base["has_battery_status"] = has_battery
    base["has_esc_status"] = has_esc
    base["rpm_source"] = rpm_source

    # --- Timestamps ---
    base["timestamp_us"] = base["timestamp"]
    base["time_s"] = (base["timestamp_us"] - base["timestamp_us"].iloc[0]) / 1e6

    # --- Sample validity ---
    base["sample_valid"] = base[[
        "velocity_n_m_s", "velocity_e_m_s", "velocity_d_m_s",
        "q_w", "q_x", "q_y", "q_z",
    ]].notna().all(axis=1)

    meta.update({
        "has_wind": has_wind,
        "has_air_data": has_air_data,
        "has_battery": has_battery,
        "has_esc": has_esc,
        "rpm_source": rpm_source,
    })

    return base, meta


# ---------------------------------------------------------------------------
# Output helpers.
# ---------------------------------------------------------------------------

# Canonical column order for the CSV.
_FIXED_COLS = [
    "timestamp_us", "time_s",
    "position_n_m", "position_e_m", "position_d_m",
    "velocity_n_m_s", "velocity_e_m_s", "velocity_d_m_s",
    "q_w", "q_x", "q_y", "q_z",
    "roll_rad", "pitch_rad", "yaw_rad",
    "roll_rate_rad_s", "pitch_rate_rad_s", "yaw_rate_rad_s",
    "wind_n_m_s", "wind_e_m_s", "wind_d_m_s",
    "air_density_kg_m3", "air_temperature_c", "dynamic_viscosity_pa_s",
    "battery_voltage_v", "battery_current_a", "battery_remaining_frac",
] + [
    col
    for i in range(MOTOR_COUNT)
    for col in (f"motor_cmd_{i}", f"motor_rpm_{i}",
                f"motor_voltage_v_{i}", f"motor_current_a_{i}")
] + [
    col
    for i in range(MOTOR_COUNT)
    for col in (f"v_inf_{i}_m_s", f"v_normal_{i}_m_s", f"v_inplane_{i}_m_s",
                f"alpha_disk_{i}_rad", f"j_{i}", f"j_n_{i}", f"j_p_{i}", f"re_07_{i}")
] + [
    "has_wind", "has_air_data", "has_battery_status", "has_esc_status",
    "rpm_source", "sample_valid",
]


def ordered_columns(df: pd.DataFrame) -> list:
    """Return _FIXED_COLS filtered to columns that exist in df."""
    present = set(df.columns)
    return [c for c in _FIXED_COLS if c in present]


def write_csv(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df[ordered_columns(df)].to_csv(output_path, index=False)


def write_report(df: pd.DataFrame, meta: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ts = df["timestamp_us"]
    monotonic = bool((ts.diff().dropna() > 0).all())

    q = df[["q_w", "q_x", "q_y", "q_z"]].to_numpy(dtype=float)
    q_norm = np.linalg.norm(q, axis=1)
    q_ok = int(np.sum(np.abs(q_norm - 1.0) <= 0.01))
    q_bad = len(df) - q_ok

    def _topic_row(name: str) -> str:
        status = "found" if name in meta["topics_found"] else "**missing**"
        return f"| `{name}` | {status} |"

    sections: list = [
        "# x500 Baseline Timeseries — Validation Report",
        "",
        "**Scope:** PX4/Gazebo x500 SITL simulation only.  "
        "This report does **not** validate aerodynamic truth, "
        "hardware equivalence, or final BEMT accuracy.",
        "",
        "## Input",
        f"- File: `{meta['input_file']}`",
        "",
        "## Topic availability",
        "| Topic | Status |",
        "|---|---|",
        *[_topic_row(t) for t in [
            "vehicle_local_position", "vehicle_attitude", "vehicle_angular_velocity",
            "vehicle_air_data", "wind", "battery_status", "actuator_motors", "esc_status",
        ]],
        "",
        "## Dataset summary",
        f"- Total rows: {len(df)}",
        f"- Timestamp range: {ts.iloc[0]:.0f} – {ts.iloc[-1]:.0f} µs "
        f"({df['time_s'].iloc[-1]:.3f} s elapsed)",
        f"- Timestamps monotonically increasing: {'yes' if monotonic else '**NO — non-monotonic timestamps detected**'}",
        f"- Rows with `sample_valid = True`: {int(df['sample_valid'].sum())}",
        "",
        "## Quaternion quality",
        f"- Samples with |‖q‖ − 1| ≤ 0.01: {q_ok} / {len(df)}",
        *([] if q_bad == 0 else [f"- **Warning:** {q_bad} sample(s) outside quaternion norm tolerance."]),
        "",
        "## RPM availability",
        f"- Source: `{meta['rpm_source']}`",
    ]

    for i in range(MOTOR_COUNT):
        col = f"motor_rpm_{i}"
        if col in df.columns:
            n_valid = int(df[col].notna().sum())
            sections.append(f"- Motor {i}: {n_valid} / {len(df)} samples with finite RPM")

    sections += ["", "## Per-rotor nondimensional variable ranges (finite samples only)"]
    for i in range(MOTOR_COUNT):
        sections.append(f"\n### Motor {i}")
        for label, col in [("J", f"j_{i}"), ("J_n", f"j_n_{i}"),
                            ("J_p", f"j_p_{i}"), ("Re_07", f"re_07_{i}")]:
            finite = df[col].dropna()
            if len(finite):
                sections.append(
                    f"- {label}: min={finite.min():.4g}, max={finite.max():.4g}, "
                    f"mean={finite.mean():.4g}  (n={len(finite)})"
                )
            else:
                sections.append(f"- {label}: no finite samples")

    sections += [
        "",
        "## Assumptions and missing-data policy",
        "- **Wind:** "
        + ("measured from `wind` topic."
           if meta["has_wind"]
           else "**topic absent — set to zero**. "
                "Relative air velocity equals vehicle velocity. "
                "Inflow computed under no-wind assumption."),
        "- **Air density:** "
        + ("measured from `vehicle_air_data.rho`."
           if meta["has_air_data"]
           else f"**topic absent — ISA sea-level fallback ({_ISA_RHO} kg/m³) used.**"),
        "- **Dynamic viscosity:** derived from air temperature via Sutherland's law "
        "(engineering assumption; not a direct measurement).",
        "- **Wind vertical component (`wind_d_m_s`):** assumed zero (not logged by PX4).",
        "- **Motor RPM:** "
        + ("measured from `esc_status`. Not estimated or fabricated."
           if meta["has_esc"]
           else "**topic absent — J, J_n, J_p, Re_07 are NaN.** "
                "RPM was not estimated or fabricated from command signals."),
        "",
        "## Limitations",
        "- Results apply only to PX4/Gazebo x500 simulation conditions.",
        "- This report does not validate aerodynamic model truth.",
        "- This report does not establish hardware equivalence with any real platform.",
        "- ESC-reported RPM reflects simulated motor dynamics, not validated real-hardware RPM.",
        "- Single-shot induced-velocity correction in `bemt_model.cpp` is an engineering "
        "approximation (see code comments); it is not used in this offline tool.",
        "",
        "---",
        "_Generated by `tools/ulog_to_csv.py`_",
    ]

    output_path.write_text("\n".join(sections), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a PX4 SITL .ulg log to a flat CSV timeseries and "
            "a markdown validation report."
        )
    )
    parser.add_argument("--ulog", type=Path, required=True,
                        help="Input PX4 .ulg file.")
    parser.add_argument("--output-dir", type=Path, default=Path("output"),
                        help="Output directory (default: output/).")
    return parser.parse_args()


def _run(ulog_path: Path, output_dir: Path) -> int:
    if not ulog_path.exists():
        print(f"Error: {ulog_path} not found.", file=sys.stderr)
        return 1
    try:
        df, meta = build_timeseries(ulog_path)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    csv_path = output_dir / "x500_baseline_timeseries.csv"
    report_path = output_dir / "x500_baseline_validation_report.md"

    write_csv(df, csv_path)
    write_report(df, meta, report_path)

    print(f"Wrote {len(df)} rows to {csv_path}")
    print(f"Wrote validation report to {report_path}")
    return 0


def main() -> int:
    args = _parse_args()
    return _run(args.ulog, args.output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
