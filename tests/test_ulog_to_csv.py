from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ulog_to_csv import (
    CHORD_07_M,
    MOTOR_COUNT,
    PROP_DIAMETER_M,
    _FIXED_COLS,
    _RADIUS_07_M,
    _SUTH_C,
    _SUTH_MU_REF,
    _SUTH_T_REF_K,
    build_timeseries,
    compute_rotor_inflow_state,
    ordered_columns,
    quat_to_euler_rad,
    rotate_ned_to_body,
    sutherland_viscosity,
    _merge_bemt_legacy,
    _merge_bemt_new_cols,
)

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

_ISA_RHO = 1.225
_ISA_MU = sutherland_viscosity(15.0)


def _inflow(motor_idx, v_ned_3, *, omega_3=None, q_4=None, rpm=6000.0, rho=None, mu=None):
    """Single-sample wrapper around compute_rotor_inflow_state.

    Returns the raw dict with shape-(1,) arrays.  Access values as state["key"][0].
    """
    v = np.array([v_ned_3], dtype=float)
    o = np.array([omega_3 if omega_3 is not None else [0.0, 0.0, 0.0]])
    q = np.array([q_4 if q_4 is not None else [1.0, 0.0, 0.0, 0.0]])
    r = np.array([rpm])
    rh = np.array([rho if rho is not None else _ISA_RHO])
    mu_ = np.array([mu if mu is not None else _ISA_MU])
    return compute_rotor_inflow_state(motor_idx, v, o, q, r, rh, mu_)


# ---------------------------------------------------------------------------
# Sutherland viscosity
# ---------------------------------------------------------------------------

class TestSutherlandViscosity:
    def test_at_reference_temperature(self):
        """At 0 °C (= _SUTH_T_REF_K in kelvin) the ratio is 1, returning mu_ref."""
        result = sutherland_viscosity(0.0)
        assert math.isclose(result, _SUTH_MU_REF, rel_tol=1e-9)

    def test_isa_sea_level_in_physical_range(self):
        """At ISA 15 °C, dynamic viscosity should be in 1.7–1.9 × 10⁻⁵ Pa·s."""
        mu = sutherland_viscosity(15.0)
        assert 1.7e-5 < mu < 1.9e-5

    def test_extreme_cold_fallback(self):
        """Temperature at or below absolute zero must return the reference value."""
        assert sutherland_viscosity(-300.0) == _SUTH_MU_REF

    def test_viscosity_increases_with_temperature(self):
        """Sutherland's law: air viscosity increases with temperature."""
        assert sutherland_viscosity(50.0) > sutherland_viscosity(15.0)

    def test_sutherland_formula_at_25c(self):
        """Spot-check against hand-computed value at 25 °C."""
        T_k = 25.0 + 273.15
        ratio = T_k / _SUTH_T_REF_K
        expected = _SUTH_MU_REF * (ratio ** 1.5) * (_SUTH_T_REF_K + _SUTH_C) / (T_k + _SUTH_C)
        assert math.isclose(sutherland_viscosity(25.0), expected, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# rotate_ned_to_body
# ---------------------------------------------------------------------------

class TestRotateNedToBody:
    def test_identity_quaternion_passthrough(self):
        """Identity q → body frame equals NED frame; v_body == v_ned."""
        q = np.array([1.0, 0.0, 0.0, 0.0])
        v = np.array([3.0, -2.0, 1.0])
        np.testing.assert_allclose(rotate_ned_to_body(q, v), v, atol=1e-10)

    def test_yaw90_north_to_body_minus_y(self):
        """Heading east (90° yaw): NED north maps to body −y (left)."""
        s = 1.0 / math.sqrt(2.0)
        q = np.array([s, 0.0, 0.0, s])       # q = [cos45, 0, 0, sin45]
        v_ned = np.array([1.0, 0.0, 0.0])    # NED north unit vector
        v_body = rotate_ned_to_body(q, v_ned)
        np.testing.assert_allclose(v_body, [0.0, -1.0, 0.0], atol=1e-10)

    def test_vectorised_shape_preserved(self):
        """Batch of N samples must return shape (N, 3)."""
        N = 12
        rng = np.random.default_rng(42)
        q = np.tile([1.0, 0.0, 0.0, 0.0], (N, 1))
        v = rng.uniform(-5.0, 5.0, (N, 3))
        result = rotate_ned_to_body(q, v)
        assert result.shape == (N, 3)
        np.testing.assert_allclose(result, v, atol=1e-10)

    def test_rotation_preserves_vector_magnitude(self):
        """Rotation must not change vector magnitude."""
        s = 1.0 / math.sqrt(2.0)
        q = np.array([s, s * 0.5, s * 0.5, 0.0])
        q /= np.linalg.norm(q)
        v = np.array([3.0, 4.0, 5.0])
        v_body = rotate_ned_to_body(q, v)
        assert math.isclose(np.linalg.norm(v_body), np.linalg.norm(v), rel_tol=1e-9)


# ---------------------------------------------------------------------------
# quat_to_euler_rad
# ---------------------------------------------------------------------------

class TestEulerAngles:
    def test_identity_all_zero(self):
        """Identity quaternion → roll = pitch = yaw = 0."""
        q = np.array([[1.0, 0.0, 0.0, 0.0]])
        roll, pitch, yaw = quat_to_euler_rad(q)
        assert roll[0] == pytest.approx(0.0, abs=1e-10)
        assert pitch[0] == pytest.approx(0.0, abs=1e-10)
        assert yaw[0] == pytest.approx(0.0, abs=1e-10)

    def test_90deg_roll(self):
        """q for 90° roll about body x-axis → roll = π/2."""
        s = 1.0 / math.sqrt(2.0)
        q = np.array([[s, s, 0.0, 0.0]])
        roll, pitch, yaw = quat_to_euler_rad(q)
        assert roll[0] == pytest.approx(math.pi / 2.0, abs=1e-6)
        assert pitch[0] == pytest.approx(0.0, abs=1e-6)
        assert yaw[0] == pytest.approx(0.0, abs=1e-6)

    def test_90deg_yaw(self):
        """q for 90° yaw about z-axis → yaw = π/2."""
        s = 1.0 / math.sqrt(2.0)
        q = np.array([[s, 0.0, 0.0, s]])
        roll, pitch, yaw = quat_to_euler_rad(q)
        assert yaw[0] == pytest.approx(math.pi / 2.0, abs=1e-6)
        assert roll[0] == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Hover (zero velocity, zero angular rates, identity attitude)
# ---------------------------------------------------------------------------

class TestHover:
    def setup_method(self):
        self.state = _inflow(0, [0.0, 0.0, 0.0])

    def test_v_inf_zero(self):
        assert self.state["v_inf"][0] == pytest.approx(0.0, abs=1e-10)

    def test_j_zero(self):
        assert self.state["j"][0] == pytest.approx(0.0, abs=1e-10)

    def test_j_n_zero(self):
        assert self.state["j_n"][0] == pytest.approx(0.0, abs=1e-10)

    def test_re07_positive(self):
        """Even in hover the blade moves, so Re_07 > 0."""
        assert self.state["re_07"][0] > 0.0

    def test_re07_matches_hover_formula(self):
        """Re_07 = rho * omega_prop * 0.7R * c_07 / mu at hover (v_inf = 0)."""
        rpm = 6000.0
        omega = rpm * (2.0 * math.pi / 60.0)
        v_tip07 = omega * _RADIUS_07_M
        expected = _ISA_RHO * v_tip07 * CHORD_07_M / _ISA_MU
        assert self.state["re_07"][0] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Edgewise forward flight (level, heading north, identity attitude)
# ---------------------------------------------------------------------------

class TestEdgewiseFlight:
    V_FWD = 10.0
    RPM = 6000.0

    def setup_method(self):
        self.state = _inflow(0, [self.V_FWD, 0.0, 0.0], rpm=self.RPM)

    def test_v_inf_equals_v_forward(self):
        assert self.state["v_inf"][0] == pytest.approx(self.V_FWD, rel=1e-6)

    def test_v_normal_zero(self):
        """Level flight has no axial velocity component."""
        assert self.state["v_normal"][0] == pytest.approx(0.0, abs=1e-10)

    def test_v_inplane_equals_v_forward(self):
        assert self.state["v_inplane"][0] == pytest.approx(self.V_FWD, rel=1e-6)

    def test_alpha_disk_pi_over_2(self):
        """Fully in-plane inflow: disk incidence = 90°."""
        assert self.state["alpha_disk"][0] == pytest.approx(math.pi / 2.0, rel=1e-6)

    def test_j_matches_thesis_formula(self):
        """J = V_inf / (n * D), the thesis baseline equation."""
        n = self.RPM / 60.0
        expected = self.V_FWD / (n * PROP_DIAMETER_M)
        assert self.state["j"][0] == pytest.approx(expected, rel=1e-6)

    def test_j_n_zero_in_level_flight(self):
        """J_n = V_normal / (nD) = 0 when no axial inflow."""
        assert self.state["j_n"][0] == pytest.approx(0.0, abs=1e-10)

    def test_j_p_equals_j_in_edgewise_flight(self):
        """In purely edgewise flight all inflow is in-plane, so J_p == J."""
        assert self.state["j_p"][0] == pytest.approx(self.state["j"][0], rel=1e-6)

    def test_re07_positive(self):
        assert self.state["re_07"][0] > 0.0


# ---------------------------------------------------------------------------
# Descent (NED vz positive = downward = axial headwind for rotor)
# ---------------------------------------------------------------------------

class TestDescentFlight:
    V_DESC = 2.0
    RPM = 6000.0

    def setup_method(self):
        # NED vz = +V_DESC → vehicle moves downward
        self.state = _inflow(0, [0.0, 0.0, self.V_DESC], rpm=self.RPM)

    def test_v_normal_positive(self):
        """Descent creates upward inflow through the rotor (positive axial component)."""
        assert self.state["v_normal"][0] > 0.0

    def test_v_normal_equals_v_descent(self):
        assert self.state["v_normal"][0] == pytest.approx(self.V_DESC, rel=1e-6)

    def test_v_inplane_zero(self):
        """No horizontal motion → no in-plane inflow."""
        assert self.state["v_inplane"][0] == pytest.approx(0.0, abs=1e-10)

    def test_alpha_disk_zero(self):
        """Purely axial inflow → disk incidence = 0."""
        assert self.state["alpha_disk"][0] == pytest.approx(0.0, abs=1e-10)

    def test_j_n_positive(self):
        """Axial headwind in descent makes J_n positive."""
        assert self.state["j_n"][0] > 0.0

    def test_j_n_matches_formula(self):
        n = self.RPM / 60.0
        expected = self.V_DESC / (n * PROP_DIAMETER_M)
        assert self.state["j_n"][0] == pytest.approx(expected, rel=1e-6)

    def test_j_equals_abs_j_n_in_pure_axial(self):
        """In purely axial flight, J == |J_n| because v_inplane == 0."""
        assert abs(self.state["j"][0]) == pytest.approx(abs(self.state["j_n"][0]), rel=1e-6)


# ---------------------------------------------------------------------------
# Signed axial convention in climb
# ---------------------------------------------------------------------------

class TestClimbFlight:
    V_CLIMB = -2.0
    RPM = 6000.0

    def test_v_normal_and_j_n_are_negative_in_climb(self):
        """Upward climb (negative NED vz) produces negative signed axial inflow."""
        state = _inflow(0, [0.0, 0.0, self.V_CLIMB], rpm=self.RPM)
        assert state["v_normal"][0] < 0.0
        assert state["j_n"][0] < 0.0


# ---------------------------------------------------------------------------
# Missing / invalid RPM → NaN propagation for RPM-dependent outputs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad_rpm", [float("nan"), 0.0, -100.0, -1.0])
class TestMissingRpm:
    def test_j_is_nan(self, bad_rpm):
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert not np.isfinite(state["j"][0])

    def test_j_n_is_nan(self, bad_rpm):
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert not np.isfinite(state["j_n"][0])

    def test_j_p_is_nan(self, bad_rpm):
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert not np.isfinite(state["j_p"][0])

    def test_re07_is_nan(self, bad_rpm):
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert not np.isfinite(state["re_07"][0])

    def test_v_inf_still_finite(self, bad_rpm):
        """v_inf is independent of RPM and must remain finite."""
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert np.isfinite(state["v_inf"][0])

    def test_alpha_disk_still_finite(self, bad_rpm):
        """alpha_disk is independent of RPM and must remain finite."""
        state = _inflow(0, [5.0, 0.0, 0.0], rpm=bad_rpm)
        assert np.isfinite(state["alpha_disk"][0])


# ---------------------------------------------------------------------------
# Quaternion input quality
# ---------------------------------------------------------------------------

class TestQuaternionInputs:
    def test_unit_quaternion_produces_finite_output(self):
        state = _inflow(0, [0.0, 0.0, 0.0])
        assert np.isfinite(state["v_inf"][0])

    def test_all_motors_symmetric_at_hover(self):
        """All four rotors should give identical results at hover with identity attitude."""
        results = [_inflow(i, [0.0, 0.0, 0.0]) for i in range(MOTOR_COUNT)]
        ref = results[0]["re_07"][0]
        for state in results[1:]:
            assert state["re_07"][0] == pytest.approx(ref, rel=1e-6)

    def test_forward_flight_all_motors_same_v_inf(self):
        """All rotors share the same free-stream speed in pure translation (no rotation)."""
        states = [_inflow(i, [10.0, 0.0, 0.0]) for i in range(MOTOR_COUNT)]
        ref = states[0]["v_inf"][0]
        for st in states[1:]:
            assert st["v_inf"][0] == pytest.approx(ref, rel=1e-6)

    def test_zero_vector_inflow_is_handled(self):
        """All-zero velocity should not raise and should return finite values."""
        state = _inflow(0, [0.0, 0.0, 0.0])
        for key in ("v_inf", "v_normal", "v_inplane", "alpha_disk"):
            assert np.isfinite(state[key][0])


# ---------------------------------------------------------------------------
# Higher-level build_timeseries coverage
# ---------------------------------------------------------------------------

class TestBuildTimeseries:
    def test_build_timeseries_merges_topics_and_preserves_missing_topic_policy(self, monkeypatch, tmp_path):
        sentinel_ulog = object()

        topics = {
            "vehicle_local_position": pd.DataFrame({
                "timestamp": [1_000_000, 1_020_000],
                "x": [0.0, 0.5],
                "y": [0.0, 0.0],
                "z": [0.0, -0.1],
                "vx": [5.0, 5.0],
                "vy": [0.0, 0.0],
                "vz": [0.0, 0.0],
            }),
            "vehicle_attitude": pd.DataFrame({
                "timestamp": [995_000, 1_015_000],
                "q[0]": [1.0, 1.0],
                "q[1]": [0.0, 0.0],
                "q[2]": [0.0, 0.0],
                "q[3]": [0.0, 0.0],
            }),
            "vehicle_air_data": pd.DataFrame({
                "timestamp": [999_000, 1_019_000],
                "rho": [1.20, 1.21],
                "baro_temp_celcius": [20.0, 21.0],
            }),
            "wind": pd.DataFrame({
                "timestamp": [1_090_000],
                "windspeed_north": [1.5],
                "windspeed_east": [-0.5],
            }),
            "actuator_motors": pd.DataFrame({
                "timestamp": [1_001_000, 1_019_000],
                "control[0]": [0.11, 0.22],
                "control[1]": [0.12, 0.23],
                "control[2]": [0.13, 0.24],
                "control[3]": [0.14, 0.25],
            }),
        }

        def fake_load(path):
            assert path == tmp_path / "sample.ulg"
            return sentinel_ulog

        def fake_topic_df(ulog, name, instance=0):
            assert ulog is sentinel_ulog
            assert instance == 0
            return topics.get(name)

        monkeypatch.setattr("ulog_to_csv._load_ulog", fake_load)
        monkeypatch.setattr("ulog_to_csv._topic_df", fake_topic_df)

        df, meta = build_timeseries(tmp_path / "sample.ulg")

        assert df["timestamp_us"].tolist() == [1_000_000, 1_020_000]
        assert df["time_s"].tolist() == pytest.approx([0.0, 0.02])

        np.testing.assert_allclose(df["q_w"], [1.0, 1.0])
        np.testing.assert_allclose(df["air_density_kg_m3"], [1.20, 1.21])
        np.testing.assert_allclose(df["wind_n_m_s"], [1.5, 1.5])
        np.testing.assert_allclose(df["wind_e_m_s"], [-0.5, -0.5])
        np.testing.assert_allclose(df["motor_cmd_0"], [0.11, 0.22])

        assert "vehicle_angular_velocity" in meta["topics_missing"]
        assert "battery_status" in meta["topics_missing"]
        assert "esc_status" in meta["topics_missing"]
        assert "wind" in meta["topics_found"]

        assert (df["has_wind"] == True).all()
        assert (df["has_air_data"] == True).all()
        assert (df["has_battery_status"] == False).all()
        assert (df["has_esc_status"] == False).all()
        assert (df["rpm_source"] == "unavailable").all()
        assert (df["rotor_state_source"] == "reconstruction").all()
        assert (df["has_bemt_topic"] == False).all()  # noqa: E712

        np.testing.assert_allclose(df["roll_rate_rad_s"], [0.0, 0.0])
        assert df["battery_voltage_v"].isna().all()
        assert df["motor_rpm_0"].isna().all()
        assert df["j_0"].isna().all()
        assert df["j_n_0"].isna().all()
        assert df["re_07_0"].isna().all()
        assert np.isfinite(df["v_inf_0_m_s"]).all()
        assert df["sample_valid"].all()


# ---------------------------------------------------------------------------
# Timestamp monotonicity (synthetic DataFrame checks)
# ---------------------------------------------------------------------------

class TestTimestampMonotonicity:
    def test_ascending_timestamps_are_monotonic(self):
        ts = np.arange(0, 10_000_000, 100_000, dtype=np.int64)
        df = pd.DataFrame({"timestamp_us": ts})
        assert bool((df["timestamp_us"].diff().dropna() > 0).all())

    def test_shuffled_timestamps_not_monotonic(self):
        ts = np.array([0, 200_000, 100_000, 300_000], dtype=np.int64)
        df = pd.DataFrame({"timestamp_us": ts})
        assert not bool((df["timestamp_us"].diff().dropna() > 0).all())

    def test_duplicate_timestamps_not_strictly_monotonic(self):
        ts = np.array([0, 100_000, 100_000, 200_000], dtype=np.int64)
        df = pd.DataFrame({"timestamp_us": ts})
        assert not bool((df["timestamp_us"].diff().dropna() > 0).all())


# ---------------------------------------------------------------------------
# Schema and column ordering
# ---------------------------------------------------------------------------

class TestSchema:
    def _minimal_df(self):
        cols = [
            "timestamp_us", "time_s",
            "velocity_n_m_s", "velocity_e_m_s", "velocity_d_m_s",
            "q_w", "q_x", "q_y", "q_z",
            "has_wind", "rpm_source", "sample_valid",
        ]
        return pd.DataFrame({c: [0.0] for c in cols})

    def test_ordered_columns_returns_list(self):
        df = self._minimal_df()
        assert isinstance(ordered_columns(df), list)

    def test_ordered_columns_only_present_cols(self):
        """ordered_columns must not include columns absent from the DataFrame."""
        df = self._minimal_df()
        cols = ordered_columns(df)
        df_cols = set(df.columns)
        assert all(c in df_cols for c in cols)

    def test_velocity_columns_have_si_suffix(self):
        vel_cols = [c for c in _FIXED_COLS if "velocity" in c and "wind" not in c]
        assert all(c.endswith("_m_s") for c in vel_cols), (
            f"Non-SI velocity columns: {[c for c in vel_cols if not c.endswith('_m_s')]}"
        )

    def test_motor_rpm_columns_present_for_all_motors(self):
        rpm_cols = [c for c in _FIXED_COLS if "motor_rpm" in c]
        assert len(rpm_cols) == MOTOR_COUNT
        for i in range(MOTOR_COUNT):
            assert f"motor_rpm_{i}" in _FIXED_COLS

    def test_advance_ratio_columns_present_for_all_motors(self):
        j_cols = [c for c in _FIXED_COLS if c.startswith("j_") and not c.startswith("j_n") and not c.startswith("j_p")]
        assert len(j_cols) == MOTOR_COUNT


# ---------------------------------------------------------------------------
# Provenance and assumptions
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_has_esc_status_flag_in_schema(self):
        assert "has_esc_status" in _FIXED_COLS

    def test_rpm_source_column_in_schema(self):
        assert "rpm_source" in _FIXED_COLS

    def test_wind_d_assumed_zero_in_schema(self):
        """Vertical wind column must be in the fixed schema (set to 0 when absent from log)."""
        assert "wind_d_m_s" in _FIXED_COLS

    def test_no_aerodynamic_coefficient_columns(self):
        """c_t / c_q / c_p must not appear — aerodynamic validation is out of scope."""
        bad = [c for c in _FIXED_COLS if c.startswith(("c_t", "c_q", "c_p", "ct_", "cq_", "cp_"))]
        assert bad == [], f"Unexpected coefficient columns: {bad}"

    def test_all_per_rotor_columns_present(self):
        for i in range(MOTOR_COUNT):
            for col in (f"j_{i}", f"j_n_{i}", f"j_p_{i}", f"re_07_{i}"):
                assert col in _FIXED_COLS, f"Missing column: {col}"

    def test_sample_valid_column_in_schema(self):
        assert "sample_valid" in _FIXED_COLS

    def test_rotor_state_source_in_schema(self):
        assert "rotor_state_source" in _FIXED_COLS

    def test_has_bemt_topic_in_schema(self):
        assert "has_bemt_topic" in _FIXED_COLS

    def test_kinematic_full_columns_in_schema(self):
        for i in range(MOTOR_COUNT):
            for col in (
                f"kinematic_signed_axial_speed_{i}_m_s",
                f"kinematic_v_normal_{i}_m_s",
                f"kinematic_v_inf_{i}_m_s",
                f"kinematic_j_{i}",
                f"kinematic_j_n_{i}",
                f"kinematic_j_p_{i}",
            ):
                assert col in _FIXED_COLS, f"Missing: {col}"

    def test_corrected_columns_in_schema(self):
        for i in range(MOTOR_COUNT):
            for col in (
                f"corrected_v_inf_{i}_m_s",
                f"corrected_j_{i}",
                f"corrected_j_n_{i}",
                f"corrected_j_p_{i}",
                f"induced_axial_velocity_{i}_m_s",
            ):
                assert col in _FIXED_COLS, f"Missing: {col}"


# ---------------------------------------------------------------------------
# BEMT topic: source preference and column presence
# ---------------------------------------------------------------------------

def _make_base_topics():
    """Minimal topics dict for build_timeseries that satisfies required topics."""
    return {
        "vehicle_local_position": pd.DataFrame({
            "timestamp": [1_000_000, 1_020_000],
            "x": [0.0, 0.5], "y": [0.0, 0.0], "z": [0.0, -0.1],
            "vx": [5.0, 5.0], "vy": [0.0, 0.0], "vz": [0.0, 0.0],
        }),
        "vehicle_attitude": pd.DataFrame({
            "timestamp": [995_000, 1_015_000],
            "q[0]": [1.0, 1.0], "q[1]": [0.0, 0.0],
            "q[2]": [0.0, 0.0], "q[3]": [0.0, 0.0],
        }),
        "vehicle_air_data": pd.DataFrame({
            "timestamp": [999_000, 1_019_000],
            "rho": [1.20, 1.21], "baro_temp_celcius": [20.0, 21.0],
        }),
    }


def _make_bemt_topic():
    """Minimal bemt_rotor_state DataFrame for two timestamps."""
    data: dict = {"timestamp": [1_000_000, 1_020_000]}
    for i in range(MOTOR_COUNT):
        data[f"motor_rpm[{i}]"] = [6000.0, 6100.0]
        data[f"kinematic_signed_axial_speed_m_s[{i}]"] = [-0.10, -0.12]
        data[f"kinematic_v_normal_m_s[{i}]"] = [0.10, 0.12]
        data[f"kinematic_v_inplane_m_s[{i}]"] = [5.0, 5.1]
        data[f"kinematic_v_inf_m_s[{i}]"] = [5.001, 5.101]
        data[f"kinematic_alpha_disk_rad[{i}]"] = [1.55, 1.54]
        data[f"kinematic_j[{i}]"] = [0.077, 0.078]
        data[f"kinematic_j_n[{i}]"] = [-0.002, -0.002]
        data[f"kinematic_j_p[{i}]"] = [0.077, 0.078]
        data[f"induced_axial_velocity_m_s[{i}]"] = [2.1, 2.2]
        data[f"corrected_signed_axial_speed_m_s[{i}]"] = [2.0, 2.08]
        data[f"corrected_v_normal_m_s[{i}]"] = [2.0, 2.08]
        data[f"corrected_v_inplane_m_s[{i}]"] = [5.0, 5.1]
        data[f"corrected_v_inf_m_s[{i}]"] = [5.39, 5.49]
        data[f"corrected_alpha_disk_rad[{i}]"] = [1.19, 1.18]
        data[f"corrected_j[{i}]"] = [0.083, 0.084]
        data[f"corrected_j_n[{i}]"] = [0.031, 0.032]
        data[f"corrected_j_p[{i}]"] = [0.077, 0.078]
        data[f"re_07[{i}]"] = [55_000.0, 56_000.0]
    return pd.DataFrame(data)


def _patch_topics(monkeypatch, topics):
    sentinel = object()

    def fake_load(path):
        return sentinel

    def fake_topic_df(ulog, name, instance=0):
        assert ulog is sentinel
        return topics.get(name)

    monkeypatch.setattr("ulog_to_csv._load_ulog", fake_load)
    monkeypatch.setattr("ulog_to_csv._topic_df", fake_topic_df)


class TestBuildTimeseriesWithBemtTopic:
    """build_timeseries preferred path: bemt_rotor_state topic is present."""

    def test_rotor_state_source_is_bemt_topic(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, meta = build_timeseries(tmp_path / "sample.ulg")
        assert (df["rotor_state_source"] == "bemt_rotor_state").all()
        assert meta["rotor_state_source"] == "bemt_rotor_state"
        assert meta["has_bemt_topic"] is True

    def test_has_bemt_topic_flag_true(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        assert (df["has_bemt_topic"] == True).all()  # noqa: E712

    def test_kinematic_columns_present(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        for i in range(MOTOR_COUNT):
            assert f"kinematic_v_inf_{i}_m_s" in df.columns
            assert f"kinematic_j_{i}" in df.columns
            assert f"kinematic_j_n_{i}" in df.columns

    def test_corrected_columns_present(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        for i in range(MOTOR_COUNT):
            assert f"corrected_v_inf_{i}_m_s" in df.columns
            assert f"corrected_j_{i}" in df.columns
            assert f"induced_axial_velocity_{i}_m_s" in df.columns

    def test_legacy_columns_populated_from_kinematic_signed_axial(self, monkeypatch, tmp_path):
        """v_normal_*_m_s must carry the SIGNED axial value from kinematic_signed_axial_speed."""
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        bemt = _make_bemt_topic()
        np.testing.assert_allclose(
            df["v_normal_0_m_s"].values,
            bemt["kinematic_signed_axial_speed_m_s[0]"].values,
            rtol=1e-5,
        )

    def test_legacy_j_populated_from_kinematic_j(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        bemt = _make_bemt_topic()
        np.testing.assert_allclose(
            df["j_0"].values,
            bemt["kinematic_j[0]"].values,
            rtol=1e-5,
        )

    def test_bemt_topic_in_topics_found(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        topics["bemt_rotor_state"] = _make_bemt_topic()
        _patch_topics(monkeypatch, topics)
        _, meta = build_timeseries(tmp_path / "sample.ulg")
        assert "bemt_rotor_state" in meta["topics_found"]


class TestBuildTimeseriesFallbackReconstruction:
    """build_timeseries fallback path: bemt_rotor_state topic is absent."""

    def test_rotor_state_source_is_reconstruction(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        _patch_topics(monkeypatch, topics)
        df, meta = build_timeseries(tmp_path / "sample.ulg")
        assert (df["rotor_state_source"] == "reconstruction").all()
        assert meta["rotor_state_source"] == "reconstruction"
        assert meta["has_bemt_topic"] is False

    def test_has_bemt_topic_flag_false(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        assert (df["has_bemt_topic"] == False).all()  # noqa: E712

    def test_corrected_columns_absent(self, monkeypatch, tmp_path):
        """Corrected columns must not appear when only reconstruction is available."""
        topics = _make_base_topics()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        for i in range(MOTOR_COUNT):
            assert f"corrected_v_inf_{i}_m_s" not in df.columns
            assert f"induced_axial_velocity_{i}_m_s" not in df.columns

    def test_legacy_columns_still_present(self, monkeypatch, tmp_path):
        """Legacy columns (v_inf_*, j_*, etc.) must be present from reconstruction."""
        topics = _make_base_topics()
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        for i in range(MOTOR_COUNT):
            assert f"v_inf_{i}_m_s" in df.columns
            assert f"j_{i}" in df.columns

    def test_bemt_topic_in_topics_missing(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        _patch_topics(monkeypatch, topics)
        _, meta = build_timeseries(tmp_path / "sample.ulg")
        assert "bemt_rotor_state" in meta["topics_missing"]


class TestBemtTopicPreference:
    """When bemt_rotor_state is present, it must be preferred over reconstruction."""

    def test_prefers_bemt_topic_values_over_reconstruction(self, monkeypatch, tmp_path):
        """v_inf_0_m_s must come from the BEMT topic, not offline reconstruction."""
        topics = _make_base_topics()
        bemt = _make_bemt_topic()
        topics["bemt_rotor_state"] = bemt
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        np.testing.assert_allclose(
            df["v_inf_0_m_s"].values,
            bemt["kinematic_v_inf_m_s[0]"].values,
            rtol=1e-5,
        )

    def test_re_07_comes_from_bemt_topic(self, monkeypatch, tmp_path):
        topics = _make_base_topics()
        bemt = _make_bemt_topic()
        topics["bemt_rotor_state"] = bemt
        _patch_topics(monkeypatch, topics)
        df, _ = build_timeseries(tmp_path / "sample.ulg")
        np.testing.assert_allclose(
            df["re_07_0"].values,
            bemt["re_07[0]"].values,
            rtol=1e-5,
        )


class TestMergeBemtHelpers:
    """Unit tests for _merge_bemt_legacy and _merge_bemt_new_cols."""

    def _base_df(self):
        return pd.DataFrame({
            "timestamp": [1_000_000, 1_020_000],
            "dummy": [0.0, 0.0],
        })

    def test_legacy_merge_populates_v_normal_from_signed_axial(self):
        base = self._base_df()
        bemt = pd.DataFrame({
            "timestamp": [1_000_000, 1_020_000],
            "kinematic_signed_axial_speed_m_s[0]": [-0.5, -0.6],
            "kinematic_v_inplane_m_s[0]": [3.0, 3.1],
            "kinematic_v_inf_m_s[0]": [3.04, 3.16],
            "kinematic_alpha_disk_rad[0]": [1.4, 1.4],
            "kinematic_j[0]": [0.046, 0.048],
            "kinematic_j_n[0]": [-0.007, -0.009],
            "kinematic_j_p[0]": [0.046, 0.048],
            "re_07[0]": [50_000.0, 51_000.0],
        })
        result = _merge_bemt_legacy(base, bemt)
        np.testing.assert_allclose(result["v_normal_0_m_s"], [-0.5, -0.6], rtol=1e-5)
        np.testing.assert_allclose(result["j_0"], [0.046, 0.048], rtol=1e-5)

    def test_new_cols_merge_adds_kinematic_full_names(self):
        base = self._base_df()
        bemt = pd.DataFrame({
            "timestamp": [1_000_000, 1_020_000],
            "kinematic_v_inf_m_s[0]": [3.0, 3.1],
            "kinematic_j[0]": [0.046, 0.048],
            "kinematic_alpha_disk_rad[0]": [1.4, 1.4],
            "corrected_v_inf_m_s[0]": [4.0, 4.1],
            "corrected_j[0]": [0.06, 0.063],
            "induced_axial_velocity_m_s[0]": [2.0, 2.1],
        })
        result = _merge_bemt_new_cols(base, bemt)
        assert "kinematic_v_inf_0_m_s" in result.columns
        assert "kinematic_j_0" in result.columns
        assert "kinematic_alpha_disk_0_rad" in result.columns
        assert "corrected_v_inf_0_m_s" in result.columns
        assert "induced_axial_velocity_0_m_s" in result.columns
        np.testing.assert_allclose(result["kinematic_v_inf_0_m_s"], [3.0, 3.1], rtol=1e-5)
        np.testing.assert_allclose(result["corrected_j_0"], [0.06, 0.063], rtol=1e-5)
