#include "bemt_model.hpp"

#include <drivers/drv_hrt.h>
#include <px4_platform_common/getopt.h>
#include <px4_platform_common/log.h>
#include <px4_platform_common/module.h>
#include <px4_platform_common/posix.h>
#include <uORB/Subscription.hpp>
#include <uORB/topics/actuator_motors.h>
#include <uORB/topics/battery_status.h>
#include <uORB/topics/esc_status.h>
#include <uORB/topics/vehicle_air_data.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/wind.h>

#include <cmath>

using namespace time_literals;

class BEMTOptimize : public ModuleBase<BEMTOptimize>
{
public:
	BEMTOptimize() = default;
	~BEMTOptimize() override = default;

	static int task_spawn(int argc, char *argv[])
	{
		_task_id = px4_task_spawn_cmd("bemt_optimize",
					      SCHED_DEFAULT,
					      SCHED_PRIORITY_DEFAULT,
					      4000,
					      (px4_main_t)&run_trampoline,
					      argv);

		if (_task_id < 0) {
			_task_id = -1;
			return PX4_ERROR;
		}

		return PX4_OK;
	}

	static BEMTOptimize *instantiate(int argc, char *argv[])
	{
		return new BEMTOptimize();
	}

	static int custom_command(int argc, char *argv[])
	{
		return print_usage("unknown command");
	}

	int print_status() override
	{
		PX4_INFO("running");
		return 0;
	}

	void run() override
	{
		while (!should_exit()) {
			step();
			px4_usleep(100_ms);
		}
	}

	static int print_usage(const char *reason = nullptr)
	{
		if (reason) {
			PX4_WARN("%s\n", reason);
		}

		PRINT_MODULE_DESCRIPTION(
			R"DESCR_STR(
### Description
Runs the current BEMT state-extraction model against PX4 simulation topics and
prints per-rotor inflow, correction, and coefficient summaries for SITL validation.

This module is for observation and validation only. Printed aerodynamic outputs
remain provisional until the thesis pipeline is fully validated.
)DESCR_STR");

		PRINT_MODULE_USAGE_NAME("bemt_optimize", "template");
		PRINT_MODULE_USAGE_COMMAND("start");
		PRINT_MODULE_USAGE_DEFAULT_COMMANDS();
		return 0;
	}

private:
	static constexpr int kQuatSize = 4;
	static constexpr float kMinUsefulRpm = 1.0F;

	uORB::Subscription _vehicle_local_position_sub{ORB_ID(vehicle_local_position)};
	uORB::Subscription _vehicle_attitude_sub{ORB_ID(vehicle_attitude)};
	uORB::Subscription _vehicle_angular_velocity_sub{ORB_ID(vehicle_angular_velocity)};
	uORB::Subscription _vehicle_air_data_sub{ORB_ID(vehicle_air_data)};
	uORB::Subscription _wind_sub{ORB_ID(wind)};
	uORB::Subscription _battery_status_sub{ORB_ID(battery_status)};
	uORB::Subscription _actuator_motors_sub{ORB_ID(actuator_motors)};
	uORB::Subscription _esc_status_sub{ORB_ID(esc_status)};

	hrt_abstime _last_print_time{0};

	static float dynamic_viscosity_from_temperature_c(float temperature_c)
	{
		// Engineering assumption: Sutherland law for air using PX4 air-data temperature.
		static constexpr float kReferenceTempK = 273.15F;
		static constexpr float kReferenceViscosityPaS = 1.716e-5F;
		static constexpr float kSutherlandConstantK = 110.4F;

		const float temperature_k = temperature_c + 273.15F;

		if (temperature_k <= 0.0F) {
			return kReferenceViscosityPaS;
		}

		const float ratio = temperature_k / kReferenceTempK;
		return kReferenceViscosityPaS *
		       std::pow(ratio, 1.5F) *
		       ((kReferenceTempK + kSutherlandConstantK) /
			(temperature_k + kSutherlandConstantK));
	}

	static int count_valid_rpm_sources(const bemt::Input &input)
	{
		int count = 0;

		for (int i = 0; i < bemt::kMotorCount; ++i) {
			if (PX4_ISFINITE(input.motor_rpm[i]) && input.motor_rpm[i] > kMinUsefulRpm) {
				++count;
			}
		}

		return count;
	}

	void step()
	{
		vehicle_local_position_s lpos{};
		vehicle_attitude_s attitude{};
		vehicle_air_data_s air_data{};
		vehicle_angular_velocity_s ang_vel{};

		if (!_vehicle_local_position_sub.copy(&lpos) ||
		    !_vehicle_attitude_sub.copy(&attitude) ||
		    !_vehicle_air_data_sub.copy(&air_data)) {
			return;
		}

		_vehicle_angular_velocity_sub.copy(&ang_vel);

		if (!lpos.v_xy_valid || !lpos.v_z_valid) {
			return;
		}

		wind_s wind{};
		const bool have_wind = _wind_sub.copy(&wind);

		battery_status_s battery{};
		const bool have_battery = _battery_status_sub.copy(&battery);

		actuator_motors_s motors{};
		const bool have_motors = _actuator_motors_sub.copy(&motors);

		esc_status_s esc_status{};
		const bool have_esc_status = _esc_status_sub.copy(&esc_status);

		bemt::Input input{};
		input.air_density_kg_m3 = PX4_ISFINITE(air_data.rho) && (air_data.rho > 0.0F) ? air_data.rho : 1.225F;
		input.dynamic_viscosity_pa_s = dynamic_viscosity_from_temperature_c(
			PX4_ISFINITE(air_data.baro_temp_celcius) ? air_data.baro_temp_celcius : 15.0F);

		input.altitude_m = PX4_ISFINITE(lpos.z) ? -lpos.z : 0.0F;
		input.velocity_north_m_s = PX4_ISFINITE(lpos.vx) ? lpos.vx : 0.0F;
		input.velocity_east_m_s = PX4_ISFINITE(lpos.vy) ? lpos.vy : 0.0F;
		input.velocity_down_m_s = PX4_ISFINITE(lpos.vz) ? lpos.vz : 0.0F;

		input.wind_north_m_s = have_wind && PX4_ISFINITE(wind.windspeed_north) ? wind.windspeed_north : 0.0F;
		input.wind_east_m_s = have_wind && PX4_ISFINITE(wind.windspeed_east) ? wind.windspeed_east : 0.0F;

		// Copy quaternion atomically: if any component is non-finite, fall back to
		// the identity quaternion rather than patching individual components, which
		// would silently produce a non-unit quaternion and corrupt body-frame rotation.
		bool q_valid = true;

		for (int i = 0; i < kQuatSize; ++i) {
			if (!PX4_ISFINITE(attitude.q[i])) {
				q_valid = false;
				break;
			}
		}

		if (q_valid) {
			for (int i = 0; i < kQuatSize; ++i) {
				input.attitude_q[i] = attitude.q[i];
			}

		} else {
			// Defensive fallback: identity quaternion (no rotation, NED == body).
			input.attitude_q[0] = 1.0F;
			input.attitude_q[1] = 0.0F;
			input.attitude_q[2] = 0.0F;
			input.attitude_q[3] = 0.0F;
		}

		input.roll_rate_rad_s = PX4_ISFINITE(ang_vel.xyz[0]) ? ang_vel.xyz[0] : 0.0F;
		input.pitch_rate_rad_s = PX4_ISFINITE(ang_vel.xyz[1]) ? ang_vel.xyz[1] : 0.0F;
		input.yaw_rate_rad_s = PX4_ISFINITE(ang_vel.xyz[2]) ? ang_vel.xyz[2] : 0.0F;

		input.battery_voltage_v = have_battery && PX4_ISFINITE(battery.voltage_v) ? battery.voltage_v : 0.0F;
		input.battery_current_a = have_battery && PX4_ISFINITE(battery.current_a) ? battery.current_a : 0.0F;
		input.battery_remaining = have_battery && PX4_ISFINITE(battery.remaining) ? battery.remaining : 0.0F;

		for (int i = 0; i < bemt::kMotorCount; ++i) {
			input.motor_command[i] = (have_motors && PX4_ISFINITE(motors.control[i])) ? motors.control[i] : 0.0F;
			input.motor_voltage_v[i] = input.battery_voltage_v;
			input.motor_current_a[i] = 0.0F;
			input.motor_rpm[i] = 0.0F;

			if (have_esc_status && (i < esc_status.esc_count)) {
				input.motor_rpm[i] = PX4_ISFINITE(static_cast<float>(esc_status.esc[i].esc_rpm)) ?
						       static_cast<float>(esc_status.esc[i].esc_rpm) : 0.0F;
				input.motor_voltage_v[i] = PX4_ISFINITE(esc_status.esc[i].esc_voltage) ?
							 esc_status.esc[i].esc_voltage : input.battery_voltage_v;
				input.motor_current_a[i] = PX4_ISFINITE(esc_status.esc[i].esc_current) ?
							 esc_status.esc[i].esc_current : 0.0F;
			}
		}

		bemt::Output output{};

		if (!bemt::calculate(input, output)) {
			if (hrt_elapsed_time(&_last_print_time) > 1_s) {
				_last_print_time = hrt_absolute_time();
				PX4_WARN("bemt_optimize: calculate() rejected current input sample");
			}

			return;
		}

		if (hrt_elapsed_time(&_last_print_time) > 1_s) {
			_last_print_time = hrt_absolute_time();

			const int rpm_valid_count = count_valid_rpm_sources(input);
			const char *rpm_source = (have_esc_status && rpm_valid_count > 0) ? "esc_status" : "missing_or_zero";
			const unsigned raw_esc_rpm_m0 =
				(have_esc_status && esc_status.esc_count > 0) ? static_cast<unsigned>(esc_status.esc[0].esc_rpm) : 0U;

			PX4_INFO("[bemt] src=%s valid=%d/%d raw_esc_rpm_m0=%u rpm_m0=%.1f cmd_m0=%.3f batt_v=%.2f",
				 rpm_source,
				 rpm_valid_count,
				 bemt::kMotorCount,
				 raw_esc_rpm_m0,
				 (double)input.motor_rpm[0],
				 (double)input.motor_command[0],
				 (double)input.battery_voltage_v);

			PX4_INFO("[bemt] vx=%.2f vy=%.2f vz=%.2f wind_n=%.2f wind_e=%.2f [m/s NED]",
				 (double)input.velocity_north_m_s,
				 (double)input.velocity_east_m_s,
				 (double)input.velocity_down_m_s,
				 (double)input.wind_north_m_s,
				 (double)input.wind_east_m_s);

			PX4_INFO("[bemt] kin V0=%.2f J=%.3f Jn=%.3f Jp=%.3f adisk=%.3f",
				 (double)output.kinematic_v_inf_m_s[0],
				 (double)output.kinematic_j[0],
				 (double)output.kinematic_j_n[0],
				 (double)output.kinematic_j_p[0],
				 (double)output.kinematic_alpha_disk_rad[0]);

			PX4_INFO("[bemt] cor V0=%.2f J=%.3f Jn=%.3f Jp=%.3f adisk=%.3f vi=%.2f",
				 (double)output.corrected_v_inf_m_s[0],
				 (double)output.corrected_j[0],
				 (double)output.corrected_j_n[0],
				 (double)output.corrected_j_p[0],
				 (double)output.corrected_alpha_disk_rad[0],
				 (double)output.induced_axial_velocity_m_s[0]);

			PX4_INFO("[bemt] Re0=%.0f Ct0=%.4f Cq0=%.4f Cp0=%.4f [PROVISIONAL]",
				 (double)output.re_07[0],
				 (double)output.c_t[0],
				 (double)output.c_q[0],
				 (double)output.c_p[0]);

			PX4_INFO("[bemt] dbg r=%.4f phi=%.3f alpha=%.3f Ftip=%.3f",
				 (double)output.section_debug_radius_m[0],
				 (double)output.section_debug_inflow_angle_rad[0],
				 (double)output.section_debug_angle_of_attack_rad[0],
				 (double)output.section_debug_tip_loss_factor[0]);
		}
	}
};

extern "C" __EXPORT int bemt_optimize_main(int argc, char *argv[])
{
	return BEMTOptimize::main(argc, argv);
}