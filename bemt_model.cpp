#include "bemt_model.hpp"

#include <cmath>

namespace bemt
{

static constexpr float kPi = 3.14159265358979323846F;
static constexpr float kMinAirDensity = 0.5F;
static constexpr float kMaxAirDensity = 1.5F;
static constexpr float kMinDynamicViscosityPaS = 1.0e-9F;
static constexpr float kMinRpm = 0.0F;
static constexpr float kMaxRpm = 30000.0F;
static constexpr float kMinVectorNorm = 1.0e-6F;
static constexpr float kMinRevolutionsPerSecond = 1.0e-6F;
static constexpr uint8_t kDebugSectionIndex = 26;

// Engineering-assumption airfoil model until validated section lookup data are wired in.
// This is acceptable for structural BEMT testing, but not for final quantitative claims.
static constexpr float kLiftCurveSlopePerRad = 2.0F * kPi;
static constexpr float kMaxLiftCoefficient = 1.2F;
static constexpr float kBaseDragCoefficient = 0.012F;
static constexpr float kInducedDragFactor = 0.02F;

static bool is_input_valid(const Input &input)
{
	bool valid = true;

	// Aligned with calculate_induced_axial_velocity and calculate_section_forces,
	// which both reject at <=kMinAirDensity; use the same inclusive lower bound here.
	if ((input.air_density_kg_m3 <= kMinAirDensity) ||
	    (input.air_density_kg_m3 > kMaxAirDensity)) {
		valid = false;
	}

	if (input.dynamic_viscosity_pa_s <= kMinDynamicViscosityPaS) {
		valid = false;
	}

	for (int i = 0; i < kMotorCount; i++) {
		if ((input.motor_rpm[i] < kMinRpm) ||
		    (input.motor_rpm[i] > kMaxRpm)) {
			valid = false;
		}
	}

	// Quaternion validity: all components must be finite and the quaternion must be
	// approximately unit-norm. Tolerance of 0.01 on |‖q‖ − 1| is a conservative
	// sanity bound; EKF output is typically within 1e-5 of unity after normalization.
	// A corrupted quaternion would silently poison the body-frame rotation.
	static constexpr float kQuatNormTolerance = 0.01F;
	bool q_finite = true;
	float q_norm_sq = 0.0F;

	for (int i = 0; i < 4; ++i) {
		if (!std::isfinite(input.attitude_q[i])) {
			q_finite = false;
			break;
		}

		q_norm_sq += input.attitude_q[i] * input.attitude_q[i];
	}

	if (!q_finite || (std::fabs(std::sqrt(q_norm_sq) - 1.0F) > kQuatNormTolerance)) {
		valid = false;
	}

	return valid;
}

static float rpm_to_rad_s(const float rpm)
{
	return rpm * 2.0F * kPi / 60.0F;
}

static float rpm_to_rev_s(const float rpm)
{
	return rpm / 60.0F;
}

static float clamp(const float value, const float min_value, const float max_value)
{
	if (value < min_value) {
		return min_value;
	}

	if (value > max_value) {
		return max_value;
	}

	return value;
}

static Vector3 make_vector3(const float x_m_s, const float y_m_s, const float z_m_s)
{
	return Vector3{x_m_s, y_m_s, z_m_s};
}

static Vector3 add(const Vector3 &a, const Vector3 &b)
{
	return make_vector3(a.x_m_s + b.x_m_s, a.y_m_s + b.y_m_s, a.z_m_s + b.z_m_s);
}

static Vector3 subtract(const Vector3 &a, const Vector3 &b)
{
	return make_vector3(a.x_m_s - b.x_m_s, a.y_m_s - b.y_m_s, a.z_m_s - b.z_m_s);
}

static Vector3 scale(const Vector3 &v, const float scalar)
{
	return make_vector3(v.x_m_s * scalar, v.y_m_s * scalar, v.z_m_s * scalar);
}

static float dot(const Vector3 &a, const Vector3 &b)
{
	return (a.x_m_s * b.x_m_s) + (a.y_m_s * b.y_m_s) + (a.z_m_s * b.z_m_s);
}

static Vector3 cross(const Vector3 &a, const Vector3 &b)
{
	return make_vector3((a.y_m_s * b.z_m_s) - (a.z_m_s * b.y_m_s),
			    (a.z_m_s * b.x_m_s) - (a.x_m_s * b.z_m_s),
			    (a.x_m_s * b.y_m_s) - (a.y_m_s * b.x_m_s));
}

static float norm(const Vector3 &v)
{
	return std::sqrt(dot(v, v));
}

static bool normalize(const Vector3 &v, Vector3 &unit_v)
{
	const float magnitude = norm(v);

	if (magnitude < kMinVectorNorm) {
		return false;
	}

	unit_v = scale(v, 1.0F / magnitude);
	return true;
}

static float section_width_from_index(const uint8_t section_index)
{
	if (section_index == 0) {
		return kSectionRadiusM[1] - kSectionRadiusM[0];
	}

	if (section_index == (kBladeSections - 1)) {
		return kSectionRadiusM[kBladeSections - 1] - kSectionRadiusM[kBladeSections - 2];
	}

	return 0.5F * (kSectionRadiusM[section_index + 1] - kSectionRadiusM[section_index - 1]);
}

static Vector3 rotate_ned_to_body(const float attitude_q[4], const Vector3 &ned_vector_m_s)
{
	const float q0 = attitude_q[0];
	const float q1 = attitude_q[1];
	const float q2 = attitude_q[2];
	const float q3 = attitude_q[3];

	const float r11 = 1.0F - 2.0F * ((q2 * q2) + (q3 * q3));
	const float r12 = 2.0F * ((q1 * q2) - (q0 * q3));
	const float r13 = 2.0F * ((q1 * q3) + (q0 * q2));
	const float r21 = 2.0F * ((q1 * q2) + (q0 * q3));
	const float r22 = 1.0F - 2.0F * ((q1 * q1) + (q3 * q3));
	const float r23 = 2.0F * ((q2 * q3) - (q0 * q1));
	const float r31 = 2.0F * ((q1 * q3) - (q0 * q2));
	const float r32 = 2.0F * ((q2 * q3) + (q0 * q1));
	const float r33 = 1.0F - 2.0F * ((q1 * q1) + (q2 * q2));

	return make_vector3((r11 * ned_vector_m_s.x_m_s) + (r21 * ned_vector_m_s.y_m_s) +
				    (r31 * ned_vector_m_s.z_m_s),
			    (r12 * ned_vector_m_s.x_m_s) + (r22 * ned_vector_m_s.y_m_s) +
				    (r32 * ned_vector_m_s.z_m_s),
			    (r13 * ned_vector_m_s.x_m_s) + (r23 * ned_vector_m_s.y_m_s) +
				    (r33 * ned_vector_m_s.z_m_s));
}

bool get_rotor_geometry(uint8_t motor_index, RotorGeometry &rotor_geometry)
{
	if (motor_index >= kMotorCount) {
		return false;
	}

	rotor_geometry.position_x_m = kRotorPos[motor_index][0];
	rotor_geometry.position_y_m = kRotorPos[motor_index][1];
	rotor_geometry.position_z_m = kRotorPos[motor_index][2];
	rotor_geometry.axis_x = kRotorAxis[0];
	rotor_geometry.axis_y = kRotorAxis[1];
	rotor_geometry.axis_z = kRotorAxis[2];
	return true;
}

bool get_blade_section(uint8_t section_index, BladeSection &section)
{
	if (section_index >= kBladeSections) {
		return false;
	}

	section.radius_m = kSectionRadiusM[section_index];
	section.radius_ratio = kSectionRadiusM[section_index] / kPropRadiusM;
	section.chord_m = kSectionChordM[section_index];
	section.pitch_rad = kSectionPitchRad[section_index];
	return true;
}

bool interpolate_blade_section(float radius_m, BladeSection &section)
{
	if (radius_m < kSectionRadiusM[0] || radius_m > kSectionRadiusM[kBladeSections - 1]) {
		return false;
	}

	if (radius_m <= kSectionRadiusM[0]) {
		return get_blade_section(0, section);
	}

	for (uint8_t i = 0; i < (kBladeSections - 1); ++i) {
		const float r0 = kSectionRadiusM[i];
		const float r1 = kSectionRadiusM[i + 1];

		if (radius_m <= r1) {
			const float t = (radius_m - r0) / (r1 - r0);
			section.radius_m = radius_m;
			section.radius_ratio = radius_m / kPropRadiusM;
			section.chord_m = kSectionChordM[i] + t * (kSectionChordM[i + 1] - kSectionChordM[i]);
			section.pitch_rad = kSectionPitchRad[i] + t * (kSectionPitchRad[i + 1] - kSectionPitchRad[i]);
			return true;
		}
	}

	return get_blade_section(kBladeSections - 1, section);
}

bool calculate_air_relative_velocity_body(const Input &input, Vector3 &air_velocity_body_m_s)
{
	const Vector3 vehicle_velocity_ned_m_s = make_vector3(input.velocity_north_m_s,
							      input.velocity_east_m_s,
							      input.velocity_down_m_s);
	const Vector3 wind_velocity_ned_m_s = make_vector3(input.wind_north_m_s,
							   input.wind_east_m_s,
							   0.0F);
	const Vector3 air_velocity_ned_m_s = subtract(vehicle_velocity_ned_m_s, wind_velocity_ned_m_s);
	air_velocity_body_m_s = rotate_ned_to_body(input.attitude_q, air_velocity_ned_m_s);
	return true;
}

bool calculate_local_air_relative_velocity_body(const Vector3 &air_velocity_body_m_s,
						const BodyRates &body_rates_rad_s,
						const RotorGeometry &rotor_geometry,
						Vector3 &local_air_velocity_body_m_s)
{
	const Vector3 omega_body_rad_s = make_vector3(body_rates_rad_s.roll_rad_s,
						      body_rates_rad_s.pitch_rad_s,
						      body_rates_rad_s.yaw_rad_s);
	const Vector3 rotor_position_m = make_vector3(rotor_geometry.position_x_m,
						      rotor_geometry.position_y_m,
						      rotor_geometry.position_z_m);
	const Vector3 rigid_body_velocity_m_s = cross(omega_body_rad_s, rotor_position_m);
	local_air_velocity_body_m_s = add(air_velocity_body_m_s, rigid_body_velocity_m_s);
	return true;
}

bool calculate_rotor_inflow_vector_body(const Vector3 &local_air_velocity_body_m_s,
					Vector3 &inflow_velocity_body_m_s)
{
	inflow_velocity_body_m_s = scale(local_air_velocity_body_m_s, -1.0F);
	return true;
}

bool decompose_inflow_components(const Vector3 &inflow_velocity_body_m_s,
				 const RotorGeometry &rotor_geometry,
				 InflowComponents &components)
{
	const Vector3 rotor_axis = make_vector3(rotor_geometry.axis_x,
						rotor_geometry.axis_y,
						rotor_geometry.axis_z);
	Vector3 rotor_axis_unit{};

	if (!normalize(rotor_axis, rotor_axis_unit)) {
		return false;
	}

	components.signed_axial_speed_m_s = dot(inflow_velocity_body_m_s, rotor_axis_unit);
	components.axial_velocity_body_m_s = scale(rotor_axis_unit, components.signed_axial_speed_m_s);
	components.inplane_velocity_body_m_s =
		subtract(inflow_velocity_body_m_s, components.axial_velocity_body_m_s);
	components.v_normal_m_s = std::fabs(components.signed_axial_speed_m_s);
	components.v_inplane_m_s = norm(components.inplane_velocity_body_m_s);
	return true;
}

bool calculate_induced_axial_velocity(float air_density_kg_m3,
				      float thrust_n,
				      float v_normal_kinematic_m_s,
				      float &induced_axial_velocity_m_s)
{
	static constexpr float kMinDiskAreaM2 = 1.0e-9F;
	const float rotor_disk_area_m2 = kPi * kPropRadiusM * kPropRadiusM;

	if ((air_density_kg_m3 <= kMinAirDensity) ||
	    (rotor_disk_area_m2 <= kMinDiskAreaM2) ||
	    (thrust_n <= 0.0F)) {
		induced_axial_velocity_m_s = 0.0F;
		return false;
	}

	const float disk_loading_term = thrust_n / (2.0F * air_density_kg_m3 * rotor_disk_area_m2);

	induced_axial_velocity_m_s =
		-0.5F * v_normal_kinematic_m_s +
		std::sqrt((0.25F * v_normal_kinematic_m_s * v_normal_kinematic_m_s) +
			  disk_loading_term);

	if (induced_axial_velocity_m_s < 0.0F) {
		induced_axial_velocity_m_s = 0.0F;
	}

	return true;
}

bool apply_induced_axial_correction(const InflowComponents &kinematic_components,
				    float induced_axial_velocity_m_s,
				    InflowComponents &corrected_components)
{
	corrected_components = kinematic_components;

	const float axial_sign = (kinematic_components.signed_axial_speed_m_s >= 0.0F) ? 1.0F : -1.0F;
	corrected_components.signed_axial_speed_m_s =
		kinematic_components.signed_axial_speed_m_s + (axial_sign * induced_axial_velocity_m_s);
	corrected_components.v_normal_m_s = std::fabs(corrected_components.signed_axial_speed_m_s);

	if (std::fabs(kinematic_components.signed_axial_speed_m_s) > kMinVectorNorm) {
		const float axial_scale =
			corrected_components.signed_axial_speed_m_s / kinematic_components.signed_axial_speed_m_s;
		corrected_components.axial_velocity_body_m_s =
			scale(kinematic_components.axial_velocity_body_m_s, axial_scale);
	}

	// The in-plane component is intentionally left unchanged by this simple
	// induced-flow correction. This is an engineering approximation.
	return true;
}

bool calculate_section_flow_state(const BladeSection &section,
				  const InflowComponents &corrected_components,
				  float motor_rpm,
				  SectionFlowState &flow_state)
{
	const float omega_rad_s = rpm_to_rad_s(motor_rpm);

	flow_state.axial_velocity_m_s = corrected_components.v_normal_m_s;
	flow_state.inplane_velocity_m_s = corrected_components.v_inplane_m_s;
	flow_state.rotational_velocity_m_s = omega_rad_s * section.radius_m;

	// Azimuth-averaged RMS tangential velocity using an axisymmetric engineering approximation.
	flow_state.tangential_velocity_m_s =
		std::sqrt((flow_state.rotational_velocity_m_s * flow_state.rotational_velocity_m_s) +
			  (0.5F * flow_state.inplane_velocity_m_s * flow_state.inplane_velocity_m_s));
	flow_state.resultant_velocity_m_s =
		std::sqrt((flow_state.axial_velocity_m_s * flow_state.axial_velocity_m_s) +
			  (flow_state.tangential_velocity_m_s * flow_state.tangential_velocity_m_s));

	if (flow_state.resultant_velocity_m_s < kMinVectorNorm) {
		flow_state.inflow_angle_rad = 0.0F;
		return false;
	}

	flow_state.inflow_angle_rad =
		std::atan2(flow_state.axial_velocity_m_s, flow_state.tangential_velocity_m_s);
	return true;
}

bool calculate_section_angle_of_attack(const BladeSection &section,
				       const SectionFlowState &flow_state,
				       SectionAerodynamicState &aero_state)
{
	aero_state.angle_of_attack_rad = section.pitch_rad - flow_state.inflow_angle_rad;
	aero_state.lift_coefficient = clamp(kLiftCurveSlopePerRad * aero_state.angle_of_attack_rad,
						-kMaxLiftCoefficient,
						kMaxLiftCoefficient);
	aero_state.drag_coefficient =
		kBaseDragCoefficient +
		(kInducedDragFactor * aero_state.lift_coefficient * aero_state.lift_coefficient);
	return true;
}

bool calculate_prandtl_tip_loss(const BladeSection &section,
				const SectionFlowState &flow_state,
				float &tip_loss_factor)
{
	tip_loss_factor = 1.0F;

	const float sin_phi = std::sin(flow_state.inflow_angle_rad);
	const float radius_gap_m = kPropRadiusM - section.radius_m;

	// Standard propeller theory: Prandtl tip-loss factor for finite blade count.
	if ((sin_phi <= kMinVectorNorm) || (radius_gap_m <= 0.0F) || (section.radius_m <= kHubRadiusM)) {
		return false;
	}

	const float f =
		(static_cast<float>(kBladeCount) / 2.0F) *
		(radius_gap_m / (section.radius_m * sin_phi));
	const float exponent = std::exp(-f);
	const float clamped_exponent = clamp(exponent, 0.0F, 1.0F);
	tip_loss_factor = (2.0F / kPi) * std::acos(clamped_exponent);
	tip_loss_factor = clamp(tip_loss_factor, 0.0F, 1.0F);
	return true;
}

bool calculate_section_forces(float air_density_kg_m3,
			      const BladeSection &section,
			      const SectionFlowState &flow_state,
			      const SectionAerodynamicState &aero_state,
			      float section_width_m,
			      SectionForces &section_forces)
{
	section_forces.thrust_n = 0.0F;
	section_forces.torque_nm = 0.0F;
	section_forces.tip_loss_factor = 1.0F;

	if ((air_density_kg_m3 <= kMinAirDensity) ||
	    (section_width_m <= 0.0F) ||
	    (flow_state.resultant_velocity_m_s < kMinVectorNorm)) {
		return false;
	}

	const float dynamic_pressure_pa =
		0.5F * air_density_kg_m3 * flow_state.resultant_velocity_m_s * flow_state.resultant_velocity_m_s;
	const float lift_n = dynamic_pressure_pa * aero_state.lift_coefficient * section.chord_m * section_width_m;
	const float drag_n = dynamic_pressure_pa * aero_state.drag_coefficient * section.chord_m * section_width_m;

	const float elemental_thrust_per_blade_n =
		(lift_n * std::cos(flow_state.inflow_angle_rad)) -
		(drag_n * std::sin(flow_state.inflow_angle_rad));
	const float elemental_torque_per_blade_nm =
		((lift_n * std::sin(flow_state.inflow_angle_rad)) +
		 (drag_n * std::cos(flow_state.inflow_angle_rad))) * section.radius_m;

	calculate_prandtl_tip_loss(section, flow_state, section_forces.tip_loss_factor);

	section_forces.thrust_n =
		static_cast<float>(kBladeCount) * section_forces.tip_loss_factor * elemental_thrust_per_blade_n;
	section_forces.torque_nm =
		static_cast<float>(kBladeCount) * section_forces.tip_loss_factor * elemental_torque_per_blade_nm;
	return true;
}

bool integrate_rotor_sections(float air_density_kg_m3,
			      const InflowComponents &components,
			      float motor_rpm,
			      float &thrust_n,
			      float &torque_nm,
			      BladeSection &debug_section,
			      SectionFlowState &debug_flow_state,
			      SectionAerodynamicState &debug_aero_state,
			      float &debug_tip_loss_factor)
{
	thrust_n = 0.0F;
	torque_nm = 0.0F;
	debug_tip_loss_factor = 0.0F;

	bool have_valid_section = false;

	for (uint8_t section_index = 0; section_index < kBladeSections; ++section_index) {
		BladeSection section{};
		SectionFlowState flow_state{};
		SectionAerodynamicState aero_state{};
		SectionForces section_forces{};

		if (!get_blade_section(section_index, section)) {
			continue;
		}

		if (!calculate_section_flow_state(section, components, motor_rpm, flow_state)) {
			continue;
		}

		if (!calculate_section_angle_of_attack(section, flow_state, aero_state)) {
			continue;
		}

		if (!calculate_section_forces(air_density_kg_m3,
					      section,
					      flow_state,
					      aero_state,
					      section_width_from_index(section_index),
					      section_forces)) {
			continue;
		}

		thrust_n += section_forces.thrust_n;
		torque_nm += section_forces.torque_nm;
		have_valid_section = true;

		if (section_index == kDebugSectionIndex) {
			debug_section = section;
			debug_flow_state = flow_state;
			debug_aero_state = aero_state;
			debug_tip_loss_factor = section_forces.tip_loss_factor;
		}
	}

	return have_valid_section;
}

bool calculate_rotor_advance_state(const InflowComponents &components,
				   float motor_rpm,
				   RotorAdvanceState &advance_state)
{
	const float revolutions_per_second = rpm_to_rev_s(motor_rpm);

	advance_state.v_inf_m_s = std::sqrt((components.v_normal_m_s * components.v_normal_m_s) +
					    (components.v_inplane_m_s * components.v_inplane_m_s));
	advance_state.alpha_disk_rad = std::atan2(components.v_inplane_m_s,
						  components.v_normal_m_s);

	if (revolutions_per_second < kMinRevolutionsPerSecond) {
		advance_state.j = 0.0F;
		advance_state.j_n = 0.0F;
		advance_state.j_p = 0.0F;
		return false;
	}

	const float denominator = revolutions_per_second * kPropDiameterM;

	if (denominator < kMinVectorNorm) {
		advance_state.j = 0.0F;
		advance_state.j_n = 0.0F;
		advance_state.j_p = 0.0F;
		return false;
	}

	advance_state.j = advance_state.v_inf_m_s / denominator;
	advance_state.j_n = components.v_normal_m_s / denominator;
	advance_state.j_p = components.v_inplane_m_s / denominator;
	return true;
}

bool calculate_reynolds_07(float air_density_kg_m3,
			   float dynamic_viscosity_pa_s,
			   float motor_rpm,
			   float v_inf_m_s,
			   float &re_07)
{
	if (dynamic_viscosity_pa_s <= kMinDynamicViscosityPaS) {
		re_07 = 0.0F;
		return false;
	}

	const float omega_rad_s = rpm_to_rad_s(motor_rpm);
	const float radius_07_m = 0.7F * kPropRadiusM;
	const float tangential_speed_07_m_s = omega_rad_s * radius_07_m;
	const float v_rel_07_m_s =
		std::sqrt((tangential_speed_07_m_s * tangential_speed_07_m_s) +
			  (v_inf_m_s * v_inf_m_s));

	BladeSection section_07{};

	if (!interpolate_blade_section(radius_07_m, section_07)) {
		re_07 = 0.0F;
		return false;
	}

	re_07 = (air_density_kg_m3 * v_rel_07_m_s * section_07.chord_m) / dynamic_viscosity_pa_s;
	return true;
}

bool calculate_rotor_coefficients(float air_density_kg_m3,
				  float motor_rpm,
				  float thrust_n,
				  float torque_nm,
				  float power_w,
				  RotorCoefficients &coefficients)
{
	const float revolutions_per_second = rpm_to_rev_s(motor_rpm);

	coefficients.c_t = 0.0F;
	coefficients.c_q = 0.0F;
	coefficients.c_p = 0.0F;

	if ((air_density_kg_m3 <= kMinAirDensity) ||
	    (revolutions_per_second < kMinRevolutionsPerSecond)) {
		return false;
	}

	const float n2 = revolutions_per_second * revolutions_per_second;
	const float n3 = n2 * revolutions_per_second;
	const float d2 = kPropDiameterM * kPropDiameterM;
	const float d4 = d2 * d2;
	const float d5 = d4 * kPropDiameterM;

	const float thrust_denominator = air_density_kg_m3 * n2 * d4;
	const float torque_denominator = air_density_kg_m3 * n2 * d5;
	const float power_denominator = air_density_kg_m3 * n3 * d5;

	if ((thrust_denominator < kMinVectorNorm) ||
	    (torque_denominator < kMinVectorNorm) ||
	    (power_denominator < kMinVectorNorm)) {
		return false;
	}

	coefficients.c_t = thrust_n / thrust_denominator;
	coefficients.c_q = torque_nm / torque_denominator;
	coefficients.c_p = power_w / power_denominator;
	return true;
}

bool calculate(const Input &input, Output &output)
{
	if (!is_input_valid(input)) {
		return false;
	}

	output.total_thrust_n = 0.0F;
	output.total_power_w = 0.0F;

	const Vector3 air_velocity_body_m_s = [&input]() {
		Vector3 value{};
		calculate_air_relative_velocity_body(input, value);
		return value;
	}();

	const BodyRates body_rates_rad_s{
		input.roll_rate_rad_s,
		input.pitch_rate_rad_s,
		input.yaw_rate_rad_s
	};

	for (int i = 0; i < kMotorCount; i++) {
		output.thrust_n[i] = 0.0F;
		output.torque_nm[i] = 0.0F;
		output.power_w[i] = 0.0F;

		output.kinematic_signed_axial_speed_m_s[i] = 0.0F;
		output.kinematic_v_normal_m_s[i] = 0.0F;
		output.kinematic_v_inplane_m_s[i] = 0.0F;
		output.kinematic_v_inf_m_s[i] = 0.0F;
		output.kinematic_alpha_disk_rad[i] = 0.0F;
		output.kinematic_j[i] = 0.0F;
		output.kinematic_j_n[i] = 0.0F;
		output.kinematic_j_p[i] = 0.0F;

		output.induced_axial_velocity_m_s[i] = 0.0F;

		output.corrected_signed_axial_speed_m_s[i] = 0.0F;
		output.corrected_v_normal_m_s[i] = 0.0F;
		output.corrected_v_inplane_m_s[i] = 0.0F;
		output.corrected_v_inf_m_s[i] = 0.0F;
		output.corrected_alpha_disk_rad[i] = 0.0F;
		output.corrected_j[i] = 0.0F;
		output.corrected_j_n[i] = 0.0F;
		output.corrected_j_p[i] = 0.0F;

		output.re_07[i] = 0.0F;
		output.c_t[i] = 0.0F;
		output.c_q[i] = 0.0F;
		output.c_p[i] = 0.0F;
		output.section_debug_radius_m[i] = 0.0F;
		output.section_debug_axial_velocity_m_s[i] = 0.0F;
		output.section_debug_tangential_velocity_m_s[i] = 0.0F;
		output.section_debug_resultant_velocity_m_s[i] = 0.0F;
		output.section_debug_inflow_angle_rad[i] = 0.0F;
		output.section_debug_angle_of_attack_rad[i] = 0.0F;
		output.section_debug_tip_loss_factor[i] = 0.0F;

		RotorGeometry rotor_geometry{};
		Vector3 local_air_velocity_body_m_s{};
		Vector3 inflow_velocity_body_m_s{};
		InflowComponents kinematic_components{};
		InflowComponents corrected_components{};
		RotorAdvanceState kinematic_advance_state{};
		RotorAdvanceState corrected_advance_state{};
		RotorCoefficients coefficients{};
		BladeSection debug_section{};
		SectionFlowState debug_flow_state{};
		SectionAerodynamicState debug_aero_state{};
		float debug_tip_loss_factor = 0.0F;

		if (!get_rotor_geometry(static_cast<uint8_t>(i), rotor_geometry)) {
			continue;
		}

		if (!calculate_local_air_relative_velocity_body(air_velocity_body_m_s,
								body_rates_rad_s,
								rotor_geometry,
								local_air_velocity_body_m_s)) {
			continue;
		}

		if (!calculate_rotor_inflow_vector_body(local_air_velocity_body_m_s,
							inflow_velocity_body_m_s)) {
			continue;
		}

		if (!decompose_inflow_components(inflow_velocity_body_m_s,
						 rotor_geometry,
						 kinematic_components)) {
			continue;
		}

		calculate_rotor_advance_state(kinematic_components,
					      input.motor_rpm[i],
					      kinematic_advance_state);

		output.kinematic_signed_axial_speed_m_s[i] = kinematic_components.signed_axial_speed_m_s;
		output.kinematic_v_normal_m_s[i] = kinematic_components.v_normal_m_s;
		output.kinematic_v_inplane_m_s[i] = kinematic_components.v_inplane_m_s;
		output.kinematic_v_inf_m_s[i] = kinematic_advance_state.v_inf_m_s;
		output.kinematic_alpha_disk_rad[i] = kinematic_advance_state.alpha_disk_rad;
		output.kinematic_j[i] = kinematic_advance_state.j;
		output.kinematic_j_n[i] = kinematic_advance_state.j_n;
		output.kinematic_j_p[i] = kinematic_advance_state.j_p;

		float first_pass_thrust_n = 0.0F;
		float first_pass_torque_nm = 0.0F;
		BladeSection unused_debug_section{};
		SectionFlowState unused_debug_flow_state{};
		SectionAerodynamicState unused_debug_aero_state{};
		float unused_debug_tip_loss_factor = 0.0F;

		integrate_rotor_sections(input.air_density_kg_m3,
					 kinematic_components,
					 input.motor_rpm[i],
					 first_pass_thrust_n,
					 first_pass_torque_nm,
					 unused_debug_section,
					 unused_debug_flow_state,
					 unused_debug_aero_state,
					 unused_debug_tip_loss_factor);

		// Single-shot momentum-theory induced-velocity correction: one pass only.
		// Engineering approximation — may underpredict converged induced velocity
		// relative to an iterative solution, especially at high disk loading (hover).
		calculate_induced_axial_velocity(input.air_density_kg_m3,
						 first_pass_thrust_n,
						 kinematic_components.v_normal_m_s,
						 output.induced_axial_velocity_m_s[i]);

		apply_induced_axial_correction(kinematic_components,
					       output.induced_axial_velocity_m_s[i],
					       corrected_components);

		float integrated_thrust_n = 0.0F;
		float integrated_torque_nm = 0.0F;

		integrate_rotor_sections(input.air_density_kg_m3,
					 corrected_components,
					 input.motor_rpm[i],
					 integrated_thrust_n,
					 integrated_torque_nm,
					 debug_section,
					 debug_flow_state,
					 debug_aero_state,
					 debug_tip_loss_factor);

		output.thrust_n[i] = integrated_thrust_n;
		output.torque_nm[i] = integrated_torque_nm;
		output.power_w[i] = output.torque_nm[i] * rpm_to_rad_s(input.motor_rpm[i]);
		output.total_thrust_n += output.thrust_n[i];
		output.total_power_w += output.power_w[i];

		calculate_rotor_advance_state(corrected_components,
					      input.motor_rpm[i],
					      corrected_advance_state);

		calculate_reynolds_07(input.air_density_kg_m3,
				      input.dynamic_viscosity_pa_s,
				      input.motor_rpm[i],
				      corrected_advance_state.v_inf_m_s,
				      output.re_07[i]);

		output.corrected_signed_axial_speed_m_s[i] = corrected_components.signed_axial_speed_m_s;
		output.corrected_v_normal_m_s[i] = corrected_components.v_normal_m_s;
		output.corrected_v_inplane_m_s[i] = corrected_components.v_inplane_m_s;
		output.corrected_v_inf_m_s[i] = corrected_advance_state.v_inf_m_s;
		output.corrected_alpha_disk_rad[i] = corrected_advance_state.alpha_disk_rad;
		output.corrected_j[i] = corrected_advance_state.j;
		output.corrected_j_n[i] = corrected_advance_state.j_n;
		output.corrected_j_p[i] = corrected_advance_state.j_p;

		calculate_rotor_coefficients(input.air_density_kg_m3,
					     input.motor_rpm[i],
					     output.thrust_n[i],
					     output.torque_nm[i],
					     output.power_w[i],
					     coefficients);

		output.c_t[i] = coefficients.c_t;
		output.c_q[i] = coefficients.c_q;
		output.c_p[i] = coefficients.c_p;

		output.section_debug_radius_m[i] = debug_section.radius_m;
		output.section_debug_axial_velocity_m_s[i] = debug_flow_state.axial_velocity_m_s;
		output.section_debug_tangential_velocity_m_s[i] = debug_flow_state.tangential_velocity_m_s;
		output.section_debug_resultant_velocity_m_s[i] = debug_flow_state.resultant_velocity_m_s;
		output.section_debug_inflow_angle_rad[i] = debug_flow_state.inflow_angle_rad;
		output.section_debug_angle_of_attack_rad[i] = debug_aero_state.angle_of_attack_rad;
		output.section_debug_tip_loss_factor[i] = debug_tip_loss_factor;
	}

	return true;
}

} // namespace bemt