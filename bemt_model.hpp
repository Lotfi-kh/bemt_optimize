#pragma once

#include <cstdint>

namespace bemt
{

static constexpr uint8_t kMotorCount = 4;
static constexpr uint8_t kBladeCount = 2;
static constexpr uint8_t kBladeSections = 39;

static constexpr float kPropDiameterM = 0.254F;
static constexpr float kPropRadiusM = 0.127F;
static constexpr float kHubRadiusM = 0.02032F;

static constexpr float kArmLengthM = 0.175F;

// Rotor z-offset in PX4 FRD body frame (z-down): rotors are 0.06 m above base_link.
// Source: x500_base Gazebo SDF (Tools/simulation/gz/models/x500_base/model.sdf),
// all four rotor links at z = +0.06 m in SDF FLU frame (z-up), which maps to -0.06 m
// in PX4 FRD (z-down). This value is used only for the omega-cross-r rigid-body velocity
// correction; its primary effect is on in-plane velocity, not axial thrust.
static constexpr float kRotorPosZBodyM = -0.06F;

static constexpr float kRotorPos[kMotorCount][3] = {
	{ kArmLengthM * 0.7071F,  kArmLengthM * 0.7071F, kRotorPosZBodyM},
	{-kArmLengthM * 0.7071F, -kArmLengthM * 0.7071F, kRotorPosZBodyM},
	{ kArmLengthM * 0.7071F, -kArmLengthM * 0.7071F, kRotorPosZBodyM},
	{-kArmLengthM * 0.7071F,  kArmLengthM * 0.7071F, kRotorPosZBodyM},
};

static constexpr float kRotorAxis[3] = {0.0F, 0.0F, -1.0F};

static constexpr float kSectionRadiusM[kBladeSections] = {
	0.02032F, 0.02184F, 0.02337F, 0.02489F, 0.02642F, 0.02794F, 0.02946F, 0.03112F,
	0.03401F, 0.03703F, 0.04006F, 0.04308F, 0.04610F, 0.04912F, 0.05215F, 0.05517F,
	0.05819F, 0.06121F, 0.06424F, 0.06726F, 0.07028F, 0.07330F, 0.07633F, 0.07935F,
	0.08237F, 0.08539F, 0.08842F, 0.09144F, 0.09446F, 0.09749F, 0.10051F, 0.10353F,
	0.10655F, 0.10958F, 0.11260F, 0.11562F, 0.11864F, 0.12165F, 0.12448F
};

static constexpr float kSectionChordM[kBladeSections] = {
	0.01867F, 0.01978F, 0.02081F, 0.02178F, 0.02268F, 0.02352F, 0.02428F, 0.02504F,
	0.02617F, 0.02710F, 0.02777F, 0.02818F, 0.02834F, 0.02825F, 0.02801F, 0.02771F,
	0.02736F, 0.02694F, 0.02648F, 0.02596F, 0.02540F, 0.02478F, 0.02412F, 0.02341F,
	0.02266F, 0.02188F, 0.02105F, 0.02018F, 0.01928F, 0.01835F, 0.01739F, 0.01639F,
	0.01537F, 0.01432F, 0.01325F, 0.01216F, 0.01104F, 0.00973F, 0.00713F
};

static constexpr float kSectionPitchRad[kBladeSections] = {
	0.59803F, 0.60902F, 0.61077F, 0.60508F, 0.59306F, 0.57549F, 0.55313F, 0.52902F,
	0.49116F, 0.45661F, 0.42632F, 0.39959F, 0.37586F, 0.35468F, 0.33568F, 0.31854F,
	0.30302F, 0.28890F, 0.27600F, 0.26418F, 0.25331F, 0.24329F, 0.23401F, 0.22540F,
	0.21736F, 0.20989F, 0.20292F, 0.19638F, 0.19025F, 0.18449F, 0.17906F, 0.17394F,
	0.16910F, 0.16452F, 0.16018F, 0.15606F, 0.15214F, 0.14844F, 0.14511F
};

struct Vector3
{
	float x_m_s;
	float y_m_s;
	float z_m_s;
};

struct BodyRates
{
	float roll_rad_s;
	float pitch_rad_s;
	float yaw_rad_s;
};

struct RotorGeometry
{
	float position_x_m;
	float position_y_m;
	float position_z_m;

	float axis_x;
	float axis_y;
	float axis_z;
};

struct BladeSection
{
	float radius_m;
	float radius_ratio;
	float chord_m;
	float pitch_rad;
};

struct InflowComponents
{
	Vector3 axial_velocity_body_m_s;
	Vector3 inplane_velocity_body_m_s;

	float signed_axial_speed_m_s;
	float v_normal_m_s;
	float v_inplane_m_s;
};

struct SectionFlowState
{
	float axial_velocity_m_s;
	float inplane_velocity_m_s;
	float rotational_velocity_m_s;
	float tangential_velocity_m_s;
	float resultant_velocity_m_s;
	float inflow_angle_rad;
};

struct SectionAerodynamicState
{
	float angle_of_attack_rad;
	float lift_coefficient;
	float drag_coefficient;
};

struct SectionForces
{
	float thrust_n;
	float torque_nm;
	float tip_loss_factor;
};

struct RotorAdvanceState
{
	float v_inf_m_s;
	float j;
	float j_n;
	float j_p;
	float alpha_disk_rad;
	float re_07;
};

struct RotorCoefficients
{
	float c_t;
	float c_q;
	float c_p;
};

struct Input
{
	float air_density_kg_m3;
	float dynamic_viscosity_pa_s;

	float altitude_m;
	float velocity_north_m_s;
	float velocity_east_m_s;
	float velocity_down_m_s;
	float wind_north_m_s;
	float wind_east_m_s;

	float attitude_q[4];
	float roll_rate_rad_s;
	float pitch_rate_rad_s;
	float yaw_rate_rad_s;

	float battery_voltage_v;
	float battery_current_a;
	float battery_remaining;

	float motor_command[kMotorCount];
	float motor_voltage_v[kMotorCount];
	float motor_current_a[kMotorCount];
	float motor_rpm[kMotorCount];
};

struct Output
{
	float thrust_n[kMotorCount];
	float torque_nm[kMotorCount];
	float power_w[kMotorCount];

	float kinematic_signed_axial_speed_m_s[kMotorCount];
	float kinematic_v_normal_m_s[kMotorCount];
	float kinematic_v_inplane_m_s[kMotorCount];
	float kinematic_v_inf_m_s[kMotorCount];
	float kinematic_alpha_disk_rad[kMotorCount];
	float kinematic_j[kMotorCount];
	float kinematic_j_n[kMotorCount];
	float kinematic_j_p[kMotorCount];

	float induced_axial_velocity_m_s[kMotorCount];

	float corrected_signed_axial_speed_m_s[kMotorCount];
	float corrected_v_normal_m_s[kMotorCount];
	float corrected_v_inplane_m_s[kMotorCount];
	float corrected_v_inf_m_s[kMotorCount];
	float corrected_alpha_disk_rad[kMotorCount];
	float corrected_j[kMotorCount];
	float corrected_j_n[kMotorCount];
	float corrected_j_p[kMotorCount];

	float re_07[kMotorCount];
	float c_t[kMotorCount];
	float c_q[kMotorCount];
	float c_p[kMotorCount];

	float section_debug_radius_m[kMotorCount];
	float section_debug_axial_velocity_m_s[kMotorCount];
	float section_debug_tangential_velocity_m_s[kMotorCount];
	float section_debug_resultant_velocity_m_s[kMotorCount];
	float section_debug_inflow_angle_rad[kMotorCount];
	float section_debug_angle_of_attack_rad[kMotorCount];
	float section_debug_tip_loss_factor[kMotorCount];

	float total_thrust_n;
	float total_power_w;
};

bool calculate(const Input &input, Output &output);

bool get_rotor_geometry(uint8_t motor_index, RotorGeometry &rotor_geometry);
bool get_blade_section(uint8_t section_index, BladeSection &section);
bool interpolate_blade_section(float radius_m, BladeSection &section);

bool calculate_air_relative_velocity_body(const Input &input, Vector3 &air_velocity_body_m_s);

bool calculate_local_air_relative_velocity_body(const Vector3 &air_velocity_body_m_s,
						const BodyRates &body_rates_rad_s,
						const RotorGeometry &rotor_geometry,
						Vector3 &local_air_velocity_body_m_s);

bool calculate_rotor_inflow_vector_body(const Vector3 &local_air_velocity_body_m_s,
					Vector3 &inflow_velocity_body_m_s);

bool decompose_inflow_components(const Vector3 &inflow_velocity_body_m_s,
				 const RotorGeometry &rotor_geometry,
				 InflowComponents &components);

bool calculate_induced_axial_velocity(float air_density_kg_m3,
				      float thrust_n,
				      float v_normal_kinematic_m_s,
				      float &induced_axial_velocity_m_s);

bool apply_induced_axial_correction(const InflowComponents &kinematic_components,
				    float induced_axial_velocity_m_s,
				    InflowComponents &corrected_components);

bool calculate_section_flow_state(const BladeSection &section,
				  const InflowComponents &corrected_components,
				  float motor_rpm,
				  SectionFlowState &flow_state);

bool calculate_section_angle_of_attack(const BladeSection &section,
				       const SectionFlowState &flow_state,
				       SectionAerodynamicState &aero_state);

bool calculate_prandtl_tip_loss(const BladeSection &section,
				const SectionFlowState &flow_state,
				float &tip_loss_factor);

bool calculate_section_forces(float air_density_kg_m3,
			      const BladeSection &section,
			      const SectionFlowState &flow_state,
			      const SectionAerodynamicState &aero_state,
			      float section_width_m,
			      SectionForces &section_forces);

bool integrate_rotor_sections(float air_density_kg_m3,
			      const InflowComponents &components,
			      float motor_rpm,
			      float &thrust_n,
			      float &torque_nm,
			      BladeSection &debug_section,
			      SectionFlowState &debug_flow_state,
			      SectionAerodynamicState &debug_aero_state,
			      float &debug_tip_loss_factor);

bool calculate_rotor_advance_state(const InflowComponents &components,
				   float motor_rpm,
				   RotorAdvanceState &advance_state);

bool calculate_reynolds_07(float air_density_kg_m3,
			   float dynamic_viscosity_pa_s,
			   float motor_rpm,
			   float v_inf_m_s,
			   float &re_07);

bool calculate_rotor_coefficients(float air_density_kg_m3,
				  float motor_rpm,
				  float thrust_n,
				  float torque_nm,
				  float power_w,
				  RotorCoefficients &coefficients);

} // namespace bemt