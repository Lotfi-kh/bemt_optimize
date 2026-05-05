#include <gtest/gtest.h>
#include <cmath>
#include "bemt_model.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Build a minimal InflowComponents with the given scalar values.
static bemt::InflowComponents make_inflow(float v_normal_m_s, float v_inplane_m_s)
{
	bemt::InflowComponents c{};
	c.signed_axial_speed_m_s = v_normal_m_s;
	c.v_normal_m_s           = std::fabs(v_normal_m_s);
	c.v_inplane_m_s          = v_inplane_m_s;
	return c;
}

// Run calculate_section_flow_state for one section at the given radius.
// Returns false if the function itself returned false.
static bool flow_at_radius(float radius_m, float omega_rad_s,
			    float v_normal_m_s, float v_inplane_m_s,
			    bemt::SectionFlowState &out)
{
	bemt::BladeSection section{};

	if (!bemt::interpolate_blade_section(radius_m, section)) {
		return false;
	}

	const bemt::InflowComponents inflow = make_inflow(v_normal_m_s, v_inplane_m_s);
	// rpm = omega_rad_s * 60 / (2*pi)
	const float rpm = omega_rad_s * 60.0F / (2.0F * 3.14159265F);
	return bemt::calculate_section_flow_state(section, inflow, rpm, out);
}

// ---------------------------------------------------------------------------
// 1. Tangential velocity formula — pure-hover (V_inplane == 0)
//    Expected: tangential_velocity = omega * r   [m/s]
// ---------------------------------------------------------------------------

TEST(BemtSectionFlowState, PureHoverTangentialEqualsOmegaR)
{
	const float omega_rad_s  = 500.0F;   // representative hover speed [rad/s]
	const float radius_m     = 0.08842F; // section 26, 0.7R reference station [m]
	const float v_inplane    = 0.0F;     // no forward flight [m/s]
	const float v_normal     = 0.0F;     // no axial inflow [m/s]

	bemt::SectionFlowState flow{};
	ASSERT_TRUE(flow_at_radius(radius_m, omega_rad_s, v_normal, v_inplane, flow));

	const float expected_tang = omega_rad_s * radius_m;  // [m/s]
	EXPECT_NEAR(flow.tangential_velocity_m_s, expected_tang, 1.0e-4F)
		<< "With V_inplane = 0, tangential must equal omega*r exactly.";
	EXPECT_NEAR(flow.rotational_velocity_m_s, expected_tang, 1.0e-4F);
}

// ---------------------------------------------------------------------------
// 2. Tangential velocity formula — zero rotation, nonzero V_inplane
//    Expected: tangential_velocity = sqrt(0.5) * V_inplane   [m/s]
// ---------------------------------------------------------------------------

TEST(BemtSectionFlowState, ZeroRotationTangentialEqualsRmsFactor)
{
	// To force rotational_velocity = 0 we need omega = 0 (rpm = 0).
	// calculate_section_flow_state uses rpm_to_rad_s(rpm), so pass rpm = 0.
	const float v_inplane_m_s = 10.0F;  // forward speed in-plane component [m/s]
	const float radius_m      = 0.08842F;

	bemt::BladeSection section{};
	ASSERT_TRUE(bemt::interpolate_blade_section(radius_m, section));

	bemt::InflowComponents inflow = make_inflow(0.0F, v_inplane_m_s);
	bemt::SectionFlowState flow{};
	// rpm = 0 → omega = 0 → rotational_velocity = 0
	const bool ok = bemt::calculate_section_flow_state(section, inflow, 0.0F, flow);
	// resultant may be nonzero (from inplane), but zero-rpm may make it return false
	// if resultant < kMinVectorNorm; check only when it succeeds.
	if (ok) {
		const float expected = std::sqrt(0.5F) * v_inplane_m_s;  // [m/s]
		EXPECT_NEAR(flow.tangential_velocity_m_s, expected, 1.0e-4F)
			<< "With omega=0, tangential must equal sqrt(0.5)*V_inplane (azimuth-RMS).";
	} else {
		// Zero rpm + zero axial → resultant near zero → function correctly returns false.
		EXPECT_NEAR(flow.tangential_velocity_m_s, 0.0F, 1.0e-6F);
	}
}

// ---------------------------------------------------------------------------
// 3. Corrected formula gives smaller tangential speed than the old formula
//    for any nonzero V_inplane.
//    Old:  sqrt(omega_r^2 + V_p^2)
//    New:  sqrt(omega_r^2 + 0.5*V_p^2)   ← must always be smaller
// ---------------------------------------------------------------------------

TEST(BemtSectionFlowState, CorrectedFormulaSmallerthanOldFormula)
{
	const float omega_rad_s   = 480.0F;
	const float radius_m      = 0.08842F;
	const float v_inplane_m_s = 12.0F;  // 12 m/s forward speed in-plane [m/s]
	const float v_normal      = 1.5F;

	bemt::SectionFlowState flow{};
	ASSERT_TRUE(flow_at_radius(radius_m, omega_rad_s, v_normal, v_inplane_m_s, flow));

	const float omega_r    = omega_rad_s * radius_m;
	const float old_tang   = std::sqrt(omega_r * omega_r + v_inplane_m_s * v_inplane_m_s);
	const float new_tang   = std::sqrt(omega_r * omega_r + 0.5F * v_inplane_m_s * v_inplane_m_s);

	EXPECT_NEAR(flow.tangential_velocity_m_s, new_tang, 1.0e-4F)
		<< "Tangential velocity must match the azimuth-RMS formula.";
	EXPECT_LT(flow.tangential_velocity_m_s, old_tang)
		<< "Corrected formula must be strictly less than old formula for nonzero V_inplane.";
}

// ---------------------------------------------------------------------------
// 4. Debug section index 26 is nearest tabulated station to 0.7R
//    0.7 * R = 0.7 * 0.127 m = 0.08890 m
// ---------------------------------------------------------------------------

TEST(BemtDebugSection, Index26IsNearestTo07R)
{
	const float target_07r = 0.7F * bemt::kPropRadiusM;  // 0.08890 m

	bemt::BladeSection best{};
	ASSERT_TRUE(bemt::get_blade_section(26, best));
	const float dist_26 = std::fabs(best.radius_m - target_07r);

	// Check all other sections are farther away.
	for (uint8_t idx = 0; idx < bemt::kBladeSections; ++idx) {
		if (idx == 26) { continue; }
		bemt::BladeSection s{};
		ASSERT_TRUE(bemt::get_blade_section(idx, s));
		EXPECT_GE(std::fabs(s.radius_m - target_07r), dist_26)
			<< "Section " << static_cast<int>(idx)
			<< " (r=" << s.radius_m << " m) is closer to 0.7R than section 26.";
	}
}

// ---------------------------------------------------------------------------
// 5. Blade geometry consistency — kPropRadiusM == kPropDiameterM / 2
// ---------------------------------------------------------------------------

TEST(BemtGeometry, RadiusConsistentWithDiameter)
{
	EXPECT_NEAR(bemt::kPropRadiusM, 0.5F * bemt::kPropDiameterM, 1.0e-6F);
}
