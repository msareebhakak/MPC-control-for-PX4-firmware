#include "mc_att_control.hpp"

#include <drivers/drv_hrt.h>
#include <mathlib/math/Limits.hpp>
#include <mathlib/math/Functions.hpp>

//#include "test_macros.hpp"

#include <matrix/math.hpp>


using namespace matrix;
//Constructor
MulticopterAttitudeControl::MulticopterAttitudeControl(bool vtol) :
	ModuleParams(nullptr),
	WorkItem(MODULE_NAME, px4::wq_configurations::att_pos_ctrl),
	_vehicle_attitude_setpoint_pub(vtol ? ORB_ID(mc_virtual_attitude_setpoint) : ORB_ID(vehicle_attitude_setpoint)),
	_loop_perf(perf_alloc(PC_ELAPSED, MODULE_NAME": cycle"))
{
	if (vtol) {
		int32_t vt_type = -1;

		if (param_get(param_find("VT_TYPE"), &vt_type) == PX4_OK) {
			_is_tailsitter = (static_cast<vtol_type>(vt_type) == vtol_type::TAILSITTER);
		}
	}

	_vehicle_status.vehicle_type = vehicle_status_s::VEHICLE_TYPE_ROTARY_WING;

	/* initialize quaternions in messages to be valid */
	_v_att.q[0] = 1.f;
	_v_att_sp.q_d[0] = 1.f;

	parameters_updated();
//

}
///Destructor
MulticopterAttitudeControl::~MulticopterAttitudeControl()
{
	perf_free(_loop_perf);
}

bool
MulticopterAttitudeControl::init()
{
	if (!_vehicle_attitude_sub.registerCallback()) {
		PX4_ERR("vehicle_attitude callback registration failed!");
		return false;
	}

	return true;
}

void
MulticopterAttitudeControl::parameters_updated()
{
	// Store some of the parameters in a more convenient way & precompute often-used values
	_attitude_control.setProportionalGain(Vector3f(_param_mc_roll_p.get(), _param_mc_pitch_p.get(), _param_mc_yaw_p.get()));

	// angular rate limits
	using math::radians;
	_attitude_control.setRateLimit(Vector3f(radians(_param_mc_rollrate_max.get()), radians(_param_mc_pitchrate_max.get()),
						radians(_param_mc_yawrate_max.get())));

	_man_tilt_max = math::radians(_param_mpc_man_tilt_max.get());
}

float
MulticopterAttitudeControl::throttle_curve(float throttle_stick_input)
{
	float throttle_min = _vehicle_land_detected.landed ? 0.0f : _param_mpc_manthr_min.get();

	// throttle_stick_input is in range [0, 1]
	switch (_param_mpc_thr_curve.get()) {
	case 1: // no rescaling to hover throttle
		return throttle_min + throttle_stick_input * (_param_mpc_thr_max.get() - throttle_min);

	default: // 0 or other: rescale to hover throttle at 0.5 stick
		if (throttle_stick_input < 0.5f) {
			return (_param_mpc_thr_hover.get() - throttle_min) / 0.5f * throttle_stick_input +
			       throttle_min;

		} else {
			return (_param_mpc_thr_max.get() - _param_mpc_thr_hover.get()) / 0.5f * (throttle_stick_input - 1.0f) +
			       _param_mpc_thr_max.get();
		}
	}
}

void
MulticopterAttitudeControl::generate_attitude_setpoint(float dt, bool reset_yaw_sp)
{
	vehicle_attitude_setpoint_s attitude_setpoint{};
	const float yaw = Eulerf(Quatf(_v_att.q)).psi();

	/* reset yaw setpoint to current position if needed */
	if (reset_yaw_sp) {
		_man_yaw_sp = yaw;

	} else if (_manual_control_sp.z > 0.05f || _param_mc_airmode.get() == (int32_t)Mixer::Airmode::roll_pitch_yaw) {

		const float yaw_rate = math::radians(_param_mpc_man_y_max.get());
		attitude_setpoint.yaw_sp_move_rate = _manual_control_sp.r * yaw_rate;
		_man_yaw_sp = wrap_pi(_man_yaw_sp + attitude_setpoint.yaw_sp_move_rate * dt);
	}

	/*
	 * Input mapping for roll & pitch setpoints
	 * ----------------------------------------
	 * We control the following 2 angles:
	 * - tilt angle, given by sqrt(x*x + y*y)
	 * - the direction of the maximum tilt in the XY-plane, which also defines the direction of the motion
	 *
	 * This allows a simple limitation of the tilt angle, the vehicle flies towards the direction that the stick
	 * points to, and changes of the stick input are linear.
	 */
	const float x = _manual_control_sp.x * _man_tilt_max;
	const float y = _manual_control_sp.y * _man_tilt_max;

	// we want to fly towards the direction of (x, y), so we use a perpendicular axis angle vector in the XY-plane
	Vector2f v = Vector2f(y, -x);
	float v_norm = v.norm(); // the norm of v defines the tilt angle

	if (v_norm > _man_tilt_max) { // limit to the configured maximum tilt angle
		v *= _man_tilt_max / v_norm;
	}

	Quatf q_sp_rpy = AxisAnglef(v(0), v(1), 0.f);
	Eulerf euler_sp = q_sp_rpy;
	attitude_setpoint.roll_body = euler_sp(0);
	attitude_setpoint.pitch_body = euler_sp(1);
	// The axis angle can change the yaw as well (noticeable at higher tilt angles).
	// This is the formula by how much the yaw changes:
	//   let a := tilt angle, b := atan(y/x) (direction of maximum tilt)
	//   yaw = atan(-2 * sin(b) * cos(b) * sin^2(a/2) / (1 - 2 * cos^2(b) * sin^2(a/2))).
	attitude_setpoint.yaw_body = _man_yaw_sp + euler_sp(2);

	/* modify roll/pitch only if we're a VTOL */
	if (_vehicle_status.is_vtol) {
		// Construct attitude setpoint rotation matrix. Modify the setpoints for roll
		// and pitch such that they reflect the user's intention even if a large yaw error
		// (yaw_sp - yaw) is present. In the presence of a yaw error constructing a rotation matrix
		// from the pure euler angle setpoints will lead to unexpected attitude behaviour from
		// the user's view as the euler angle sequence uses the  yaw setpoint and not the current
		// heading of the vehicle.
		// However there's also a coupling effect that causes oscillations for fast roll/pitch changes
		// at higher tilt angles, so we want to avoid using this on multicopters.
		// The effect of that can be seen with:
		// - roll/pitch into one direction, keep it fixed (at high angle)
		// - apply a fast yaw rotation
		// - look at the roll and pitch angles: they should stay pretty much the same as when not yawing

		// calculate our current yaw error
		float yaw_error = wrap_pi(attitude_setpoint.yaw_body - yaw);

		// compute the vector obtained by rotating a z unit vector by the rotation
		// given by the roll and pitch commands of the user
		Vector3f zB = {0.0f, 0.0f, 1.0f};
		Dcmf R_sp_roll_pitch = Eulerf(attitude_setpoint.roll_body, attitude_setpoint.pitch_body, 0.0f);
		Vector3f z_roll_pitch_sp = R_sp_roll_pitch * zB;

		// transform the vector into a new frame which is rotated around the z axis
		// by the current yaw error. this vector defines the desired tilt when we look
		// into the direction of the desired heading
		Dcmf R_yaw_correction = Eulerf(0.0f, 0.0f, -yaw_error);
		z_roll_pitch_sp = R_yaw_correction * z_roll_pitch_sp;

		// use the formula z_roll_pitch_sp = R_tilt * [0;0;1]
		// R_tilt is computed from_euler; only true if cos(roll) not equal zero
		// -> valid if roll is not +-pi/2;
		attitude_setpoint.roll_body = -asinf(z_roll_pitch_sp(1));
		attitude_setpoint.pitch_body = atan2f(z_roll_pitch_sp(0), z_roll_pitch_sp(2));
	}

	/* copy quaternion setpoint to attitude setpoint topic */
	Quatf q_sp = Eulerf(attitude_setpoint.roll_body, attitude_setpoint.pitch_body, attitude_setpoint.yaw_body);
	q_sp.copyTo(attitude_setpoint.q_d);
	attitude_setpoint.q_d_valid = true;

	attitude_setpoint.thrust_body[2] = -throttle_curve(_manual_control_sp.z);
	attitude_setpoint.timestamp = hrt_absolute_time();

	_vehicle_attitude_setpoint_pub.publish(attitude_setpoint);
}

/**
 * Attitude controller.
 * Input: 'vehicle_attitude_setpoint' topics (depending on mode)
 * Output: '_rates_sp' vector
 */
void
MulticopterAttitudeControl::control_attitude()
{
	_v_att_sp_sub.update(&_v_att_sp);
	_rates_sp = _attitude_control.update(Quatf(_v_att.q), Quatf(_v_att_sp.q_d), _v_att_sp.yaw_sp_move_rate);
}

void
MulticopterAttitudeControl::publish_rates_setpoint()
{
	vehicle_rates_setpoint_s v_rates_sp{};

	v_rates_sp.roll = _rates_sp(0);
	v_rates_sp.pitch = _rates_sp(1);
	v_rates_sp.yaw = _rates_sp(2);
	v_rates_sp.thrust_body[0] = _v_att_sp.thrust_body[0];
	v_rates_sp.thrust_body[1] = _v_att_sp.thrust_body[1];
	v_rates_sp.thrust_body[2] = _v_att_sp.thrust_body[2];
	v_rates_sp.timestamp = hrt_absolute_time();

	_v_rates_sp_pub.publish(v_rates_sp);
}

void
MulticopterAttitudeControl::Run()
{
	if (should_exit()) {
		_vehicle_attitude_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}

	perf_begin(_loop_perf);

	// Check if parameters have changed
	if (_params_sub.updated()) {
		// clear update
		parameter_update_s param_update;
		_params_sub.copy(&param_update);

		updateParams();
		parameters_updated();
	}

	// run controller on attitude updates
	const uint8_t prev_quat_reset_counter = _v_att.quat_reset_counter;

	if (_vehicle_attitude_sub.update(&_v_att)) {

		// Check for a heading reset
		if (prev_quat_reset_counter != _v_att.quat_reset_counter) {
			// we only extract the heading change from the delta quaternion
			_man_yaw_sp += Eulerf(Quatf(_v_att.delta_q_reset)).psi();
		}

		const hrt_abstime now = hrt_absolute_time();

		// Guard against too small (< 0.2ms) and too large (> 20ms) dt's.
		const float dt = math::constrain(((now - _last_run) / 1e6f), 0.0002f, 0.02f);
		_last_run = now;

		/* check for updates in other topics */
		_manual_control_sp_sub.update(&_manual_control_sp);
		_v_control_mode_sub.update(&_v_control_mode);
		_vehicle_land_detected_sub.update(&_vehicle_land_detected);
		_vehicle_status_sub.update(&_vehicle_status);

		/* Check if we are in rattitude mode and the pilot is above the threshold on pitch
		* or roll (yaw can rotate 360 in normal att control). If both are true don't
		* even bother running the attitude controllers */
		if (_v_control_mode.flag_control_rattitude_enabled) {
			_v_control_mode.flag_control_attitude_enabled =
				fabsf(_manual_control_sp.y) <= _param_mc_ratt_th.get() &&
				fabsf(_manual_control_sp.x) <= _param_mc_ratt_th.get();
		}

		bool attitude_setpoint_generated = false;

		const bool is_hovering = _vehicle_status.vehicle_type == vehicle_status_s::VEHICLE_TYPE_ROTARY_WING
					 && !_vehicle_status.in_transition_mode;

		// vehicle is a tailsitter in transition mode
		const bool is_tailsitter_transition = _vehicle_status.in_transition_mode && _is_tailsitter;

		bool run_att_ctrl = _v_control_mode.flag_control_attitude_enabled && (is_hovering || is_tailsitter_transition);

		if (run_att_ctrl) {
			// Generate the attitude setpoint from stick inputs if we are in Manual/Stabilized mode
			if (_v_control_mode.flag_control_manual_enabled &&
			    !_v_control_mode.flag_control_altitude_enabled &&
			    !_v_control_mode.flag_control_velocity_enabled &&
			    !_v_control_mode.flag_control_position_enabled) {

				generate_attitude_setpoint(dt, _reset_yaw_sp);
				attitude_setpoint_generated = true;
			}

			control_attitude();

			if (_v_control_mode.flag_control_yawrate_override_enabled) {
				/* Yaw rate override enabled, overwrite the yaw setpoint */
				_v_rates_sp_sub.update(&_v_rates_sp);
				const auto yawrate_reference = _v_rates_sp.yaw;
				_rates_sp(2) = yawrate_reference;
			}

			publish_rates_setpoint();
		}

		// reset yaw setpoint during transitions, tailsitter.cpp generates
		// attitude setpoint for the transition
		_reset_yaw_sp = (!attitude_setpoint_generated && !_v_control_mode.flag_control_rattitude_enabled) ||
				_vehicle_land_detected.landed ||
				(_vehicle_status.is_vtol && _vehicle_status.in_transition_mode);

	}

	perf_end(_loop_perf);
}

int MulticopterAttitudeControl::task_spawn(int argc, char *argv[])
{
	bool vtol = false;

	if (argc > 1) {
		if (strcmp(argv[1], "vtol") == 0) {
			vtol = true;
		}
	}

	MulticopterAttitudeControl *instance = new MulticopterAttitudeControl(vtol);

	if (instance) {
		_object.store(instance);
		_task_id = task_id_is_work_queue;

		if (instance->init()) {
			return PX4_OK;
		}

	} else {
		PX4_ERR("alloc failed");
	}

	delete instance;
	_object.store(nullptr);
	_task_id = -1;

	return PX4_ERROR;
}

int MulticopterAttitudeControl::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int MulticopterAttitudeControl::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
This implements the multicopter attitude controller. It takes attitude
setpoints (`vehicle_attitude_setpoint`) as inputs and outputs a rate setpoint.
The controller has a P loop for angular error
Publication documenting the implemented Quaternion Attitude Control:
Nonlinear Quadrocopter Attitude Control (2013)
by Dario Brescianini, Markus Hehn and Raffaello D'Andrea
Institute for Dynamic Systems and Control (IDSC), ETH Zurich
https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/154099/eth-7387-01.pdf
)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("mc_att_control", "controller");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_ARG("vtol", "VTOL mode", true);
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

/// Own function
void MulticopterAttitudeControl::control_attitude_rates(float dt)
{
	Matrix<float,6, 1> QPhild(const Matrix<float, 6, 6> &E,
	Matrix<float,6, 1> &F, const Matrix<float, 24, 6> &CC,
	Matrix<float,24, 1> &d);
	Matrix<float,6, 1> x;
	Matrix<float,3, 1> uu;
  Matrix<float,9, 1> Xf;
	Matrix<float,3, 1> y;
	SquareMatrix<float,6> E;

// Initialisations
x.setZero();
uu.setZero();
y.setZero();
Xf.setZero();

float dataP[135]= {
	0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,
	2.0f,0.0f,0.0f,3.0f,0.0f,0.0f,4.0f,0.0f,0.0f,
	5.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,2.0f,0.0f,0.0f,3.0f,0.0f,
	0.0f,4.0f,0.0f,0.0f,5.0f,0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,2.0f,
	0.0f,0.0f,3.0f,0.0f,0.0f,4.0f,0.0f,0.0f,5.0f,
	1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,
	1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,1.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,1.0f,
	0.0f,0.0f,1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,1.0f};
Matrix<float,15, 9> P(dataP);

// Output prediction matrix H
float dataH[90] = {
	9.3897f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 9.0212f, 0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 7.0922f, 0.0f, 0.0f, 0.0f,
	18.7793f, 0.0f, 0.0f, 9.3897f, 0.0f, 0.0f,
	0.0f, 18.0424f, 0.0f, 0.0f, 9.0212f, 0.0f,
	0.0f, 0.0f, 14.1844f, 0.0f, 0.0f, 7.0922f,
	28.169f, 0.0f, 0.0f, 18.7793f, 0.0f, 0.0f,
	0.0f, 27.0636f, 0.0f, 0.0f, 18.0424f, 0.0f,
	0.0f, 0.0f, 21.2766f, 0.0f, 0.0f, 14.1844f,
	37.5587f, 0.0f, 0.0f, 28.1690f, 0.0f, 0.0f,
	0.0f, 36.0848f, 0.0f, 0.0f, 27.0636f, 0.0f,
	0.0f, 0.0f, 28.3688f, 0.0f, 0.0f, 21.2766f,
	46.9484f, 0.0f, 0.0f, 37.5587f, 0.0f, 0.0f,
	0.0f, 45.1060f, 0.0f, 0.0f, 36.0848f, 0.0f,
	0.0f, 0.0f, 35.4610f, 0.0f, 0.0f, 28.3688f};
Matrix<float,15, 6> H(dataH);

// Input weight
float dataW[36] = {
0.065f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.065f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.085f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.065f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 0.065f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.085f};
Matrix<float,6, 6> W(dataW);

// Transpose of H output matrix
//H_t = H.transpose();
Matrix<float,6, 15> H_trans = H.transpose();
// quadratic programming variable, E
E =  (H_trans*H + W );
E = E * 2;
// Constraint matrix
float dataCC[144] = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, -1.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, -1.0f, 0.0f, 0.0f, 0.0f,
-1.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,
0.0f, -1.0f, 0.0f, 0.0f, -1.0f, 0.0f,
0.0f, 0.0f, -1.0f, 0.0f, 0.0f, -1.0f};
Matrix<float,24, 6> CC(dataCC);
// Constraint vector of input and rate of input change
float datadd[24] ={
0.7548f,
0.7548f,
0.1288f,
0.7548f,
0.7548f,
0.1288f,
0.7548f,
0.7548f,
0.1288f,
0.7548f,
0.7548f,
0.1288f,
1.258f,
1.258f,
0.2147f,
1.258f,
1.258f,
0.2147f,
1.258f,
1.258f,
0.2147f,
1.258f,
1.258f,
0.2147f};
Matrix<float,24, 1> dd(datadd);
// Past rate of input change matrix
float datadupast[72] ={0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f,
-1.0f, 0.0f, 0.0f,
0.0f, -1.0f, 0.0,
0.0f, 0.0f, -1.0f,
-1.0f, 0.0f, 0.0f,
0.0f, -1.0f, 0.0f,
0.0f, 0.0f, -1.0f,
1.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f,
0.0f, 0.0f, 1.0f,
1.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f,
0.0f, 0.0f, 1.0f};

Matrix<float,24, 3> dupast(datadupast);
float dataAd[36] = {1.0f, 0.2f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 1.0f, 0.2f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.2f,
0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
// State space matrices
static Matrix<float, 6, 6> Ad(dataAd);

float dataBd[18] = {0.9390f, 0.0f, 0.0f,
9.3897f, 0.0f, 0.0f,
0.0f, 0.9021f, 0.0f,
0.0f, 9.0212f, 0.0f,
0.0f, 0.0f, 0.7092f,
0.0f, 0.0f, 7.0922f};
static Matrix<float, 6, 3> Bd(dataBd);


float dataCd[18] = {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
static Matrix<float, 3, 6> Cd(dataCd);

// Desired attitude rates/ rates setpoint
float datades[3] = {_rates_sp(0),_rates_sp(1),_rates_sp(2)};
// Reference adjusting matrix
Matrix<float,3, 1> des(datades);

float dataRs[45] ={ 1, 0, 0,
0, 1, 0,
0, 0, 1,
1, 0, 0,
0, 1, 0,
0, 0, 1,
1, 0, 0,
0, 1, 0,
0, 0, 1,
1, 0, 0,
0, 1, 0,
0, 0, 1,
1, 0, 0,
0, 1, 0,
0, 0, 1};
static Matrix<float, 15, 3> Rs(dataRs);
// disturbance

float datadist[3] = {1,1.1,0.9};
Matrix<float,3, 1> dist(datadist);
dist = dist*0.06;
// current rates

float datax_curr[3]={_rates_sp(0),
_rates_sp(1),
_rates_sp(2)};

Matrix<float,3, 1> x_curr(datax_curr);
// Defining previous state vector

float datax_prev[6] = {0,
	_rates_sp(0),
	0,
	_rates_sp(1),
	0,
	_rates_sp(2)};
Matrix<float,6, 1> x_prev(datax_prev);

int i = 0;
do{
Matrix<float, 6, 1> F =  (H_trans)*(Rs*des - P*Xf);
F = F * -2;
//Matrix<float, 24, 1> d = dd + dupast*uu;
///////////////////Function starts here////////////////////////////////

static SquareMatrix<float, 6> E_cholesky = cholesky(E);
//static Matrix<float, 6, 24> CC_transd = CC.transpose();
LeastSquaresSolver<float, 6, 6> qrd = LeastSquaresSolver<float, 6, 6>(E_cholesky);
///error line
///
//Matrix<float, 24, 1> K = (CC*(qrd.solve(F)) + d);

//int k_row = 24;
Matrix<float, 24, 1> lambda;
lambda.setZero();
//float al = 3.0f;
//int km = 0;
/***
do
{
Matrix<float, 24, 1> lambda_p = lambda;
// loop to determine lambda values for respective iterations
int j = 0;
do
{
/// solved upto this
float Tjj = T(j, j);
T(j, j) = 0;
float la = -(T.col(j).dot(lambda) + K(j)) / Tjj;
T(j, j) = Tjj;
if (la < 0.0f) lambda(j) = 0.0f;
else lambda(j) = la;
j++;
} while (j < k_row);
al = (lambda - lambda_p).norm_squared();

if (al < 0.001f) break;
km++;
} while (km < 15);
***/
Matrix<float, 6, 1> DeltaU = -qrd.solve(F); //-(qrd.solve(CC_transd))*lambda;

//////////////////Function ends here//////////////////////////
float DeltaU_1data[6] = {DeltaU(0, 0), DeltaU(1, 0), DeltaU(2, 0),
DeltaU(3, 0), DeltaU(4, 0), DeltaU(5, 0)};
Matrix<float, 2, 3> DeltaU_1(DeltaU_1data);
Matrix<float, 1, 3> newDeltaU_1(DeltaU_1.row(0));
Matrix<float, 3, 1>deltau_tran=newDeltaU_1.transpose();
//Matrix<float, 3, 1>deltau_tran=(DeltaU_1.row(0)).transpose();
uu = uu + deltau_tran;
x = Ad*x_prev + Bd*uu;
y = Cd*x;//+ dist;
Matrix<float, 6, 1> xx = (x-x_prev);
Xf(1,0) = xx(1,0);Xf(2,0) = xx(0,0);Xf(3,0) = xx(3,0);
Xf(4,0) = xx(4,0);Xf(5,0) = xx(0,0);Xf(6,0) = xx(6,0);
Xf(7,0) = y(1,0); Xf(8,0) = y(2,0);Xf(9,0) = y(3,0);
x_prev = x;
i++;
} while(i < 2);
Matrix<float,3,1> _att_control;
_att_control(0,0) = uu(0,0);
_att_control(1,0) = uu(1,0);
_att_control(2,0) = uu(2,0);
float umaxx_data[3] = {1.258f, 1.258f, 0.2147f};
Vector<float,3> umaxx(umaxx_data);

Vector<float,3> uminn = -umaxx;

for (int k = 0; k < 3; k++)
{
_att_control(k,0) = (2 * ((_att_control(k,0) - uminn(k))
/(umaxx(k) - uminn(k)))-1);
}
}
////Solved upto this point
//////////////////////////////////////////



int mc_att_control_main(int argc, char *argv[])
{
	return MulticopterAttitudeControl::main(argc, argv);
}
