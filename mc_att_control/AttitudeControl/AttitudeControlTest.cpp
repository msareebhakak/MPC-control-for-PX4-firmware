/****************************************************************************
 *
 *   Copyright (C) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#include <gtest/gtest.h>
#include <AttitudeControl.hpp>
#include <mathlib/math/Functions.hpp>

using namespace matrix;

TEST(AttitudeControlTest, AllZeroCase)
{
	AttitudeControl attitude_control;
	Vector3f rate_setpoint = attitude_control.update(Quatf(), Quatf(), 0.f);
	EXPECT_EQ(rate_setpoint, Vector3f());
}

class AttitudeControlConvergenceTest : public ::testing::Test
{
public:
	AttitudeControlConvergenceTest()
	{
		_attitude_control.setProportionalGain(Vector3f(.5f, .6f, .3f));
		_attitude_control.setRateLimit(Vector3f(100, 100, 100));
	}

	void checkConvergence()
	{
		int i; // need function scope to check how many steps
		Vector3f rate_setpoint(1000, 1000, 1000);

		for (i = 100; i > 0; i--) {
			// run attitude control to get rate setpoints
			const Vector3f rate_setpoint_new = _attitude_control.update(_quat_state, _quat_goal, 0.f);
			// rotate the simulated state quaternion according to the rate setpoint
			_quat_state = _quat_state * Quatf(AxisAnglef(rate_setpoint_new));

			// expect the error and hence also the output to get smaller with each iteration
			if (rate_setpoint_new.norm() >= rate_setpoint.norm()) {
				break;
			}

			rate_setpoint = rate_setpoint_new;
		}

		EXPECT_EQ(adaptAntipodal(_quat_state), adaptAntipodal(_quat_goal));
		// it shouldn't have taken longer than an iteration timeout to converge
		EXPECT_GT(i, 0);
	}

	Quatf adaptAntipodal(const Quatf q)
	{
		for (int i = 0; i < 4; i++) {
			if (fabs(q(i)) > FLT_EPSILON) {
				return q * math::sign(q(i));
			}
		}

		return q;
	}

	AttitudeControl _attitude_control;
	Quatf _quat_state;
	Quatf _quat_goal;
};

TEST_F(AttitudeControlConvergenceTest, AttitudeControlConvergence)
{
	const int inputs = 8;

	const Quatf QArray[inputs] = {
		Quatf(),
		Quatf(0, 1, 0, 0),
		Quatf(0, 0, 1, 0),
		Quatf(0, 0, 0, 1),
		Quatf(0.698f, 0.024f, -0.681f, -0.220f),
		Quatf(-0.820f, -0.313f, 0.225f, -0.423f),
		Quatf(0.599f, -0.172f, 0.755f, -0.204f),
		Quatf(0.216f, -0.662f, 0.290f, -0.656f)
	};

	for (int i = 0; i < inputs; i++) {
		for (int j = 0; j < inputs; j++) {
			printf("--- Input combination: %d %d\n", i, j);
			_quat_state = QArray[i];
			_quat_goal = QArray[j];
			_quat_state.normalize();
			_quat_goal.normalize();
			checkConvergence();
		}
	}
}
