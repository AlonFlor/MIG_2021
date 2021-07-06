#include "MIGBasicSim.h"


class Person{
public:
	static const int dim = 2;

	float X[dim];
	float V[dim];
	float destination[dim];
	float desired_speed;

	Person()
	{
		
	}

	Person(float* X_in, float* V_in, float* destination_in, float desired_speed_in)
	{
		X[0] = X_in[0];
		X[1] = X_in[1];
		V[0] = V_in[0];
		V[1] = V_in[1];
		destination[0] = destination_in[0];
		destination[1] = destination_in[1];
		desired_speed = desired_speed_in;
	}

	float* get_desired_velocity()
	{
		float* ans = diff(X, destination);
		float distance = dist(ans);
		if (distance <= 0.001) {
			ans[0] = 0.0;
			ans[1] = 0.0;
			return ans;
		}
		ans[0] = ans[0] / distance * desired_speed;
		ans[1] = ans[1] / distance * desired_speed;
		return ans;
	}

	void update(float dt, float acceleration[]) {
		X[0] = X[0] + V[0] * dt;
		X[1] = X[1] + V[1] * dt;
		V[0] = V[0] + acceleration[0] * dt;
		V[1] = V[1] + acceleration[1] * dt;
	}
};