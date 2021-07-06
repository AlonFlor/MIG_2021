// MIGBasicSim.cpp : Defines the entry point for the application.
//

#include "MIGBasicSim.h"
#include "Person.h"
#include <math.h>

using namespace std;

float interaction_strength = 1.0;
float interaction_range = 1.0;
float relaxation_time = 1.0;
float desired_speed = 1.0;

float interaction_radius = 0.5;
//float max_interation_force = 20;

const float dt = 0.001;
const float total_time = 100;
const int num_people = 50;


Person initialize(float initial_X[], float initial_V[], float desired_displacement[])
{
	//Places
	float destination[2];
	destination[0] = initial_X[0] + desired_displacement[0];
	destination[1] = initial_X[1] + desired_displacement[1];
	return Person(initial_X, initial_V, destination, desired_speed);
}


float* interaction_force(Person p1, Person p2)
{
	float* ans = diff(p1.X, p2.X);
	float distance = dist(ans);
	float mult_factor = interaction_strength * exp(-1 * distance / interaction_range);
	if (distance > interaction_radius) {
		mult_factor = 0.0;
	}
	ans[0] = mult_factor * ans[0] / distance;
	ans[1] = mult_factor * ans[1] / distance;

	/*//cap force
	float amount = dist(ans);
	if (amount > max_interation_force) {
		ans[0] = ans[0] / amount * max_interation_force;
		ans[1] = ans[1] / amount * max_interation_force;
	}*/
	return ans;
}

float* net_interaction_force(int p1_index, Person people[])
{
	float force[2];
	force[0] = 0.0;
	force[1] = 0.0;
	for (int i = 0; i < num_people; ++i) {
		if (i != p1_index) {
			float* this_force = interaction_force(people[i], people[p1_index]);
			force[0] += this_force[0];
			force[1] += this_force[1];
		}
	}
	return force;
}

int main()
{
	//centers
	float initial_center1[2];
	initial_center1[0] = 1.0;
	initial_center1[1] = 1.0;
	float initial_center2[2];
	initial_center2[0] = 4.0;
	initial_center2[1] = 0.5;

	//displacements
	float displacement1[2];
	displacement1[0] = 1.0;
	displacement1[1] = 0.0;
	float displacement2[2];
	displacement2[0] = -8.0;
	displacement2[1] = 0.0;

	//initial velocity
	float initial_V[2];
	initial_V[0] = 0.;
	initial_V[1] = 0.;

	//initialize particles
	Person people[num_people];
	int count = 0;
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			float x[2];
			x[0] = initial_center1[0] - j;
			x[1] = initial_center1[1] + i;
			people[count] = initialize(x, initial_V, displacement1);
			count += 1;
		}
	}
	count = 0;
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			float x[2];
			x[0] = initial_center2[0] + j;
			x[1] = initial_center2[1] + i;
			people[count + 25] = initialize(x, initial_V, displacement2);
			count += 1;
		}
	}
	//people[0] = initialize(initial_center1, initial_V, displacement1);
	//people[1] = initialize(initial_center2, initial_V, displacement2);

	int num_steps = ceil(total_time / dt);
	
	//main loop
	FILE* fp;
	for (int i = 0; i < num_steps; ++i) {
		if (i % 100 == 0) {
			float current_time = (float(i)) * dt;
			char t[50];
			sprintf(t, "%f.csv", current_time);

			//open file
			fp = fopen(t, "w");

			//write header
			fprintf(fp, "pos_x");
			fprintf(fp, ",pos_y");
			fprintf(fp, ",v_x");
			fprintf(fp, ",v_y");
			fprintf(fp, "\n");
		}

		for (int j = 0; j < num_people; ++j) {
			if(i%100==0)fprintf(fp, "%f,%f,%f,%f\n", people[j].X[0], people[j].X[1], people[j].V[0], people[j].V[1]); //print person position and velocity

			float* desired_velocity = people[j].get_desired_velocity();

			float desired_velocity_force[2];
			desired_velocity_force[0] = (desired_velocity[0] - people[j].V[0]) / relaxation_time;
			desired_velocity_force[1] = (desired_velocity[1] - people[j].V[1]) / relaxation_time;

			float* this_net_interaction_force = net_interaction_force(j, people);

			float acceleration[2];
			acceleration[0] = desired_velocity_force[0] +this_net_interaction_force[0];
			acceleration[1] = desired_velocity_force[1] +this_net_interaction_force[1];

			people[j].update(dt, acceleration);
		}
		//close file
		if(i%100==0)fclose(fp);
	}
	float current_time = (float(num_steps)) * dt;
	char t[50];
	sprintf(t, "%f.csv", current_time);

	//open file
	fp = fopen(t, "w");

	//write header
	fprintf(fp, "pos_x");
	fprintf(fp, ",pos_y");
	fprintf(fp, ",v_x");
	fprintf(fp, ",v_y");
	fprintf(fp, "\n");
	fprintf(fp, "%f", current_time);

	//print people positions and velocities
	for (int j = 0; j < num_people; ++j) {
		fprintf(fp, "%f,%f,%f,%f\n", people[j].X[0], people[j].X[1], people[j].V[0], people[j].V[1]);
	}
	//close file
	fclose(fp);

	cout << "Done." << endl;
	return 0;
}