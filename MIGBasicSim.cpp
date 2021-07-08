// MIGBasicSim.cpp : Defines the entry point for the application.
//

#include "MIGBasicSim.h"
#include "Person.h"
#include <math.h>
#include <stdio.h>
#include <vector>
#include <random>

using namespace std;

float interaction_strength;
float interaction_range;
float relaxation_time;
float desired_speed;

float interaction_radius;
//float max_interation_force = 20;

float dt;
float total_time;
int num_groups;
int people_per_group;
int num_people;

vector<float> group_inits;
vector<float> group_displacements;

void read_parameters(char* parameters_file) {
	FILE* params_file;
	params_file = fopen(parameters_file, "r");
	char nothing[100];
	int also_nothing = 0;

	also_nothing = fscanf(params_file, "%f %s", &interaction_strength, nothing);
	printf("interaction_strength: %f\n", interaction_strength);
	also_nothing = fscanf(params_file, "%f %s", &interaction_range, nothing);
	printf("interaction_range: %f\n", interaction_range);
	also_nothing = fscanf(params_file, "%f %s", &relaxation_time, nothing);
	printf("relaxation_time: %f\n", relaxation_time);
	also_nothing = fscanf(params_file, "%f %s", &desired_speed, nothing);
	printf("desired_speed: %f\n", desired_speed);

	also_nothing = fscanf(params_file, "%f %s", &interaction_radius, nothing);
	printf("interaction_radius: %f\n", interaction_radius);
	also_nothing = fscanf(params_file, "%f %s", &dt, nothing);
	printf("dt: %f\n", dt);
	also_nothing = fscanf(params_file, "%f %s", &total_time, nothing);
	printf("total_time: %f\n", total_time);
	also_nothing = fscanf(params_file, "%d %s", &people_per_group, nothing);
	printf("people_per_group: %d\n", people_per_group);
	also_nothing = fscanf(params_file, "%d %s", &num_groups, nothing);
	printf("num_groups: %d\n", num_groups);
	num_people = num_groups * people_per_group;

	float x_init;
	float y_init;
	float x_disp;
	float y_disp;
	for (int i = 0; i < num_groups; ++i) {
		also_nothing = fscanf(params_file, "%f %f %s", &x_init, &y_init, nothing);
		group_inits.push_back(x_init);
		group_inits.push_back(y_init);
		printf("group %d init: (%f, %f)\n", i, group_inits[2*i], group_inits[2*i + 1]);
		
		also_nothing = fscanf(params_file, "%f %f %s", &x_disp, &y_disp, nothing);
		group_displacements.push_back(x_disp);
		group_displacements.push_back(y_disp);
		printf("group %d displacement: (%f, %f)\n", i, group_displacements[2*i], group_displacements[2*i + 1]);
	}

	fclose(params_file);
}

Person initialize(float initial_X[], float initial_V[], float desired_displacement[], int group_ID)
{
	//Places
	float destination[2];
	destination[0] = initial_X[0] + desired_displacement[0];
	destination[1] = initial_X[1] + desired_displacement[1];
	return Person(initial_X, initial_V, destination, desired_speed, group_ID);
}


float** interaction_force_and_disease_spread(Person p1, Person p2)
{
	float* difference = diff(p2.X, p1.X);
	float distance = dist(difference);
	float mult_factor = interaction_strength * exp(-1 * distance / interaction_range);
	float disease_mult_factor = 1.0;
	if (distance > interaction_radius) {
		mult_factor = 0.0;
		disease_mult_factor = 0.0;
	}
	float force_ans[2];
	force_ans[0] = mult_factor * difference[0] / distance;
	force_ans[1] = mult_factor * difference[1] / distance;

	float disease_ans[1];
	float disease_difference = max(p2.disease - p1.disease, 0.0f);	//if p1 is more more infected than p2, p2 cannot infect p1.
	float distance_sq = distance * distance;
	disease_ans[0] = exp(-1*distance) * disease_difference * disease_mult_factor / distance_sq;

	float* ans[2];
	ans[0] = force_ans;
	ans[1] = disease_ans;

	/*//cap force
	float amount = dist(ans);
	if (amount > max_interation_force) {
		ans[0] = ans[0] / amount * max_interation_force;
		ans[1] = ans[1] / amount * max_interation_force;
	}*/
	return ans;
}

float** net_interaction_force_and_disease_spread(int p1_index, vector<Person>people)
{
	float force[2];
	force[0] = 0.0;
	force[1] = 0.0;
	float disease_change[1];
	disease_change[0] = 0.0;

	for (int i = 0; i < num_people; ++i) {
		if (i != p1_index) {
			float** this_force_and_disease_spread = interaction_force_and_disease_spread(people[p1_index], people[i]);

			float* this_force = this_force_and_disease_spread[0];
			force[0] += this_force[0];
			force[1] += this_force[1];

			float* disease_spread = this_force_and_disease_spread[1];
			disease_change[0] += disease_spread[0];
		}
	}

	float* ans[2];
	ans[0] = force;
	ans[1] = disease_change;
	return ans;
}

int main()
{
	read_parameters("description.txt");

	//initial velocity
	float initial_V[2];
	initial_V[0] = 0.;
	initial_V[1] = 0.;

	//initialize particles
	vector<Person> people;
	random_device r;
	default_random_engine e1(r());
	uniform_real_distribution<float> distribution(0,5);
	printf("num_groups %d\n",num_groups);
	printf("people_per_group %d\n",people_per_group);
	for (int i = 0; i < num_groups; ++i) {
		float init_center[2];
		init_center[0] = group_inits[2*i];
		init_center[1] = group_inits[2*i + 1];
		float disp_center[2];
		disp_center[0] = group_displacements[2*i];
		disp_center[1] = group_displacements[2*i + 1];
		printf("(%f, %f)\t(%f, %f)\n",init_center[0], init_center[1], disp_center[0], disp_center[1]);
		for (int j = 0; j < people_per_group; ++j) {
			float rand_coords[2];
			rand_coords[0] = distribution(e1);
			rand_coords[1] = distribution(e1);

			//set x
			while (dist(rand_coords) > 5) {
				rand_coords[0] = distribution(e1);
				rand_coords[1] = distribution(e1);
			}
			float x[2];
			x[0] = rand_coords[0] + init_center[0];
			x[1] = rand_coords[1] + init_center[1];

			//set displacement
			while (dist(rand_coords) > 5) {
				rand_coords[0] = distribution(e1);
				rand_coords[1] = distribution(e1);
			}
			float displacement[2];
			displacement[0] = rand_coords[0] + disp_center[0];
			displacement[1] = rand_coords[1] + disp_center[1];

			//write particle
			people.push_back(initialize(x, initial_V, displacement, i));
		}
	}

	people[42].disease = 1.0; //one person is infected
	
	//people[0] = initialize(initial_center1, initial_V, displacement1);
	//people[1] = initialize(initial_center2, initial_V, displacement2);

	int num_steps = ceil(total_time / dt);
	
	//main loop
	FILE* fp;
	int img_count = 1;
	for (int i = 0; i < num_steps; ++i) {
		if (i % 100 == 0) {
			//float current_time = (float(i)) * dt;
			char t[50];
			//sprintf(t, "%f.csv", current_time);
			sprintf(t, "%04d.csv",img_count);
			img_count += 1;

			//open file
			fp = fopen(t, "w");

			//write header
			fprintf(fp, "pos_x");
			fprintf(fp, ",pos_y");
			fprintf(fp, ",v_x");
			fprintf(fp, ",v_y");
			fprintf(fp, ",group_ID");
			fprintf(fp, ",disease");
			fprintf(fp, "\n");
		}

		for (int j = 0; j < num_people; ++j) {
			if(i%100==0)fprintf(fp, "%f,%f,%f,%f,%d,%f\n", people[j].X[0], people[j].X[1], people[j].V[0], people[j].V[1], people[j].group_ID, people[j].disease); //print person position and velocity

			float* desired_velocity = people[j].get_desired_velocity();

			float desired_velocity_force[2];
			desired_velocity_force[0] = (desired_velocity[0] - people[j].V[0]) / relaxation_time;
			desired_velocity_force[1] = (desired_velocity[1] - people[j].V[1]) / relaxation_time;

			float** force_and_disease = net_interaction_force_and_disease_spread(j, people);
			float* this_net_interaction_force = force_and_disease[0];
			float this_disease_change = force_and_disease[1][0];

			float acceleration[2];
			acceleration[0] = desired_velocity_force[0] +this_net_interaction_force[0];
			acceleration[1] = desired_velocity_force[1] +this_net_interaction_force[1];

			people[j].update_pos_vel(dt, acceleration);
			people[j].update_disease(dt, this_disease_change);
		}
		//close file
		if(i%100==0)fclose(fp);
	}
	//float current_time = (float(num_steps)) * dt;
	char t[50];
	//sprintf(t, "%f.csv", current_time);
	sprintf(t, "%04d.csv", img_count);

	//open file
	fp = fopen(t, "w");

	//write header
	fprintf(fp, "pos_x");
	fprintf(fp, ",pos_y");
	fprintf(fp, ",v_x");
	fprintf(fp, ",v_y");
	fprintf(fp, ",group_ID");
	fprintf(fp, ",disease");
	fprintf(fp, "\n");

	//print people positions and velocities
	for (int j = 0; j < num_people; ++j) {
		fprintf(fp, "%f,%f,%f,%f,%d,%f\n", people[j].X[0], people[j].X[1], people[j].V[0], people[j].V[1], people[j].group_ID, people[j].disease);
	}
	//close file
	fclose(fp);

	cout << "Done." << endl;
	return 0;
}