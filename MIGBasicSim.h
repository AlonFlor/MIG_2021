// MIGBasicSim.h : Include file for standard system include files,
// or project specific include files.

#ifndef MIGBASICSIM_H
#define MIGBASICSIM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <random>
#include <Windows.h>
#include "Person.h"
//#include "view.h"

//float interaction_strength = 5.0;
//float interaction_range = 1.0;
//float relaxation_time = 1.0;
//float desired_speed = 1.0;

//float interaction_radius = 1.5;
//float max_interation_force = 20;

const float dt = 0.001;
const float total_time = 100;
const int num_people = 100;

class MIGBasicSim
{
public:
    MIGBasicSim();
    ~MIGBasicSim();
    void readDescriptionTxt();
    void initializeCrowd();
    void reset();
    void initializeCrowdNoMask();
    void update();
    //void updateOnePerson(Person*);
    //Person initialize(float initial_X[], float initial_V[], float desired_displacement[], int group_ID);
    float** interaction_force_and_disease_spread(Person p1, Person p2);
    float** net_interaction_force_and_disease_spread(int p1_index, Person people[]);
    std::vector<std::vector<float>> interaction_force_and_disease_spread(Person*,Person*);
    std::vector<std::vector<float>> net_interaction_force_and_disease_spread(Person*);
    std::vector<std::vector<float>> concurrent_net_interaction_force_and_disease_spread(Person* one);
    void outputCSV(int);

    //setters
    void setNumPeople(int);
    void setInfection(int,int);
    void setInteractionStrength(float);
    void setInteractionRange(float);
    void setRelaxationTime(float);
    void setDesiredSpeed(float);
    void setMaxInteractionForce(float);

    //getters
    int getNumPeople();
    float getInteractionStrength();
    float getInteractionRange();
    float getRelaxationTime();
    float getDesiredSpeed();
    float getMaxInteractionForce();
    status getStatus(int);
private:
    float dt;
    float total_time;
    float time_till_now;
    int time_steps_per_frame;

    std::vector<Person*> people;
    int num_people;
    int num_groups;
    int row_in_group;
    int column_in_group;
    std::vector<float> initial_velocity;

    float interaction_strength;
    float interaction_range;
    float relaxation_time;
    float desired_speed;
    float max_interation_force;

    float social_force_strength;
    float social_force_range_b;
    float social_force_radius;
    float contact_force_strength;

    int disease_a;
    float disease_s;
    float disease_radius;
    float check_radius;
    float disease_range_rs;
    float mask_prob;
    float healthy_life_prob;
    float vaccine_factor;

    float particle_radius;
    float group_init_radius;
    int people_per_group;
    std::vector<std::vector<float>> initial_centers;
    std::vector<std::vector<float>> displacements;

    struct simData sim_data;
};

#endif
