// MIGBasicSim.h : Include file for standard system include files,
// or project specific include files.

#ifndef MIGBASICSIM_H
#define MIGBASICSIM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <Windows.h>
#include "Person.h"
#include "view.h"

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
    void update(View* v=nullptr);
    void updateOnePerson(Person*,View*);
    //Person initialize(float initial_X[], float initial_V[], float desired_displacement[], int group_ID);
    float** interaction_force_and_disease_spread(Person p1, Person p2);
    float** net_interaction_force_and_disease_spread(int p1_index, Person people[]);
    std::vector<std::vector<float>> interaction_force_and_disease_spread(Person*,Person*);
    std::vector<std::vector<float>> net_interaction_force_and_disease_spread(Person*);
    void outputCSV(int);
    int getNumPeople();
    void setNumPeople(int);
    void setInfection(int,int);
    float getInteractionStrength();
    void setInteractionStrength(float);
    float getInteractionRange();
    void setInteractionRange(float);
    float getRelaxationTime();
    void setRelaxationTime(float);
    float getDesiredSpeed();
    void setDesiredSpeed(float);
    float getInteractionRadius();
    void setInteractionRadius(float);
    float getMaxInteractionForce();
    void setMaxInteractionForce(float);
private:
    std::vector<Person*> people;
    int numPeople;
    int group;
    int row_in_group;
    int column_in_group;
    std::vector<std::vector<float>> initial_centers;
    std::vector<std::vector<float>> displacements;
    std::vector<float> initial_velocity;

    float interaction_strength;
    float interaction_range;
    float relaxation_time;
    float desired_speed;
    float interaction_radius;
    float max_interation_force;
};

#endif
