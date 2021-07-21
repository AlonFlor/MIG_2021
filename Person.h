#ifndef PERSON_H
#define PERSON_H
#include <math.h>
#include <vector>
#include <string>
#include <stdio.h>
#include "particle.h"
class Person{
public:
    static const int dim = 2;

    //float X[dim];
    //float V[dim];
    //float destination[dim];
    std::vector<float> X;
    std::vector<float> V;
    std::vector<float> destination;
    float desired_speed;
    std::vector<float> desired_velocity;
    float disease = 0.0;
    float mask = 0.0;

    Person();
    //Person(float* X_in, float* V_in, float* destination_in, float desired_speed_in, int group_ID_in);
    Person(int,std::vector<float>,std::vector<float>,std::vector<float>,float,int);
    float getDesiredSpeed();
    void setDesiredSpeed(float);
    //float* get_desired_velocity();
    std::vector<float> getDesiredVelocity();
    void setDesiredVelocity();
    //void update_pos_vel(float dt, float acceleration[]);
    void update(float,std::vector<std::vector<float>>);
    void updateXV(float,std::vector<float>);
    void updateDisease(float dt,float disease_in);
    float getDisease();
    void setDisease(int);
    //std::vector<float> getStatus();
    status getStatus();
    int getID();
    static std::vector<float> diff(std::vector<float> a,std::vector<float> b);
    float dist(Person*);
    float dist(std::vector<float>,std::vector<float>);
private:
    int people_ID;
    int group_ID;
};

#endif // PERSON_H
