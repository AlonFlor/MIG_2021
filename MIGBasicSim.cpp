#include "MIGBasicSim.h"
#include <math.h>

using namespace std;

MIGBasicSim::MIGBasicSim()
{
//    //centers
//    float initial_center1[2];
//    initial_center1[0] = -15.0;
//    initial_center1[1] = 0.0;
//    float initial_center2[2];
//    initial_center2[0] = 15.0;
//    initial_center2[1] = 0.0;

//    //displacements
//    float displacement1[2];
//    displacement1[0] = 30.0;
//    displacement1[1] = 0.0;
//    float displacement2[2];
//    displacement2[0] = -30.0;
//    displacement2[1] = 0.0;

//    //initial velocity
//    float initial_V[2];
//    initial_V[0] = 0.;
//    initial_V[1] = 0.;

//    //initialize particles
//    Person people[num_people];
//    int count = 0;
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 10; ++j) {
//            float x[2];
//            x[0] = initial_center1[0] - j;
//            x[1] = initial_center1[1] + i;
//            people[count] = MIGBasicSim::initialize(x, initial_V, displacement1,0);
//            count += 1;
//        }
//    }
//    count = 0;
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 10; ++j) {
//            float x[2];
//            x[0] = initial_center2[0] + j;
//            x[1] = initial_center2[1] + i;
//            people[count + 50] = MIGBasicSim::initialize(x, initial_V, displacement2,1);
//            count += 1;
//        }
//    }

    interaction_strength = 5.0;
    interaction_range = 1.0;
    relaxation_time = 1.0;
    desired_speed = 1.0;

    interaction_radius = 1.5;
    //max_interation_force = 20.0;

    vector<float> center1 = {-15,0};
    vector<float> center2 = {15,0};
    initial_centers.push_back(center1);
    initial_centers.push_back(center2);

    vector<float> disp1 = {30,0};
    vector<float> disp2 = {-30,0};
    displacements.push_back(disp1);
    displacements.push_back(disp2);

    initial_velocity.push_back(0.0);
    initial_velocity.push_back(0.0);

    numPeople = num_people;
    group = 2;
    row_in_group = 5;
    column_in_group = numPeople/(group*row_in_group);
    int p_index = 0;
    for(int i=0;i<row_in_group;++i){
        for(int j=0;j<column_in_group;++j){
            vector<float> X;
            X.push_back(initial_centers[0][0]-j);
            X.push_back(initial_centers[0][1]+i);
            Person* p;
            p = new Person(p_index,X,initial_velocity,displacements[0],desired_speed,0);
            ++p_index;
            people.push_back(p);
        }
    }
    for(int i=0;i<row_in_group;++i){
        for(int j=0;j<column_in_group;++j){
            vector<float> X;
            X.push_back(initial_centers[1][0]-j);
            X.push_back(initial_centers[1][1]+i);
            Person* p;
            p = new Person(p_index,X,initial_velocity,displacements[1],desired_speed,1);
            ++p_index;
            people.push_back(p);
        }
    }
}

MIGBasicSim::~MIGBasicSim()
{
    for(auto p:people){
        delete p;
    }
    people.clear();
}

//float** MIGBasicSim::interaction_force_and_disease_spread(Person p1, Person p2)
//{
//    float* difference = diff(p2.X, p1.X);
//    float distance = dist(difference);
//    float mult_factor = interaction_strength * exp(-1 * distance / interaction_range);
//    float disease_mult_factor = 1.0;
//    if (distance > interaction_radius) {
//        mult_factor = 0.0;
//        disease_mult_factor = 0.0;
//    }
//    float force_ans[2];
//    force_ans[0] = mult_factor * difference[0] / distance;
//    force_ans[1] = mult_factor * difference[1] / distance;

//    float disease_ans[1];
//    float disease_difference = max(p2.disease - p1.disease, 0.0f);	//if p1 is more more infected than p2, p2 cannot infect p1.
//    float distance_sq = distance * distance;
//    disease_ans[0] = exp(-1*distance) * disease_difference * disease_mult_factor / distance_sq;

//    float* ans[2];
//    ans[0] = force_ans;
//    ans[1] = disease_ans;

//    /*//cap force
//    float amount = dist(ans);
//    if (amount > max_interation_force) {
//        ans[0] = ans[0] / amount * max_interation_force;
//        ans[1] = ans[1] / amount * max_interation_force;
//    }*/
//    return ans;
//}

vector<vector<float>> MIGBasicSim::interaction_force_and_disease_spread(Person* p1, Person* p2){
    vector<float> difference = Person::diff(p1->X,p2->X);
    float distance = p1->dist(p2);
    float mult_factor = interaction_strength * exp(-1 * distance / interaction_range);
    float disease_mult_factor = 1.0;
    if(distance>interaction_radius){
//        mult_factor = 0;
//        disease_mult_factor = 0.0;
        return vector<vector<float>> {{0,0},{0,0}};
    }

    vector<float> force = {mult_factor*difference[0]/distance, mult_factor*difference[1]/distance};
    float d_disease = max(p2->getDisease()-p1->getDisease(),0.0f); //if p1 is more more infected than p2, p2 cannot infect p1.
    float distance_sq = distance*distance;
    vector<float> disease = {(exp(-1*distance)*d_disease*disease_mult_factor)/distance_sq};

    vector<vector<float>> force_disease;
    force_disease.push_back(force);
    force_disease.push_back(disease);
    return force_disease;
}

//float** MIGBasicSim::net_interaction_force_and_disease_spread(int p1_index, Person people[])
//{
//    float force[2];
//    force[0] = 0.0;
//    force[1] = 0.0;
//    float disease_change[1];
//    disease_change[0] = 0.0;

//    for (int i = 0; i < num_people; ++i) {
//        if (i != p1_index) {
//            float** this_force_and_disease_spread = interaction_force_and_disease_spread(people[p1_index], people[i]);

//            float* this_force = this_force_and_disease_spread[0];
//            force[0] += this_force[0];
//            force[1] += this_force[1];

//            float* disease_spread = this_force_and_disease_spread[1];
//            disease_change[0] += disease_spread[0];
//        }
//    }

//    float* ans[2];
//    ans[0] = force;
//    ans[1] = disease_change;
//    return ans;
//}

vector<vector<float>> MIGBasicSim::net_interaction_force_and_disease_spread(Person* one){
    vector<float> force = {0.0,0.0};
    vector<float> disease_change = {0.0};

    for(auto p: people){
        if(p!=one){
            vector<vector<float>> current_spread = interaction_force_and_disease_spread(one,p);
            vector<float> current_force = current_spread[0];
            force[0]+=current_force[0];
            force[1]+=current_force[1];
            disease_change[0]+=current_spread[1][0];
        }
    }
    vector<vector<float>> force_disease;
    force_disease.push_back(force);
    force_disease.push_back(disease_change);
    return force_disease;
}

vector<vector<float>> MIGBasicSim::concurrent_net_interaction_force_and_disease_spread(Person* one){
    vector<float> force = {0.0,0.0};
    vector<float> disease_change = {0.0};

    for(auto p: people){
        if(p!=one){
            vector<vector<float>> current_spread = interaction_force_and_disease_spread(one,p);
            vector<float> current_force = current_spread[0];
            force[0]+=current_force[0];
            force[1]+=current_force[1];
            disease_change[0]+=current_spread[1][0];
        }
    }
    vector<vector<float>> force_disease;
    force_disease.push_back(force);
    force_disease.push_back(disease_change);
    return force_disease;
}

void MIGBasicSim::update(){
    for(auto p:people){
        updateOnePerson(p);
//        if(v!=nullptr)v->setParticle(p->getID(),p->getStatus());
    }
}

void MIGBasicSim::updateOnePerson(Person* one){
    vector<float> desired_velocity = one->getDesiredVelocity();
    vector<float> desired_velocity_force = Person::diff(desired_velocity,one->V);
    vector<vector<float>> force_disease = net_interaction_force_and_disease_spread(one);
    one->update(dt,force_disease);
}

void MIGBasicSim::outputCSV(int i){
    //zero padding
    stringstream ss;
    ss<<setw(4)<<setfill('0')<<i;
    string filename(".\\output\\"+ss.str()+".csv");

    fstream fout;
    fout.open(filename,fstream::out);
    fout<<"pos_x,pos_y,v_x,v_y,group_ID,disease"<<endl<<flush;
    for(auto p:people){
        status status = p->getStatus();
        fout<<Particle::statusString(status)<<flush;
        fout<<endl<<flush;
    }
    fout.close();
}

int MIGBasicSim::getNumPeople()
{
    return numPeople;
}

void MIGBasicSim::setNumPeople(int n)
{
    numPeople = n;
}

void MIGBasicSim::setInfection(int p_index, int disease)
{
    if(p_index<int(people.size())){
        Person* p = people[p_index];
        p->setDisease(disease);
    }
}

float MIGBasicSim::getInteractionStrength(){
    return interaction_strength;
}

void MIGBasicSim::setInteractionStrength(float s){
    interaction_strength = s;
}

float MIGBasicSim::getInteractionRange(){
    return interaction_range;
}

void MIGBasicSim::setInteractionRange(float r){
    interaction_range = r;
}

float MIGBasicSim::getRelaxationTime(){
    return relaxation_time;
}

void MIGBasicSim::setRelaxationTime(float t){
    relaxation_time = t;
}

float MIGBasicSim::getDesiredSpeed(){
    return desired_speed;
}

void MIGBasicSim::setDesiredSpeed(float s){
    desired_speed = s;
}

float MIGBasicSim::getInteractionRadius(){
    return interaction_radius;
}

void MIGBasicSim::setInteractionRadius(float r){
    interaction_radius = r;
}

float MIGBasicSim::getMaxInteractionForce(){
    return max_interation_force;
}

void MIGBasicSim::setMaxInteractionForce(float f){
    max_interation_force = f;
}

status MIGBasicSim::getStatus(int p_index){
    return people[p_index]->getStatus();
}
