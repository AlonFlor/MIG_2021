#include "Person.h"

Person::Person(int i,std::vector<float> X,std::vector<float> V,std::vector<float> displacement,float speed,int group){
    people_ID=i;
    this->X=X;
    this->V=V;
    group_ID=group;
    destination.push_back(X[0]+displacement[0]);
    destination.push_back(X[1]+displacement[1]);
    setDesiredSpeed(speed);
    setDesiredVelocity();
}

Person::Person(int i,std::vector<float> X,std::vector<float> V,std::vector<float> dest,float speed,int group,float rs,bool wear_mask,bool healthy){
    people_ID=i;
    this->X=X;
    this->V=V;
    X0=X;
    V0=V;
    group_ID=group;
    destination=dest;
    setDesiredSpeed(speed);
    setDesiredVelocity();
    mask=(wear_mask)?(0.9*rs):0.0;
    immunity=0.5;
    immunity_strength=(healthy)?0.1:-0.1;
}

void Person::update(float dt, std::vector<std::vector<float>> net,float relaxation_time){
    std::vector<float> net_interaction_force=net[0];
    std::vector<float> desired_velocity_force=getDesiredVelocityForce(relaxation_time);
    std::vector<float> a={\
        desired_velocity_force[0]+net_interaction_force[0],\
        desired_velocity_force[1]+net_interaction_force[1]};
    updateXV(dt,a);
    float disease_change=net[1][0];
    updateDisease(dt,disease_change);
}

void Person::update(std::vector<std::vector<float>> net,float t,struct simData d){
    std::vector<float> net_interaction_force=net[0];
    std::vector<float> desired_velocity_force=getDesiredVelocityForce(d.relaxation_time);
    std::vector<float> a={\
        desired_velocity_force[0]+net_interaction_force[0],\
        desired_velocity_force[1]+net_interaction_force[1]};
    updateXV(d.dt,a);
    float disease_change=net[1][0];
    updateDisease(d,disease_change);
    updateImmunity(d.dt,t,d.vaccine_factor);
    //printStatus();
}

/**
 * @brief update position and velocity
 * @param dt
 * @param a
 */
void Person::updateXV(float dt,std::vector<float> a){
    X[0]=X[0]+V[0]*dt;
    X[1]=X[1]+V[1]*dt;
    V[0]=V[0]+a[0]*dt;
    V[1]=V[1]+a[1]*dt;
}

void Person::updateDisease(float dt,float disease_change){
    disease+=disease_change*dt;
    disease=clamp(disease,0,1);
}

//with second disease change
void Person::updateDisease(struct simData d,float disease_change){
    float disease_sq=disease*disease;
    float second_disease_change=0.;
    if(d.disease_a!=0){second_disease_change=disease*(d.disease_s-disease)-d.disease_a*disease_sq*immunity/(disease_sq+1);}
    disease_change+=second_disease_change;
    disease+=disease_change*d.dt;
    disease=clamp(disease,0,1);
}

void Person::updateImmunity(float dt,float current_time,float vaccine_factor){
    float disease_sq=getDiseaseSq();
    float immunity_change=-0.5*disease_sq/(1.0+disease_sq)*immunity+immunity_strength*immunity+vaccine_factor*tanh(current_time);
    immunity+=immunity_change*dt;
    immunity=clamp(immunity,0,1);
}

void Person::reset(){
    X=X0;
    V=V0;
    disease=0;
}

void Person::printStatus(){
    printf("p[%d],X=[%f,%f],V=[%f,%f],disease=%f\n",people_ID,X[0],X[1],V[0],V[1],disease);
}

float Person::getDesiredSpeed(){
    return desired_speed;
}

std::vector<float> Person::getDesiredVelocity(){
    return desired_velocity;
}

std::vector<float> Person::getDesiredVelocityForce(float relaxation_time){
    std::vector<float> f={(desired_velocity[0]-V[0])/relaxation_time,
                          (desired_velocity[1]-V[1])/relaxation_time};
    return f;
}

float Person::getDisease(){
    return disease;
}

float Person::getDiseaseSq(){
    return disease*disease;
}

status Person::getStatus(){
    status s{X[0],X[1],V[0],V[1],group_ID,disease};
    return s;
}

int Person::getID(){
    return people_ID;
}

float Person::getMask(){
    return mask;
}

void Person::setDesiredSpeed(float s){
    desired_speed=s;
}

void Person::setDesiredVelocity(){
    std::vector<float> df=diff(X,destination);
    float dst=dist(X,destination);
    std::vector<float> velocity;
    if(dst<=0.001){
        velocity.push_back(0.0);
        velocity.push_back(0.0);
    }else{
        velocity.push_back(df[0]/dst*getDesiredSpeed());
        velocity.push_back(df[1]/dst*getDesiredSpeed());
    }
    desired_velocity=velocity;
}

void Person::setDisease(int d){
    disease=d;
}

std::vector<float> Person::diff(std::vector<float> a, std::vector<float> b){
    std::vector<float> d={b[0]-a[0],b[1]-a[1]};
    return d;
}

float Person::dist(std::vector<float> a, std::vector<float> b){
    float x1=a[0];
    float y1=a[1];
    float x2=b[0];
    float y2=b[1];
    return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

float Person::dist(Person* b){
    float x1=X[0];
    float y1=X[1];
    float x2=b->X[0];
    float y2=b->X[1];
    return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

float Person::clamp(float input,float lo,float hi){
    if(input<lo)return lo;
    if(input>hi)return hi;
    return input;
}
