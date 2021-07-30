#include "MIGBasicSim.h"
#include <math.h>

MIGBasicSim::MIGBasicSim(){
    time_till_now=0;
    readDescriptionTxt();
    initializeCrowd();
}

MIGBasicSim::~MIGBasicSim(){
    for(auto p:people){
        delete p;
    }
    people.clear();
}

void MIGBasicSim::readDescriptionTxt(){
    social_force_strength=5.0;
    social_force_range_b=1.0;
    social_force_radius=1.0;
    contact_force_strength=10.0;
    particle_radius=0.05;
    disease_a=0;
    disease_s=1.0;
    disease_radius=2.0;
    check_radius=std::max(social_force_radius,disease_radius)+0.1;
    disease_range_rs=1.0;
    mask_prob=0;
    healthy_life_prob=0.5;
    vaccine_factor=0;
    relaxation_time=1.0;
    desired_speed=1.0;
    dt=0.01;
    total_time=100;
    time_steps_per_frame=10;
    group_init_radius=10;
    people_per_group=50;
    num_groups=2;
    num_people=num_groups*people_per_group;

    std::vector<float> init_coords_1={-30,0};
    std::vector<float> displacement_coords_1={600,0};
    std::vector<float> init_coords_2={30,0};
    std::vector<float> displacement_coords_2={-600,0};
    initial_centers={init_coords_1,init_coords_2};
    displacements={displacement_coords_1,displacement_coords_2};

    sim_data={dt,relaxation_time,vaccine_factor,disease_a,disease_s};
}

void MIGBasicSim::initializeCrowd(){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> distribution(-1*group_init_radius,group_init_radius);
    std::uniform_real_distribution<float> probability(0,1);
    for(int g=0;g<num_groups;++g){
        std::vector<float> init_center=initial_centers[g];
        std::vector<float> disp_center=displacements[g];
        //printf("(%f, %f)\t(%f, %f)\n",init_center[0], init_center[1], disp_center[0], disp_center[1]);
        for(int i=0;i<people_per_group;++i){
            std::vector<float> rand_coords={distribution(mt),distribution(mt)};
            std::vector<float> origin={0,0};
            while(Person::dist(rand_coords,origin)>group_init_radius)
                rand_coords={distribution(mt),distribution(mt)};

            std::vector<float> X0={rand_coords[0]+init_center[0],rand_coords[1]+init_center[1]}; //position
            std::vector<float> V0={0.0,0.0}; //velocity
            //std::vector<float> destination={X0[0]+rand_coords[0]+disp_center[0],X0[1]+rand_coords[1]+disp_center[1]};
            std::vector<float> destination={X0[0]+disp_center[0],X0[1]+disp_center[1]};
            bool masked=probability(mt)<mask_prob;
            bool healthy_life=probability(mt)<healthy_life_prob;
            int p_index=g*people_per_group+i;

            Person* p;
            p=new Person(p_index,X0,V0,destination,desired_speed,g,disease_range_rs,masked,healthy_life);
            //printf("p[%d],(%f,%f),group=%d,destination=(%f,%f)\n",p_index,X0[0],X0[1],g,destination[0],destination[1]);
            //printf("desired velocity=(%f,%f)\n",p->getDesiredVelocity()[0],p->getDesiredVelocity()[1]);
            people.push_back(p);
        }
    }
}

void MIGBasicSim::reset(){
    time_till_now=0;
    Person* p;
    for(int i=0;i<int(people.size());++i){
        p=people[i];
        p->reset();
    }
}

void MIGBasicSim::initializeCrowdNoMask(){
    row_in_group=5;
    column_in_group=num_people/(num_groups*row_in_group);
    int p_index=0;
    for(int i=0;i<row_in_group;++i){
        for(int j=0;j<column_in_group;++j){
            std::vector<float> X;
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
            std::vector<float> X;
            X.push_back(initial_centers[1][0]-j);
            X.push_back(initial_centers[1][1]+i);
            Person* p;
            p = new Person(p_index,X,initial_velocity,displacements[1],desired_speed,1);
            ++p_index;
            people.push_back(p);
        }
    }
}

std::vector<std::vector<float>> MIGBasicSim::interaction_force_and_disease_spread(Person *p1, Person *p2){
    std::vector<float> difference=Person::diff(p2->X,p1->X);
    float distance=p1->dist(p2);
    //skip if particles are distant
    if(distance>check_radius)
        return std::vector<std::vector<float>>{{0,0},{0,0}};

    std::vector<float> vel_difference={(p2->X[0]-p1->X[0])*dt,(p2->X[1]-p1->X[1])*dt};
    float vel_distance=Person::dist(vel_difference,std::vector<float>{0,0});
    std::vector<float> vel_pos_difference={difference[0]-vel_difference[0],difference[1]-vel_difference[1]};
    float vel_pos_distance=Person::dist(vel_pos_difference,std::vector<float>{0,0});

    //social force
    float distance_plus_vel_pos_distance=distance+vel_pos_distance;
    float b_ij=0.5*std::sqrt(distance_plus_vel_pos_distance*distance_plus_vel_pos_distance-vel_distance);
    float middle_factor=0.5*distance_plus_vel_pos_distance/b_ij;
    float mult_factor=social_force_strength*middle_factor*std::exp(-b_ij/social_force_range_b);
    if(distance>social_force_radius)
        mult_factor=0.0;
    std::vector<float> force={
        mult_factor*((difference[0]/distance)+(vel_pos_difference[0]/vel_pos_distance)),
        mult_factor*((difference[1]/distance)+(vel_pos_difference[1]/vel_pos_distance))};

    //contact force
    float contact_mult_factor=contact_force_strength*std::max(0.0f,2*particle_radius-distance);
    force[0]+=contact_mult_factor*difference[0]/distance;
    force[1]+=contact_mult_factor*difference[1]/distance;

    //disease
    float disease_mult_factor=1.0;
    if(distance>disease_radius)
        disease_mult_factor=0.0;
    float disease_difference=std::max(p2->getDisease()-p1->getDisease(),0.0f);
    float distance_sq=distance*distance;
    std::vector<float> disease={exp(-1*(distance*p2->getMask())/disease_range_rs)*disease_difference*disease_mult_factor/distance_sq};

    std::vector<std::vector<float>> force_disease={force,disease};
    return force_disease;
}

std::vector<std::vector<float>> MIGBasicSim::net_interaction_force_and_disease_spread(Person* one){
    std::vector<float> force={0.0,0.0};
    std::vector<float> disease_change={0.0};
    for(auto p: people){
        if(p!=one){
            std::vector<std::vector<float>> current=interaction_force_and_disease_spread(one,p);
            std::vector<float> current_force=current[0];
            force[0]+=current_force[0];
            force[1]+=current_force[1];
            disease_change[0]+=current[1][0];
        }
    }
    std::vector<std::vector<float>> force_disease={force,disease_change};
    return force_disease;
}

void MIGBasicSim::update(){
    for(auto p:people){
        std::vector<std::vector<float>> force_disease = net_interaction_force_and_disease_spread(p);
        time_till_now+=dt;
        p->update(force_disease,time_till_now,sim_data);
    }
}

void MIGBasicSim::outputCSV(int i){
    //zero padding
    std::stringstream ss;
    ss<<std::setw(4)<<std::setfill('0')<<i;
    std::string filename(".\\output\\"+ss.str()+".csv");

    std::fstream fout;
    fout.open(filename,std::fstream::out);
    fout<<"pos_x,pos_y,v_x,v_y,group_ID,disease"<<std::endl<<std::flush;
    for(auto p:people){
        status status = p->getStatus();
        fout<<Particle::statusString(status)<<std::flush;
        fout<<std::endl<<std::flush;
    }
    fout.close();
}

int MIGBasicSim::getNumPeople(){
    return num_people;
}

void MIGBasicSim::setNumPeople(int n){
    num_people = n;
}

void MIGBasicSim::setInfection(int p_index, int disease){
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

float MIGBasicSim::getMaxInteractionForce(){
    return max_interation_force;
}

void MIGBasicSim::setMaxInteractionForce(float f){
    max_interation_force = f;
}

status MIGBasicSim::getStatus(int p_index){
    return people[p_index]->getStatus();
}
