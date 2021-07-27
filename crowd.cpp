#include <QPainter>
#include "crowd.h"

Crowd::Crowd(){
    float radius=0.3;
    setCount(2);
    addParticle(0,590,radius,0);
    addParticle(590,1,radius,1);
}

Crowd::Crowd(int count){
    float radius=0.5;
    setCount(count);
    Particle* particle;
    for(int i=0;i<getCount();++i){
        float x=rand()%80-40;
        float y=rand()%80-40;
        particle = new Particle(x,y,radius);
        particle->setID(i);
        if(i<getCount()/2){
            particle->setGroup(0);
        }else{
            particle->setGroup(1);
        }
        float disease = rand()%100*0.01;
        particle->setDisease(disease);
        //cout<<"particle disease="<<particle->getDisease()<<"=>"<<disease<<endl;
        particles.push_back(particle);
        //std::cout<<"Crowd-particle["<<i<<"]="<<particle<<", x="<<particle->getX()<<", y="<<particle->getY()<<std::endl;
    }
}

Crowd::~Crowd(){
    Particle* particle;
    foreach(particle,particles)
        delete particle;
}

void Crowd::addParticle(float x,float y,float radius,int group,float disease){
    Particle* particlePtr;
    particlePtr=new Particle(x,y,radius);
    particlePtr->setGroup(group);
    particlePtr->setDisease(disease);
    particles.push_back(particlePtr);
}

void Crowd::toggleDisplayMode(){
    Particle* particle;
    foreach(particle,particles)
        particle->toggleDisplayMode();
}

void Crowd::setCount(int n){
    count = n;
}

void Crowd::setScale(float s){
    Particle* particle;
    foreach(particle,particles)
        particle->setScale(s);
}

int Crowd::getCount(){
    return count;
}

void Crowd::setParticleStatus(int i,struct status s){
    particles[i]->setStatus(s);
}

QVector<Particle*> Crowd::getParticles(){
    return particles;
}

void Crowd::setTextPointer(QGraphicsTextItem* text){
    for(int i=0;i<particles.size();++i){
        particles[i]->setTextPointer(text);
    }
}
