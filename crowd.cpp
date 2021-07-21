#include <QPainter>
#include "crowd.h"
using namespace std;

Crowd::Crowd(){
    float radius=0.3;
    setCount(2);
    addParticle(this,0,590,radius,0);
    addParticle(this,590,1,radius,1);
}

Crowd::Crowd(int count){
    float radius=0.5;
    setCount(count);
    Particle* particle;
    for(int i=0;i<getCount();++i){
        float x=rand()%600;
        float y=rand()%600;
        particle = new Particle(this,x,y,radius);
        if(i<getCount()/2){
            particle->setGroup(0);
        }else{
            particle->setGroup(1);
        }
        float disease = rand()%100*0.01;
        particle->setDisease(disease);
        //cout<<"particle disease="<<particle->getDisease()<<"=>"<<disease<<endl;
        particles.push_back(particle);
        cout<<"constructor 2 particle["<<i<<"]="<<particle<<", x="<<particle->getX()<<", y="<<particle->getY()<<endl;
    }
}

Crowd::~Crowd(){
    Particle* particle;
    foreach(particle,particles)
        delete particle;
}

QRectF Crowd::boundingRect() const{
    return QRectF(0,0,10,10);
}

void Crowd::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
    Q_UNUSED(painter);
    Q_UNUSED(option);
    Q_UNUSED(widget);
}

void Crowd::addParticle(QGraphicsItem *parent,float x,float y,float radius,int group,float disease){
    Particle* particlePtr;
    particlePtr=new Particle(parent,x,y,radius);
    particlePtr->setGroup(group);
    particlePtr->setDisease(disease);
    particles.push_back(particlePtr);
}

void Crowd::setCount(int n){
    count = n;
}

int Crowd::getCount(){
    return count;
}

//void Crowd::setParticleStatus(int index,vector<float> status){
//    particles[index]->setStatus(status);
//}

void Crowd::setParticleStatus(int i,struct status s){
    particles[i]->setStatus(s);
}

//int Crowd::getWindowH(){
//    return window_h;
//}

//int Crowd::getEllipseSize(){
//    return ellipse_size;
//}

//void Crowd::setEllipseSize(int s){
//    ellipse_size = s;
//}