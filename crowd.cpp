#include <QPainter>
#include "crowd.h"
using namespace std;

Crowd::Crowd()
{
    window_h=640;
    ellipse_size=5;
    setCount(10);
    Particle* particle;
    for(int i=0;i<getCount();++i){
        particle = new Particle(this,rand()%600,rand()%600);
        if(i<getCount()/2){
            particle->setGroup(0);
        }else{
            particle->setGroup(1);
        }
        float disease = rand()%100*0.01;
        particle->setDisease(disease);
        //cout<<"particle disease="<<particle->getDisease()<<"=>"<<disease<<endl;
        particles.push_back(particle);
        //cout<<"particle="<<particle<<", x="<<particle->getX()<<", y="<<particle->getY()<<endl;
    }
}

Crowd::Crowd(int count)
{
    window_h=640;
    ellipse_size=5;
    setCount(count);
    Particle* particle;
    for(int i=0;i<getCount();++i){
        particle = new Particle(this,rand()%600,rand()%600);
        if(i<getCount()/2){
            particle->setGroup(0);
        }else{
            particle->setGroup(1);
        }
        float disease = rand()%100*0.01;
        particle->setDisease(disease);
        //cout<<"particle disease="<<particle->getDisease()<<"=>"<<disease<<endl;
        particles.push_back(particle);
        //cout<<"particle="<<particle<<", x="<<particle->getX()<<", y="<<particle->getY()<<endl;
    }
}

Crowd::~Crowd()
{
    Particle* particle;
    foreach(particle,particles)
        delete particle;
}

QRectF Crowd::boundingRect() const
{
    return QRectF(0,0,10,10);
}

void Crowd::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
//    QPen pen;
//    pen.setCosmetic(true);
//    painter->setPen(pen);
//    for(Particle* particle: qAsConst(particles)){
//        painter->setBrush(QColor(255,0,0));
//        int group=particle->getGroup();
//        if(group==1){
//            painter->setBrush(QColor(0,0,255));
//        }
//        //cout<<"particle="<<particle<<"x="<<particle->getX()<<", y="<<particle->getY()<<endl;
//        painter->drawEllipse(QRect(particle->getX(),particle->getY(),ellipse_size,ellipse_size));

//        // plot disease
//        qreal disease=particle->getDisease();
//        painter->setBrush(QColor(255*(1-disease),255*(1-disease),255*(1-disease)));
//        painter->drawEllipse(QRect(particle->getX()+640,particle->getY(),ellipse_size,ellipse_size));
//    }
}

void Crowd::setCount(int n){
    count = n;
}

int Crowd::getCount(){
    return count;
}

void Crowd::setParticleStatus(int index,vector<float> status){
    particles[index]->setStatus(status);
}

int Crowd::getWindowH(){
    return window_h;
}

int Crowd::getEllipseSize(){
    return ellipse_size;
}

void Crowd::setEllipseSize(int s){
    ellipse_size = s;
}
