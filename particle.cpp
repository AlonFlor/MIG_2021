#include "particle.h"
#include <QPainter>
using namespace std;
Particle::Particle()
{
    setX(0);
    setY(0);
}

Particle::Particle(QGraphicsItem *parent,float x,float y)
{
    setParentItem(parent);
    setPos(x,y);
    setX(x);
    setY(y);
}

Particle::~Particle()
{
}

QRectF Particle::boundingRect() const{
    return QRectF(0,0,10,10);
}

void Particle::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
//    QPen pen;
//    pen.setCosmetic(true);
//    painter->setPen(pen);
//    painter->setBrush(QColor(255,0,0));
//    if(group==1){
//        painter->setBrush(QColor(0,0,255));
//    }
//    painter->drawEllipse(QRect(x,y,5,5));
//    painter->setBrush(QColor(255*(1-disease),255*(1-disease),255*(1-disease)));
//    painter->drawEllipse(QRect(x+640,y,5,5));
    QPen pen;
    pen.setCosmetic(true);
    painter->setPen(pen);
    painter->setBrush(QColor(255,0,0));
    if(group==1){
        painter->setBrush(QColor(0,0,255));
    }
    painter->drawEllipse(QRect(0,0,5,5));
    painter->setBrush(QColor(255*(1-disease),255*(1-disease),255*(1-disease)));
    painter->drawEllipse(QRect(640,0,5,5));
}

float Particle::getX()
{
    return x;
}

void Particle::setX(float x)
{
    this->x = x;
}

float Particle::getY()
{
    return y;
}

void Particle::setY(float y)
{
    this->y = y;
}

float Particle::getDirection()
{
    return 0;
}

void Particle::setDirection(float dx,float dy)
{
}

void Particle::setGroup(int g)
{
    group = g;
}

int Particle::getGroup()
{
    return group;
}

void Particle::setDisease(float d)
{
    disease = d;
}

float Particle::getDisease()
{
    return disease;
}

void Particle::setStatus(float x,float y,float dx,float dy,int group,float dss)
{
    setPos(x,y);
    setDirection(dx,dy);
    setGroup(group);
    setDisease(dss);
    appendHistory();
    //cout<<"Particle::setStatus_v1 history.size()="<<history.size()<<endl;
}

void Particle::setStatus(float* status)
{
    float *f;
    f=status;
    setX(*f);
    setY(*(f+1));
    setDirection(*(f+2),*(f+3));
    setGroup(int(*(f+4)));
    setDisease(*(f+5));
    appendHistory();
    //cout<<"Particle::setStatus_v2 history.size()="<<history.size()<<endl;
}

void Particle::setStatus(vector<float> s)
{
    setPos(s[0],s[1]);
    setDirection(s[2],s[3]);
    setGroup(s[4]);
    setDisease(s[5]);
    appendHistory();
    //cout<<"Particle::setStatus_v3 x="<<s[0]<<", y="<<s[1]<<", group="<<s[4]<<", step="<<history.size()<<endl;
}

void Particle::appendHistory()
{
    float* fp;
    float status[5]={getX(),getY(),getDirection(),float(getGroup()),getDisease()};
    fp = status;
    history.push_back(fp);
}
