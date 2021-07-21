#include "particle.h"
#include <QPainter>
using namespace std;

Particle::Particle(QGraphicsItem *parent,float x,float y,float radius)
    : QGraphicsEllipseItem(parent)
{
    s.x=x;
    s.y=y;
    this->radius=radius;

    setFlags(ItemIsSelectable);
    setAcceptHoverEvents(true);
    //TODO: display info with hover event
}

QRectF Particle::boundingRect() const{
    return QRectF(0,0,radius*2,radius*2);
}

void Particle::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
    Q_UNUSED(widget);
    Q_UNUSED(option);

    QPen pen;
    qreal width = radius*0.2;
    pen.setWidthF(width);
    painter->setPen(pen);
    QBrush b = painter->brush();
        painter->setBrush(b);

    if(s.group==0){
        painter->setBrush(QColor(255,0,0));
    }else if(s.group==1){
        painter->setBrush(QColor(0,0,255));
    }
    painter->drawEllipse(QRect(s.x,s.y,radius*2,radius*2));
//    painter->setBrush(QColor(255*(1-disease),255*(1-disease),255*(1-disease)));
//    painter->drawEllipse(QRect(640,0,radius*2,radius*2));
}

//getters
float Particle::getX(){
    return s.x;
}

float Particle::getY(){
    return s.y;
}

float Particle::getDirection(){
    return 0;
}

int Particle::getGroup(){
    return s.group;
}

void Particle::setDisease(float d){
    s.disease = d;
}

float Particle::getDisease(){
    return s.disease;
}

//setters
void Particle::setX(float x){
    this->s.x = x;
}

void Particle::setY(float y){
    this->s.y = y;
}

void Particle::setXY(float x,float y){
    this->s.x = x;
    this->s.y = y;
}

void Particle::setDirection(float dx,float dy){
    this->s.dx=dx;
    this->s.dy=dy;
}

void Particle::setGroup(int g){
    s.group = g;
}

//void Particle::setStatus(float x,float y,float dx,float dy,int group,float dss){
//    setPos(x,y);
//    setDirection(dx,dy);
//    setGroup(group);
//    setDisease(dss);
//    //appendHistory();
//    //cout<<"Particle::setStatus_v1 history.size()="<<history.size()<<endl;
//}

//void Particle::setStatus(float* status){
//    float *f;
//    f=status;
//    setX(*f);
//    setY(*(f+1));
//    setDirection(*(f+2),*(f+3));
//    setGroup(int(*(f+4)));
//    setDisease(*(f+5));
//    //appendHistory();
//    //cout<<"Particle::setStatus_v2 history.size()="<<history.size()<<endl;
//}

//void Particle::setStatus(vector<float> s){
//    setPos(s[0],s[1]);
//    setDirection(s[2],s[3]);
//    setGroup(s[4]);
//    setDisease(s[5]);
//    //appendHistory();
//    //cout<<"Particle::setStatus_v3 x="<<s[0]<<", y="<<s[1]<<", group="<<s[4]<<", step="<<history.size()<<endl;
//}

// set multiple properties of a particle at a time
//void Particle::setStatus(float x,float y,float dx,float dy,int group,float disease){
//    if(x){this->x=x;}
//    if(y){this->y=y;}
//    if(dx){this->direction[0]=dx;}
//    if(dy){this->direction[1]=dy;}
//    if(group){this->group=group;}
//    if(disease){this->disease=disease;}
//}

void Particle::setStatus(struct status s){
    this->s=s;
    update();
}

//void Particle::appendHistory(){
//    float* fp;
//    float status[5]={getX(),getY(),getDirection(),float(getGroup()),getDisease()};
//    fp = status;
//    history.push_back(fp);
//}

struct status Particle::getStatus(){
    return s;
}

// all possible status of a particle
//vector<float> Particle::getStatusFormat(float const x,float const y,float const dx,float const dy,int const group,float const disease){
//    vector<float> status={x,y,dx,dy,float(group),disease};
//    return status;
//}

std::string Particle::statusString(struct status s){
    std::string str=
            std::to_string(s.x)+','+
            std::to_string(s.y)+','+
            std::to_string(s.dx)+','+
            std::to_string(s.dy)+','+
            std::to_string(s.group)+','+
            std::to_string(s.disease);
    return str;
}
