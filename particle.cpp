#include "particle.h"
#include <QPainter>
using namespace std;

Particle::Particle(float x,float y,float radius){
    s.x=x;
    s.y=y;
    this->radius=radius;
    display_mode=0;
    hovered=false;
    setFlag(ItemIsSelectable,true);
    setAcceptHoverEvents(true);
    //TODO: display info with hover event
}

Particle::~Particle(){}

QRectF Particle::boundingRect() const{
    return QRectF(0,0,getR()*2,getR()*2);
}

void Particle::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
    Q_UNUSED(widget);
    Q_UNUSED(option);

    QPen pen;
    float r=getR();
    qreal width = r*0.2;
    pen.setWidthF(width);
    painter->setPen(pen);
    QBrush b = painter->brush();
    painter->setBrush(b);

    //std::cout<<"Particle::paint x="<<x<<", y="<<y<<std::endl;
    if(hovered){
        painter->setBrush(QColor(255,220,140));
    }else if(display_mode==0){
        if(s.group==0){
            painter->setBrush(QColor(255,120,120));
        }else if(s.group==1){
            painter->setBrush(QColor(120,120,255));
        }
    }else if(display_mode==1){
        //painter->setBrush(QColor(255*(1-s.disease),255*(1-s.disease),255*(1-s.disease)));
        painter->setBrush(getDiseaseColor());
    }
    painter->drawEllipse(QRect(0,0,r*2,r*2));
}

void Particle::toggleDisplayMode(){
    display_mode=(display_mode+1)%2;
    //std::cout<<"toggleDisplayMode() display_mode="<<display_mode;
    update();
}

//getters
int Particle::getID(){
    return id;
}

float Particle::getX(){
    return s.x*scale;
}

float Particle::getY(){
    return s.y*scale;
}

float Particle::getR() const{
    return radius*scale;
}

float Particle::getDirection(){
    return 0;
}

int Particle::getGroup(){
    return s.group;
}


float Particle::getDisease(){
    return s.disease;
}

struct status Particle::getStatus(){
    return s;
}

QColor Particle::getDiseaseColor(){
    float d=s.disease;
    float red=255*d;
    float green=255*(1-d);
    float blue=0;
    QColor c(red,green,blue);
    return c;
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

void Particle::setDisease(float d){
    s.disease = d;
}

void Particle::setStatus(struct status s){
    this->s=s;
    setPos(getX(),getY());
    setStatusDisplayString();
    update();
}

void Particle::setScale(float s){
    scale=s;
}

void Particle::setID(int i){
    id=i;
}

//void Particle::appendHistory(){
//    float* fp;
//    float status[5]={getX(),getY(),getDirection(),float(getGroup()),getDisease()};
//    fp = status;
//    history.push_back(fp);
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

void Particle::setTextPointer(QGraphicsTextItem *text){
    this->text=text;
}

void Particle::setStatusDisplayString(){
    QString status="Particle Info\n";
    status+="Index="+QString::number(id)+'\n';
    status+='['+QString::number(s.x)+','+QString::number(s.y)+']'+'\n';
    statusDisplayString=status;
}

void Particle::hoverEnterEvent(QGraphicsSceneHoverEvent *event){
    Q_UNUSED(event);
    hovered=true;
    text->setPlainText(statusDisplayString);
    update();
}

void Particle::hoverLeaveEvent(QGraphicsSceneHoverEvent *event){
    Q_UNUSED(event);
    hovered=false;
    text->setPlainText("Particle Info\n");
    update();
}

void Particle::mousePressEvent(QGraphicsSceneMouseEvent *event){
    Q_UNUSED(event);
    std::cout<<"Particle::press ["<<getX()<<","<<getY()<<"]"<<std::endl;
    update();
}

void Particle::mouseReleaseEvent(QGraphicsSceneMouseEvent *event){
    Q_UNUSED(event);
    update();
}
