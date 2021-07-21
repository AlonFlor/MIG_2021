#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <QGraphicsEllipseItem>
#include <QGraphicsItem>
#include <QList>
#include <vector>
#include <string>

struct status{
    float x{NAN};
    float y{NAN};
    float dx{NAN};
    float dy{NAN};
    int group{int(NAN)};
    float disease{NAN};
};

class Particle : public QGraphicsEllipseItem
{
public:
    Particle();
    Particle(QGraphicsItem*,float,float,float);
    ~Particle();
    // QGraphicsItem interface
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

    float getX();
    float getY();
    float getDirection();
    int getGroup();
    float getDisease();

    void setX(float);
    void setY(float);
    void setXY(float,float);
    void setDirection(float,float);
    void setGroup(int);
    void setDisease(float);

//    void setStatus(float,float,float,float,int,float);
//    void setStatus(float*);
//    void setStatus(std::vector<float>);
    void setStatus(struct status);
//    std::vector<float> getStatus();
    struct status getStatus();
//    static std::vector<float> getStatusFormat(float const x=NULL,float const y=NULL,float const dx=NULL,float const dy=NULL,int const group=NULL,float const disease=NULL);
//    void appendHistory();
    static std::string statusString(struct status s);
private:
//    float x=0;
//    float y=0;
//    int group=0;
//    float disease=0;
//    QList<float*> history;
//    std::vector<float> direction={0,0};
    float radius;
    struct status s;
};

#endif // PARTICLE_H
