#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <QGraphicsItem>
#include <QList>
#include <vector>
class Particle : public QGraphicsItem
{
public:
    Particle();
    Particle(QGraphicsItem*,float,float);
    // QGraphicsItem interface
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    ~Particle();
    float getX();
    void setX(float);
    float getY();
    void setY(float);
    float getDirection();
    void setDirection(float,float);
    void setGroup(int);
    int getGroup();
    void setDisease(float);
    float getDisease();
    void setStatus(float,float,float,float,int,float);
    void setStatus(float*);
    void setStatus(std::vector<float>);
    void appendHistory();
private:
    float x;
    float y;
    int group;
    float disease;
    QList<float*> history;
};

#endif // PARTICLE_H
