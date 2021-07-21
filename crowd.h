#ifndef CROWD_H
#define CROWD_H

#include <QGraphicsItem>
#include "particle.h"
class Crowd : public QGraphicsItem
{
public:
    Crowd();
    Crowd(int);
    ~Crowd();
    // QGraphicsItem interface
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void addParticle(QGraphicsItem*,float,float,float,int,float disease=0);
    void setCount(int);
    int getCount();
//    void setParticleStatus(int i,std::vector<float>);
    void setParticleStatus(int i,struct status);
//    int getWindowH();
//    int getEllipseSize();
//    void setEllipseSize(int);
private:
    int count;
//    int ellipse_size;
    QVector<Particle*> particles;
//    int window_h;
};

#endif // CROWD_H
