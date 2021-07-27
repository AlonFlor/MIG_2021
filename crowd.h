#ifndef CROWD_H
#define CROWD_H

#include "particle.h"
class Crowd
{
public:
    Crowd();
    Crowd(int);
    ~Crowd();
    void addParticle(float,float,float,int,float disease=0);
    void toggleDisplayMode();
    void setCount(int);
    void setScale(float);
    int getCount();
    void setParticleStatus(int i,struct status);
    QVector<Particle*> getParticles();
    void setTextPointer(QGraphicsTextItem*);
private:
    int count;
    float scale;
    QVector<Particle*> particles;
};

#endif // CROWD_H
