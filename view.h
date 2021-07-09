#ifndef VIEW_H
#define VIEW_H

#include <QWidget>
#include <QGraphicsView>
#include <QWheelEvent>
#include <QGraphicsEllipseItem>
#include "crowd.h"

class View : public QGraphicsView
{
    Q_OBJECT

public:
    View(QWidget *parent=nullptr);
    ~View();
    void wheelEvent(QWheelEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void setPauseCondition(int);
    void togglePause();
    int getPauseCondition();
    void setStepCount(int);
    void incrementStep();
    int getStepCount();
    void restart();
    void setCrowdCount(int);
//    void setParticle(int,float,float,float,float,int,float);
    void setParticle(int,std::vector<float>);
    void setScaleX(int);
    void setOffsetX(int);
    void setScaleY(int);
    void setOffsetY(int);
private:
    int particle_count;
    int step_count;
    int pause;
    Crowd* crowd;
    int window_h;
    int scaleX;
    int offsetX;
    int scaleY;
    int offsetY;
};

#endif // VIEW_H
