#ifndef VIEW_H
#define VIEW_H

#include <QWidget>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QWheelEvent>
#include <QGraphicsEllipseItem>
#include "crowd.h"
#include "MIGBasicSim.h"

class View : public QGraphicsView
{
    Q_OBJECT

public:
//    View(QWidget *parent=nullptr);
    View(QGraphicsScene *scene=nullptr,QWidget *parent=nullptr);
    ~View();
//    void wheelEvent(QWheelEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void togglePause();
    void incrementStep();
    void update();
    void restart();

    // getters
    int getPauseCondition();
    int getStepCount();

    // setters
    void setSimulation(MIGBasicSim* sim);
    void setPauseCondition(int);
    void setStepCount(int);
    //void setParticle(int,float,float,float,float,int,float);
    void setCrowdCount(int);
    void setParticle(int,std::vector<float>);
    void setScaleX(int);
    void setOffsetX(int);
    void setScaleY(int);
    void setOffsetY(int);
private:
    QGraphicsScene* graphicsScene;
    Crowd* crowd;
    MIGBasicSim* simulation;
    int particle_count;
    int step_count;
    int pause;
    int window_h;
    int scaleX;
    int offsetX;
    int scaleY;
    int offsetY;
};

#endif // VIEW_H
