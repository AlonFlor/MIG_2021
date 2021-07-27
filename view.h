#ifndef VIEW_H
#define VIEW_H

#include <QWidget>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QWheelEvent>
#include <QGraphicsEllipseItem>
#include <QGraphicsSimpleTextItem>
#include <QPointF>
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
    virtual void keyPressEvent(QKeyEvent *event) override;
    void toggleDisplayMode();
    void togglePause();
    void incrementStep();
    void bigStep();
    void display();
    void restart();

    // getters
    int getPauseCondition();
    int getStepCount();

    // setters
    void setSimulation(MIGBasicSim* sim);
    void setPauseCondition(int);
    void setStepCount(int);
    void setCrowdCount(int);
    void setStepSize(int);
    void setScale(float);

    void addParticles();
    void printItems();
protected:
    virtual void mousePressEvent(QMouseEvent *event) override;
private:
    QGraphicsScene* graphicsScene;
    Crowd* crowd;
    MIGBasicSim* simulation;
    QGraphicsTextItem* text;
    int particle_count;
    int step_count;
    int pause;
    int window_h;
    int step_size;
};

#endif // VIEW_H
