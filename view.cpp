#include "view.h"

/*
View::View(QWidget *parent)
    : QGraphicsView(parent)
{
    setRenderHint(QPainter::Antialiasing);
    setScene(new QGraphicsScene);
    setStepCount(0);
    setPauseCondition(1);
//    for(int particleNum=0;particleNum<getParticleCount();++particleNum) {
//        QGraphicsEllipseItem *circleItem = new QGraphicsEllipseItem(QRect(rand()%600,rand()%600,10,10));
//        circleItem->setBrush(QColor(255,0,0));
//        scene()->addItem(circleItem);
//    }

    crowd = new Crowd();
    scene()->addItem(crowd);
//    setAlignment(Qt::AlignTop|Qt::AlignLeft);
    scene()->setSceneRect(QRect(0,0,600,600));
//    QAbstractScrollArea::setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
//    QAbstractScrollArea::setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setScaleX(10);
    setOffsetX(360);
    setScaleY(10);
    setOffsetY(300);
    setWindowTitle("Crowd Simulation");
}
*/

View::View(QGraphicsScene *scene,QWidget *parent): graphicsScene(scene)
{
    Q_UNUSED(parent);
    setScene(scene);
    setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    crowd = new Crowd();
    scene->addItem(crowd);
    float scale=8.1;
    scene->setSceneRect(QRect(-50*scale,-40*scale,100*scale,80*scale));
    setMouseTracking(true);
    viewport()->setMouseTracking(true);
}

View::~View()
{
    delete scene();
    delete crowd;
}

//void View::wheelEvent(QWheelEvent *event)
//{
//    //Q_UNUSED(event);
//    QGraphicsView::wheelEvent(event);
//    if(event->isAccepted()){
//        return;
//    }
//    const float factor = 1.1;
//    if(event->angleDelta().y()>0){
//        scale(factor,factor);
//    }else{
//        scale(1/factor,1/factor);
//    }
//    event->accept();
//}

void View::keyPressEvent(QKeyEvent *event){
    if(event->isAutoRepeat()) {
        return;
    }
    switch(event->key()) {
        case Qt::Key_Space:
            togglePause();
            std::cout<<"View::keyPressEvent space pressed. pause="<<getPauseCondition()<<std::endl;
            break;
        case Qt::Key_S:
            incrementStep();
            std::cout<<"View::keyPressEvent S pressed. step="<<getStepCount()<<std::endl;
            break;
        case Qt::Key_R:
            restart();
            std::cout<<"View::keyPressEvent R pressed. step="<<getStepCount()<<std::endl;
            break;
        case Qt::Key_D:
            toggleDisplayMode();
            std::cout<<"View::keyPressEvent D pressed. step="<<getStepCount()<<std::endl;
    }
}

void View::toggleDisplayMode(){
    crowd->toggleDisplayMode();
}

void View::togglePause(){
    pause = -pause;
}

void View::incrementStep(){
    for(int i=0;i<step_size;++i){
        step_count++;
        simulation->update();
    }
    display();
}

/**
 * @brief display simulation
 */
void View::display(){
    int numPeople=simulation->getNumPeople();
    for(int i=0;i<numPeople;++i){
        struct status s=simulation->getStatus(i);
        crowd->setParticleStatus(i,s);
    }
}

void View::restart(){
    setStepCount(0);
}

int View::getStepCount(){
    return step_count;
}

int View::getPauseCondition(){
    return pause;
}

void View::setPauseCondition(int p){
    this->pause = p;
}

void View::setSimulation(MIGBasicSim* sim){
    simulation=sim;
    display();
}

void View::setStepCount(int s)
{
    this->step_count = s;
}

void View::setCrowdCount(int count){
    delete crowd;
    crowd = new Crowd(count);
    scene()->addItem(crowd);
}

void View::setStepSize(int s){
    step_size=s;
}

void View::setScale(float s){
    crowd->setScale(s);
}
