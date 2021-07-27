#include "view.h"

View::View(QGraphicsScene *scene,QWidget *parent): graphicsScene(scene)
{
    Q_UNUSED(parent);
    setScene(scene);
    setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    crowd = new Crowd();
    float scale=8.1;
    scene->setSceneRect(QRect(-50*scale,-40*scale,100*scale,80*scale));
    setMouseTracking(true);
    viewport()->setMouseTracking(true);
    //Top left info display
    QString string="Particle Info\n";
    text=scene->addText(string);
    text->setPos(-400,-320);
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
            break;
        case Qt::Key_B:
            bigStep();
            std::cout<<"View::keyPressEvent B pressed. step="<<getStepCount()<<std::endl;
            break;
    }
}

void View::toggleDisplayMode(){
    crowd->toggleDisplayMode();
}

void View::togglePause(){
    pause=-pause;
}

void View::incrementStep(){
    for(int i=0;i<step_size;++i){
        step_count++;
        simulation->update();
    }
    display();
}

void View::bigStep(){
    for(int i=0;i<630;++i){
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
    //scene()->addItem(crowd);
    addParticles();
    crowd->setTextPointer(text);
}

void View::setStepSize(int s){
    step_size=s;
}

void View::setScale(float s){
    //scene()->setSceneRect(QRect(-50*s,-40*s,100*s,80*s));
    crowd->setScale(s);
}

void View::mousePressEvent(QMouseEvent *event){
    if(!itemAt(event->pos())){
        qDebug("Nothing");
    }else{
        qDebug("Item Selected");
    }
    QGraphicsView::mousePressEvent(event);
}

void View::addParticles(){
    QVector<Particle*> vector=crowd->getParticles();
    for(QVector<Particle*>::iterator i=vector.begin();i!=vector.end();++i){
        scene()->addItem(*i);
    }
}

void View::printItems(){
    QList<QGraphicsItem*> item_list=items();
    QList<QGraphicsItem*>::iterator i;
    for(i=item_list.begin();i!=item_list.end();++i){
        std::cout<<(*i)->pos().x()<<" "<<(*i)->pos().y()<<std::endl;
    }
}
