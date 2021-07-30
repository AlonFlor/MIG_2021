#include "view.h"

View::View(QGraphicsScene *scene,QWidget *parent): graphicsScene(scene)
{
    Q_UNUSED(parent);
    graphicsScene=scene;
    setScene(scene);
    setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    crowd = new Crowd();
    float scale=8.5;
    scene->setSceneRect(QRect(-50*scale,-40*scale,100*scale,80*scale));
    setMouseTracking(true);
    viewport()->setMouseTracking(true);

    //Info Display at Top of GUI
    text=addText("Particle Info\n",-400,-230,20);
    addText("(Hover to view)\n",-400,-250,15);
    addText("S - Step Forward    D - Switch Display Mode    R - Reset\n",-400,200,25);
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
    pause=!pause;
}

void View::incrementStep(){
    for(int i=0;i<step_size;++i){
        step_count++;
        simulation->update();
    }
    display();
}

void View::bigStep(){
    for(int i=0;i<800;++i){
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
    simulation->reset();
    display();
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
    addParticles();
    crowd->setTextPointer(text);
}

void View::setStepSize(int s){
    step_size=s;
}

void View::setScale(float s){
    crowd->setScale(s);
}

void View::mousePressEvent(QMouseEvent *event){
    Particle* selected=(Particle*)itemAt(event->pos());
    if(!selected){
        qDebug("Nothing");
    }else{
        int index=selected->getID();
        float disease=1;
        simulation->setInfection(index,disease);
        qDebug("Item Selected");
        struct status s=simulation->getStatus(index);
        crowd->setParticleStatus(index,s);
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

QGraphicsTextItem* View::addText(QString content,float x,float y,int size){
    QGraphicsTextItem* text=graphicsScene->addText(content);
    text->setPos(x,y);
    if(size>0){
        QFont font=text->font();
        font.setPointSize(size);
        text->setFont(font);
    }
    return text;
}
