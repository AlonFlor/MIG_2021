#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsTextItem>
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

class Particle : public QGraphicsItem
{
public:
    Particle();
    Particle(float,float,float);
    ~Particle();
    // QGraphicsItem interface
    virtual QRectF boundingRect() const override;
    virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;
    void toggleDisplayMode();

    float getX();
    float getY();
    float getR() const;
    float getDirection();
    int getGroup();
    float getDisease();
    struct status getStatus();

    void setX(float);
    void setY(float);
    void setXY(float,float);
    void setDirection(float,float);
    void setGroup(int);
    void setDisease(float);
    void setStatus(struct status);
    void setScale(float);
    void setID(int);
//    void appendHistory();
    static std::string statusString(struct status s);
    void setTextPointer(QGraphicsTextItem *);
    void setStatusDisplayString();
protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *) override;
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *) override;
    void mousePressEvent(QGraphicsSceneMouseEvent*) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *) override;
private:
//    float x=0;
//    float y=0;
//    int group=0;
//    float disease=0;
//    QList<float*> history;
//    std::vector<float> direction={0,0};
    int id;
    float radius;
    struct status s;
    int display_mode;
    float scale;
    bool hovered;
    bool selected;
    QGraphicsTextItem *text;
    QString statusDisplayString;
};

#endif // PARTICLE_H
