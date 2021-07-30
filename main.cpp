#include <QApplication>
#include "view.h"
#include "MIGBasicSim.h"
#include "Person.h"
#include <QThread>

int main(int argc, char *argv[]){
    QApplication a(argc,argv);
    QGraphicsScene scene;
    View w(&scene);
    w.show();

    w.setCrowdCount(num_people);
    w.setScale(11);
    w.setStepSize(15);
    //create output directory
    std::string outputDir = ".\\output";
    CreateDirectoryA(outputDir.c_str(),NULL);
    MIGBasicSim sim;
    //one person is infected, people[42].disease = 1.0;
    //sim.setInfection(42,1.0);
    w.setSimulation(&sim);
    return a.exec();;
}
