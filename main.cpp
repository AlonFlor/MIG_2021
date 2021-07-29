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
    //int num_steps = ceil(total_time / dt);
    //int test_steps = 1;
    w.setSimulation(&sim);

//    //main loop
//    int img_count = 1;
//    for (int i = 0; i < num_steps; ++i) {
//        if (i % 100 == 0) {
//            if(test_steps<0) break;
//            std::cout<<"main::outputCSV at step "<<i<<std::endl;
//            sim.outputCSV(img_count);
//            img_count++;
//            test_steps--;
//            sim.update(&w);
//        }
//        sim.update();
//    }
//    std::cout << "Done." << std::endl;
    return a.exec();;
}
