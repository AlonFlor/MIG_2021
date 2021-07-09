#include <QApplication>
#include "view.h"
#include "MIGBasicSim.h"
#include "Person.h"

//int main(int argc, char *argv[])
//{
//    QApplication a(argc,argv);
//    View w;
//    w.setCrowdCount(20);
//    float status[6]={0,0,0,0,1,0};
//    w.setParticle(19,status);
//    float status2[6]={3,3,0,0,0,1};
//    w.setParticle(0,status2);
//    w.show();
//    return a.exec();
//}

int main(int argc, char *argv[])
{
    //view with Qt
    QApplication a(argc,argv);
    View w;
    w.setCrowdCount(20);
    std::vector<float> status={0,0,0,0,1,0};
    w.setParticle(19,status);
    std::vector<float> status2={3,3,0,0,0,1};
    w.setParticle(0,status2);
    w.show();
    w.setCrowdCount(num_people);

    //create output directory
    std::string outputDir = ".\\output";
    CreateDirectoryA(outputDir.c_str(),NULL);
    MIGBasicSim sim;
    //one person is infected, people[42].disease = 1.0;
    sim.setInfection(42,1.0);
    int num_steps = ceil(total_time / dt);
    int test_steps = 3;
    //main loop
    int img_count = 1;
    for (int i = 0; i < num_steps; ++i) {
        if (i % 100 == 0) {
            if(test_steps<0) break;
            std::cout<<"main::outputCSV at step "<<i<<std::endl;
            sim.outputCSV(img_count);
            test_steps--;
            img_count++;
            sim.update(&w);
        }
        sim.update();
    }

    std::cout << "Done." << std::endl;

    return a.exec();
}
