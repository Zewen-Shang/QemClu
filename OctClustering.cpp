// VertexClustering.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <vtkDelaunay3D.h>
#include <iostream>
#include <chrono>
#include "Model.h"


using namespace std;

#define START_TIME startTime = chrono::system_clock::now();

#define END_TIME(str) endTime = chrono::system_clock::now();\
dur = endTime - startTime;\
second = chrono::duration<double>(dur);\
cout << str <<"Use time : " << second.count() << endl

using namespace chrono;

//70——10%
//25——20%
//8——40%
//6——50%

int main()
{
    chrono::system_clock::time_point startTime = chrono::system_clock::now();
    chrono::system_clock::time_point endTime = chrono::system_clock::now();
    auto dur = endTime - startTime;
    chrono::duration<double> second(dur);

    
    Model m("cylinder1m");

    START_TIME;
    m.buildOct(8);
    m.addQ();
    END_TIME("build");

    //m.selectBorder(5);
    START_TIME;
    m.getBestPos();
    END_TIME("simplification");

    m.outputVtk("cylinder1m_Oct_40");
}
