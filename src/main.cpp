#include <string>
#include <iostream>
#include <fstream>

#include "main.h"
#include "Camera.h"
#include "Runge Kutta.h"
#include "Runge Kutta Data.h"
#include "Black Hole.h"
#include "Files.h"

using namespace std;

int main(){
    //Tester();
    RunUniverse(((string)"C:\\dosgames\\kerr schild\\" + INPUT).c_str());
	return 0;
}

void RunUniverse(const char *file)
{
    ifstream in;
    in.open(file);
    int cameras, i, j, maxiterations;
    double initstep, maxdistance, errscale;
    double position[4];
    double cameralook[16];
    string error;
    string output;
    int skip = 0;
    int coordinate;
    double M, a;
    Camera3D *cam;
    RungeKuttaForm *solver;
    BlackHole *bh;

    int pixh, pixv;
    double fovh, fovv, zoomh1, zoomh2, zoomv1, zoomv2;

    solver = RungeKutta::DormandPrince;
    /*solver = RungeKutta::Midpoint;
    solver = RungeKutta::Kutta;
    solver = RungeKutta::BogackiSchampine;
    solver = RungeKutta::CashKarp;
    solver = RungeKutta::Fehlberg;/**/

    in >> error;
    in >> output;

    in >> cameras;
    in >> initstep;
    in >> errscale;
    errscale *= 100;
    cout << "cameras: " << cameras << ". error: " << error << ". output: " << output << "." << endl;

    in >> M;
    in >> a;

    bh = new BlackHole(M, a);

    for(i = 0; i < cameras; i++)
    {
        in >> maxiterations;
        in >> maxdistance;

        in >> pixh;
        in >> pixv;
        in >> fovh;
        in >> fovv;
        in >> zoomh1;
        in >> zoomh2;
        in >> zoomv1;
        in >> zoomv2;
        in >> coordinate;

        for(j = 0; j < 4; j++) in >> position[j];

        for(j = 0; j < 16; j++) in >> cameralook[j];

        cam = new Camera3D(bh, position, cameralook, coordinate, CAMERA_FLAT, pixh, pixv, fovh, fovv, zoomh1, zoomh2, zoomv1, zoomv2);

        if(skip<1) {
            cam->Snapshot(solver, ((string)"c:\\dosgames\\kerr schild\\" + output +  stringify(i) + ".txt").c_str(), maxiterations>MAXSTEPS?MAXSTEPS:maxiterations, maxdistance, initstep, errscale, coordinate, CoordinateTest);
            cout << "snapshot complete!" << endl << flush;
        } else {skip--;}

        cout << "Image " << i << " complete!" << endl << flush;
        delete cam;
    }

    in.close();
}
