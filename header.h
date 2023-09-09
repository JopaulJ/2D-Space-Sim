#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;

//Defintions of objects to be used in program:

//An object for storing parameter values to be used in the program
struct parameters {
    //Initialise variables in parameters struct 
    double G = 0.0;
    double T = 0.0;
    double delta_t = 0.0;
};

//An object for storing initial conditions and mass of a body
struct Body {
    //Initialise variables in body struct
    double xi = 0.0;
    double yi = 0.0;
    double xdoti = 0.0;
    double ydoti = 0.0;
    double mi = 0.0;
};
