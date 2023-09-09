/*
 * Jopaul Jobi  CID: 02023052
 * Coursework Part 2 Q5)
 * 
 */
#include "header.h"

//Declarations of functions to be used in program (defined after main funtion):
double vel(const double& ti, const double& xi, const double& xdoti);
double acc_x(const double& xi, const double& yi, const vector<Body>& bodies, const double& G, const int& N_Bod, int i);
double acc_y(const double& xi, const double& yi, const vector<Body>& bodies, const double& G, const int& N_Bod, int i);
void rungekutta4(vector<Body>& bodies, const parameters& param, const int& N_Bod, const double& ti);


int main()
{
    //initialise variables
    string line;                    //Used as temporary variable to store each line of input file
    parameters param1;              //Creates parameters struct to store values
    vector<Body> bodies;            //Used to store info about each body, length is adjustable for any number of bodies
    unsigned int N_Iterations = 0;  //Stores number of total time steps to run the simulation
    unsigned int N_Bodies = 0;      //Stores number of total bodies in the given system
    double ti = 0.0;                //Stores current time step of the simulation
    
    //reads parameters file as input
    ifstream InputFile("parameters.txt", ios::in);
    //Ensures the input file is readable and will output an error if not
    if (InputFile.good()) {
        
        //First line will always be T, G and delta_t values 
        getline(InputFile, line);
        stringstream ssline(line);
        ssline >> param1.G >> param1.T >> param1.delta_t;
        
        //iterate until end of file is reached
        while (true) {
            
            /* This section of code gets the line of initiial conditions data for each body, adds it 
             * into the struct 'body' and then adds it to the vector of bodies*/
            getline(InputFile, line);
            stringstream ssline(line);
            Body body;
            ssline >> body.xi >> body.yi >> body.xdoti >> body.ydoti >> body.mi;
            bodies.push_back(body);
            
            if (InputFile.eof()) {
                break;
                }
            }
        }
        
    else {
        cout << "Error: Failed to open file" << endl;
        return 1;
    }
    InputFile.close();
    
    //Calculates number of iterations for the sim and evaulates number of bodies in the system
    N_Iterations = param1.T / param1.delta_t;
    N_Bodies = bodies.size();
    
    //Opens output file to write outputs to, will clear any previous output file that already exists
    ofstream OutputFile("output.txt", ios::out | ios::trunc);
    
    /*This loop executes for each iteration of the sim where the current values of data (i.e current time step, 
     * positions etc) for each body are written to the file and then it executes the 4th order runge kutta 
     * function to evaulate the data at the next time step and repeats until desired time */
    for (unsigned int n = 0; n <= N_Iterations; n++) {
        
        for (unsigned int i = 0; i < N_Bodies; i++) {
            
            OutputFile << (i + 1) << " " << ti
            << " " << bodies[i].xi << " " << bodies[i].yi
            << " " << bodies[i].xdoti << " " << bodies[i].ydoti << endl;
        }
        
        rungekutta4(bodies, param1, N_Bodies, ti);
        ti = ti + param1.delta_t;
    }
    OutputFile.close();
    
    return 0;
}


//Definitions of funtions used:

//Funtion to calculate the velocity of a body from the given ODE
double vel(const double& ti, const double& xi, const double& xdoti) {
    return xdoti;
}

//Function to calculate the x acceleration of a body due to the gravitational force of all other bodies
double acc_x(const double& xi, const double& yi, const vector<Body>& bodies, const double& G, const int& N_Bod, int i) {
    //Initialise variables to be used in function
    double a = 0.0;
    double a_ij = 0.0;
    
    //Executes formula of Newtons Law of graviation
    for (int j = 0; j < N_Bod; j++) {
        if (i == j) continue;
        double dx = (bodies[j].xi - xi);
        double dy = (bodies[j].yi - yi);
        double r = sqrt(dx*dx + dy*dy);
        a_ij = G * (bodies[j].mi * (bodies[j].xi - xi)) / (r*r*r);
        a += a_ij;
    }
    return a;
}

//Function to calculate the y acceleration of a body due to the gravitational force of all other bodies
double acc_y(const double& xi, const double& yi, const vector<Body>& bodies, const double& G, const int& N_Bod, int i) {
    //Initialise variables to be used in function
    double a = 0.0;
    double a_ij = 0.0;
    
    //Executes formula of Newtons Law of graviation
    for (int j = 0; j < N_Bod; j++) {
        if (i == j) continue;
        double dx = (bodies[j].xi - xi);
        double dy = (bodies[j].yi - yi);
        double r = sqrt(dx*dx + dy*dy);
        a_ij = G * (bodies[j].mi * (bodies[j].yi - yi)) / (r*r*r);
        a += a_ij;
    }
    return a;
}

/*Funtion that uses 4th order runge kutta method to evaluate data (i.e position and velcoty) for
 *next time step for each body in the system */
void rungekutta4(vector<Body>& bodies, const parameters& param, const int& N_Bod, const double& ti) {
    
    //Extracts data from param object to be used in program
    const double h = param.delta_t;
    const double G = param.G;
    
    /*Initialise arrays that will store the values at each stage of 4th order runge kutta method,
    * they have a length of 4 as there are 4 data points to be evaluated for each body (x position,
    * y position, x velocity and y velocity) */
    double k1s[4];
    double k2s[4];
    double k3s[4];
    double k4s[4];
    
    //Initialising temporary variables to be used in calculations throughout the stages of RK4
    double xi = 0.0;
    double yi = 0.0;
    double xdoti = 0.0;
    double ydoti = 0.0;
    double x_half = 0.0;
    double y_half = 0.0;
    double xdot_half = 0.0;
    double ydot_half = 0.0;
    
    //Loops for each body in system
    for (int i = 0; i < N_Bod; i++) {
        
        //Gets curent data of body and stores in temporary variables
        xi = bodies[i].xi;
        yi = bodies[i].yi;
        xdoti = bodies[i].xdoti;
        ydoti = bodies[i].ydoti;
        
        //Evaluates first stage values of RK4 method
        k1s[1] = h * vel(ti, xi, xdoti);
        k1s[2] = h * vel(ti, yi, ydoti);
        k1s[3] = h * acc_x(xi, yi, bodies, G, N_Bod, i);
        k1s[4] = h * acc_y(xi, yi, bodies, G, N_Bod, i);
        
        //Calculates inputs to be used in second stage of RK4
        x_half = xi + 0.5 * k1s[1];
        y_half = yi + 0.5 * k1s[2];
        xdot_half = xdoti + 0.5 * k1s[3];
        ydot_half = ydoti + 0.5 * k1s[4];
        
        //Evaluates second stage values of RK4 method
        k2s[1] = h * vel(ti + 0.5 * h, x_half, xdot_half);
        k2s[2] = h * vel(ti + 0.5 * h, y_half, ydot_half);
        k2s[3] = h * acc_x(x_half, y_half, bodies, G, N_Bod, i);
        k2s[4] = h * acc_y(x_half, y_half, bodies, G, N_Bod, i);
        
        //Calculates inputs to be used in third stage of RK4
        x_half = xi + 0.5 * k2s[1];
        y_half = yi + 0.5 * k2s[2];
        xdot_half = xdoti + 0.5 * k2s[3];
        ydot_half = ydoti + 0.5 * k2s[4];

        //Evaluates third stage values of RK4 method
        k3s[1] = h * vel(ti + 0.5 * h, x_half, xdot_half);
        k3s[2] = h * vel(ti + 0.5 * h, y_half, ydot_half);
        k3s[3] = h * acc_x(x_half, y_half, bodies, G, N_Bod, i);
        k3s[4] = h * acc_y(x_half, y_half, bodies, G, N_Bod, i);

        //Calculates inputs to be used in fourth stage of RK4
        x_half = xi + k3s[1];
        y_half = yi + k3s[2];
        xdot_half = xdoti + k3s[3];
        ydot_half = ydoti + k3s[4];

        //Evaluates fourth stage values of RK4 method
        k4s[1] = h * vel(ti + h, x_half, xdot_half);
        k4s[2] = h * vel(ti + h, y_half, ydot_half);
        k4s[3] = h * acc_x(x_half, y_half, bodies, G, N_Bod, i);
        k4s[4] = h * acc_y(x_half, y_half, bodies, G, N_Bod, i);

        /*Calculates new data values for the body (final stage of RK4) and updates the 
         *values in the 'bodies' vecotor which has been passeed by refrence*/
        bodies[i].xi = xi + (k1s[1] + (2 * k2s[1]) + (2 * k3s[1]) + k4s[1]) / 6.0;
        bodies[i].yi = yi + (k1s[2] + (2 * k2s[2]) + (2 * k3s[2]) + k4s[2]) / 6.0;
        bodies[i].xdoti = xdoti + (k1s[3] + (2 * k2s[3]) + (2 * k3s[3]) + k4s[3]) / 6.0;
        bodies[i].ydoti = ydoti + (k1s[4] + (2 * k2s[4]) + (2 * k3s[4]) + k4s[4]) / 6.0;
    }
}