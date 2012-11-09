//  main.cpp
//  Runge_Kutta
//
//  Uses fourth order Runge-Kutta algorithm to calculate the trajectory of a derivative given the derivative as a C++ function, initial values y0 and t0, time step dt, and the number of steps n. Returns vector of the t values vector and the y values vector
//
//  Created by Kevin McLaughlin on 10/30/12.
//  Copyright (c) 2012 Kevin McLaughlin. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <vector>
using namespace std;


//iterates through and prints values for a vector of vectors
void vector_of_vectors_iterate(vector<vector<float>>vector_of_vectors){
    
    //To access values
    for(vector<vector<float> >::iterator it = vector_of_vectors.begin(); it != vector_of_vectors.end(); ++it)
    {
        //it is now a pointer to a vector<int>
        for(vector<float>::iterator jt = it->begin(); jt != it->end(); ++jt)
        {
            // jt is now a pointer to an integer.
            cout << *jt<<" ";
        }
        cout << endl;
    }
}//end vector of vector iteration function


//Here is the rk4 algorithm with paramters:the derivative as a C++ function, initial values y0 and t0, time step dt, and the number of steps n. Returns vector of the t values vector and the y values vector
vector<vector<float>> rk4(float (*the_deriv)(float,float),float t0, float y0, float dt, int n){
    vector<float> y_values;
    y_values.push_back(y0);
    
    vector<float> t_values;
    t_values.push_back(t0);
    
    for (int i=0; i<n; i++){
        
        //calculate the k-values for rk4
        float k1 = dt*the_deriv(t0, y0);
        float k2 = dt*the_deriv((t0+dt/2),(y0+k1/2));
        float k3 = dt*the_deriv((t0+dt/2),(y0+k2/2));
        float k4 = dt*the_deriv((t0+dt),(y0+k3));
        
        //calculate y1 based on k-values, add dt to t, append y and t to the respective vectors
        float y1 = y0 +k1/6+k2/3+k3/3+k4/6;
        float t1 = t0+dt;
        y_values.push_back(y1);
        t_values.push_back(t1);
        
        //reset values for next cycle
        y0 = y1;
        t0 = t1;
        
    }//end for loop
    
    vector<vector<float>> t_y_values;
    t_y_values.push_back(t_values);
    t_y_values.push_back(y_values);
    return t_y_values;
    
};//end RK4 Function

//defines the derivate function for use in RK4
float the_derivative(float t, float y){
    
    float out;
    //place the derivate here. Example out = y*t for y'=y*t
    out = -2.3*y;
    return out;
}//end derivative


int main(int argc, const char * argv[])
{
    
    //Set initial values
    float y0=1;
    float t0=0;
    float dt=.0125;
    float tmax = 15;
    int n = int(tmax/dt);
    
    //run the RK4 function
    vector<vector<float>> t_y_values;
    t_y_values = rk4(the_derivative, t0, y0, dt, n);
    
    //print the t, y values
    vector_of_vectors_iterate(t_y_values);
    
    
    return 0;
}//end main


