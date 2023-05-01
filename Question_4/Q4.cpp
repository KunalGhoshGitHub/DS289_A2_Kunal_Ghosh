#include <iostream>
#include <cassert>
#include<algorithm>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>
#include "DS289.h"

using namespace std;

typedef double real;

// This function calculates the analytical solution at the grid points
void analytical_solution(vector<real> &x, real t, real alpha, vector <real> &Analytical_solution)
{
    for (int i = 0;i < x.size();i++)
    {
        Analytical_solution[i] = ((exp(-16.0*alpha*t))*sin(4.0*x[i])) + ((exp(-alpha*t))*sin(x[i]));
    }
}

int main()
{
    vector<real> Value_of_N_s,Value_of_rd_s,Inputs;
    vector<string> Output_file_names;
    real dx,X_l,X_u,pi;

    pi = 4.0*atan(1.0);

    // Reading the inputs from the input file
    Inputs = input_parameters("Input.txt");
    X_l = Inputs[0];
    X_u = (Inputs[1])*pi;
    real alpha = Inputs[2];
    real t_end = Inputs[3];

    // Reading the values of N from the file
    Value_of_N_s = input_parameters("Value_of_Ns.txt");

    // Reading the values of rd from the file
    Value_of_rd_s = input_parameters("Value_of_rd_s.txt");

    // This loop iterates over the different values of N
    for (int n = 0;n < Value_of_N_s.size();n++)
    {
        int N = Value_of_N_s[n];

        vector<real> x;

        // Discretize the domain
        dx = (X_u-X_l)/(N-1.0);
        Discretize_1D(x,dx,X_l,X_u);

        // This loop iterates over the different values of rd
        for (int r = 0;r < Value_of_rd_s.size();r++)
        {
            real r_d = Value_of_rd_s[r];

            // Calculating the value of dt
            real dt;
            dt = r_d*dx*dx/alpha;

            // Calculating the number of the time steps required
            int num_steps = t_end/(dt);

            vector<real> u_n;

            // Adjustments for the last time step
            real dt_new = t_end - (dt*num_steps);
            real r_d_new = dt_new*alpha/(dx*dx);

            // Initial Conditions
            for (int i = 0;i < x.size();i++)
            {
                u_n.push_back(sin(4.0*x[i])+ sin(x[i]));
            }

            // Time marching loop
            for (int  k = 0;k < num_steps;k++)
            {
                vector<real> u_n_1(N,0.0);
                u_n_1[0] = (r_d*u_n[N-2]) + ((1.0-(2.0*r_d))*u_n[0])+ (r_d*u_n[1]);
                for (int i = 1;i < N-1-1;i++)
                {
                    u_n_1[i] = (r_d*u_n[i-1]) + ((1.0-(2.0*r_d))*u_n[i])+ (r_d*u_n[i+1]);
                }
                u_n_1[N-2] = (r_d*u_n[N-3]) + ((1.0-(2.0*r_d))*u_n[N-2])+ (r_d*u_n[0]);
                u_n_1[N-1] = u_n_1[0];
                u_n = u_n_1;
            }

            // Calculating the A matrix (Adjustments for the last time step)
            vector<real> u_n_1(N,0.0);

            // Solving for the last time step
            u_n_1[0] = (r_d_new*u_n[N-2]) + ((1.0-(2.0*r_d_new))*u_n[0])+ (r_d_new*u_n[1]);
            for (int i = 1;i < N-1-1;i++)
            {
                u_n_1[i] = (r_d_new*u_n[i-1]) + ((1.0-(2.0*r_d_new))*u_n[i])+ (r_d_new*u_n[i+1]);
            }
            u_n_1[N-2] = (r_d_new*u_n[N-3]) + ((1.0-(2.0*r_d_new))*u_n[N-2])+ (r_d_new*u_n[0]);
            u_n_1[N-1] = u_n_1[0];

            // Writing the numerical solution to a file
            write_to_file(u_n_1,"Question_4_Explicit_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Question_4_Explicit_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

            // Calculating the analytical solution
            vector<real> Analytical_solution(N,0.0);
            analytical_solution(x,t_end,alpha,Analytical_solution);

            // Writing the analytical solution to a file
            write_to_file(Analytical_solution,"Analytical_Solution_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Analytical_Solution_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

            // Calculating the average of absolute error
            vector<real> average_error;
            average_error.push_back(avg_error_vector(u_n_1,Analytical_solution));

            // Writing the average of absolute error to a file
            write_to_file(average_error,"Question_4_Explicit_Average_Error_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Question_4_Explicit_Average_Error_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
        }
    }

    // Writing the names of all the generated files to a file
    write_to_file(Output_file_names,"Output_file_names.csv");
}
