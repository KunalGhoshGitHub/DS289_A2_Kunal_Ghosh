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

// This function calculates the analytical solution
void analytical_solution(vector<real> &x, real t, real alpha, vector <real> &Analytical_solution)
{
    for (int i = 0;i < x.size();i++)
    {
        Analytical_solution[i] = ((exp(-16.0*alpha*t))*sin(4.0*x[i])) + ((exp(-alpha*t))*sin(x[i]));
    }
}

// This function solves the system linear of equations using the Jacobi Iterative Method
void jacobi_solver(vector<vector<real>> &A, vector<real> &b, real tol, vector<real> &x)
{
    // Calculating the number of rows in A matrix
    int n = A.size();

    // Declaring a new vector to store the updated value of the unkowns during the iterations
    vector<real> xnew(n, 0.0);

    // Tolerance achieved (diff)
    real diff = tol + 1.0;

    // This loop will continue as long as the tolerance achieved is greater than the required tolerance
    while (diff > tol)
    {
        for (int i = 0; i < n; ++i)
        {
            real sum = b[i];
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    sum -= A[i][j] * x[j];
                }
            }
            xnew[i] = sum / A[i][i];
        }

        // Calculating the achieved tolerance
        diff = 0.0;
        for (int i = 0; i < n; ++i)
        {
            diff = diff + abs(xnew[i]-x[i]);
        }

        // Updating the value of unknowns
        x = xnew;
    }
}

void A_matrix(vector<vector<real>> &A,real r_d,int N)
{
    real rows = N-1;

    A[0][0] = (1.0+(2.0*r_d));
    A[0][1] = -r_d;
    A[0][N-2] = -r_d;
    for (int i = 1;i < rows-1;i++)
    {
        A[i][i] = (1.0+(2.0*r_d));
        A[i][i-1] = -r_d;
        A[i][i+1] = -r_d;
    }
    A[N-2][N-2] = (1.0+(2.0*r_d));
    A[N-2][N-2-1] = -r_d;
    A[N-2][0] = -r_d;
}

int main()
{
    vector<real> Value_of_N_s,Value_of_rd_s,Inputs;
    vector<string> Output_file_names;
    real dx,X_l,X_u,pi,Tolerance;

    // Reading the inputs from the input file
    Inputs = input_parameters("Input.txt");
    pi = 4.0*atan(1.0);
    X_l = Inputs[0];
    X_u = (Inputs[1])*pi;
    real alpha = Inputs[2];
    real t_end = Inputs[3];
    Tolerance = Inputs[4];

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
            real rows,cols;

            // Calculating the value of dt
            real dt;
            dt = r_d*dx*dx/alpha;

            // Calculating the number of the time steps required
            int num_steps = t_end/(dt);

            rows = N-1;
            cols = N-1;
            vector<real> u_n;
            vector<real> u_n_1;
            vector<vector <real> > A(rows,vector<real> (cols,0.0));
            vector<real> b;

            // Adjustments for the last time step
            real dt_new = t_end - (dt*num_steps);
            real r_d_new = dt_new*alpha/(dx*dx);

            // Initial Conditions
            for (int i = 0;i < x.size()-1;i++)
            {
                u_n.push_back(sin(4.0*x[i])+ sin(x[i]));
            }

            // Calculating the RHS
            b = u_n;

            // Calculating the A matrix
            A_matrix(A,r_d,N);

            // Time marching loop
            for (int t_steps = 0;t_steps < num_steps;t_steps ++)
            {
                vector<real> x_sol(b.size(),0.0);
                // Using Jacobi Iterative Method to solve the system of equations
                jacobi_solver(A, b,Tolerance,x_sol);
                b = x_sol;
            }

            // Calculating the A matrix (Adjustments for the last time step)
            vector<real> x_sol(b.size(),0.0);
            A_matrix(A,r_d_new,N);

            // Solving for the last time step
            jacobi_solver(A, b,Tolerance,x_sol);
            b = x_sol;


            vector<real> U_sol;
            U_sol = b;

            // Ghost point approach is used
            U_sol.push_back(U_sol[0]);

            // Writing the numerical solution to a file
            write_to_file(U_sol,"Question_5_Implicit_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Question_5_Implicit_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

            // Calculating the analytical solution
            vector<real> Analytical_solution(N,0.0);
            analytical_solution(x,t_end,alpha,Analytical_solution);

            // Writing the analytical solution to a file
            write_to_file(Analytical_solution,"Analytical_Solution_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Analytical_Solution_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

            // Calculating the average of absolute error
            vector <real> average_error;
            average_error.push_back(avg_error_vector(U_sol,Analytical_solution));

            // Writing the average of absolute error to a file
            write_to_file(average_error,"Question_5_Implicit_Average_Error_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
            Output_file_names.push_back("Question_5_Implicit_Average_Error_u_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
        }
    }

    // Writing the names of all the generated files to a file
    write_to_file(Output_file_names,"Output_file_names.csv");
}
