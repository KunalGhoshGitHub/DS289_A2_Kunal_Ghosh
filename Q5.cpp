#include <iostream>
#include <cassert>
#include<algorithm>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>

using namespace std;

typedef double real;

// This function calculates the analytical solution
void analytical_solution(vector<real> &x, real t, real Alpha, vector <real> &Analytical_solution)
{
    for (int i = 0;i < x.size();i++)
    {
        Analytical_solution[i] = ((exp(-16.0*Alpha*t))*sin(4.0*x[i])) + ((exp(-Alpha*t))*sin(x[i]));
    }
}

real avg_error_vector(vector<real> &sol1,vector<real> &sol2)
{
    // Checking if the vectors are of the same size
    assert(sol1.size() == sol2.size());

    // avg_error = Average error
    // Error = Sum of absolute errors at all griid points
    real avg_error,Error;

    // Intializing the Error to 0.0
    Error = 0.0;

    // This loop iterates over the entire loop
    for(int i = 0;i<sol1.size();i++)
    {
        // Adding up the absolute errors at each of the grid points
        Error = Error + (abs(sol1[i]-sol2[i]));
    }
    // Calculating the average absolute error
    avg_error = Error/sol1.size();

    return avg_error;
}

// Function to save a vector<double> to a file of given name
void Writing_to_file(vector<double> &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<scientific<<endl;
    }

    // Closing the file
    file.close();
    cout<<endl;
}

void Discretization_in_1d(vector<real> &x, real dx, real X0, real X)
{
    // This loop iterates over the entire vector (t)
    for (int i = 0; (X0 + (i*dx)) <= X ;i++)
    {
        // Assigning values to different elements of t
        x.push_back(X0 + (i*dx));
        //t[i] = T0 + (i*dt);
    }
}

void Jacobi_Solver(vector<vector<real>> &A, vector<real> &b, real tol, vector<real> &x)
{
    int n = A.size();
    vector<real> xnew(n, 0.0);

    real diff = tol + 1.0;

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

        diff = 0.0;
        for (int i = 0; i < n; ++i)
        {
            diff = diff + abs(xnew[i]-x[i]);
        }
        x = xnew;
    }
}

int main()
{
    vector<real> Value_of_N_s,Value_of_rd_s,Inputs;
    vector<string> Output_file_names;
    real dx,Xl,Xu,pi,Tolerance;

    // Reading the inputs from the input file
    pi = 4.0*atan(1.0);
    Xl = 0;
    Xu = (2)*pi;
    real Alpha = 1.5;
    real t_end = 0.4;
    Tolerance = 1e-4;


    int N = 128;

    vector<real> x;

    dx = (Xu-Xl)/(N-1.0);
    Discretization_in_1d(x,dx,Xl,Xu);
    real r_d = 0.5;
    real rows,cols;

    // Calculating the value of dt
    real dt;
    dt = r_d*dx*dx/Alpha;

    // Calculating the number of the time steps required
    int Number_of_Steps = t_end/(dt);

    rows = N-1;
    cols = N-1;
    vector<real> U_n;
    vector<real> U_n_1;
    vector<vector <real> > A(rows,vector<real> (cols,0.0));
    vector<real> b;

    // Adjustments for the last time step
    real dt_new = t_end - (dt*Number_of_Steps);
    real r_d_new = dt_new*Alpha/(dx*dx);

    // Initial Conditions
    for (int i = 0;i < x.size()-1;i++)
    {
        U_n.push_back(sin(4.0*x[i])+ sin(x[i]));
    }

    // Computing the RHS
    b = U_n;

    // Computing the A matrix
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

    // Time marching loop
    for (int t_steps = 0;t_steps < Number_of_Steps;t_steps ++)
    {
        vector<real> x_sol(b.size(),0.0);
        Jacobi_Solver(A, b,Tolerance,x_sol);
        b = x_sol;
    }

    // Calculating the A matrix (Adjustments for the last time step)
    vector<real> x_sol(b.size(),0.0);
    A[0][0] = (1.0+(2.0*r_d_new));
    A[0][1] = -r_d_new;
    A[0][N-2] = -r_d_new;
    for (int i = 1;i < rows-1;i++)
    {
        A[i][i] = (1.0+(2.0*r_d_new));
        A[i][i-1] = -r_d_new;
        A[i][i+1] = -r_d_new;
    }
    A[N-2][N-2] = (1.0+(2.0*r_d));
    A[N-2][N-2-1] = -r_d_new;
    A[N-2][0] = -r_d_new;

    // Solving for the last time step
    Jacobi_Solver(A, b,Tolerance,x_sol);
    b = x_sol;


    vector<real> U_sol;
    U_sol = b;

    // Ghost point approach is used
    U_sol.push_back(U_sol[0]);

    Writing_to_file(U_sol,"Question_5_Implicit_U_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

    vector<real> Analytical_solution(N,0.0);
    analytical_solution(x,t_end,Alpha,Analytical_solution);

    Writing_to_file(Analytical_solution,"Analytical_Solution_U_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");

    vector <real> average_error;
    average_error.push_back(avg_error_vector(U_sol,Analytical_solution));

    Writing_to_file(average_error,"Question_5_Implicit_Average_Error_U_n_t_"+to_string(t_end)+"_N_"+to_string(N)+"_rd_"+to_string(r_d)+"_.csv");
}
