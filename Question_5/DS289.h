#ifndef DS289

#define DS289
#include<algorithm>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>
#include<cassert>

using namespace std;

// Using the double as a real number
typedef double real;

// t_interval(vector<real> &t, real dt, real T0, real T) populates the vector t, starting at T0 and ending at T with step size dt
void t_interval(vector<real> &t, real dt, real T0, real T)
{
    // This loop iterates over the entire vector (t)
    for (int i = 0; (T0 + (i*dt)) <= T ;i++)
    {
        // Assigning values to different elements of t
        t[i] = T0 + (i*dt);
    }
}


// Discretize_1D(vector<real> &x, real dx, real X0, real X) populates the vector x, starting at X0 and ending at X with step size dx
void Discretize_1D(vector<real> &x, real dx, real X0, real X)
{
    // This loop iterates over the entire vector (t)
    for (int i = 0; (X0 + (i*dx)) <= X ;i++)
    {
        // Assigning values to different elements of t
        x.push_back(X0 + (i*dx));
        //t[i] = T0 + (i*dt);
    }
}

// This function computes the product of two given matrix
void matrix_multiply(vector<vector<real>> &A, vector<vector<double>> &B, vector<vector<double>> &C)
{
    int m = A.size();
    int n = B.size();
    int p = B[0].size();

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            for (int k = 0; k < n; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// This function calculates the error between the analytical solution and numerical solution
// void Error(vector<real> &Analytical_Solution,vector<real> &Numerical_Solution,vector<real> &Error)
void Error(vector<real> &Analytical_T,vector<real> &Numerical_Solution,vector<real> &Error_T)
{
    // Just checking if the size of the analytical solution and numerical solution vectors are same or not
    // assert(Analytical_T.size() == Numerical_Solution.size());

    for(int i = 0;i < Analytical_T.size();i++)
    {
        // Appending the Error vector
        Error_T.push_back(abs(Analytical_T[i] - Numerical_Solution[i]));
    }
}

// adam(vector<real> &sol, vector<real> &t, function<double(real, real)> func) solves 1st order ODE
// sol = solution vector
// t = vector corresponding to the grid points
// func = function pointer to the function evaluating y'
void adam(vector<real> &sol, vector<real> &t, function<double(real, real)> func)
{
    // dt = Step size
    real dt;

    // This for loop iterates over the vector corresponding to the grid points (t)
    for(int i = 2; i<t.size();i++)
    {
        // Evaluating the step size
        // Generalized for non-uniform grid size
        dt = (t[i] - t[i-1]);

        // Scheme of the second-order Adams-Bashforth method
        // y_(i) = y(i-1) + (3*dt*(y'(i-1))/2) - (dt*(y'(i-2))/2)
        sol[i] = sol[i-1] + (3.0*dt*func(sol[i-1],t[i-1])/2.0) - (dt*func(sol[i-2],t[i-2])/2.0);
    }
}

// euler_explicit(vector<real> &sol, vector<real> &t, function<double(real, real)> func) solves 1st order ODE
// sol = solution vector
// t = vector corresponding to the grid points
// func = function pointer to the function evaluating y'
void euler_explicit(vector<real> &sol, vector<real> &t, function<double(real, real)> func)
{
    // This for loop iterates over the vector corresponding to the grid points (t)
    for (int i = 1;i < t.size();i++)
    {
        // Scheme of the explicit euler method
        // y_(i) = y(i-1) + (((t(i) - t(i-1))*y'(i-1))
        sol[i] = sol[i-1] + ((t[i] - t[i-1])*func(sol[i-1],t[i-1]));
    }
}

// euler_explicit(real sol,int i, vector<real> &t, function<double(real, real)> func) works same as
// euler_explicit(vector<real> &sol, vector<real> &t, function<double(real, real)> func) but here we are just returning the
// the value at the next step
real euler_explicit(real sol,int i, vector<real> &t, function<double(real, real)> func)
{
    // Scheme of the explicit euler method
    // y_(i) = y(i-1) + (((t(i) - t(i-1))*y'(i-1))
    return sol + ((t[i] - t[i-1])*func(sol,t[i-1]));
}

// RK4(vector<real> &sol, vector<real> &t, function<double(real, real)> func) solves 1st order ODE using RK4 method
// sol = solution vector
// t = vector corresponding to the grid points
// func = function pointer to the function evaluating y'
void RK4(vector<real> &sol, vector<real> &t, function<double(real, real)> func)
{
    // dt = step size
    double dt;

    // K1,K2,K3 and K4 are the "representative" slopes
    real K1,K2,K3,K4;

    // Iterates over the entire grid point vectors
    for(int i = 1; i<t.size();i++)
    {
        // dt = step size
        dt = (t[i] - t[i-1]);

        // K1 = dt*(y'(y(i-1),t(i-1)))
        K1 = dt*func(sol[i-1],t[i-1]);

        // K2 = dt*(y'(y(i-1) + (dt*K1/2),t(i-1) + (dt/2)))
        K2 = dt*func(sol[i-1] + (K1*(1.0/2.0)),t[i-1] + (dt*(1.0/2.0)));

        // K3 = dt*(y'(y(i-1) + (dt*K2/2),t(i-1) + (dt/2)))
        K3 = dt*func(sol[i-1] + (K2*(1.0/2.0)),t[i-1] + (dt*(1.0/2.0)));

        // K4 = dt*(y'(y(i-1) + (dt*K3),t(i-1) + (dt)))
        K4 = dt*func(sol[i-1] + K3,t[i-1] + (dt));

        // y(i) = y(i-1) + ((1.0/6.0)*(K1 + (2*K2) + (2*K3) + K4))
        sol[i] = sol[i-1] + (K1/6.0)+ (K2/3.0)+ (K3/3.0)+ (K4/6.0);
    }
}

// RK4(vector<real> &sol1, vector<real> &sol2, vector<real> &t, function<real(real)> func1, function<real(real, real, real, real, real)> func2, real c, real m, real k)
// solves 2nd order ODE using RK4 method as a system of 2 first order ODE
// sol1 and sol2 are solution vectors
// t = vector corresponding to the grid points
// func1 and func2 are function pointers to the functions evaluating the first order derivatives of the 2 first order ODE
void RK4(vector<real> &sol1, vector<real> &sol2, vector<real> &t, function<real(real)> func1, function<real(real, real, real, real, real)> func2, real c, real m, real k)
{
    // dt = step size
    double dt;

    // K11,K21,K31,K41,K12,K22,K32 and K42 are the "representative" slopes
    real K11,K21,K31,K41,K12,K22,K32,K42;

    // Iterates over the entire grid point vectors
    for(int i = 1; i<t.size();i++)
    {
        // dt = step size
        dt = (t[i] - t[i-1]);

        // Here in this particular case y1'(y1(i-1),y2(i-1),t(i-1)) = y1'(y2(i-1))

        // K11 = dt*(y1'(y1(i-1),y2(i-1),t(i-1)))
        K11 = dt*func1(sol2[i-1]);

        // We are passing the the values of c,m,k just to genneralize the function
        // K12 = dt*(y2'(y1(i-1),y2(i-1),t(i-1)))
        K12 = dt*func2(sol1[i-1],sol2[i-1],c,m,k);

        // K21 = dt*(y1'(y1(i-1) + (K11/2), y2(i-1) + (K12/2), t(i-1)+(dt/2)))
        K21 = dt*func1(sol2[i-1] + (K12*(1.0/2.0)));

        // K22 = dt*(y2'(y1(i-1) + (K11/2), y2(i-1) + (K12/2), t(i-1)+(dt/2)))
        K22 = dt*func2(sol1[i-1] + (K11*(1.0/2.0)), sol2[i-1] + (K12*(1.0/2.0)), c, m, k);

        // K31 = dt*(y1'(y1(i-1) + (K21/2), y2(i-1) + (K22/2), t(i-1)+(dt/2)))
        K31 = dt*func1(sol2[i-1] + (K22*(1.0/2.0)));

        // K32 = dt*(y2'(y1(i-1) + (K21/2), y2(i-1) + (K22/2), t(i-1)+(dt/2)))
        K32 = dt*func2(sol1[i-1] + (K21*(1.0/2.0)), sol2[i-1] + (K22*(1.0/2.0)), c, m, k);

        // K41 = dt*(y1'(y1(i-1) + (K31), y2(i-1) + (K32), t(i-1)+(dt)))
        K41 = dt*func1(sol2[i-1] + K32);

        // K42 = dt*(y2'(y1(i-1) + (K31), y2(i-1) + (K32), t(i-1)+(dt)))
        K42 = dt*func2(sol1[i-1] + K31,sol2[i-1] + K32,c,m,k);

        // y1(i) = y1(i-1) + (K11/6)+ (K21/3)+ (K31/3)+ (K41/6)
        sol1[i] = sol1[i-1] + (K11/6.0)+ (K21/3.0)+ (K31/3.0)+ (K41/6.0);

        // y2(i) = y2(i-1) + (K12/6)+ (K22/3)+ (K32/3)+ (K42/6)
        sol2[i] = sol2[i-1] + (K12/6)+ (K22/3)+ (K32/3)+ (K42/6);
    }
}

// avg_error_vector(vector<real> &sol1,vector<real> &sol2) calculates the average absolute error
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

// max_error_vector(vector<real> &sol1,vector<real> &sol2) calculates the maximum absolute error
real max_error_vector(vector<real> &sol1,vector<real> &sol2)
{
    // Checking if the vectors are of the same size
    assert(sol1.size() == sol2.size());

    // Vector to store the absolute error at each of the grid points
    vector<real> Error;

    // This loop iterate over the elements of the vectors
    for(int i = 0;i<sol1.size();i++)
    {
        // Calculating and storing the absolute errors at each of the grid points
        Error.push_back(abs(sol1[i]-sol2[i]));
    }

    // Getting the address of the maximum element in the vector Error
    auto int_result = std::max_element(Error.begin(),Error.end());

    // return the value of the maximum element in the vector Error
    return *int_result;
}

// Thomas_Algorithm_Matix_Solver(vector<vector <real> > &A, vector<real> &B, vector<real> &X)
// solves a system of equations: AX = B, where A is a Tridiagonal matrix
void Thomas_Algorithm_Matix_Solver(vector<vector <real> > &A, vector<real> &B, vector<real> &X)
{
    // Forward Sweep
    for (int row = 1;row < B.size(); row++)
    {
        A[row][row] = A[row][row] - (A[row-1][row]*(A[row][row-1]/A[row-1][row-1]));
        B[row] = B[row] - (B[row-1]*(A[row][row-1]/A[row-1][row-1]));
        A[row][row-1] = 0.0;
    }

    // Backward Sweep
    for (int row = B.size()-1; row >= 1; row--)
    {
        B[row-1] = B[row-1] - (B[row]*(A[row-1][row]/A[row][row]));
        A[row-1][row] = 0.0;
    }

    for (int row = 0;row < B.size(); row++)
    {
        X[row] = B[row]/A[row][row];
    }
}

// Function to save a vector<double> to a file of given name
void write_to_file(vector<double> &u, string str)
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

    cout<<"Output file "<< str <<" is opened."<<endl;

    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<scientific<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// Function to save a vector<string> to a file of given name
void write_to_file(vector<string> &u, string str)
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

    cout<<"Output file "<< str <<" is opened."<<endl;

    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// Function to save a vector<vector <real> > (MATRIX) to a file of given name
void write_to_file(vector<vector <real> > &u, string str)
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

    cout<<"Output file "<< str <<" is opened."<<endl;

    int rows = u.size();
    int cols = u[0].size();

    // Writing the vector values to the file in scientific notation

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols-1; col++)
        {
            file<<u[row][col]<<scientific<<",";
        }
        file<<u[row][cols-1]<<scientific<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

vector<real> input_parameters(string file_name)
{
    // Vector to read the input parameters from a input file
    vector<real> input;
    string item_name;
    int i = 0;

    // Open the file
    ifstream file;
    file.open(file_name);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: Input file could not be opened" << endl;
        exit(1);
    }

    cout<<"Input file "<<file_name<<" is opened."<<endl;

    string line;
    while (getline(file, line))
    {
        // Classifying a string as double if the first character is a numeric
        if(line[0]!= '/')
        {
            // To ensure that we are not reading white space characters
            if(isdigit(line[0]))
            {
                input.push_back(stod(line));
            }
        }
    }

    // Closing the input file
    file.close();

    cout<<"Input file "<<file_name<<" is closed."<<endl;
    cout<<endl;

    return input;
}

#endif

