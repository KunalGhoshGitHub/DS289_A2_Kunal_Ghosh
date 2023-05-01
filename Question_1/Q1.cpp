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

extern "C" void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );

int main()
{
    real x_l,x_u,y_l,y_u,T_L,T_R,d_T_D,d_T_U_a;
    int Nx,Ny;

    // Vector to store input from the input file
    vector<real> Inputs;
    Inputs = input_parameters("Input.txt");

    // Computational domain:
    // x_l <= x <= x_u
    // y_l <= y <= y_u
    x_l = Inputs[0];
    x_u = Inputs[1];
    y_l = Inputs[2];
    y_u = Inputs[3];

    // Number of grid points:
    // Number of grid points along x-direction: (Nx)
    // Number of grid points along y-direction: (Ny)
    Nx =  Inputs[4];
    Ny = Inputs[5];

    // Boundary Conditions:
    // T_L = T (x_l, y)
    T_L = Inputs[6];
    // T_R = T (x_u, y)
    T_R = Inputs[7];
    // d_T_D = dT/dy (x, y_l) (This is a partial derivative with respect to y)
    d_T_D = Inputs[8];
    // d_T_U_a = dT /dy(x, y_u) = (T (x, y_u) − a), where a is (This is a partial derivative with respect to y)
    d_T_U_a = Inputs[9];

    real m,dx,dy;

    // Calculating the value of dx, dy and m

    dx = (x_u-x_l)/(Nx-1.0);
    dy = (y_u-y_l)/(Ny-1.0);
	m = (dx/dy)*(dx/dy);

    int size_of_A;
    // Calculating the size of A
    size_of_A = Nx*Ny*Nx*Ny;

    // Allocating the memory for A array
    real *A_array = (real *)malloc(size_of_A * sizeof(real));

    int size_of_b;

    // Calculating the size of b array
    size_of_b = Nx*Ny;

    // Allocating the memory for A array
    real *b_array = (real *)malloc(size_of_b * sizeof(real));
    // real b_array[size_of_b];

    int rows,cols;
    rows = Nx*Ny;
    cols = Nx*Ny;

    vector<vector <real> > A(rows,vector<real> (cols,0.0));
    vector<real> b(rows, 0.0);

    // Creating the A matrix  and setting up the boundary conditions
    for (int i = 0;i < Ny;i++)
    {
        A[i][i] = 1.0;
        b[i] = T_L;
    }


    for (int i = 1;i <= Ny;i++)
    {
        A[(Nx*Ny)-i][(Nx*Ny)-i] = 1.0;
        b[(Nx*Ny)-i] = T_R;
    }


    for (int j = 1;j < Nx-1;j++)
    {
        for (int  i = 1; i < Ny-1;i++)
        {
            A[(j*Ny)+i][(j*Ny)+i] =-2.0*(m+1.0);
            A[(j*Ny)+i][(j*Ny)+i-1] =m;
            A[(j*Ny)+i][(j*Ny)+i-Ny] =1.0;
            A[(j*Ny)+i][(j*Ny)+i+1] =m;
            A[(j*Ny)+i][(j*Ny)+i+Ny] =1.0;
        }
    }


    for (int j = 1;j < Nx-1;j++)
    {
        A[(j*Ny)][(j*Ny)] =-3;
        A[(j*Ny)][(j*Ny)+1] =4.0;
        A[(j*Ny)][(j*Ny)+2] =-1.0;
    }

        for (int j = 1;j < Nx-1;j++)
    {
        A[(j*Ny)+Ny-1][(j*Ny)+Ny-1] =(3.0 -(2.0*dy));
        A[(j*Ny)+Ny-1][(j*Ny)+Ny-1-1] =-4.0;
        A[(j*Ny)+Ny-1][(j*Ny)+Ny-1-2] =1;
        b[(j*Ny)+Ny-1] =(-d_T_U_a*2.0*dy);
    }

    int ctr = 0;
    for (int j = 0; j < cols;j++)
    {
        for (int i=0; i < rows;i++)
        {
            A_array[ctr++] =  A[i][j];
        }
    }

    ctr = 0;
    for (int i=0; i < rows;i++)
    {
        b_array[ctr++] = b[i];
    }

    // Seeting the parameter for Lapack dgesv_ function
	char trans = 'N';
	int dim = Nx*Ny;
	int nrhs = 1;
	int LDA =dim;
	int LDB = dim;
	int info;
	int n = dim;
	int lda = LDA;
	int ldb = LDB;

    int ipiv[Nx*Ny];

	dgesv_( &n, &nrhs, A_array, &lda, ipiv, b_array, &ldb, &info );

	cout<<"Output from the dgesv_ function: "<<"Info  = " <<info<<endl;
	cout<<endl;

	vector<real> sol;
	vector<vector <real> > Solution(Ny,vector<real> (Nx,0.0));
	ctr = 0;


	// Converting the solution to a Matrix
	for (int i = 0; i < Nx; i++)
	{

        for (int j = 0; j < Ny; j++)
        {
            Solution[j][i] = b_array[ctr++];
        }
	}

	// Writing the solution to a file
	write_to_file(Solution,"Question_1_Solution.csv");

}
