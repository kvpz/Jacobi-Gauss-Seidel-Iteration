#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <iomanip>
#include <cstdlib>
using namespace std;

typedef vector<double> doubleVec;

int n = 3; //number of rows

void printZeroes(int i)
{
	int temp = 0;
	while (temp < i)
	{
		cout << 0 << ' ';
		temp++;
	}
}

//Prints a three band matrix of size nxn
void printMatrix(const doubleVec& U, const doubleVec& L, const doubleVec& D, const doubleVec& b)
{
	for (int i = 0; i < n ; i++)
	{
		if (i > 1) //when to start printing 0s below diagonal
			printZeroes(i-1); //prints 0s before the three bands in A

		if (i == 0) //prints first row
		{
			cout << D[i] <<U[i]<<' ';			
		}
		else if (i == 1) //prints second row
		{
			cout << L[i-1] << D[i] << U[i]<<' ';
		}
		else if (i>1 && i < n - 1) //prints rows 3->n-1
		{
			cout << L[i-1] << D[i] << U[i]<<' ';
		}
		else if (i >= n - 1) //prints the last row n
		{
			cout << L[i - 1] << D[i] << ' ';
		}

		printZeroes(n - 2 - i); //prints 0s after the three bands in A
		cout << b[i] << endl; //prints b[i] after every row of A
	}
}

/*
void inputDLU(vector<double> &D, vector<double> &L, vector<double> &U, double d_val, double l_val, double u_val)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) //diagonal values
				matrix[i][j] = d_val;
			else if (j == i + 1) //upper diagonal values
				matrix[i][j] = u_val;
			else if (j == i - 1) //lower diagonal values
				matrix[i][j] = l_val;
			else //all other matrix values
				matrix[i][j] = 0;
		}//inner for

		if (i % 2 == 0)
			matrix[i][m - 1] = 1;
		else
			matrix[i][m - 1] = 0;
	}//outer for
}
*/

void outputSolution(const doubleVec& x)
{
	//output
	cout << "solution:" << endl;
	for (auto i = x.begin(); i != x.end(); i++)
	{
		//i is of type iterator 
		cout << "x" << i - x.begin() << ": " << *i << endl;
	}
}

void findAx(doubleVec& Ax,const doubleVec& x,const doubleVec& D,const doubleVec& L,const doubleVec& U)
{
	Ax = vector < double >(n, 0.0); //reinitialize to an empty vector
	//find the value of diagonally dominant Ax per every iterative of x
	for (auto itr = Ax.begin(); itr != Ax.end(); itr++)
	{
		//auto position = distance(itr, Ax.begin());
		auto position = itr - Ax.begin(); //returns position in Ax

		if (position == 0)
			Ax[position] += x[position] * D[position] + x[1] * U[position];
		else if (position > 0 && position < n - 1)//middle rows
			Ax[position] += x[position - 1] * L[position - 1] + x[position] * D[position] + x[position + 1] * U[position];
		else
			Ax[position] += x[position - 1] * L[position - 1] + x[position] * D[position];

		//	cout << "AX" << position<<": "<<*itr << endl;
	}
}

int main()
{
	
	//The vectors U,L and D are used to form a three band diagonal matrix
	vector<double> U(n-1,-1.0);
	vector<double> L(n-1,-1.0); 
	vector<double> D(n,4.0);
	vector<double> b(n);
	vector<double> xKnot(n,0.0);//will hold previous iteration approximation
	vector<double> x(n, 0.0);//holds current iteration approximation
	vector<double> Ax(n, 0.0); //Ax matrix initialized to 0
	double tolerance = 0.0000000001;
	cout << setprecision(10)<<fixed;
	//input b vector
	for (int i = 0; i < n; i++)
	{
		if (i < n / 2)
			b[i] = 1.0;
		else
			b[i] = 0.0;
	}
	printMatrix(U,L,D,b);

	int k = 0;
	double sum;
	int N = 100000;

	int tolCount = 0;
	while (k < N)
	{//n:# of rows
		

		for (int i = 0; i < n; i++) //i represents the x component being estimated
		{
			sum = 0.0; //summation for first row
			if (i == 0)// summation for first row
					sum += U[i] * xKnot[i+1]; 
			else if (i>0 && i < n - 1) //summation for central rows
				sum += L[i - 1] * xKnot[i-1] + U[i] * xKnot[i + 1];
			else if (i == n-1) //summation for last row
				sum += L[i - 1] * xKnot[i-1];

			//x[component] ith approximation 
			x[i] = (-1.0*sum + b[i]) / D[i];
	
				if (abs(Ax[i] - b[i]) < tolerance && abs(Ax[i]-b[i]) != 0)
				{	
					cout << "Ax[" << i << "] has reached tolerance at iteration " <<k<< endl;
					tolCount++;
				}

		} //for loop
		
		//All components of the iterated x have reached the tolerance difference with b
		if (tolCount == 3)
			break;

		//find Ax for a diagonally dominant matrix, A
		findAx(Ax, x, D, L, U);

		//xKnot will hold this iteration of x approximation for the next iteration
		for (int j = 0; j < n; j++)
			xKnot[j] = x[j]; 
		
		k++;//k is compared to max number of desired iterations (N)

	}//while


	outputSolution(x); //output solution x


	
	return 0;

}

/*
U and L store n-1 elements. D, x and b store n elements.

Iteration does not always converge.  A sufficient condition for iteration to Jacobi iteration to converge is that A is strictly diagonally  dominant (the magnitude of the diagonal entry in a row is larger than or equal to the sum of the magnitudes of all the other (non-diagonal) entries in that row).
*/
/*
sum = 0.0;
if (k == 0) //first row
sum += -1.0*x[k] * U[i];
else if (k == n - 1) //last row
sum += -1.0*L[i - 1] * x[k];
else//rows in between
{
	sum += -1.0* L[i - 1] * x[k];
	sum += -1.0*x[k] * U[i];
}

x[k] = (sum + b[k]) / D[k];
cout << "x(" << k << "): " << x[k] << endl;
k++;
*/
