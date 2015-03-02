#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>/*fabs()*/

using namespace std;


//global variables
const int n = 3; //rows
const int m = n+1; //columns

void inputMatrixRows(double(&mat)[n][m])
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cin >> mat[i][j];
		}
	}
}

void inputMatrixColumns(double(&mat)[n][m])
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cin >> mat[j][i];
		}
	}
}

//print nxm matrix
void printMatrix(double(&mat)[n][m])
{
	cout << setprecision(1) << fixed;
	for (int i = 0; i<n; i++)    //This loops on the rows
	{
		for (int j = 0; j<m; j++) //This loops on the columns
		{
			cout << mat[i][j] << "  "; //prints a row
		}
		cout << endl;
	}
}

//prints nxn matrix
void printMatrix(double(&mat)[n][n])
{
	cout << setprecision(1) << fixed;
	for (int i = 0; i<n; i++)    //This loops on the rows
	{
		for (int j = 0; j<m; j++) //This loops on the columns
		{
			cout << mat[i][j] << "  "; //prints a row
		}
		cout << endl;
	}
}

void swapRows(double (&mat)[n][m], int rowBelow, int rowOnTop)
{
	double temp;
	for (int c = 0; c < m; c++)
	{
		temp = mat[rowBelow][c];
		mat[rowBelow][c] = mat[rowOnTop][c];
		mat[rowOnTop][c] = temp;
	}
}

void extractAb(double (&mat)[n][m], double (&A)[n][n], double *b)
{
	//Extract A and b from augmented matrix
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) //A is nxn matrix
		{
			A[i][j] = mat[i][j];
		}
		b[i] = mat[i][m-1]; //gets values(b) in last column of augmented matrix
	}
}

void inputDLU(double matrix[n][m],double d_val, double l_val, double u_val)
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

		/*
		if (i % 2 == 0)
			matrix[i][m - 1] = 1;
		else
			matrix[i][m - 1] = 0;
		*/
		if (i < n / 2)
			matrix[i][m-1] = 1;
		else
			matrix[i][m-1] = 0;
	}//outer for
}

int main()
{   
	
	double matrix[n][m];
//	inputMatrixColumns(matrix);
//	printMatrix(matrix);
	// = { { 0, 1, 0, 1 }, { -1, 0, 0, 1 }, { 0, 0, 1, 1 } };
	/*{{ 4, 1, 2, 9 },
	 { 2, 4, -1, -5 },
	 { 1, 1, -3, -9 } };
	
	 = { { 4, -1, 1, 8 },
	    { 2, 5, 2, 3 },
	    { 1, 2, 4, 11 } };
	*/
	
	inputDLU(matrix,4, -1, -1); //matrix,D,L,U
	/*
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j) //diagonal values
				matrix[i][j] = 2;
			else if (j == i + 1 || j == i - 1) //off diagonal values
				matrix[i][j] = -1;
			else //all other matrix values
				matrix[i][j] = 0;
		}//inner for
		if (i % 2 == 0)
			matrix[i][m-1] = 1;
		else
			matrix[i][m-1] = 0;
	}//outer for
	*/
	printMatrix(matrix);
	
	/**************/
	//this is performed here because "matrix" is modified in the rest of the program
	double A[n][n];
	double b[n];
	extractAb(matrix, A, b); //A and b are to be used in 2nd subroutine at bottom
	/**********/
//	printMatrix(matrix);

	double largestCoefTemp; //largest pivot value for a swap
	double largestCoefficient;
	int largestRow=0; //collects value of row to swap 
	bool swap=false; //true if swapping condition met(found larger coefficient in column below pivot)
	double multiple; //=matrix[z][i]/matrix[i][i], z != 0 used during elimination

	/*STEP 1*/
	for (int i = 0; i < n-1; i++)    //This loops on the columns  
	{ 
		largestCoefTemp = 0; //stores largest column value to compare with others
		/*STEP 2*/													
		//search column for largest value in order to know what row to swap				
		for (int p=i; p < n; p++) //p= 0 : number of rows
		{
			//cout << "asdfa" << endl;
			if (fabs(matrix[p][i]) > fabs(matrix[i][i])) //find larger coefficient in the column than current pivot
			{
				//cout << "fabs(matrix[p][i])>fabs(matrix[i][i])" << endl;
				largestCoefTemp = matrix[p][i];
				//check if larger than the previous large value found
				if (fabs(matrix[p][i]) > largestCoefTemp)
				{
				//	cout << "fabs(matrix[p][i])>largestCoefTemp" << endl;
					largestCoefficient = matrix[p][i]; //collects largest column value to compare to others
					largestRow = p; //row with largest potential pivot
					swap = true; //must swap when a larger number than current pivot is found
					multiple = matrix[p][i] / matrix[i][i]; //where p>i
				}//inner if
			}//outer if
		}//for loop 2

		if (largestCoefficient == 0)
		{
			cout << "Pivot point is 0. End of algorithm" << endl;
			exit(0);
		}

		/*STEP 3*/
		//For increased numerical stability largest possible pivot is used.*
		if (swap == true && largestRow!=i) //don't swap if current pivot is largest value in column
		{
			swapRows(matrix, largestRow,i);
			swap = false; //reset
		}
		
		if (largestRow == i) //assign row immediately below current pivot if no larger coefficient in column was found
			largestRow = i + 1;
		
		/*STEP 4*/
		for (int z = largestRow; z < n; z++)
		{	/*STEP 5*/
			multiple = matrix[z][i] / matrix[i][i]; //used to multiply through row that will eliminate one below
			/*STEP 6*/
			for (int j = 0; j < m; j++) //This loops through a row
				matrix[z][j] = matrix[z][j] - multiple*matrix[i][j]; //elimination/ change of coefficients
		}//for loop 3
	}//outer for, for loop 1


	if (matrix[n - 1][n - 1] == 0)
	{
		cout << "No Solution" << endl;
		exit(0);
	}

	/*Finding solution X*/
	double x[m-1]; //x values to be found 
	double sum; //collects sum of product with existing x values and their coefficients
	x[n-1] = matrix[n - 1][m - 1] / matrix[n - 1][n - 1];//value of last x (x[n-1])
	for (int i = n - 2; i >= 0; i--)
	{ //backward substitution
		sum = 0;
		for (int j = i + 1; j < n; j++)
		{
			sum += matrix[i][j] * x[j]; //sum of rest of equation when other x are plugged in
		}
		x[i] = (matrix[i][m - 1] - sum) / matrix[i][i];
	}

	//output solution X
	cout << setprecision(10);
	cout << "X = (";
	for(int i = 0; i < n; i++)
	{
		cout << x[i];
		if(i<=n-2)
			cout<< ", ";
	}
	cout << ")"<<endl;


	/****2nd subroutine****/


	double r[n]; //r=b-Ax. Gets r for each row
	double Ax[n][n]; 
	double SUM=0;// [n] = { 0, 0, 0 }; //causes problems in arithmetic(+=) if not initialized to 0 before

	for (int i = 0; i < n; i++)//iterates through the rows
	{
		SUM = 0;
		for (int j = 0; j < n; j++)//iterates through each column for a row in order to get SUM
		{
			Ax[i][j] = A[i][j] * x[j];
		//	cout << "Ax: " << Ax[i][j] << endl;
			SUM += Ax[i][j];
			
		}
		r[i] = SUM - b[i];

		//to avoid the output of -0.0 (for neatness)
		cout << "r= "<<(r[i] <= 0.0 ? 0.0 : r[i]) << endl;
	}

	

	return 0;
}






/*
	NOTES
The largest pivot point is used because it results in most accurate numerical output
*/
