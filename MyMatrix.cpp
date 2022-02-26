// MyMatrix.cpp

#include <iostream>
#include"Matrix.h"
using namespace std;

int main()
{
	/*Matrix<int> a(5,5,1);
	a(1, 2) = 5;
	a(2, 2) = 6;
	a(3, 2) = 7;
	a(4, 2) = 8;
	a(5, 2) = 9;
	a(5, 3) = 9;
	a.print();
	a.Delete_col(3);
	a.Delete_row(2);
	a.print();*/

	Matrix<int> a(2, 2, 1);
	a(2, 2) = 3;
	a(2, 1) = 7;
	a.print();
	cout << Cofactor(1, 2, a);
	//cout <<a.Det();

	/*Matrix<double> c(2, 2, 3.2);
	Matrix<double> d(2, 2, 2.1);
	d(1, 1) = 9;
	Matrix<double> f = d/2.0;
	f.print();*/
	/*Matrix<int> b=(a*3);
	b.print();
	(a + b).print();*/
}


