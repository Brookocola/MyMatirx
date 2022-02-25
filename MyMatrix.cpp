// MyMatrix.cpp

#include <iostream>
#include"Matrix.h"
using namespace std;

int main()
{
	Matrix<int> a(3,4,1);
	a(1, 2) = 5;
	a.print();
	Matrix<double> c(2, 2, 3.2);
	Matrix<double> d(2, 2, 2.1);
	d(1, 1) = 9;
	Matrix<double> f = d/2.0;
	f.print();
	/*Matrix<int> b=(a*3);
	b.print();
	(a + b).print();*/
}


