// MyMatrix.cpp

#include <iostream>
#include"Matrix.h"
using namespace std;

int main()
{
	Matrix<double> a(4,5,1);
	a(1, 2) = 2;
	a(1, 3) = 3;
	a(1, 4) = 4;
	a(2, 1) = 5;
	a(2, 2) = 6;
	a(2, 3) = 7;
	a(2, 4) = 8;
	a(3, 1) = 7;
	a(3, 2) = 7;
	a(3, 3) = 8;
	a(3, 4) = 9;
	a(4, 1) = 2;
	a(4, 2) = 4;
	a(4, 3) = 5;
	a.print();
	/*a.Company().print();
	cout << a.Det()<<endl;*/
	a.Swap_col(1, 5);
	a.print();
	/*Matrix<double> c(2, 2, 3.2);
	Matrix<double> d(2, 2, 2.1);
	d(1, 1) = 9;
	Matrix<double> f = d/2.0;
	f.print();*/
	/*Matrix<int> b=(a*3);
	b.print();
	(a + b).print();*/
}


