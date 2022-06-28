// MyMatrix.cpp

#include <iostream>
#include"Matrix.h"
#include<string>
using namespace std;


int main()
{
	Matrix<double> a;
	a.Read_file("D:/mine4ever/1.txt");
	a.print();
	//cout<<a.isSymmetric();
	cout << a.isPositiveDefinite();
	/*Matrix<double> b = a.Gaussian_elimination();
	b.print();*/

	/*Matrix<double> c = a.Rref();
	c.print();
	cout << c.Rank();*/

	/*Matrix<double> c(2, 2, 3.2);
	Matrix<double> d(2, 2, 2.1);
	d(1, 1) = 9;
	Matrix<double> f = d/2.0;
	f.print();*/
	/*Matrix<int> b=(a*3);
	b.print();
	(a + b).print();*/

}


