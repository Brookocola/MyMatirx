// MyMatrix.cpp

#include <iostream>
#include"Matrix.h"
using namespace std;

int main()
{
	Matrix<int> a(3,4);
	a(1, 2) = 5;
	a.print();
	Matrix<int> b=(a-21);
	b.print();
	(a + b).print();
}


