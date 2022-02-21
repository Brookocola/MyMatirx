#pragma once
#include <iostream>
using namespace std;

//矩阵类
template<class T>
class Matrix
{
public:
	Matrix();//默认构造函数
	Matrix(const Matrix&M);//拷贝类构造函数,类构造
	Matrix(int row,int col);//构造函数_行列数
	Matrix(int row, int col, T value);//构造函数_行列数,值
	~Matrix();//析构函数
	int Getrows() const {return rows;}//获取矩阵行数
	int Getcols() const {return cols;}//获取矩阵列数
	int Getsize() const {return size;}//获取矩阵大小
	void print() const;//输出矩阵
	Matrix<T> Transpose();//矩阵转置
	Matrix<double> Inverse();//矩阵求逆(求逆结果无法保证元素数据类型不变)
	T Det() const;//矩阵求行列式
	int Rank() const;//求矩阵的秩


	//操作符重载
	T& operator ()(int row, int col);//()操作符重载,用于访问矩阵元素
	void operator +=(T value);//操作符重载,矩阵数加
	void operator -=(T value);//操作符重载,矩阵数减
	void operator *=(T value);//操作符重载,矩阵数乘
	void operator /=(T value);//操作符重载,矩阵数除

	template<class Type>//友元函数声明时需要加template<class Type>
	friend Matrix<Type> operator+(Matrix<Type> M, Type value);//操作符重载,矩阵数加1
	template<class Type>
	friend Matrix<Type> operator+(Type value, Matrix<Type> M);//操作符重载,矩阵数加2
	template<class Type>
	friend Matrix<Type> operator+(Matrix<Type> M1, Matrix<Type> M2);//操作符重载,矩阵加法
	template<class Type>
	friend Matrix<Type> operator-(Matrix<Type> M, Type value);//操作符重载,矩阵数减
	template<class Type>
	friend Matrix<Type> operator-(Matrix<Type> M1, Matrix<Type> M2);//操作符重载,矩阵减法
	template<class Type>
	friend Matrix<Type> operator*(Matrix<Type> M, Type value);//操作符重载,矩阵数乘1
	template<class Type>
	friend Matrix<Type> operator*(Type value, Matrix<Type> M);//操作符重载,矩阵数乘2
	template<class Type>
	friend Matrix<Type> operator*(Matrix<Type> M1, Matrix<Type> M2);//操作符重载,矩阵乘法
	template<class Type>
	friend Matrix<Type> operator/(Matrix<Type> M, Type value);//操作符重载,矩阵数除
	template<class Type>
	friend Matrix<Type> operator/(Matrix<Type> M1, Matrix<Type> M2);//操作符重载,矩阵除法

private:
	int rows, cols;//矩阵行列数
	int size;//矩阵大小
	T* data;//存储矩阵元素
};


//函数实现
template<class T>
Matrix<T>::Matrix() {
	rows = 0;
	cols = 0;
	size = 0;
	data = NULL;
}
template<class T>
Matrix<T>::Matrix(const Matrix&M) {
	rows = M.rows;
	cols = M.cols;
	size = M.size;
	data = new T[size];
	for (int i = 0; i < size; i++) {
		data[i] = M.data[i];
	}
}
template<class T>
Matrix<T>::Matrix(int row,int col)
{
	rows = row;
	cols = col;
	size = rows*cols;
	data = new T[size];
	for (int i = 0; i < size; i++) {
		data[i] = 0;//全赋值为0
	}
}
template<class T>
Matrix<T>::Matrix(int row, int col,T value)
{
	rows = row;
	cols = col;
	size = rows * cols;
	data = new T[size];
	for (int i = 0; i < size; i++) {
		data[i] = value;//全赋值为value
	}
}
template<class T>
Matrix<T>::~Matrix()
{
	delete []data;
}
template<class T>
void Matrix<T>::print() const {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			cout << data[i*cols + j] << " ";
		}
		cout << endl;
	}
}
template<class T>
Matrix<T> Matrix<T>::Transpose() {
	Matrix<T> tmp(cols,rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			tmp.data[j*this->rows + i] = data[i*this->cols + j];
		}
	}
	return tmp;
}
template<class T>
T& Matrix<T>::operator()(int row, int col) {
	return data[row*cols + col];
}
template<class T>
void Matrix<T>::operator+=(T value) {
	for (int i = 0; i < size; i++) {
		data[i] += value;
	}
}
template<class T>
void Matrix<T>::operator-=(T value) {
	for (int i = 0; i < size; i++) {
		data[i] -= value;
	}
}
template<class T>
void Matrix<T>::operator*=(T value) {
	for (int i = 0; i < size; i++) {
		data[i] *= value;
	}
}
template<class T>
void Matrix<T>::operator/=(T value) {
	if (value == 0)
		return;
	for (int i = 0; i < size; i++) {
		data[i] /= value;
	}
}

template<class Type>
Matrix<Type> operator+(Matrix<Type> M, Type value) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] += value;
	}
	return tmp;
}
template<class Type>
Matrix<Type> operator+(Type value, Matrix<Type> M) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] += value;
	}
	return tmp;
}
template<class Type>
Matrix<Type> operator+(Matrix<Type> M1, Matrix<Type> M2) {
	if (M1.cols != M2.cols || M1.rows != M2.rows) {
		//矩阵大小不同无法相加
	}
	else
	{
		Matrix<Type> tmp = M1;
		for (int i = 0; i < tmp.size; i++) {
			tmp.data[i] = M1.data[i] + M2.data[i];
		}
		return tmp;
	}
}
template<class Type>
Matrix<Type> operator-(Matrix<Type> M, Type value) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] -= value;
	}
	return tmp;
}
template<class Type>
Matrix<Type> operator-(Matrix<Type> M1, Matrix<Type> M2) {
	if (M1.cols != M2.cols || M1.rows != M2.rows) {
		//矩阵大小不同无法相减
	}
	else
	{
		Matrix<Type> tmp = M1;
		for (int i = 0; i < tmp.size; i++) {
			tmp.data[i] = M1.data[i] - M2.data[i];
		}
		return tmp;
	}
}