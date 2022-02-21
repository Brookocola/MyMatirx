#pragma once
#include <iostream>
using namespace std;

//������
template<class T>
class Matrix
{
public:
	Matrix();//Ĭ�Ϲ��캯��
	Matrix(const Matrix&M);//�����๹�캯��,�๹��
	Matrix(int row,int col);//���캯��_������
	Matrix(int row, int col, T value);//���캯��_������,ֵ
	~Matrix();//��������
	int Getrows() const {return rows;}//��ȡ��������
	int Getcols() const {return cols;}//��ȡ��������
	int Getsize() const {return size;}//��ȡ�����С
	void print() const;//�������
	Matrix<T> Transpose();//����ת��
	Matrix<double> Inverse();//��������(�������޷���֤Ԫ���������Ͳ���)
	T Det() const;//����������ʽ
	int Rank() const;//��������


	//����������
	T& operator ()(int row, int col);//()����������,���ڷ��ʾ���Ԫ��
	void operator +=(T value);//����������,��������
	void operator -=(T value);//����������,��������
	void operator *=(T value);//����������,��������
	void operator /=(T value);//����������,��������

	template<class Type>//��Ԫ��������ʱ��Ҫ��template<class Type>
	friend Matrix<Type> operator+(Matrix<Type> M, Type value);//����������,��������1
	template<class Type>
	friend Matrix<Type> operator+(Type value, Matrix<Type> M);//����������,��������2
	template<class Type>
	friend Matrix<Type> operator+(Matrix<Type> M1, Matrix<Type> M2);//����������,����ӷ�
	template<class Type>
	friend Matrix<Type> operator-(Matrix<Type> M, Type value);//����������,��������
	template<class Type>
	friend Matrix<Type> operator-(Matrix<Type> M1, Matrix<Type> M2);//����������,�������
	template<class Type>
	friend Matrix<Type> operator*(Matrix<Type> M, Type value);//����������,��������1
	template<class Type>
	friend Matrix<Type> operator*(Type value, Matrix<Type> M);//����������,��������2
	template<class Type>
	friend Matrix<Type> operator*(Matrix<Type> M1, Matrix<Type> M2);//����������,����˷�
	template<class Type>
	friend Matrix<Type> operator/(Matrix<Type> M, Type value);//����������,��������
	template<class Type>
	friend Matrix<Type> operator/(Matrix<Type> M1, Matrix<Type> M2);//����������,�������

private:
	int rows, cols;//����������
	int size;//�����С
	T* data;//�洢����Ԫ��
};


//����ʵ��
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
		data[i] = 0;//ȫ��ֵΪ0
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
		data[i] = value;//ȫ��ֵΪvalue
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
		//�����С��ͬ�޷����
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
		//�����С��ͬ�޷����
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