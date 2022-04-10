#pragma once
#include <iostream>
#include <fstream>
#include<string>
using namespace std;

//������
template<class T>
class Matrix
{
public:
	Matrix();//Ĭ�Ϲ��캯��
	Matrix(const Matrix&M);//�����๹�캯��
	Matrix(int row,int col);//���캯��_������
	Matrix(int row, int col, T value);//���캯��_������,ֵ
	~Matrix();//��������

	bool Read_file(string filename);//���ļ���ȡ�������
	int Getrows() const {return rows;}//��ȡ��������
	int Getcols() const {return cols;}//��ȡ��������
	int Getsize() const {return size;}//��ȡ�����С
	void print() const;//�������
	Matrix<T> Transpose() const;//����ת��
	Matrix<double> Inverse() const;//��������(�������޷���֤Ԫ���������Ͳ���)
	T Det() const;//����������ʽ
	int Rank() const;//��������
	Matrix<double> Rref() const;//��򻯵��н����ξ���->��˹��Ԫ��
	void Delete_row(int row);//ɾ������ĳ��
	void Delete_col(int col);//ɾ������ĳ��
	void Swap_row(int r1, int r2);//����������
	void Swap_col(int c1, int c2);//����������
	Matrix<T> Delete_row_col(int row, int col) const;//ɾ������ĳ��ĳ��
	Matrix<T> Company() const;//��������
	Matrix<double> Convert2double() const;//��������ת��
	template<class Type> 
	friend Type Cofactor(int i, int j, Matrix<Type> M);//���������ʽ
	template<class Type>
	friend Type Det(Matrix<Type> M);//���������ʽ

	//����������
	T& operator ()(int row, int col);//()����������,���ڷ��ʾ���Ԫ��
	Matrix<T>& operator =(const Matrix<T> &M);//����������,ͬ�ڸ�ֵ
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
	friend Matrix<double> operator/(Matrix<Type> M1, Matrix<Type> M2);//����������,�������

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
bool Matrix<T>::Read_file(string filename) {
	//��ȡtxt�ļ��������
	ifstream infile(filename);
	string line;
	int row = 0;
	int col = 0;
	char pk;
	T tmp;
	if (infile.is_open() == false) {
		cout << "��ȡ�ļ�ʧ��" << endl;
		return false;
	}
	else {
		//��ȡ������
		while (true)
		{
			infile >> tmp;
			col++;
			pk = infile.peek();
			if (pk == '\n')
				break;
		}
		infile.clear();//ifstream���ؿ�ͷ
		infile.seekg(ios::beg);
		while (true)
		{
			getline(infile, line);
			row++;
			if (infile.eof()) break;
		}
		infile.clear();//ifstream���ؿ�ͷ
		infile.seekg(ios::beg);
		this->rows = row;
		this->cols = col;
		size = rows * cols;
		data = new T[size];
		for (int i = 0; i < size; i++) {
			infile >> data[i];
		}
		return true;
	}
}
template<class T>
void Matrix<T>::print() const {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			cout << data[i*cols + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
template<class T>
Matrix<T> Matrix<T>::Transpose() const {
	Matrix<T> tmp(cols,rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			tmp.data[j*this->rows + i] = data[i*this->cols + j];
		}
	}
	return tmp;
}
template<class T>
Matrix<double> Matrix<T>::Inverse() const {
	if (this->cols != this->rows) {
		//ֻ�з����������
		cout << "ֻ�з����������" << endl;
	}
	else {
		double det =(double)this->Det();
		Matrix<T> acc_T = this->Company();
		Matrix<double> acc_double = acc_T.Convert2double();
		Matrix<double> inv = acc_double / det;
		return inv;
	}
}
template<class T>
T Matrix<T>::Det() const {
	if (this->size == 1) {
		return this->data[0];
	}
	else {
		T det = 0;
		//����һ��չ��
		for (int j = 0; j < this->cols; j++) {
			Matrix<T> tmp = this->Delete_row_col(1, j + 1);
			det +=this->data[0 * this->cols + j] * tmp.Det()*(j % 2 ? -1 : 1);
		}
		return det;
	}
}
template<class T>
int Matrix<T>::Rank() const {
	//��������
	int rank = 0;
	Matrix<double> M = this->Rref();
	for (int i = 1, j = 1; i <= M.rows&&j <= M.cols;) {
		if (M(i, j) != 0) {
			rank++;
			i++;
		}
		j++;
	}
	return rank;
}
template<class T>
Matrix<double> Matrix<T>::Rref() const {
	//�����������;���
	Matrix<double> M_double = this->Convert2double();
	bool ifbreak=false;
	for (int i = 1,j = 1; i <= rows&&j <= cols; i++, j++) {
		int flag = i;
		while (abs(M_double(i,j))<1e-6)
		{
			if (j == cols&&flag==rows) {
				//û����һ����Ԫ��
				ifbreak=true;
				break;
			}
			if (flag == rows) {
				j++;
				flag=i;
			}
			else if (abs(M_double(flag + 1, j)) >= 1e-6) {
				M_double.Swap_row(i, flag + 1);
			}
			else {
				flag++;
			}
		}
		if (ifbreak) {
			break;
		}
		for (int k = i + 1; k <= rows; k++) {
			double coef = M_double(k, j) / M_double(i, j);
			for (int m = 1; m <= cols; m++) {
				M_double(k, m) = M_double(k, m) - coef * M_double(i, m);
			}
		}
	}
	for (int i = this->rows; i >= 1; i--) {
		for (int j = 1; j <= this->cols;j++) {
			//Ѱ����Ԫ
			if (abs(M_double(i, j)) >= 1e-6) {
				for (int k = i - 1; k >= 1; k--) {
					double coef = M_double(k, j) / M_double(i, j);
					//row k 
					for (int m = 1; m <= this->cols; m++) {
						M_double(k, m) = M_double(k, m) - coef * M_double(i, m);
					}
				}
				if (M_double(i, j) != 1) {
					//��Ԫ��λ��
					double pivot = M_double(i, j);
					for (int n = 1; n <= this->cols; n++) {
						M_double(i, n) = M_double(i, n) / pivot;
					}
				}
				break;
			}
		}
	}
	for (int i = 0; i < size; i++) {
		if (abs(M_double.data[i]) < 1e-6) {
			M_double.data[i] = 0;
		}
	}
	return M_double;
}
template<class T>
void Matrix<T>::Delete_row(int row) {
	if (row<1 || row>this->rows) {
		return;//��������Χ
	}
	else {
		Matrix<T> tmp=*this;
		delete[] data;
		this->rows = tmp.rows - 1;
		this->size = this->rows*this->cols;
		data = new T[this->size];
		row--;
		for (int i = 0; i < this->size; i++) {
			if (i < row * this->cols) {
				this->data[i] = tmp.data[i];
			}
			else {
				this->data[i] = tmp.data[i + this->cols];
			}
		}
	}
}
template<class T>
void Matrix<T>::Delete_col(int col) {
	if (col<1 || col>this->cols) {
		return;//��������Χ
	}
	else {
		Matrix<T> tmp = *this;
		delete[] data;
		this->cols = tmp.cols - 1;
		this->size = this->rows*this->cols;
		data = new T[this->size];
		col--;
		int flag = 0;
		for (int i = 0; i < this->size; i++) {
			if (col == this->cols) {
				if (i != 0 && i%this->cols == 0) {
					flag++;
				}
			}
			else {
				if (i%this->cols == col) {
					flag++;
				}
			}
			
			this->data[i] = tmp.data[i + flag];
		}
	}
}
template<class T>
void Matrix<T>::Swap_row(int r1, int r2) {
	//���󽻻���
	if (r1 > this->rows || r2 > this->rows||r1<1||r2<1) {
		//��������Χ
		cout << "�г�������Χ,�޷�����" << endl;
	}
	else
	{
		T* row1 = new T[this->cols];
		for (int j = 0; j < this->cols; j++) {
			row1[j] = this->data[(r1 - 1)*this->cols + j];
		}
		for (int j = 0; j < this->cols; j++) {
			this->data[(r1 - 1)*this->cols + j] = this->data[(r2 - 1)*this->cols + j];
			this->data[(r2 - 1)*this->cols + j] = row1[j];
		}
		delete[] row1;
	}
}
template<class T>
void Matrix<T>::Swap_col(int c1, int c2) {
	//���󽻻���
	if (c1 > this->cols || c2 > this->cols||c1<1||c2<1) {
		//��������Χ
		cout << "�г�������Χ,�޷�����" << endl;
	}
	else
	{
		T* col1 = new T[this->rows];
		for (int i = 0; i < this->rows; i++) {
			col1[i] = this->data[i*this->cols + c1-1];
		}
		for (int i = 0; i < this->rows; i++) {
			this->data[i*this->cols + c1 - 1] = this->data[i*this->cols + c2 - 1];
			this->data[i*this->cols + c2 - 1] = col1[i];
		}
		delete[] col1;
	}
}
template<class T>
Matrix<T> Matrix<T>::Delete_row_col(int row, int col) const {
	Matrix<T> tmp = *this;
	tmp.Delete_row(row);
	tmp.Delete_col(col);
	return tmp;
}
template<class T>
Matrix<T> Matrix<T>::Company() const {
	Matrix<T> acc = *this;
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			acc.data[i*this->cols + j] = Cofactor(i + 1, j + 1,*this);
		}
	}
	return acc.Transpose();
}
template<class T>
Matrix<double> Matrix<T>::Convert2double() const {
	Matrix<double> tmp(this->rows, this->cols,0);
	for(int i = 0; i<this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			tmp(i + 1, j + 1) = (double)this->data[i*this->cols + j];
		}
	}
	return tmp;
}
template<class Type>
Type Cofactor(int i, int j, Matrix<Type> M) {
	Matrix<Type> tmp = M;
	tmp.Delete_row(i);
	tmp.Delete_col(j);
	return ((i+j) % 2 ? -1 : 1)*Det(tmp);
}
template<class Type>
Type Det(Matrix<Type> M) {
	if (M.size== 1) {
		return M.data[0];
	}
	else {
		Type det = 0;
		//����һ��չ��
		for (int j = 0; j < M.cols; j++) {
			Matrix<Type> tmp=M.Delete_row_col(1, j + 1);
			det += M.data[0 * M.cols + j]*Det(tmp)*(j % 2 ? -1 : 1);
		}
		return det;
	}
}
//���������غ���
template<class T>
T& Matrix<T>::operator()(int row, int col) {
	row--;//����������1��ʼ
	col--;
	return data[row*cols + col];
}
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &M) {
	if (this == &M) {
		return *this;
	}
	else {
		delete[] data;
		this->rows = M.rows;
		this->cols = M.cols;
		this->size = M.size;
		data = new T[this->size];
		for (int i = 0; i < this->size; i++) {
			this->data[i] = M.data[i];
		}
		return *this;
	}
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
		cout << "�����С��ͬ�޷����" << endl;
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
		cout << "�����С��ͬ�޷����" << endl;
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
template<class Type>
Matrix<Type> operator*(Matrix<Type> M, Type value) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] *= value;
	}
	return tmp;
}
template<class Type>
Matrix<Type> operator*(Type value,Matrix<Type> M ) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] *= value;
	}
	return tmp;
}
template<class Type>
Matrix<Type> operator*(Matrix<Type> M1, Matrix<Type> M2) {
	if (M1.cols!=M2.rows) {
		//�����޷����
		cout << "M1������������M2������,�޷����" << endl;
	}
	else {
		int tmp_rows, tmp_cols;
		tmp_rows = M1.rows;
		tmp_cols = M2.cols;
		Matrix<Type> tmp(tmp_rows, tmp_cols);
		for (int i = 0; i < tmp_rows; i++) {
			for (int j = 0; j < tmp_cols; j++) {
				Type value = 0;
				for (int m = 0; m <M1.cols; m++) {
					value += M1.data[i*M1.cols + m] * M2.data[m*M2.cols + j];
				}
				tmp.data[i*tmp_cols + j] = value;
			}
		}
		return tmp;
	}
}
template<class Type>
Matrix<Type> operator/(Matrix<Type> M, Type value) {
	Matrix<Type> tmp = M;
	for (int i = 0; i < tmp.size; i++) {
		tmp.data[i] /= value;
	}
	return tmp;
}
template<class Type>
Matrix<double> operator/(Matrix<Type> M1, Matrix<Type> M2) {
	//M1/M2=M1*M2_inverse
	Matrix<double> M2_inverse = M2.Inverse();
	Matrix<double> tmp = M1 * M2_inverse;
	return tmp;
}