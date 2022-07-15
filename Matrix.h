#pragma once
#include <iostream>
#include <fstream>
#include<string>
using namespace std;

//矩阵类
template<class T>
class Matrix
{
public:
	Matrix();//默认构造函数
	Matrix(const Matrix&M);//拷贝类构造函数
	Matrix(int row,int col);//构造函数_行列数
	Matrix(int row, int col, T value);//构造函数_行列数,值
	~Matrix();//析构函数

	bool Read_file(string filename);//从文件读取构造矩阵
	int Getrows() const {return rows;}//获取矩阵行数
	int Getcols() const {return cols;}//获取矩阵列数
	int Getsize() const {return size;}//获取矩阵大小
	void print() const;//输出矩阵
	Matrix<T> Transpose() const;//矩阵转置
	Matrix<double> Inverse() const;//矩阵求逆(求逆结果无法保证元素数据类型不变)
	T Det() const;//矩阵求行列式
	int Rank() const;//求矩阵的秩
	Matrix<double> Rref() const;//求简化的行阶梯形矩阵->高斯消元法
	void Delete_row(int row);//删除矩阵某行
	void Delete_col(int col);//删除矩阵某列
	void Swap_row(int r1, int r2);//交换矩阵行
	void Swap_col(int c1, int c2);//交换矩阵列
	Matrix<T> Delete_row_col(int row, int col) const;//删除矩阵某行某列
	Matrix<T> Company() const;//求伴随矩阵
	Matrix<double> Convert2double() const;//矩阵类型转换
	bool isPositiveDefinite() const;//矩阵是否正定
	bool isSymmetric() const;//矩阵是否对称
	bool equal(Matrix<T> M) const;//矩阵是否相等
	template<class Type> 
	friend Type Cofactor(int i, int j, Matrix<Type> M);//求代数余子式
	template<class Type>
	friend Type Det(Matrix<Type> M);//求矩阵行列式
	template<class Type>
	friend bool equal(Matrix<Type> M1, Matrix<Type> M2);//矩阵是否相等

	//操作符重载
	T& operator ()(int row, int col);//()操作符重载,用于访问矩阵元素
	Matrix<T>& operator =(const Matrix<T> &M);//操作符重载,同于赋值
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
	friend Matrix<double> operator/(Matrix<Type> M1, Matrix<Type> M2);//操作符重载,矩阵除法
	
	//TODO:
	//1.判断两个矩阵是否相似
	//2.矩阵相似对角化
	//3.矩阵求伪逆
	//4.矩阵特征值和特征向量
	//5.奇异值SVD分解

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
bool Matrix<T>::Read_file(string filename) {
	//读取txt文件构造矩阵
	ifstream infile(filename);
	string line;
	int row = 0;
	int col = 0;
	char pk;
	T tmp;
	if (infile.is_open() == false) {
		cout << "读取文件失败" << endl;
		return false;
	}
	else {
		//获取行列数
		while (true)
		{
			infile >> tmp;
			col++;
			pk = infile.peek();
			if (pk == '\n') break;
			if (infile.eof()) break;

		}
		infile.clear();//ifstream返回开头
		infile.seekg(ios::beg);
		while (true)
		{
			getline(infile, line);
			row++;
			if (infile.eof()) break;
		}
		infile.clear();//ifstream返回开头
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
		//只有方阵才能求逆
		cout << "只有方阵才能求逆" << endl;
	}
	else {
		//逆矩阵等于伴随矩阵除以行列式
		double det =(double)this->Det();
		Matrix<T> acc_T = this->Company();
		Matrix<double> acc_double = acc_T.Convert2double();
		Matrix<double> inv = acc_double / det;
		//检查-0
		for (int i = 0; i < inv.size; i++) {
			if (abs(inv.data[i]) < 1e-6) {
				inv.data[i] = 0;
			}
		}
		return inv;
	}
}
template<class T>
bool Matrix<T>::isPositiveDefinite() const {
	if (this->cols != this->rows) {
		cout << "只有方阵才能求正定性" << endl;
	}
	else if(!this->isSymmetric()) {
		//正定矩阵必须为对称矩阵
		//cout << "正定矩阵必须为对称矩阵";
		return false;
	}
	else {
		bool flag = false;
		int n = this->cols;
		Matrix<T> tmp = *this;
		while (n>=1) {
			if (tmp.Det() > 0) {
				tmp.Delete_row_col(n, n);
				n--;
			}
			else
				break;
		}
		if (n == 0)
			flag = true;
		return flag;
	}
}
template<class T>
bool Matrix<T>::isSymmetric() const {
	Matrix<T> m_transpose = this->Transpose();
	if (this->equal(m_transpose)) {
		return true;
	}
	else
		return false;
}
template<class T>
bool Matrix<T>::equal(Matrix<T> M) const {
	bool flag=false;
	int i;
	if (this->rows == M.rows && this->cols == M.cols) {
		for (i = 0; i < M.size; i++) {
			if (this->data[i] != M.data[i])
				break;
		}
		if (i == M.size)
			flag = true;
	}
	return flag;
}
template<class T>
T Matrix<T>::Det() const {
	if (this->size == 1) {
		return this->data[0];
	}
	else {
		T det = 0;
		//按第一行展开
		for (int j = 0; j < this->cols; j++) {
			Matrix<T> tmp = this->Delete_row_col(1, j + 1);
			det +=this->data[0 * this->cols + j] * tmp.Det()*(j % 2 ? -1 : 1);
		}
		return det;
	}
}
template<class T>
int Matrix<T>::Rank() const {
	//矩阵求秩
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
	//求行最简阶梯型矩阵
	Matrix<double> M_double = this->Convert2double();
	bool ifbreak=false;
	for (int i = 1,j = 1; i <= rows&&j <= cols; i++, j++) {
		int flag = i;
		while (abs(M_double(i,j))<1e-6)
		{
			if (j == cols&&flag==rows) {
				//没有下一个主元了
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
			//寻找主元
			if (abs(M_double(i, j)) >= 1e-6) {
				for (int k = i - 1; k >= 1; k--) {
					double coef = M_double(k, j) / M_double(i, j);
					//row k 
					for (int m = 1; m <= this->cols; m++) {
						M_double(k, m) = M_double(k, m) - coef * M_double(i, m);
					}
				}
				if (M_double(i, j) != 1) {
					//主元单位化
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
		return;//超过矩阵范围
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
		return;//超过矩阵范围
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
	//矩阵交换行
	if (r1 > this->rows || r2 > this->rows||r1<1||r2<1) {
		//超出矩阵范围
		cout << "行超出矩阵范围,无法交换" << endl;
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
	//矩阵交换列
	if (c1 > this->cols || c2 > this->cols||c1<1||c2<1) {
		//超出矩阵范围
		cout << "列超出矩阵范围,无法交换" << endl;
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
		//按第一行展开
		for (int j = 0; j < M.cols; j++) {
			Matrix<Type> tmp=M.Delete_row_col(1, j + 1);
			det += M.data[0 * M.cols + j]*Det(tmp)*(j % 2 ? -1 : 1);
		}
		return det;
	}
}
//操作符重载函数
template<class T>
T& Matrix<T>::operator()(int row, int col) {
	row--;//矩阵索引从1开始
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
		//矩阵大小不同无法相加
		cout << "矩阵大小不同无法相加" << endl;
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
		cout << "矩阵大小不同无法相减" << endl;
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
		//矩阵无法相乘
		cout << "M1的列数不等于M2的行数,无法相乘" << endl;
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