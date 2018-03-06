#include"matrix.hpp"

bool canAdd(const Matrix & left, const Matrix & right) {
	return left.h_ == right.h_ && left.w_ == right.w_;
}

bool canMul(const Matrix & left, const Matrix & right) {
	return right.h == left.w;
}

Matrix::Matrix(float* in, unsigned char w, unsigned char h) : h_(h), w_(w) {
	for (unsigned char i = 0; i < w*h; ++i) {
		elements[i] = in[i];
	}
}

Matrix::Matrix(unsigned char w, unsigned char h) : h_(h), w_(w) {
	if (h*w == 0) {
		elements = nullptr;
	} else {
		elements = new float[h*w];
	}
};

Matrix::Matrix() : Matrix(0, 0) {};

Matrix::Matrix(const Matrix& other) {
	h_ = other.h_;
	w_ = other.w_;
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] = other.elements[i];
	}
}

Matrix & Matrix::operator=(Matrix other) {
	swap(*this, other);
}

Matrix::~Matrix() {
	delete[] elements;
}

float & Matrix::operator[](loc coord) {
	return elements[coord.x, coord.y];
}

float Matrix::operator[](loc coord) const {
	return elements[coord.x, coord.y];
}

Matrix Matrix::operator+(const Matrix& other) const {
	Matrix result(w_, h_);
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_;++i) {
			result.elements[i] = other.elements[i]+elements[i];
		}
	} else {
		result.NaNify();
	}
}

Matrix & Matrix::operator+=(const Matrix& other) {
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			elements[i] = other.elements[i] + elements[i];
		}
	}
	else {
		NaNify();
	}
	return *this;
}

Matrix Matrix::operator-(const Matrix& other) const {
	Matrix result(w_, h_);
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			result.elements[i] = other.elements[i] - elements[i];
		}
	}
	else {
		result.NaNify();
	}
}

Matrix & Matrix::operator-=(const Matrix& other) {
	if (canAdd(*this, other)) {
		for (unsigned char i = 0; i < h_*w_; ++i) {
			elements[i] = other.elements[i] - elements[i];
		}
	}
	else {
		NaNify();
	}
	return *this;
}

Matrix Matrix::operator*(const float &coeff) const {
	Matrix result(w_, h_);
	for (unsigned char i = 0; i < h_*w_; ++i) {
		result.elements[i] = elements[i] * coeff;
	}
}

Matrix & Matrix::operator*=(const float& coeff) {
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] *= coeff;
	}
}

Matrix Matrix::operator*(const Matrix & other) const {
<<<<<<< HEAD
	// toDo: max it "work" for two row vectors or two colum vectors, or two vectors in the wrong order
=======
	if(!canMul(*this,other)){
		if (/*column, column*/){
			return transp(*this)*other;
		} else if (/*row, row*/){
			return (*this)*transp(other);
		} else if (/*column, row*/){
			return transp(*this)*transp(other);
		} else {
			Matrix result(1,1);
			result.NaNify;
			return result;
		}
	}
>>>>>>> 325ba1e071fe52eb2a2f03b9d75cf310936d23cf
	// from https://github.com/eecharlie/MatrixMath/blob/master/MatrixMath.cpp
	Matrix result(h_, other.w_);
	unsigned char i, j, k;
	for (i = 0; i < h_; i++) {
		for (j = 0; j < other.w_; j++)
		{
			result.elements[other.w_ * i + j] = 0;
			for (k = 0; k < w_; k++)
				result.elements[other.w_ * i + j] = result.elements[other.w_ * i + j] + elements[w_ * i + k] * other.elements[other.w_ * k + j];
		}
	}
}

Matrix & Matrix::operator*=(const Matrix & other) {
	(*this) = (*this)*other;
	return *this;
}

Matrix & Matrix::operator/=(const float& f) {
	(*this) = (*this) / f;
	return *this;
}

Matrix transp(const Matrix & in){
	Matrix result(in.h_,in.w_);
	for(unsigned char i=0; i<in.w_;++i) for(unsigned char j=0; j<in.h_;++j) result.elements[j][i]=result.elements[i][j];
}

Matrix invert()const Matrix in){
	if(in.w_!=in.h_){
		in.NaNify;
	} else {
		// A = input matrix AND result matrix
		// n = number of rows = number of columns in A (n x n)
		unsigned char pivrow;		// keeps track of current pivot row
		unsigned char k, i, j;		// k: overall index along diagonal; i: row index; j: col index
		unsigned char pivrows[n]; // keeps track of rows swaps to undo at end
		float tmp;		// used for finding max value and making column swaps

		for (k = 0; k < n; k++)
		{
			// find pivot row, the row with biggest entry in current column
			tmp = 0;
			for (i = k; i < n; i++)
			{
				if (abs(in.elements[i * in.w_ + k]) >= tmp)	// 'Avoid using other functions inside abs()?'
				{
					tmp = abs(in.elements[i * in.w_ + k]);
					pivrow = i;
				}
			}

			// check for singular matrix
			if (in.elements[pivrow * in.w_ + k] == 0.0f)
			{
				in.NaNify;
				goto doneCalc;
			}

			// Execute pivot (row swap) if needed
			if (pivrow != k)
			{
				// swap row k with pivrow
				for (j = 0; j < n; j++)
				{
					tmp = in.elements[k * in.w_ + j];
					in.elements[k * in.w_ + j] = in.elements[pivrow * in.w_ + j];
					in.elements[pivrow * in.w_ + j] = tmp;
				}
			}
			pivrows[k] = pivrow;	// record row swap (even if no swap happened)

			tmp = 1.0f / in.elements[k * in.w_ + k];	// invert pivot element
			in.elements[k * in.w_ + k] = 1.0f;		// This element of input matrix becomes result matrix

			// Perform row reduction (divide every element by pivot)
			for (j = 0; j < n; j++)
			{
				in.elements[k * in.w_ + j] = in.elements[k * in.w_ + j] * tmp;
			}

			// Now eliminate all other entries in this column
			for (i = 0; i < n; i++)
			{
				if (i != k)
				{
					tmp = in.elements[i * in.w_ + k];
					in.elements[i * in.w_ + k] = 0.0f; // The other place where in matrix becomes result mat
					for (j = 0; j < n; j++)
					{
						in.elements[i * in.w_ + j] = in.elements[i * in.w_ + j] - in.elements[k * in.w_ + j] * tmp;
					}
				}
			}
		}

		// Done, now need to undo pivot row swaps by doing column swaps in reverse order
		for (k = in.w_ - 1; k >= 0; k--)
		{
			if (pivrows[k] != k)
			{
				for (i = 0; i < n; i++)
				{
					tmp = in.elements[i * in.w_ + k];
					in.elements[i * in.w_ + k] = in.elements[i * in.w_ + pivrows[k]];
					in.elements[i * in.w_ + pivrows[k]] = tmp;
				}
			}
		}
	}
	doneCalc:
	return in;
}

float Matrix::mag(){
	float resultSquared=0;
	for(unsigned char i=0; i<h_*w_;++i) resultSquared+=elements[i]*elements[i];
	return sqrt(resultSquared);
}