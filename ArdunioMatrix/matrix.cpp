#include<math.h>

#include"matrix.hpp"

#define SQ(x) x*x
#define abs(x) (x<0) ? -x : x

bool canAdd(const Matrix & left, const Matrix & right) {
	return left.h_ == right.h_ && left.w_ == right.w_;
}

bool canMul(const Matrix & left, const Matrix & right) {
	return right.h_ == left.w_;
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

/*Matrix::Matrix(const imu::Quaternion & Q) : Matrix(0,0) {
	float q = SQ(Q.x())+SQ(Q.y())+SQ(Q.z())+SQ(Q.w());
	float s = (q == 0) ? 1 : (1.0 / q);
	elements[0] = 1 - 2 * s * (SQ(Q.y()) + SQ(Q.z())); elements[1]= 2 * s * (Q.x()*Q.y() - Q.z() * Q.w()); elements[2] = 2 * s*(Q.x()*Q.z() + Q.y() * Q.w());
	elements[3] = 2 * s * (Q.x()*Q.y() + Q.z() * Q.w()); elements[4] = 1 - 2 * s*(SQ(Q.x()) + SQ(Q.z())); elements[5] = 2 * s*(Q.y()*Q.z() - Q.x() * Q.w());
	elements[6] = 2 * s*(Q.x()*Q.z() + Q.y() * Q.w()); elements[7] = 2 * s*(Q.y()*Q.z() + Q.x() * Q.w()); elements[8] = 1 - 2 * s*(SQ(Q.x()) + SQ(Q.y()));

}


Matrix::Matrix(const imu::Vector<3>& v) {
	elements[0] = v[0];
	elements[1] = v[1];
	elements[2] = v[2];
}*/

Matrix::Matrix(const Matrix& other) {
	h_ = other.h_;
	w_ = other.w_;
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] = other.elements[i];
	}
}

Matrix & Matrix::operator=(Matrix other) {
	swap(*this, other);
	return *this;
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
	return result;
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
	return result;
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
	return result;
}

Matrix Matrix::operator/(float f) const{ 
	return (1.0 / f)*(*this); 
}

Matrix & Matrix::operator*=(const float& coeff) {
	for (unsigned char i = 0; i < h_*w_; ++i) {
		elements[i] *= coeff;
	}
	return *this;
}

Matrix Matrix::operator*(const Matrix & other) const {
	// toDo: max it "work" for two row vectors or two colum vectors, or two vectors in the wrong order

	if(!canMul(*this,other)){
		if (h_==1 && w_!=1 && other.h_==1 && other.w_!=1){
			return transp(*this)*other;
		} else if (w_==1 && h_!=1 && other.h_==1 && other.w_!=1){
			return (*this)*transp(other);
		} else if (h_==1 && w_!=1 && other.h_==1 && other.h_==1 && other.w_!=1){
			return transp(*this)*transp(other);
		} else {
			Matrix result(1,1);
			result.NaNify();
			return result;
		}
	}
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
	return result;
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
	for (unsigned char i = 0; i < in.w_; ++i) {
		for (unsigned char j = 0; j < in.h_; ++j) {
			result.elements[in.h_*j + i] = in.elements[in.w_*i+j];
		}
	}
	return result;
}


Matrix invert(const Matrix in){
	Matrix result(in);
	if (result.w_ != result.h_) {
		result.NaNify();
	}
	else {
		// from https://github.com/eecharlie/MatrixMath/blob/master/MatrixMath.cpp

		unsigned char pivrow;		// keeps track of current pivot row
		unsigned char k, i, j;		// k: overall index along diagonal; i: row index; j: col index
		unsigned char* pivrows = new unsigned char[result.w_]; // keeps track of rows swaps to undo at end
		float tmp;		// used for finding max value and making column swaps

		for (k = 0; k < result.w_; k++) {
			// find pivot row, the row with biggest entry in current column
			tmp = 0;
			for (i = k; i < result.w_; i++) {
				if (abs(result.elements[i * result.w_ + k]) >= tmp) {
					tmp = abs(result.elements[i * result.w_ + k]);
					pivrow = i;
				}
			}

			// check for singular matrix
			if (result.elements[pivrow * result.w_ + k] == 0.0f) {
				result.NaNify();
				goto doneCalc;
			}

			// Execute pivot (row swap) if needed
			if (pivrow != k) {
				// swap row k with pivrow
				for (j = 0; j < result.w_; j++) {
					tmp = result.elements[k * result.w_ + j];
					result.elements[k * result.w_ + j] = result.elements[pivrow * result.w_ + j];
					result.elements[pivrow * result.w_ + j] = tmp;
				}
			}
			pivrows[k] = pivrow;	// record row swap (even if no swap happened)

			tmp = 1.0f / result.elements[k * result.w_ + k];	// invert pivot element
			result.elements[k * result.w_ + k] = 1.0f;		// This element of input matrix becomes result matrix

															// Perform row reduction (divide every element by pivot)
			for (j = 0; j < result.w_; j++) {
				result.elements[k * result.w_ + j] = result.elements[k * result.w_ + j] * tmp;
			}

			// Now eliminate all other entries in this column
			for (i = 0; i < result.w_; i++) {
				if (i != k) {
					tmp = result.elements[i * result.w_ + k];
					result.elements[i * result.w_ + k] = 0.0f; // The other place where in matrix becomes result mat
					for (j = 0; j < result.w_; j++) {
						result.elements[i * result.w_ + j] = result.elements[i * result.w_ + j] - result.elements[k * result.w_ + j] * tmp;
					}
				}
			}
		}

		// Done, now need to undo pivot row swaps by doing column swaps in reverse order
		for (k = result.w_ - 1; k >= 0; k--) {
			if (pivrows[k] != k) {
				for (i = 0; i < result.w_; i++) {
					tmp = result.elements[i * result.w_ + k];
					result.elements[i * result.w_ + k] = result.elements[i * result.w_ + pivrows[k]];
					result.elements[i * result.w_ + pivrows[k]] = tmp;
				}
			}
		}
		delete[] pivrows;
	}
doneCalc:
	return result;
}


Matrix Matrix::cross(const Matrix & m) const {
	Matrix result(1,3);
	if (!((m.h_==1 && m.w_==3) || (m.w_==1 && m.h_==3)) && ((h_==1 && w_==3) || (w_==1 && h_==3))){
		result.NaNify();
		return result;
	}
	result.elements[0]=elements[1]*m.elements[2]-elements[2]*m.elements[1];
	result.elements[1]=elements[2]*m.elements[0]-elements[0]*m.elements[2];
	result.elements[3]=elements[0]*m.elements[1]-elements[1]*m.elements[0];
	return result;
}


float Matrix::mag() const {

	float resultSquared=0;
	for(unsigned char i=0; i<h_*w_;++i) resultSquared+=elements[i]*elements[i];
	return sqrt(resultSquared);
}

void swap(Matrix &a, Matrix &b){

	unsigned char tempW = a.w_;
	unsigned char tempH = a.h_;
	float * tempElements = a.elements;
	a.w_ = b.w_;
	a.h_ = b.h_;
	a.elements = b.elements;

	b.w_ = tempW;
	b.h_ = tempH;
	b.elements = tempElements;
}

void Matrix::NaNify() {
	for (unsigned char i = 0; i < w_*h_; ++i) {
		elements[i] = NAN;
	}
}

bool operator==(const Matrix &left, const Matrix & right) {
	if (!(left.h_ == right.h_ && left.w_ == right.w_)) return false;
	for (unsigned char i = 0; i < left.h_*left.w_; ++i) if (left.elements[i] != right.elements[i]) return false;
	return true;
}

Matrix rotBetweenVec(const Matrix & orig, const Matrix & target) {
	Matrix result(3, 3);
	if (!(orig.width() == 1 && orig.height() == 3 && target.width() == 1 && target.height() == 3)) result.NaNify();
	else {
		Matrix axis = orig.cross(target);

		float angle2vert = asin(mag(axis) / (mag(orig)*mag(target)));

		float Qw = cos(angle2vert / 2);
		float Qx = sin(angle2vert / 2)*axis.elements[0];
		float Qy = sin(angle2vert / 2)*axis.elements[1];
		float Qz = sin(angle2vert / 2)*axis.elements[2];

		float * R2 = result.elements;

		R2[0] = 1 - 2 * (SQ(Qy) + SQ(Qz)); R2[1] = 2 * (Qx*Qy - Qz * Qw); R2[2] = 2 * (Qx*Qz + Qy * Qw);
		R2[3] = 2 * (Qx*Qy + Qz * Qw); R2[4] = 1 - 2 * (SQ(Qx) + SQ(Qz)); R2[5] = 2 * (Qy*Qz - Qx * Qw);
		R2[6] = 2 * (Qx*Qz + Qy * Qw); R2[7] = 2 * (Qy*Qz + Qx * Qw); R2[8] = 1 - 2 * (SQ(Qx) + SQ(Qy));
	}
	return result;
}

Matrix identity(unsigned char size) {
	Matrix result(size, size);
	for (unsigned char i = 0; i < size; ++i) for (unsigned char j = 0; j < size; ++j) result.elements[i*size + j] = (i == j) ? 1.0 : 0.0;
	return result;
}