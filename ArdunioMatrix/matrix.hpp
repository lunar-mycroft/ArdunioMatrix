#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<iostream>;
using std::istream;
//#include"adafruitVector.hpp"
//#include"adafruitQuaternion.hpp"

class loc {
public:
	loc(unsigned char X, unsigned char Y) : x(X), y(Y) {};
	unsigned char x;
	unsigned char y;
};

class Matrix {
private:
	unsigned char w_;
	unsigned char h_;
	float* elements;

	void NaNify();
public:
	Matrix();
	Matrix(unsigned char, unsigned char);
	Matrix(float*, unsigned char, unsigned char);
	/*Matrix(const imu::Vector<3>&);
	Matrix(const imu::Quaternion &);  //Generates a rotation matrix from a quaternion*/


	//Big three
	Matrix(const Matrix&);
	Matrix & operator=(Matrix);
	~Matrix();

	//Index operator;
	float & operator[](loc);
	float   operator[](loc) const;

	//Dimentions
	unsigned char height() const {return h_;};
	unsigned char width() const { return w_; };

	//Math operators

	Matrix operator+(const Matrix&) const;
	Matrix & operator+=(const Matrix&);

	Matrix operator-(const Matrix&) const;
	Matrix & operator-=(const Matrix&);

	Matrix operator*(const float &) const;
	Matrix & operator*=(const float&);

	Matrix operator*(const Matrix &) const;
	Matrix & operator*=(const Matrix&);

	Matrix cross(const Matrix &) const; //Cross product

	Matrix operator/(float) const;
	Matrix & operator/=(const float&);

	Matrix operator/(const Matrix & other) const { return (*this)*invert(other); }

	//Math functions

	friend Matrix invert(const Matrix);
	friend Matrix transp(const Matrix &);

	float mag() const;

	friend bool canAdd(const Matrix &, const Matrix &);
	friend bool canMul(const Matrix &, const Matrix &);
	friend void swap(Matrix &, Matrix &);

	friend float toFloat(const Matrix & m) { return m.elements[0]; };

	friend bool operator==(const Matrix &, const Matrix &);
	friend Matrix rotBetweenVec(const Matrix & orig, const Matrix & target); //Returns a rotation matrix which will turn 

	friend Matrix identity(unsigned char);

};

Matrix operator*(float f, const Matrix & m);
Matrix operator/(float f, const Matrix & m);

bool operator!=(const Matrix & left, const Matrix & right);

float mag(const Matrix & m);

#endif
