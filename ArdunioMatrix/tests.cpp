#include<cassert>
#include<iostream>

using std::cout;

#include"matrix.hpp"

int main() {
	float f[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	float id3[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
	float z[9] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };

	float xp[3] = { 1.0,0.0,0.0 };
	float yp[3] = { 0.0,1.0,0.0 };
	float zp[3] = { 0.0,0.0,1.0 };
	float xn[3] = { -1.0,0.0,0.0 };
	float yn[3] = { 0.0,-1.0,0.0 };
	float zn[3] = { 0.0,0.0,-1.0 };

	{
		Matrix a;
		Matrix b;
		assert(a == b);
	}
	{
		Matrix a = identity(3);
		Matrix b((float *)id3, 3, 3);
		assert(a == b);
	}
	{
		Matrix a = identity(3);
		Matrix b((float *)f, 3, 3);
		assert(a != b);
	}
	{
		Matrix up((float *)zp, 3, 1);
		Matrix down((float *)zn, 3, 1);
		Matrix zero3((float *)z, 3, 1);
		Matrix zero1((float *)z, 1, 1);


		cout << int(up.height()) << ',' << int(up.width()) << '\n';

		assert(up + down == zero3);
		assert(-1 * up == down);
		Matrix prod = up*down;
		cout << prod[loc(0, 0)];
		cout << prod[loc(0, 1)];
		cout << prod[loc(0, 2)];
		assert(up*down == zero1);
	}
}