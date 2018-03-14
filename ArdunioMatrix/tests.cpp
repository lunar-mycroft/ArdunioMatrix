#include<cassert>
#include<iostream>

#include"matrix.hpp"

using std::cout;

int main() {
	{
		Matrix a;
		Matrix b;
		assert(a == b);
	}
}