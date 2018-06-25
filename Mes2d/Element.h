#pragma once
#include <vector>

using namespace std;

class Element {
public:
	vector <int> ID;
	vector <bool> pow;
	double localMatrixH[4][4];
	double localMatrixC[4][4];
	double localVectorP[4];
	double NpoXmatrix[4][4];
	double NpoYmatrix[4][4];

	Element();
	//Element(int, int, int, int);
};