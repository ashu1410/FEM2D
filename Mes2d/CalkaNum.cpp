#include <iostream>
#include "CalkaNum.h"

using namespace std;

CalkaNum::CalkaNum(int j) {
	if (j == 3) {
		tabPkt = new float[3];
		tabWag = new float[3];
		tabPkt[0] = 0, 7745;
		tabPkt[1] = 0;
		tabPkt[2] = -0, 7745;
		tabWag[0] = (float)5 / (float)9;
		tabWag[1] = (float)8 / (float)9;
		tabWag[2] = (float)5 / (float)9;
		this->wartoscFun = 0;
		this->suma = 0;
	}
	else if (j == 2) {
		tabPkt = new float[2];
		tabWag = new float[2];
		tabPkt[0] = 0, 5773;
		tabPkt[1] = -0, 5773;
		tabWag[0] = 1;
		tabWag[1] = 1;
		this->wartoscFun = 0;
		this->suma = 0;
	}
}

CalkaNum::CalkaNum() {
	cout << "Z³y konstruktor" << endl;
}


float CalkaNum::fun2d(float x, float y) {

	this->wartoscFun = (5 * x*x * y* y) + 6 * x*y + 2 * x + 2 * y;
	return wartoscFun;
}

float CalkaNum::fun3d(float x, float y, float z) {
	this->wartoscFun = (2 * x*y*z) + x + 5 * z + 2 * y;
	return wartoscFun;
}

float CalkaNum::fun1d(float x) {
	this->wartoscFun = x + 2;
	return wartoscFun;
}

void CalkaNum::calkuj2d() {
	if (this->tabPkt[1] == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				this->suma += fun2d(this->tabPkt[i], this->tabPkt[j])*this->tabWag[i] * this->tabWag[j];
			}
		}
	}
	else {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				this->suma += fun2d(this->tabPkt[i], this->tabPkt[j])*this->tabWag[i] * this->tabWag[j];
			}
		}
	}
}

void CalkaNum::calkuj3d() {
	if (this->tabPkt[1] == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; i < 3; i++) {
					this->suma += fun3d(this->tabPkt[i], this->tabPkt[j], this->tabPkt[k])*this->tabWag[i] * this->tabWag[j] * this->tabWag[k];
				}
			}
		}
	}
	else {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; i < 3; i++) {
					this->suma += fun3d(this->tabPkt[i], this->tabPkt[j], this->tabPkt[k])*this->tabWag[i] * this->tabWag[j] * this->tabWag[k];
				}
			}
		}
	}
}

void CalkaNum::calkuj1d() {
	if (this->tabPkt[1] == 0) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				this->suma += fun1d(this->tabPkt[i])*this->tabWag[i];
			}
		}
	}
	else {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				this->suma += fun1d(this->tabPkt[i])*this->tabWag[i];
			}
		}
	}
}


float CalkaNum::getSuma() {
	return this->suma;
}