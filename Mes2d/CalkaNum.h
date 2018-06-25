#pragma once

class CalkaNum {
public:
	float *tabPkt;
	float *tabWag;
	float wartoscFun;
	float suma;
	CalkaNum(int j);
	CalkaNum();
	float fun1d(float x);
	float fun2d(float x, float y);
	float fun3d(float x, float y, float z);
	void calkuj2d();
	void calkuj3d();
	void calkuj1d();
	float getSuma();
};