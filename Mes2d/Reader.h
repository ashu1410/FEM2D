#pragma once
#include <fstream>
using namespace std;

class Reader
{
public:
	double wysokoscR, dlugoscR;
	double condKR, cpR, gestoscR, tempOtoczeniaR, alfaR, tempPoczatkowaR;
	double krokCzasuR, czasSymulacjiR;
	int ileWezlWysR, ileWezlDlR;
	Reader();
	fstream data;
};
