#pragma once
#pragma once
#include "Element.h"
#include "Node.h"
#include "Reader.h"
class Siatka {
public:
	vector<Node> vectWezl;
	vector<Element> vectElement;
	int ileWezlDl;
	int ileWezlWys;
	int ileEl;
	double wysokosc;
	double dlugosc;
	int ileWez;
	double krokWys;
	double krokDl;
	double condK, cp, gestosc;
	double tempOtoczenia, alfa, tempPoczatkowa;
	double krokCzasu, czasSymulacji;
	void zeroMatrix(double**, int, int);
	bool ludist(int, double**);
	bool lusolve(int, double **, double *, double *&);
	double max(double*, int);
	double min(double*, int);
	Siatka(Reader*);
	void tworzNode();
	void tworzElementy();

	//koniec siatki, Jakobian

	double matrixN[4][4];
	double **hGlobal;
	double **cGlobal;
	double *pGlobal;
	double N1(double, double);
	double N2(double, double);
	double N3(double, double);
	double N4(double, double);
	void funkcjeKsztaltu(double*, double*);
	double Jakobian();
	void Mes();
};