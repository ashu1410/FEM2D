#include <iostream>
#include "Siatka.h"
#include <math.h>
#include<stdio.h>
#include "Element.h"
using namespace std;

/*Siatka::Siatka(double dl,double wys,int wezWys,int wezDl) {
this->dlugosc = dl;
this->wysokosc = wys;
this->ileWezlDl = wezDl;
this->ileWezlWys = wezWys;
this->krokDl = (double)this->dlugosc / (double)this->ileWezlDl;
this->krokWys = (double)this->wysokosc / (double) this->ileWezlWys;
this->ileWez = this->ileWezlDl*this->ileWezlWys;
this->ileEl = (this->ileWezlDl - 1)*(this->ileWezlWys - 1);
this->vectWezl = new Node[ileWez];
this->vectElement = new Element[ileEl];
funkcjeKsztaltu();
//cout << "ile wezlow" << this->ileWez << "ile elementow" << this->ileEl << "ile wezlow dl"<< this->ileWezlDl<<"ile wysokosc"<<this->ileWezlWys<<endl;
}*/


Siatka::Siatka(Reader *reader) {
	this->wysokosc = reader->wysokoscR;
	this->dlugosc = reader->dlugoscR;
	this->ileWezlWys = reader->ileWezlWysR;
	this->ileWezlDl = reader->ileWezlDlR;
	this->condK = reader->condKR;
	this->cp = reader->cpR;
	this->gestosc = reader->gestoscR;
	this->tempOtoczenia = reader->tempOtoczeniaR;
	this->alfa = reader->alfaR;
	this->czasSymulacji = reader->czasSymulacjiR;
	this->krokCzasu = reader->krokCzasuR;
	this->tempPoczatkowa = reader->tempPoczatkowaR;
	this->ileWez = (this->ileWezlDl*this->ileWezlWys);
	this->ileEl = (this->ileWezlDl - 1)*(this->ileWezlWys - 1);
	this->vectWezl.resize(ileWez);
	this->vectElement.resize(ileEl);
	this->tworzNode();
	this->tworzElementy();
	this->hGlobal = new double*[ileWez];
	this->cGlobal = new double*[ileWez];
	this->pGlobal = new double[ileWez];
	for (int i = 0; i < ileWez; i++) {
		this->hGlobal[i] = new double[ileWez];
		this->cGlobal[i] = new double[ileWez];
		for (int j = 0; j < ileWez; j++) {
			this->hGlobal[i][j] = 0.0;
			this->cGlobal[i][j] = 0.0;

		}
	}
	for (int i = 0; i < ileWez; i++) {
		this->pGlobal[i] = 0.0;
	}
	this->Jakobian();
	this->Mes();
}




void Siatka::tworzNode()
{
	int column = 0;
	for (int i = 0; i < ileWez; i++) 
	{
		if (i % ileWezlWys == 0 && i != 0)
			column++;

		if (i % ileWezlWys == 0) {
			vectWezl[i].y = 0;
			vectWezl[i].status = true;

		}
		else {
			vectWezl[i].y = (wysokosc / ((double)ileWezlWys - 1)) * (i % ileWezlWys);
		}
		if (i == 0) {
			vectWezl[i].x = 0;
		}
		else {
			vectWezl[i].x = (dlugosc / ((double)ileWezlDl - 1)) * column;
		}

		if (column == 0 || column == (ileWezlWys - 1) || (i %ileWezlWys == ileWezlWys - 1))
			vectWezl[i].status = true;



	}
}

void Siatka::tworzElementy()
{
	for (int i = 0; i < ileEl; i++)
	{

		for (int k = 0; k < 4; k++) {
			vectElement[i].localVectorP[k] = 0.0;
			for (int j = 0; j < 4; j++) {
				vectElement[i].localMatrixH[k][j] = 0.0;
				vectElement[i].localMatrixC[k][j] = 0.0;

				vectElement[i].NpoXmatrix[k][j] = 0.0;
				vectElement[i].NpoYmatrix[k][j] = 0.0;
			}
		}

		if (i == 0)
			vectElement[i].ID[0] = 0;

		else
			if (i % (ileWezlWys- 1) != 0)
				vectElement[i].ID[0] =vectElement[i - 1].ID[3];
			else
				vectElement[i].ID[0] = vectElement[i - 1].ID[3] + 1;

		vectElement[i].ID[1] = vectElement[i].ID[0] + ileWezlWys;
		vectElement[i].ID[2] = vectElement[i].ID[1] + 1;
		vectElement[i].ID[3] = vectElement[i].ID[0] + 1;

	}


	for (int i = 0; i < ileEl; i++)
	{
		vectElement[i].pow[0] = vectWezl[vectElement[i].ID[0]].status && vectWezl[vectElement[i].ID[3]].status ? 1 : 0;
		vectElement[i].pow[1] = vectWezl[vectElement[i].ID[0]].status && vectWezl[vectElement[i].ID[1]].status ? 1 : 0;
		vectElement[i].pow[2] = vectWezl[vectElement[i].ID[1]].status && vectWezl[vectElement[i].ID[2]].status ? 1 : 0;
		vectElement[i].pow[3] = vectWezl[vectElement[i].ID[2]].status && vectWezl[vectElement[i].ID[3]].status ? 1 : 0;
	}


}

double Siatka::N1(double ksi, double eta) {
	return ((1 - ksi)*(1 - eta)*0.25);
}
double Siatka::N2(double ksi, double eta) {
	return ((1 + ksi)*(1 - eta)*0.25);
}
double Siatka::N3(double ksi, double eta) {
	return ((1 + ksi)*(1 + eta)*0.25);
}
double Siatka::N4(double ksi, double eta) {
	return ((1 - ksi)*(1 + eta)*0.25);
}

void Siatka::funkcjeKsztaltu(double *ksi, double *eta) {
	this->matrixN[0][0] = N1(ksi[0], eta[0]);
	this->matrixN[0][1] = N2(ksi[0], eta[0]);
	this->matrixN[0][2] = N3(ksi[0], eta[0]);
	this->matrixN[0][3] = N4(ksi[0], eta[0]);

	this->matrixN[1][0] = N1(ksi[1], eta[1]);
	this->matrixN[1][1] = N2(ksi[1], eta[1]);
	this->matrixN[1][2] = N3(ksi[1], eta[1]);
	this->matrixN[1][3] = N4(ksi[1], eta[1]);

	this->matrixN[2][0] = N1(ksi[2], eta[2]);
	this->matrixN[2][1] = N2(ksi[2], eta[2]);
	this->matrixN[2][2] = N3(ksi[2], eta[2]);
	this->matrixN[2][3] = N4(ksi[2], eta[2]);

	this->matrixN[3][0] = N1(ksi[3], eta[3]);
	this->matrixN[3][1] = N2(ksi[3], eta[3]);
	this->matrixN[3][2] = N3(ksi[3], eta[3]);
	this->matrixN[3][3] = N4(ksi[3], eta[3]);

}


double Siatka::Jakobian() {
	double ksi[4] = { -0.5773,0.5773,0.5773,-0.5773 };
	double eta[4] = { -0.5773,-0.5773,0.5773,0.5773 };

	funkcjeKsztaltu(ksi, eta);
	/*for (int a = 0; a < 4; a++) {
	for (int b = 0; b < 4; b++)
	cout << "Funkcje ksztaltu" << matrixN[a][b] << endl;
	}*/
	double NpoKsi[4][4];
	double NpoEta[4][4];
	double JakxDetJ[2][2];
	double XpoKsi = 0.0, XpoEta = 0.0, YpoKsi = 0.0, YpoEta = 0.0, detJ = 0.0, detOrig = 0.0, jakob = 0.0;
	bool warunek = false;
	double npc[2][2];
	double tmpLocalVectorP[4][4];
	double tmpLocal[4][4];
	double calka[4][4];
	double tmpNP[2][4];

	NpoKsi[0][0] = (-0.25)*(1 - eta[0]); //pochodne f.k po lokalnych
	NpoKsi[0][1] = (0.25)*(1 - eta[0]);
	NpoKsi[0][2] = (0.25)*(1 + eta[0]);
	NpoKsi[0][3] = (-0.25)*(1 + eta[0]);

	NpoKsi[1][0] = (-0.25)*(1 - eta[1]);
	NpoKsi[1][1] = (0.25)*(1 - eta[1]);
	NpoKsi[1][2] = (0.25)*(1 + eta[1]);
	NpoKsi[1][3] = (-0.25)*(1 + eta[1]);

	NpoKsi[2][0] = (-0.25)*(1 - eta[2]);
	NpoKsi[2][1] = (0.25)*(1 - eta[2]);
	NpoKsi[2][2] = (0.25)*(1 + eta[2]);
	NpoKsi[2][3] = (-0.25)*(1 + eta[2]);

	NpoKsi[3][0] = (-0.25)*(1 - eta[3]);
	NpoKsi[3][1] = (0.25)*(1 - eta[3]);
	NpoKsi[3][2] = (0.25)*(1 + eta[3]);
	NpoKsi[3][3] = (-0.25)*(1 + eta[3]);

	NpoEta[0][0] = (-0.25)*(1 - ksi[0]);
	NpoEta[0][1] = (-0.25)*(1 + ksi[0]);
	NpoEta[0][2] = (0.25)*(1 + ksi[0]);
	NpoEta[0][3] = (0.25)*(1 - ksi[0]);

	NpoEta[1][0] = (-0.25)*(1 - ksi[1]);
	NpoEta[1][1] = (-0.25)*(1 + ksi[1]);
	NpoEta[1][2] = (0.25)*(1 + ksi[1]);
	NpoEta[1][3] = (0.25)*(1 - ksi[1]);

	NpoEta[2][0] = (-0.25)*(1 - ksi[2]);
	NpoEta[2][1] = (-0.25)*(1 + ksi[2]);
	NpoEta[2][2] = (0.25)*(1 + ksi[2]);
	NpoEta[2][3] = (0.25)*(1 - ksi[2]);

	NpoEta[3][0] = (-0.25)*(1 - ksi[3]);
	NpoEta[3][1] = (-0.25)*(1 + ksi[3]);
	NpoEta[3][2] = (0.25)*(1 + ksi[3]);
	NpoEta[3][3] = (0.25)*(1 - ksi[3]);


	for (vector<Element>::iterator element = this->vectElement.begin(); element != this->vectElement.end(); element++) {
		XpoKsi = (0.5773 / 4.0)*(this->vectWezl[element->ID[0]].x - (this->vectWezl[element->ID[1]].x) + (this->vectWezl[element->ID[2]].x) - this->vectWezl[element->ID[3]].x) +
			0.25*(-this->vectWezl[element->ID[0]].x + (this->vectWezl[element->ID[1]].x) + (this->vectWezl[element->ID[2]].x) - this->vectWezl[element->ID[3]].x);

		XpoEta = (0.5773 / 4.0)*(this->vectWezl[element->ID[0]].x - (this->vectWezl[element->ID[1]].x) + (this->vectWezl[element->ID[2]].x) - this->vectWezl[element->ID[3]].x) +
			0.25*(-this->vectWezl[element->ID[0]].x - (this->vectWezl[element->ID[1]].x) + (this->vectWezl[element->ID[2]].x) + this->vectWezl[element->ID[3]].x);

		YpoKsi = (0.5773 / 4.0)*(this->vectWezl[element->ID[0]].y - this->vectWezl[element->ID[1]].y + (this->vectWezl[element->ID[2]].y) - (this->vectWezl[element->ID[3]].y)) +
			0.25*(-this->vectWezl[element->ID[0]].y + this->vectWezl[element->ID[1]].y + (this->vectWezl[element->ID[2]].y) - (this->vectWezl[element->ID[3]].y));

		YpoEta = (0.5773 / 4.0)*(this->vectWezl[element->ID[0]].y - this->vectWezl[element->ID[1]].y + (this->vectWezl[element->ID[2]].y) - (this->vectWezl[element->ID[3]].y)) +
			0.25*(-this->vectWezl[element->ID[0]].y - this->vectWezl[element->ID[1]].y + (this->vectWezl[element->ID[2]].y) + (this->vectWezl[element->ID[3]].y));

		detOrig = ((XpoKsi*YpoEta) - (XpoEta*YpoKsi));
		detJ = 1 / detOrig;


		JakxDetJ[0][0] = detJ*YpoEta;
		JakxDetJ[0][1] = -detJ*YpoKsi;
		JakxDetJ[1][0] = -detJ*XpoEta;
		JakxDetJ[1][1] = detJ*XpoKsi;


		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				element->NpoXmatrix[j][k] += JakxDetJ[0][0] * NpoKsi[j][k] + JakxDetJ[0][1] * NpoEta[j][k];
				element->NpoYmatrix[j][k] += JakxDetJ[1][0] * NpoKsi[j][k] + JakxDetJ[1][1] * NpoEta[j][k];
			}
		}


		for (int k = 0; k < 4; k++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {//pierwsza calka po objetosci
					element->localMatrixH[i][j] += detOrig*this->condK*((element->NpoXmatrix[k][j] * element->NpoXmatrix[k][i]) + (element->NpoYmatrix[k][j] * element->NpoYmatrix[k][i]));
				}
			}


			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {//pierwsza czesc do macierzy C z f kszt
					element->localMatrixC[i][j] += detOrig*(matrixN[k][j] * matrixN[k][i]);
				}
				//cout << endl;
			}

		}
		//cout << "Macierz cLokalna:\n";
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				//this->element->localMatrixH[i][j] *= detOrig*this->_condK;
				element->localMatrixC[i][j] *= (this->cp*this->gestosc);
				//printf("%lf\t", this->element->localMatrixC[i][j]);
			}
			//	cout << endl;
		}



		jakob = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				tmpLocal[i][j] = tmpLocalVectorP[i][j] = calka[i][j] = 0.0;
			}
		}

		for (int i = 0; i < 4; i++) {
			if (element->pow[i] == 1) {
				warunek = true;
				switch (i) {
				case 0:
					jakob = (this->vectWezl[element->ID[3]].y - this->vectWezl[element->ID[0]].y) / 2;
					npc[0][0] = -1.0;
					npc[0][1] = -0.5773502691896258;
					npc[1][0] = -1.0;
					npc[1][1] = 0.5773502691896258;
					break;
				case 1:
					jakob = (this->vectWezl[element->ID[1]].x - this->vectWezl[element->ID[0]].x) / 2;
					npc[0][0] = -0.5773502691896258;
					npc[0][1] = -1.0;
					npc[1][0] = 0.5773502691896258;
					npc[1][1] = -1.0;
					break;
				case 2:
					jakob = (this->vectWezl[element->ID[2]].y - this->vectWezl[element->ID[1]].y) / 2;
					npc[0][0] = 1.0;
					npc[0][1] = -0.5773502691896258;
					npc[1][0] = 1.0;
					npc[1][1] = 0.5773502691896258;
					break;
				case 3:
					jakob = (this->vectWezl[element->ID[2]].x - this->vectWezl[element->ID[3]].x) / 2;
					npc[0][0] = -0.5773502691896258;
					npc[0][1] = 1.0;
					npc[1][0] = 0.5773502691896258;
					npc[1][1] = 1.0;
					break;
				}
				for (int n = 0; n < 2; n++) {
					tmpNP[n][0] = this->N1(npc[n][0], npc[n][1]);
					tmpNP[n][1] = this->N2(npc[n][0], npc[n][1]);
					tmpNP[n][2] = this->N3(npc[n][0], npc[n][1]);
					tmpNP[n][3] = this->N4(npc[n][0], npc[n][1]);
				}

				for (int j = 0; j < 4; j++) {//same funkje ksztaltu na razie
					tmpLocalVectorP[i][j] = jakob*(tmpNP[0][j] + tmpNP[1][j]);
				}

				for (int k = 0; k < 4; k++) {
					for (int j = 0; j < 4; j++) {// {N}{N}T
						calka[j][k] = calka[j][k] + ((tmpNP[0][k] * tmpNP[0][j])*jakob);
						calka[j][k] = calka[j][k] + ((tmpNP[1][k] * tmpNP[1][j])*jakob);
					}
				}
			}
		}



		//	lokal P
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				element->localVectorP[i] += tmpLocalVectorP[j][i];
			}
		}

		if (warunek) {//dorzucamy 2 calke do H gdy na brzegu
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					element->localMatrixH[i][j] += this->alfa*(calka[i][j]);
				}
			}
		}


		//   A G R E G A C J E
		for (int i = 0; i < 4; i++) {
			element->localVectorP[i] = element->localVectorP[i] * this->alfa * this->tempOtoczenia;
		}
	}

	for (int k = 0; k < this->vectElement.size(); k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->hGlobal[vectElement[k].ID[i]][vectElement[k].ID[j]] += vectElement[k].localMatrixH[i][j];
				this->cGlobal[vectElement[k].ID[i]][vectElement[k].ID[j]] += vectElement[k].localMatrixC[i][j];
			}
		}
	}

	for (int k = 0; k < this->vectElement.size(); k++) {
		for (int j = 0; j < 4; j++) {
			this->pGlobal[this->vectElement[k].ID[j]] += this->vectElement[k].localVectorP[j];
		}
	}

	return 0;
}

void Siatka::Mes() {
	double **h = new double*[this->ileWez];
	double **cDivided = new double*[this->ileWez];
	double *p = new double[this->ileWez];
	double *t0 = new double[this->ileWez];
	double *x = new double[this->ileWez];
	for (int i = 0; i < this->ileWez; i++) {
		h[i] = new double[this->ileWez];
		t0[i] = this->tempPoczatkowa;
		cDivided[i] = new double[this->ileWez];
		x[i] = 0;
		for (int j = 0; j < this->ileWez; j++) {
			h[i][j] = this->hGlobal[i][j] + (this->cGlobal[i][j] / this->krokCzasu);//pierwsza czesc rownania niestacjonarnego H+C/dTau
			cDivided[i][j] = this->cGlobal[i][j] / this->krokCzasu;//do drugiej czesci C/dTau
		}
	}

	double minTemp = 0;
	double maxTemp = 0;
	cout << "Time[s]\tMinTemp[s]\tMaxTemp[s]\n";
	for (int i = this->krokCzasu; i <= this->czasSymulacji; i += this->krokCzasu) {
		for (int l = 0; l < this->ileWez; l++) {
			p[l] = x[l] = 0.0;
		}
		for (int j = 0; j < this->ileWez; j++) {

			for (int k = 0; k < this->ileWez; k++) {
				p[j] += cDivided[j][k] * t0[k];
			}
			p[j] += this->pGlobal[j];
		}

		if (lusolve(this->ileWez, h, p, x)) {
			minTemp = min(x, this->ileWez);
			maxTemp = max(x, this->ileWez);
			for (int l = 0; l < this->ileWez; l++) {
				t0[l] = x[l];
			}
		}

		cout << i << "\t" << minTemp << "\t\t" << maxTemp << endl;
	}

}

void Siatka::zeroMatrix(double **matrix, int w, int h) {
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			matrix[i][j] = 0;
		}
	}
}

const double eps = 1e-12;
bool Siatka::ludist(int n, double ** A)
{
	int i, j, k;
	for (k = 0; k < n - 1; k++) {
		if (fabs(A[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}
	return true;
}

bool Siatka::lusolve(int n, double ** A, double * B, double *& X)
{
	int    i, j;
	double s;
	double **tab = new double*[n];
	for (int z = 0; z < n; z++) {
		tab[z] = new double[n];
		for (int x = 0; x < n; x++) {
			tab[z][x] = A[z][x];
		}
	}
	ludist(n, tab);
	X[0] = B[0];

	for (i = 1; i < n; i++) {
		s = 0;

		for (j = 0; j < i; j++) s += tab[i][j] * X[j];

		X[i] = B[i] - s;
	}

	if (fabs(tab[n - 1][n - 1]) < eps) return false;

	X[n - 1] /= tab[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--) {
		s = 0;

		for (j = i + 1; j < n; j++) s += tab[i][j] * X[j];

		if (fabs(tab[i][i]) < eps) return false;

		X[i] = (X[i] - s) / tab[i][i];
	}
	for (int i = 0; i < n; i++) {
		delete tab[i];
	}
	delete tab;
	return true;
}

double Siatka::max(double *vec, int n) {
	double max = vec[0];
	for (int i = 1; i < n - 1; i++) {
		if (vec[i] > max) {
			max = vec[i];
		}
	}
	return max;
}
double Siatka::min(double *vec, int n) {
	double min = vec[0];
	for (int i = 1; i < n - 1; i++) {
		if (vec[i] < min) {
			min = vec[i];
		}
	}
	return min;
}
