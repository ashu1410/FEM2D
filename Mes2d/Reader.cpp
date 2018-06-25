#include "Reader.h"
#include <iostream>

Reader::Reader()
{
	data.open("data.txt", ios::in);
	if (data.good())
	{
		data >> wysokoscR;
		data >> dlugoscR;
		data >> ileWezlWysR;
		data >> ileWezlDlR;
		data >> condKR;
		data >> cpR;
		data >> gestoscR;
		data >> tempOtoczeniaR;
		data >> alfaR;
		data >> czasSymulacjiR;
		data >> krokCzasuR;
		data >> tempPoczatkowaR;
		data.clear();
		data.close();
	}

	else
		cout << "Error reading file: data.txt" << endl;

}

