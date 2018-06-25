#include <iostream>
#include "Element.h"
#include "Node.h"
#include "Siatka.h"
#include "CalkaNum.h"
#include "Reader.h"
#include <cstdio>
using namespace std;




void main() {
	setlocale(LC_ALL, "");
	Reader *czytam = new Reader();
	Siatka siata(czytam);

	system("pause");
}

