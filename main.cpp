#include <iostream>
#include <fstream>
#include "Solve.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	Solver s;
	SolverModif sm;
	SolverNew sn;
	SolverRoot sr;
	File f;

	Practice p;

	int n = 1000;
	vector <int> K(n);

	p.BA(K,n);
	//f.InputToFile_percent(K, n, 1);

	p.BA_Modif(K, n);
	//f.InputToFile_percent(K, n, 2);

	//p.NewA(K, n);
	//f.InputToFile_percent(K, n, 3);

	//p.RootA(K,n);
	//f.InputToFile_percent(K, n, 4);
}
