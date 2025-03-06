#include <iostream>
#include <fstream>
#include "Solve.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	Solver s;
	SolverNew sn;
	SolverRoot sr;
	File f;

	THillProblem hp;
	THillProblemFamily hpf;
	TShekelProblem shp;
	TShekelProblemFamily shf;

	int n = 1;
	vector <int> K(n);

	for (int i = 0; i < n; i++) {
		s.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = s.GetK();
	}
	f.InputToFile_percent(K, n, 1);
	//f.InputToFile_x_solver(s.trials, n-1, 1);

	for (int i = 0; i < n; i++) {
		sn.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = sn.GetK();
	}
	f.InputToFile_percent(K, n, 2);
	//f.InputToFile_x_solver(sn.trials, n-1, 2);

	for (int i = 0; i < n; i++) {
		sr.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = sr.GetK();
	}
	f.InputToFile_percent(K, n, 3);
	//f.InputToFile_x_solver(sr.trials, n - 1, 3);

}
