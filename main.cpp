#include <iostream>
#include "Solve.h"

using namespace std;
int main()
{
	setlocale(LC_ALL, "Russian");
	Solver s;
	THillProblem hp;
	THillProblemFamily hpf;

	for (int i = 0; i < 10; i++) {
		s.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
	}

	/*s.Solve(0);
	cout << "\nOptimumValue = " << hp.GetOptimumValue() << endl;
	cout << "OptimumPoint = " << hp.GetOptimumPoint()[0] << endl;*/

	//Trial t = s.GetSolution();
}