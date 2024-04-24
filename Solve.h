#ifndef GLOBAL_OPTIMIZATION_ALGORITHMS_SOLVE_H
#define GLOBAL_OPTIMIZATION_ALGORITHMS_SOLVE_H

#include <iostream>
#include <vector>

using namespace std;

struct Trial {
	double x, z;
	int k;
};

class Solver
{
	double eps = 0.001;
	double r=2;

	double a=-0.5, b=1.5; 
	int Kmax = 200;

	Trial BestTrial;
	Trial NewTrial;
public:
	void Solve();
	void FirstTrial(double* minFunc, double* argmin, double* l, vector<Trial>& trials);
	void EstimateConstant(int *K, vector<Trial>& trials, double *m);
	void FindMaxR_BA(int* K, vector<Trial>& trials, int* t, double* m); //Base Algorithm
	void FindMaxR_SS(int* K, vector<Trial>& trials, int* t); //Sequential Scan
	void MakeNewTrial_BA(int* K, vector<Trial>& trials, int* t, double* m); //Base Algorithm
	void MakeNewTrial_SS(int* K, vector<Trial>& trials, int* t); //Sequential Scan
	void UpdateOptimum(double* minFunc, double* argmin, double* l, vector<Trial>& trials);
	void CheckStopCondition(bool *flag, vector<Trial>& trials, int* t);
	void MakeBestTrial(double* minFunc, double* argmin, double* l);

	Trial GetSolution();

	double obj_func(double x) {
		return x * x - cos(18 * x);
		//return (x - 1) * (x - 5) * x + sin(18 * x);
	}

	/*double eps, r;//               
	double a, b;
	int Kmax;
	
	void SetParams();
	*/
};
#endif//GLOBAL_OPTIMIZATION_ALGORITHMS_SOLVE_H