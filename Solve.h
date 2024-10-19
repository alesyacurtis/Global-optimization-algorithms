#include <iostream>
#include <vector>
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblem.hpp"
#include "ShekelProblemFamily.hpp"

using namespace std;

struct Trial {
	double x, z;
	int k;
};

class Solver
{
	double eps = 0.001;
	double r=2;

	double a=0, b=1; 
	int Kmax = 200;

	Trial BestTrial;
	Trial NewTrial;
	THillProblem hp;
	THillProblemFamily hpf;
	TShekelProblem shp;
	TShekelProblemFamily shf;
public:
	void Solve(int i);
	//void BaseAlgorithm();
	//void SequentialScan();
	void FirstTrial(int i, double* minFunc, double* argmin, double* l, vector<Trial>& trials);
	void EstimateConstant(int *K, vector<Trial>& trials, double *m);
	void FindMaxR_BA(int* K, vector<Trial>& trials, int* t, double* m); //Base Algorithm
	void FindMaxR_SS(int* K, vector<Trial>& trials, int* t); //Sequential Scan
	void MakeNewTrial_BA(int i, int* K, vector<Trial>& trials, int* t, double* m); //Base Algorithm
	void MakeNewTrial_SS(int i, int* K, vector<Trial>& trials, int* t); //Sequential Scan
	void UpdateOptimum(double* minFunc, double* argmin, double* l, vector<Trial>& trials);
	void CheckStopCondition(bool *flag, vector<Trial>& trials, int* t);
	void MakeBestTrial(double* minFunc, double* argmin, double* l);

	Trial GetSolution();
	int mProblemIndex;
	
	double obj_func(double x, int i) {
		//return x * x - cos(18 * x);
		//return (x - 1) * (x - 5) * x + sin(18 * x);

		//return hp.ComputeFunction({ x });
		return hpf[i]->ComputeFunction({ x });
		//return shp.ComputeFunction({ x });
		//return shf[i]->ComputeFunction({ x }); //берет десятую функцию из семейства
	}

	/*double eps, r;//набор с консоли
	double a, b;
	int Kmax;
	
	void SetParams();
	*/
};
