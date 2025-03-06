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
protected:
	double eps = 0.0001;
	double r = 2;

	double a = 0, b = 1;
	int Kmax = 5000;

	Trial BestTrial;
	Trial NewTrial;
	THillProblem hp;
	THillProblemFamily hpf;
	TShekelProblem shp;
	TShekelProblemFamily shf;

public:
	vector<Trial> trials;

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
	void CheckStopConditionWithKnownOptPoint(int i, bool* flag, vector<Trial>& trials, int* t);
	void CheckStopConditionCurrentZ(bool* flag, vector<Trial>& trials, int* t);
	void MakeBestTrial(double* minFunc, double* argmin, double* l);

	int GetK();
	int h;
	Trial GetSolution();
	int mProblemIndex;
	
	double obj_func(double x, int i) {
		//return (x - 0.3) * (x - 0.3);
		//return (x - 1) * x * x - sin(20 * x);

		//return hp.ComputeFunction({ x });
		double known_min = hpf[i]->GetOptimumValue();
		return hpf[i]->ComputeFunction({ x }) - known_min;
		//return shp.ComputeFunction({ x });
		//double known_min = shf[i]->GetOptimumValue();
		//return shf[i]->ComputeFunction({ x }) - known_min; //берет десятую функцию из семейства
	}
	/*double eps, r;//набор с консоли
	double a, b;
	int Kmax;
	
	void SetParams();
	*/
};

class SolverNew : public Solver
{
public:
	void Solve(int i);
	void FindMaxR(int* K, vector<Trial>& trials, int* t);
	void MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t);
};
class SolverRoot : public Solver
{
public:
	void Solve(int i);
	void FindMaxR(int* K, vector<Trial>& trials, int* t);
	void MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t);
};

class File {
public:
	void InputToFile_percent(vector<int>K, int n, int a);
	void InputToFile_x_solver(vector<Trial> trials, int n, int a);
};
