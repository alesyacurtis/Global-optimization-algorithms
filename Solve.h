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
	double eps = 0.001;
	double r = 3.6;

	double a = 0, b = 1;
	int Kmax = 1000;

	Trial BestTrial;
	Trial NewTrial;
	THillProblem hp;
	THillProblemFamily hpf;
	TShekelProblem shp;
	TShekelProblemFamily shf;

public:
	vector<Trial> trials;

	void Solve(int i);

	void FirstTrial(int i, double* minFunc, double* argmin, double* l, static vector<Trial>& trials);
	void EstimateConstant(int *K, static vector<Trial>& trials, double *m);
	void FindMaxR_BA(int* K, static vector<Trial>& trials, int* t, double* m); //Base Algorithm

	void FindMaxR_BA_UPD(int* K, double* minFunc, static vector<Trial>& trials, int* t, double* m); //Base Algorithm

	void FindMaxR_SS(int* K, static vector<Trial>& trials, int* t); //Sequential Scan
	void MakeNewTrial_BA(int i, int* K, static vector<Trial>& trials, int* t, double* m); //Base Algorithm
	void MakeNewTrial_SS(int i, int* K, static vector<Trial>& trials, int* t); //Sequential Scan
	void UpdateOptimum(double* minFunc, double* argmin, double* l, static vector<Trial>& trials);
	void CheckStopCondition(bool *flag, static vector<Trial>& trials, int* t);
	void CheckStopConditionWithKnownOptPoint(int i, bool* flag, static vector<Trial>& trials, int* t);
	void CheckStopConditionCurrentZ(int i, bool* flag, static vector<Trial>& trials, int* t);
	void MakeBestTrial(double* minFunc, double* argmin, double* l);

	int GetK();
	int GetKcur();
	int h, hh;
	Trial GetSolution();
	int mProblemIndex;
	
	double obj_func(double x, int i) {
		//return (x - 0.3) * (x - 0.3);
		//return (x - 1) * x * x - sin(20 * x);

		//return x * x - cos(18 * x) + sin(50 * x) + 1.878638;

		//return hp.ComputeFunction({ x });

		double known_min = hpf[i]->GetOptimumValue();
		return hpf[i]->ComputeFunction({ x }) - known_min;

		//return shp.ComputeFunction({ x });

		//double known_min = shf[i]->GetOptimumValue();
		//return shf[i]->ComputeFunction({ x }) - known_min; 
	}
};

class SolverNew : public Solver
{
public:
	void Solve(int i);
	void FindMinR(int* K, vector<Trial>& trials, int* t);
	void MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t);
};
class SolverRoot : public Solver
{
public:
	void Solve(int i);
	void FindMaxR(int* K, vector<Trial>& trials, int* t);
	void MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t);
};
class SolverModif : public Solver
{
protected:
	double r_glob = 3.6;
	double r_loc = 1.3;
	double ro = ((1 - 1 / r_glob) / (1 - 1 / r_loc)) * ((1 - 1 / r_glob) / (1 - 1 / r_loc));
public:
	void Solve(int i);
	void FindMaxR_BA_Modif(int* K, double* minFunc, vector<Trial>& trials, int* t, double* m, double* r_upd);
	void MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t, double* m, double* r_upd);
};

class File {
public:
	void InputToFile_percent(vector<int>&K, int n, int a);
	void InputToFile_x_solver(vector<Trial>& trials, int n, int a);
};

class Practice {
public:
	void BA(vector<int>& K, int n);
	void BA_Modif(vector<int>& K, int n);
	void NewA(vector<int>& K, int n);
	void RootA(vector<int>& K, int n);
};
