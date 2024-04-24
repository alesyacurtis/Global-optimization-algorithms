#include <iostream>
#include "Solve.h"
#include <math.h>
#include <vector>

using namespace std;

void Solver::Solve() {
	int K = 2;
	bool flag = 1;
	vector<Trial> trials(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(&minFunc, &argmin, &l, trials);

	while (flag && (K < Kmax)) {
		double m;
		int t;
		EstimateConstant(&K, trials, &m);

		FindMaxR_BA(&K, trials, &t, &m);
		MakeNewTrial_BA(&K, trials, &t, &m);

		//FindMaxR_SS(&K, trials, &t);
		//MakeNewTrial_SS(&K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		CheckStopCondition(&flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}

/*void Solver::SetParams() {//console input
	cout << "Input interval for function: ";
	cin >> a >> b;
	cout << "\nInput method parameter 2<=r<=4: ";
	cin >> r;
	cout << "\nInput accuracy: ";
	cin >> eps;
	cout << "\nInput max iterations quantity: ";
	cin >> Kmax;
}*/

Trial Solver:: GetSolution() {
	cout << "Current trial" << endl;
	cout << "z = " << NewTrial.z << endl;
	cout << "x = " << NewTrial.x << endl;
	cout << "k = " << NewTrial.k << endl;

	cout << "Best trial" << endl;
	cout << "z = " << BestTrial.z << endl;
	cout << "x = " << BestTrial.x << endl;
	cout << "k = " << BestTrial.k;
	return BestTrial;
}

void Solver::FirstTrial(double* minFunc, double* argmin, double* l, vector<Trial>& trials) {
	trials[0].x = a; trials[0].z = obj_func(trials[0].x); trials[0].k = 1;
	trials[1].x = b; trials[1].z = obj_func(trials[1].x); trials[1].k = 2;

	if (*minFunc > trials[1].z) {
		*minFunc = trials[1].z;
		*argmin = trials[1].x;
		*l = trials[1].k;
	}
}
void Solver::EstimateConstant(int *K, vector<Trial>& trials, double *m) {
	double M = 0, res;
	for (int i = 1; i < *K; i++) {
		res = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
		if (res > M)
			M = res;
	}
	if (M == 0) *m = 1;
	else *m = r * M;
}
void Solver::FindMaxR_BA(int* K, vector<Trial>& trials, int* t, double *m) { //Base Algorithm
	double max = -10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		R[i] = *m * (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (*m * (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z);
		if (R[i] > max) {
			max = R[i];
			*t = i;
		}
	}
}
void Solver::FindMaxR_SS(int* K, vector<Trial>& trials, int* t) { //Sequential Scan
	double max = -10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		R[i] = trials[i].x - trials[i - 1].x;
		if (R[i] > max) {
			max = R[i];
			*t = i;
		}
	}
}
void Solver::MakeNewTrial_BA(int* K, vector<Trial>& trials, int* t, double* m) {
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2 - ((trials[*t].z - trials[*t - 1].z) / (2 * *m));
	NewTrial.z = obj_func(NewTrial.x);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void Solver::MakeNewTrial_SS(int* K, vector<Trial>& trials, int* t) {
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2;
	NewTrial.z = obj_func(NewTrial.x);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void Solver::UpdateOptimum(double* minFunc, double* argmin, double* l, vector<Trial>& trials) {
	if (*minFunc > NewTrial.z) {
		*minFunc = NewTrial.z;
		*argmin = NewTrial.x;
		*l = NewTrial.k;
	}
}
void Solver::CheckStopCondition(bool *flag, vector<Trial>& trials, int* t) {
	if (trials[*t].x - trials[*t - 1].x <= eps * (b - a))
		*flag = 0;
	trials.insert(trials.begin() + *t, NewTrial);
}
void Solver::MakeBestTrial(double* minFunc, double* argmin, double* l) {
	BestTrial.z = *minFunc;
	BestTrial.x = *argmin;
	BestTrial.k = *l;
}