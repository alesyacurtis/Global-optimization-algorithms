#include <iostream>
#include "Solve.h"
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

void Solver::Solve(int i) {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	trials.resize(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(i, &minFunc, &argmin, &l, trials);

	while (flag && (K < Kmax)) {
		double m;
		int t;
		EstimateConstant(&K, trials, &m);

		FindMaxR_BA(&K, trials, &t, &m);
		MakeNewTrial_BA(i, &K, trials, &t, &m);

		//FindMaxR_SS(&K, trials, &t);
		//MakeNewTrial_SS(i, &K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		//CheckStopCondition(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
		CheckStopConditionCurrentZ(&flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();

	/*cout << "Значения базового алгоритма:" << endl;
	BaseAlgorithm();
	cout << "\nЗначения алгоритма перебора:" << endl;
	SequentialScan();*/
}
void SolverNew::Solve(int i) {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	trials.resize(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(i, &minFunc, &argmin, &l, trials);

	while (flag && (K < Kmax)) {
		int t;

		FindMaxR(&K, trials, &t);
		MakeNewTrial(i, &K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		//CheckStopCondition(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
		CheckStopConditionCurrentZ(&flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}
void SolverRoot::Solve(int i) {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	trials.resize(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(i, &minFunc, &argmin, &l, trials);

	while (flag && (K < Kmax)) {
		int t;

		FindMaxR(&K, trials, &t);
		MakeNewTrial(i, &K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		CheckStopConditionCurrentZ(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}

/*void Solver::SetParams() {//набор с консоли
	cout << "Введите интервал, на котором будет рассматриваться функция: ";
	cin >> a >> b;
	cout << "\nВведите параметр метода 2<=r<=4: ";
	cin >> r;
	cout << "\nВведите точность: ";
	cin >> eps;
	cout << "\nВведите максимальное количество итераций: ";
	cin >> Kmax;
}*/

Trial Solver:: GetSolution() {
	cout << "\nТекущий результат" << endl;
	cout << "z = " << NewTrial.z << endl;
	cout << "x = " << NewTrial.x << endl;
	cout << "k = " << NewTrial.k << endl;

	cout << "Лучший результат" << endl;
	cout << "z = " << BestTrial.z << endl;
	cout << "x = " << BestTrial.x << endl;
	cout << "k = " << BestTrial.k;
	return BestTrial;
}
int Solver::GetK() {
	Solver::h;
	h = BestTrial.k;
	return h;
}

/*void Solver::BaseAlgorithm() {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	vector<Trial> trials(Kmax);
	trials[0].x = a; trials[0].z = obj_func(trials[0].x); trials[0].k = 1;
	trials[1].x = b; trials[1].z = obj_func(trials[1].x); trials[1].k = 2;

	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	if (minFunc > trials[1].z) {
		minFunc = trials[1].z;
		argmin = trials[1].x;
		l = trials[1].k;
	}

	while (flag && (K < Kmax)) {

		double M = 0, res, m;
		for (int i = 1; i < K; i++) {//оценим максимальное абсолютное значение относительной первой разности
			res = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
			if (res > M) //res = max(M,fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x)))
				M = res;
		}
		if (M == 0) m = 1;
		else m = r * M;

		double max = -10000;
		int t;
		vector <double> R(K);
		for (int i = 1; i < K; i++) {
			R[i] = m * (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (m * (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z);
			if (R[i] > max) {
				max = R[i];
				t = i;
			}
		}
		NewTrial.x = (trials[t].x + trials[t - 1].x) / 2 - ((trials[t].z - trials[t - 1].z) / (2 * m));
		NewTrial.z = obj_func(NewTrial.x);
		NewTrial.k = K + 1;
		K = K + 1;

		if (minFunc > NewTrial.z) {
			minFunc = NewTrial.z;
			argmin = NewTrial.x;
			l = NewTrial.k;
		}

		if (trials[t].x - trials[t - 1].x <= eps * (b - a))//условие остановки
			flag = 0;

		trials.insert(trials.begin() + t, NewTrial);
	}

	BestTrial.z = minFunc;
	BestTrial.x = argmin;
	BestTrial.k = l;
	GetSolution();
}

void Solver::SequentialScan() {
	int K = 2; //изначально 2 точки
	bool flag = 1;

	vector<Trial> trials(Kmax);
	trials[0].x = a; trials[0].z = obj_func(trials[0].x); trials[0].k = 1;
	trials[1].x = b; trials[1].z = obj_func(trials[1].x); trials[1].k = 2;

	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	if (minFunc > trials[1].z) {
		minFunc = trials[1].z;
		argmin = trials[1].x;
		l = trials[1].k;
	}

	while (flag && (K < Kmax)) {
		double max = -10000;
		int t;
		vector <double> R(K);
		for (int i = 1; i < K; i++) {
			R[i] = trials[i].x - trials[i - 1].x;
			if (R[i] > max) {
				max = R[i];
				t = i;
			}
		}
		NewTrial.x = (trials[t].x + trials[t - 1].x) / 2;
		NewTrial.z = obj_func(NewTrial.x);
		NewTrial.k = K + 1;
		K = K + 1;

		if (minFunc > NewTrial.z) {
			minFunc = NewTrial.z;
			argmin = NewTrial.x;
			l = NewTrial.k;
		}

		if (trials[t].x - trials[t - 1].x <= eps * (b - a))//условие остановки
			flag = 0;

		trials.insert(trials.begin() + t, NewTrial);
	}
	BestTrial.z = minFunc;
	BestTrial.x = argmin;
	BestTrial.k = l;
	GetSolution();
}*/

void Solver::FirstTrial(int i, double* minFunc, double* argmin, double* l, vector<Trial>& trials) {
	trials[0].x = a; trials[0].z = obj_func(trials[0].x,i); trials[0].k = 1;
	trials[1].x = b; trials[1].z = obj_func(trials[1].x,i); trials[1].k = 2;
	*minFunc = trials[0].z, * argmin = trials[0].x, * l = trials[0].k;

	if (*minFunc > trials[1].z) {
		*minFunc = trials[1].z;
		*argmin = trials[1].x;
		*l = trials[1].k;
	}
}
void Solver::EstimateConstant(int *K, vector<Trial>& trials, double *m) {
	double M = 0, res;
	for (int i = 1; i < *K; i++) {//оценим максимальное абсолютное значение относительной первой разности
		res = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
		if (res > M) //res = max(M,fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x)))
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
void SolverNew::FindMaxR(int* K, vector<Trial>& trials, int* t) { ////////// new
	double min = 10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		R[i] = (trials[i].z + trials[i - 1].z)/(trials[i].x - trials[i - 1].x);
		if (R[i] < min) {
			min = R[i];
			*t = i;
		}
	}
}
void SolverRoot::FindMaxR(int* K, vector<Trial>& trials, int* t) { ////////// root
	double min = 10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		if (trials[i].z * trials[i - 1].z >= 0) {
			R[i] = (trials[i].z * trials[i - 1].z) / (trials[i].x - trials[i - 1].x);
		}
		else R[i] = 0;
		if (R[i] < min) {
			min = R[i];
			*t = i;
		}
	}
}
void Solver::MakeNewTrial_BA(int i, int* K, vector<Trial>& trials, int* t, double* m) {//добавлено i для семейства
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2 - ((trials[*t].z - trials[*t - 1].z) / (2 * *m));
	NewTrial.z = obj_func(NewTrial.x,i);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void Solver::MakeNewTrial_SS(int i, int* K, vector<Trial>& trials, int* t) {
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2;
	NewTrial.z = obj_func(NewTrial.x,i);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void SolverNew::MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t) { ////////// new
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2 - ((trials[*t].z - trials[*t - 1].z) / (trials[*t].z + trials[*t - 1].z)) * (trials[*t].x - trials[*t - 1].x) * 1 / (2 * r);
	NewTrial.z = obj_func(NewTrial.x, i);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void SolverRoot::MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t) { ////////// root
	NewTrial.x = ((abs(trials[*t].z) * trials[*t - 1].x) + (abs(trials[*t - 1].z) * trials[*t].x)) / (abs(trials[*t].z) + abs(trials[*t - 1].z));
	NewTrial.z = obj_func(NewTrial.x, i);
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
	if (trials[*t].x - trials[*t - 1].x <= eps * (b - a))//условие остановки
		*flag = 0;
	trials.insert(trials.begin() + *t, NewTrial);
}
void Solver::CheckStopConditionWithKnownOptPoint(int i, bool* flag, vector<Trial>& trials, int* t) {
	if (abs(trials[*t].x - hpf[i]->GetOptimumPoint()[0]) <= eps * (b - a))//условие остановки
		*flag = 0;
	//if (abs(trials[*t].x - shf[i]->GetOptimumPoint()[0]) <= eps * (b - a))//условие остановки
	//	*flag = 0;
	trials.insert(trials.begin() + *t, NewTrial);
}
void Solver::CheckStopConditionCurrentZ(bool* flag, vector<Trial>& trials, int* t) {
	if (trials[*t].z <= eps * (b - a))//условие остановки
		*flag = 0;
	trials.insert(trials.begin() + *t, NewTrial);
}
void Solver::MakeBestTrial(double* minFunc, double* argmin, double* l) {
	BestTrial.z = *minFunc;
	BestTrial.x = *argmin;
	BestTrial.k = *l;
}

///////////////////////////////////////////////////////

void File::InputToFile_percent(vector<int>K, int n, int a) {
	fstream file;
	int p = 71;
	vector <double> P(p);
	for (int i = 1; i < n; i++) {
		int j = i - 1;
		while (j >= 0 && K[j] > K[j + 1])
		{
			swap(K[j], K[j + 1]);
			j--;
		}
	}
	int m = 0;
	int t = 0;
	bool f = 1;
	for (int i = 0; i < p; i++) {
		f = 1;
		m = 0;
		while (m < n && K[m] <= i * 10) {
			P[i] = P[i] + 1;
			m++;
			f = 0;
		}
		if (f) P[i] = 0;
	}
	if (a == 1)	file.open("percBA.txt", fstream::in | fstream::out);
	else if (a==2) file.open("percNewA.txt", fstream::in | fstream::out);
	else if (a==3) file.open("percRoot.txt", fstream::in | fstream::out);

	file << "K: ";
	for (int i = 0; i < p; i++) {
		if (i == 0) file << i * 10;
		else file << ", " << i * 10;
	}
	file << "\n" << "P(%): ";
	for (int i = 0; i < p; i++) {
		if (i == 0) file << P[i] * 100 / n;
		else file << ", " << P[i] * 100 / n;
	}
	file.close();
}
void File::InputToFile_x_solver (vector<Trial> trials, int n, int a) {
	fstream file;
	THillProblemFamily hpf;

	if (a == 1)	file.open("x_solverBA.txt", fstream::in | fstream::out);
	else if (a == 2) file.open("x_solverNewA.txt", fstream::in | fstream::out);
	else if (a == 3) file.open("x_solverRoot.txt", fstream::in | fstream::out);

	vector <int> b(201);
	int g = 0;
	for (int i = 0; i < 201; i++) {
		if (i == 0) file << "x: " << g * 0.01 / 2;
		else file << ", " << g * 0.01 / 2;
		g++;
	}
	g = 0;
	file << "\n\n";
	for (int i = 0; i < 201; i++) {
		if (i == 0) file << "phi: " << hpf[n]->ComputeFunction({g * 0.01 / 2}) - hpf[n]->GetOptimumValue();
		else file << ", " << hpf[n]->ComputeFunction({g * 0.01 / 2}) - hpf[n]->GetOptimumValue();
		g++;
	}
	file << "\n\n";
	for (int i = 0; i < trials.size(); i++) {
		if (i == 0) file << "x solver: " << trials[i].x;
		else file << ", " << trials[i].x;
	}
	file.close();
}
