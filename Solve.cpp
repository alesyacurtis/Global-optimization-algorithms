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

		//FindMaxR_BA(&K, trials, &t, &m);
		FindMaxR_BA_UPD(&K, &minFunc, trials, &t, &m);
		MakeNewTrial_BA(i, &K, trials, &t, &m);

		//FindMaxR_SS(&K, trials, &t);
		//MakeNewTrial_SS(i, &K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		//CheckStopCondition(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
		CheckStopConditionCurrentZ(i, &flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}
void SolverModif::Solve(int i) {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	trials.resize(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(i, &minFunc, &argmin, &l, trials);
	while (flag && (K < Kmax)) {
		double m;
		double r_upd;
		int t;
		EstimateConstant(&K, trials, &m);

		FindMaxR_BA_Modif(&K, &minFunc, trials, &t, &m,&r_upd);
		MakeNewTrial(i, &K, trials, &t, &m, &r_upd);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		//CheckStopCondition(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
		CheckStopConditionCurrentZ(i, &flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}
void SolverNew::Solve(int i) {
	int K = 2; //изначально 2 точки
	bool flag = 1;
	trials.resize(Kmax);
	double minFunc = trials[0].z, argmin = trials[0].x, l = trials[0].k;
	FirstTrial(i, &minFunc, &argmin, &l, trials);

	while (flag && (K < Kmax)) {
		int t;

		FindMinR(&K, trials, &t);
		MakeNewTrial(i, &K, trials, &t);

		UpdateOptimum(&minFunc, &argmin, &l, trials);
		//CheckStopCondition(&flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
		CheckStopConditionCurrentZ(i, &flag, trials, &t);
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
		CheckStopConditionCurrentZ(i, &flag, trials, &t);
		//CheckStopConditionWithKnownOptPoint(i, &flag, trials, &t);
	}
	MakeBestTrial(&minFunc, &argmin, &l);
	GetSolution();
}

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
int Solver::GetKcur() {
	Solver::hh;
	hh = NewTrial.k;
	return hh;
}

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
	else *m = M;
}
void Solver::FindMaxR_BA(int* K, vector<Trial>& trials, int* t, double *m) { //Base Algorithm
	double max = -10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		R[i] = *m * r * (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (*m * r * (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z);
		if (R[i] > max) {
			max = R[i];
			*t = i;
		}
	}
}
void Solver::FindMaxR_BA_UPD(int* K, double* minFunc, vector<Trial>& trials, int* t, double* m) { //Base Algorithm
	double max = -10000;
	vector <double> R(*K);
	for (int i = 1; i < *K; i++) {
		R[i] = (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (*m * *m *r*r* (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z - 2 * *minFunc)/ (*m*r);
		if (R[i] > max) {
			max = R[i];
			*t = i;
		}
	}
}
void SolverModif::FindMaxR_BA_Modif(int* K, double* minFunc, vector<Trial>& trials, int* t, double* m, double* r_upd) {
	double max = -10000;
	vector <double> R(*K);
	vector <double> R_glob(*K);
	vector <double> R_loc(*K);
	for (int i = 1; i < *K; i++) {
		R_glob[i] = (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (*m * *m * r_glob * r_glob * (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z - 2 * *minFunc) / (*m * r_glob);
		R_loc[i] = (trials[i].x - trials[i - 1].x) + (((trials[i].z - trials[i - 1].z) * (trials[i].z - trials[i - 1].z)) / (*m * *m * r_loc * r_loc * (trials[i].x - trials[i - 1].x))) - 2 * (trials[i].z + trials[i - 1].z - 2 * *minFunc) / (*m * r_loc);
		if (ro * R_loc[i] > R_glob[i]) {
			R[i] = ro * R_loc[i];
			*r_upd = r_loc;
		}
		else { R[i] = R_glob[i]; *r_upd = r_glob; }
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
void SolverNew::FindMinR(int* K, vector<Trial>& trials, int* t) { ////////// new
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
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2 - ((trials[*t].z - trials[*t - 1].z) / (2 * *m * r));
	NewTrial.z = obj_func(NewTrial.x,i);
	NewTrial.k = *K + 1;
	*K = *K + 1;
}
void SolverModif::MakeNewTrial(int i, int* K, vector<Trial>& trials, int* t, double* m, double* r_upd) {
	NewTrial.x = (trials[*t].x + trials[*t - 1].x) / 2 - ((trials[*t].z - trials[*t - 1].z) / (2 * *m * *r_upd));
	NewTrial.z = obj_func(NewTrial.x, i);
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
void Solver::CheckStopConditionCurrentZ(int i, bool* flag, vector<Trial>& trials, int* t) {
	if ((trials[*t].z <= eps))//условие остановки
		*flag = 0;
	trials.insert(trials.begin() + *t, NewTrial);
}
void Solver::MakeBestTrial(double* minFunc, double* argmin, double* l) {
	BestTrial.z = *minFunc;
	BestTrial.x = *argmin;
	BestTrial.k = *l;
}

///////////////////////////////////////////////////////

void File::InputToFile_percent(vector<int>&K, int n, int a) {
	fstream file;
	int p = 801;//Для Kmax=1000 100*8+1, если 1.25(для шекеля) / 101, если 10(для хилла)
	int Median = 0;
	vector <double> P(p);
	for (int i = 1; i < n; i++) {
		int j = i - 1;
		while (j >= 0 && K[j] > K[j + 1])
		{
			swap(K[j], K[j + 1]);
			j--;
		}
	}
	if (K.size() % 2 == 0) {
		Median = ((K[K.size() / 2 + 1]) + ((K[K.size() / 2 - 1]))) / 2;
	}
	else {
		Median = K[K.size() / 2 + 0.5];
	}
	int m = 0;
	int t = 0;
	bool f = 1;
	for (int i = 0; i < p; i++) {
		f = 1;
		m = 0;
		while (m < n && K[m] <= i * 1.25) {//1.25 //10
			P[i] = P[i] + 1;
			m++;
			f = 0;
		}
		if (f) P[i] = 0;
	}
	if (a == 1)	file.open("percBA.txt", fstream::in | fstream::out);
	else if (a == 2) file.open("percBA_M.txt", fstream::in | fstream::out);
	else if (a==3) file.open("percNewA.txt", fstream::in | fstream::out);
	else if (a==4) file.open("percRoot.txt", fstream::in | fstream::out);

	file << "K: ";
	for (int i = 0; i < p; i++) {
		if (i == 0) file << i * 1.25;//1.25 //10
		else file << ", " << i * 1.25;//1.25 //10
	}
	file << "\n" << "Median: " << Median;
	file << "\n" << "P(%): ";
	for (int i = 0; i < p; i++) {
		if (i == 0) file << P[i] * 100 / n;
		else file << ", " << P[i] * 100 / n;
	}
	file.close();
}
void File::InputToFile_x_solver (vector<Trial>& trials, int n, int a) {
	fstream file;
	THillProblemFamily hpf;
	TShekelProblemFamily shf;

	if (a == 1)	file.open("x_solverBA.txt", fstream::in | fstream::out);
	else if (a == 2) file.open("x_solverBA_M.txt", fstream::in | fstream::out);
	else if (a == 3) file.open("x_solverNewA.txt", fstream::in | fstream::out);
	else if (a == 4) file.open("x_solverRoot.txt", fstream::in | fstream::out);

	//vector <int> b(401);
	int g = 0;//-200
	for (int i = 0; i < 401; i++) {//801//401
		if (i == 0) file << "x: " << g * 0.05 / 2;
		else file << ", " << g * 0.05 / 2;
		g++;
	}
	g = 0;
	file << "\n\n";
	for (int i = 0; i < 401; i++) {
		if (i == 0) file << "phi: " << shf[n]->ComputeFunction({ g * 0.05 / 2 }) - shf[n]->GetOptimumValue();//(g * 0.005 / 2) * (g * 0.005 / 2) - cos(18 * (g * 0.005 / 2)) + sin(50 * (g * 0.005 / 2)) + 1.878638; //\\hpf[n]->ComputeFunction({g * 0.005 / 2}) - hpf[n]->GetOptimumValue() //\\shf[n]->ComputeFunction({g * 0.05 / 2}) - shf[n]->GetOptimumValue()
		else file << ", " << shf[n]->ComputeFunction({ g * 0.05 / 2 }) - shf[n]->GetOptimumValue();//(g * 0.005 / 2) * (g * 0.005 / 2) - cos(18 * (g * 0.005 / 2)) + sin(50 * (g * 0.005 / 2)) + 1.878638;
		g++;
	}
	file << "\n\n";
	for (int i = 0; i < trials.size(); i++) {
		if (i == 0) file << "x solver: " << trials[i].x;
		else file << ", " << trials[i].x;
	}
	file.close();
}
/* for specific func
 int g = -200;//0
	for (int i = 0; i < 801; i++) {//401
		if (i == 0) file << "x: " << g * 0.005 / 2;
		else file << ", " << g * 0.005 / 2;
		g++;
	}
	g =-200;
	file << "\n\n";
	for (int i = 0; i < 801; i++) {
		if (i == 0) file << "phi: " << (g * 0.005 / 2) * (g * 0.005 / 2) - cos(18 * (g * 0.005 / 2)) + sin(50 * (g * 0.005 / 2)) + 1.878638; //hpf[n]->ComputeFunction({g * 0.005 / 2}) - hpf[n]->GetOptimumValue();
		else file << ", " << (g * 0.005 / 2) * (g * 0.005 / 2) - cos(18 * (g * 0.005 / 2)) + sin(50 * (g * 0.005 / 2)) + 1.878638; //hpf[n]->ComputeFunction({g * 0.005 / 2}) - hpf[n]->GetOptimumValue();
		g++;
	}
*/

///////////

void Practice::BA(vector<int>& K, int n) {
	Solver s;
	THillProblemFamily hpf;
	TShekelProblemFamily shf;
	File f;
	vector<int> p(50);
	vector <int> l(n);
	static int u, t;
	double sumK = 0;
	for (int i = 0; i < n; i++) {
		s.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = s.GetKcur();
		sumK = sumK + K[i];
		l[i] = s.GetKcur();
		if (l[i] == 1000) {
			u++;
			p[t] = i;
			t++;
		}
	}
	cout << u << endl;
	for (int i = 0; i < t; i++) {
		cout << p[i] << endl;
	}
	cout << "K:" << sumK / 1000;
	f.InputToFile_x_solver(s.trials, n - 1, 1);
}
void Practice::BA_Modif(vector<int>& K, int n) {
	SolverModif sm;
	THillProblemFamily hpf;
	TShekelProblemFamily shf;
	File f;
	vector<int> p(20);
	vector <int> l(n);
	int u=0;
	int t = 0;
	double sumK = 0;
	for (int i = 0; i < n; i++) {
		sm.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = sm.GetKcur();
		sumK = sumK + K[i];
		l[i] = sm.GetKcur();
		if (l[i] == 1000) {
			u++;
			p[t] = i;
			t++;
		}
	}
	cout << u << endl;
	for (int i = 0; i < t; i++) {
		cout << p[i] << endl;
	}
	cout << "K:" << sumK/1000;

	f.InputToFile_x_solver(sm.trials, n - 1, 2);
}
void Practice::NewA(vector<int>& K, int n) {
	SolverNew sn;
	THillProblemFamily hpf;
	TShekelProblemFamily shf;
	File f;
	vector<int> p(10);
	vector <int> l(n);
	int u=0;
	int t = 0;
	double sumK = 0;
	for (int i = 0; i < n; i++) {
		sn.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = sn.GetKcur();
		sumK = sumK + K[i];
		l[i] = sn.GetKcur();
		if (l[i] == 1000) {
			u++;
			p[t] = i;
			t++;
		}
	}
	cout << u << endl;
	for (int i = 0; i < t; i++) {
		cout << p[i] << endl;
	}
	cout << "K:" << sumK / 1000;

	f.InputToFile_x_solver(sn.trials, n - 1, 3);
}
void Practice::RootA(vector<int>& K, int n) {
	SolverRoot sr;
	THillProblemFamily hpf;
	File f;
	vector<int> p(52);
	int u;
	int t = 0;
	for (int i = 0; i < n; i++) {
		sr.Solve(i);
		cout << "\nOptimumValue = " << hpf[i]->GetOptimumValue() << endl;
		cout << "OptimumPoint = " << hpf[i]->GetOptimumPoint()[0] << endl;
		//cout << "\nOptimumValue = " << shf[i]->GetOptimumValue() << endl;
		//cout << "OptimumPoint = " << shf[i]->GetOptimumPoint()[0] << endl;
		K[i] = sr.GetKcur();
	}

	f.InputToFile_x_solver(sr.trials, n - 1, 4);
}
