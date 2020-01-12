#ifndef _IInstance_H
#define _IInstance_H

#include <vector>
#include <iostream>
using namespace std;

class IInstance_Euclidean
{
public:	
	virtual int Dimension() const = 0;
	virtual int VecDim() const = 0;

	virtual double Dis(int node_1, int node_2) const = 0;
	virtual double Dis(const double *vec_1, const double *vec_2) const = 0;
	virtual double Dis(const vector<double> *vec_1, const vector<double> *vec_2) const = 0;
	virtual double Dis(int, const vector<double> *vec_2) const = 0;
	virtual double Dis(int, const double *vec_2) const = 0;
	virtual const double* get(int node) const = 0;
	virtual const vector<double>* getVec(int node) const = 0;
	virtual double get(int node, int kthCoord) const = 0;
	virtual inline const double* operator[](int node) const = 0;

	virtual double Center(int index) const = 0;
	virtual double Max(int index) const = 0;
	virtual double Min(int index) const = 0;
};

#endif