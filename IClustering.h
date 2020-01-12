#ifndef _IClustering_H
#define _IClustering_H



#include "IInstance-Euclidean.h"

class IClustering
{
public:
	virtual int getNumberOfClusters() = 0;
	virtual void setInstance(IInstance_Euclidean* IInstance) = 0;
	virtual const vector<int>* Apply() = 0;
	virtual double Center(int centerId, int j) = 0;
	virtual int getClusterIdOf(int) = 0;
	virtual const int* getClusterIds() = 0;

	virtual void setNumberOfClusters(int) = 0;
	virtual void setParm(int index, double value) = 0;

	virtual string Name() = 0;
};

#endif