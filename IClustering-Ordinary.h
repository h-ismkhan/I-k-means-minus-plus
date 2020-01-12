#ifndef _IClustering_Ordinary_H
#define _IClustering_Ordinary_H



#include <iostream>
using namespace std;

#include "IInstance-Euclidean.h"
#include "IClustering.h"
#include "Time_Random.h"


bool incOr(pair<int, double> i, pair<int, double> j) { return (i.second < j.second); }
bool decOr(pair<int, double> i, pair<int, double> j) { return (i.second > j.second); }

int BinaryLargestSmall(const vector<pair<int, double>> &data, double c)
{
	int f = 0, e = data.size() - 1;
	while (e > f)
	{
		int mid = (f + e) / 2;
		if (c == data[mid].second)
			return mid - 1;
		else if (c < data[mid].second)
			e = mid - 1;
		else f = mid + 1;
	}
	if (e < 0) return -1;
	if (c > data[e].second)
		return e;
	if (c < data[f].second)
		return f - 1;
	return f;
}

class IClustering_Ordinary : public IClustering
{
protected:
	int vecDim, dim, num_of_clusters;
	IInstance_Euclidean *inst;
	int recent_num_of_clusters;
	vector<vector<double>> ccenters;

	vector<int> cluster_id;
	//vector<int> firstCenterCandIndex;
	vector<double> firstCenterCandDis_or_clusterDis;

	vector<int> SizeOf;
public:
	IClustering_Ordinary()
	{
		dim = vecDim = num_of_clusters = 0;
		inst = 0;
	}
	virtual void setInstance(IInstance_Euclidean* IInstance)
	{
		this->inst = IInstance;
		this->dim = inst->Dimension();
		this->vecDim = inst->VecDim();

		this->cluster_id.assign(dim, -1);
		this->firstCenterCandDis_or_clusterDis.assign(dim, DBL_MAX);
	}
	virtual double Center(int centerId, int j)
	{
		return ccenters[centerId][j];
	}
	virtual void setNumberOfClusters(int K)
	{
		SizeOf.assign(K, 0);
		num_of_clusters = K;
	}
	virtual int getNumberOfClusters()
	{
		return num_of_clusters;
	}
	virtual int getClusterIdOf(int n)
	{
		return cluster_id[n];
	}
	virtual const int* getClusterIds(){ return cluster_id.data(); }

	
	const double*getDisFromItsCenter()
	{
		return this->firstCenterCandDis_or_clusterDis.data();
	}
	virtual double DisFromItsCluster(int id)
	{
		return this->firstCenterCandDis_or_clusterDis[id];
	}

	virtual const int* getSizeOfClusters()
	{
		return SizeOf.data();
	}
	virtual const double* getCenter(int i)
	{
		return ccenters[i].data();
	}
	virtual void setParm(int index, double value){/*do no thing*/ }
};

#endif