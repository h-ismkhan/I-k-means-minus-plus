#ifndef _IInstance_Euclidean_Ordinary_H
#define _IInstance_Euclidean_Ordinary_H

#include <vector>
#include "IInstance-Euclidean.h"

class IInstance_Euclidean_Ordinary : public IInstance_Euclidean
{
protected:
	int dim, vecDim, actualK;
	vector<double> maxvector, minvector, centervector;
	mutable vector<vector<double>> vectors;
	vector<int> actualClusIds;
	
public:
	virtual int Dimension() const
	{
		return dim;
	}
	virtual int VecDim() const
	{
		return vecDim;
	}
	
	virtual double Dis(int v1index, int v2index) const
	{
		double dis = 0;
		for (int i = 0; i < vecDim; i++)
		{
			double sub = vectors[v1index][i] - vectors[v2index][i];
			dis = dis + sub * sub;
		}
		return sqrt(dis);
	}
	virtual double Dis(const double *vec_1, const double *vec_2) const
	{
		double dis = 0;
		for (int i = 0; i < vecDim; i++)
		{
			double sub = vec_1[i] - vec_2[i];
			dis = dis + sub * sub;
		}
		return sqrt(dis);
	}
	virtual double Dis(const vector<double> *vec_1, const vector<double> *vec_2) const
	{
		double dis = 0;
		for (int i = 0; i < vecDim; i++)
		{
			double sub = (*vec_1)[i] - (*vec_2)[i];
			dis = dis + sub * sub;
		}
		return sqrt(dis);
	}
	virtual double Dis(int o, const vector<double> *vec_2) const
	{
		return Dis(&vectors[o], vec_2);
	}
	virtual double Dis(int o, const double *vec_2) const
	{
		return Dis(vectors[o].data(), vec_2);
	}

	virtual bool HasEdge(int x, int y) const { return 1; }
	virtual double get(int node, int kthCoord) const
	{
		return vectors[node][kthCoord];
	}
	inline const double* operator[](int node) const
	{
		return vectors[node].data();
	}
	virtual const double* get(int node) const
	{
		return vectors[node].data();
	}
	virtual const vector<double>* getVec(int node) const
	{
		return &vectors[node];
	}

	virtual const int* getActualClustersIds() const{ return actualClusIds.data(); }
	virtual int getActualNumOfClusters() const{ return actualK; }

	virtual double Center(int index) const
	{
		return centervector[index];
	}
	virtual double Max(int index) const
	{
		return maxvector[index];
	}
	virtual double Min(int index) const
	{
		return minvector[index];
	}
	virtual void UpdateMinMax()
	{		
		maxvector.assign(vecDim, 0.0);
		minvector.assign(vecDim, 0.0);

		for (int i = 0; i < vecDim; i++)
		{
			maxvector[i] = minvector[i] = vectors[0][i];
		}
		for (int i = 1; i < dim; i++)
		for (int j = 0; j < vecDim; j++)
		{
			if (maxvector[j] < vectors[i][j])
				maxvector[j] = vectors[i][j];
			if (minvector[j] > vectors[i][j])
				minvector[j] = vectors[i][j];
		}
	}
	virtual void UpdateCenter()
	{
		centervector.assign(vecDim, 0.0);

		for (int i = 0; i < vecDim; i++)
		{
			centervector[i] = vectors[0][i];
		}
		for (int i = 1; i < dim; i++)
		for (int j = 0; j < vecDim; j++)
		{			
			centervector[j] = centervector[j] + vectors[i][j];
		}
		for (int i = 0; i < vecDim; i++)
		{
			centervector[i] = centervector[i] / (double)dim;
		}
	}
	void Normalize()
	{
		if (minvector.size() == 0)
			UpdateMinMax();

		double min_new = 0.0, max_new = 2.0;
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < vecDim; j++)
			{
				double min_i = minvector[j], max_i = maxvector[j];
				if (min_i == max_i)
				{
					vectors[i][j] = (min_new + max_new) / 2;
					continue;
				}
				vectors[i][j] = min_new + (((vectors[i][j] - min_i) / (max_i - min_i)) * (max_new - min_new));
			}
		}		
	}
	IInstance_Euclidean_Ordinary* Copy()
	{
		IInstance_Euclidean_Ordinary *toret = new IInstance_Euclidean_Ordinary();
		toret->actualClusIds = this->actualClusIds;
		toret->actualK = this->actualK;
		toret->centervector = this->centervector;
		toret->dim = this->dim;
		toret->maxvector = this->maxvector;
		toret->minvector = this->minvector;
		toret->vecDim = this->vecDim;
		toret->vectors = this->vectors;
		return toret;
	}
};

#endif