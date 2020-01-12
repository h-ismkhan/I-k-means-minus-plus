#ifndef _IGlobalKMeans_User_IMP_H
#define _IGlobalKMeans_User_IMP_H

#include <algorithm>
#include "I-KMeans-IMP.h"

//bool incOr(pair<int, double> i, pair<int, double> j) { return (i.second < j.second); }
//bool decOr(pair<int, double> i, pair<int, double> j) { return (i.second > j.second); }

class IGlobalKMeans_User_IMP : public IKMeans_IMP
{
protected:
	KMeans_OneReplacement_ManualSeeds *kmeans_global;
protected:
	void UpdateIsAdjacent(vector<vector<bool>> &IsInAdjacent)
	{
		int dim = this->inst->Dimension();
		int thek = kmeans_global->getNumberOfClusters();
		const int* clusIds = kmeans_global->getClusterIds();
		const int* secondClusIds = kmeans_global->getSecondCenterIndex();

		if (IsInAdjacent.size() < thek)
		{
			vector<bool> isInVec(thek, false);
			IsInAdjacent.assign(thek, isInVec);
		}
		else
		{
			for (int i = 0; i < thek; i++)
			for (int j = 0; j < thek; j++)
				IsInAdjacent[i][j] = false;
		}

		for (int i = 0; i < this->dim; i++)
			IsInAdjacent[clusIds[i]][secondClusIds[i]] = true;//NOTE: NOT IsInAdjacent[secondClusIds[i]][clusIds[i]] = true -->please see next line
		// ++++++++++++Ck++++++++++++            ++Cj++Ci++ ==>   Cj is adjacent of Ck but Ck is not adjacent Cj
	}
	void UpdateGainCost(vector<double> &Gain, vector<double> &Cost)
	{
		int dim = this->inst->Dimension();
		int thek = kmeans_global->getNumberOfClusters();
		const int* clusIds = kmeans_global->getClusterIds();

		if (Cost.size() < thek)
		{
			Cost.assign(thek, 0);
			Gain.assign(thek, 0);
		}
		else
		{
			for (int i = 0; i < thek; i++)
				Cost[i] = Gain[i] = 0;
		}

		const double *secondCandCenterDis = kmeans_global->getSecondCenterDis();
		const double *disFromCenter = kmeans_global->getDisFromItsCenter();
		for (int i = 0; i < dim; i++)
		{
			int clusId = clusIds[i];

			//we compute the current SSEDM of what potentially is to be added.
			//we then, only estimate what we gain using 3*SSEDM/4 (SSEDM - (1/4)SSEDM)
			double dis = disFromCenter[i];
			double disP2 = dis * dis;
			Gain[clusId] = Gain[clusId] + disP2;


			//dis =is= disFromCenter[i];
			//comute what we gain by deleting:
			Cost[clusId] = Cost[clusId] - disP2;
			//comute what we *MISS* by deleting:
			dis = secondCandCenterDis[i];
			Cost[clusId] = Cost[clusId] + dis * dis;
		}
	}
	void UpdateGainCost(vector<pair<int, double>> &Gain, vector<pair<int, double>> &Cost)
	{
		int dim = this->inst->Dimension();
		int thek = kmeans_global->getNumberOfClusters();
		const int* clusIds = kmeans_global->getClusterIds();

		if (Cost.size() < thek)
		{
			pair<int, double> empty(0, 0.0);
			Cost.assign(thek, empty);
			Gain.assign(thek, empty);
		}
		else
		{
			for (int i = 0; i < thek; i++)
			{
				Gain[i].first = Cost[i].first = i;
				Gain[i].second = Cost[i].second = 0;
			}
		}

		const double *secondCandCenterDis = kmeans_global->getSecondCenterDis();
		const double *disFromCenter = kmeans_global->getDisFromItsCenter();
		for (int i = 0; i < dim; i++)
		{
			int clusId = clusIds[i];

			//we compute the current SSEDM of what potentially is to be added.
			//we then, only estimate what we gain using 3*SSEDM/4 (SSEDM - (1/4)SSEDM)
			double dis = disFromCenter[i];
			double disP2 = dis * dis;
			Gain[clusId].second = Gain[clusId].second + disP2;


			//dis =is= disFromCenter[i];
			//comute what we gain by deleting:
			Cost[clusId].second = Cost[clusId].second - disP2;
			//comute what we *MISS* by deleting:
			dis = secondCandCenterDis[i];
			Cost[clusId].second = Cost[clusId].second + dis * dis;
		}
	}
protected:
	double SSEDM(IInstance_Euclidean *inst, IClustering_Ordinary *kmmanu, int theK)
	{
		double SSEDMV = 0;
		for (int i = 0; i < dim; i++)
		{
			double dis = kmmanu->DisFromItsCluster(i);
			SSEDMV = SSEDMV + dis * dis;
		}
		return SSEDMV;
	}
	double SEDM(IInstance_Euclidean *inst, IClustering_Ordinary *kmmanu, int theK)
	{
		double SEDMV = 0;
		for (int i = 0; i < dim; i++)
		{
			SEDMV = SEDMV + kmmanu->DisFromItsCluster(i);
		}
		return SEDMV;
	}	
public:
	IGlobalKMeans_User_IMP(){ kmeans_global = 0; }
	~IGlobalKMeans_User_IMP(){ if (kmeans_global) delete kmeans_global; }
	virtual double Center(int centerId, int j)
	{
		return kmeans_global->Center(centerId, j);
	}
	virtual int getClusterIdOf(int n)
	{
		return kmeans_global->getClusterIdOf(n);
	}
	virtual const int* getSizeOfClusters()
	{
		return kmeans_global->getSizeOfClusters();
	}
	virtual const double* getCenter(int i)
	{
		return kmeans_global->getCenter(i);
	}
	virtual const int* getClusterIds(){ return kmeans_global->getClusterIds(); }
};

#endif