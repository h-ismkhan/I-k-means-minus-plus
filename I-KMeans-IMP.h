#ifndef _ikmean_imp_H
#define _ikmean_imp_H

#include <algorithm>
#include "KM-OneReplcmnt-ManualSeeds.h"
#include "ISeed-Initializer.h"


class IKMeans_IMP : public IClustering_Ordinary
{
protected:
	ISeed_Initializer *seedInit = 0;
protected:
	class ISolution_KPoints
	{
	public:
		virtual int theK() = 0;
		virtual const int* Points() = 0;
		virtual double fitness() = 0;
		virtual void SetPoint(int kthPoint, int id) = 0;
		virtual IInstance_Euclidean* GetInstance() = 0;
	};
	class ISolution_KPoints_Ordinary : public ISolution_KPoints
	{
	protected:
		int num_of_clusters;
		vector<int> points;
		double fit;
		IInstance_Euclidean *inst;
		bool is_sol_updated;
		void updateFitness()
		{
			fit = 0;

			for (int i = 0; i < num_of_clusters; i++)
			{
				vector<double>disfrom;
				int r = points[i];
				for (int j = 0; j < num_of_clusters; j++)
				{
					if (j == i) continue;
					int s = points[j];
					disfrom.push_back(inst->Dis(r, s));
				}
				std::sort(disfrom.begin(), disfrom.end());
				fit = fit + disfrom[0];
				for (int j = 0; j < num_of_clusters - 1; j++)
					fit = fit + disfrom[j];
			}

			is_sol_updated = false;
		}
	public:
		ISolution_KPoints_Ordinary(IInstance_Euclidean *inst, int K)
		{
			this->inst = inst;
			this->num_of_clusters = K;
			is_sol_updated = true;
			points.assign(K, -1);
		}
		int theK(){ return num_of_clusters; }
		const int* Points(){ return points.data(); }
		void SetPoint(int kthPoint, int id){ points[kthPoint] = id; is_sol_updated = true; }
		IInstance_Euclidean* GetInstance(){ return inst; }
		double fitness(){ if (is_sol_updated)updateFitness(); return fit; }
	};
protected:
	virtual ISolution_KPoints_Ordinary* getNext(IInstance_Euclidean *inst, int thek, int start)
	{

		int dim = inst->Dimension();

		vector<int> nodes(dim);
		MakeVectorRandom(&nodes);

		ISolution_KPoints_Ordinary* outpt = new ISolution_KPoints_Ordinary(inst, thek);
		const int *k_1_selectedCenters = outpt->Points();
		vector<bool> isin(dim, false);
		outpt->SetPoint(0, start);
		isin[start] = true;
		
		int L = thek * 30;
		int nodescount = min(L, dim);
		
		vector<double>sumDisFrom_K_1_PreviousCenter(nodescount, 0.0);

		for (int i = 1; i < thek; i++)
		{
			int nextnode = -1;
			double nextweight = DBL_MIN;
			for (int z = 0; z < nodescount; z++)
			{
				int node = nodes[z];
				if (isin[node]) continue;/*{ L++;  }*/
				
				double newdisFromLastCenter = inst->Dis(node, k_1_selectedCenters[i - 1]);
				sumDisFrom_K_1_PreviousCenter[z] = sumDisFrom_K_1_PreviousCenter[z] + newdisFromLastCenter;
				double diff = abs(newdisFromLastCenter - sumDisFrom_K_1_PreviousCenter[z] / (double)i);
				
				if (diff == 0)diff = 1;
				double curweight = sumDisFrom_K_1_PreviousCenter[z] / diff;
				if (curweight > nextweight)
				{
					nextweight = curweight;
					nextnode = node;
				}
			}
			outpt->SetPoint(i, nextnode);
			isin[nextnode] = true;
		}
		return outpt;
	}
public:
	virtual void setSeedInitializer(ISeed_Initializer *inputSeedInit){ this->seedInit = inputSeedInit; }
	virtual void setNumberOfClusters(int K)
	{
		num_of_clusters = K;
	}
	virtual void setInstance(IInstance_Euclidean* IInstance)
	{
		this->inst = IInstance;
		this->dim = inst->Dimension();
		this->vecDim = inst->VecDim();
	}
	virtual int getNumberOfClusters()
	{
		return num_of_clusters;
	}
	virtual const vector<int>* Apply() = 0;
	virtual string Name() = 0;
};

#endif