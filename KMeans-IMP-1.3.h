#ifndef _kmean_imp_1v3_H
#define _kmean_imp_1v3_H

#include "I-GlobalKMeans-User-IMP.h"


class KMeans_IMP_1v3 : public IGlobalKMeans_User_IMP
{
public:
	virtual const vector<int>* Apply()
	{
		int thek = num_of_clusters;
		int dim = inst->Dimension();
		int vecdim = inst->VecDim();


		if (kmeans_global)
			delete kmeans_global;

		kmeans_global = new KMeans_OneReplacement_ManualSeeds();

		kmeans_global->setInstance(inst);
		kmeans_global->setNumberOfClusters(thek);

		if (this->seedInit == 0)
		{
			ISolution_KPoints_Ordinary *InitialSol = getNext(inst, thek, Random(0, dim - 1));
			kmeans_global->setSeeds(inst, InitialSol->Points());

			delete InitialSol;
		}
		else
		{
			this->seedInit->setInstance(inst);
			this->seedInit->setTheK(thek);
			kmeans_global->setSeeds(inst, this->seedInit->nextKPoints());
		}


		kmeans_global->Apply();
		double curSSEDM_Global = SSEDM(inst, kmeans_global, thek);

		vector<bool> ToBeAvoided_Delete(thek, false);
		vector<bool> ToBeAvoided_Add(thek, false);

		vector<double> Cost, Gain;
		this->UpdateGainCost(Gain, Cost);

		vector<vector<bool>> isInAdjacent;
		this->UpdateIsAdjacent(isInAdjacent);

		vector<bool> toBePairVec(thek, true);
		vector<vector<bool>> canBePair(thek, toBePairVec);

		vector<int> failAdd(thek, 0);
		vector<int> failDel(thek, 0);

		int numOfSuc = 0;
		while (true)
		{
			int Ci = -1, Cj = -1;
			double largestGain = 0;
			for (int i = 0; i < thek; i++)
			{
				if (Gain[i] > largestGain && !ToBeAvoided_Add[i])
				{
					Ci = i;
					largestGain = Gain[i];
				}
			}
			if (Ci < 0) break;
			int rank = 0;
			for (int i = 0; i < thek; i++)
			{
				if (Gain[i] > Gain[Ci])
					rank++;
			}
			if (rank > thek / 2) break;
			double smallestCost = DBL_MAX;
			for (int j = 0; j < thek; j++)
			{
				if (Cost[j] < smallestCost && !ToBeAvoided_Delete[j] && !isInAdjacent[j][Ci] && !isInAdjacent[Ci][j] && j != Ci && canBePair[Ci][j])
				{
					Cj = j;
					smallestCost = Cost[j];
				}
			}
			if (Cj < 0) break;
			rank = 0;
			for (int j = 0; j < thek; j++)
			{
				if (Cost[j] < Cost[Cj] && !isInAdjacent[j][Ci] && !isInAdjacent[Ci][j] && j != Ci && canBePair[Ci][j])
					rank++;
			}
			if (rank > thek / 2) break;
			rank = 0;
			for (int i = 0; i < thek; i++)
			{
				if (Cost[i] < Cost[Cj])
					rank++;
			}
			if (rank > thek / 2) { ToBeAvoided_Add[Ci] = true; continue; }

			if (3 * Gain[Ci] / 4 - Cost[Cj] > 0)
			{
				vector<int> fromSi;
				const int *clusIds = this->kmeans_global->getClusterIds();
				for (int i = 0; i < dim; i++)
				{
					if (clusIds[i] == Ci)
						fromSi.push_back(i);
				}
				int nodeToAdd = fromSi[Random(0, fromSi.size() - 1)];
				KMeans_OneReplacement_ManualSeeds *candOf_kmeans_global = kmeans_global->getCopy();
				candOf_kmeans_global->setCenter(Cj, nodeToAdd);

				candOf_kmeans_global->Apply();
				double cand_SSEDM_Global = SSEDM(inst, candOf_kmeans_global, thek);
				if (cand_SSEDM_Global < curSSEDM_Global)
				{
					curSSEDM_Global = cand_SSEDM_Global;
					delete kmeans_global;
					kmeans_global = candOf_kmeans_global;
					clusIds = kmeans_global->getClusterIds();

					ToBeAvoided_Delete[Ci] = ToBeAvoided_Delete[Cj] = true;

					//previous (before updating) adjacents of Cj shoud be mareked NOT to be Added
					for (int i = 0; i < thek; i++)
					{
						if (isInAdjacent[Cj][i] && isInAdjacent[i][Cj])
							ToBeAvoided_Add[i] = true;
					}

					this->UpdateGainCost(Gain, Cost);
					this->UpdateIsAdjacent(isInAdjacent);

					//current (after updating) adjacents of Ci and Cj shoud be mareked NOT to be Deleted
					for (int i = 0; i < thek; i++)
					{
						if ((isInAdjacent[Ci][i] && isInAdjacent[i][Ci]) || (isInAdjacent[Cj][i] && isInAdjacent[i][Cj]))
							ToBeAvoided_Delete[i] = true;
					}

					numOfSuc++;
				}
				else
				{					
					canBePair[Ci][Cj] = canBePair[Cj][Ci] = false;

					failAdd[Ci]++;
					failDel[Cj]++;
					if (failAdd[Ci] > thek / 2)
						ToBeAvoided_Add[Ci] = true;
					if (failDel[Cj] > thek / 2)
						ToBeAvoided_Delete[Cj] = true;

					delete candOf_kmeans_global;
				}
			}
			else
				break;
			if (numOfSuc > thek / 2)
				break;
		}
		const vector<int> toret(kmeans_global->getClusterIds(), kmeans_global->getClusterIds() + thek);
		return &toret;
	}
	virtual string Name()
	{
		string name = "KM-IMP-1v3";
		if (this->seedInit != 0)
			name = name + "-Init: " + this->seedInit->name();
		return name;
	}
};

#endif