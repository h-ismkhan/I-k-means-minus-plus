#ifndef _kmean_OneReplcmnt_ManualSeeds_H
#define _kmean_OneReplcmnt_ManualSeeds_H

#include "IClustering-Ordinary.h"
#include "ISeed-Initializer.h"
#include <assert.h>

class KMeans_OneReplacement_ManualSeeds : public IClustering_Ordinary
{
private:
	ISeed_Initializer *seedInitializer;

	vector <bool> isCenterChanged;//it is used after apply function, so->
	//it can be initialized in the first of each other function, ofcourcse after its->
	//first utilization (its previous state may be used in the first of each function)

	//vector<int> changedCenterList;

	vector<int> secondCenterCandIndex;
	vector<double> secondCenterCandDis;

	int addedCenter = -1, replacedCenter = -1;
private:
	virtual void AssignCentersMem()
	{
		ccenters.assign(num_of_clusters, vector<double>(vecDim));
	}
	virtual void InitializeCcenters()
	{
		for (int i = 0; i < num_of_clusters; i++)
		{
			for (int j = 0; j < vecDim; j++)
				ccenters[i][j] = Random(inst->Min(j), inst->Max(j));
		}
	}

	virtual void UpdateCenters()
	{
		//this function uses this->isCenterChanged
		//if a center is changed in UpdateId, then it is participated here.
		for (int i = 0; i < num_of_clusters; i++)
		{
			if (!this->isCenterChanged[i])
				continue;
			SizeOf[i] = 0;
			for (int j = 0; j < vecDim; j++)
				ccenters[i][j] = 0;
		}
		for (int i = 0; i < dim; i++)
		{
			if (!this->isCenterChanged[this->cluster_id[i]])
				continue;

			SizeOf[cluster_id[i]]++;
			int c_id = cluster_id[i];
			const double* row = inst->get(i);
			for (int j = 0; j < vecDim; j++)
				ccenters[c_id][j] = ccenters[c_id][j] + row[j];
		}
		for (int i = 0; i < num_of_clusters; i++)
		{
			if (!this->isCenterChanged[i])
				continue;

			if (SizeOf[i])
			{
				for (int j = 0; j < vecDim; j++)
					ccenters[i][j] = ccenters[i][j] / SizeOf[i];
			}
		}
	}


	virtual void DeterminePotentialNearestCentersOf(int centerX, vector<int>& listOf_targetCenters, vector<int>& listOf_affectedPoints)
	{
		vector<bool> isCenterConsidered(this->num_of_clusters, false);
		isCenterConsidered[centerX] = true;
		int size = listOf_targetCenters.size();
		for (int i = 0; i < size; i++)
			isCenterConsidered[listOf_targetCenters[i]] = true;

		vector<bool> isPointAdded(this->dim, false);
		size = listOf_affectedPoints.size();
		for (int i = 0; i < size; i++)
			isPointAdded[listOf_affectedPoints[i]] = true;

		for (int i = 0; i < dim; i++)
		{
			if ((secondCenterCandIndex[i] == centerX || cluster_id[i] == centerX) && !isPointAdded[i])
			{
				listOf_affectedPoints.push_back(i);
				isPointAdded[i] = true;
			}
		}

		size = listOf_affectedPoints.size();
		for (int i = 0; i < size; i++)
		{
			int point = listOf_affectedPoints[i];

			int firstCenterPoint = cluster_id[point];
			int secondCenterPoint = this->secondCenterCandIndex[point];

			if (!isCenterConsidered[firstCenterPoint])
			{
				listOf_targetCenters.push_back(firstCenterPoint);
				isCenterConsidered[firstCenterPoint] = true;

			}
			if (!isCenterConsidered[secondCenterPoint])
			{
				listOf_targetCenters.push_back(secondCenterPoint);
				isCenterConsidered[secondCenterPoint] = true;

			}
		}
	}
	virtual short UpdateClustersID_UseActiveCenter()
	{	
		//determine affected points. affected points={points with active first and second NC}
		vector<int> affectePoints;
		vector<int> activeCenters;
		for (int i = 0; i < this->dim; i++)
		{
			if (this->isCenterChanged[cluster_id[i]]
				||
				this->isCenterChanged[secondCenterCandIndex[i]])
				affectePoints.push_back(i);
		}
		for (int k = 0; k < this->num_of_clusters; k++){
			if (!this->isCenterChanged[k]) continue;
			activeCenters.push_back(k);
		}

		vector<bool> backup_isCenterChanged = isCenterChanged;
		this->isCenterChanged.assign(this->num_of_clusters, false);//it sholud be initialize by false

		short IsChanged = 0;
		int size = affectePoints.size();
		for (int pi = 0; pi < size; pi++)
		{
			int point = affectePoints[pi];
			int pointOldId = this->cluster_id[point];
			
			if (backup_isCenterChanged[this->secondCenterCandIndex[point]])
			{
				this->secondCenterCandIndex[point] = -1;
				this->secondCenterCandDis[point] = DBL_MAX;
			}
			if (backup_isCenterChanged[this->cluster_id[point]])
			{
				this->cluster_id[point] = this->secondCenterCandIndex[point];
				this->firstCenterCandDis_or_clusterDis[point] = this->secondCenterCandDis[point];

				this->secondCenterCandIndex[point] = -1;
				this->secondCenterCandDis[point] = DBL_MAX;
			}

			double min_dis = this->firstCenterCandDis_or_clusterDis[point];
			int index = this->cluster_id[point];
			int sizeActiveCenters = activeCenters.size();
			for (int cj = 0; cj < sizeActiveCenters; cj++)
			{
				int centerActive = activeCenters[cj];
				double cur_dis = inst->Dis(point, &(ccenters[centerActive]));

				if (cur_dis < min_dis)
				{
					secondCenterCandDis[point] = min_dis;
					secondCenterCandIndex[point] = index;

					min_dis = cur_dis;
					index = centerActive;
				}
				else if (cur_dis < secondCenterCandDis[point])
				{
					secondCenterCandDis[point] = cur_dis;
					secondCenterCandIndex[point] = centerActive;
				}
			}

			firstCenterCandDis_or_clusterDis[point] = min_dis;
			cluster_id[point] = index;

			if (index != pointOldId)
			{
				IsChanged = 1;
				this->isCenterChanged[pointOldId] = this->isCenterChanged[index] = true;
			}
		}

		//==================================================================================
		//center c of a cluster is changed, so it may become first center of other points (points with second center c and unchanged cluster)
		for (int i = 0; i < dim; i++)
		{
			if (this->cluster_id[i] == this->addedCenter)
				this->isCenterChanged[this->secondCenterCandIndex[i]] = true;
		}

		return IsChanged;
	}
	virtual short UpdateClustersID_OneAddRemove()
	{	

		//===========================================================================================================================
		//the first step update first and second centers of affected points ****with deleted*** corresponding center.
		vector<int> deletedCenter_affectePoints;
		vector<int> deletedCenter_potentiaNearestCenters;
		DeterminePotentialNearestCentersOf(this->replacedCenter, deletedCenter_potentiaNearestCenters, deletedCenter_affectePoints);

		int theSize = deletedCenter_affectePoints.size();
		for (int pi = 0; pi < theSize; pi++)
		{
			int point = deletedCenter_affectePoints[pi];
			if (this->cluster_id[point] == replacedCenter)
			{
				this->isCenterChanged[this->secondCenterCandIndex[point]] = true;
				this->cluster_id[point] = this->secondCenterCandIndex[point];
				this->firstCenterCandDis_or_clusterDis[point] = this->secondCenterCandDis[point];
			}

			this->secondCenterCandDis[point] = DBL_MAX;
			this->secondCenterCandIndex[point] = -1;
		}

		deletedCenter_potentiaNearestCenters.push_back(replacedCenter);//->next line
		//some points of cluster with deleted centers may be affected by updated center-replacedCenter-.
		theSize = deletedCenter_affectePoints.size();
		for (int pi = 0; pi < theSize; pi++)
		{
			int point = deletedCenter_affectePoints[pi];
			int pointOldId = this->cluster_id[point];
			//firstCenterCandDis_or_clusterDis[point] = secondCenterCandDis[point] = DBL_MAX;
						
			int theSize2 = deletedCenter_potentiaNearestCenters.size();
			for (int cj = 0; cj < theSize2; cj++)
			{
				int potentialCenter = deletedCenter_potentiaNearestCenters[cj];
				double cur_dis = inst->Dis(point, &(ccenters[potentialCenter]));

				if (cur_dis < this->firstCenterCandDis_or_clusterDis[point])
				{
					this->secondCenterCandDis[point] = this->firstCenterCandDis_or_clusterDis[point];
					this->secondCenterCandIndex[point] = this->cluster_id[point];

					this->firstCenterCandDis_or_clusterDis[point] = cur_dis;
					this->cluster_id[point] = potentialCenter;
				}
				else if (cur_dis < secondCenterCandDis[point])
				{
					secondCenterCandDis[point] = cur_dis;
					secondCenterCandIndex[point] = potentialCenter;
				}
			}			

			if (pointOldId != cluster_id[point])
			{
				this->isCenterChanged[pointOldId] = this->isCenterChanged[cluster_id[point]] = true;
			}
		}

		//===========================================================================================================================
		//the second step: update first and second centers of affected points ****with deleted*** corresponding center.
		vector<int> addedCenter_affectePoints;
		vector<int> addedCenter_potentiaNearestCenters;
		DeterminePotentialNearestCentersOf(this->addedCenter, addedCenter_potentiaNearestCenters, addedCenter_affectePoints);

		addedCenter_potentiaNearestCenters.push_back(this->replacedCenter);//->next line
		//some points of points with added center may be affected by updated center-replacedCenter-.
		theSize = addedCenter_affectePoints.size();
		for (int pi = 0; pi < theSize; pi++)
		{
			int point = addedCenter_affectePoints[pi];
			int pointOldId = this->cluster_id[point];

			int theSize2 = addedCenter_potentiaNearestCenters.size();
			for (int cj = 0; cj < theSize2; cj++)
			{
				int potentialCenter = addedCenter_potentiaNearestCenters[cj];
				double cur_dis = inst->Dis(point, &(ccenters[potentialCenter]));

				if (cur_dis < this->firstCenterCandDis_or_clusterDis[point])
				{
					this->secondCenterCandDis[point] = this->firstCenterCandDis_or_clusterDis[point];
					this->secondCenterCandIndex[point] = this->cluster_id[point];

					this->firstCenterCandDis_or_clusterDis[point] = cur_dis;
					this->cluster_id[point] = potentialCenter;
				}
				else if (cur_dis < secondCenterCandDis[point])
				{
					secondCenterCandDis[point] = cur_dis;
					secondCenterCandIndex[point] = potentialCenter;
				}
			}
			if (pointOldId != cluster_id[point])
			{
				this->isCenterChanged[pointOldId] = this->isCenterChanged[cluster_id[point]] = true;
			}
		}

		//==================================================================================
		//center c of split cluster will changed, so it may become first center of other points (points whch thier second was c)
		for (int i = 0; i < dim; i++)
		{
			if (this->cluster_id[i] == this->addedCenter)
				this->isCenterChanged[this->secondCenterCandIndex[i]] = true;
		}

		//========================================================================================
				

		this->addedCenter = -1; this->replacedCenter = -1;
		//========================================================================================

		short ischanged = 0;
		for (int i = 0; i < this->num_of_clusters; i++)
		{
			if (this->isCenterChanged[i])
			{
				ischanged = true;
				break;
			}
		}
		
		return ischanged;
	}
	virtual short UpdateClustersID()
	{
		if (addedCenter >= 0 && replacedCenter >= 0)
			return UpdateClustersID_OneAddRemove();

		if (this->isCenterChanged.size() == this->num_of_clusters)
		{
			for (int j = 0; j < num_of_clusters; j++)
			{
				if (!isCenterChanged[j])
					return UpdateClustersID_UseActiveCenter();
			}
		}
		this->isCenterChanged.assign(this->num_of_clusters, false);//it sholud be initialize by false

		short IsChanged = 0;
		for (int i = 0; i < dim; i++)
		{
			double min_dis = DBL_MAX;
			int index = -1;
			for (int j = 0; j < num_of_clusters; j++)
			{
				double cur_dis = inst->Dis(i, &(ccenters[j]));
				
				if (cur_dis < min_dis)
				{
					secondCenterCandDis[i] = min_dis;
					secondCenterCandIndex[i] = index;

					min_dis = cur_dis;
					index = j;
				}
				else if (cur_dis < secondCenterCandDis[i])
				{
					secondCenterCandDis[i] = cur_dis;
					secondCenterCandIndex[i] = j;
				}
			}

			firstCenterCandDis_or_clusterDis[i] = min_dis;

			if (index != cluster_id[i])
			{
				IsChanged = 1;
				this->isCenterChanged[index] = true;
				if (cluster_id[i] >= 0)
					this->isCenterChanged[cluster_id[i]] = true;
				cluster_id[i] = index;
			}
		}
		return IsChanged;
	}	
public:
	KMeans_OneReplacement_ManualSeeds()
	{
		addedCenter = -1; replacedCenter = -1;

		seedInitializer = 0;
		num_of_clusters = 0;
	}
	~KMeans_OneReplacement_ManualSeeds()
	{
	}
	virtual void setInstance(IInstance_Euclidean* IInstance)
	{
		IClustering_Ordinary::setInstance(IInstance);
		secondCenterCandDis.assign(dim, DBL_MAX);
		secondCenterCandIndex.assign(dim, -1);
		this->isCenterChanged.clear();
	}

	int getSecondCenterIndex(int id)
	{
		return this->secondCenterCandIndex[id];
	}
	double getSecondCenterDis(int id)
	{
		return this->secondCenterCandDis[id];
	}
	const int*getSecondCenterIndex()
	{
		return this->secondCenterCandIndex.data();
	}
	const double*getSecondCenterDis()
	{
		return this->secondCenterCandDis.data();
	}

	void setSeedInitializer(ISeed_Initializer *seedInitializer)
	{
		this->seedInitializer = seedInitializer;
	}
	void setSeeds(IInstance_Euclidean* inst, const int* points)
	{
		vecDim = inst->VecDim();
		ccenters.assign(num_of_clusters, vector<double>(vecDim));

		this->isCenterChanged.assign(this->num_of_clusters, true);
		
		for (int i = 0; i < num_of_clusters; i++)
		{
			const double* vec = inst->get(points[i]);
			for (int j = 0; j < vecDim; j++)
				ccenters[i][j] = vec[j];
		}
	}
	virtual void setNumberOfClusters(int K)
	{
		IClustering_Ordinary::setNumberOfClusters(K);
		isCenterChanged.assign(K, true);
	}
	void setCenterAsChanged_KM_ManualSeed(int cid)
	{
		isCenterChanged[cid] = true;
	}
	void addCenter_KM_ManualSeed(const double* point, int vecdim)
	{
		if (ccenters.size() > 0)
			assert(this->vecDim == vecdim);
		else this->vecDim = vecdim;

		ccenters.push_back(vector<double>(vecDim));
		num_of_clusters++;
		for (int i = 0; i < vecDim; i++)
		{
			ccenters[num_of_clusters - 1][i] = point[i];
		}
		isCenterChanged.push_back(true);
	}
	void setCenter(int cindex, int pointId)
	{
		assert(replacedCenter < 0 && addedCenter < 0);

		replacedCenter = cindex;
		addedCenter = cluster_id[pointId];

		/*this->changedCenterList.clear();
		this->changedCenterList.push_back(replacedCenter);
		this->changedCenterList.push_back(addedCenter);*/
		/*this->isCenterChanged.assign(this->num_of_clusters, false);
		isCenterChanged[replacedCenter] = isCenterChanged[addedCenter] = true;*/

		const double*point = inst->get(pointId);
		if (ccenters.size() != num_of_clusters)
			ccenters.assign(num_of_clusters, vector<double>(vecDim));

		for (int i = 0; i < vecDim; i++)
		{
			ccenters[cindex][i] = point[i];
		}		
	}
	void setCenter_KM_ManualSeed(int cindex, const double* point)
	{
		if (ccenters.size() != num_of_clusters)
			ccenters.assign(num_of_clusters, vector<double>(vecDim));

		for (int i = 0; i < vecDim; i++)
		{
			ccenters[cindex][i] = point[i];
		}
		isCenterChanged[cindex] = true;
	}

	virtual vector<int>* Apply()
	{
		if (this->seedInitializer)
		{
			this->seedInitializer->setInstance(this->inst);
			this->seedInitializer->setTheK(this->num_of_clusters);
			this->setSeeds(this->inst, this->seedInitializer->nextKPoints());
		}

		while (1)
		{	
			if (!UpdateClustersID())
				break;
			UpdateCenters();
		}		

		replacedCenter = -1;
		addedCenter = -1;
		return &cluster_id;
	}
	virtual string Name()
	{
		if (this->seedInitializer)
			return "KMeans-" + seedInitializer->name();
		return "KMeans-ManualSeeds";
	}
	
	KMeans_OneReplacement_ManualSeeds* getCopy()
	{
		KMeans_OneReplacement_ManualSeeds *toret = new KMeans_OneReplacement_ManualSeeds();

		toret->dim = this->dim;
		toret->vecDim = this->vecDim;
		toret->num_of_clusters = this->num_of_clusters;
		toret->recent_num_of_clusters = this->recent_num_of_clusters;
		toret->inst = this->inst;
		
		toret->ccenters = this->ccenters;
		toret->cluster_id = this->cluster_id;
		toret->SizeOf = this->SizeOf;
		toret->isCenterChanged = this->isCenterChanged;
		
		toret->firstCenterCandDis_or_clusterDis = this->firstCenterCandDis_or_clusterDis;

		toret->secondCenterCandDis = this->secondCenterCandDis;
		toret->secondCenterCandIndex = this->secondCenterCandIndex;

		toret->addedCenter = this->addedCenter;
		toret->replacedCenter = this->replacedCenter;
		
		return toret;
	}
};

#endif