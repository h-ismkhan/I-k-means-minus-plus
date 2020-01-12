#include "My-Instance.h"
#include "KMeans-IMP-1.3.h"
#include "Mult-UNCs-ByRateAveMax-FFC.h"

#include <conio.h>

#include <math.h>
class SSEDM_Measure
{
protected:
	IInstance_Euclidean *inst;
	const int*clusid;
public:
	virtual void setInstance(IInstance_Euclidean *input)
	{
		this->inst = input;
	}
	virtual void setClusterIds(const vector<int>* clusterIds)
	{
		this->clusid = clusterIds->data();
	}
	virtual void setClusterIds(const int* clusterIds)
	{
		this->clusid = clusterIds;
	}
	virtual void setClassIds(const vector<int>*)
	{
	}
	virtual double measure()
	{
		int dim = inst->Dimension();
		int vecDim = inst->VecDim();

		int num_of_clusters = 0;
		for (int i = 0; i < dim; i++)
		{
			if (clusid[i] > num_of_clusters)
				num_of_clusters = clusid[i];
		}
		num_of_clusters++;

		vector<vector<double>> ccenters;
		ccenters.assign(num_of_clusters, vector<double>(vecDim));
		for (int i = 0; i < num_of_clusters; i++)
		{
			for (int j = 0; j < vecDim; j++)
				ccenters[i][j] = 0;
		}
		vector<int>SizeOf(dim, 0);
		for (int i = 0; i < dim; i++)
		{
			SizeOf[clusid[i]]++;
			int c_id = clusid[i];
			const double* row = inst->get(i);
			for (int j = 0; j < vecDim; j++)
				ccenters[c_id][j] = ccenters[c_id][j] + row[j];
		}
		for (int i = 0; i < num_of_clusters; i++)
		{
			if (SizeOf[i])
			{
				for (int j = 0; j < vecDim; j++)
					ccenters[i][j] = ccenters[i][j] / SizeOf[i];
			}
		}

		double SSEDM = 0;
		for (int i = 0; i < dim; i++)
		{
			const double *vec = inst->get(i);
			double dis = inst->Dis(i, ccenters[clusid[i]].data());
			SSEDM = SSEDM + dis * dis;
		}
		return SSEDM;
	}
};

void main()
{
	ifstream in("path.txt");
	string path;

	in >> path;
	IInstance_Euclidean *inst = new My_Instance(path);

	cout << "Read Inst!" << endl;

	int thek;
	in >> thek;
	cout << "the k is read:" << thek << endl;

	KMeans_IMP_1v3 km;
	km.setInstance(inst);
	km.setNumberOfClusters(thek);
	LogMultUNCDis_by_RateUNCDisAveToMax_FFC_seedInit multiplyLogFFC_aveToMax;
	km.setSeedInitializer(&multiplyLogFFC_aveToMax);

	const vector<int>* x = km.Apply();
	
	SSEDM_Measure ssedmComputer;
	ssedmComputer.setClusterIds(km.getClusterIds());
	ssedmComputer.setInstance(inst);
	double ssedm = ssedmComputer.measure();
	cout << "ssedm: " << ssedm << endl;

	_getch();
}