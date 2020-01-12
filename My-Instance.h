#ifndef _My_Instance_H
#define _My_Instance_H

#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;
#include "IInstance-Euclidean-Ordinary.h"

class My_Instance : public IInstance_Euclidean_Ordinary
{
public:
	My_Instance(){}
	My_Instance(string path) 
	{
		ReadfromFile(path); 
		//UpdateStatistics(); 
	}
	virtual void ReadfromFile(string path)
	{
		ifstream in(path);
		in >> dim;
		in >> vecDim;
				
		vectors.assign(dim, vector<double>(vecDim));
		
		for (int i = 0; i < dim; i++)
		{
			for(int j = 0; j < vecDim; j++)
				in >> vectors[i][j];
		}
		if (!in.eof())
		{
			int minK = 1, maxK = 0;
			this->actualK = 0;
			actualClusIds.assign(dim, 0);
			for (int j = 0; j < dim; j++)
			{
				if (in.eof())
				{
					this->actualK = -1;
					maxK = -1;
					minK = 0;
					break;
				}
				try
				{					
					in >> actualClusIds[j];
					if (actualClusIds[j] > maxK)
						maxK = actualClusIds[j];
					if (actualClusIds[j] < minK)
						minK = actualClusIds[j];
				}
				catch (exception e)
				{
					maxK = -1;
					minK = 0;
					break;
				}
			}
			if (minK == 0)
				actualK = maxK + 1;
			else
				actualK = maxK;
		}
		in.close();
		//UpdateStatistics();
	}
};

#endif // !_My_Instance