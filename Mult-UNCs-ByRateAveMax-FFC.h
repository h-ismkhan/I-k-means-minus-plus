#ifndef _UsefulNearestCenters_LogMultUNCDis_by_RateUNCDisAveToMax_FFC_seedInit_H
#define _UsefulNearestCenters_LogMultUNCDis_by_RateUNCDisAveToMax_FFC_seedInit_H


#include"ISeed-Initializer.h"
#include "IInstance-Euclidean.h"
#include "Time_Random.h"

#include<string>

//[log(dis-1)+...log(dis-1)]*(AveDis/MaxDis), FFC:First Fartesht from Center center of data set
class LogMultUNCDis_by_RateUNCDisAveToMax_FFC_seedInit : public ISeed_Initializer
{
private:
	int thek;
	IInstance_Euclidean* inst;
	int L = INT_MAX;
	vector<int> solution;
private:

public:
	virtual void setMaxProbedPoints(int m){ this->L = m; }
	virtual void setTheK(int theK) { this->thek = theK; /*this->L = theK * 30;*/ }
	virtual void setInstance(IInstance_Euclidean *inst){ this->inst = inst; }
	virtual const int* nextKPoints()
	{
		int dim = this->inst->Dimension();

		vector<int> nodes(dim);
		MakeVectorRandom(&nodes);

		vector<double> initDisFromOtherC(thek, DBL_MAX);
		vector<vector<double>> centers_dis(thek, initDisFromOtherC);

		int nodescount = min(this->L, dim);

		int theFirstCenterIndex = 0;
		for (int i = 1; i < nodescount; i++)
		{
			int node = nodes[i];
			if (inst->get(node, 0) > inst->get(nodes[theFirstCenterIndex], 0))
				theFirstCenterIndex = i;
		}

		solution.assign(thek, -1);
		vector<bool> isin(dim, false);
		solution[0] = nodes[theFirstCenterIndex];
		isin[solution[0]] = true;
		centers_dis[0][0] = 0;

		vector<int> emptyInt;
		vector<vector<int>> usefulNC(dim, emptyInt);
		vector<vector<int>> NOTusefulNC(dim, emptyInt);

		vector<double> emptyDouble;
		vector<vector<double>> usefulNC_dis(dim, emptyDouble);
		vector<vector<double>> NOTusefulNC_dis(dim, emptyDouble);

		for (int i = 1; i < thek; i++)
		{
			int nextnode = -1;
			double nextweight = 0;
			//MakeVectorRandom(&nodes);
			for (int z = 0; z < nodescount; z++)
			{
				int node = nodes[z];
				if (isin[node]) continue;

				double node_nearUsefulC = 0;

				if (usefulNC[node].size() == 0)
				{
					usefulNC[node].push_back(i - 1); // NOT push_back[solution[i - 1]]
					usefulNC_dis[node].push_back(this->inst->Dis(solution[i - 1], node));//NOT Dis(i - 1, node)

					node_nearUsefulC = log(usefulNC_dis[node][0]);
				}
				else
				{
					// i - 1 th center is recnt center added in previous iteration (i - 1) 
					double disBetween_node_Cneter_i_1 = this->inst->Dis(solution[i - 1], node);//solution[i] is not known, yet.

					//we should consider whether the i - 1 th center is useful or not.
					bool the_i_1_is_useful = true;
					for (int useful_index = 0; useful_index < usefulNC[node].size()/*some usefuls may be removed*/; useful_index++)
					{
						double disBetween_PrevUseful_Cneter_i_1 = centers_dis[i - 1][usefulNC[node][useful_index]];
						if (disBetween_PrevUseful_Cneter_i_1 == DBL_MAX)
							centers_dis[i - 1][usefulNC[node][useful_index]] = disBetween_PrevUseful_Cneter_i_1
							= inst->Dis(solution[i - 1], solution[usefulNC[node][useful_index]]);
						double disBetween_node_PrevUseful = usefulNC_dis[node][useful_index];

						//if it is NOT useful, then we should mark it.
						//node==PrevUseful
						//            |
						//              |
						//              center_i-1
						if (disBetween_node_PrevUseful < disBetween_node_Cneter_i_1
							&& disBetween_PrevUseful_Cneter_i_1 < disBetween_node_Cneter_i_1)
						{
							the_i_1_is_useful = false;
							break;
						}
					}
					//the i-1 th center is useful, then by adding it to useful list, some previous usefules become NOT useful
					if (the_i_1_is_useful)
					{
						for (int useful_index = 0; useful_index < usefulNC[node].size()/*some usefuls may be removed*/; useful_index++)
						{
							double disBetween_PrevUseful_Cneter_i_1 = centers_dis[i - 1][usefulNC[node][useful_index]];
							if (disBetween_PrevUseful_Cneter_i_1 == DBL_MAX)
								centers_dis[i - 1][usefulNC[node][useful_index]] = disBetween_PrevUseful_Cneter_i_1
								= inst->Dis(solution[i - 1], solution[usefulNC[node][useful_index]]);
							double disBetween_node_PrevUseful = usefulNC_dis[node][useful_index];

							//node==center_i-1
							//          |
							//            |
							//             PrevUseful
							if (disBetween_node_Cneter_i_1 < disBetween_node_PrevUseful
								&& disBetween_PrevUseful_Cneter_i_1 < disBetween_node_PrevUseful)
							{
								//NOTE:
								//if i-1 th center makes some other un-useful, they should be removed from useful list of node.

								//In addtion, the removed from useful, should be added to NOT useful.
								NOTusefulNC[node].push_back(usefulNC[node][useful_index]);
								NOTusefulNC_dis[node].push_back(usefulNC_dis[node][useful_index]);


								usefulNC[node][useful_index] = usefulNC[node][usefulNC[node].size() - 1];
								usefulNC_dis[node][useful_index] = usefulNC_dis[node][usefulNC[node].size() - 1];

								usefulNC[node].pop_back();
								usefulNC_dis[node].pop_back();
								useful_index--;

							}
						}
					}
					if (the_i_1_is_useful)
					{
						usefulNC[node].push_back(i - 1); // NOT push_back[solution[i - 1]]
						usefulNC_dis[node].push_back(disBetween_node_Cneter_i_1);//NOT Dis(i - 1, node)
					}

					double usefulSizeOf_node = usefulNC[node].size();

					double max_dis = 0, min_dis = DBL_MAX;

					node_nearUsefulC = 0;

					double sum = 0;

					for (int useful_index = 0; useful_index < usefulSizeOf_node; useful_index++)
					{
						node_nearUsefulC = node_nearUsefulC + log(usefulNC_dis[node][useful_index]);
						sum = sum + usefulNC_dis[node][useful_index];
						if (usefulNC_dis[node][useful_index] > max_dis)
							max_dis = usefulNC_dis[node][useful_index];
						if (usefulNC_dis[node][useful_index] < min_dis)
							min_dis = usefulNC_dis[node][useful_index];
					}
					double ave = sum / usefulSizeOf_node;
					node_nearUsefulC = node_nearUsefulC * (ave / max_dis);

				} // end of **else** of "if (usefulNC[node].size() == 0)"


				if (node_nearUsefulC > nextweight)
				{
					nextnode = node;
					nextweight = node_nearUsefulC;
				}
			}
			solution[i] = nextnode;
			isin[nextnode] = true;
			centers_dis[i][i] = 0;

			//to update dstance of thios from other centers, but useful centers
			//?if x is usefu y => y is useful x : centers_dis[x][y]=centers_dis[y][x]=newval?
			int node_usefulSize = usefulNC[nextnode].size();
			for (int k = 0; k < node_usefulSize; k++)
				centers_dis[i][usefulNC[nextnode][k]] = centers_dis[usefulNC[nextnode][k]][i] =
				usefulNC_dis[nextnode][k];
		}
		L = INT_MAX;
		return solution.data();
	}
	virtual string name()
	{
		string s = "_withINT_MAX";
		if (this->L < INT_MAX)
			s = std::to_string(L / thek);

		return "LogMulti-UNC-by-Ave-toMax-FFC-" + s;
	}
};

#endif