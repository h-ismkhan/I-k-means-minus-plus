#ifndef _ISeed_Initializer_H
#define _ISeed_Initializer_H


#include "IInstance-Euclidean.h"

class ISeed_Initializer
{
public:
	virtual void setTheK(int theK) = 0;
	virtual void setInstance(IInstance_Euclidean *) = 0;
	virtual const int* nextKPoints() = 0;
	virtual string name() = 0;
};

#endif