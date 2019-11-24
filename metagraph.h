#ifndef _METAGRAPH_H
#define _METAGRAPH_H

#include <string>
#include "readinput.h"

class metaGraph
{
public:

	Graph m_graphMeta;
	int m_iAlignmentSize;

	metaGraph(const int in_iAlignmentSize)
	{
		m_iAlignmentSize = in_iAlignmentSize;
	}

	
};



#endif