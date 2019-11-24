#ifndef _READINPUT_H
#define _READINPUT_H

#include <string>
#include <Eigen/Dense>

#include "minla.h"

using namespace std;

class readInput{
	public:	
	Graph m_graphInput;
	Graph m_graphMeta;
	MatrixXi m_matrixEdgeWeights;
	EdgeArray<int> m_arrayMetaEdgeWeights;

	int m_iAlignmentCount;
	vector<list<node> > m_vectorAlignment;
	NodeArray<list<int> > m_arrayAlignmentMap;
	
	readInput(Graph &in_graphInp)
	{
		m_graphInput = in_graphInp;
	}

	void readAlignmentFileOption1(const char *in_strFilename);
	int findAlignmentCount(const char *in_strFilename);
	char * safeCopyCharArray(const char *in_strToCopy);
	node findNode(const int in_iIndex);
	void printAlignments();
	void convertVectorToNodeArray();
	void formMetaGraph();
};

void readInput::convertVectorToNodeArray()
{
	m_arrayAlignmentMap.init(m_graphInput);

	for ( int tmp_iCount = 0; tmp_iCount < m_vectorAlignment.size(); ++tmp_iCount)
	{
		list<node>::iterator iterate = m_vectorAlignment[tmp_iCount].begin();
		for (; iterate != m_vectorAlignment[tmp_iCount].end(); iterate++ )
		{
			m_arrayAlignmentMap[(*iterate)].push_back(tmp_iCount);
		}
	}

	node n;

	forall_nodes(n, m_graphInput)
	{
		cout << m_arrayAlignmentMap[n].size() << endl;
	}
}

void readInput::printAlignments()
{
	for (int tmp_iCount = 0; tmp_iCount < m_vectorAlignment.size(); ++tmp_iCount)
	{
		cout << "Alignment -" << tmp_iCount << ":";
		list<node>::iterator iterate = m_vectorAlignment[tmp_iCount].begin();
		for ( ; iterate != m_vectorAlignment[tmp_iCount].end(); ++iterate )
		{
			cout << *iterate << "-";
		}
		cout << endl;
	}
}

char * readInput::safeCopyCharArray(const char *in_strToCopy)
{
	char *tmp_strSave;
	tmp_strSave = (char*)malloc(sizeof(strlen(in_strToCopy)));
	strncpy(tmp_strSave, in_strToCopy, strlen(in_strToCopy));
	tmp_strSave[ strlen(in_strToCopy) ] = '\0';
	return tmp_strSave;
}

node readInput::findNode(const int in_iIndex)
{
	node n;
	int tmp_iCount = 0;
	forall_nodes(n,m_graphInput)
	{
		if (in_iIndex == tmp_iCount)
		{
			return n;
		}
		tmp_iCount++;
	}
}

int readInput::findAlignmentCount(const char *in_strFilename)
{
	string tmp_strLine;
	ifstream tmp_oFile (in_strFilename);
	int tmp_iCount = 0;

	if (tmp_oFile.is_open())
	{
		while ( tmp_oFile.good() )
		{
			getline (tmp_oFile,tmp_strLine);			
			tmp_iCount++;
		}
		tmp_oFile.close();
	}
	else
	{
		cout << "Unable to open file";
		return -1;
	}

	return tmp_iCount;
}

void readInput::readAlignmentFileOption1(const char *in_strFilename)
{
	node n,tmp_nodeSource, tmp_nodeTarget;
	int tmp_iCount = 0;
	string tmp_strLine;
	ifstream tmp_oFile (in_strFilename);
	m_iAlignmentCount = findAlignmentCount(in_strFilename);
	m_vectorAlignment.resize(m_iAlignmentCount);

	if (tmp_oFile.is_open())
	{
		while ( tmp_oFile.good() )
		{
			getline (tmp_oFile,tmp_strLine);
			cout << tmp_strLine << endl;

			int tmp_intPosition = tmp_strLine.find_first_of(" \t");
			string tmp_strPartString = tmp_strLine.substr(0,tmp_intPosition);
			tmp_nodeSource = findNode(atoi(tmp_strPartString.c_str()));

			tmp_strLine = tmp_strLine.substr(tmp_intPosition+1);
			tmp_intPosition = tmp_strLine.find_first_of(" \t\n");
			tmp_strPartString = tmp_strLine.substr(0,tmp_intPosition);
			tmp_nodeTarget = findNode(atoi(tmp_strPartString.c_str()));
			tmp_strLine = tmp_strLine.substr(tmp_intPosition+1);

			list<node> tmp_listNodes;
			tmp_listNodes.push_back(tmp_nodeTarget);
			while ( tmp_strLine.empty() != true)
			{
				tmp_intPosition = tmp_strLine.find_first_of(" \t\n");
				if (-1 == tmp_intPosition)
				{
					tmp_strPartString =  tmp_strLine;
					tmp_intPosition = 0;
				}
				else
				{
					tmp_strPartString = tmp_strLine.substr(0,tmp_intPosition);
				}
				tmp_nodeTarget = findNode(atoi(tmp_strPartString.c_str()));
				tmp_strLine = tmp_strLine.substr(tmp_intPosition+1);
				tmp_listNodes.push_back(tmp_nodeTarget);
			}
			m_vectorAlignment[tmp_iCount] = tmp_listNodes;
			tmp_iCount++;
		}
		tmp_oFile.close();
	}
	else
	{
		cout << "Unable to open file";
	}
}

void readInput::formMetaGraph()
{
	node source, target, n, m;
	list<int>::iterator iterateSource, iterateTarget;
	list<node> tmp_listSourceAligns, tmp_listTargetAligns;
	list<node>::iterator iterateSourceNodes, iterateTargetNodes;
	

	m_iAlignmentCount = m_vectorAlignment.size();

	vector<node> tmp_vectorNodes(m_iAlignmentCount);
	for (int tmp_iCount = 0; tmp_iCount < m_iAlignmentCount; ++tmp_iCount)
	{
		node tmp_nodeCreated = m_graphMeta.newNode();
		tmp_vectorNodes[tmp_iCount] = tmp_nodeCreated;
	}

	MatrixXi tmp_matrix(m_graphMeta.numberOfNodes(),m_graphMeta.numberOfNodes());

	int j_count = 0, i_count = 0;
	forall_nodes(n, m_graphMeta)
	{
		j_count = 0;
		forall_nodes(m, m_graphMeta)
		{
			tmp_matrix(i_count,j_count) = 0;
			j_count++;
		}
		i_count++;
	}

	NodeArray<int> m_narrayIndex(m_graphMeta, 0);

	int tmp_iCount = 0;
	forall_nodes(n, m_graphMeta)
	{
		m_narrayIndex[n] = tmp_iCount;
		tmp_iCount++;
	}

	edge tmp_edgeE;
	forall_edges(tmp_edgeE,m_graphInput)
	{
		source = tmp_edgeE->source();
	    target = tmp_edgeE->target();
		
		if ( m_arrayAlignmentMap[source].size() > 0 &&  m_arrayAlignmentMap[target].size() > 0)
		{

			iterateSource = m_arrayAlignmentMap[source].begin();
			iterateTarget = m_arrayAlignmentMap[target].begin();

			tmp_listSourceAligns.clear();
			tmp_listTargetAligns.clear();

			for ( ; iterateSource != m_arrayAlignmentMap[source].end(); ++iterateSource )
			{
				tmp_listSourceAligns.push_back(tmp_vectorNodes[*iterateSource]);
			}

			for ( ; iterateTarget != m_arrayAlignmentMap[target].end(); ++iterateTarget )
			{
				tmp_listTargetAligns.push_back(tmp_vectorNodes[*iterateTarget]);
			}

			if ( tmp_listTargetAligns.empty() != true && tmp_listSourceAligns.empty() != true )
			{
				iterateSourceNodes = tmp_listSourceAligns.begin();
				iterateTargetNodes = tmp_listTargetAligns.begin();
				for ( ; iterateSourceNodes != tmp_listSourceAligns.end(); ++iterateSourceNodes )
				{
					for ( ; iterateTargetNodes != tmp_listTargetAligns.end(); ++iterateTargetNodes )
					{
						if ( *iterateSourceNodes != *iterateTargetNodes )
						{
							tmp_matrix(m_narrayIndex[*iterateSourceNodes],m_narrayIndex[*iterateTargetNodes])++; 
							tmp_matrix(m_narrayIndex[*iterateTargetNodes],m_narrayIndex[*iterateSourceNodes])++;
							m_graphMeta.newEdge(*iterateSourceNodes, *iterateTargetNodes);
						}
					}
				}
			}
		}
	}

	m_arrayMetaEdgeWeights.init(m_graphMeta);

	edge e;
	forall_edges(e,m_graphMeta)
	{
		m_arrayMetaEdgeWeights[e] = tmp_matrix(m_narrayIndex[e->source()],m_narrayIndex[e->target()]);
	}

	// TODO use this matrix during edge creation
	cout << tmp_matrix;
	m_graphMeta.writeGML("metagraph.txt");
}

#endif