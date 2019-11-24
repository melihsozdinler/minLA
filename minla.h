#ifndef _MINLA_H
#define _MINLA_H

//#include <ogdf/layered/SugiyamaLayout.h>
//#include <ogdf/layered/OptimalRanking.h>
//#include <ogdf/layered/MedianHeuristic.h>
//#include <ogdf/layered/FastHierarchyLayout.h>
//#include <ogdf/planarity/PlanarizationLayout.h>
//#include <ogdf/planarity/VariableEmbeddingInserter.h>
//#include <ogdf/planarity/FastPlanarSubgraph.h>
//#include <ogdf/orthogonal/OrthoLayout.h>
//#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
//#include <ogdf/cluster/ClusterPlanarizationLayout.h>
//#include <ogdf/cluster/ClusterOrthoLayout.h>
//#include <ogdf/cluster/ClusterOrthoShaper.h>
//#include <ogdf/energybased/GEMLayout.h>
//#include <ogdf/energybased/FMMMLayout.h>
//#include <ogdf/basic/geometry.h>
//#include <ogdf/graphalg/ShortestPathWithBFM.h>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/basic/List.h>

#include "eigenobject.h"

#include <stdio.h>
#include <iostream> 
#include <fstream> 
#include <vector>
#include <list>
#include <map>

using namespace ogdf;

struct Entities{
	char entityName[256];
};

typedef struct Entities ENT;

class minLA: public eigenObject{
	public:

		Graph m_graphMinLA;
		GraphAttributes m_oGraphAttributes;
		NodeArray< int > m_arrayComponent;
		NodeArray< int > m_arrayNodeIndex;

		bool m_bComponentSet;
		bool m_bNodeIndexSet;

		vector<list<node> > m_vecListOfNodes;

		int nodeNumber;
		int edgeNumber;

		int m_iNumberOfComponents;

		minLA(const char * in_strFilename)
		{
			m_bComponentSet = false;
			m_bNodeIndexSet = false;
			GraphAttributes tmp_oGraphAttributes(m_graphMinLA,
			GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
			GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
			GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
			GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
			
			tmp_oGraphAttributes.readGML(m_graphMinLA,in_strFilename);
			m_oGraphAttributes = tmp_oGraphAttributes;
		}

		minLA(Graph in_graphMinLA)
		{
			m_bComponentSet = false;
			m_bNodeIndexSet = false;

			GraphAttributes tmp_oGraphAttributes(m_graphMinLA,
			GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
			GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
			GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
			GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
			
			m_graphMinLA = in_graphMinLA;
		}

		void createSimpleGraph(const int in_iNumberOfEdges, const int in_iNumberOfNodes);
		void makeAcyclic();
		void makeAcyclic(Graph *in_graphInput);
		void setComponentNodeArray(NodeArray< int > &in_nodeArray);
		void setComponentNodeArray(NodeArray< int > &in_nodeArray, Graph *in_graphInput);
		void setNodeIndexNodeArray(NodeArray< int > &in_nodeArray);
		void setNodeIndexNodeArray(NodeArray< int > &in_nodeArray, Graph *in_graphInput);
		void saveToFile(const char *in_strFileName);
		void saveToFileGraph(const char *in_strFileName);
		void initializeNodesWithIndexes();
		void initializeNodesWithIndexes(Graph *in_graphInput);
		void findConnectedComponents();
		void findConnectedComponents(Graph *in_graphInput);
		void initializeSolverMatrix();
		void initializeSolverLapMatrices(Graph *in_graphInput);
		void initializeSolverAdjMatrices();
		void initializeSolverLapMatrices();
		void separateComponents();
		void separateComponents(Graph *in_graphInput);

		friend class eigenObject;
};

void minLA::createSimpleGraph(const int in_iNumberOfEdges, const int in_iNumberOfNodes)
{
	randomSimpleGraph(m_graphMinLA, in_iNumberOfNodes, in_iNumberOfEdges);
}

void minLA::makeAcyclic()
{
	DfsAcyclicSubgraph m_oDAS;
	m_oDAS.callAndReverse(m_graphMinLA);
}

void minLA::makeAcyclic(Graph *in_graphInput)
{
	DfsAcyclicSubgraph m_oDAS;
	m_oDAS.callAndReverse(*in_graphInput);
}

void minLA::setComponentNodeArray(NodeArray< int > &in_nodeArray)
{
	m_bComponentSet = true;
	in_nodeArray.init(m_graphMinLA);
}

void minLA::setComponentNodeArray(NodeArray< int > &in_nodeArray, Graph *in_graphInput)
{
	m_bComponentSet = true;
	in_nodeArray.init(*in_graphInput);
}

void minLA::setNodeIndexNodeArray(NodeArray< int > &in_nodeArray)
{
	m_bNodeIndexSet = true;
	in_nodeArray.init(m_graphMinLA);
}

void minLA::setNodeIndexNodeArray(NodeArray< int > &in_nodeArray, Graph *in_graphInput)
{
	m_bNodeIndexSet = true;
	in_nodeArray.init(*in_graphInput);
}

void minLA::saveToFile(const char *in_strFileName)
{
	m_oGraphAttributes.writeGML(in_strFileName);
}

void minLA::saveToFileGraph(const char *in_strFileName)
{
	m_graphMinLA.writeGML(in_strFileName);
}


/*
	Fills m_arrayComponent with component ids and node array can be traced for components
*/
void minLA::findConnectedComponents()
{
	if (true == m_bComponentSet)
		connectedComponents( m_graphMinLA, m_arrayComponent);
	else
		cout << "Node array of components is not set" << endl;
}

void minLA::findConnectedComponents(Graph *in_graphInput)
{
	m_arrayComponent.init(*in_graphInput);
	if (true == m_bComponentSet)
		connectedComponents( *in_graphInput, m_arrayComponent);
	else
		cout << "Node array of components is not set" << endl;
}

void minLA::initializeNodesWithIndexes()
{
	node n;
	int i_count = 0;

	forall_nodes(n, m_graphMinLA)
	{
		m_arrayNodeIndex[ n ] = i_count++;
	}
}

void minLA::initializeNodesWithIndexes(Graph *in_graphInput)
{
	node n;
	int i_count = 0;
	m_arrayNodeIndex.init(*in_graphInput);

	forall_nodes(n, *in_graphInput)
	{
		m_arrayNodeIndex[ n ] = i_count++;
	}
}

void minLA::initializeSolverMatrix()
{
	int i_count = 0, j_count = 0;
	node n,m;

	initializeMatrix(m_graphMinLA.numberOfNodes());
	forall_nodes(n, m_graphMinLA)
	{
		j_count = 0;
		forall_nodes(m, m_graphMinLA)
		{
			m_matrixEigen(i_count,j_count) = 0;
			j_count++;
		}
		i_count++;
	}
	
	edge e;
	forall_edges(e, m_graphMinLA)
	{
		m_matrixEigen(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->target()]) = 1;
	}
}

// prereq is separateComponents
void minLA::initializeSolverAdjMatrices()
{
	int tmp_kCount = 0;
	node n,m;
	edge e;
	if (m_iNumberOfComponents > 0)
	{
		vector<int> tmp_vectorIndexes(m_iNumberOfComponents+1);
		m_vectorMatricesEigen.resize(m_vecListOfNodes.size());

		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
		{
			initializeMatrix(m_vecListOfNodes.at(tmp_kCount).size(),tmp_kCount);
			for (int tmp_iCount = 0; tmp_iCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_iCount++)
			{
				for (int tmp_jCount = 0; tmp_jCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_jCount++)
				{
					m_vectorMatricesEigen.at(tmp_kCount)(tmp_iCount,tmp_jCount) = 0;
				}
			}
			tmp_vectorIndexes.at(tmp_kCount) = 0;
		}

		/*for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);*/
		
		//Fill the map node, index
		forall_nodes(n, m_graphMinLA)
		{
			m_arrayNodeIndex[ n ] = tmp_vectorIndexes.at(m_arrayComponent[ n ]);
			cout << m_arrayComponent[ n ] << " - " << m_arrayNodeIndex[ n ] << endl;
			tmp_vectorIndexes.at(m_arrayComponent[ n ])++;
		}

		forall_edges(e, m_graphMinLA)
		{
			if (m_arrayComponent[e->source()] == m_arrayComponent[e->target()])
			{
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->target()]) = 1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->target()],m_arrayNodeIndex[e->source()]) = 1;
			}
		}
		
		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);

		reserveSolverMemory(m_iNumberOfComponents);
	}
}


// prereq is separateComponents
void minLA::initializeSolverLapMatrices()
{
	int tmp_kCount = 0;
	node n,m;
	edge e;
	if (m_iNumberOfComponents > 0)
	{
		vector<int> tmp_vectorIndexes(m_iNumberOfComponents+1);
		m_vectorMatricesEigen.resize(m_vecListOfNodes.size());

		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
		{
			initializeMatrix(m_vecListOfNodes.at(tmp_kCount).size(),tmp_kCount);
			for (int tmp_iCount = 0; tmp_iCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_iCount++)
			{
				for (int tmp_jCount = 0; tmp_jCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_jCount++)
				{
					m_vectorMatricesEigen.at(tmp_kCount)(tmp_iCount,tmp_jCount) = 0;
				}
			}
			tmp_vectorIndexes.at(tmp_kCount) = 0;
		}

		/*for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);*/
		
		//Fill the map node, index
		forall_nodes(n, m_graphMinLA)
		{
			m_arrayNodeIndex[ n ] = tmp_vectorIndexes.at(m_arrayComponent[ n ]);
			cout << m_arrayComponent[ n ] << " - " << m_arrayNodeIndex[ n ] << endl;
			tmp_vectorIndexes.at(m_arrayComponent[ n ])++;
		}

		forall_edges(e, m_graphMinLA)
		{
			if (m_arrayComponent[e->source()] == m_arrayComponent[e->target()])
			{
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->target()]) = -1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->target()],m_arrayNodeIndex[e->source()]) = -1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->source()]) += 1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->target()],m_arrayNodeIndex[e->target()]) += 1;
			}
		}
		
		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);

		reserveSolverMemory(m_iNumberOfComponents);
	}
}

// prereq is separateComponents
void minLA::initializeSolverLapMatrices(Graph *in_graphInput)
{
	int tmp_kCount = 0;
	node n,m;
	edge e;
	if (m_iNumberOfComponents >= 0)
	{
		vector<int> tmp_vectorIndexes(m_iNumberOfComponents+1);
		m_vectorMatricesEigen.resize(m_vecListOfNodes.size());

		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
		{
			initializeMatrix(m_vecListOfNodes.at(tmp_kCount).size(),tmp_kCount);
			for (int tmp_iCount = 0; tmp_iCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_iCount++)
			{
				for (int tmp_jCount = 0; tmp_jCount < m_vecListOfNodes.at(tmp_kCount).size(); tmp_jCount++)
				{
					m_vectorMatricesEigen.at(tmp_kCount)(tmp_iCount,tmp_jCount) = 0;
				}
			}
			tmp_vectorIndexes.at(tmp_kCount) = 0;
		}

		/*for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);*/
		
		//Fill the map node, index
		forall_nodes(n, *in_graphInput)
		{
			m_arrayNodeIndex[ n ] = tmp_vectorIndexes.at(m_arrayComponent[ n ]);
			cout << m_arrayComponent[ n ] << " - " << m_arrayNodeIndex[ n ] << endl;
			tmp_vectorIndexes.at(m_arrayComponent[ n ])++;
		}

		forall_edges(e, *in_graphInput)
		{
			if (m_arrayComponent[e->source()] == m_arrayComponent[e->target()])
			{
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->target()]) = -1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->target()],m_arrayNodeIndex[e->source()]) = -1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->source()],m_arrayNodeIndex[e->source()]) += 1;
				m_vectorMatricesEigen.at(m_arrayComponent[e->source()])(m_arrayNodeIndex[e->target()],m_arrayNodeIndex[e->target()]) += 1;
			}
		}
		
		for (tmp_kCount = 0; tmp_kCount < m_vecListOfNodes.size(); tmp_kCount++)
			printMatrix(tmp_kCount);

		reserveSolverMemory(m_iNumberOfComponents);
	}
}

void minLA::separateComponents()
{
	m_iNumberOfComponents = 0;
	node n;
	int i_count = 0;

	forall_nodes(n, m_graphMinLA)
	{
		if( m_arrayComponent[ n ] > m_iNumberOfComponents )
			m_iNumberOfComponents = m_arrayComponent[ n ];
	}

	m_vecListOfNodes.resize(m_iNumberOfComponents+1);

	forall_nodes(n, m_graphMinLA)
	{
		m_vecListOfNodes.at( m_arrayComponent[ n ] ).push_back(n);
	}
}

void minLA::separateComponents(Graph *in_graphInput)
{
	m_iNumberOfComponents = 0;
	node n;
	int i_count = 0;

	forall_nodes(n, *in_graphInput)
	{
		if( m_arrayComponent[ n ] > m_iNumberOfComponents )
			m_iNumberOfComponents = m_arrayComponent[ n ];
	}

	m_vecListOfNodes.resize(m_iNumberOfComponents+1);

	forall_nodes(n, *in_graphInput)
	{
		m_vecListOfNodes.at( m_arrayComponent[ n ] ).push_back(n);
	}
}

#endif