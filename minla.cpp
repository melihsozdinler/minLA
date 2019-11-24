#include "minla.h"
#include "readInput.h"
#include "metagraph.h"
#include "fileoperations.h"

void solveMetaGraph();
readInput* readAlignmentInput(Graph *in_graphMinLA);
void readNetworkGraphAndSolve(minLA *m_pMinLA);

int main()
{	
	solveMetaGraph();
	return 0;
}


void solveMetaGraph()
{
	int tmp_iCount = 0;
	node n;

	// Use minLA class for remaining operations
	minLA *m_pMinLA = new minLA("graph.gml");	
	// Save graph as protovis arc format
	saveProtovis("saveprot.js", "graphmy", m_pMinLA->m_graphMinLA);
	// Use readInput class for alignment read operations
	readInput *m_pReadInput = readAlignmentInput(&m_pMinLA->m_graphMinLA);
	
	minLA *m_pMinLAM = new minLA(m_pReadInput->m_graphMeta);

	m_pMinLAM->setComponentNodeArray(m_pMinLAM->m_arrayComponent,&m_pReadInput->m_graphMeta);
	m_pMinLAM->setNodeIndexNodeArray(m_pMinLAM->m_arrayNodeIndex,&m_pReadInput->m_graphMeta);
	m_pMinLAM->findConnectedComponents(&m_pReadInput->m_graphMeta);
	m_pMinLAM->makeAcyclic(&m_pReadInput->m_graphMeta);
	m_pMinLAM->initializeNodesWithIndexes(&m_pReadInput->m_graphMeta);	
	m_pMinLAM->separateComponents(&m_pReadInput->m_graphMeta);
	m_pMinLAM->initializeSolverLapMatrices(&m_pReadInput->m_graphMeta);
	for (tmp_iCount = 0; tmp_iCount <= m_pMinLAM->m_iNumberOfComponents; tmp_iCount++)
	{
		m_pMinLAM->initializeSolver(tmp_iCount);
	}
	list<vector<int> > m_listVectorSolutions = m_pMinLAM->solveAll();
	
	m_pMinLAM->saveToFileGraph("metagraph.gml");
	saveProtovisWithOrderWeighted("saveportmeta.js", "graphmy", m_pReadInput->m_graphMeta, m_listVectorSolutions.front(), m_pReadInput->m_arrayMetaEdgeWeights);

	delete m_pMinLAM;
	delete m_pMinLA;
	delete m_pReadInput;
}

readInput* readAlignmentInput(Graph *in_graphMinLA)
{
	readInput *m_pReadInput = new readInput(*in_graphMinLA);
	m_pReadInput->readAlignmentFileOption1("inputalignment.txt");
	m_pReadInput->printAlignments();
	m_pReadInput->convertVectorToNodeArray();
	m_pReadInput->formMetaGraph();
	return m_pReadInput;
}

void readNetworkGraphAndSolve(minLA *m_pMinLA)
{	
	// Use minLA class for remaining operations
	m_pMinLA = new minLA("graph.gml");
	//m_pMinLA->createSimpleGraph(50,50);
	//m_pMinLA->makeAcyclic();
	m_pMinLA->setComponentNodeArray(m_pMinLA->m_arrayComponent);
	m_pMinLA->setNodeIndexNodeArray(m_pMinLA->m_arrayNodeIndex);
	m_pMinLA->findConnectedComponents();
	m_pMinLA->makeAcyclic();
	//m_pMinLAM->saveToFile("test.gml");
	m_pMinLA->initializeNodesWithIndexes();	
	m_pMinLA->separateComponents();
	m_pMinLA->initializeSolverLapMatrices();
	for (int tmp_iCount = 0; tmp_iCount < m_pMinLA->m_iNumberOfComponents; tmp_iCount++)
	{
		m_pMinLA->initializeSolver(tmp_iCount);
	}
	//m_pMinLAM->printMatrix();
	m_pMinLA->solveAll();
	//if (m_pMinLA->m_oEigenSolver.info() != Success) 
	//	abort();	
	//m_pMinLA->printEigenValues();
	//m_pMinLA->printEigenVectors();
	cout << m_pMinLA->m_graphMinLA.numberOfEdges() << " - " <<  m_pMinLA->m_graphMinLA.numberOfNodes() << endl;
	
	delete m_pMinLA;
}