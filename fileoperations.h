#ifndef _FILEOPERATIONS_H
#define _FILEOPERATIONS_H

#include <iostream> 
#include <fstream>
#include <string>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/GraphAttributes.h>
#include <Eigen/Dense>

void saveProtovis( const char *in_strFilename, 
			  const char *in_strGraphName, 
			  Graph in_graphInput){

	ofstream tmp_oFile (in_strFilename);
	NodeArray<int> tmp_narrayGraph(in_graphInput);

	node n;
	edge e;

	if (!tmp_oFile.is_open()){
		cout << "File can not be opened" << endl;
	}
	else{
		tmp_oFile << "var graph" << in_strGraphName << " = {" << endl;

		// Begin Nodes of a graph
		tmp_oFile << "\tnodes:[";
		
		int tmp_iCount = 0;
		forall_nodes(n,in_graphInput)
		{
			if (tmp_iCount != in_graphInput.numberOfNodes() - 1){
				tmp_oFile << "\n\t\t{nodeName:\"node" << tmp_iCount << "\", group:1" << "},";
			}
			else{
				tmp_oFile << "\n\t\t{nodeName:\"node" << tmp_iCount << "\", group:1" << "}" << endl;
			}
			tmp_narrayGraph[n] = tmp_iCount;
			tmp_iCount++;
		}
	
		tmp_oFile << "\t]," << endl;
		// End Nodes of a graph
		// Begin Edges of a graph
		tmp_oFile << "\tlinks:[";

		tmp_iCount = 0;
		forall_edges(e,in_graphInput)
		{
			if (tmp_iCount != in_graphInput.numberOfEdges() - 1){
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:1" << "},";
			}
			else{
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:1" << "}" << endl;
			}
			tmp_iCount++;
		}

		tmp_oFile << "\t]\n}" << endl;
		// End Edges of a graph
		tmp_oFile.close();
	}
}


void saveProtovisWithOrderUnweighted( const char *in_strFilename, 
			  const char *in_strGraphName, 
			  Graph in_graphInput, vector<int> in_vectorNodesIndex){

	ofstream tmp_oFile (in_strFilename);
	NodeArray<int> tmp_narrayGraph(in_graphInput);
	EdgeArray<int> tmp_earrayGraph(in_graphInput);

	node n;
	edge e;

	if (!tmp_oFile.is_open()){
		cout << "File can not be opened" << endl;
	}
	else{
		tmp_oFile << "var graph" << in_strGraphName << " = {" << endl;

		// Begin Nodes of a graph
		tmp_oFile << "\tnodes:[";
		
		int tmp_iCount = 0;

		vector<node> m_vectorNodes(in_graphInput.numberOfNodes());

		tmp_iCount = 0;
		forall_nodes(n,in_graphInput)
		{
			m_vectorNodes[tmp_iCount] = n;
			tmp_iCount++;
		}

		list<node> m_listNodes;
		for(tmp_iCount = 0; tmp_iCount < in_vectorNodesIndex.size(); ++tmp_iCount)
		{
			m_listNodes.push_back(m_vectorNodes[in_vectorNodesIndex[tmp_iCount]]);
		}

		list<node>::iterator iterate = m_listNodes.begin();

		tmp_iCount = 0;
		for( ; iterate != m_listNodes.end(); ++iterate)
		{
			if (tmp_iCount != in_graphInput.numberOfNodes() - 1){
				tmp_oFile << "\n\t\t{nodeName:\"node" << tmp_iCount << "\", group:1" << "},";
			}
			else{
				tmp_oFile << "\n\t\t{nodeName:\"node" << tmp_iCount << "\", group:1" << "}" << endl;
			}
			tmp_narrayGraph[*iterate] = tmp_iCount;
			tmp_iCount++;
		}
	
		tmp_oFile << "\t]," << endl;
		// End Nodes of a graph
		// Begin Edges of a graph
		tmp_oFile << "\tlinks:[";

		tmp_iCount = 0;
		forall_edges(e,in_graphInput)
		{
			if (tmp_iCount != in_graphInput.numberOfEdges() - 1){
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:1" << "},";
			}
			else{
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:1" << "}" << endl;
			}
			tmp_iCount++;
		}

		tmp_oFile << "\t]\n}" << endl;
		// End Edges of a graph
		tmp_oFile.close();
	}
}

void saveProtovisWithOrderWeighted( const char *in_strFilename, 
			  const char *in_strGraphName, 
			  Graph &in_graphInput, vector<int> in_vectorNodesIndex, EdgeArray<int> in_arrayAdjWeights){
				  
	ofstream tmp_oFile (in_strFilename);
	NodeArray<int> tmp_narrayGraph(in_graphInput);
	EdgeArray<int> tmp_earrayGraph(in_graphInput);

	node n;
	edge e;

	if (!tmp_oFile.is_open()){
		cout << "File can not be opened" << endl;
	}
	else{
		tmp_oFile << "var graph" << in_strGraphName << " = {" << endl;

		// Begin Nodes of a graph
		tmp_oFile << "\tnodes:[";
		
		int tmp_iCount = 0;

		vector<node> m_vectorNodes(in_graphInput.numberOfNodes());

		tmp_iCount = 0;
		forall_nodes(n,in_graphInput)
		{
			m_vectorNodes[tmp_iCount] = n;
			tmp_iCount++;
		}

		list<node> m_listNodes;
		for(tmp_iCount = 0; tmp_iCount < in_vectorNodesIndex.size(); ++tmp_iCount)
		{
			m_listNodes.push_back(m_vectorNodes[in_vectorNodesIndex[tmp_iCount]]);
		}

		list<node>::iterator iterate = m_listNodes.begin();

		tmp_iCount = 0;
		for( ; iterate != m_listNodes.end(); ++iterate)
		{
			if (tmp_iCount != in_graphInput.numberOfNodes() - 1){
				tmp_oFile << "\n\t\t{nodeName:\"Alignment " << in_vectorNodesIndex[tmp_iCount] << "\", group:1" << "},";
			}
			else{
				tmp_oFile << "\n\t\t{nodeName:\"Alignment" << in_vectorNodesIndex[tmp_iCount] << "\", group:1" << "}" << endl;
			}
			tmp_narrayGraph[*iterate] = tmp_iCount;
			tmp_iCount++;
		}
	
		tmp_oFile << "\t]," << endl;
		// End Nodes of a graph
		// Begin Edges of a graph
		tmp_oFile << "\tlinks:[";

		tmp_iCount = 0;
		forall_edges(e,in_graphInput)
		{
			if (tmp_iCount != in_graphInput.numberOfEdges() - 1){
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:" << in_arrayAdjWeights[e] << "},";
			}
			else{
				tmp_oFile << "\n\t\t{source:" << tmp_narrayGraph[e->source()] << " , target:" << tmp_narrayGraph[e->target()] << " , value:" << in_arrayAdjWeights[e] << "}" << endl;
			}
			tmp_iCount++;
		}

		tmp_oFile << "\t]\n}" << endl;
		// End Edges of a graph
		tmp_oFile.close();
	}
}


#endif