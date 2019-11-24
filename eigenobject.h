#ifndef _EIGENOBJECT_H
#define _EIGENOBJECT_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <list>
#include <map> 

using namespace std;
using namespace Eigen;

class eigenObject
{
	public:
		SelfAdjointEigenSolver<MatrixXf> m_oEigenSolver;
		vector<SelfAdjointEigenSolver<MatrixXf> > m_vectorEigenSolvers;	
		void printEigenValues();
		void printMatrix();
		void printMatrix(const int in_iIndex);
		void printEigenVectors();
		void initializeSolver();
		void initializeSolver(const int index);
		list<vector<int> > solveAll();

	protected:
		MatrixXf m_matrixEigen;
		vector<MatrixXf> m_vectorMatricesEigen;
		void initializeMatrix(const int in_iNumberOfNodes);
		void initializeMatrix(const int in_iNumberOfNodes, const int in_iIndex);
		void reserveSolverMemory(const int in_iNumberOfSolvers);
};


void eigenObject::initializeSolver()
{
	SelfAdjointEigenSolver<MatrixXf> tmp_eigensolver(m_matrixEigen);
	m_oEigenSolver = tmp_eigensolver;
}

void eigenObject::initializeSolver(const int in_iIndex)
{
	SelfAdjointEigenSolver<MatrixXf> tmp_eigensolver(m_vectorMatricesEigen.at(in_iIndex));
	m_vectorEigenSolvers.at(in_iIndex) = tmp_eigensolver;
}

void eigenObject::reserveSolverMemory(const int in_iNumberOfSolvers)
{
	m_vectorEigenSolvers.resize(in_iNumberOfSolvers+1);
}

void eigenObject::printMatrix()
{
	cout << "Here is the matrix A:\n" 
		 << m_matrixEigen 
		 << endl;
}

void eigenObject::printMatrix(const int in_iIndex)
{
	cout << "Here is the matrix A:\n" 
		 << m_vectorMatricesEigen.at(in_iIndex)
		 << endl;
}

void eigenObject::printEigenValues()
{
	cout << "The eigenvalues of matrix are:\n" 
		 << m_oEigenSolver.eigenvalues() 
		 << endl;
}

void eigenObject::printEigenVectors()
{
	cout << "Here'are eigenvectors of a Matrix corresponding to eigenvalues:\n"
		 << m_oEigenSolver.eigenvectors() 
		 << endl;
}

void eigenObject::initializeMatrix(const int in_iNumberOfNodes)
{
	MatrixXf tmp_matrix(in_iNumberOfNodes,in_iNumberOfNodes);
	m_matrixEigen = tmp_matrix;
}

void eigenObject::initializeMatrix(const int in_iNumberOfNodes, const int in_iIndex)
{
	MatrixXf tmp_matrix(in_iNumberOfNodes,in_iNumberOfNodes);
	m_vectorMatricesEigen.at(in_iIndex) = tmp_matrix;
}

list<vector<int> > eigenObject::solveAll()
{
	list<vector<int> > tmp_listVectorOrders;
	int tmp_iEigenVecCount1 = 0, tmp_intVectorSize = 0;

	for (int tmp_iCount = 0; tmp_iCount < m_vectorEigenSolvers.size(); tmp_iCount++)
	{
		if (m_vectorEigenSolvers.at(tmp_iCount).info() != Success) 
			abort();
		else
		{
			tmp_intVectorSize =  m_vectorMatricesEigen.at(tmp_iCount).rows();

			/*cout << "The eigenvalues of matrix are:\n" 
				 << m_vectorEigenSolvers.at(tmp_iCount).eigenvalues()
				 << endl;*/
			
			// to learn matrix row number: m_vectorMatricesEigen.at(tmp_iCount).rows()
			// to learn eigen vector matrix second row and column m_vectorEigenSolvers.at(tmp_iCount).eigenvectors().col(1).row(1)

			list<double> tmp_listSortedValues;

			for (tmp_iEigenVecCount1 = 0; tmp_iEigenVecCount1 < tmp_intVectorSize ; ++tmp_iEigenVecCount1)
			{
				tmp_listSortedValues.push_back(m_vectorEigenSolvers.at(tmp_iCount).eigenvectors().col(1).row(tmp_iEigenVecCount1).value());
			}
			tmp_listSortedValues.sort();

			vector<int> tmp_vectorSortedOrders(tmp_listSortedValues.size());
			list<double>::iterator iterate = tmp_listSortedValues.begin();
			vector<bool> tmp_boolMarked(tmp_listSortedValues.size(), false);

			for ( tmp_iEigenVecCount1 = 0; iterate != tmp_listSortedValues.end() ; ++tmp_iEigenVecCount1,++iterate)
			{
				for ( int tmp_iEigenVecCount2 = 0; tmp_iEigenVecCount2 < tmp_intVectorSize ; ++tmp_iEigenVecCount2)
				{
					if( false == tmp_boolMarked[tmp_iEigenVecCount2] && *(iterate) == m_vectorEigenSolvers.at(tmp_iCount).eigenvectors().col(1).row(tmp_iEigenVecCount2).value() )
					{
						tmp_vectorSortedOrders[tmp_iEigenVecCount1] = tmp_iEigenVecCount2;
						tmp_boolMarked[tmp_iEigenVecCount2] = true;
						break;
					}
				}
			}
			tmp_listVectorOrders.push_back(tmp_vectorSortedOrders);
		}
	}
	return tmp_listVectorOrders;
}

#endif 