#include "cliqueCreator.h"

#include<iostream>
#include<fstream>


List<edge> simpleMIMMLParserContingency( char filename[64], Graph &G, NodeArray<int> &Gn, NodeArray<int> &Type, NodeArray<Entities> &ENTITIES, List<CONTINGENCY> &addedNodes, List<CONTINGENCY> &contingencyList ){
	FILE *fptr,*fptr2;
	List<edge> added;
	fptr2 = std::fopen( "outputs/reactionTypes.txt", "a" );
	//GRAPH<int,int> G;	/* Introducing new graph */
	node n,n2,walkOnNodes;
	edge e = NULL;
	NodeArray<int> Type2;
	ifstream myReadFile;
	myReadFile.open( filename );

	if( myReadFile.is_open() == NULL ){
		cout << "\n The File named " << filename << " does not exist\n" << endl;
	}
	else{
		cout <<  "/* File opened */\n";
		char parse[128],parse2[128],parse3[128],parse4[128],parse5[128],rname[128],rname2[128];
		myReadFile >>  parse ;
		List<node> NODES;
		List<Entities> ENTY;
		int satir = 1;
		int typeNow = 1;
		while( !myReadFile.eof() ){
			CONTINGENCY ctemp;
			myReadFile >>  parse ;
			if( strcmp( parse, "reactions:" ) == 0 ){
				myReadFile >>  parse ;
				typeNow = 1;
			}
			if( strcmp( parse, "contingencies:" ) == 0 ){
				myReadFile >>  parse ;
				typeNow = 2;
			}
			if( strcmp( parse, "-" ) == 0 ){
				myReadFile >>  rname;
				myReadFile >>  parse2;
				myReadFile >>  parse3;
				myReadFile >>  parse4;
				myReadFile >>  rname2;
				//fscanf( fptr, "%s%s%s%s%s",rname, parse2, parse3, parse4, rname2 );
				for( int k1 = 0; rname[ k1 ] != '\0'; k1++ ){
					if( rname[ k1 ] == ':' )
						rname[ k1 ] = '\0';
				}
				node ent_s1,ent_t1;
				bool found1 = false,found2 = false;
				forall_nodes( n, G ){
					if( ( strcmp( ENTITIES[ n ].entityName, parse2) == 0 || strcmp( ENTITIES[ n ].reactionName, parse2) == 0 ) ){
						ent_s1 = n;
						found1 = true;
					}
					if( ( strcmp( ENTITIES[ n ].entityName, parse4) == 0 || strcmp( ENTITIES[ n ].reactionName, parse4) == 0 ) ){
						ent_t1 = n;
						found2 = true;
					}
				}
				if( found1 != true ){
					ent_s1 = G.newNode();
					Type2.init(G);
					Type2[ ent_s1 ] = 0;
					forall_nodes( walkOnNodes, G){
						if( walkOnNodes != ent_s1 )
							Type2[ walkOnNodes ] = Type[ walkOnNodes ];
					}
					Type = Type2;
					NodeArray<Entities> ENTITIEST = ENTITIES;
					ENTITIES.init( G );
					forall_nodes( walkOnNodes, G ){
						if( walkOnNodes != ent_s1 ){
							ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
						}
					}
					sprintf( ENTITIES[ ent_s1 ].entityName, "%s", parse2 );
				}
				if( found2 != true ){
					ent_t1 = G.newNode();
					Type2.init(G);
					Type2[ ent_t1 ] = 0;
					forall_nodes( walkOnNodes, G){
						if( walkOnNodes != ent_t1 )
							Type2[ walkOnNodes ] = Type[ walkOnNodes ];
					}
					Type = Type2;
					NodeArray<Entities> ENTITIEST = ENTITIES;
					ENTITIES.init( G );
					forall_nodes( walkOnNodes, G ){
						if( walkOnNodes != ent_t1 ){
							ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
						}
					}
					sprintf( ENTITIES[ ent_t1 ].entityName, "%s", parse4 );
				}
				/*if( found1 && found2 ){*/
					NodeArray<int> Type2 = Type;
					n = G.newNode();
					Type.init(G);
					Type[ n ] = 9;
					forall_nodes( walkOnNodes, G){
						if( walkOnNodes != n )
							Type[ walkOnNodes ] = Type2[ walkOnNodes ];
					}
					NodeArray<Entities> ENTITIEST = ENTITIES;
					ENTITIES.init( G );
					forall_nodes( walkOnNodes, G ){
						if( walkOnNodes != n ){
							ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
						}
					}
					sprintf( ENTITIES[ n ].reactionName, "%s", rname );
					sprintf( ENTITIES[ n ].entityName, "%s", rname );
					CONTINGENCY X;
					e = G.newEdge( ent_s1, n );
					ctemp.source1ToComplex = e;
					X.source1ToComplex = e;
					added.pushBack( e );
					e = G.newEdge( n, ent_t1 );
					ctemp.complexToTarget1 = e;
					X.complexToTarget1 = e;
					added.pushBack( e );
					X.c1 = n;					
					addedNodes.pushBack( X );
					ctemp.c1 = n;
					ctemp.s1 = ent_s1;
					ctemp.t1 = ent_t1;
					sprintf( ctemp.contigName, "%s", ENTITIES[ n ].entityName );
					sprintf( ctemp.source1, "%s", ENTITIES[ ent_s1 ].entityName );
					sprintf( ctemp.target1, "%s", ENTITIES[ ent_t1 ].entityName );
					ctemp.typeInt = getContingecyType( rname2 );
					sprintf( ctemp.typeName, "%s", rname2 );
				//}
			}
			contingencyList.pushBack( ctemp );
		}
	}
	myReadFile.close();
	return added;
}

NodeArray<Entities> simpleMIMMLParserReactionAndContingency( char filename[64], Graph &G, NodeArray<int> &Gn, NodeArray<int> &Type, List<node> &reactionAmbigous, List<REACTION> &RList ){
	FILE *fptr,*fptr2;
	fptr2 = std::fopen( "outputs/reactionTypes.txt", "w" );
	//GRAPH<int,int> G;	/* Introducing new graph */
	node n,n2,walkOnNodes;
	NodeArray<Entities> ENTITIES( G );
	NodeArray<Entities> ENTITIEST( G );
	edge e = NULL;
	NodeArray<int> Type2;
ifstream myReadFile;
	myReadFile.open( filename );

	if( myReadFile.is_open() == NULL ){
		cout << "\n The File named " << filename << " does not exist\n" << endl;
	}
	else{
		cout <<  "/* File opened */\n";
		char parse[128],parse2[128],parse3[128],parse4[128],parse5[128],rname[128],rname2[128];
		myReadFile >>  parse ;
		List<node> NODES;
		List<Entities> ENTY;
		int satir = 1;
		int typeNow = 1;
		while( !myReadFile.eof() ){	
			REACTION rtemp;
			myReadFile >>  parse ;
			if( strcmp( parse, "reactions:" ) == 0 ){
			/*  'ncb' : Non-covalent binding (SBO:0000177)
				'cvm' : Covalent modification (SBO:0000210 ! addition of a chemical group)
				'cvb' : Covalent bond (SBO:0000210 ! addition of a chemical group, ?)
				'stc' : Stoichiometric conversion (SBO:0000182 ! conversion)
				'syn' : Production without loss of reactants (SBO:0000205 - biochemical process, ?)
				'tra' : Transcription (SBO:0000183)
				'cle' : Cleavage (SBO:0000178)
				'deg' : Degradation (SBO:0000179)	*/
				myReadFile >>  parse ;
				typeNow = 1;
			}
			if( strcmp( parse, "contingencies:" ) == 0 ){
				myReadFile >>  parse ;
				typeNow = 2;
			}
			if( strcmp( parse, "-" ) == 0 ){
				myReadFile >>  rname ;
				myReadFile >>  parse2 ;
				myReadFile >>  parse3 ;
				//fscanf( fptr, "%s%s%s",rname, parse2, parse3);
				int k2 = 0;
				/*for( int k1 = 0; rname2[ k1 ] != '\0'; k1++ ){
					rname2[ k1 ] = ' ';
				}
				for( int k1 = 0; rname[ k1 ] != '\0'; k1++ ){
					rname2[ k2 ] = rname[ k1 ];
					if( rname[ k1 ] == ':' )
						rname2[ k2 ] = '\0';
					k2++;
				}

				sprintf( rname, "%s", rname2 );*/
				for( int k1 = 0; rname[ k1 ] != '\0'; k1++ ){
					if( rname[ k1 ] == ':' )
						rname[ k1 ] = '\0';
				}
				if( strcmp( parse3, "+" ) == 0 ){
					myReadFile >>  parse4 ;
					myReadFile >>  parse5 ;
					//fscanf( fptr, "%s%s", parse4, parse5 );
					
					/* parse2 is one node, parse3 is adding operator,
					parse4 is the second node of entities and parse 5 is
					assignment opertor to determine reaction type. */ 
					node ent_s1,ent_s2,ent_t1,ent_t2;
					bool found1 = false,found2 = false;
					forall_nodes( n, G ){
						if( ( strcmp( ENTITIES[ n ].entityName, parse2) == 0 || strcmp( ENTITIES[ n ].reactionName, parse2) == 0 ) && strcmp( "P", parse2 ) != 0 ){
							ent_s1 = n;
							found1 = true;
						}
						if( ( strcmp( ENTITIES[ n ].entityName, parse4) == 0 || strcmp( ENTITIES[ n ].reactionName, parse4) == 0 ) &&  strcmp( "P", parse4 ) != 0 ){
							ent_s2 = n;
							found2 = true;
						}
					}
					if( found1 == false ){
						Type = Type2;
						ent_s1 = G.newNode();

						Type2.init(G);
						Type2[ ent_s1 ] = 0;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != ent_s1 )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != ent_s1 ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_s1 ].entityName, "%s", parse2 );
// 						cout << ENTITIES[ ent_s1 ].entityName << "\t" <<  parse2 << "\n";
					}

					if( found2 == false ){
						Type = Type2;
						ent_s2 = G.newNode();
						Type2.init(G);
						Type2[ ent_s2 ] = 0;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != ent_s2 )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						NodeArray<Entities> ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != ent_s2 ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_s2 ].entityName, "%s", parse4 );
					}
 				
					myReadFile >>  parse2 ;
					myReadFile >>  parse3 ;
					//fscanf( fptr, "%s%s", parse2, parse3 );
					found2 = false;
					forall_nodes( n, G ){
						if( ( strcmp( ENTITIES[ n ].entityName, parse2 ) == 0 || strcmp( ENTITIES[ n ].reactionName, parse2) == 0 ) && strcmp( "P", parse2 ) != 0){
							ent_t1 = n;
							found2 = true;
						}
					}
					if( found2 == false ){
						Type = Type2;
						ent_t1 = G.newNode();
						Type2.init(G);
						Type2[ ent_t1 ] = 0;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != ent_t1 )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						NodeArray<Entities> ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != ent_t1 ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_t1 ].entityName, "%s", parse2 );
					}
					int moreWeHave = 0;
					if( strcmp( parse3, "+" ) == 0 ){
						myReadFile >>  parse4 ;
						myReadFile >>  parse5 ;
						//fscanf( fptr, "%s%s", parse4, parse5 );
// 						cout << parse4 << "\t" <<  parse5 << "\n";
						found2 = false;
						forall_nodes( n, G ){
							if( ( strcmp( ENTITIES[ n ].entityName, parse4) == 0 || strcmp( ENTITIES[ n ].reactionName, parse4) == 0 )&& strcmp( "P", parse4 ) != 0 != 0 ){
								ent_t2 = n;
								found2 = true;
							}
						}
						if( found2 == false ){
							Type = Type2;
							ent_t2 = G.newNode();
							Type2.init(G);
							Type2[ ent_t2 ] = 0;
							forall_nodes( walkOnNodes, G){
								if( walkOnNodes != ent_t2 )
									Type2[ walkOnNodes ] = Type[ walkOnNodes ];
							}
							Type = Type2;
							NodeArray<Entities> ENTITIEST = ENTITIES;
							ENTITIES.init( G );
							forall_nodes( walkOnNodes, G ){
								if( walkOnNodes != ent_t2 ){
									ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
								}
							}
							sprintf( ENTITIES[ ent_t2 ].entityName, "%s", parse4 );
						}
						moreWeHave = 1;
					}

					if( moreWeHave != 1 ){
						sprintf( ENTITIES[ ent_t1 ].reactionName, "%s", rname );
						e = G.newEdge( ent_s1, ent_t1 );
						if( typeNow == 1 )
							rtemp.source1ToComplex = e;
						e = G.newEdge( ent_s2, ent_t1 );
						if( typeNow == 1 )
							rtemp.source2ToComplex = e;
						reactionAmbigous.pushBack( ent_t1 );
						if( typeNow == 1 ){
								rtemp.s1 = ent_s1;
								rtemp.s2 = ent_s2;
								rtemp.r1 = ent_t1;
								sprintf( rtemp.complex, "%s", ENTITIES[ ent_t1 ].entityName );
								sprintf( rtemp.source1, "%s", ENTITIES[ ent_s1 ].entityName );
								sprintf( rtemp.source2, "%s", ENTITIES[ ent_s2 ].entityName );
								sprintf( rtemp.typeName, "%s", parse5 );
								sprintf( rtemp.type, "%s", "2+1" );
								rtemp.typeInt = getReactionType( parse5 );
						}
					}
					else{
						Type = Type2;
						n = G.newNode();
						Type2.init(G);
						Type2[ n ] = typeNow;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != n )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != n ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_t1 ].reactionName, "%s", rname );
		
						if( moreWeHave == 1 ){
							Gn[ n ] = getReactionType( parse5 );
							fprintf( fptr2, "%d\n", getReactionType( parse5 ) );
							rtemp.typeInt = getReactionType( parse5 );
						}
						else{
							Gn[ n ] = getReactionType( parse3 );
							fprintf( fptr2, "%d\n", getReactionType( parse3 ) );							
							rtemp.typeInt = getReactionType( parse3 );
						}
						e = G.newEdge( ent_s1, ent_t1 );
						if( typeNow == 1 )
							rtemp.source1ToComplex = e;
						e = G.newEdge( ent_s2, ent_t1 );
						if( typeNow == 1 )
							rtemp.source2ToComplex = e;
						e = G.newEdge( n, ent_t1 );
						if( typeNow == 1 )
							rtemp.complexToTarget1 = e;
						e = G.newEdge( n, ent_t2 );		
						if( typeNow == 1 )
							rtemp.complexToTarget2 = e;
						if( typeNow == 1 ){
								rtemp.s1 = ent_s1;
								rtemp.s2 = ent_s2;
								rtemp.t1 = ent_t1;
								rtemp.t2 = ent_t2;
								rtemp.r1 = n;
								sprintf( rtemp.complex, "%s", ENTITIES[ n ].entityName );
								sprintf( rtemp.source1, "%s", ENTITIES[ ent_s1 ].entityName );
								sprintf( rtemp.source2, "%s", ENTITIES[ ent_s2 ].entityName );
								sprintf( rtemp.target1, "%s", ENTITIES[ ent_t1 ].entityName );
								sprintf( rtemp.target2, "%s", ENTITIES[ ent_t2 ].entityName );
								sprintf( rtemp.typeName, "%s", parse5 );
								sprintf( rtemp.type, "%s", "2+2" );
						}
					}
				}
				else{
					myReadFile >>  parse4 ;
					myReadFile >>  parse5 ;
					//fscanf( fptr, "%s%s", parse4, parse5 );
					node ent_s1,ent_s2,ent_t1,ent_t2;
					bool found2 = false,found1 = false;
					forall_nodes( n, G ){
						if( strcmp( ENTITIES[ n ].entityName, parse2) == 0 || strcmp( ENTITIES[ n ].reactionName, parse2) == 0  ){
							ent_s1 = n;
							found1 = true;
						}
						if( strcmp( ENTITIES[ n ].entityName, parse4) == 0 || strcmp( ENTITIES[ n ].reactionName, parse4) == 0 ){
							ent_t1 = n;
							found2 = true;
						}
					}
					
					if( found1 == false ){
						Type = Type2;
						ent_s1 = G.newNode();
						Type2.init(G);
						Type2[ ent_s1 ] = 0;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != ent_s1 )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						NodeArray<Entities> ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != ent_s1 ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_s1 ].entityName, "%s", parse2 );
					}
	
					if( found2 == false ){
						Type = Type2;
						ent_t1 = G.newNode();
						Type2.init(G);
						Type2[ ent_t1 ] = 0;
						forall_nodes( walkOnNodes, G){
							if( walkOnNodes != ent_t1 )
								Type2[ walkOnNodes ] = Type[ walkOnNodes ];
						}
						Type = Type2;
						NodeArray<Entities> ENTITIEST = ENTITIES;
						ENTITIES.init( G );
						forall_nodes( walkOnNodes, G ){
							if( walkOnNodes != ent_t1 ){
								ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
							}
						}
						sprintf( ENTITIES[ ent_t1 ].entityName, "%s", parse4 );
					}
	
					fprintf( fptr2, "%d\n", getReactionType( parse5 ) );
					Type = Type2;
					n = G.newNode();
					Type2.init(G);
					Type2[ n ] = typeNow;
					forall_nodes( walkOnNodes, G){
						if( walkOnNodes != n )
							Type2[ walkOnNodes ] = Type[ walkOnNodes ];
					}
					Type = Type2;
					NodeArray<Entities> ENTITIEST = ENTITIES;
					ENTITIES.init( G );
					forall_nodes( walkOnNodes, G ){
						if( walkOnNodes != n ){
							ENTITIES[ walkOnNodes] = ENTITIEST[ walkOnNodes ];
						}
					}
					sprintf( ENTITIES[ n ].entityName, "%s", rname );
					e = G.newEdge( ent_s1, n );
					if( typeNow == 1 )
						rtemp.source1ToComplex = e;
					e = G.newEdge( n, ent_t1 );
					if( typeNow == 1 )
						rtemp.complexToTarget1 = e;
					if( typeNow == 1 ){
						rtemp.s1 = ent_s1;
						rtemp.t1 = ent_t1;
						rtemp.r1 = n;
						sprintf( rtemp.complex, "%s", ENTITIES[ n ].entityName );
						sprintf( rtemp.source1, "%s", ENTITIES[ ent_s1 ].entityName );
						sprintf( rtemp.target1, "%s", ENTITIES[ ent_t1 ].entityName );
						sprintf( rtemp.type, "%s", "1+1" );
						rtemp.typeInt = getReactionType( parse5 );
					}
				}
				satir++;
			}
			RList.pushBack( rtemp );
		}
		sprintf( rname2, "X" );
	}
	myReadFile.close();
	fclose( fptr2 );
	G.writeGML("outputs/output.gml");
	return ENTITIES;
}

void produceExperiment1Inputs( int S1, int S2 ){

	int dim1 = 20;
	int dim2 = 20;
	double ratio = 0.0;
	for( int i = 0; i < 10; i++ ){
		for( int j = 0; j <= 5; j++ ){
			for( int k = 0; k < 50; k++ ){
				char arry[128];
				sprintf( arry, "my/experiment_%dx%d_%dx%d_%lf_%d", dim1, dim2, S1, S2, ratio, k );
				cliqueCreator clq;
				clq.createClique( arry, dim1, dim2, S1, S2, ratio );
				//bimaxMain( "my_bimax.txt", 100 );
			}
			ratio += 0.01;
		}
		ratio = 0.0;
		dim1 += 10;
		dim2 += 10;
	}
}

void produceExperiment2Inputs( int S1, int S2 ){
	int dim1 = 10;
	int dim2 = 10;
	double ratio = 0.0;
	for( int i = 0; i < 5; i++ ){
		for( int j = 0; j <= 5; j++ ){
			for( int k = 0; k < 50; k++ ){
				char arry[128];
				sprintf( arry, "my/experiment_%dx%d_%dx%d_%lf_%d", dim1, dim2, S1, S2, ratio, k );
				cliqueCreator clq;
				clq.createClique( arry, dim1, dim2, S1, S2, ratio );
				//bimaxMain( "my_bimax.txt", 100 );
			}
			ratio += 0.01;
		}
		ratio = 0.0;
		dim1 += 5;
		dim2 += 5;
	}
}

int main()
{		
    produceExperiment2Inputs( 50, 50 );
	//bimaxMain( "my/experiment_100x100_200x200_0.050000_46_bimax.txt", 100 );
	Graph G;
	node n,v,source,target;
	edge e;
	NodeArray<int> Gn( G );
	NodeArray<Entities> ENTITIES;
	NodeArray<int> Type;
	List<node> reactAmbig;
	List<REACTION> reactionsList;
	List<CONTINGENCY> contingencyList;
	List<node> Layer1, Layer2, TempN;
	List<int> Layer1Index, Layer2Index;

	ENTITIES = simpleMIMMLParserReactionAndContingency( "inputs/input.txt", G, Gn, Type, reactAmbig, reactionsList );
	List<CONTINGENCY> addedNodes;
	List<CONTINGENCY> contingencyList2;
	List<edge> addedContingencies = simpleMIMMLParserContingency( "inputs/input2.txt", G, Gn, Type, ENTITIES, addedNodes, contingencyList );
	/*forall_nodes( n, G ){
		cout << ENTITIES[ n ].entityName << endl;
	}*/
	contingencyList2 = contingencyList;
	int limit = contingencyList2.size();
	FILE *fptr;
	/*for( int i = 0; i < limit; i++ ){
		CONTINGENCY c = contingencyList2.popFrontRet();
		cout << ENTITIES[ c.s1 ].entityName << "\t" <<  ENTITIES[ c.t1 ].entityName << endl;
	}*/
	NodeArray<int> contingencyNumber( G );
	forall_nodes( n, G )
		contingencyNumber[ n ] = 0;
	contingencyList2 = contingencyList;
	int count1 = 0, count2 = 0;
	for( int i = 0; i < limit; i++ ){
		CONTINGENCY c = contingencyList2.popFrontRet();
		TempN = Layer1;
		bool found = false;
		int size1 = TempN.size();
		for( int j = 0; j < size1; j++ ){
			n = TempN.popFrontRet();
			if( c.s1 == n ){
				found = true;
			}
		}			
		if( found == false ){
			Layer1.pushBack( c.s1 );
			Layer1Index.pushBack( count1 );
			count1++;
		}
		found = false;			
		TempN = Layer2;
		int size2 = TempN.size();
		for( int j = 0; j < size2; j++ ){
			n = TempN.popFrontRet();
			if( c.t1 == n ){
				found = true;
			}
		}
		if( found == false ){
			Layer2.pushBack( c.t1  );
			Layer2Index.pushBack( count2 );
			count2++;
		}
		contingencyNumber[ c.s1 ]++;
		contingencyNumber[ c.t1 ]++;
	}
	/*forall_nodes( n, G ){
		cout << contingencyNumber[ n ] << "\t";
	}*/
	Array2D<int> adjMatrix( 0, Layer1.size(), 0, Layer2.size(), 0 );
	contingencyList2 = contingencyList;
	int cplexEdgeCount = 0;
	for( int i = 0; i < limit; i++ ){
		count1 = -1, count2 = -1;
		CONTINGENCY c = contingencyList2.popFrontRet();
		TempN = Layer1;
		int size1 = TempN.size();
		for( int j = 0; j < size1; j++ ){
			n = TempN.popFrontRet();
			if( n == c.s1 ){
				source = n;
				count1 = j;		
			}
		}
		TempN = Layer2;
		int size2 = TempN.size();
		for( int j = 0; j < size2; j++ ){
			n = TempN.popFrontRet();
			if( n == c.t1 ){
				target = n;
				count2 = j;
			}
		}
		//cout << count1 << "\t" << count2 << " | ";
		adjMatrix( count1, count2 ) = 1;
		cplexEdgeCount++;
		/*count1 = 0, count2 = 0;
		TempN = Layer2;
		for( int j = 0; j < TempN.size(); j++ ){
			n = TempN.popFrontRet();
			if( n == c.s1 ){
				source = n;
				count1 = i;		
			}
		}
		TempN = Layer1;
		for( int j = 0; j < TempN.size(); j++ ){
			n = TempN.popFrontRet();
			if( n == c.t1 ){
				target = n;
				count2 = j;
			}
		}
		adjMatrix( count2, count1 ) = 1;
		cout << count2 << "\t" << count1 << endl;*/
	}
	fptr = std::fopen( "adjConting.txt", "w" );
	TempN = Layer2;
	int size2 = TempN.size();
	if( FILETYPE ){
		fprintf( fptr, "probset\t" );
		for( int j = 0; j < size2; j++ ){
			n = TempN.popFrontRet();
			if( j != size2 - 1 )
				fprintf( fptr, "%s\t", ENTITIES[ n ].entityName );
			else
				fprintf( fptr, "%s\n", ENTITIES[ n ].entityName );
		}
	}
	else{
		fprintf( fptr, "%d %d 2 2\n", Layer1.size(), Layer2.size() );
	}
	TempN = Layer1;
	for( int i = 0; i < Layer1.size(); i++ ){
		n = TempN.popFrontRet();
		if( FILETYPE )
			fprintf( fptr, "%s\t", ENTITIES[ n ].entityName );
		for( int j = 0; j < Layer2.size(); j++ ){
			if( j != Layer2.size() - 1 )
				fprintf( fptr, "%d\t", adjMatrix( i, j ) );
			else{
				if( i != Layer1.size() - 1 ){
					fprintf( fptr, "%d\n", adjMatrix( i, j ) );
				}
				else{
					fprintf( fptr, "%d", adjMatrix( i, j ) );
				}
			}
		}
	}
	bimaxMain( "adjConting.txt", 100 );
	fclose( fptr );
	fptr = std::fopen( "cplexInput.txt", "w" );
	fprintf( fptr, "%d %d %d 1 1\n", Layer1.size(), Layer2.size(), cplexEdgeCount );
	for( int i = 0; i < Layer1.size(); i++ ){
		for( int j = 0; j < Layer2.size(); j++ ){
			if( adjMatrix( i, j ) == 1 ){
				fprintf( fptr, "%d %d\n", i, j );
			}
		}
	}
	fclose( fptr );


	cout << "****************************************\n";
	cout << "Orthogonal Begins\n";
	cout << "****************************************\n";
	//node n;
	//edge e;
	GraphAttributes GA(G,
		GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
		GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
		GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
		GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate
	);
	int emax = 0;
	forall_nodes( n, G ){
		List<edge> elist;
		G.outEdges( n, elist );
		if( elist.size() > emax ){
			emax = elist.size();
			v = n;
		}
	}
	forall_nodes( n, G ){		
		int count_adj = 0;		
		GA.labelNode(n) = String( ENTITIES[ n ].entityName );
		forall_adj_edges( e, n)
			count_adj++;
		if( n == v ){
			GA.height(n) = 100.0;
			GA.width(n) = 1400.0;
		}
		else{
			//GA.height(n) = 5.0 * count_adj;
			GA.height(n) = 25.0;
			if( count_adj < 10 )
				GA.width(n) = 20.0 * count_adj;
			else
				GA.width(n) = 300.0;
		}
	}
	GA.writeGML("outputs/orto_org.gml");
	PlanarizationLayout pl;
 
	FastPlanarSubgraph *ps = new FastPlanarSubgraph;
	ps->runs(500);
	VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;
	ves->removeReinsert(EdgeInsertionModule::rrNone);
	pl.setSubgraph(ps);
	pl.setInserter(ves);

	OrthoLayout *ol = new OrthoLayout;
	//ol->separation(20.0);
	//ol->cOverhang(5.0);
	//ol->setOptions(2+2);
	ol->separation(15.0);
    ol->cOverhang(1.0);
    ol->setOptions(2+2);

	pl.setPlanarLayouter(ol);
	pl.call(GA);
	GA.writeGML("outputs/orto.gml");
	return 0;
}