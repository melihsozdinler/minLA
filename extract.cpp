#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/cluster/ClusterPlanarizationLayout.h>
#include <ogdf/cluster/ClusterOrthoLayout.h>
#include <ogdf/cluster/ClusterOrthoShaper.h>
#include <ogdf/energybased/SpringEmbedderFR.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/graphalg/ShortestPathWithBFM.h>

#include <stdio.h>


using namespace ogdf;
 
struct entityNames{
	char entityProp[256];
};
typedef entityNames Type;
typedef entityNames ReactionType;
typedef entityNames ContingecyType;

struct Entities{
	char entityName[256];
	char reactionName[256];
	List<Type> Types;
	List<entityNames> entityProp;
};

struct contingency{
	char source1[256];
	char target1[256];
	char contigName[256];
	node s1;
	node t1;
	node c1;
	edge source1ToComplex;
	edge complexToTarget1;
	char typeName[256];
	int typeInt;
};

struct reaction{
	char source1[256];
	char source2[256];
	edge source1ToComplex;
	edge source2ToComplex;
	edge complexToTarget1;
	edge complexToTarget2;
	node s1;
	node t1;
	node s2;
	node t2;
	node r1;
	char complex[256];
	char target1[256];
	char target2[256];
	char type[10]; // Possible; 1+1, 1+2, 2+1, 2+2;
	char typeName[256];
	int typeInt;
};

typedef contingency CONTINGENCY;
typedef reaction REACTION;

int getReactionType( char *parse ){
						/*  'ncb' : Non-covalent binding (SBO:0000177)
							'cvm' : Covalent modification (SBO:0000210 ! addition of a chemical group)
							'cvb' : Covalent bond (SBO:0000210 ! addition of a chemical group, ?)
							'stc' : Stoichiometric conversion (SBO:0000182 ! conversion)
							'syn' : Production without loss of reactants (SBO:0000205 - biochemical process, ?)
							'tra' : Transcription (SBO:0000183)
							'cle' : Cleavage (SBO:0000178)
							'deg' : Degradation (SBO:0000179)	*/
	ReactionType arrayOfTypes[ 8 ] = { {"(ncb)"},{"(cvm)"},{"(cvb)"},{"(stc)"},{"(syn)"},{"(tra)"},{"(cle)"},{"(deg)"}};
	for( int i = 0; i <= 7; i++ ){
		if( strcmp( arrayOfTypes[ i ].entityProp, parse ) == 0 )
			return i;
	}
	return 0;
}

int getContingecyType( char *parse ){
						/*  'sti' : Stimulation
							'cat' : Catalysis
							'stireq' : Stimulation and requirement
							'inh' : Inhibition 
							'req' : Requirement  
							'enzcat_trans' : Enz Cat Trans.*/
	ContingecyType arrayOfTypes[ 6 ] = { {"(cat)"},{"(enzcat_trans)"},{"(req)"},{"(sti)"},{"(inh)"},{"stireq"}};
	for( int i = 8; i <= 13; i++ ){
		if( strcmp( arrayOfTypes[ i - 8 ].entityProp, parse ) == 0 )
			return i;
	}
	return 0;
}

List<edge> simpleMIMMLParserContingency( char filename[64], Graph &G, NodeArray<int> &Gn, NodeArray<int> &Type, NodeArray<Entities> &ENTITIES, List<CONTINGENCY> &addedNodes, List<CONTINGENCY> &contingencyList ){
	FILE *fptr,*fptr2;
	List<edge> added;
	fptr2 = fopen( "outputs/reactionTypes.txt", "a" );
	//GRAPH<int,int> G;	/* Introducing new graph */
	node n,n2,walkOnNodes;
	edge e = NULL;
	NodeArray<int> Type2;

	if( ( fptr = fopen( filename, "r" ) ) == NULL ){
		cout << "\n The File named " << filename << " does not exist\n" << endl;
	}
	else{
		cout <<  "/* File opened */\n";
		char parse[128],parse2[128],parse3[128],parse4[128],parse5[128],rname[128],rname2[128];
		fscanf( fptr, "%s", parse );
		List<node> NODES;
		List<Entities> ENTY;
		int satir = 1;
		int typeNow = 1;
		while( !feof( fptr ) ){
			CONTINGENCY ctemp;
			fscanf( fptr, "%s", parse );
			if( strcmp( parse, "reactions:" ) == 0 ){
				fscanf( fptr, "%s", parse );
				typeNow = 1;
			}
			if( strcmp( parse, "contingencies:" ) == 0 ){
				fscanf( fptr, "%s", parse );
				typeNow = 2;
			}
			if( strcmp( parse, "-" ) == 0 ){
				fscanf( fptr, "%s%s%s%s%s",rname, parse2, parse3, parse4, rname2 );
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
	return added;
}

NodeArray<Entities> simpleMIMMLParserReactionAndContingency( char filename[64], Graph &G, NodeArray<int> &Gn, NodeArray<int> &Type, List<node> &reactionAmbigous, List<REACTION> &RList ){
	FILE *fptr,*fptr2;
	fptr2 = fopen( "outputs/reactionTypes.txt", "w" );
	//GRAPH<int,int> G;	/* Introducing new graph */
	node n,n2,walkOnNodes;
	NodeArray<Entities> ENTITIES( G );
	NodeArray<Entities> ENTITIEST( G );
	edge e = NULL;
	NodeArray<int> Type2;
// 	ListIterator() it = NULL;

	if( ( fptr = fopen( filename, "r" ) ) == NULL ){
		cout << "\n The File named " << filename << " does not exist\n" << endl;
	}
	else{
		cout <<  "/* File opened */\n";
		char parse[128],parse2[128],parse3[128],parse4[128],parse5[128],rname[128],rname2[128];
		fscanf( fptr, "%s", parse );
		List<node> NODES;
		List<Entities> ENTY;
		int satir = 1;
		int typeNow = 1;
		while( !feof( fptr ) ){	
			REACTION rtemp;
			fscanf( fptr, "%s", parse );
			if( strcmp( parse, "reactions:" ) == 0 ){
			/*  'ncb' : Non-covalent binding (SBO:0000177)
				'cvm' : Covalent modification (SBO:0000210 ! addition of a chemical group)
				'cvb' : Covalent bond (SBO:0000210 ! addition of a chemical group, ?)
				'stc' : Stoichiometric conversion (SBO:0000182 ! conversion)
				'syn' : Production without loss of reactants (SBO:0000205 - biochemical process, ?)
				'tra' : Transcription (SBO:0000183)
				'cle' : Cleavage (SBO:0000178)
				'deg' : Degradation (SBO:0000179)	*/
				fscanf( fptr, "%s", parse );
				typeNow = 1;
			}
			if( strcmp( parse, "contingencies:" ) == 0 ){
				fscanf( fptr, "%s", parse );
				typeNow = 2;
			}
			if( strcmp( parse, "-" ) == 0 ){
				fscanf( fptr, "%s%s%s",rname, parse2, parse3);
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
					fscanf( fptr, "%s%s", parse4, parse5 );
					
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
 				
					fscanf( fptr, "%s%s", parse2, parse3 );
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
						fscanf( fptr, "%s%s", parse4, parse5 );
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
					fscanf( fptr, "%s%s", parse4, parse5 );
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
	fclose( fptr );
	fclose( fptr2 );
	G.writeGML("outputs/output.gml");
	return ENTITIES;
}


int main()
{
    Graph G;
    NodeArray<int> Gn( G );
    NodeArray<Entities> ENTITIES;
    NodeArray<int> Type;
	List<node> reactAmbig;
	List<REACTION> reactionsList;
	List<CONTINGENCY> contingencyList;
    ENTITIES = simpleMIMMLParserReactionAndContingency( "inputs/input.txt", G, Gn, Type, reactAmbig, reactionsList );

    cout << "****************************************\n";
    Array<String> listOfTypes( 10 );
    listOfTypes[ 0 ] = String( "#00FF00" );
    listOfTypes[ 1 ] = String( "#FFD700" );
    listOfTypes[ 2 ] = String( "#6A5ACD" );
    listOfTypes[ 3 ] = String( "#0000FF" );
    listOfTypes[ 4 ] = String( "#FFE4E1" );
    listOfTypes[ 5 ] = String( "#7FFFD4" );
    listOfTypes[ 6 ] = String( "#FA8072" );
    listOfTypes[ 7 ] = String( "#FF5511" );
    listOfTypes[ 8 ] = String( "#CCCCCC" );
    listOfTypes[ 9 ] = String( "#AA00FF" );

    node n,v,x;
    edge e,e2;
    NodeArray<String> labels( G );
	int count = 0;

    GraphAttributes GA(G,
    GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
    GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
    GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
    GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
    cout << "****************************************\n";
    ClusterGraph CG = ClusterGraph(G);
    ClusterGraphAttributes CA(CG,
    GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
    GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
    GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
    GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
    cout << "****************************************\n";
    GA.setAllWidth( 15.0 );
    CA.setAllWidth( 15.0 );
    GA.setAllHeight( 8.0 );
    CA.setAllHeight( 8.0 );

    FILE *fptr;
    fptr = fopen( "outputs/reactionTypes.txt", "r" ); 
    forall_nodes( n, G ){
// 	    cout << " LOOP BEGIN " << endl;
 	    if( Type[ n ] == 1 ){
//  		  GA.shapeNode(n) = 0;
		  int type;
		  fscanf( fptr, "%d", &type );
// 		  cout << type << endl;
		  GA.colorNode(n) = listOfTypes[ type ];
		  CA.colorNode(n) = listOfTypes[ type ];
// 		  GA.colorNode(n) = String( "#FFAA11" );
		  GA.height(n) = 6.0;
		  CA.height(n) = 6.0;
// 		  cout << " type1 " << endl;
	    }
	    else{
		  if( Type[ n ] == 0 ){
			GA.colorNode(n) = listOfTypes[ 8 ];
			GA.width(n) = 7.0;
			GA.height(n) = 7.0;
			CA.colorNode(n) = listOfTypes[ 8 ];
			CA.width(n) = 7.0;
			CA.height(n) = 7.0;
// 			cout << " type2 " << endl;
		  }
		  else{
		      if( Type[ n ] == 2 ){  
				GA.colorNode(n) = listOfTypes[ 9 ];
				GA.width(n) = 7.0;
				GA.height(n) = 7.0;
				CA.colorNode(n) = listOfTypes[ 9 ];
				CA.width(n) = 7.0;
				CA.height(n) = 7.0;
// 				cout << " type3 " << endl;
		      }
		  }
	    }
// 	    cout << " CHECK DONE " << endl;
	    int count_adj = 0;
	    forall_adj_edges( e, n)
		count_adj++;

		if( Type[ n ] == 0 ){
			GA.width(n) = 60.0;
			GA.height(n) = 30.0 * count_adj;
		}
		if( Type[ n ] == 1 ){
			GA.width(n) = 10.0 * count_adj;
			GA.height(n) = 1.0;
		}
		if( Type[ n ] == 2 ){
			GA.width(n) = 10.0 * count_adj;
			GA.height(n) = 1.0;
		}

	    GA.labelNode(n) = String( ENTITIES[ n ].entityName );
	    CA.width(n) = 5.0 * count_adj;
	    CA.labelNode(n) = String( ENTITIES[ n ].entityName );
// 	    cout <<  GA.labelNode(n)  << endl;
// 	    cout << " LOOP DONE " << endl;
//  	    cout << ENTITIES[ n ].entityName << endl;
// 	    cout << GA.labelNode(n) << endl;
    }

    fclose( fptr );
    cout << "****************************************\n";
    cout << "Create AND/OR GRAPH Begins\n";
    cout << "****************************************\n";
    FILE *fptrandor;
    fptrandor = fopen( "outputs/AndOrGraph.txt", "w" ); 
    forall_nodes( n, G ){
	int out = 0, in = 0;
	if( Type[ n ] == 2 ){
		forall_adj_edges( e, n){
			if( strcmp( ENTITIES[ e->target() ].entityName ,  ENTITIES[ n ].entityName ) == 0 )
				out++;
			if( strcmp( ENTITIES[ e->source() ].entityName ,  ENTITIES[ n ].entityName ) == 0 )
				in++;
		}
		if( out > 0 && in > 0 ){
			{
				int out_ = 0;
				forall_adj_edges( e, n){
					if( strcmp( ENTITIES[ e->target() ].entityName ,  ENTITIES[ n ].entityName ) == 0 ){
						if( out_ < out - 1 )
							fprintf( fptrandor, "%s\t+\t", ENTITIES[ e->source() ].entityName );
						else
							fprintf( fptrandor, "%s\t", ENTITIES[ e->source() ].entityName );
						out_++;
					}
				}
			}
			{
				fprintf( fptrandor, "=>\t" );
				int in_ = 0;
				forall_adj_edges( e, n){
					if( strcmp( ENTITIES[ e->source() ].entityName ,  ENTITIES[ n ].entityName ) == 0 ){
						if( in_ < in - 1 )
							fprintf( fptrandor, "%s\t+\t", ENTITIES[ e->target() ].entityName );
						else
							fprintf( fptrandor, "%s", ENTITIES[ e->target() ].entityName );
						in_++;
					}
				}
			}
			fprintf( fptrandor, "\n" );
		}
	}
    }
    fclose( fptrandor );
    

//    cout << "****************************************\n";
//    cout << "Sugiyama Begins\n";
//    cout << "****************************************\n";
//
//
//	forall_nodes( n, G ){
//		GA.width(n) = 60.0;
//		GA.height(n) = 60.0;
//	}
//
//    SugiyamaLayout SL;
//
//    SL.setRanking(new OptimalRanking);
//    SL.setCrossMin(new MedianHeuristic);
// 
//    FastHierarchyLayout *ohl = new FastHierarchyLayout();
//	
//    ohl->layerDistance(300.0);
//    ohl->nodeDistance(150.0);
//
////     ohl->weightBalancing(0.1);
//    SL.setLayout(ohl);
//    SL.call(GA);
//    
//    GA.writeGML("outputs/sugiyama.gml");
// 
//
//    fptr = fopen( "outputs/sugi_layout.txt", "w" );
//
//    count = 1;
//    fprintf( fptr, "%d\n", G.numberOfEdges() );
//    forall_edges( e, G ){
//// 	  cout << e << "\n";
//	  DPolyline lines = GA.bends( e );
//	  if( Type[ e->source() ] == 2 ){
//		for( int i = 0; i < listOfTypes.size(); i++ ){
//			if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//				fprintf( fptr, "%d\t%d\t", lines.size(), i );
//			}
//		}
//	  }
//	  else{
// 		if( Type[ e->source() ] == 1 ){
//			for( int i = 0; i < listOfTypes.size(); i++ ){
//				if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//					fprintf( fptr, "%d\t%d\t", lines.size(), i );
//				}
//			}
//		}
//		else{
//			if( Type[ e->target() ] == 2 ){
//				for( int i = 0; i < listOfTypes.size(); i++ ){
//					if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//						fprintf( fptr, "%d\t%d\t", lines.size(), i );
//					}
//				}
//			}
//			else{
//				if( Type[ e->target() ] == 1 ){
//					for( int i = 0; i < listOfTypes.size(); i++ ){
//						if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//							fprintf( fptr, "%d\t%d\t", lines.size(), i );
//						}
//					}
//				}
//				else{
//					fprintf( fptr, "%d\t%d\t", lines.size(), 1 );
//				}
//			}
//		}
//	  }
//
//	  fprintf( fptr, "%s\t%s\t%lf\t%lf\t", ENTITIES[ e->source() ].entityName, ENTITIES[ e->target() ].entityName, GA.x( e->source() ), GA.y( e->source() ) );
//	  while( lines.empty() != true  ){
//		 DPoint n2 = lines.popFrontRet();
//		 fprintf( fptr, "%lf\t%lf\t", n2.m_x, n2.m_y );
//	  }
//	  fprintf( fptr, "\n" );
//    }
//
//    count = 1;
//    fprintf( fptr, "%d\n", G.numberOfNodes() );
//    forall_nodes( n, G ){
//	  fprintf( fptr, "%d\t%d\t%lf\t%lf\t%s\t%lf\t%lf\n", count, Type[ n ], GA.width(n),
//							     GA.height(n), ENTITIES[ n ].entityName, GA.x( n ), GA.y( n ) );
//	  count++;
//    }
//
//    fclose( fptr );

//     cout << "****************************************\n";
//     cout << "Clustered Orthogonal Begins\n";
//     cout << "****************************************\n";
// 
//     NodeArray<cluster> node2cluster;
// 
//     cluster root = CG.firstCluster();
// 
//     SList<node> cl1,cl2;
// 
// 
//     forall_nodes(v,G){
// // 	cout << CA.labelNode(v) << "  "  << Type[ v ] << endl;
// 	if(Type[ v ]== 0){
// 	  CG.reassignNode(v,root);
// 	}
// 	else{
// 		if(Type[ v ]== 1){
// 		      cl1.pushBack( v );
// 		}
// 		else{
// 			if(Type[ v ]== 2){
// 			    cl2.pushBack( v );
// 			}
// 		}
// 	}
//     }
//     cluster c1 = CG.createCluster(cl1,root);
//     cluster c2 = CG.createCluster(cl2,root);
//     cout << CG.numberOfClusters() << endl;


//     ClusterPlanarizationLayout plan;
//     ClusterOrthoLayout *col = new ClusterOrthoLayout;
//     col->separation(100.0);
//     col->cOverhang(0.5);
//     col->setOptions(2+4);
//     CCLayoutPackModule *ccpm;
//     plan.pageRatio( 1.2 );
//     plan.setPlanarLayouter( col );
//     plan.setPacker( ccpm );
//     cout << "****************************************\n";
//     plan.call(G,CA,CG,true);
// 
//     CA.writeGML("out4.gml");

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
	    forall_adj_edges( e, n)
			count_adj++;
		if( Type[ n ] == 0 ){
			if( n == v ){
				GA.height(n) = 3500;
				GA.width(n) = 50.0;
			}
			else{
				GA.height(n) = 25.0;
				if( count_adj < 10 )
					GA.width(n) = 20.0 * count_adj;
				else
					GA.width(n) = 250.0;
			}
		}
		if( Type[ n ] == 1 ){
			GA.height(n) = 5.0 * count_adj;
			GA.width(n) = 2.0;
		}
		if( Type[ n ] == 2 ){
			GA.height(n) = 5.0 * count_adj;
			GA.width(n) = 2.0;
		}
	}
	int size = reactAmbig.size();
	for( int i = 0; i < size; i++ ){
		n = reactAmbig.popFrontRet();
		GA.height(n) = 20.0;
		GA.width(n) = 20.0;
	}

	List<DRect> boundingRectangles;
	List<DLine> boundingLines;
	forall_nodes( n, G ){
		DRect rect( GA.y( n ) + GA.height( n ) / 2, GA.x( n ) - GA.width( n ) / 2, GA.y( n ) - GA.height( n ) / 2, GA.x( n ) + GA.width( n ) / 2 );
		boundingRectangles.pushBack( rect );
		boundingLines.pushBack( rect.bottomLine() );
		boundingLines.pushBack( rect.leftLine() );
		boundingLines.pushBack( rect.rightLine() );
		boundingLines.pushBack( rect.topLine() );
	}
	List<CONTINGENCY> addedNodes;
	List<edge> addedContingencies = simpleMIMMLParserContingency( "inputs/input2.txt", G, Gn, Type, ENTITIES, addedNodes, contingencyList );

	size = addedNodes.size();
	for( int i = 0; i < size; i++ ){
		CONTINGENCY X = addedNodes.popFrontRet();
		if( X.source1ToComplex->source() != NULL && X.complexToTarget1->target() != NULL ){
			GA.x( X.c1 ) = (double)( GA.x( X.source1ToComplex->source() ) + GA.x( X.complexToTarget1->target() ) ) / 2.0;
			GA.y( X.c1 ) = (double)( GA.y( X.source1ToComplex->source() ) + GA.y( X.complexToTarget1->target() ) ) / 2.0;
			GA.height( X.c1 ) = 5.0;
			GA.width( X.c1 ) = 5.0;
		}
	}

	//size = addedContingencies.size();
	//for( int i = 0; i < size; i++ ){
	//	e = addedContingencies.popFrontRet();
	//	DPolyline lines;
	//	DPoint src;
	//	DPoint tgt;
	//	DPoint tmp;
	//	src.m_x = GA.x( e->source() );
	//	src.m_y = GA.y( e->source() );
	//	tgt.m_x = GA.x( e->target() );
	//	tgt.m_y = GA.y( e->target() );
	//	tmp.m_y = src.m_y >= tgt.m_y ? tgt.m_y + abs( src.m_y - tgt.m_y ) / 2.0 : src.m_y + abs( src.m_y - tgt.m_y ) / 2.0;
	//	tmp.m_x = src.m_x;
	//	lines.pushBack( tmp );
	//	tmp.m_x = tgt.m_x;
	//	lines.pushBack( tmp );
	//	lines.pushBack( tgt );
	//	GA.bends( e ) = lines;
	//	//lines = GA.bends( e );
	//	/*while( lines.empty() != true  ){
	//		 DPoint n2 = lines.popFrontRet();
	//		 printf( "%lf\t%lf\t", n2.m_x, n2.m_y );
	//	}*/
	//	//cout << " done \n";
	//	GA.colorEdge( e ) = listOfTypes[ 1 ];
	//}
	//GA.writeGML("outputs/orthogonal2.gml");

	FILE *connectivity;
	//connectivity = fopen( "connectivity.txt", "w" );
	//fprintf( connectivity, "%s\t", "connectivity" );
	//forall_nodes( n, G ){
	//	fprintf( connectivity, "%s\t",ENTITIES[ n ].entityName );
	//}
	//fprintf( connectivity, "\n" );
	//forall_nodes( n, G ){
	//	fprintf( connectivity, "%s\t", ENTITIES[ n ].entityName );
	//	forall_nodes( v, G ){
	//		bool flag = false;
	//		forall_edges( e, G ){
	//			if( (e->source() == n && e->target() == v ) || (e->source() == v && e->target() == n ) ){
	//				flag = true;
	//				fprintf( connectivity, "%d\t", 1 );
	//				break;
	//			}
	//		}
	//		if( flag == false )
	//			fprintf( connectivity, "%d\t", 0 );
	//	}
	//	fprintf( connectivity, "\n" );
	//}
	//fclose( connectivity );

	int maxNeighborCount = 0;
	int bagSizeForDegree = 5;
	forall_nodes( n, G ){
		int n_count = 0;
		forall_adj_edges( e, n ){
			n_count++;
		}
		if( n_count >= maxNeighborCount )
			maxNeighborCount = n_count;
	}

	Array<List<node>> bagAccordingToDegree( maxNeighborCount / bagSizeForDegree + 1 );
	connectivity = fopen( "connectivity.txt", "w" );
	forall_nodes( n, G ){
		fprintf( connectivity, "%s\t", ENTITIES[ n ].entityName );
		int n_count = 0;
		forall_adj_edges( e, n ){
			n_count++;
		}
		bagAccordingToDegree[ n_count / bagSizeForDegree ].pushBack( n );
		fprintf( connectivity, "%d\t", n_count );
		fprintf( connectivity, "\n" );
	}
	fclose( connectivity );
	
	Array<List<CONTINGENCY> > contingencyBagsDegree( maxNeighborCount / bagSizeForDegree + 1 );
	List<CONTINGENCY> suspects_tmp = contingencyList;
	// Now form contingecy bags according to type
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		List<edge> elist,elist2;
		G.adjEdges( c1.s1, elist );
		G.adjEdges( c1.t1, elist2 );
		//if( elist.size() / bagSizeForDegree == elist2.size() / bagSizeForDegree){
			contingencyBagsDegree[ elist.size() / bagSizeForDegree ].pushBack( c1 );
			contingencyBagsDegree[ elist2.size() / bagSizeForDegree ].pushBack( c1 );
		/*}*/
	}

	suspects_tmp = contingencyList;
	// Now form contingecy bags according to type
	List<node> sources,sourcestemp;
	List<node> targets,targetstemp;
	int maxDegree = 0,maxDegree2 = 0;;
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		List<edge> elist,elist2;
		sources.pushBack( c1.s1 );
		G.adjEdges( c1.s1 , elist );
		if( elist.size() > maxDegree )
			maxDegree = elist.size();
		targets.pushBack( c1.t1 );
		G.adjEdges( c1.t1 , elist2 );
		if( elist2.size() > maxDegree2 )
			maxDegree2 = elist2.size();
	}

	FILE *wptr;
	wptr = fopen( "degreeDist.txt", "w" );
	targetstemp = targets;
	Array2D<int> degreeMatrix( 0, maxDegree, 0, maxDegree2, 0 );
	suspects_tmp = contingencyList;
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		List<edge> elist,elist2;
		G.adjEdges( c1.s1, elist );
		G.adjEdges( c1.t1, elist2 );
		//cout << elist.size() << " - " << elist2.size() << endl;
		degreeMatrix( elist.size(), elist2.size() )++;
	}

	suspects_tmp.clear();
	fprintf( wptr, "\t" );
	for( int j = 0; j < degreeMatrix.high2(); j++ ){
		fprintf( wptr, "%d\t", j );
		cout << j << " ";
	}
	fprintf( wptr, "\n" );
	cout << endl;
	for( int i = 0; i < degreeMatrix.high1(); i++ ){
		fprintf( wptr, "%d\t", i );
		cout << i << " ";
		for( int j = 0; j < degreeMatrix.high2(); j++ ){ 
			fprintf( wptr, "%d\t", degreeMatrix( i, j ) );
			cout << degreeMatrix( i, j ) << " ";
		}
		fprintf( wptr, "\n" );
		cout << endl;
	}
	fclose( wptr );
	wptr = fopen( "degreeDist2.txt", "w" );
	suspects_tmp = contingencyList;
	// Now form contingecy bags according to type
	sources.clear();
	targets.clear();
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		List<edge> elist,elist2;
		sources.pushBack( c1.s1 );
		targets.pushBack( c1.t1 );
	}

	/*************************************************/
	/*   Produce CPLEX input						 */
	/*************************************************/
	
	FILE *cplex;
	cplex = fopen( "cplex.txt", "w" );
	suspects_tmp = contingencyList;
	NodeArray<int> nodeId( G );
	NodeArray<bool> nodeMarked( G );
	forall_nodes( n, G )
		nodeMarked[ n ] = false;
	List<node> stemp = sources, ttemp = targets;
	int count1 = 0;
	while( stemp.empty() != true ){
		n = stemp.popFrontRet();
		if( nodeMarked[ n ] == false ){
			nodeId[ n ] = count1;
			count1++;
			nodeMarked[ n ] = true;
		}
	}
	int count2 = 0;
	while( ttemp.empty() != true ){
		n = ttemp.popFrontRet();
		if( nodeMarked[ n ] == false && Type[ n ] != 9 ){
			nodeId[ ttemp.popFrontRet() ] = count2;
			count2++;			
			nodeMarked[ n ] = true;
		}
	}
	
	fprintf( cplex, "%d %d %d\n", contingencyList.size(), count1, count2 );
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		if( Type[ c1.t1 ] != 9 )
			fprintf( cplex, "%d %d\n", nodeId[ c1.s1 ], nodeId[ c1.t1 ] );
	}
	fclose( cplex );
	/*************************************************/
	
	NodeArray<int> COUNTS( G, 0 );
	sourcestemp = sources;
	targetstemp = targets;
	maxDegree = 0;
	while( sourcestemp.empty() != true ){
		node source = sourcestemp.popFrontRet();
		suspects_tmp = contingencyList;
		int count = 0;
		while( suspects_tmp.empty() != true ){
			CONTINGENCY c1 = suspects_tmp.popFrontRet();
			if( source == c1.t1 || source == c1.s1 )
				count++;
		}
		COUNTS[ source ] = count;
		if( count > maxDegree )
			maxDegree = count;
	}
	maxDegree2 = 0;
	while( targetstemp.empty() != true ){
		node target = targetstemp.popFrontRet();
		suspects_tmp = contingencyList;
		int count = 0;
		while( suspects_tmp.empty() != true ){
			CONTINGENCY c1 = suspects_tmp.popFrontRet();
			if( target == c1.s1 || target == c1.t1 )
				count++;
		}
		COUNTS[ target ] = count;
		if( count > maxDegree2 )
			maxDegree2 = count;
	}
	Array2D<int> degreeMatrix2( 0, maxDegree, 0, maxDegree2, 0 );
	suspects_tmp = contingencyList;
	while( suspects_tmp.empty() != true ){
		CONTINGENCY c1 = suspects_tmp.popFrontRet();
		//cout << elist.size() << " - " << elist2.size() << endl;
		degreeMatrix2( COUNTS[ c1.s1 ], COUNTS[ c1.t1 ] )++;
	}

	suspects_tmp.clear();
	fprintf( wptr, "\t" );
	for( int j = 0; j < degreeMatrix2.high2(); j++ ){
		fprintf( wptr, "%d\t", j );
		cout << j << " ";
	}
	fprintf( wptr, "\n" );
	cout << endl;
	for( int i = 0; i < degreeMatrix2.high1(); i++ ){
		fprintf( wptr, "%d\t", i );
		cout << i << " ";
		for( int j = 0; j < degreeMatrix2.high2(); j++ ){ 
			fprintf( wptr, "%d\t", degreeMatrix2( i, j ) );
			cout << degreeMatrix2( i, j ) << " ";
		}
		fprintf( wptr, "\n" );
		cout << endl;
	}
	fclose( wptr );
	//wptr = fopen( "degreeDistR.txt", "w" );
	//suspects_tmp = contingencyList;
	//// Now form contingecy bags according to type
	//sources.clear();
	//targets.clear();
	//while( suspects_tmp.empty() != true ){
	//	CONTINGENCY c1 = suspects_tmp.popFrontRet();
	//	List<edge> elist,elist2;
	//	if( Type[ c1.t1 ] == 1 ){
	//		sources.pushBack( c1.s1 );
	//		targets.pushBack( c1.t1 );
	//	}
	//}
	//NodeArray<int> COUNTS2( G, 0 );
	//sourcestemp = sources;
	//targetstemp = targets;
	//maxDegree = 0;
	//while( sourcestemp.empty() != true ){
	//	node source = sourcestemp.popFrontRet();
	//	suspects_tmp = contingencyList;
	//	int count = 0;
	//	while( suspects_tmp.empty() != true ){
	//		CONTINGENCY c1 = suspects_tmp.popFrontRet();
	//		if( source == c1.t1 || source == c1.s1 )
	//				count++;
	//	}
	//	COUNTS2[ source ] = count;
	//	if( count > maxDegree )
	//		maxDegree = count;
	//}
	//maxDegree2 = 0;
	//while( targetstemp.empty() != true ){
	//	node target = targetstemp.popFrontRet();
	//	suspects_tmp = contingencyList;
	//	int count = 0;
	//	while( suspects_tmp.empty() != true ){
	//		CONTINGENCY c1 = suspects_tmp.popFrontRet();
	//		if( target == c1.s1 || target == c1.t1 )
	//				count++;
	//	}
	//	COUNTS2[ target ] = count;
	//	if( count > maxDegree2 )
	//		maxDegree2 = count;
	//}
	//Array2D<int> degreeMatrix3( 0, maxDegree, 0, maxDegree2, 0 );
	//suspects_tmp = contingencyList;
	//while( suspects_tmp.empty() != true ){
	//	CONTINGENCY c1 = suspects_tmp.popFrontRet();
	//	//cout << elist.size() << " - " << elist2.size() << endl;
	//	if( Type[ c1.t1 ] == 1 )
	//		degreeMatrix3( COUNTS2[ c1.s1 ], COUNTS2[ c1.t1 ] )++;
	//}

	//suspects_tmp.clear();
	//fprintf( wptr, "\t" );
	//for( int j = 0; j < degreeMatrix3.high2(); j++ ){
	//	fprintf( wptr, "%d\t", j );
	//	cout << j << " ";
	//}
	//fprintf( wptr, "\n" );
	//cout << endl;
	//for( int i = 0; i < degreeMatrix3.high1(); i++ ){
	//	fprintf( wptr, "%d\t", i );
	//	cout << i << " ";
	//	for( int j = 0; j < degreeMatrix3.high2(); j++ ){ 
	//		fprintf( wptr, "%d\t", degreeMatrix3( i, j ) );
	//		cout << degreeMatrix3( i, j ) << " ";
	//	}
	//	fprintf( wptr, "\n" );
	//	cout << endl;
	//}
	//fclose( wptr );


	// Loop over each contingency bag
	for( int i = 0; i < contingencyBagsDegree.size(); i++ ){
		if( contingencyBagsDegree[ i ].empty() != true ){
			suspects_tmp = contingencyBagsDegree[ i ];
			List<node> contingencyNodes;
			ListIterator<CONTINGENCY> coniterate;
			ListIterator<node> nodeiterate;
			Array<int> markedContingencies( 0, suspects_tmp.size(), 0 );
			int k = 0;
			forall_listiterators( CONTINGENCY, coniterate, suspects_tmp ){
				node source = (*coniterate).s1;
				node target = (*coniterate).t1;
				int strue = 0;
				int ttrue = 0;
				forall_listiterators( node, nodeiterate, contingencyNodes ){
					if( (*nodeiterate) == source )
						strue = 1;
					if( (*nodeiterate) == target )
						ttrue = 1;
				}
				if( ttrue == 0 ){
					contingencyNodes.pushBack( source );
					markedContingencies[ k ] = 1;
				}
				if( strue == 0 ){
					contingencyNodes.pushBack( target );
					markedContingencies[ k ] = 1;
				}	
				k++;
			}
			Array2D<int> similarity( 0, contingencyNodes.size(), 0, contingencyNodes.size(), 0 ); // Similarity matrix in order to determine similar contingencies between each bag nodes
			k = 0;
			forall_listiterators( CONTINGENCY, coniterate, suspects_tmp ){
				if( markedContingencies[ k ] == 1 ){
					node source = (*coniterate).s1;
					node target = (*coniterate).t1;
					int l = 0;
					int indexX = 0;
					int indexY = 0;
					forall_listiterators( node, nodeiterate, contingencyNodes ){
						if( (*nodeiterate) == source ){
							indexX = l;
						}
						if( (*nodeiterate) == target ){
							indexY = l;
						}
						l++;
					}
					similarity( indexX, indexY ) = (*coniterate).typeInt;
				}
				k++;
			}
			cout << endl;
			for( int k = 0; k < similarity.high1(); k++ ){
				for( int l = 0; l < similarity.high2(); l++ ){
					cout << similarity( k, l ) << "-";
				}
				cout << endl;
			}
		}
	}







	List<List<node> > toBeMerged;
	 //For each bag of similar degreed nodes do this loop
	for( int i = 0; i < bagAccordingToDegree.size(); i++ ){
		List<node> mergeCandidates; // for storing merge candidates.
		List<node> t_list_n = bagAccordingToDegree[ i ]; // Used for walking on bag nodes
		List<node> t_list_n2 = t_list_n; // Used for walking on bag nodes
		List<node> t_list_n3 = t_list_n; // Used for walking on bag nodes
		//while( t_list_n.empty() != true ){
		//	cout << ENTITIES[ t_list_n.popBackRet() ].entityName << "\t";
		//}
		//cout << "\n****************************************\n";
		Array<List<CONTINGENCY> > contingecies( t_list_n2.size() + 1 ); // For each node inside the bag store its all contingencies.
		Array<List<CONTINGENCY> > contingecies2( t_list_n2.size() + 1 ); // For each node inside the bag store its all contingencies.

		// Extract these contingencies
		for( int k = 0; k < t_list_n2.size(); k++ ){
			List<CONTINGENCY> contingencyList2 = contingencyList;
			if( t_list_n.empty() == true )
				break;
			n = t_list_n.popFrontRet();
			for( int j = 0; j < contingencyList.size(); j++ ){
				CONTINGENCY c1 = contingencyList2.popFrontRet();
				if( c1.t1 == n ){
					contingecies[ k ].pushBack( c1 );
				}
				if( c1.s1 == n ){
					contingecies2[ k ].pushBack( c1 );
				}
			}
		}

		Array2D<int> similarity( 0, contingecies.size(), 0, contingecies.size(), 0 ); // Similarity matrix in order to determine similar contingencies between each bag nodes
		Array2D<int> marked( 0, contingecies.size(), 0, contingecies.size(), 0 ); 

		// Produce similarity matrix
		for( int k = 0; k < t_list_n2.size(); k++ ){
			for( int j = 0; j < t_list_n2.size(); j++ ){
				if( k != j ){
					List<CONTINGENCY> tmp = contingecies[ k ];
					int size1 = tmp.size();
					for( int a = 0; a < size1; a++ ){
						CONTINGENCY c1 = tmp.popFrontRet();
						List<CONTINGENCY> tmp2 = contingecies[ j ];
						int size2 = tmp2.size();
						for( int b = 0; b < size2; b++ ){
							CONTINGENCY c2 = tmp2.popFrontRet();
							if( strcmp( c1.target1, c2.target1 ) == 0 && c1.typeInt == c2.typeInt ){
								similarity( k, j ) += 1;
							}							
						}
						tmp2 = contingecies2[ j ];
						size2 = tmp2.size();
						for( int b = 0; b < size2; b++ ){
							CONTINGENCY c2 = tmp2.popFrontRet();
							if( strcmp( c1.source1, c2.source1 ) == 0 && c1.typeInt == c2.typeInt ){
								similarity( k, j ) += 1;
							}
						}
					}
				}				
			}
		}

		//connectivity = fopen( "similarity.txt", "w" );
		//t_list_n = t_list_n2;
		//cout << endl;
		//fprintf( connectivity, "%s\t", "Similarity" );
		//for( int k = 0; k < t_list_n.size(); k++ ){
		//	n = t_list_n.popFrontRet();
		//	fprintf( connectivity, "%s\t", ENTITIES[ n ].entityName );
		//}
		//t_list_n = t_list_n2;
		//fprintf( connectivity, "\n" );
		//for( int k = 0; k < t_list_n.size(); k++ ){
		//	n = t_list_n2.popFrontRet();
		//	fprintf( connectivity, "%s\t", ENTITIES[ n ].entityName );
		//	for( int j = 0; j < t_list_n.size(); j++ ){
		//		fprintf( connectivity, "%d\t", similarity( k, j ) );
		//		//cout << similarity( k, j ) << "\t";
		//	}
		//	fprintf( connectivity, "\n" );
		//	//cout << endl;
		//	//cout << endl;
		//}
		//fclose(connectivity);
		//cout << endl;

		int max_similarity_i = 0;
		t_list_n = t_list_n2;
		t_list_n3 = t_list_n2;
		node n_t, v_t;
		// Determine max similarity values where it is larger than 3
		for( int k = 0; k < t_list_n.size(); k++ ){
			n = t_list_n2.popFrontRet();
			for( int j = 0; j < t_list_n.size(); j++ ){
				v = t_list_n3.popFrontRet();
				if( similarity( k, j ) >= max_similarity_i && similarity( k, j ) > 3 ){
					max_similarity_i = similarity( k, j );
					n_t = n;
					v_t = v;
				}
			}
			t_list_n3 = t_list_n;
		}	
		mergeCandidates.pushBack( n_t );
		mergeCandidates.pushBack( v_t );
		t_list_n2 = t_list_n;
		t_list_n3 = t_list_n;
		// Determine merge candidates if it is similar to max_similarity value with bagSizeForDegree
		for( int k = 0; k < t_list_n.size(); k++ ){
			n = t_list_n2.popFrontRet();
			for( int j = k; j < t_list_n.size(); j++ ){
				v = t_list_n3.popFrontRet();
				if( (similarity( k, j ) > max_similarity_i - bagSizeForDegree && similarity( k, j ) < max_similarity_i + bagSizeForDegree ) && max_similarity_i > 3 && marked( k, j ) == 0 ){
					mergeCandidates.pushBack( n );
					mergeCandidates.pushBack( v );
					marked( k, j ) = 1;
				}
			}
			t_list_n3 = t_list_n;
		}	
		//toBeMerged.pushBack( mergeCandidates );
		
		t_list_n = mergeCandidates;
		
		//List<node> allNeighbors;
		//while( mergeCandidates.empty() != true ){
		//	n = mergeCandidates.popBackRet();
		//	forall_adj_edges( e, n ){
		//		if( e->source() == n )
		//			allNeighbors.pushBack( e->target() );
		//		else
		//			allNeighbors.pushBack( e->source() );
		//	}
		//}
		//printf( "\n*********************************\n" );
		//while( allNeighbors.empty() != true ){
		//	printf( "%s\t", ENTITIES[ allNeighbors.popBackRet() ].entityName );
		//}
		//printf( "\n*********************************\n" );


		// Find Suspect list from merged candidate list
		List<CONTINGENCY> suspects;
		suspects_tmp.clear();
		Array<int> marked2( contingencyList.size() + 1 );
		for( int i = 0; i < marked2.size(); i++ ){
			marked2[ i ] = 0;
		}
		// Main loop Find Suspect list from merged candidate list
		while( mergeCandidates.empty() != true ){
			n = mergeCandidates.popBackRet();
			List<CONTINGENCY> contingencyList2 = contingencyList;
			// Find Suspect list inner loop to provide uniqueness check that is it appended before.
			for( int j = 0; j < contingencyList.size(); j++ ){
				CONTINGENCY c1 = contingencyList2.popFrontRet();
				if( c1.t1 == n && marked2[ j ] == 0 ){
					//cout << "1" << "\t";
					suspects.pushBack( c1 );
					marked2[ j ]  = 1;
				}
				else{
					//cout << marked2[ j ] << "\t" << "2" << "\t";
					if( c1.s1 == n && marked2[ j ] == 0 ){
						suspects.pushBack( c1 );
						marked2[ j ] = 1;
					}
				}
				//cout << j << "\n";
			}
		}

		// Now we have suspect List, but we do not have idea about contingecy types. 
		suspects_tmp = suspects;
		int max_i = 0;
		// Determine the max size of type of contingecy in suspect list
		while( suspects_tmp.empty() != true ){
			CONTINGENCY c1 = suspects_tmp.popFrontRet();
			if( max_i <= c1.typeInt )
				max_i = c1.typeInt;
		}
		
		Array<List<CONTINGENCY> > contingencyBags( max_i + 1 );
		suspects_tmp = suspects;

		// Now form contingecy bags according to type
		while( suspects_tmp.empty() != true ){
			CONTINGENCY c1 = suspects_tmp.popFrontRet();
			contingencyBags[ c1.typeInt ].pushBack( c1 );
		}
	}


	//				List<node> sources, sourcesComp;
	//				List<node> targets, targetsComp;
	//				while( suspects_tmp.empty() != true ){
	//					CONTINGENCY c1 = suspects_tmp.popFrontRet();
	//					sourcesComp = sources;
	//					int flag = 0;
	//					for( int k = 0; k < sources.size(); k++ ){
	//						n = sourcesComp.popBackRet();
	//						if( strcmp( ENTITIES[c1.s1].entityName, ENTITIES[n].entityName ) == 0 ){
	//							flag = 1;
	//						}
	//					}
	//					if( flag == 0 ){
	//						sources.pushBack( c1.s1 );
	//						cout << ENTITIES[c1.s1].entityName << " |s| ";
	//					}
	//					targetsComp = targets;
	//					flag = 0;
	//					for( int k = 0; k < targets.size(); k++ ){
	//						n = targetsComp.popBackRet();
	//						if( strcmp( ENTITIES[c1.t1].entityName, ENTITIES[n].entityName ) == 0 ){
	//							flag = 1;
	//						}
	//					}
	//					if( flag == 0 ){
	//						targets.pushBack( c1.t1 );
	//						cout << ENTITIES[c1.t1].entityName << " |t| ";
	//					}
	//				}
	//				t_list_n3 = sources;
	//				t_list_n2 = targets;
	//				/*  
	//					For all source nodes and its incoming edges 
	//					update s.t. Combine Node is connected to all 
	//				    and its out edges goes to contingencies that 
	//					are reduced due to that simplification.
	//				*/

	//				// Form grouping node
	//				NodeArray<int> Type2;
	//				node combined = G.newNode();
	//				Type2.init(G);
	//				Type2[ combined ] = 99;
	//				forall_nodes( v, G){
	//					if( v != combined )
	//						Type2[ v ] = Type[ v ];
	//				}
	//				Type = Type2;
	//				NodeArray<Entities> ENTITIEST = ENTITIES;
	//				ENTITIES.init( G );
	//				forall_nodes( v, G ){
	//					if( v != combined ){
	//						ENTITIES[ v] = ENTITIEST[ v ];
	//					}
	//				}
	//				sprintf( ENTITIES[ combined ].reactionName, "%s", "Combined Group" );
	//				sprintf( ENTITIES[ combined ].entityName, "%s", "Combined Group" );

	//				int sourceInedges = 0, targetInedges = 0;
	//				int sourceOutedges = 0, targetOutedges = 0;
	//				for( int k = 0; k < sources.size(); k++ ){
	//					n = t_list_n3.popBackRet();
	//					List<edge> elist;
	//					G.inEdges( n, elist );
	//					sourceInedges += elist.size();
	//					G.outEdges( n, elist );
	//					sourceOutedges += elist.size();
	//				}
	//				for( int k = 0; k < targets.size(); k++ ){
	//					n = t_list_n2.popBackRet();
	//					List<edge> elist;
	//					G.inEdges( n, elist );
	//					targetInedges += elist.size();
	//					G.outEdges( n, elist );
	//					targetOutedges += elist.size();
	//				}
	//				int count = 1;
	//				if( targetOutedges < targetInedges ){
	//					if( sourceInedges < sourceOutedges ){
	//						if( sourceInedges > targetOutedges ){ // Simplify source size
	//							t_list_n3 = sources;					
	//							for( int k1 = 0; k1 < sources.size(); k1++ ){
	//								CONTINGENCY c1;
	//								node c_s = G.newNode();
	//								c1.s1 = combined;
	//								c1.c1 = c_s;
	//								Type2.init(G);
	//								Type2[ c_s ] = 99;
	//								forall_nodes( v, G){
	//									if( v != c_s )
	//										Type2[ v ] = Type[ v ];
	//								}
	//								Type = Type2;
	//								NodeArray<Entities> ENTITIEST = ENTITIES;
	//								ENTITIES.init( G );
	//								forall_nodes( v, G ){
	//									if( v != c_s ){
	//										ENTITIES[ v] = ENTITIEST[ v ];
	//									}
	//								}
	//								sprintf( ENTITIES[ c_s ].reactionName, "%s%d", "C_", count );
	//								sprintf( ENTITIES[ c_s ].entityName, "%s%d", "C_", count );
	//								sprintf( c1.typeName, "%s%d", "C_", count );
	//								sprintf( c1.source1, "%s",  ENTITIES[ combined ].entityName );
	//								count++;
	//								t_list_n2 = targets;
	//								for( int k2 = k1; k2 < targets.size(); k2++ ){
	//									v = t_list_n2.popBackRet();
	//									sprintf( c1.target1, "%s",  ENTITIES[ v ].entityName );
	//									c1.t1 = v;
	//									e = G.newEdge( combined, c_s  );
	//									c1.complexToTarget1 = e;
	//									e = G.newEdge( c_s , v );
	//									c1.source1ToComplex = e;
	//									c1.typeInt = i;
	//									contingencyList.pushBack( c1 );
	//									break;
	//								}
	//								v = t_list_n3.popBackRet();
	//								List<edge> elist;
	//								G.outEdges( v, elist );
	//								int size = elist.size();
	//								for( int k2 = 0; k2 < size; k2++ ){
	//									e = elist.popBackRet();
	//									List<node> nlist = targets;
	//									int flag = 0;
	//									for( int k3 = 0; k3 < sources.size(); k3++ ){
	//										node z = nlist.popBackRet();
	//										List<CONTINGENCY> contingencyList2 = contingencyList;
	//										for( int k4 = 0; k4 < contingencyList.size(); k4++ ){
	//											CONTINGENCY c1 = contingencyList2.popBackRet();
	//											if( c1.t1 == z && c1.c1 == e->target() ){
	//												G.delNode( e->target() );
	//												flag = 1;
	//												break;
	//											}
	//										}
	//									}
	//									if( flag == 0 ){
	//										n = e->target();
	//										G.delEdge( e );
	//										G.newEdge( combined, n );
	//									}
	//								}
	//								G.inEdges( v, elist );
	//								size = elist.size();
	//								for( int k2 = 0; k2 < size; k2++ ){
	//									e = elist.popBackRet();
	//									n = e->source();
	//									G.delEdge( e );
	//									G.newEdge( n, combined );
	//								}
	//							}								
	//						}
	//						else{
	//							t_list_n3 = targets;					
	//							for( int k1 = 0; k1 < targets.size(); k1++ ){
	//								CONTINGENCY c1;
	//								node c_s = G.newNode();
	//								c1.t1 = combined;
	//								c1.c1 = c_s;
	//								Type2.init(G);
	//								Type2[ c_s ] = 2;
	//								forall_nodes( v, G){
	//									if( v != c_s )
	//										Type2[ v ] = Type[ v ];
	//								}
	//								Type = Type2;
	//								NodeArray<Entities> ENTITIEST = ENTITIES;
	//								ENTITIES.init( G );
	//								forall_nodes( v, G ){
	//									if( v != c_s ){
	//										ENTITIES[ v] = ENTITIEST[ v ];
	//									}
	//								}
	//								sprintf( ENTITIES[ c_s ].reactionName, "%s%d", "C_", count );
	//								sprintf( ENTITIES[ c_s ].entityName, "%s%d", "C_", count );
	//								sprintf( c1.typeName, "%s%d", "C_", count );
	//								sprintf( c1.target1, "%s",  ENTITIES[ combined ].entityName );
	//								count++;
	//								t_list_n2 = sources;
	//								for( int k2 = k1; k2 < sources.size(); k2++ ){
	//									v = t_list_n2.popBackRet();
	//									sprintf( c1.target1, "%s",  ENTITIES[ v ].entityName );
	//									c1.s1 = v;
	//									e = G.newEdge( c_s, combined  );
	//									c1.complexToTarget1 = e;
	//									e = G.newEdge( v , c_s );
	//									c1.source1ToComplex = e;
	//									c1.typeInt = i;
	//									contingencyList.pushBack( c1 );
	//									break;
	//								}
	//								v = t_list_n3.popBackRet();
	//								List<edge> elist;
	//								G.inEdges( v, elist );
	//								int size = elist.size();
	//								for( int k2 = 0; k2 < size; k2++ ){
	//									e = elist.popBackRet();
	//									List<node> nlist = sources;
	//									int flag = 0;
	//									for( int k3 = 0; k3 < sources.size(); k3++ ){
	//										node z = nlist.popBackRet();
	//										List<CONTINGENCY> contingencyList2 = contingencyList;
	//										for( int k4 = 0; k4 < contingencyList.size(); k4++ ){
	//											CONTINGENCY c1 = contingencyList2.popBackRet();
	//											if( c1.s1 == z && c1.c1 == e->source() ){
	//												G.delNode( e->source() );
	//												flag = 1;
	//												break;
	//											}
	//										}
	//									}
	//									if( flag == 0 ){
	//										n = e->source();
	//										G.delEdge( e );
	//										G.newEdge( n, combined );
	//									}
	//								}
	//								G.outEdges( v, elist );
	//								size = elist.size();
	//								for( int k2 = 0; k2 < size; k2++ ){
	//									e = elist.popBackRet();
	//									n = e->target();
	//									G.delEdge( e );
	//									G.newEdge( combined, n );
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	forall_nodes( n, G ){		
		int count_adj = 0;
		forall_adj_edges( e, n)
			count_adj++;
		if( Type[ n ] == 0 ){
			if( n == v ){
				GA.height(n) = 2000;
				GA.width(n) = 50.0;
			}
			else{
				GA.height(n) = 25.0;
				if( count_adj < 10 )
					GA.width(n) = 20.0 * count_adj;
				else
					GA.width(n) = 250.0;
			}
		}
		if( Type[ n ] == 1 ){
			GA.height(n) = 5.0 * count_adj;
			GA.width(n) = 2.0;
		}
		if( Type[ n ] == 2 ){
			GA.height(n) = 5.0 * count_adj;
			GA.width(n) = 2.0;
		}
		if( Type[ n ] == 99 ){
			GA.height(n) = 5.0 * count_adj;
			GA.width(n) = 4.0;
		}
		GA.labelNode(n) = String( ENTITIES[ n ].entityName );
	}

	GA.writeGML( "output2.gml" );

	cout << "****************************************\n";
    cout << "Orthogonal Begins\n";
    cout << "****************************************\n";

    PlanarizationLayout pl;
 
    FastPlanarSubgraph *ps = new FastPlanarSubgraph;
cout << " DONE \n";
    ps->runs(1000);
    VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;
cout << " DONE \n";
	ves->removeReinsert(EdgeInsertionModule::rrAll);
    pl.setSubgraph(ps);
    pl.setInserter(ves);
//     EmbedderMinDepthMaxFaceLayers *emb = new EmbedderMinDepthMaxFaceLayers;
//     pl.setEmbedder(emb); 
cout << " DONE \n";

    OrthoLayout *ol = new OrthoLayout;

    ol->separation(10.0);
    ol->cOverhang(10.0);
    ol->setOptions(2+4);

    pl.setPlanarLayouter(ol);
    pl.call(GA);

	GA.writeGML("outputs/orthogonal.gml");
cout << " DONE \n";



//	GraphAttributes GA2(G,
//    GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
//    GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
//    GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
//    GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
//	GA = GA2;
//	
//	fptr = fopen( "outputs/reactionTypes.txt", "r" ); 
//    forall_nodes( n, G ){
//// 	    cout << " LOOP BEGIN " << endl;
// 	    if( Type[ n ] == 1 ){
////  		  GA.shapeNode(n) = 0;
//		  int type;
//		  fscanf( fptr, "%d", &type );
//// 		  cout << type << endl;
//		  GA.colorNode(n) = listOfTypes[ type ];
//// 		  GA.colorNode(n) = String( "#FFAA11" );
//// 		  cout << " type1 " << endl;
//	    }
//	    else{
//		  if( Type[ n ] == 0 ){
//			GA.colorNode(n) = listOfTypes[ 8 ];
//// 			cout << " type2 " << endl;
//		  }
//		  else{
//		      if( Type[ n ] == 9 ){  
//				GA.colorNode(n) = listOfTypes[ 9 ];
//// 				cout << " type3 " << endl;
//		      }
//		  }
//	    }
//	    GA.labelNode(n) = String( ENTITIES[ n ].entityName );
//    }
//    fclose( fptr );





    fptr = fopen( "outputs/ortho_layout.txt", "w" );

    fprintf( fptr, "%d\n", G.numberOfEdges() );
    count = 1;
    forall_edges( e, G ){
// 	  cout << e << "\n";
	  DPolyline lines = GA.bends( e );
	  if( Type[ e->source() ] == 2 ){
		for( int i = 0; i < listOfTypes.size(); i++ ){
			if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
				fprintf( fptr, "%d\t%d\t", lines.size(), i );
			}
		}
	  }
	  else{
 		if( Type[ e->source() ] == 1 ){
			for( int i = 0; i < listOfTypes.size(); i++ ){
				if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
					fprintf( fptr, "%d\t%d\t", lines.size(), i );
				}
			}
		}
		else{
			if( Type[ e->target() ] == 2 ){
				for( int i = 0; i < listOfTypes.size(); i++ ){
					if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
						fprintf( fptr, "%d\t%d\t", lines.size(), i );
					}
				}
			}
			else{
				if( Type[ e->target() ] == 1 ){
					for( int i = 0; i < listOfTypes.size(); i++ ){
						if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
							fprintf( fptr, "%d\t%d\t", lines.size(), i );
						}
					}
				}
				else{
					fprintf( fptr, "%d\t%d\t", lines.size(), 1 );
				}
			}
		}
	  }

	  fprintf( fptr, "%s\t%s\t%lf\t%lf\t", ENTITIES[ e->source() ].entityName, ENTITIES[ e->target() ].entityName, GA.x( e->source() ), GA.y( e->source() ) );
	  while( lines.empty() != true  ){
		 DPoint n2 = lines.popFrontRet();
		 fprintf( fptr, "%lf\t%lf\t", n2.m_x, n2.m_y );
	  }
	  fprintf( fptr, "\n" );
    }

    count = 1;
    fprintf( fptr, "%d\n", G.numberOfNodes() );
    forall_nodes( n, G ){
	  fprintf( fptr, "%d\t%d\t%lf\t%lf\t%s\t%lf\t%lf\n", count, Type[ n ], GA.width(n),
	  GA.height(n), ENTITIES[ n ].entityName, GA.x( n ), GA.y( n ) );
	  count++;
	  //cout << ENTITIES[ n ].entityName << " - " << ENTITIES[ n ].reactionName << endl;
    }
    fclose( fptr );

//    cout << "****************************************\n";
//    cout << "Spring Begins\n";
//    cout << "****************************************\n";
//
//    SpringEmbedderFR FMMMLayout;
//    FMMMLayout.iterations(100);
//    FMMMLayout.minDistCC( 500.0 );
//    FMMMLayout.call( GA );
//
//    GA.writeGML("outputs/spring.gml");
//
//    fptr = fopen( "outputs/spring_layout.txt", "w" );
//
//    fprintf( fptr, "%d\n", G.numberOfEdges() );
//    count = 1;
//    forall_edges( e, G ){
//// 	  cout << e << "\n";
//	  DPolyline lines = GA.bends( e );
//	  if( Type[ e->source() ] == 2 ){
//		for( int i = 0; i < listOfTypes.size(); i++ ){
//			if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//				fprintf( fptr, "%d\t%d\t", lines.size(), i );
//			}
//		}
//	  }
//	  else{
// 		if( Type[ e->source() ] == 1 ){
//			for( int i = 0; i < listOfTypes.size(); i++ ){
//				if( strcmp( GA.colorNode( e->source() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//					fprintf( fptr, "%d\t%d\t", lines.size(), i );
//				}
//			}
//		}
//		else{
//			if( Type[ e->target() ] == 2 ){
//				for( int i = 0; i < listOfTypes.size(); i++ ){
//					if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//						fprintf( fptr, "%d\t%d\t", lines.size(), i );
//					}
//				}
//			}
//			else{
//				if( Type[ e->target() ] == 1 ){
//					for( int i = 0; i < listOfTypes.size(); i++ ){
//						if( strcmp( GA.colorNode( e->target() ).cstr() , listOfTypes[ i ].cstr() ) == 0 ){
//							fprintf( fptr, "%d\t%d\t", lines.size(), i );
//						}
//					}
//				}
//				else{
//					fprintf( fptr, "%d\t%d\t", lines.size(), 1 );
//				}
//			}
//		}
//	  }
//
//	  fprintf( fptr, "%s\t%s\t%lf\t%lf\t", ENTITIES[ e->source() ].entityName, ENTITIES[ e->target() ].entityName, GA.x( e->source() ), GA.y( e->source() ) );
//	  while( lines.empty() != true  ){
//		 DPoint n2 = lines.popFrontRet();
//		 fprintf( fptr, "%lf\t%lf\t", n2.m_x, n2.m_y );
//	  }
//	  fprintf( fptr, "\n" );
//    }
//
//    count = 1;
//    fprintf( fptr, "%d\n", G.numberOfNodes() );
//    forall_nodes( n, G ){
//	  fprintf( fptr, "%d\t%d\t%lf\t%lf\t%s\t%lf\t%lf\n", count, Type[ n ], GA.width(n),
//							     GA.height(n), ENTITIES[ n ].entityName, GA.x( n ), GA.y( n ) );
//	  count++;
//    }
//
//    fclose( fptr );


	//cout << "****************************************\n";
 //   cout << "Produce Interaction Table \n";
 //   cout << "****************************************\n";
 //   FILE *fptrtable;
 //   char chr;
 //   int ecount = 0;
 //   fptrtable = fopen( "outputs/interaction_table.html", "w" );
 //   fptr = fopen( "inputs/table.txt", "r" );
 //   while( !feof( fptr ) ){
 //   	fscanf( fptr, "%c", &chr );
	//fprintf( fptrtable, "%c", chr );
 //   }
 //   fclose( fptr );

 //   forall_nodes( n, G ){
	//if( Type[ n ] == 0 ){
	//	ecount = 0;
	//	{
	//		forall_adj_edges( e, n){
	//			if( Type[ e->target() ] == 1 || Type[ e->target() ] == 2 )
	//				ecount++;	
	//		}
	//	}
	//	fprintf( fptrtable, "<tr> <td rowspan=%d> %s </td>", ecount, ENTITIES[ n ].entityName );
	//	int loop1 = 0;
	//	forall_adj_edges( e, n){
	//		if( loop1 != 0 ){
	//			fprintf( fptrtable, "<tr> " ); 
	//			loop1++;
	//		}
	//		if( Type[ e->target() ] == 1 ){
	//			fprintf( fptrtable, " <td> %s </td> <td>", ENTITIES[ e->target() ].entityName ); 
	//			forall_adj_edges( e2, e->target() ){
	//				if( Type[ e2->target() ] == 0 ){
	//					fprintf( fptrtable, "%s", ENTITIES[ e2->target() ].entityName );
	//				}
	//			}
	//			fprintf( fptrtable, " </td> <td> </td> <td> </td>" );
	//			fprintf( fptrtable, " </tr>\n" );
	//		}
	//		if( Type[ e->target() ] == 2 ){
	//			fprintf( fptrtable, " <td> </td> <td> </td> <td> %s </td> <td>", ENTITIES[ e->target() ].entityName ); 
	//			forall_adj_edges( e2, e->target() ){
 //					if( strcmp( ENTITIES[ e2->target() ].entityName , ENTITIES[ e->target() ].entityName ) != 0 ){
	//					fprintf( fptrtable, "%s", ENTITIES[ e2->target() ].entityName );
 //					}
	//			}
	//			fprintf( fptrtable, " </td>" );
	//			fprintf( fptrtable, " </tr>\n" );
	//		}
	//	}
	//	if( loop1 == 0 )
	//		fprintf( fptrtable, " <td> </td> <td> </td> <td> </td> <td> </td> </tr>\n" );
	//}
 //   }
 //   fprintf( fptrtable, "</table>\n</body>\n</html>" );
 //   fclose( fptrtable );
 //   cout << "****************************************\n";
 //   cout << "Calculating Shortest Paths\n";
 //   cout << "****************************************\n";

 //   NodeArray<int> D(G);
 //   NodeArray<edge> DE(G);
 //   EdgeArray<int> E(G);
 //   forall_edges( e, G ){
	//  E[ e ] = 1;
 //   }

 //   fptrtable = fopen( "outputs/path_table.html", "w" );
 //   FILE *fptr2 = fopen( "inputs/table2.txt", "r" );
 //   fptr = fopen( "outputs/shortestPaths.txt", "w" );
 //   fprintf( fptr, "\t" );
 //   forall_nodes( n, G ){
	//  fprintf( fptr, "%s\t", ENTITIES[ n ].entityName );
 //   }
 //   fprintf( fptr, "\n" );
 //   while( !feof( fptr2 ) ){
	//fscanf( fptr2, "%c", &chr );
	//fprintf( fptrtable, "%c", chr );
 //   }

 //   count = 0;
 //   forall_nodes( n, G ){
	//ShortestPathWithBFM shp = ShortestPathWithBFM();
	//shp.call( G, n, E, D, DE );
	//forall_nodes( v, G ){
	//	if( D[ v ] < G.numberOfEdges() && v != n && Type[ v ] == 0 && strcmp( ENTITIES[ n ].entityName, ENTITIES[ v ].entityName ) != 0 )
	//		count++;
	//}
	//if( Type[ n ] == 0 && count != 0 && strcmp( ENTITIES[ n ].entityName, "P" ) != 0 ){
	//	fprintf( fptrtable, "\n<table id=\"interactions\"> <tr> <th rowspan=2 width=300>Entity Name</th>" );
	//	for( int i = 0; i < count; i++ ){
	//		fprintf( fptrtable, " <th colspan=2 width=400> Path %d </th>", i + 1 );
	//	}
	//	fprintf( fptrtable, " </tr>\n<tr> " );
	//	for( int i = 0; i < count; i++ ){
	//		fprintf( fptrtable, "<td> To </td> <td> Length </td> " );
	//	}
	//	fprintf( fptrtable, " </tr>\n<tr> " );
	//	fprintf( fptrtable, " <td> %s </td> ", ENTITIES[ n ].entityName );
	//	forall_nodes( v, G ){
	//		if( D[ v ] < G.numberOfEdges() && v != n && Type[ v ] == 0 && strcmp( ENTITIES[ n ].entityName, ENTITIES[ v ].entityName ) != 0 )
	//			fprintf( fptrtable, " <td> %s </td> <td> %d </td>", ENTITIES[ v ].entityName, D[ v ] );
	//	}
	//	fprintf( fptrtable, " </tr>\n</table><hr>\n" );
	//}
	//fprintf( fptr, "%s\t", ENTITIES[ n ].entityName );
	//forall_nodes( v, G ){
	//	if( D[ v ] > G.numberOfEdges() )
	//		fprintf( fptr, "%d\t", -1 ); 
	//	else
	//		fprintf( fptr, "%d\t", D[ v ] ); 
	//}  
	//fprintf( fptr, "\n" );
	//count = 0;
 //   } 
 //   fprintf( fptrtable, "\n</body>\n</html>" );
 //   fclose( fptr );
 //   fclose( fptr2 );
 //   fclose( fptrtable );

 //   cout << "****************************************\n";
    return 0;
}