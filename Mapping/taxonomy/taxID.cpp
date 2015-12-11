#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <iostream>
using namespace std;
// accessioID_taxID file
void buildSet( unordered_set<string> & my_set, char * fileName ){
	ifstream fin;
	fin.open( fileName );
	string str, taxID, record;
	int pos;
	while( fin.good() ){
		getline( fin, str );
		if( str == "" )
			break;
		pos = str.find_last_of( " " );
		taxID = str.substr( pos + 1 );
		if( taxID[ taxID.length() - 1 ] == '\r' )
			taxID = taxID.substr( 0, taxID.length() - 1 );
		if( taxID == record )
			continue;
		record = taxID;
		my_set.insert( taxID );
	}
	fin.close();
}
//list of taxID file
void buildSet2( unordered_set<string> & my_set, char * fileName ){
	ifstream fin;
	fin.open( fileName );
	string taxID, record;
	while( fin.good() ){
		getline( fin, taxID );
		if( taxID == "" )
			break;
		if( taxID[ taxID.length() - 1 ] == '\r' )
			taxID = taxID.substr( 0, taxID.length() - 1 );
		if( taxID == record )
			continue;
		record = taxID;
		my_set.insert( taxID );
	}
	fin.close();
}

// add old taxID into list
/*
General information.
Field terminator is "\t|\t"
Row terminator is "\t|\n"

12	|	74109	|
30	|	29	|
36	|	184914	|
37	|	42	|
46	|	39	|
*/
void addToSet( unordered_set<string> & my_set, char * fileName ){
	ifstream fin;
	fin.open( fileName );
	string str, old_taxID, new_taxID;
	int pos, start;
	while( fin.good() ){
		getline( fin, str );
		if( str == "" )
			break;
		pos = str.find_first_of ( '\t' );
		old_taxID = str.substr( 0, pos );
		start = pos + 3;
		pos = str.find_first_of( '\t', start );
		new_taxID = str.substr( start, pos - start );
		//cout << "new:" << new_taxID << ", old:" << old_taxID << endl;
		if( my_set.find( new_taxID ) != my_set.end() ){  // then this new taxID is bacteria, add old taxID in
			my_set.insert( old_taxID );
			//cout << "entered" << endl;
		}
	}
	fin.close();
}
// argv[1] is file name for accessionID_taxID.txt
// argv[2] is file name for accessionID_taxID_bac_draft.txt
// argv[3] is file name for taxonomy_result.txt  for ncbi
// argv[3] is file name for merged.dmp
//how to run
// ./taxID accessionID_taxID.txt accessionID_taxID_bac_draft.txt taxonomy_result.txt merged.dmp
int main( int argc, char *argv[] ){

	unordered_set<string> bac_set;  // string is taxID
	unordered_set<string> bac_draft_set;
	unordered_set<string> bac_ncbi_set;
	

	buildSet( bac_set, argv[1] );  
	buildSet( bac_draft_set, argv[2] );
	buildSet2( bac_ncbi_set, argv[3] );  // ncbi taxID list
	addToSet( bac_ncbi_set, argv[4] );   // add old taxID into list

	ofstream fout;
	fout.open( "bac_genome_but_not_ncbi.txt" );

	cout << endl << "#of taxID in bac_genomes: " << bac_set.size() << endl;
	cout << "#of taxID in bac_draft_genomes: " << bac_draft_set.size() << endl;
	cout << "#of taxID in bac_ncbi: " << bac_ncbi_set.size() << endl;

	int intsect_bac_and_draft = 0;
	int intsect_bac_and_ncbi = 0;
	int intsect_draft_and_ncbi = 0;
	int intsect_all = 0;
	for( auto it = bac_draft_set.begin(); it != bac_draft_set.end(); ++it ){
		if( bac_ncbi_set.find( *it ) != bac_ncbi_set.end() ){
			++intsect_draft_and_ncbi;
		}
	}
	
	bool bac_ncbi_found, bac_draft_found;
	for( auto it = bac_set.begin(); it != bac_set.end(); ++it ){
		bac_ncbi_found = false;
		bac_draft_found = false;
		if( bac_ncbi_set.find( *it ) != bac_ncbi_set.end() ){
			++intsect_bac_and_ncbi;
			bac_ncbi_found = true;
		}
		else
			fout<< *it << endl;
		if( bac_draft_set.find( *it ) != bac_draft_set.end() ){
			++intsect_bac_and_draft;
			bac_draft_found = true;
		}
		if( bac_ncbi_found && bac_draft_found )
			++intsect_all;
	}

	cout << endl << "intersect of bac genomes and bac_draft_genomes: " <<  intsect_bac_and_draft << endl;
	cout << "intersect of bac genomes and bac_ncbi: " <<  intsect_bac_and_ncbi << endl;
	cout << "intersect of bac_draft genomes and bac_ncbi: " <<  intsect_draft_and_ncbi << endl;
	cout << "intersect of all: " <<  intsect_all << endl;
}