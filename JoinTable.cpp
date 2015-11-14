#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include "Utils.h"
using namespace std;

// inputFile1 is PHostUMicrobe.host_gene.P1P2P3P4.txt, nmae is in the third column , inputFile2 is PHost.blastout.P21.txt, name is the first column
void modify( char *inputFile1, char *inputFile2, char* outputFile){
	ifstream input1, input2;
	ofstream out;

	input1.open( inputFile1 );
	input2.open( inputFile2 );

	out.open( outputFile );

	string readLineA, readLineB, titleA, titleB;
	string name, outputTitle;
	int pos, i, index_host_read;

	unordered_map<string, string> mp;
	// read titles
	getline( input1, titleA );
	getline( input2, titleB );

	outputTitle = titleA + "\t" + "host_read" + "\t#Bac_hits" + "\t#Unique_taxID" + "\tList_Bac"
								"\t#Env_hits" + "\t#Unique_taxID" + "\tList_Env"
								"\t#Human_hits" + "\t#Unique_taxID" + "\tList_Human"
								"\t#Other_hits" + "\t#Unique_taxID" + "\tList_Other" ;

	// output title line

	out << outputTitle << endl;

	vector<string> vecColumns;
	vecColumns = split( titleA, "\t" );
	for( i = 0; i < vecColumns.size() ; ++i )
		if( vecColumns[i] == "host_read" ){
			index_host_read = i;
			break;
		}
	cout << vecColumns[ index_host_read ] << "  index is :" << index_host_read << endl;
	// read in table 2 into a map
	while( getline( input2, readLineB ) ){
		pos = readLineB.find_first_of( "\t" );
		name = readLineB.substr( 0, pos );
		mp[ name ] = readLineB;
	}

	while( getline( input1, readLineA) ){	
		out << readLineA ;	
		
		vecColumns = split( readLineA, "\t" );
		name = vecColumns[ index_host_read ];
		if( mp.find( name ) != mp.end() ){
			out << "\t" << mp[ name ] << endl;
		}
		else {
			for( i = 0; i < 13; ++i )
				out << "\t\\" ;
			out << endl;		
		}
	}
	out.close();
	input1.close();
	input2.close();
}

// g++ -std=c++0x -o JoinTable JoinTable.cpp

// ./JoinTable PHostUMicrobe.host_gene.P1P2P3P4.txt PHost.blastout.P21.txt PHostUMicrobe.combine_table.P1P2P3P4P5.txt
int main(int argc, char *argv[]){
	//string inputFileName( "sample.txt" );
	char *inputFile1, *inputFile2;
	inputFile1 = argv[1];
	inputFile2 = argv[2];
	modify( inputFile1, inputFile2, argv[3] ) ;
	return 0;
}

