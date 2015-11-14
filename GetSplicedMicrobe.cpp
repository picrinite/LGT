#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
using namespace std;

bool isPaired ( string &titleA, string &titleB ){
	int i;
	char ch_A, ch_B;
	i = 0;
	ch_A = titleA[ i ];
	ch_B = titleB[ i ];
	while( ch_A != '_' || ch_B != '_' ){
		if( ch_A != ch_B )
			return false;
		++i;
		ch_A = titleA[ i ];
		ch_B = titleB[ i ];
	}
	return true;
}

// As blast ouput file has no query sequence, this function find its raw sequence given read title 
// @param inputFileFastq1 and inputFileFastq2 are raw fastq files
// @param inputFile is a table whose first column is Fastq title
void modify( char *inputFileFastq1, char *inputFileFastq2, char *inputFile, char* outputFile){
	ifstream inFastq1, inFastq2, in;
	ofstream out;

	inFastq1.open( inputFileFastq1 );
	inFastq2.open( inputFileFastq2 );
	in.open( inputFile );
	out.open( outputFile );

	string readLine, readFastq, line2, line3, line4;
	string name, str;
	bool isMatch;
	int pos;
	stringstream ss;

	unordered_map<string, string> mp;
	// read titles
	getline( in, readLine );
	out  << "microbe_read_length\tmicrobe_read\t" << readLine << endl;
	while( getline( in, readLine ) ){
		pos = readLine.find_first_of( "\t" );
		name = "@" + readLine.substr( 0, pos );
		mp[ name ] = readLine;
	}
	while( getline( inFastq1, readFastq ) ){		
		getline( inFastq1, line2 );
		getline( inFastq1, line3);
		getline( inFastq1, line4 );
		if( mp.find( readFastq ) != mp.end() ){
			ss.str("");
			ss << line2.size() << "\t" << line2 << "\t" << mp[ readFastq ];
			mp[ readFastq ] = ss.str();
		}
	}
	while( getline( inFastq2, readFastq ) ){		
		getline( inFastq2, line2 );
		getline( inFastq2, line3);
		getline( inFastq2, line4 );
		if( mp.find( readFastq ) != mp.end() ){
			ss.str("");
			ss << line2.size() << "\t" << line2 << "\t" << mp[ readFastq ];
			mp[ readFastq ] = ss.str();
		}
	}
	for( auto it = mp.begin(); it != mp.end(); ++it )
		out << it->second << endl;
	in.close();
	out.close();
	inFastq1.close();
	inFastq2.close();
}
//
// g++ -std=c++0x -o GetSplicedMicrobe GetSplicedMicrobe.cpp

// ./GetSplicedMicrobe ../7_LJ_RNA_Human_CAGATC_F_1.fastq ../7_LJ_RNA_Human_CAGATC_F_2.fastq PHostUMicrobe.add_host.P1P2.txt PHostUMicrobe.add_microbe.P1P2P3.txt
int main(int argc, char *argv[]){
	//string inputFileName( "sample.txt" );
	char *inputFileFastq1, *inputFileFastq2,  *inputFile;
	inputFileFastq1 = argv[1];
	inputFileFastq2 = argv[2];
	inputFile = argv[3];
	modify( inputFileFastq1, inputFileFastq2, inputFile, argv[4] ) ;
	return 0;
}

