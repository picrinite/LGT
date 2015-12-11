#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <unordered_map>
using namespace std;

// title line in mapping file
//Sample Name;	Animal Type; Barcode Sequence

int lenghBarcode;
int lengthAfterSpace;
const int POS_BARCODE = 3; // barcode is in the 3rd column of mapping file
char chEndType;   // '1' or '2' :   _1.fastq or _2.fastq
void globalVarInit( char * input ){
	int len = strlen( input );
	//@MG00HS08:584:C5YUBACXX:6:1101:6407:32413 1:N:0:AAAAAA
	int pos;
	string strLen( input ); // which end of RNAseq: 1 or 2
	pos = strLen.find_last_of( '.' );
	chEndType = input[ pos - 1 ];
	ifstream fin( input );
	getline( fin, strLen );
	pos = strLen.find_last_of ( ':' );
	lenghBarcode = strLen.size() - 1 - pos;
	pos = strLen.find_last_of ( ' ', pos );
	lengthAfterSpace =  strLen.size() - 1 - pos;
	fin.close();
}

vector<string> split( string & str, string separator ){
	int len_separator = separator.length();
	int i,j;
	vector<string> res;
	bool startFound = false;
	bool isSeparator = false;
	string field;
	for( i = 0; i < str.size() ; ++i ){
		char & ch_c = str[i];
		isSeparator = false;
		for( j = 0; j < len_separator; ++j ){
			if(  ch_c == separator[j] ){
				isSeparator = true;
				break;
			}
		}
		if( !isSeparator )
			field += ch_c;
		else if ( field.length() ){
			res.push_back( field );
			field.clear();
		}
	}
	if( field.length() )
		res.push_back( field );
	return res;
}

void demultiplex( char *input, unordered_map<string, string> & mp ){
	string strLine;
	string strBarcode;
	unordered_map<string, ofstream*>  mpFout;
	ifstream fin( input );
	ofstream foutNotInList( "not_extracted_" + string( input ) ) ;
	int i = 0;
	int pos = 0;
	
	for( auto it = mp.begin(); it != mp.end() ; ++it ){
		//ofstream fout ( it -> second.c_str() );
		//mpFout[ it->first ] = move(fout);
		//mpFout.insert( make_pair( it->first, ofstream ( it -> second.c_str() ) ) );
		mpFout[ it->first ] = new ofstream ( it -> second.c_str() ) ;
	}
	while( fin.good() ){
		getline( fin, strLine );
		if( strLine == "" )
			break;
		pos = strLine.length() - 1 - lengthAfterSpace ;
		if( strLine[ pos ] == ' ' )
			strLine[ pos ] = '_';
		pos = strLine.length() - lenghBarcode ; // index position of barcode
		strBarcode = strLine.substr( pos );
		if( mpFout.find( strBarcode ) != mpFout.end() ) {
			ofstream &fout = *mpFout[ strBarcode ];
			fout << strLine << endl;
			for( i = 0; i < 3 ; ++i ){
				getline( fin, strLine );
				fout << strLine << endl;
			}
		}
		else {
			ofstream & fout = foutNotInList;
			fout << strLine << endl;
			for( i = 0; i < 3 ; ++i ){
				getline( fin, strLine );
				fout << strLine << endl;
			}
		}
	}
	
	for( auto it = mpFout.begin() ; it != mpFout.end() ; ++it )
		(*it->second).close();
	foutNotInList.close();
}

void getBarcodes( char * mapping_file, unordered_map<string, string> & mp ){
	ifstream fin( mapping_file );
	string strLine;
	string fileName;
	string strBarcode;
	int i;
	getline( fin, strLine ); // read off title line Sample Name	Animal Type	Barcode Sequence
	vector<string> vecColumns;
	while ( fin.good() ){
		getline( fin, strLine );
		if( strLine == "" )
			break;
		vecColumns = split( strLine, " \t" );
		fileName = vecColumns[ 0 ];
		strBarcode = vecColumns[ POS_BARCODE - 1 ];
		for( i = 1; i < vecColumns.size(); ++i )
			fileName = fileName + "_" + vecColumns[ i ] ;
		fileName = fileName + "_" + chEndType + ".fastq";
		mp[ strBarcode ] = fileName;
	}
	fin.close();
}
// argv[1] is the file argv[2] is mapping file 
// test0502.txt
//	JB41_Mapping_File_1-11-15.xlsx
int main( int argc, char *argv[]){
	//string inputFileName( "sample.txt" );
	//string inputFileName("flowcell253_lane1_pair1_ACAGTG.fastq");
	//string outputFileName("flowcell253_lane1_pair1_ACAGTG_1.fastq");
	//modify( argv[1], argv[2] );
	globalVarInit( argv[ 1 ] );
	vector<string> res;
	unordered_map<string, string>  mp_barcode_fileName;   // barcode,  output file name
	getBarcodes( argv[ 2 ], mp_barcode_fileName );
	demultiplex( argv[ 1 ], mp_barcode_fileName );

	/*
	for( auto it = mp_barcode_fileName.begin() ; it != mp_barcode_fileName.end(); ++it )
		cout << it -> first << "," << it -> second << endl;
		*/
	return 0;
}