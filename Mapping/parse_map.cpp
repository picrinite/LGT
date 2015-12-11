// First Name:        Xing
// Last Name:         Zhang


#include <sys/stat.h>
#include <dirent.h>    // DIR *dirp = opendir( char filename );  dirent * direntp = readdir( dirp );
#include <iostream>
#include <string>
#include <errno.h>
#include <cstdlib>
#include <fstream>
using namespace std;
// example of .rmp file
/*
Accession: NC_015671.1
GI: 336319134
DNA  length = 3526441
Taxname: [Cellvibrio] gilvus ATCC 13127
Taxid: 593907
Genetic Code: 11
Protein count: 3164
CDS count: 3164
Pseudo CDS count: 0
RNA count: 54
Gene count: 3263
Pseudo gene count: 45
Others: 3263
Total: 6486
*/
// parse .rpt file , output to mapping file
void parseFile( string &fileName, ostream &fout ){
	ifstream fin;
	fin.open( fileName.c_str() );
	string str;
	bool accessionFound = false;
	bool taxidFound = false;
	string accessionID;
	string taxID;
	int pos;
	while( fin.good() && ( !accessionFound || !taxidFound ) ){
		getline( fin, str );
		if( str == "" )
			break;
		if( !accessionFound && str.length() > 10 && str.substr( 0, 9 ) == "Accession" ){
			accessionFound = true;
			pos = str.find_last_of( " " );
			accessionID = str.substr( pos + 1 );
			continue;
		}
		if( !taxidFound && str.length() > 6 && str.substr( 0, 5 ) == "Taxid" ){
			taxidFound = true;
			pos = str.find_last_of( " " );
			taxID = str.substr( pos + 1 );
			continue;
		}
	}
	fin.close();
	if( accessionFound && taxidFound )
		fout << accessionID << " " << taxID << endl;
}

// parse directory recursively
void parseFolder( string dirName, ostream & fout ){
	DIR * dirp;
	dirp = opendir( dirName.c_str() );
	if( !dirp ){
		cerr << "Error(" << errno << ") opening " << dirName << endl;
		exit (1);
	}
	dirent *direntp;
	struct stat statbuf;
	while( ( direntp = readdir( dirp ) ) ){
		const string& localName = string( direntp->d_name );
		if(  localName == "." || localName == ".." ){
			direntp = readdir( dirp );
			continue;
		}
		string fileName = dirName + localName;
		stat( fileName.c_str(), &statbuf );
		// directory
		if( S_ISDIR( statbuf.st_mode ) )		
			parseFolder( fileName + "/", fout );
		else 
			parseFile( fileName, fout );
	}
	closedir( dirp );
}
int main( int argc, char * argv[] ){
	if( argc != 3 ){    // no arguments
		cerr << "Please provide input folder name and output mapping file name" << endl;
		exit( 1 );
	}
	struct stat statbuf;
	stat( argv[1], &statbuf );
	ofstream fout;
	fout.open( argv[2] );
	if( S_ISDIR(statbuf.st_mode) ){
		parseFolder( string( argv[1] ) + "/", fout );
		fout.close();
	}
	else {
		cerr << "Input( argv[1] is not a folder name" << endl;
		exit( 1 );
	}
}	
