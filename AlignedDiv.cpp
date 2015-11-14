#include <iostream>
#include <fstream>
#include <string>
using namespace std;
// Divide aligned reads from tophat into
// PairedHost & SingleHost


// return 0 for two fragments, return 1 for paired reads, return 2 for the same read
int compare ( string &titleA, string &titleB ){
	int index_A, index_B, res;
	char ch_A, ch_B;
	bool A_found, B_found; // '_' has been found
	index_A = 0;
	index_B = 0;
	A_found = false;
	B_found = false;
	while( !A_found || !B_found){
		ch_A = titleA[ index_A ];
		ch_B = titleB[ index_B ];
		if( ch_A != ch_B )
			return 0;   // reads from two fragments
		if( ch_A != '_' )
			++index_A;
		else
			A_found = true;
		if( ch_B != '_' )
			++index_B;
		else
			B_found = true;
	}
	if ( titleA[ ++index_A ] == titleB[ ++index_B ] )
		return 2;    // same read
	else
		return 1;    // paired read
}

void modify( char *input, char* barcode ){
	string readA,readB;
	ifstream in;
	ofstream outSingle, outPair;
	in.open( input );
	string str_bar( barcode );
	string paired_host;
	string single_host;
	paired_host = str_bar + "_PairedHost.sam";
	single_host = str_bar + "_SingleHost.sam";
	outPair.open( &paired_host[0] );
	outSingle.open( &single_host[0] );
	int res;
	bool readA_is_paired;
	while( true ){
		getline(in,readA);
		if (readA[0] != '@' )
			break;
	}
	readA_is_paired = false;
	while( in.good() ){
		getline( in,readB );
		if( readB == "" )
			break;
		res = compare( readA, readB );
		if ( res == 2 )
			readA = readA + '\n' + readB ;
		else if ( res == 1 ){
			outPair << readA << endl;
			readA = readB;
			readA_is_paired = true;
		}
		else if ( readA_is_paired ) {
			outPair << readA << endl ;
			readA = readB;
			readA_is_paired = false;
		}
		else {
			outSingle << readA << endl;
			readA = readB;
		}
	}
	if ( readA_is_paired )
		outPair << readA << endl ;
	else
		outSingle << readA << endl ;
}

// first argument is mapped_sorted.sam  , argv[2] is barcode
int main(int argc, char *argv[]){  
	//string inputFileName( "sample.txt" );
	modify( argv[1], argv[2] );
	return 0;
}