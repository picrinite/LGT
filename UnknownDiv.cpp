#include <iostream>
#include <fstream>
#include <string>
#include "Util.h"
using namespace std;
// Divide unknown reads after tophat into
// PairedUnknown & Single Unknown


void modify( char *input, char* barcode ){
	string titleA,seqA,titleB,seqB;
	ifstream in;
	ofstream outSingle, outPair;
	in.open( input );
	string str_bar( barcode );
	string paired_unknown;
	string single_unknown;
	paired_unknown = str_bar + "_PairedUnknown.fasta";
	single_unknown = str_bar + "_SingleUnknown.fasta";
	outPair.open( &paired_unknown[0] );
	outSingle.open( &single_unknown[0] );

	bool paired = true;
	while( in.good() ){
		if( paired ){
			getline(in,titleA);
			if ( titleA.find_first_of(":N:") == string::npos )  // header line, not a title line
				continue;
			getline(in,seqA);
			paired = false;
			continue;
		}
		if( !paired && in.good() ){
			getline(in,titleB);
			if (titleB == "")
				break;
			getline(in,seqB);
		}
		paired = compare( titleA, titleB );
		if ( paired ) {
			outPair << titleA << endl << seqA << endl;
			outPair << titleB << endl << seqB << endl;
		}
		else {
			outSingle << titleA << endl << seqA << endl;
			titleA = titleB;
			seqA = seqB;
		}
	}
	if ( !paired )
		outSingle << titleA << endl << seqA << endl;
}
// argument  argv[1] is unmapped_sorted.fasta  , argv[2] is barcode
int main(int argc, char *argv[]){
	//string inputFileName( "sample.txt" );
	modify( argv[1], argv[2] );
	return 0;
}