#include <iostream>
#include <fstream>
using namespace std;
/*

example of file "accessionId_human.txt"

>gi|568815597|ref|NC_000001.11| Homo sapiens chromosome 1, GRCh38.p2 Primary Assembly

*/

int main(){
	
	ifstream fin( "accessionId_human.txt" );
	ofstream fout( "accessionID_taxID_human.txt" );
	string line;
	string taxID ("9606");	// human taxID 9606
	string accessionID;
	int pos_start;
	int pos_end;
	while( fin.good() ){
		getline( fin, line );
		if( line == "" )
			break;
		pos_start = line.find_first_of( "NC" );
		pos_end = line.find_first_of ( '|', pos_start );
		accessionID = line.substr ( pos_start, pos_end - pos_start );
		fout << accessionID << " " << taxID << endl;
	}
	fin.close();
	fout.close();
	return 0;
}