#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "Utils.h"
using namespace std;

unordered_map<string, string> mpParameters;

//   query id, subject id, subject tax ids, % identity, alignment length, mismatches, evalue, bit score, subject title, s. start, s. end, subject seq

//MG00HS08:584:C5YUBACXX:6:1101:1238:13451_2:N:0:CAGATC   gi|704366570|gb|KM224290.1|     9606    99.01   101     1       4e-43    183    Homo sapiens isolate FRT66 12S ribosomal RNA gene, partial sequence; mitochondrial      125     225     CTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGG

const int QUERY_ID = 0;
const int SUBJECT_ID = 1;
const int SUBJECT_TAX_ID = 2;
const int PERCENT_IDENTITY = 3;
const int ALIGNMENT_LENGTH = 4;
const int MISMATCH = 5;
const int EVALUE = 6;
const int BIT_SCORE = 7;
const int SUBJECT_TITLE = 8;
const int START = 9;
const int END = 10;
const int SUBJECT_SEQUENCE= 11;

struct BlastHit{
	string read_id;
	string gi_no;
	string accession_no;
	string taxon_id;
	string identity;
	string align_length;
	string mismatch;
	string evalue;
	string bitscore;
	string subject_title;
	string subject_start;
	string subject_end;
	string subject_sequence;;
	void setGiAndAccession(string & subject_id){  // gi|665506052|gb|AC255368.1| 
		if(count(subject_id.begin(), subject_id.end(), '|') == 4 && subject_id.substr(0, 3) == "gi|"){
			int pos = subject_id.find_first_of('|', 3);
			if(pos != string::npos){
				gi_no = subject_id.substr(3, pos - 3);
				accession_no = subject_id.substr(pos + 1, subject_id.length() - 1 - (pos + 1));
				return;
			}
		}
		gi_no = "";
		accession_no = "";
	}
	BlastHit(){};
	BlastHit( vector<string> & vec_hit ){
		read_id = vec_hit[QUERY_ID];
		setGiAndAccession(vec_hit[SUBJECT_ID]);
		taxon_id = vec_hit[SUBJECT_TAX_ID];
		identity = vec_hit[PERCENT_IDENTITY];
		align_length = vec_hit[ALIGNMENT_LENGTH];
		mismatch = vec_hit[MISMATCH];
		evalue = vec_hit[EVALUE];
		bitscore = vec_hit[BIT_SCORE];
		subject_title = vec_hit[SUBJECT_TITLE];
		subject_start = vec_hit[START];
		subject_end = vec_hit[END];
		subject_sequence = vec_hit[SUBJECT_SEQUENCE];	
	}
	string getDatabaseRow(){
		return read_id + "\t" + gi_no + "\t" + accession_no + "\t" + taxon_id + "\t"
		 + identity + "\t" + align_length + "\t" + mismatch + "\t"
		 + evalue + "\t" + bitscore + "\t" + subject_title+ "\t" + subject_start + "\t"
		 + subject_end + "\t" + subject_sequence;
	}
};

void getOutput( string input, string output ){
	ifstream in(input);
	ofstream out(output);
	vector<string> vecColumns;
	string readLine;
	string title = "Read_ID\tGI_No\tAccession_No\tTaxon_ID\tIdentity\tAlign_Length\tMismatches\tEvalue\tBitscore\tSubject_Title\tS_Start\tS_End\tS_Seq";
	out << title << endl;
	while( getline( in, readLine ) ) {
		if(readLine.empty())
			continue;
		vecColumns = split(readLine, "\t");
		if(vecColumns.size() > SUBJECT_SEQUENCE){
			BlastHit blasthit(vecColumns);
			out << mpParameters["sample_id"] << ":" << blasthit.getDatabaseRow() << endl;
		}
	}
	in.close();
	out.close();
}

// g++ -std=c++0x -o Blast_CreateDatabaseTable Blast_CreateDatabaseTable.cpp
// ../Blast_CreateDatabaseTable ../parameters.txt
int main(int argc, char *argv[] ){
	initParameters(argv[1], mpParameters);
	getOutput("PairedBlast.txt", "DatabasePairedBlast.txt");

	return 0;
}
