#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>

using namespace std;


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
const int NUM_BLAST_HITS = 150;

int countReadsNotComplete = 0; // reads with at least 150 best hits
int countNotCompleteWithoutBac = 0;  // human reads that don't complete
int countNotCompleteWithBac = 0;  // bac reads that don't complete
int COUNT[ NUM_BLAST_HITS ] ;


void buildSet( unordered_set<string> & my_set, char * fileName ){
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


struct Hit{
	string subject_id;
	string subject_tax_id;
	string identity;
	string alignment_length;
	string mismatches;
	string evalue;
	string subject_title;
	string s_start;
	string s_end;
	Hit( vector<string> & vec_hit ){
		subject_id = vec_hit[ SUBJECT_ID ];
		subject_tax_id = vec_hit[ SUBJECT_TAX_ID ];
		identity = vec_hit[ PERCENT_IDENTITY ];
		alignment_length = vec_hit[ ALIGNMENT_LENGTH ];
		mismatches = vec_hit[ MISMATCH ];
		evalue = vec_hit[ EVALUE ];
		subject_title = vec_hit[ SUBJECT_TITLE ];
		s_start = vec_hit[ START ];
		s_end = vec_hit[ END ];
	}
	string getOutput(){
		string res;
		res = subject_id + "," + subject_tax_id + "," + identity + "," + alignment_length + "," 
			+ mismatches + "," + evalue + "," + subject_title + "," + s_start + "," + s_end + "#" ;
		return res;
	}
};

// all best_hits assigned to one group
struct Group{
	string query_id;
	string groupName;
	vector<Hit> best_hits;
	Group( string name ): groupName ( name ) {
	}
	string getAllHits() {
		string output = query_id ;
		if( !best_hits.size() )
			return string( "\\" );
		for( auto it = best_hits.begin(); it != best_hits.end(); ++it )
					output = output + (*it).getOutput() ;	
		return output;
	}
	string getOutput() {
		stringstream ss;
		ss << getSeqCount() << "\t" ;
		ss << getTaxCount() << "\t" ;
		if( groupName != "human" || best_hits.size() <= 10 )
			ss << getAllHits() ;
		else 
			ss << "# of best_hits > 10 and not shown";
		return ss.str();
	}
	int getSeqCount() {
		return best_hits.size();
	}	
	int getTaxCount() {
		unordered_set<string> set_taxID;
		for( auto it = best_hits.begin(); it != best_hits.end(); ++it )
			set_taxID.insert( (*it).subject_tax_id );
		return set_taxID.size();
	}	
	void cleanUp() {
		best_hits.clear();
	}
};



// Keep only best hits, only keep reads with best hits from bacteria

void modify( char *inputFile, unordered_set<string> & bac_set, unordered_set<string> & envi_set, unordered_set<string> & incomplete_reads_set  ){
	ifstream in;
	ofstream out; //_human, out_bac, out_env_non_bac, out_bac_env, out_non_bac_non_env_non_human;  // out_bac is bacterial_not_environmental from ncbi + bac_genome + bac_draft_genome
	// out_bac_env is bacterial_and_environmental
	string str, taxID, bit_score, max_bit_score, query_id;
	bool shouldDiscard = false;
	bool shouldOutput = false;
	int numBestHits = 0;
	vector<string> vecColumns;	

	Group gr_human( "human" );
	Group gr_other( "other" );
	Group gr_bac( "bac" );
	Group gr_env( "env" );
	
	string outputTitle;
	
	outputTitle = outputTitle + "bacterial_read" + "\t#Bac_hits" + "\t#Unique_taxID" + "\tList_Bac"
								"\t#Env_hits" + "\t#Unique_taxID" + "\tList_Env"
								"\t#Human_hits" + "\t#Unique_taxID" + "\tList_Human"
								"\t#Other_hits" + "\t#Unique_taxID" + "\tList_Other" ;

	in.open( inputFile );

	out.open( "PHost.blastout.P21.txt");
	/*
	out_bac.open( "PHostPMicrobe.blastout.P1.txt" );  
	out_env_non_bac.open( "PHostPEnv.blastout.P1.txt" );
	out_non_bac_non_env_non_human.open( "PHostOthers.blastout.P1.txt" );*/

	out << outputTitle << endl;
	/*
	out_non_bac_non_env_non_human << outputTitle << endl;
	out_bac << outputTitle << endl;
	out_env_non_bac << outputTitle << endl;*/
	

	int n = 0;
	while( getline( in, str ) ){
		if( str[0] == '#' ){
			if( shouldOutput ){
				++COUNT[ numBestHits - 1 ] ;
				string output ;

				output = vecColumns[ QUERY_ID ] + "\t" 
				       + gr_bac.getOutput() + "\t" 
					   + gr_env.getOutput() + "\t"
					   + gr_human.getOutput() + "\t" 
					   + gr_other.getOutput(); 

				out <<  output << endl;

				shouldOutput = false;
				shouldDiscard = false;
				max_bit_score = "";
				numBestHits = 0;
				gr_human.cleanUp();
				gr_other.cleanUp();
				gr_bac.cleanUp();
				gr_env.cleanUp();
			}
			continue;
		}
		if( shouldDiscard )
			continue;
		shouldOutput = true;
		// parsing
		vecColumns = split( str, "\t" );


		Hit hit( vecColumns );

		//cout << n++ << endl;

		if( max_bit_score == "" ){
			max_bit_score = vecColumns[ BIT_SCORE ];
			query_id = vecColumns[ QUERY_ID ];
		}
		else if( vecColumns[ BIT_SCORE ] != max_bit_score ){
			shouldDiscard = true;
			continue;
		}
		// one of the best hits, add to map
		++numBestHits;
		if( numBestHits == NUM_BLAST_HITS  )
			shouldDiscard = true;
		taxID = vecColumns[ SUBJECT_TAX_ID ] ;
		if ( taxID  == "9606" )
			gr_human.best_hits.push_back( hit );
		else {
			bool found = false;
			if( bac_set.find( taxID ) != bac_set.end() ) {
				gr_bac.best_hits.push_back( hit );
				found = true;
			}
			if ( envi_set.find( taxID ) != envi_set.end() ) {
				gr_env.best_hits.push_back( hit );
				found = true;
			}
			if ( !found )
				gr_other.best_hits.push_back( hit );
		}
	}
	// end of file clean up
	if( shouldOutput ){
		++COUNT[ numBestHits - 1 ] ;
		string output ;

		output = vecColumns[ QUERY_ID ] + "\t" 
			+ gr_bac.getOutput() + "\t" 
			+ gr_env.getOutput() + "\t"
			+ gr_human.getOutput() + "\t" 
			+ gr_other.getOutput(); 
		
		out <<  output << endl;		
	}
	in.close();
	out.close();
}

// # Fields: query id, subject id, subject tax ids, % identity, alignment length, mismatches, evalue, bit score, subject title, s. start, s. end, subject seq
// argv[1] is CGATGT_SingleUnknown.fasta.blastout, argv[2] is all_bac_taxID.txt, argv[3] is taxonomy_environmental.txt, argv[4] is 
//taxonomy_bac_and_environmental.txt

//  g++ -std=c++0x -o ParseHostBlastOutput ParseHostBlastOutput.cpp

// ./ParseHostBlastOutput outHostFasta.fasta.blastout.0522.mixDB all_bac_taxID.txt taxonomy_environmental.txt
// 



int main(int argc, char *argv[]){
	unordered_set<string> bac_set, envi_set;
	unordered_set<string> incomplete_reads_set;
	buildSet( bac_set, argv[2] );
	buildSet( envi_set, argv[3] );

	ofstream out_result;
	out_result.open( "statistic_analysis.txt");
	out_result << endl;
	modify( argv[1], bac_set, envi_set, incomplete_reads_set );
	
	ofstream out_incomplete( "incomplete_reads.txt" );
	for( auto it = incomplete_reads_set.begin(); it != incomplete_reads_set.end(); ++it )
		out_incomplete << *it << endl;

	for( int i = 0; i < 15; ++i ){
		for( int j = 0; j < 10; ++j ){
			out_result << COUNT[ i * 10 + j ] << "\t" ;
			cout << COUNT[ i * 10 + j ] << "\t" ;
		}
		out_result << endl ;
	}

	cout << endl << "incomplete_reads size is :  " << incomplete_reads_set.size() << endl;
	//out_result << "num of read with at least 150 hits from human : " << countNotCompleteWithoutBac  << endl;
	out_result << "num of read with 150 best hits including bac or env : " << countNotCompleteWithBac << endl;
	//cout << "num of read with at least 150 hits from human : " << countNotCompleteWithoutBac  << endl;
	cout << "num of read with 150 best hits including bac or env : " << countNotCompleteWithBac << endl;
	out_result.close();
	out_incomplete.close();

	return 0;
}

