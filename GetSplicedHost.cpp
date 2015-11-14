#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

const int TOPHAT_MAPPING_THRESHOLD = 10;     // MAPQ = -10 log10 p      , p is an estimate of the probability that the alignment does not correspond to the read's true point of origin
										    //  1 - p is the probability the alignment is unique.   Mapping quality is related to "uniqueness." 
											//  MAPQ = 10 corrresponds to p = 0.1
vector<string> split( string & str, string separator ); 

// Take argv[1] as reads aligned to host and argv[2] as reads aligned to microbe
// generate output combined file

// @param titleA and titleB must contain a '_' character
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

/*
The fields are in the order as below:
1. readname  2. sum of flags  3. reference chromosome name  4. position 5. Mapping quality score( a low score means the read can be aligned to multiple position on ref sequence 
6. CIGAR string( e.g., 51M means 51 matched, 41M677N10M means 677 bases are skipped region from the reference( like intron? ) )  7.Mate read's reference chromsome name, * means mate is not aligned 8. mate position  
9. Inferred fragment size ( 0 if $7 is * ) 10. read sequence 11. read quality 12. Optional fields
*/
const int HOST_READ = 0;
const int SUM_FLGS = 1;
const int REF_CHROMOSOME = 2;
const int MAPPING_POSITION = 3;
const int MAPPING_SCORE = 4;
const int CIGAR = 5;
const int RAW_SEQUENCE = 9;

struct SamHit{
	string sum_flags;
	string ref_chromosome;
	string mapping_position;
	string mapping_score;
	string cigar;
	SamHit(){};
	SamHit( vector<string> & vec_hit ){
		sum_flags = vec_hit [ SUM_FLGS ];
		ref_chromosome = vec_hit [ REF_CHROMOSOME ];
		mapping_position = vec_hit [ MAPPING_POSITION ];
		mapping_score = vec_hit [ MAPPING_SCORE ]; 
		cigar = vec_hit[ CIGAR ];
	}
	string getOutput(){
		return sum_flags + "," + ref_chromosome + "," + mapping_position + "," + mapping_score + "," 
			+ cigar + "#" ;
	}
};

int modify( char *inputFileHost, char *inputFileMicrobe, char* outputFile){
	ifstream inMicrobe, inHost;
	ofstream out, outHostFasta;
	inHost.open( inputFileHost );
	inMicrobe.open( inputFileMicrobe );
	out.open( outputFile );
	outHostFasta.open( "outHostFasta.fasta" );

	string readHost, readMicrobe, raw_sequence, host_read, title;
	vector<string> vec_hits;
	bool paired;
	bool found;
	vector<string> vecColumns;	
	SamHit hit;
	readHost = "";
	getline( inMicrobe, title );
	title = title + "\thost_read\thost_read_length\thost_read_sequence\t#Tophat_Mapping\tList";
	out << title << endl;

	while( getline( inMicrobe, readMicrobe ) ){
		found = false;
		vec_hits.clear();
		paired = isPaired( readHost, readMicrobe );   // readHost is initialized as empty, but is set in inner loop as the row after last matched host read
		if( paired ){  //  score >= 30
			found = true;
			host_read = vecColumns[ HOST_READ ];
			raw_sequence = vecColumns[ RAW_SEQUENCE ];	
			if ( atoi( vecColumns[ MAPPING_SCORE ].c_str() ) >= TOPHAT_MAPPING_THRESHOLD )
				vec_hits.push_back( hit.getOutput() );				
		}
		// it is always possible to find readMicrobe's paired end	
		while( getline( inHost, readHost ) ){
			vecColumns = split( readHost, "\t" );
			hit = SamHit( vecColumns );
			paired = isPaired( readHost, readMicrobe );
			if( !paired ){
				if( found )  // readHost is not paired with readMicrobe and readMicrobe's pair end has been found
					break;
				else    // else keep searching until paired is true
					continue;
			}
			found = true;
			host_read = vecColumns[ HOST_READ ];
			raw_sequence = vecColumns[ RAW_SEQUENCE ];	
			if( atoi( vecColumns[ MAPPING_SCORE ].c_str() ) >= TOPHAT_MAPPING_THRESHOLD )  //  score >= 10
				vec_hits.push_back( hit.getOutput() );		
		}
		// readMicrobe's pairs have been found
		out << readMicrobe << "\t" << host_read <<  "\t" << raw_sequence.size() << "\t" << raw_sequence << "\t"  
			<< vec_hits.size() << "\t" ;
		if( vec_hits.size() == 1 )
			out << vec_hits[0] ;
		else
			out << "\\" ;
		out << endl;
		// output as fasta file for further blast
		outHostFasta << ">" << host_read << endl;
		outHostFasta << raw_sequence << endl; 
	}
	out.close();
	outHostFasta.close();
	inMicrobe.close();
	inHost.close();
	return count;
}

// Given a paired end read from microbe identified by blast, 
//  (1) find its paired host read identified by Tophat and add host read sequence into table as one additional collumn. 
//      Filter with threshold MAPQ >= 10
//  (2) generate a fasta file for these host reads

// g++ -std=c++0x -o GetSplicedHost GetSplicedHost.cpp
// argv[1] is SingleHost.sam, argv[2] is inputfile, argv[3] is outputfile
// ../GetSplicedHost CAGATC_SingleHost.sam PHostPMicrobe.blastout.P1.txt PHostUMicrobe.add_host.P1P2.txt
int main(int argc, char *argv[]){
	//string inputFileName( "sample.txt" );
	char *inputFileHost, *inputFileMicrobe;
	inputFileHost = argv[1];
	inputFileMicrobe = argv[2];
	cout << "The number of spliced DNA fragments with unique host mapping: " << endl;
	cout << modify( inputFileHost, inputFileMicrobe, argv[3] ) << endl;
	return 0;
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