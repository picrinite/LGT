#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <map>
using namespace std;

int getLength( const string &str_quality, const int Q, const int max_low_Q, const int max_consec ){
	int count_low; // record num of bases with low quality
	int consec_low; // record num of consecutive bases with low quality
	bool prev;  // record if previous base is low quality
	int i,len;
	count_low = 0;
	consec_low = 0;
	len = str_quality.length();
	for( i = 0; i < len; ++i ){
		if( str_quality[i] < Q ){	
			++consec_low;	
			++count_low;
			if( consec_low == max_consec )
				return i + 1 - consec_low;
			if( count_low > max_low_Q )
				return 0;
		}
		else
			consec_low = 0;
	}
	return len;
}

void filter( char *inputfile_pair_A, char *inputfile_pair_B, char *outputfile_pair_A, char *outputfile_pair_B, int Q, int max_low_Q, int max_consec, int min_length  ){
	string A_title, A_secquence, A_third, A_quality;
	string B_title, B_secquence, B_third, B_quality;
	int A_len, B_len;
	ifstream in_A, in_B;
	ofstream out_A, out_B;
	in_A.open( inputfile_pair_A );
	in_B.open( inputfile_pair_B );
	out_A.open( outputfile_pair_A );
	out_B.open( outputfile_pair_B );
	Q = Q + 33; // ascii value
	while( in_A.good() ){
		getline(in_A,  A_title);
		if( A_title == "" )
			break;
		getline(in_A,  A_secquence);
		getline(in_A,  A_third);
		getline(in_A,  A_quality);
		getline(in_B,  B_title);
		getline(in_B,  B_secquence);
		getline(in_B,  B_third);
		getline(in_B,  B_quality);
		A_len = getLength ( A_quality , Q, max_low_Q, max_consec );
		B_len = getLength ( B_quality , Q, max_low_Q, max_consec );
		if ( !A_len || !B_len )
			continue;
		if( A_len < min_length || B_len < min_length )
			continue;
		out_A << A_title << endl; 
		out_A << A_secquence.substr( 0, A_len ) << endl;
		out_A << A_third << endl;
		out_A << A_quality.substr( 0, A_len ) << endl;
		out_B << B_title << endl; 
		out_B << B_secquence.substr( 0, B_len ) << endl;
		out_B << B_third << endl;
		out_B << B_quality.substr( 0, B_len ) << endl;
	}
	in_A.close();
	in_B.close();
	out_A.close();
	out_B.close();
}
// program   filename_pair1.fastq    filename_pair2.fastq  -Q  20  -max_low_Q   5   -max_consecutive 3  -min_length   75
int main( int argc, char *argv[]){
	int i;
	map<string, char*> mp_argv;
	map<string, char*>::iterator it;
	i = 5;
	while( argv[i] ){
		mp_argv.insert( make_pair( argv[i], argv[i + 1] ) );
		i += 2;
	}
	int Q; // Q-score less than Q will be counted as low quality
	int max_low_Q;  // a sequence with more than max_low_Q low quality bases will be discarded
	int max_consec;  // a sequence with max_consecutive low quality bases will be truncated
	int min_length; // a truncated read smaller than min_length will be discarded
	it = mp_argv.find( string("-Q") );
	if( it != mp_argv.end() )
		Q = atoi( it->second );
	else 
		Q = 20;
	it = mp_argv.find( string("-max_low_Q") );
	if( it != mp_argv.end() )
		max_low_Q = atoi( it->second );
	else 
		max_low_Q = 3;
	it = mp_argv.find( string("-max_consec") );
	if( it != mp_argv.end() )
		max_consec = atoi( it->second );
	else 
		max_consec = 3;
	it = mp_argv.find( string("-min_length") );
	if( it != mp_argv.end() )
		min_length = atoi( it->second );
	else 
		min_length = 75;
	filter( argv[1], argv[2], argv[3], argv[4], Q, max_low_Q, max_consec, min_length );
	/*for( map<string,char*>::iterator it = mp_argv.begin(); it != mp_argv.end(); ++it ){
	cout << it->first << " , " << it->second << endl;
	}*/
	return 0;
}