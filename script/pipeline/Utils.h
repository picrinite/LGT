#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
#include <unordered_map>
#include <cassert>
using namespace std;

// return NC_000008.11 from gi|568815590|ref|NC_000008.11|
string getSeqName( string & str ){
	int start = str.find( "|N") + 1 ;
    if (start == string::npos)  
        return "";	
    int after_end = str.find_first_of( "|", start );
	if (after_end == string::npos)  
        return "";
    return str.substr( start, after_end - start );
}

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

const int SAM_QUERY_ID = 0;  // the index of query id in sam file
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
const int SAM_OPTIONAL_FIELDS = 11;

struct SamHit{
	string read_id;
	string sum_flags;
	string ref_chromosome;
	string chromosome_num;
	string mapping_position;
	string mapping_score;
	string cigar;
	string optional_fields;
	string getChromosomeNum(string & ref_chromosome){  //  gi|568815587|ref|NC_000011.10|
		string sAccession = getSeqName(ref_chromosome);
		int pos = sAccession.find_first_of('.');
		string sNum;
		if(pos != string::npos && pos >= 2)
			sNum = sAccession.substr(pos - 2, 2);
		else
			return "";
		return sNum[0] == '0' ? string(1, sNum[1]) : sNum;
	}
	SamHit(){};
	SamHit( vector<string> & vec_hit ){
        assert(vec_hit.size() >= SAM_OPTIONAL_FIELDS + 1);
		read_id = vec_hit[ HOST_READ ];
		sum_flags = vec_hit [ SUM_FLGS ];
		ref_chromosome = vec_hit [ REF_CHROMOSOME ];
		chromosome_num = getChromosomeNum(ref_chromosome);
		mapping_position = vec_hit [ MAPPING_POSITION ];
		mapping_score = vec_hit [ MAPPING_SCORE ]; 
		cigar = vec_hit[ CIGAR ];
		optional_fields = vec_hit[SAM_OPTIONAL_FIELDS];
	}
	string getDatabaseRow(){
		return read_id + "\t" + sum_flags + "\t" + ref_chromosome + "\t" + chromosome_num + "\t"
		 + mapping_position + "\t" + mapping_score + "\t" + cigar;
	}
};

void initParameters(char * file_parameters, unordered_map<string, string>  &mpParameters){
	ifstream fin(file_parameters);
	string sLine;
	while(getline(fin, sLine)){
		if( sLine.back() == '\r' )
			sLine.pop_back();
		int pos = sLine.find_first_of('=');
		if(pos != string::npos)
			mpParameters[sLine.substr(0, pos)] = sLine.substr(pos + 1);
	}
	fin.close();
}


vector<string> split( string & str, string separator ){
	int len_separator = separator.length();
	int i,j;
	vector<string> res;
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
// @param find index of string in vecotr<string>
int getIndexOfColumn( vector<string> &vecColumns, string sColumnTitle ){
    for ( int i = 0; i < vecColumns.size(); ++i ){
        if ( vecColumns[i] == sColumnTitle )
        	return i;
    }
    return -1;
}


// @ titleA and titleB are ended with _1 or _2
bool compareFastqTitle ( string &titleA, string &titleB ){
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

//return score of NM field in sam output, number of mismatch
int getNM(string & readLine){
    int start = readLine.find("NM:i:");
    if(start != string::npos){
        start += 5;
        int end = readLine.find_first_not_of("0123456789", start);
        string sNM;
        if(end != string::npos)
            sNM = readLine.substr(start, end - start);
        else
            sNM = readLine.substr(start);
        if(sNM.empty())
            return 0;
        return stoi(sNM);
    }
    return 0;
}

//whether a sam hit is labled as mapped to host, SumFlags 0x4 position is set to 1 when unaligned
//when a low --score-min is used to control mapping rate, check nMisMatch in NM field
bool isMappedSam(SamHit &samhit, int nMaxMisMatch) {
    int nSumFlags = stoi(samhit.sum_flags);
    int nNM = getNM(samhit.optional_fields);
    return !(nSumFlags & 1 << 2) && nNM <= nMaxMisMatch; // mapped and #mismatch <= max_allowed
}


#endif
