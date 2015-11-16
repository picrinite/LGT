#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include "Utils.h"

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

//int countReadsNotComplete = 0; // reads with at least 150 best hits
//int countNotCompleteWithoutBac = 0;  // human reads that don't complete
//int countNotCompleteWithBac = 0;  // bac reads that don't complete
int COUNT[ NUM_BLAST_HITS ] ;

unordered_map<string, string>  mpParameters;   // load parameters in parameters.txt file 
unordered_set<string> stCloseHost, stMicrobeNotEnv, stEnvNotMicrobe, stMicrobeAndEnv;    //list of tax id for different categories
vector<string> vecCategory(64, "Mixture");

void initParameters(char * file_parameters){
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

void setCategory(vector<string> &vecCategory, string assign_category_file){
	ifstream fin(assign_category_file);
	string sLine;
	vector<string> vecColumns;
	getline(fin, sLine);
	for(int i = 0; i < 64; ++i){
		getline(fin, sLine);
		if( sLine.back() == '\r' )
			sLine.pop_back();
		vecColumns = split(sLine, "\t");
		if(vecColumns.size() >= 7)
			vecCategory[i] = vecColumns[6];
	}
	fin.close();
}

void buildSet( unordered_set<string> & my_set, string fileName ){
	ifstream fin( fileName );
	string taxID, record;
	while( getline( fin, taxID ) ){
		if( taxID.empty() )
			break;
		if( taxID.back() == '\r' )
			taxID.pop_back();
		if( taxID == record )
			continue;
		record = taxID;
		my_set.insert( taxID );
	}
	fin.close();
}

string assignCategory(vector<string> & vecHits){
	int nCategory = 0;
	int size = vecHits.size();
	if(!size)
		return vecCategory[nCategory];
	vector<int> vecCategoryCnt(6);
	string sMaxBitScore, sCurrBitScore, sTaxID;
	string sHostTaxID = mpParameters["HostTaxonomyID"];
	vector<string> vecColumns;
	vecColumns = split(vecHits[0], "\t");
	sMaxBitScore = vecColumns[BIT_SCORE];
	for(int i = 0; i < size; ++i){
		vecColumns = split(vecHits[i], "\t");
		sCurrBitScore = vecColumns[BIT_SCORE];
		if(sCurrBitScore != sMaxBitScore){
			vecHits.resize(i);
			break;
		}
		sTaxID = vecColumns[SUBJECT_TAX_ID];
		if(sTaxID == sHostTaxID)  // e.g., 9606
			vecCategoryCnt[0] += 1;
		else if(stCloseHost.find(sTaxID) != stCloseHost.end())
			vecCategoryCnt[1] += 1;
		else if(stMicrobeNotEnv.find(sTaxID) != stMicrobeNotEnv.end())
			vecCategoryCnt[2] += 1;
		else if(stEnvNotMicrobe.find(sTaxID) != stEnvNotMicrobe.end())
			vecCategoryCnt[3] += 1;
		else if(stMicrobeAndEnv.find(sTaxID) != stMicrobeAndEnv.end())
			vecCategoryCnt[4] += 1;
		else
			vecCategoryCnt[5] += 1;
	}
	int score = 32;
	for(int i = 0; i < vecCategoryCnt.size(); ++i){
		if(vecCategoryCnt[i])
			nCategory += score;
		score /= 2;
	}
	return vecCategory[nCategory];
}

void modify( char *inputFile ){
	ifstream fin(inputFile);
	ofstream out_category, out_compressed; 	
	out_category.open("SingleUnknownCategory.txt");
	out_compressed.open("SingleUnknownBlastHit.txt");
	string sTitle = "ReadID\tCategory"; // \tHost\tCloseHost\tMicrobeNotEnv\tEnvNotMicrobe\tMicrobeAndEnv\tNotTarget""
	out_category << sTitle << endl;
	string sLine;
	vector<string> vecHits, vecColumns;
	int cnt = 0;
	while(getline(fin, sLine)){
		if(sLine.empty())
			continue;
		if(sLine[0] == '#'){
			if(vecHits.empty())
				continue;
			string sCategory = assignCategory(vecHits);
			if(sCategory != "Host" && sCategory != "CloseHost" && vecHits.size() <= 150){				
				vecColumns =  split(vecHits[0], "\t");
				out_category << vecColumns[QUERY_ID] << "\t"  << sCategory << endl;
				for(string &sHit : vecHits)
					out_compressed << sHit << endl;
			}
			vecHits.clear();
		}
		else
			vecHits.push_back(sLine);
	}
	if(!vecHits.empty()){
	string sCategory = assignCategory(vecHits);
		if(sCategory != "Host" && sCategory != "CloseHost" && vecHits.size() <= 150){				
			vecColumns =  split(vecHits[0], "\t");
			out_category << vecColumns[QUERY_ID] << "\t"  << sCategory << endl;
			for(string &sHit : vecHits)
				out_compressed << sHit << endl;
		}
	}
	fin.close();
	out_compressed.close();
	out_category.close();
}

// # Fields: query id, subject id, subject tax ids, % identity, alignment length, mismatches, evalue, bit score, subject title, s. start, s. end, subject seq
// argv[1] is CGATGT_SingleUnknown.fasta.blastout

// g++ -std=c++0x -o SingleUnknownAssignCategory SingleUnknownAssignCategory.cpp
// argv[1] is parameter files
// argv[2] is blast output file:   ../SingleUnknownAssignCategory ../parameters.txt CAGATC_SingleUnknown.fasta.blastout.0522.mixDB
int main(int argc, char *argv[]){


	initParameters(argv[1]);

	setCategory(vecCategory, mpParameters["assign_blast_category"]);

	buildSet( stCloseHost, mpParameters["taxonomy_CloseHost"] );
	buildSet( stMicrobeNotEnv, mpParameters["taxonomy_EnvNotMicrobe"] );
	buildSet( stEnvNotMicrobe, mpParameters["taxonomy_MicrobeAndEnv"] );
	buildSet( stMicrobeAndEnv, mpParameters["taxonomy_MicrobeNotEnv"] );

	modify( argv[2] );

	/*

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
	*/
	return 0;
}

