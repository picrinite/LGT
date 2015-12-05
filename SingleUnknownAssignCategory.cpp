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



const int NUM_BLAST_HITS = 150;

unordered_map<string, string>  mpParameters;   // load parameters in parameters.txt file 
unordered_set<string> stCloseHost, stMicrobeNotEnv, stEnvNotMicrobe, stMicrobeAndEnv;    //list of tax id for different categories
vector<string> vecCategory(64, "Mixture");


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


	initParameters(argv[1], mpParameters);

	setCategory(vecCategory, mpParameters["assign_blast_category"]);

	buildSet( stCloseHost, mpParameters["taxonomy_CloseHost"] );
	buildSet( stMicrobeNotEnv, mpParameters["taxonomy_EnvNotMicrobe"] );
	buildSet( stEnvNotMicrobe, mpParameters["taxonomy_MicrobeAndEnv"] );
	buildSet( stMicrobeAndEnv, mpParameters["taxonomy_MicrobeNotEnv"] );

	modify( argv[2] );

	return 0;
}

