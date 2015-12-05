#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "Utils.h"
using namespace std;


unordered_set<string> stReads;
unordered_map<string, string> mpParameters;

void buildSet(unordered_set<string> & my_set, string fileName){
	ifstream fin(fileName);
	string sLine;
	while(getline(fin, sLine)){
		if(sLine.empty())
			continue;
		int pos = sLine.find_first_of('_');
		if(pos != string::npos)
			my_set.insert(sLine.substr(0, pos));
	}
	fin.close();
}

// Retrieves hits from input file if that hit is in the unordered_set<string> stReads
int getHit(string inputFileHost, string outputFile){
	ifstream inHost(inputFileHost);
	ofstream out(outputFile);
	string sLine, sRead;
	int pos;
	while(getline(inHost, sLine)){
		if(sLine.empty())
			continue;
		pos = sLine.find_first_of('_');
		if(pos != string::npos)
			sRead = sLine.substr(0, pos);
		if(stReads.find(sRead) != stReads.end())
			out << sLine << endl;
	}
	inHost.close();
	out.close();
}

// Given a paired end read from some categories identified by blast, 
//  (1) find its paired host read identified by Tophat 

// g++ -std=c++0x -o GetHitsForPairedReads GetHitsForPairedReads.cpp
// ../GetHitsForPairedReads 
int main(int argc, char *argv[]){
	initParameters(argv[1], mpParameters);

	buildSet(stReads, mpParameters["reads_tobe_loaded_into_database"] );

	getHit("SingleHost.sam", "PairedHostTophat.txt");
	getHit("SingleUnknownBlastHit.txt", "PairedBlast.txt");

	return 0;
}
