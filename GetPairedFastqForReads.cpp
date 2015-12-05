#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "Utils.h"
using namespace std;

unordered_map<string, string> mpParameters;   // load parameters in parameters.txt file 
unordered_map<string, string> mpTargetReads;

void buildMap(unordered_map<string, string> & my_map, string fileName){   // 
//  key :   MG00HS08:584:C5YUBACXX:6:1101:1471:94901_2:N:0:CAGATC     value:  MicrobeNotEnv
	ifstream fin(fileName);
	string sLine, sRead, sCategory;
	while(getline(fin, sLine)){
		if(sLine.empty())
			continue;
		int pos = sLine.find_first_of('\t');
		if(pos != string::npos){
			sRead = sLine.substr(0, pos);
			sCategory = sLine.substr(pos + 1);
			my_map[sRead] = sCategory;
		}
	}
	fin.close();
	cout << my_map.size() << endl;
}

string getDatabaseRow(vector<string> &fastqRead,  string pairLabel, string & sCategory){  // pairLabel = 1  or 2
	int pos = fastqRead[0].find_first_of('_');
	string sRead = fastqRead[0].substr(1, pos - 1);
	vector<string> vecColumns = split(sRead, ":");
	string sLine,sReadID;
	sReadID = mpParameters["sample_id"] + ":" + fastqRead[0].substr(1);
	sLine = sLine + sReadID + "\t" + mpParameters["sample_id"] + "\t";
	for(int i = 0; i < 5; ++i)
		sLine = sLine + vecColumns[i] + "\t";
	sLine = sLine + vecColumns[5] + ":" + vecColumns[6] + "\t" + pairLabel + "\t" + fastqRead[1] + "\t"
		+ to_string(fastqRead[1].length()) + "\t" + fastqRead[3] + "\t" + sCategory;
	return sLine;
}

// Retrieves raw fastq info from input files if that sequence is in the unordered_set<string> stReads
int getFastqReads(string inputFile1, string inputFile2, string outputFile){
	ifstream inFastq1(inputFile1), inFastq2(inputFile2);
	ofstream out(outputFile);
	string sLine1, sLine2, sRead1, sRead2, sCategory1, sCategory2;
	int pos;
	bool bTargetRead1, bTargetRead2;
	if(mpTargetReads.find("MG00HS08:584:C5YUBACXX:6:1101:1471:94901") != mpTargetReads.end())
		cout << "exist in map"  << endl;
	while(getline(inFastq1, sLine1)){
		if(sLine1.empty())
			break;
		vector<string> vecRead1, vecRead2;
		getline(inFastq2, sLine2);
		sRead1 = sLine1.substr(1);
		sRead2 = sLine2.substr(1);	
		bTargetRead1 = mpTargetReads.find(sRead1) != mpTargetReads.end();
		bTargetRead2 = mpTargetReads.find(sRead2) != mpTargetReads.end();
		if(bTargetRead1 || bTargetRead2){	
			vecRead1.push_back(sLine1);
			vecRead2.push_back(sLine2);
			int cnt = 3;
			while(cnt--){
				getline(inFastq1, sLine1);
				vecRead1.push_back(sLine1);
				getline(inFastq2, sLine2);
				vecRead2.push_back(sLine2);
			}
			sCategory1 = bTargetRead1 ? mpTargetReads[sRead1] : "Host";
			sCategory2 = bTargetRead2 ? mpTargetReads[sRead2] : "Host";
			out << getDatabaseRow(vecRead1, "1", sCategory1) << endl;		
			out << getDatabaseRow(vecRead2, "2", sCategory2) << endl;
			continue;
		}
		int cnt = 3;
		while(cnt--){
			getline(inFastq1, sLine1);
			getline(inFastq2, sLine2);
		}			
	}
	inFastq1.close();
	inFastq2.close();
	out.close();
}

// Given a paired end read from some categories identified by blast, 
//  (1) find its paired host read identified by Tophat 

// g++ -std=c++0x -o GetPairedFastqForReads GetPairedFastqForReads.cpp
// ../GetPairedFastqForReads ../parameters.txt
int main(int argc, char *argv[]){
	initParameters(argv[1], mpParameters);
	buildMap(mpTargetReads, mpParameters["reads_tobe_loaded_into_database"] );

	getFastqReads(mpParameters["fastq_input_1"], mpParameters["fastq_input_2"], "DatabaseFastq.txt");

	return 0;
}
