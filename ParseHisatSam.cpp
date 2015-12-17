#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include "Utils.h"
using namespace std;


int getNM(string & readLine){
	int start = readLine.find("NM:i:");
	if(start != string::npos){
		cout << readLine.substr(start) << endl;
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
void parseHisatSam(){
	ifstream fin("./CAGATC_out_hisat/CAGATC_hisat_output.sam");
	if(!fin)
		cout << "cannot open file" << endl;
	string readLine;
	vector<string> vecSamColumns;
	vector<int>  vecMappingQ(256);
	vector<int>  vecEditDist(256);
	int nSumFlags;
	string recordName, readName, readID;
	int nNM, nMappingQ;
	int count = 0;
	while(getline(fin, readLine)){
		if(readLine.empty())
			continue;
		if(readLine[0] == '@')
			continue;
		vecSamColumns = split(readLine , "\t");
		SamHit samhit(vecSamColumns);		
		readID = samhit.read_id;
		if(readID == recordName) // multi-read{
			continue;
		else {  // a new read
			nSumFlags = stoi(samhit.sum_flags);
			if(!(nSumFlags & 1 << 2)){   // mapped read
				nMappingQ = stoi(samhit.mapping_score);
				vecMappingQ[nMappingQ]++;
				nNM = getNM(readLine);
				vecEditDist[nNM]++;
				cout << nMappingQ << "  "  << nNM << endl;
			}
			recordName = readID;
		}
		if(count++ == 2000)
			break;
		//int pos = readID.find_first_of('_');
		//readName = pos != string::npos ? readID.substr(0, pos) : readID;	
	}
	for(int i = 0; i < 256; ++i){
		if(vecMappingQ[i] || vecEditDist[i])
			cout << i << ":" << vecMappingQ[i] << ":" << vecEditDist[i] << ",";
		if(i % 10 == 0)
			cout << endl;
	}
	fin.close();
}

// g++ -std=c++0x -o ParseHisatSam ParseHisatSam.cpp
// ../ParseHisatSam ../parameters.txt
int main(int argc, char *argv[] ){
	parseHisatSam();


	return 0;
}
