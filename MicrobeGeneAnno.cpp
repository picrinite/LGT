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

unordered_map<string, string>  mpParameters;   // load parameters in parameters.txt file 

// NC_000001.11    BestRefSeq      gene    11874   14409   .       +       .       ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;pseudo=true
const int SEQNAME = 0;
const int SOURCE = 1;
const int FEATURE = 2;
const int GFF_START = 3;
const int GFF_END = 4;
const int ATTRIBUTE = 8;
//Read_ID\tGI_No\tAccession_No\tTaxon_ID\tIdentity\tAlign_Length\tMismatches\tEvalue\tBitscore\tSubject_Title\tS_Start\tS_End\tS_Seq"
const int ACCESSION_NO = 2;
const int S_START = 10;

// gb|EU811445.1  ref|NC_016829.1
// return EU811445     NC_016829
string getAccession(string &str){
	int  pos = str.find_first_of('|');
	if(pos == string::npos) 
		return "";
	int start = pos + 1;
	pos = str.find_first_of('.', pos);
	return pos == string::npos ? "" : str.substr(start, pos - start);
}


string getGeneInfo(ifstream &fin, string & s_start){
	int nSeqStart = stoi(s_start);
	string sLine;
	vector<string> vecColumns;	
	while(getline(fin, sLine)){
		if(sLine[0] == '#')
			continue;
		vecColumns = split(sLine, "\t");
		if(vecColumns[ FEATURE ] != "gene")
			continue;
		int nGeneStart = stoi(vecColumns[GFF_START]);
		int nGeneEnd = stoi(vecColumns[GFF_END]);
		if(nSeqStart >= nGeneStart && nSeqStart <= nGeneEnd)
			return vecColumns[ATTRIBUTE];
	}
	return "NotFound";
}


//input is "DatabasePairedBlast.txt" in format
// "Read_ID\tGI_No\tAccession_No\tTaxon_ID\tIdentity\tAlign_Length\tMismatches\tEvalue\tBitscore\tSubject_Title\tS_Start\tS_End\tS_Seq"
void geneAnno(string input, string output){
	ifstream fin(input), finMicrobeGff;
	ofstream fout(output);
	string sLine;
	vector<string> vecColumns;	
	getline(fin, sLine);  // read off title line of DatabasePairedBlast.txt
	fout << sLine << "\tAnnotate\tGene_Info" << endl;
	string s_microbe_gff_folder = mpParameters["microbe_gff_folder"];
	while(getline(fin, sLine)){
		vecColumns = split(sLine, "\t");
		string sAccessionNo = getAccession(vecColumns[ACCESSION_NO]);
		string sAnnotated = "NO", sGene_Info = "-", sMicrobeGffFile;
		if(!sAccessionNo.empty()){
			sMicrobeGffFile = s_microbe_gff_folder + sAccessionNo + ".gff";
			finMicrobeGff.open(sMicrobeGffFile);
			if(finMicrobeGff.is_open()){
				sAnnotated = "YES";
				sGene_Info = getGeneInfo(finMicrobeGff, vecColumns[S_START]);
				finMicrobeGff.close();
			}	
		}
		fout << sLine<< "\t" << sAnnotated << "\t"<< sGene_Info << endl;
	}
}

// g++ -std=c++0x -o MicrobeGeneAnno MicrobeGeneAnno.cpp
// ../MicrobeGeneAnno ../parameters.txt
int main(int argc, char *argv[] ){
	initParameters(argv[1], mpParameters);
	geneAnno("DatabasePairedBlast.txt", "MicrobeGeneAnno_DatabasePairedBlast.txt");


	return 0;
}
