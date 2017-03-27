#include <cassert>
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

unordered_map<string, string>  mpParameters;   // load parameters in parameters.txt file 

void parseHisatSam(string input){
    //e.g.,"./CAGATC_out_hisat/CAGATC_hisat_output_-1.sam"
    ifstream fin(input);
    if(!fin)
        cout << "cannot open file" << endl;
    ofstream fOutSingleUnknown("SingleUnknown.fasta"), fOutPairedUnknown("PairedUnknown.fasta");
    string readLine;
    vector<string> vecSamColumns;
    int nSumFlags;
    string recordReadID, recordPairID, readID, pairID;
    // example of readID:
    // MG00HS08:584:C5YUBACXX:6:2316:20502:100568_2:N:0:CAGATC
    // @SL-XAS:1:100:1000:1054_1
    // example of pairID:
    // MG00HS08:584:C5YUBACXX:6:2316:20502:100568
    // @SL-XAS:1:100:1000:1054
    int nNM;
    int nMaxMismatch = stoi(mpParameters["num_max_allowed_mismatch_hisat"]);
    bool bPairedHost = false;
    unordered_map<string, vector<string> > mpSequenceForReadID; //info. for any pair, key is readID, value is vecSamColumns
    unordered_set<string> stMappedReadID; // description for each reported pair of reads
    unordered_set<string> stSingleHostReadID;  //this set could have one member, the readID mapped to host and the other unknown
                                    // or this set has two members, when (1) L-R : Mapped-Unmapped (2) L-R : Unmapped-Mapped
    bool bSecondReadFound = true; // for each read in a pair, both have been processed. Initial value true
    size_t nPairReads = 0;
    while(getline(fin, readLine)){
        if(readLine.empty())
            continue;
        if(readLine[0] == '@')
            continue;
        vecSamColumns = split(readLine , "\t");
        SamHit samhit(vecSamColumns);       
        readID = samhit.read_id;
        size_t pos = readID.find_first_of(*mpParameters["char_before_mate_num_in_fastq_title"].begin());
        pairID = readID.substr(0, pos);
        // sam data format must not have report for the same read in two lines
        // sam data may have
        // readA's align
        // readB's align
        // readA's align 
        if(readID != recordReadID){ 
            bSecondReadFound = !bSecondReadFound;
            recordReadID = readID;
        }

        if(pairID != recordPairID){
            ++nPairReads;   //for test purpose
            recordPairID = pairID;
            bSecondReadFound = false;
            if(!mpSequenceForReadID.empty()){
                assert(mpSequenceForReadID.size() == 2);
                // this read is from a new pair, output possible spliced read and clear data structure

                // not interested in paried host
                if(!bPairedHost){
                    if(stSingleHostReadID.empty()){  // paired Unknown
                        for(auto it : mpSequenceForReadID)
                            fOutPairedUnknown << ">" << it.first << endl << it.second[9] << endl;  // it.second[9] is sequence
                    }
                    else {
                        // (1) if stSingleHostReadID has only one member, this paired read is Host-Unknown or Unknown-Host
                        // (2) if stSingleHostReadID has two members, 
                        //     it means this read pair has two possible alignments: 1. Unknown-Host and 2. Host-Unknown
                        //     pick *begin() in stSingleHostReadID as host read corresponding to either 1 or 2
                        //     ignoring the other alignment case, where *end() be the readID from host
                        assert(stSingleHostReadID.size() <= 2); 
                        for(auto it : mpSequenceForReadID)
                            if(it.first != *stSingleHostReadID.begin())
                                fOutSingleUnknown << ">" << it.first << endl << it.second[9] << endl;
                    }
                }
                bPairedHost = false;
                mpSequenceForReadID.clear();
                stSingleHostReadID.clear();
            }
        }
        else if(bPairedHost == true)  // Done with this paired sequence
            continue;
        //in hisat output sam format, sometimes a pair reads has multiple possible alignment 
        // readA's align
        // readB's align
        // readA's align
        // readB's align
        // if any pair is Host-Unknown, add the read mapped to Host to stSingleHostReadID
        
        // a possible new read, save sequence information
        if(mpSequenceForReadID.find(readID) == mpSequenceForReadID.end())
            mpSequenceForReadID[readID] = vecSamColumns;
        // custom method to determine if this read is from host
        if(stMappedReadID.find(readID) == stMappedReadID.end()){ // this read is not a valid mapped host yet
            nSumFlags = stoi(samhit.sum_flags);
            if(!(nSumFlags & 1 << 2)){   // mapped host  condition 1
                nNM = getNM(readLine);
                if(nNM <= nMaxMismatch){ // ondition 2 , edit distance <= nMaxMismatch
                    stMappedReadID.insert(readID);
                    if(stMappedReadID.size() == 2){
                        bPairedHost = true;
                    }
                }
            }   
        }

        if(bSecondReadFound){ 
            // if this pair is Host-Unknown, add host readID to stSingleHostReadID
            if(stMappedReadID.size() == 1)
                stSingleHostReadID.insert(*stMappedReadID.begin());
            //start with an empty stMappedReadID for next new reported pair
            stMappedReadID.clear();
        }
        //if(++count == 2000)
        //  break;
    }
    // take care of last paired reads in file
    if(!bPairedHost && !mpSequenceForReadID.empty()){
        assert(mpSequenceForReadID.size() == 2);
        if(stSingleHostReadID.empty()){  // paired Unknown
            for(auto it : mpSequenceForReadID)
                fOutPairedUnknown << ">" << it.first << endl << it.second[9] << endl;  // it.second[9] is sequence
        }
        else {
            // (1) if stSingleHostReadID has only one member, this paired read is Host-Unknown or Unknown-Host
            // (2) if stSingleHostReadID has two members, 
            //     it means this read pair has two possible alignments: 1. Unknown-Host and 2. Host-Unknown
            //     pick *begin() in stSingleHostReadID as host read corresponding to either 1 or 2
            //     ignoring the other alignment case, where *end() be the readID from host
            assert(stSingleHostReadID.size() == 2); 
            for(auto it : mpSequenceForReadID)
                if(it.first != *stSingleHostReadID.begin())
                    fOutSingleUnknown << ">" << it.first << endl << it.second[9] << endl;
        }
    }
    cout << nPairReads << " paired reads are processed." << endl; 
    fin.close();
    fOutSingleUnknown.close();
    fOutPairedUnknown.close();
}

// g++ -std=c++0x -o ParseHisatSam ParseHisatSam.cpp
// ./ParseHisatSam ./parameters.txt
int main(int argc, char *argv[]){
    
    initParameters(argv[1], mpParameters);
    parseHisatSam(mpParameters["hisat_output_sam"]);


    return 0;
}
