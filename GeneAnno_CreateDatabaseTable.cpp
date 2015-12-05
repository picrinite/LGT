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

// annotate 
const int START_INDEX = 1;
const int END_INDEX = 2;


// NC_000001.11    BestRefSeq      gene    11874   14409   .       +       .       ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;pseudo=true
const int SEQNAME = 0;
const int SOURCE = 1;
const int FEATURE = 2;
const int START = 3;
const int END = 4;
const int ATTRIBUTE = 8;

struct Gene{
	string Genbank;
	string GeneID;
	string product;
	string attribute;
	string description;
	static string getInfo( string attr, string info ){
		int start, after_end;
		start = attr.find( info );
		string ret;
		if( start != string::npos ){
			start = attr.find_first_of("=:", start);	
			if(start != string::npos){
				start += 1;
				after_end = attr.find_first_of( ",;", start );
				if( after_end != string::npos )
					ret = attr.substr( start, after_end - start );
				else
					ret = attr.substr( start );
			}
		}
		return ret;
	}
	friend bool operator<( const Gene & left, const Gene & right ){
		return left.GeneID < right.GeneID;
	}
	Gene( string attr ){
		GeneID = getInfo( attr, "GeneID" );
		Genbank = getInfo( attr, "Genbank" );
		product = getInfo( attr, "product" );
		description = getInfo( attr, "description" );
		attribute = attr;
	}
	string getDetails() const {
		return description + "," + GeneID + "," + product + "," + Genbank + "#";
	}
};

typedef pair<string, double> PAIR;

int getStart( const string & str ){
	int sep;
	sep = str.find_first_of ( '$' );
	return stoi( str.substr( 0, sep ) );
}

int getEnd( const string & str ){
	int sep;
	sep = str.find_first_of ( '$' );
	return stoi( str.substr( sep + 1 ) );
}

class MyLess{
public:
	bool operator()( const string & strA, const string & strB){
		int firstA, firstB, secondA, secondB;	
		firstA = getStart( strA );
		secondA = getEnd( strA );
		firstB = getStart( strB );
		secondB = getEnd( strB );
		if( firstA < firstB )
			return true;
		else if ( firstA > firstB )
			return false;
		else
			return secondA < secondB;
	}
};

typedef map<string, set<Gene>, MyLess > MSVS;

class GeneCountLess{
public:
	bool operator()( const PAIR & pairA, const PAIR & pairB) {
		return pairA.second > pairB.second;
	}
};

// @param, map key sequence accession number , map value is MSVS, a map(START$END, set<Gene> )
void buildMap(const string input, unordered_map<string, MSVS> & mp){
	string readLine, key, seqname;
	vector<string> vecColumns;	
	ifstream in(input);
	while(getline(in, readLine)){
		if(readLine[0] == '#')
			continue;
		// NC_000001.11    BestRefSeq      gene    11874   14409   .       +       .       ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;pseudo=true
		vecColumns = split(readLine, "\t");		
		if(vecColumns[ FEATURE ] != "gene")
			continue;
		seqname = vecColumns[SEQNAME];
		key = vecColumns[START] + "$" + vecColumns[END];
		Gene gene(vecColumns[ATTRIBUTE]);
		set<Gene> set_genes;
		set_genes.insert(gene);
		if(mp.find(seqname) == mp.end()){
			MSVS msvs;
			mp[seqname] = msvs;
		}
		mp[seqname].insert(make_pair(key, set_genes));
	}
}

// MSVS:map(START$END, set<Gene> ) ,  there is only one Gene in set
void block(MSVS &mp){
	if(mp.size() < 2)
		return;
	MSVS::iterator p, pNext, pTmp;
	int pStart, pEnd, nextStart, nextEnd ; // two old 
	int leftEnd, midStart, midEnd, rightStart;    // three new
	string leftKey, midKey, rightKey;
	stringstream ss;
	p = mp.begin();
	pNext = p;
	++pNext; 
	while( pNext != mp.end() ) {
		pStart = getStart( p->first );
		pEnd = getEnd( p->first );
		nextStart = getStart( pNext->first );
		nextEnd = getEnd( pNext->first );
		//cout << "Enter loop: p is ( " << pStart << ",  " << pEnd << "), pNext is ( " << nextStart << ", " << nextEnd << " )." << endl;
		if( nextStart > pEnd ){
			//cout << "no overlap " << endl;  //   (3,8 ) (9, 12 )   9 > 8
			p = pNext;
			++pNext;
			continue;
		}		

		if( nextStart > pStart ){
			//cout << "insert left " << endl;			
			leftEnd = nextStart - 1;
			ss.str("");
			ss << pStart << "$" << leftEnd;
			leftKey = ss.str();
			//cout << "leftKey is : " << leftKey << endl;
			mp.insert( MSVS::value_type ( leftKey, p->second ) );
		}
		if( nextEnd > pEnd ){   //  ( 3,10) (3,7)  10 > 7
			//cout << "insert right" << endl;
			rightStart = pEnd + 1;
			ss.str("");
			ss << rightStart << "$" << nextEnd ;
			rightKey = ss.str();
			//cout << "rightKey is : " << rightKey << endl;
			pTmp = mp.find( rightKey );
			if ( pTmp != mp.end() )
				pTmp->second.insert( pNext->second.begin(), pNext->second.end() );
			else
				mp.insert( MSVS::value_type( rightKey, pNext->second ) );
			if( nextStart > pStart ){
				//cout << "insert middle" << endl;
				midStart = nextStart;
				midEnd = pEnd;
				ss.str("");
				ss << midStart << "$" << midEnd ;
				midKey = ss.str();
				//cout << "midKey is : " << midKey << endl;
				mp.insert( MSVS::value_type( midKey, p->second ) );
				mp.erase( p );
				pTmp = mp.find( midKey );
				pTmp->second.insert( pNext->second.begin(), pNext->second.end() ); 
				p = pTmp;
			}
			else if ( nextStart == pStart ) 
				p->second.insert( pNext->second.begin(), pNext->second.end() );  
			mp.erase( pNext );
			pNext = p;
			++pNext;
		}
		// when nextEnd <= pEnd
		else {				
			pNext->second.insert( p->second.begin(), p->second.end() );	  // nextStart > pStart , nextEnd <= pEnd,  pNext is middle
			if ( nextEnd < pEnd ){   //   (4,6 ) (3,8)  
				//cout << "insert right " << endl;
				rightStart = nextEnd + 1;
				ss.str("");
				ss << rightStart << "$" << pEnd ;
				rightKey = ss.str();
				//cout << "rightKey is : " << rightKey << endl;
				pTmp = mp.find( rightKey );
				if( pTmp != mp.end() )
					pTmp->second.insert( p->second.begin(),p->second.end() );
				else
					mp.insert( MSVS::value_type ( rightKey, p->second ) );
			}
			mp.erase ( p ) ;
			p = pNext;
			++pNext;
		}
	}
}

void fillArr( vector<int> &arr, MSVS &mp ){
	arr.push_back( 0 );
	for( MSVS::iterator p = mp.begin(); p != mp.end() ; ++p ){
		string key = p->first;
		arr.push_back( getStart( key ) );
		arr.push_back( getEnd( key ) );
	}
}

// find start$end that overlaps with fragment from start to end
void getKeys ( int start, int end, vector<int> &arr, set<string> &setKeys ){
	int i, start_index, end_index;
	vector<int>::iterator low, up;
	start_index = 0;
	end_index = 0;
	if( end < arr[1] || start > arr.back() )
		return;
	low = lower_bound(arr.begin(), arr.end(), start);
	start_index = low - arr.begin();
	if(start_index % 2 == 0)
		--start_index;
	up = upper_bound(arr.begin(), arr.end(), end);
	end_index = up - arr.begin();

	if(end_index % 2)
		++end_index;
	stringstream ss;
	string key;
	for( i = start_index; i < end_index; i += 2 ){
		ss.str("");
		ss << arr[i] << "$" << arr[ i + 1 ] ;
		key = ss.str();
		setKeys.insert( key );
	}
}

void getNearestKeys ( int start, int end, vector<int> &arr, set<string> &setKeys ){
	vector<int>::iterator low= lower_bound(arr.begin(), arr.end(), start);
	int start_index = low - arr.begin();
	if(start_index == arr.size())
		start_index -= 2;
	stringstream ss;
	string key;
	ss.str("");
	ss << arr[start_index] << "$" << arr[ start_index + 1 ] ;
	key = ss.str();
	setKeys.insert( key );
}

void getOutput( string input, string output, unordered_map<string, MSVS> & mp ){
	unordered_map<string, vector<int> > mp_arr;  // key is seqname(accession number), map is block indices  
	set<string> setKeys;
	ifstream in(input);
	ofstream out(output);
	int i, size;
	string genome_seqname, readLine; 
	// 11M101N40M , len_M1 = 11, len_N = 101, len_M2 = 40
	vector<string> vecColumns, vecSamColumns;

	int start, end, found, pos, pos_M1;
	
	for( auto it = mp.begin() ; it != mp.end() ; ++it ){
		vector<int> vec_empty;
		mp_arr.insert( make_pair( it->first, vec_empty ) );
		fillArr( mp_arr[ it->first ], it->second );
		//cout << "seqname is :" << it->first << endl;
		//cout << "arr size is :" << mp_arr[ it->first ].size() << endl;
	}

	string title = "Read_ID\tSum_of_Flags\tAccession_No\tChromosome_No\tStart_Position\tMap_Q_Score\tCIGAR\tOverlap_or_Nearest\tGene_Name\tGene_ID";
	out << title << endl;
	while( getline( in, readLine ) ) {
		if(readLine.empty())
			continue;
		string cigar, len_M1, len_N, len_M2;
		setKeys.clear();
		vecSamColumns = split( readLine , "\t" ) ;
		SamHit samhit(vecSamColumns);
		cigar = samhit.cigar;
		found = cigar.find_first_not_of("0123456789MN");
		if( found == string::npos ){
			pos = cigar.find_first_of ( 'M' );
			len_M1 = cigar.substr( 0, pos );			
			pos_M1 = pos;
			pos = cigar.find_first_of ( 'N', pos_M1 + 1 );
			if( pos == string::npos )
				len_N = "";
			else{
				len_N = cigar.substr( pos_M1 + 1, pos - pos_M1 );			
				if ( cigar.find_first_of( 'N', pos + 1 ) != string::npos ) { // more than one N in CIGAR
					out << mpParameters["sample_id"] << ":" << samhit.getDatabaseRow() << "\t\\\t\\\t\\" << endl;
					continue;
				}	
				else
					len_M2 = cigar.substr( pos + 1, cigar.size() - 1 - pos );
			}
		}
		else { // not annotate
			out << mpParameters["sample_id"] << ":" << samhit.getDatabaseRow() << "\t\\\t\\\t\\" << endl;
			continue;
		}		
		start = stoi( samhit.mapping_position );
		end = start + stoi( len_M1 ) - 1;
		genome_seqname = getSeqName(samhit.ref_chromosome);
		
		getKeys ( start, end, mp_arr[ genome_seqname ], setKeys );

		if ( len_N != "" ){
			start = end + stoi( len_N ) + 1;
			end = start + stoi( len_M2 ) - 1;
			getKeys ( start, end, mp_arr[ genome_seqname ], setKeys );
		}	
		string overlap_status("overlap");

		if(setKeys.empty()){
			start = stoi( samhit.mapping_position );
			getNearestKeys ( start, end, mp_arr[ genome_seqname ], setKeys );
			overlap_status = "nearest";
		}
		string key = *setKeys.begin();   // randomly select one start$end
		if(mp[genome_seqname].find(key) == mp[genome_seqname].end()){
			cout << "key not found: "  << key << endl;
			continue;
		}
		auto gene = mp[ genome_seqname ][key].begin();  // randomly select one gene in range start$end
		out << mpParameters["sample_id"] << ":" << samhit.getDatabaseRow() << "\t" << overlap_status << "\t" << gene->description << "\t" << gene->GeneID  << endl;		
	}
	in.close();
	out.close();
}

// g++ -std=c++0x -o GeneAnno_CreateDatabaseTable GeneAnno_CreateDatabaseTable.cpp
// ../GeneAnno_CreateDatabaseTable ../parameters.txt
int main(int argc, char *argv[] ){
	initParameters(argv[1], mpParameters);

	unordered_map<string, MSVS> mp;   // key is accession number, value is MSVS
	buildMap( mpParameters["host_gff_annotation_file"], mp );
	for( auto it = mp.begin(); it != mp.end(); ++it )
		block( it->second );

	//test block()
/*	for( auto it = mp.begin(); it != mp.end(); ++it ) {
		cout << it->first << endl;
		for( auto p = it->second.begin(); p!= it->second.end(); ++p ){
			cout << p->first << ", " ;
			for( auto pp = p->second.begin(); pp != p->second.end(); ++pp )
				cout << pp->GeneID << "_";
			cout << endl;
		}
	}*/

	getOutput( "PairedHostTophat.txt", "DatabasePairedHostTophat.txt", mp );

	return 0;
}
