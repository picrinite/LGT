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
#include "Utils.h"
using namespace std;

// annotate and count bases numbers in assigned gene for reads from a single chrosome
// .gtf format:
// chr1	mm10_refFlat	CDS	3216025	3216968	0.000000	-	2	gene_id "Xkr4"; transcript_id "Xkr4"; 

// chr1	3216025	3216968	gene_id "Xkr4"; transcript_id "Xkr4"; 


const int START_INDEX = 1;
const int END_INDEX = 2;
const int READ_LENGTH = 101;

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
		string res = info + "=";
		start = attr.find( info );
		if( start != string::npos ){	
			after_end = attr.find_first_of( ",;", start );
			if( after_end != string::npos )
				res = attr.substr( start, after_end - start );
			else
				res = attr.substr( start );
		}
		return res;
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
	return atoi( &str.substr( 0, sep )[0] );
}

int getEnd( const string & str ){
	int sep;
	sep = str.find_first_of ( '$' );
	return atoi( &str.substr( sep + 1 )[0] );
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

// @param, map key sequence accession number , map value is MSVS:map(START$END, set<Gene> )
void buildMap( const char* input, unordered_map<string, MSVS> & mp ){
	string readLine, key, seqname;
	vector<string> vecColumns;	
	ifstream in;

	in.open( input );
	while( getline( in, readLine ) ){
		if( readLine[0] == '#' )
			continue;
		// NC_000001.11    BestRefSeq      gene    11874   14409   .       +       .       ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;pseudo=true
		vecColumns = split( readLine, "\t" );		
		if( vecColumns[ FEATURE ] != "gene" )
			continue;
		seqname = vecColumns[ SEQNAME ];
		key = vecColumns[ START ] + "$" + vecColumns[ END ];
		Gene gene( vecColumns[ATTRIBUTE] );
		set<Gene> set_genes;
		set_genes.insert( gene );
		if( mp.find( seqname ) != mp.end() )
			mp[ seqname ].insert( make_pair( key, set_genes ) );
		else {
			MSVS msvs;
			mp[ seqname ] = msvs;
		}
	}
}

// MSVS:map(START$END, set<Gene> ) ,  there is only one Gene in set
void block ( MSVS &mp ){
	MSVS::iterator p, pNext, pTmp;
	int pStart, pEnd, nextStart, nextEnd ; // two old 
	int leftEnd, midStart, midEnd, rightStart;    // three new
	string leftKey, midKey, rightKey;
	stringstream ss;
	if( mp.size() < 2 )
		return;
	p = mp.begin();
	pNext = p;
	++pNext ; 
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
		string key;
		key = p->first;
		arr.push_back( getStart( key ) );
		arr.push_back( getEnd( key ) );
	}
}
// return NC_000008.11 from gi|568815590|ref|NC_000008.11|
string getSeqname( string & str ){
	int start = str.find( "|N") + 1 ;
	int after_end = str.find_first_of( "|", start );
	return str.substr( start, after_end - start );
}


map<string, double>  mp_gene_stat; 
map<string, string>  mp_gene_fullname; 
set<string> setReadID; // a temporary set recording read title

void getKeys ( int start, int end, vector<int> &arr, set<string> &setKeys ){
	int i, arr_len, imin, imax, imid, start_index, end_index;
	arr_len = arr.size();
	imin = 1; // the first element in vector<int> arr is a useless zero
	imax = arr_len - 1;
	start_index = 0;
	end_index = 0;
	if ( end < arr[ imin ] || start > arr[ imax ] )
		return;
	if( start <= arr[ imin ] )   
		start_index = imin;
	else if ( start == arr[ imax ] ){
		start_index = imax - 1;
		end_index = imax;
		return;
	}
	else {           // arr[imin] < start < arr[imax] 
		while( imax - imin != 1 ){	
			imid = ( imax + imin ) / 2;
			if( start == arr[ imid ] ){
				if ( imid % 2 )  // imid is odd number
					start_index = imid;
				else 
					start_index = imid - 1;
				break;
			}
			else if ( start > arr[ imid ] )
				imin = imid ;
			else
				imax =  imid;
		}
		if( start_index == 0 ){   // if not found in while loop, then imax - imin = 1
			if ( imin % 2 )
				start_index = imin;
			else
				start_index = imax;
		}
	}
	//cout << "start_index is: " << start_index << endl;
	imin = 1;
	imax = arr_len - 1;
	if( end >= arr[ imax ] )
		end_index = imax;
	else if ( end == arr[ imin ] )
		end_index = imin + 1;
	else {      //  arr[imin] < end_index < arr[imax]
		while( imax - imin != 1 ){	
			imid = ( imax + imin ) / 2;
			if( end == arr[ imid ] ){
				if ( imid % 2 )  // imid is odd number, left of unit
					end_index = imid + 1;
				else
					end_index = imid;
				break;
			}
			else if ( end > arr[ imid ] )
				imin = imid ;
			else
				imax =  imid;
		}
		if( end_index == 0 ){   // imax - imin = 1,  arr[imin] < end < arr[imax]
			if ( imin % 2 )  // imid is odd number, left of unit, imax is right of unit
				end_index = imax;
			else
				end_index = imin;
		}
	}
	//cout << "end_index is: " << end_index << endl;
	stringstream ss;
	string key;
	for( i = start_index; end_index - i >= 1; i += 2 ){
		ss.str("");
		ss << arr[i] << "$" << arr[ i + 1 ] ;
		key = ss.str();
		setKeys.insert( key );
	}
}

int gene_count = 0;
void getGeneCount( int start, int end, set<string> &setKeys, MSVS &mp, map<string, double> &mp_gene_count , set<string> &set_gene_details, string readID ){
	int unitLeft, unitRight, overlapStart, overlapEnd;
	set<Gene> mpValue;
	set<string> setGeneNames;
	map<string, double>::iterator tmpP;
	string geneName;
	for( auto p = setKeys.begin(); p != setKeys.end(); ++p ){
		//setGeneNames.clear();
		unitLeft = getStart( *p );
		unitRight = getEnd( *p );
		overlapStart = max( unitLeft, start );
		overlapEnd = min( unitRight, end );
		double basesCount = overlapEnd - overlapStart + 1;
		mpValue = mp[ *p ];  // set of Gene
		for( auto pp = mpValue.begin(); pp != mpValue.end(); ++pp ){
			set_gene_details.insert( pp->getDetails() );   // add to set of gene_transcript
			geneName = pp->GeneID;  
			// mp_gene_fullname
			if( mp_gene_fullname.find( geneName ) != mp_gene_fullname.end() )
				++gene_count;
			mp_gene_fullname[ geneName ] = pp->attribute;

			// mp_gene_count
			tmpP = mp_gene_count.find( geneName );
			if( tmpP == mp_gene_count.end() )
				mp_gene_count.insert( make_pair( geneName, basesCount / READ_LENGTH ) );
			else
				tmpP->second += basesCount / READ_LENGTH ;
		}
		if( mpValue.size() == 1 ){    // a unique gene for this start-end range
			auto pp = mpValue.begin();
			geneName = pp->GeneID;  
			tmpP = mp_gene_stat.find( geneName );
			if( tmpP == mp_gene_stat.end() )
				mp_gene_stat.insert( make_pair( geneName, basesCount / READ_LENGTH ) );
			else
				tmpP->second += basesCount / READ_LENGTH ;

			 // Add this sequence to set: setReadID,  which is a temporary set recording read title
            setReadID.insert( readID );
		}		
	}
}  
void getOutput( char* input, char* output, unordered_map<string, MSVS> & mp ){
	unordered_map<string, vector<int> > mp_arr;  // key is seqname(accession number), map is block indices
    
    
	set<string> setKeys;
	set<string> set_gene_details;
	map<string, double > mp_gene_count;  // key is geneID + transciptID , value is count

	ifstream in;
	ofstream out;
	int i, index_num_mapping, num_mapping, index_mapping_list, index_host_read ;
	string titles, cigar, len_M1, len_N, len_M2, genome_seqname, readLine, gene_count_output, gene_details_output; 
	// 11M101N40M , len_M1 = 11, len_N = 101, len_M2 = 40
	vector<string> vecColumns, vecSamColumns;

	int start, end, found, pos, pos_M1, index_read_length, read_length;
	
	stringstream ss;

	in.open( input );
	out.open( output );

	for( auto it = mp.begin() ; it != mp.end() ; ++it ){
		vector<int> vec_empty;
		mp_arr.insert( make_pair( it->first, vec_empty ) );
		fillArr( mp_arr[ it->first ], it->second );
		//cout << "seqname is :" << it->first << endl;
		//cout << "arr size is :" << mp_arr[ it->first ].size() << endl;
	}
	// read titles
	getline( in, titles );
	
	vecColumns = split( titles, "\t");
	index_host_read = getIndexOfColumn( vecColumns, "host_read" );  // MG00HS08:584:C5YUBACXX:6:1309:7776:84081_1:N:0:CAGATC
	index_num_mapping = getIndexOfColumn( vecColumns, "#Tophat_Mapping" );
	index_read_length = getIndexOfColumn( vecColumns, "host_read_length" );
	index_mapping_list = index_num_mapping + 1;
    
	out << titles << "\t#Gene_match\t#host_gene\thost_genes_details" << endl;
	int count = 0;
	while( getline( in, readLine ) ) {
		out << readLine ;
		vecColumns = split( readLine, "\t" );
		num_mapping = atoi ( vecColumns[ index_num_mapping ].c_str() );
		if( num_mapping != 1 ){ // not unique and thus annotation not needed
			out << "\t\\\t\\\t\\" << endl;
			continue;
		}
	
		setKeys.clear();
		mp_gene_count.clear();
		set_gene_details.clear();
		gene_count_output = "\\";
		gene_details_output = "\\"; 

		//  vecColumns[ index_mapping_list ] :      73,gi|568815578|ref|NC_000020.11|,8174332,50,101M#
		vecSamColumns = split( vecColumns[ index_mapping_list ] , "," ) ;
		int size = vecSamColumns.size();
		cigar = vecSamColumns[ size - 1 ];
		size = cigar.length();
		cigar = cigar.substr( 0, size - 1 );  // remove # from 101M
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
						out << "\t\\\t\\\t\\" << endl;
						continue;
					}	
					else
						len_M2 = cigar.substr( pos + 1, cigar.size() - 1 - pos );
				}
			}
		else { // not annotate
			out << "\t\\\t\\\t\\" << endl;
			continue;
		}		
		start = atoi( &vecSamColumns[2][0] );
		end = start + atoi( &len_M1[0] ) - 1;
		genome_seqname = getSeqname( vecSamColumns[ 1 ] );
		read_length = atoi( vecColumns[ index_read_length ].c_str() ) ;
		
		getKeys ( start, end, mp_arr[ genome_seqname ], setKeys );
		getGeneCount( start, end, setKeys, mp[ genome_seqname ], mp_gene_count, set_gene_details, vecColumns[index_host_read] );

		if ( len_N != "" ){
			start = end + atoi( &len_N[0] ) + 1;
			end = start + atoi( &len_M2[0] ) -1;
			setKeys.clear();
			getKeys ( start, end, mp_arr[ genome_seqname ], setKeys );
			getGeneCount( start, end, setKeys, mp[ genome_seqname ], mp_gene_count, set_gene_details, vecColumns[index_host_read] );
		}		
		if( mp_gene_count.size() ){
			vector<PAIR> vec_gene_count( mp_gene_count.begin(), mp_gene_count.end() );
			sort( vec_gene_count.begin(), vec_gene_count.end(), GeneCountLess() );
					//print gene names and count
			ss.str("");
			for( vector<PAIR>::iterator p = vec_gene_count.begin(); p != vec_gene_count.end(); ++p )
				ss << p->first << "(" << p->second << ")$" ;
			gene_count_output = ss.str();
			ss.str("");
			for( set<string>::iterator p = set_gene_details.begin(); p != set_gene_details.end(); ++p )
				ss << *p ;
			gene_details_output = ss.str();
		}
		out << "\t" << mp_gene_count.size() << "\t" + gene_count_output + "\t" + gene_details_output << endl;
	}
	in.close();
	out.close();
}

void getStatOutput( const char* outputFile ){
	ofstream out;
	out.open( outputFile );
	out << "Count\tGene_ID\tGene_description\tFull_list" << endl;
	cout << "Total number of genes: " << mp_gene_stat.size() << endl;
	out << "Total number of genes: " << mp_gene_stat.size() << endl;
	vector<PAIR> vec_gene_count( mp_gene_stat.begin(), mp_gene_stat.end() );
	sort( vec_gene_count.begin(), vec_gene_count.end(), GeneCountLess() );

	int countWithoutDescrip = 0;
	//print gene names and count
	for( vector<PAIR>::iterator p = vec_gene_count.begin(); p != vec_gene_count.end(); ++p ) {
		string geneName = p->first;
		string attribute = mp_gene_fullname[ geneName ];
		string description = Gene::getInfo( attribute, "description" );
		if( description == "description=" )
			++countWithoutDescrip;

		out << p->second << "\t" << geneName << "\t" << description << "\t" << attribute << endl ;
	}
	out << "The number of genes without description is : " << countWithoutDescrip << endl;
	cout << "The number of genes without description is : " << countWithoutDescrip << endl;
}
//argv[1] .gff annotation file , argv[2] inputfile with host read sequence
 // g++ -std=c++0x -o GeneAnno_human GeneAnno_human_0608_gene.cpp
// ../GeneAnno_human ../ref_GRCh38.p2_top_level.gff3 PHostUMicrobe.add_microbe.P1P2P3.txt
int main(int argc, char *argv[] ){
	
	//string gffFile("../ref_GRCh38.p2_top_level.gff3");
	string gffFile( argv[1] );
	string inputFile( argv[2] );
	//string barcode( argv[3] );
	string outputFile ( "PHostUMicrobe.host_gene.P1P2P3P4.txt" );
	string outputFile2 ( "PHostUMicrobe.host_gene_summary.P1P2P3P4.txt" );

	unordered_map<string, MSVS> mp;   // key is accession number, value is MSVS
	buildMap( gffFile.c_str(), mp );
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

	getOutput( &inputFile[0], &outputFile[0], mp );
	getStatOutput( outputFile2.c_str() );
	cout << "Check repeat by gene_count: " << gene_count << endl;
	

	// output host reads that can be align to genes
    ofstream outReadID;
    outReadID.open( "outReadIDFromSingleHostGene.txt" );
    for (auto p : setReadID){
    	outReadID << p << endl;
    }
    outReadID.close();
    // Extract sam lines for readId in setReadID
    ifstream inSam;
    ofstream outHostGeneSam;
    inSam.open( "CAGATC_SingleHost.sam" );
    outHostGeneSam.open( "outHostGene.sam" );
    string readLine;
    vector<string> vecColumns;
    while( getline(inSam, readLine) ){
        vecColumns = split( readLine, "\t" );
        if ( setReadID.find( vecColumns[SAM_QUERY_ID] ) != setReadID.end() )
        	outHostGeneSam << readLine << endl;
    }
    inSam.close();
    outHostGeneSam.close();

	/*
	// begin test 
	set<string> first, second, third, fourth, vec_str;
	first.insert("geneA");
	second.insert("geneB");
	third.insert("geneC");
	fourth.insert("geneD");
	vec_str.insert("gene");
	mp.insert( MSVS::value_type ( "3$7", first ) );
	mp.insert( MSVS::value_type ( "3$10", second ) );
	mp.insert( MSVS::value_type ( "8$10", third ) );
	mp.insert( MSVS::value_type ( "10$15", fourth ) );
	mp.insert( MSVS::value_type ( "13$18", vec_str ) );
	mp.insert( MSVS::value_type ( "20$30", vec_str ) );
	mp.insert( MSVS::value_type ( "24$40", vec_str ) );
	mp.insert( MSVS::value_type ( "28$35", vec_str ) );
	mp.insert( MSVS::value_type ( "50$75", vec_str ) );
	mp.insert( MSVS::value_type ( "50$100", vec_str ) );
	mp.insert( MSVS::value_type ( "75$111", vec_str ) );
	
	set<string> setKeys;
	map<string, double> mp_gene_count;  // key is geneID + transciptID , value is count
	set<string> set_gene_transcript;
	vector<int> arr;

	fillArr( arr, mp );

	int start, end;
	start = 1;
	end = 7;
	getKeys ( start, end, arr, setKeys );
	cout << endl << " set of keys" << endl;
	for( set<string>::iterator p = setKeys.begin(); p != setKeys.end(); ++p )
		cout << *p << endl;
	getGeneCount( start, end, setKeys, mp, set_gene_transcript, mp_gene_count );

	start = 8;
	end = 25;
	setKeys.clear();
	getKeys ( start, end, arr, setKeys );
	cout << endl << " set of keys" << endl;
	for( set<string>::iterator p = setKeys.begin(); p != setKeys.end(); ++p )
		cout << *p << endl;
	getGeneCount( start, end, setKeys, mp, set_gene_transcript, mp_gene_count );

	vector<PAIR> vec_gene_count( mp_gene_count.begin(), mp_gene_count.end() );
	sort( vec_gene_count.begin(), vec_gene_count.end(), GeneCountLess() );
	//print gene names and count
	cout << endl << "Gene names and bases_counts after sorting:" << endl;
	for( vector<PAIR>::iterator p = vec_gene_count.begin(); p != vec_gene_count.end(); ++p )
		cout << p->first << "_" << p->second << ";" ;

	//print gene_transcript ID
	cout << endl << "Set of gene and transcript:" << endl;
	for( set<string>::iterator p = set_gene_transcript.begin(); p != set_gene_transcript.end(); ++p )
		cout << *p << ";" ;
	*/
	// end test


	return 0;
}
