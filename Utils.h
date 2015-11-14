#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
using namespace std;

const int SAM_QUERY_ID = 0;  // the index of query id in sam file

vector<string> split( string & str, string separator ){
	int len_separator = separator.length();
	int i,j;
	vector<string> res;
	bool startFound = false;
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


#endif