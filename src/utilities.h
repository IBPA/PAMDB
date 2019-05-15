/*
 * utilities.h
 *
 *  Created on: Oct 27, 2014
 *      Author: linh
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <map>
#include <boost/property_tree/ptree.hpp>

using namespace std;

#define EPSILON 1e-8

enum DetailPrinteLevel {
	SUMMARY,
	ALL
};

enum CommonValue {
	STR_EQ	= 0,
	UNKNOWN	= -1
};

// Pair
typedef std::pair<int,int> IdPair;
typedef std::pair<double,double> ValuePair;
typedef std::pair<string,string> StringPair;
// List
typedef std::vector<int> IdList;
typedef std::vector<double> ValueList;
typedef std::vector<string> StringList;
// Pair List
typedef std::vector<IdPair> IdPairList;
typedef std::vector<ValuePair> ValuePairList;
// Matrix
typedef std::vector<IdList> IdMatrix;
typedef std::vector<ValueList> ValueMatrix;
typedef std::vector<StringList> StringMatrix;
// Map
typedef std::map<int,int> IdConverter;
typedef std::map<string,string> StringConverter;
typedef std::map<string,int> IdMap;

// Interface for JSON Entity
class JSONEntityInferace {
public:
	virtual ~JSONEntityInferace() {
		// Do nothing
	}
	virtual void print(ostream*) = 0;
	virtual void read(boost::property_tree::ptree*) = 0;
};
void sbrome_print(ostream*, IdList);
void sbrome_print(ostream*, StringList);
void sbrome_print(ostream*, IdMap);
void sbrome_error_print(string);
// For string
string correct_to_a_var_name(string);
void add_tab_to_each_sentence(string&);
void split_string(string str, char c, StringList &word_list);
// Converters
string Int2Str(int);
string Float2Str(double);
bool delete_first_pos(IdList&, int);
void frequency_counting(string, IdMap&);
void summarize(const ValueList& value_list, double &mean, double &std_dev, double &std_error);
// For files
void combine_all_json_file(string, string);
void write_all_json_file(string, ostream*);

#endif /* UTILITIES_H_ */
