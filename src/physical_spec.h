/*
 * physical_spec.h
 *
 *  Created on: May 2, 2012
 *      Author: linh, UC Davis
 */

#ifndef PHYSICAL_SPEC_H_
#define PHYSICAL_SPEC_H_

#include <vector>
#include <boost/property_tree/ptree_fwd.hpp>

#include "utilities.h"

class DnaFeature: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string name;
	int start, end;
	string direction;
};

class DnaComponent: public JSONEntityInferace {		// this part
public:
	// Id
	static const int OPERATOR				= 1;
	static const int PROMOTER				= 2;
	static const int RBS					= 3;
	static const int CDS					= 4;
	static const int DNA					= 5;
	static const int TERMINATOR				= 6;
	static const int VECTOR_ORIGIN			= 7;
	static const int ANTIBIOTIC_RESISTANCE 	= 8;
	static const int PLASMID				= 9;
	// Converter
	static int Name2Id(string);
	static string Id2Name(int);

	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string name;
	vector<string> alias;
	string family;
	string type;
	string sequence;
	vector<DnaFeature> annotation;
private:
	static IdMap initialize_name_map();
	static const IdMap name_map;
};

class Strain: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string name;
	string parent;
	vector<string> gene_modification_list;
};

class Genotype: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string host_strain_name;
	vector<string> plasmid_name_list;
};

#endif /* PHYSICAL_SPEC_H_ */
