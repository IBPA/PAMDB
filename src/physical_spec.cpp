/*
 * physical_spec.cpp
 *
 *  Created on: May 2, 2012
 *      Author: linh, UC Davis
 */

#include <boost/foreach.hpp>
#include "physical_spec.h"

void DnaFeature::print(ostream* f) {
	*f << "Feature: " << name << "\t";
	*f << "Start: " << start << "\t";
	*f << "End: " << end << "\t";
	*f << "Direction: " << direction << endl;
}

void DnaFeature::read(boost::property_tree::ptree* pt) {
	name = pt->get<string>("feature");
	start = pt->get<int>("start", UNKNOWN);
	end = pt->get<int>("end", UNKNOWN);
	direction = pt->get<string>("direction", "forward");
}

void DnaComponent::print(ostream* f) {
	*f << "name: " << name << "\n";
	*f << "alias: ";
	sbrome_print(f, alias);
	*f << "family: " << type << "\t";
	*f << "type: " << type << "\t";
	*f << "sequence: " << sequence << endl;
	*f << "Annotation: " << endl;
	for (unsigned int i = 0; i < annotation.size(); i++) {
		annotation[i].print(f);
	}
}

void DnaComponent::read(boost::property_tree::ptree* pt) {
	name = pt->get<string>("name");
	family = pt->get<string>("family", "");
	type = pt->get<string>("type", "DNA"); // DNA is the general DNA
	sequence = pt->get<string>("sequence","");
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("alias", empty_ptree)) {
		alias.push_back(v.second.get_value<string>());
	}
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("annotation", empty_ptree)) {
		DnaFeature dna_feature;
		dna_feature.read(&v.second);
		annotation.push_back(dna_feature);
	}
}

IdMap DnaComponent::initialize_name_map() {
	IdMap name_map_tmp;
	name_map_tmp ["DNA"] 				= DNA;
	name_map_tmp ["operator"] 			= OPERATOR;
	name_map_tmp ["promoter"] 			= PROMOTER;
	name_map_tmp ["RBS"]				= RBS;
	name_map_tmp ["CDS"] 				= CDS;
	name_map_tmp ["terminator"] 		= TERMINATOR;
	name_map_tmp ["origin"]					= VECTOR_ORIGIN;
	name_map_tmp ["Antibiotic Resistance"] 	= ANTIBIOTIC_RESISTANCE;
	return name_map_tmp;
}

const IdMap DnaComponent::name_map = DnaComponent::initialize_name_map();

int DnaComponent::Name2Id(string name) {
	return name_map.at(name);
}

string DnaComponent::Id2Name(int id) {
	for (IdMap::const_iterator it = name_map.begin(); it != name_map.end(); it++)
		if (it->second == id)
			return it->first;
	sbrome_error_print("There is no DNA component type Id = " + Int2Str(id));
	return "";
}

void Strain::print(ostream* f) {
	*f << "name: " << name << "\t";
	*f << "parent: " << parent << endl;
	*f << "gene_modification: ";
	sbrome_print(f, gene_modification_list);
}

void Strain::read(boost::property_tree::ptree* pt) {
	name = pt->get<string>("name");
	parent = pt->get<string>("parent","");
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("modification")) {
		gene_modification_list.push_back(v.second.get_value<string>());
	}
}

void Genotype::print(ostream* f) {
	*f << "Genotype: " << endl;
	*f << "Host strain: " << host_strain_name << endl;
	*f << "Plasmids: ";
	sbrome_print(f, plasmid_name_list);
}

void Genotype::read(boost::property_tree::ptree* pt) {
	host_strain_name = pt->get<string>("host");	// query to fill other information later
	plasmid_name_list.clear();
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("plasmid", empty_ptree)) {
		plasmid_name_list.push_back(v.second.get_value<string>());
	}
}

