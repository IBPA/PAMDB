/*
 * PAMLib.h
 *
 *  Created on: Aug 8, 2012
 *      Author: linh
 */

#ifndef PAM_LIB_H_
#define PAM_LIB_H_

#include "utilities.h"
#include "functional_spec.h"
#include "biological_spec.h"
#include "physical_spec.h"

class PublicationReference: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string get_str();
	string pub_id;
	string title;
	StringList author;
	string journal;
	string volume;
	string number;
	string pages;
	string year;
};

class Publication: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	PublicationReference reference;
	vector<DnaComponent> dna_component_list;
	vector<Strain> strain_list;
	vector<MolecularSpecies> molecular_species_list;
	vector<Reaction> reaction_list;
	vector<Experiment> experiment_list;
};

class DnaComponentTable {
public:
	void print(ostream*);
	void add_component(DnaComponent*);
	bool get_component(string, DnaComponent&);
	string get_primary_name(string);
	void get_full_feature(string, vector<DnaFeature>&);
	double get_copy_number (vector<DnaFeature>*);
	void propagate_family_information(map<string,string>&); // Update the DNA component & reaction tables with the family member
private:
	map<string,string> alias_map;
	map<string, DnaComponent> dna_component_map;
};

class ReactionTable {
public:
	void print(ostream*);
	string find_interaction(string, string);
	string find_product(string, string);
	vector<Reaction> reaction_list;
	map<string,string> variant_to_family;
private:
	void find_reaction (Reaction &, string, string);
};

struct ModelAndExpData {
	PublicationReference paper_ref;
	string experiment_id;
	KineticModel kinetic_model;
	MeasurementData exp_dataset;
	void print(ostream* f);
	IdMap extract_parameter_name_map(int model_type);
	string write_simulation_code(int model_type, bool is_MATLAB_code = false);
	string export_MATLAB_simulation_code(IdMap, int model_type);
	string get_id();
};

class PAMLib {
public:
	PAMLib(string);
	void print(ostream*, int detail);
	void Export_MATLAB_code(string, int model_type);
	void Export_SQL_code(string sql_filename, string parameter_filename);
	static const int _PRINT_OPTION_ALL = 1;
	static const int _PRINT_OPTION_SUMMARY = 2;

private:
	// This is for reading the document database
	void import_from_a_folder(string);
	void import_from_a_file(string);
	// This is to convert form circuit design to the biological network
	void convert_experiment_to_model_and_data(Experiment, ModelAndExpData&);
	void convert_genotype_to_biological_network(Genotype, BiologicalNetwork&);
	void add_IO_node_to_biological_network(Environment, Measurement, BiologicalNetwork&);
	void add_molecular_species_interaction_to_biological_network(BiologicalNetwork&);
	// This is for parameter
	void get_parameter_value_bound(map<string,int>, StringList&, ValueList&, ValueList&);
	// This is to translate the network to the fitting problem in MATLAB code
	string get_math_model_type_name(int model_type);
	string Export_MATLAB_code_for_simultaneous_fitting(vector<ModelAndExpData>*, IdMap*, int model_type);
	string Export_MATLAB_code_for_simultaneous_fitting_LOOCV(vector<ModelAndExpData>*, IdMap*, int model_type);
	string Export_MATLAB_code_for_printing_simultaneous_fitting_result(vector<ModelAndExpData>*, IdMap* combined_parameter_name_map, int model_type);
	string Export_MATLAB_code_for_sequential_fitting(vector<ModelAndExpData>*, IdMap*, int model_type);
	string Export_MATLAB_code_for_sequential_fitting_LOOCV(vector<ModelAndExpData>*, IdMap*, int model_type);
	// Auxiliary functions
	void get_strain_information(string strain_name, StringList& modification_info);
	vector<Publication> publication_list;
	map<string,Strain> strain_table;
	DnaComponentTable dna_component_table;
	ReactionTable reaction_table;
};

#endif /* PAM_LIB_H_ */
