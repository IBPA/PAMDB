/*
 * PAMLib.cpp
 *
 *  Created on: Aug 8, 2012
 *      Author: linh
 */

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include "PAMLib.h"

void PublicationReference::print(ostream* f) {
	*f << pub_id << "\t" << author[0] << "\t" << title << "\t" << journal << "\t" << volume << "\t" << number << "\t" << pages << "\t" << year << endl;
	//*f << author << "\t" << journal << "\t" << year << endl;
}

void PublicationReference::read(boost::property_tree::ptree* pt) {
	boost::property_tree::ptree empty_ptree;
	boost::property_tree::ptree tmp = pt->get_child("reference", empty_ptree);
	if (!tmp.empty()) {
		pub_id = "";
		title = tmp.get<string>("title");
		author.push_back(tmp.get<string>("author", ""));
		journal = tmp.get<string>("journal", "");
		volume = tmp.get<string>("volume", "");
		number = tmp.get<string>("number", "");
		pages = tmp.get<string>("pages", "");
		year = tmp.get<string>("year");
	}
}

string PublicationReference::get_str() {
	return title + ". " + author[0];
}

void Publication::print(ostream* f) {
	reference.print(f);
	*f << "Dna components: " << endl;
	for (unsigned int i = 0; i < dna_component_list.size(); i++)
		dna_component_list[i].print(f);
	*f << "Strains: " << endl;
	for (unsigned int i = 0; i < strain_list.size(); i++)
		strain_list[i].print(f);
	*f << "Molecular species: " << endl;
	for (unsigned int i = 0; i < molecular_species_list.size(); i++)
		molecular_species_list[i].print(f);
	*f << "Reactions: " << endl;
	for (unsigned int i = 0; i < reaction_list.size(); i++)
		reaction_list[i].print(f);
	*f << "Experiments: " << endl;
	for (unsigned int i = 0; i < experiment_list.size(); i++) {
		*f << "--------------" << endl;
		experiment_list[i].print(f);
	}
}

void Publication::read(boost::property_tree::ptree* pt) {
	reference.read(pt);
	reference.pub_id = pt->get<string>("uri");
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("DnaComponent", empty_ptree)) {
		DnaComponent dna_component;
		dna_component.read(&v.second);
		dna_component_list.push_back(dna_component);
	}
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("Strain", empty_ptree)) {
		Strain s;
		s.read(&v.second);
		strain_list.push_back(s);
	}
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("MolecularSpecies", empty_ptree)) {
		MolecularSpecies m;
		m.read(&v.second);
		molecular_species_list.push_back(m);
	}
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("Reaction", empty_ptree)) {
		Reaction r;
		r.read(&v.second);
		reaction_list.push_back(r);
	}
	map<string,double> reporter_name_to_max_value;
	IdMap exp_ref_check_list;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("Experiment", empty_ptree)) {
		Experiment e;
		e.read(&v.second);
		if (exp_ref_check_list.find(e.reference) != exp_ref_check_list.end())
			sbrome_error_print("Experiment id " + e.reference + " is repeated: " + reference.get_str() + " !\n");
		else
			exp_ref_check_list[e.reference] = 1;
		for(unsigned int reporter_id = 0; reporter_id < e.measurement.molecular_species_list.size(); reporter_id++) {
			string unit = e.measurement.molecular_species_list[reporter_id].second;
			if (unit.compare("au") == STR_EQ || unit.compare("RPU") == STR_EQ || unit.compare("REU") == STR_EQ) {
				// TODO: this is only correct for the case with one output only
				double max_output = UNKNOWN;
				for (unsigned int i = 0; i < e.dataset.output_measurement.size(); i++)
					for (unsigned int j = 0; j < e.dataset.output_measurement[i].size(); j++)
						if (max_output < e.dataset.output_measurement[i][j].mean)
							max_output = e.dataset.output_measurement[i][j].mean;
				if (reporter_name_to_max_value.find(e.measurement.molecular_species_list[reporter_id].first) == reporter_name_to_max_value.end() || reporter_name_to_max_value[e.measurement.molecular_species_list[reporter_id].first] < max_output)
					reporter_name_to_max_value[e.measurement.molecular_species_list[reporter_id].first] = max_output;
			}
			else
				sbrome_error_print("Measurement unit is not defined: " + unit + reference.get_str() + "!\n");
		}
		experiment_list.push_back(e);
	}
	// Normalize the output
	for (unsigned int exp_id = 0; exp_id < experiment_list.size(); exp_id++) {
		// TODO: this is only correct for the case with one output only
		if (reporter_name_to_max_value.find(experiment_list[exp_id].measurement.molecular_species_list[0].first) != reporter_name_to_max_value.end()) {
			for (unsigned int i = 0; i < experiment_list[exp_id].dataset.output_measurement.size(); i++)
				for (unsigned int j = 0; j < experiment_list[exp_id].dataset.output_measurement[i].size(); j++)
					experiment_list[exp_id].dataset.output_measurement[i][j].mean /= reporter_name_to_max_value[experiment_list[exp_id].measurement.molecular_species_list[0].first];
		}
	}
}

void ReactionTable::print(ostream* f) {
	for (vector<Reaction>::iterator it = reaction_list.begin(); it != reaction_list.end(); it++) {
		it->print(f);
		*f << endl;
	}
	*f << "Variant to family" << endl;
	for (map<string,string>::iterator it = variant_to_family.begin(); it != variant_to_family.end(); it++)
		*f << it->first << "\t" << it->second << endl;
}

void ReactionTable::find_reaction (Reaction &r, string ms1, string ms2) {
	string family_1 = (variant_to_family.find(ms1) != variant_to_family.end())? variant_to_family[ms1]:ms1;
	string family_2 = (variant_to_family.find(ms2) != variant_to_family.end())? variant_to_family[ms2]:ms2;
	for (vector<Reaction>::iterator reaction = reaction_list.begin(); reaction != reaction_list.end(); reaction++) {
		if (reaction->reactant_list.size() == 2 && ((reaction->reactant_list[0].compare(family_1) == STR_EQ && reaction->reactant_list[1].compare(family_2) == STR_EQ)
			|| (reaction->reactant_list[0].compare(family_2) == STR_EQ && reaction->reactant_list[1].compare(family_1) == STR_EQ))) {
			r = *reaction;
			break;
		}
	}
}

string ReactionTable::find_interaction(string ms1, string ms2) {
	Reaction r;
	r.mechanism_type.clear();
	find_reaction(r, ms1, ms2);
	return r.mechanism_type;
}

string ReactionTable::find_product(string ms1, string ms2) {
	Reaction r;
	r.mechanism_type.clear();
	find_reaction(r, ms1, ms2);
	if (r.mechanism_type.empty())
		sbrome_error_print("Can not find a product for an wrong reaction between " + ms1 + " and " + ms2 + "\n");
	return r.product_list[0];
}

void DnaComponentTable::print(ostream* f) {
	for (map<string,DnaComponent>::iterator it = dna_component_map.begin(); it != dna_component_map.end(); it++) {
		*f << it->first << "\t" << it->second.type << endl;
		*f << "Alias: ";
		sbrome_print(f, it->second.alias);
	}
}

void DnaComponentTable::add_component(DnaComponent* dna_comp) {
	if (alias_map.find(dna_comp->name) != alias_map.end())
		sbrome_error_print("DNA component " + dna_comp->name + " is re-declared!");
	else {
		dna_component_map[dna_comp->name] = *dna_comp;
		alias_map[dna_comp->name] = dna_comp->name;
	}
	for (StringList::iterator it = dna_comp->alias.begin(); it != dna_comp->alias.end(); it++) {
		if (alias_map.find(*it) != alias_map.end())
			sbrome_error_print("Alias " + *it + " is repeated!");
		else
		alias_map[*it] = dna_comp->name;
	}
}

bool DnaComponentTable::get_component(string component_name, DnaComponent& dna_comp) {
	if (alias_map.find(component_name) == alias_map.end() || dna_component_map.find(alias_map[component_name]) == dna_component_map.end())
		return false;
	else
		dna_comp = dna_component_map[alias_map[component_name]];
	return true;
}

string DnaComponentTable::get_primary_name(string name) {
	if (alias_map.find(name) == alias_map.end())
		sbrome_error_print("DNA component " + name + " not found!");
	return alias_map[name];
}

double DnaComponentTable::get_copy_number(vector<DnaFeature>* dna_feature_list) {
	double copy_number = UNKNOWN;
	for(vector<DnaFeature>::iterator it = dna_feature_list->begin(); it != dna_feature_list->end(); it++) {
		DnaComponent tmp;
		get_component(it->name, tmp);
		if (tmp.type.compare("origin") == STR_EQ) {
			if (copy_number >= 0)
				sbrome_error_print("There are more than one orgin in a DNA sequence!\n");
			if (it->name.compare("pMB1") == STR_EQ || it->name.compare("pUC8") == STR_EQ)
				copy_number = 200;
			else if (it->name.compare("BBR1") == STR_EQ  || it->name.compare("pBR322") == STR_EQ)
				copy_number = 20;
			else if (it->name.compare("ColE1") == STR_EQ || it->name.compare("ColE1-mut") == STR_EQ)
				copy_number = 15;
			else if (it->name.compare("p15A") == STR_EQ)
				copy_number = 10;
			else if (it->name.compare("pSC101") == STR_EQ)
				copy_number = 5;
			else if (it->name.compare("Chromosome") == STR_EQ || it->name.compare("BAC") == STR_EQ || it->name.compare("pPHI80") == STR_EQ)
				copy_number = 1;
			else
				sbrome_error_print("Copy number of the origin " + it->name + " not found from the code!\n");
		}
	}
	if (copy_number < 0)
		sbrome_error_print("Origin not found from the plasmid DNA sequence \n");
	return copy_number;
}

void DnaComponentTable::get_full_feature(string dna_component_name, vector<DnaFeature>& dna_feature_list) {
	DnaComponent tmp;
	if (!get_component(dna_component_name,tmp))
		sbrome_error_print("Can not find the DNA component " + dna_component_name);
	else {
		for (unsigned int i = 0; i < tmp.annotation.size(); i++) {
			vector<DnaFeature> dna_list_tmp;
			get_full_feature(tmp.annotation[i].name, dna_list_tmp);
			if (dna_list_tmp.empty()) {
				tmp.annotation[i].name = get_primary_name(tmp.annotation[i].name);
				dna_feature_list.push_back(tmp.annotation[i]);
			}
			else
				for (unsigned int j = 0; j < dna_list_tmp.size(); j++)
					dna_feature_list.push_back(dna_list_tmp[j]);
		}
	}
}

void DnaComponentTable::propagate_family_information(map<string,string>& variant_to_family) {
	variant_to_family.clear();
	for (map<string,DnaComponent>::iterator it = dna_component_map.begin(); it != dna_component_map.end(); it++)
		if (!it->second.family.empty()) {
			variant_to_family[it->first] = it->second.family;
			if (dna_component_map.find(it->second.family) == dna_component_map.end())
				sbrome_error_print("Can not find the family " + it->second.family + " of the DNA component " + it->first);
			else
				it->second.type = dna_component_map[it->second.family].type;
		}
		else
			variant_to_family[it->first] = it->first;
}

void ModelAndExpData::print(ostream* f) {
	kinetic_model.print(f);
	*f << "==========" << endl;
	exp_dataset.print(f);
}

IdMap ModelAndExpData::extract_parameter_name_map(int model_type) {
	return kinetic_model.ExtractParameterMap(paper_ref.pub_id, model_type);
}

string add_parameter_marker(bool is_MATLAB_code, string parameter_name) {
	return is_MATLAB_code? parameter_name : "{" + parameter_name + "}";
}

string ModelAndExpData::write_simulation_code(int model_type, bool is_MATLAB_code) {
	IdList topo_order_node_id_list;
	kinetic_model.Topo_Sort(&topo_order_node_id_list);
	string simulation_code = "";
	for (unsigned int topo_node_id = 0; topo_node_id < topo_order_node_id_list.size(); topo_node_id++) {			// compute the concentration of each node
		simulation_code += "\t";
		int node_id = topo_order_node_id_list[topo_node_id];
		int node_type = kinetic_model.GetNodeInfo(node_id)->type;
		int factor_num = kinetic_model.InEdgeCount(node_id);
		string node_name = kinetic_model.GetNodeInfo(node_id)->get_var_name();
		double copy_number = kinetic_model.GetNodeInfo(node_id)->copy_number;
		switch (node_type) {
			case MolecularSpecies::mRNA: {
				string promoter_name = node_name;
				string alpha = add_parameter_marker(is_MATLAB_code, "alpha_" + promoter_name);
				if (factor_num == 0) {
					simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*" + alpha + ";\n";
				}
				else if (factor_num == 1) {
					bool is_hill_eq = false;
					if (model_type == KineticModel::_MODEL_OPTION_ELLIS) {
						if (promoter_name.compare("pLAC") == STR_EQ || promoter_name.compare("pLlacO_1") == STR_EQ || promoter_name.compare("placUV5") == STR_EQ || promoter_name.compare("pTAC") == STR_EQ) {
							int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
							string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
							string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
							string KE1 = add_parameter_marker(is_MATLAB_code, "KE1_" + promoter_name);
							string KE2 = add_parameter_marker(is_MATLAB_code, "KE2_" + promoter_name);
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*(" + beta + " + (" + alpha + " - " + beta + ")/(1 + " + tf_name + "^4/" + KE1 + " + (" + tf_name + "_free^4 - " + tf_name + "^4)/" + KE2 + "));\n";
						}
						else if (promoter_name.compare("pTET") == STR_EQ || promoter_name.compare("pLtetO_1") == STR_EQ || promoter_name.compare("pTET_star_") == STR_EQ) {
							int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
							string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
							string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
							string KE1 = add_parameter_marker(is_MATLAB_code, "KE1_" + promoter_name);
							string KE2 = add_parameter_marker(is_MATLAB_code, "KE2_" + promoter_name);
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*(" + beta + " + (" + alpha + " - " + beta + ")/(1 + " + tf_name + "^4/" + KE1 + " + " + tf_name + "^2/" + KE2 + "));\n";
						}
						else if (promoter_name.compare("pBAD") == STR_EQ) {
							is_hill_eq = true;
						}
						else {
							sbrome_error_print("Ellis model does not include " + promoter_name);
						}
					}
					if (model_type == KineticModel::_MODEL_OPTION_MOON) {
						if (promoter_name.compare("pLAC") == STR_EQ || promoter_name.compare("pLlacO_1") == STR_EQ || promoter_name.compare("placUV5") == STR_EQ || promoter_name.compare("pTAC") == STR_EQ) {
							int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
							string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
							string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
							string KM1 = add_parameter_marker(is_MATLAB_code, "KM1_" + promoter_name);
							string KM2 = add_parameter_marker(is_MATLAB_code, "KM2_" + promoter_name);
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*" + alpha + "*" + KM1 + "/(1 + " + KM1 + " + " + KM2 + "*" + tf_name + ");\n";
						}
						else if (promoter_name.compare("pTET") == STR_EQ || promoter_name.compare("pLtetO_1") == STR_EQ || promoter_name.compare("pTET_star_") == STR_EQ) {
							int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
							string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
							string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
							string KM1 = add_parameter_marker(is_MATLAB_code, "KM1_" + promoter_name);
							string KM2 = add_parameter_marker(is_MATLAB_code, "KM2_" + promoter_name);
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*" + alpha + "*" + KM1 + "/(1 + " + KM1 + " + 2*" + KM2 + "*" + tf_name + " + " + KM2 + "^2*" + tf_name + "^2" + ");\n";
						}
						else if (promoter_name.compare("pBAD") == STR_EQ) {
							int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
							string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
							string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
							string KM1 = add_parameter_marker(is_MATLAB_code, "KM1_" + promoter_name);
							string KM2 = add_parameter_marker(is_MATLAB_code, "KM2_" + promoter_name);
							string KM3 = add_parameter_marker(is_MATLAB_code, "KM3_" + promoter_name);
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*" + alpha + "*(" + KM1 + " + " + KM2 + "*(" + tf_name + "_free - " + tf_name + "))/(1 + " + KM1 + " + " + KM2 + "*(" + tf_name + "_free - " + tf_name + ")" + " + " + KM3 + "*" + tf_name + ");\n";
						}
						else {
							sbrome_error_print("Moon model does not include " + promoter_name);
						}
					}
					if (model_type == KineticModel::_MODEL_OPTION_HILL || is_hill_eq) {
						int pre_node_id = kinetic_model.GetInNodeList(node_id)[0];
						string tf_name = kinetic_model.GetNodeInfo(pre_node_id)->get_var_name();
						int regulation_type =  kinetic_model.GetEdgeInfo(pre_node_id, node_id)->type;
						string beta = add_parameter_marker(is_MATLAB_code, "beta_" + promoter_name);
						string K = add_parameter_marker(is_MATLAB_code, "K_" + promoter_name);
						string n = add_parameter_marker(is_MATLAB_code, "n_" + promoter_name);
						string pow_1 = K + "^" + n;
						string pow_2 = tf_name + "^" + n;
						if (regulation_type == BioNetEdgeInfo::ACTIVATORY)
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*(" + beta + " + (" + alpha + " - " + beta + ")*" + pow_2 + "/(" + pow_1 + " + " + pow_2 + "));\n";
						else
							simulation_code += promoter_name + " = " +  Float2Str(copy_number) + "*(" + beta + " + (" + alpha + " - " + beta + ")*" + pow_1 + "/(" + pow_1 + " + " + pow_2 + "));\n";
					}
				}
				break;
			}
			case MolecularSpecies::PROTEIN: {
				string protein_name = node_name;
				string RBS_strength = add_parameter_marker(is_MATLAB_code, correct_to_a_var_name("ALPHA_" + kinetic_model.GetNodeInfo(node_id)->aux_node_list[0].name));
				if (factor_num == 1) {
					int pre_node_id_1 = kinetic_model.GetInNodeList(node_id)[0];
					string promoter_name = kinetic_model.GetNodeInfo(pre_node_id_1)->get_var_name();
					string scale = (kinetic_model.GetOutputList()[0] == node_id)? (add_parameter_marker(is_MATLAB_code, "scale_" + node_name + "_" + correct_to_a_var_name(paper_ref.pub_id)) + "*") : "";
					simulation_code += protein_name + " = " + scale + RBS_strength + "*" + promoter_name + ";\n";
					if (kinetic_model.GetOutputList()[0] == node_id && is_MATLAB_code) {
						if (kinetic_model.GetInputList().empty())
							simulation_code += "\toutput = " + protein_name + ";\n";
						else
							simulation_code += "\toutput(i) = " + protein_name + ";\n";
					}
				}
				else { // this must be the case where a ligand binds and stops the function of a protein
					int pre_node_id_1 = kinetic_model.GetInNodeList(node_id)[0];
					int pre_node_id_2 = kinetic_model.GetInNodeList(node_id)[1];
					string promoter_name, ligand_name;
					if (kinetic_model.GetEdgeInfo(pre_node_id_1, node_id)->type == BioNetEdgeInfo::ACTIVATORY) {
						promoter_name = kinetic_model.GetNodeInfo(pre_node_id_1)->get_var_name();
						ligand_name = kinetic_model.GetNodeInfo(pre_node_id_2)->get_var_name();
					}
					else {
						promoter_name = kinetic_model.GetNodeInfo(pre_node_id_2)->get_var_name();
						ligand_name = kinetic_model.GetNodeInfo(pre_node_id_1)->get_var_name();
					}
					string n_ligand = add_parameter_marker(is_MATLAB_code, "n_" + ligand_name);
					string k_ligand = add_parameter_marker(is_MATLAB_code, "K_" + ligand_name);
					string pow_ligand = k_ligand + "^" + n_ligand;
					simulation_code += protein_name + " = " + RBS_strength + "*" + promoter_name + "*" + pow_ligand + "/(" + pow_ligand + " + " + ligand_name + "^" + n_ligand + ");\n";
					if (model_type == KineticModel::_MODEL_OPTION_ELLIS) {
						if (protein_name.compare("lacI") == STR_EQ)
							simulation_code += "\t" + protein_name + "_free = " + RBS_strength + "*" + promoter_name +  + ";\n";
					}
					if (model_type == KineticModel::_MODEL_OPTION_MOON) {
						if (protein_name.compare("araC") == STR_EQ)
							simulation_code += "\t" + protein_name + "_free = " + RBS_strength + "*" + promoter_name + ";\n";
					}
				}
				break;
			}
			case MolecularSpecies::LIGAND_PROTEIN_COMPLEX: {
				string ligand_protein = node_name;
				if (factor_num != 2) {
					sbrome_error_print("A ligand-protein binding reaction must have two reactants!\n");
				}
				else {
					int pre_node_id_1 = kinetic_model.GetInNodeList(node_id)[0];
					int pre_node_id_2 = kinetic_model.GetInNodeList(node_id)[1];
					string protein_name;
					string ligand_name;
					if (kinetic_model.GetNodeInfo(pre_node_id_1)->type == MolecularSpecies::PROTEIN) {
						protein_name = kinetic_model.GetNodeInfo(pre_node_id_1)->get_var_name();
						ligand_name = kinetic_model.GetNodeInfo(pre_node_id_2)->get_var_name();
					}
					else {
						protein_name = kinetic_model.GetNodeInfo(pre_node_id_2)->get_var_name();
						ligand_name = kinetic_model.GetNodeInfo(pre_node_id_1)->get_var_name();
					}
					string n_ligand = add_parameter_marker(is_MATLAB_code, "n_" + ligand_name);
					string k_ligand = add_parameter_marker(is_MATLAB_code, "K_" + ligand_name);
					string pow_ligand = k_ligand + "^" + n_ligand;
					simulation_code += ligand_protein + " = " + protein_name + "*" + pow_ligand + "/(" + pow_ligand + " + " + ligand_name + "^" + n_ligand + ");\n";
				}
				break;
			}
			case MolecularSpecies::LIGAND:
				// Replace ligand with the input value
				if (is_MATLAB_code)
					simulation_code += node_name + " = input(i);\n";
				// This needs to be extended for the case protein ---> ligand
				break;
			default:
				sbrome_error_print("Error: There is no type " + Int2Str(node_type) + "\n");
				break;
		}
	}
	if (!is_MATLAB_code) {
		StringList stmt_list;
		split_string(simulation_code, ';', stmt_list);
		simulation_code = "{";
		for (unsigned int i = 0; i < stmt_list.size(); i++) {
			stmt_list[i].erase(std::remove(stmt_list[i].begin(), stmt_list[i].end(), '\n'), stmt_list[i].end());
			stmt_list[i].erase(std::remove(stmt_list[i].begin(), stmt_list[i].end(), '\t'), stmt_list[i].end());
			if (!stmt_list[i].empty())
				simulation_code += ((simulation_code.length() > 1)? "; '": "'") + stmt_list[i] + "'";
		}
		simulation_code += "}";
	}
	return simulation_code;
}

string ModelAndExpData::export_MATLAB_simulation_code(IdMap combined_parameter_name_map, int model_type) {
	IdMap parameter_name_map = extract_parameter_name_map(model_type);
	string parameter_assignment_str = "";
	// Parameter assignment
	for (map<string,int>::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++) {
		string parameter_name = correct_to_a_var_name(it->first);
		parameter_assignment_str += parameter_name + " = x(" + Int2Str(combined_parameter_name_map[it->first] + 1) + ");\n";
	}
	// Simulation
	string simulation_code_str = write_simulation_code(model_type, true);
	// function description
	add_tab_to_each_sentence(parameter_assignment_str);
	string function_name = "simf_" + correct_to_a_var_name(get_id());
	string MATLAB_code;
	if (kinetic_model.GetInputList().empty()) {
		MATLAB_code = "function output = " + function_name + "(x)\n";
		MATLAB_code += parameter_assignment_str + simulation_code_str;
		MATLAB_code += "end\n";
		// error function
		MATLAB_code += "function output = " + function_name + "_error(x)\n";
		MATLAB_code += "\toutput = repmat(" + function_name + "(x) - " + Float2Str(exp_dataset.output_measurement[0][0].mean) + ",1,5);\n";
		MATLAB_code += "end\n";
	}
	else {
		MATLAB_code = "function output = " + function_name + "(x, input)\n";
		MATLAB_code += "\toutput = zeros(1,length(input));\n";
		MATLAB_code += parameter_assignment_str;
		MATLAB_code += "\tfor i = 1:length(input)\n";
		add_tab_to_each_sentence(simulation_code_str);
		MATLAB_code += simulation_code_str + "\tend\nend\n";
		// error function
		string input_str_tmp = "[" + Float2Str(exp_dataset.input_supplement[0][0]);
		for (unsigned int input_row = 1; input_row < exp_dataset.input_supplement.size(); input_row++)
			input_str_tmp += " " + Float2Str(exp_dataset.input_supplement[input_row][0]);
		input_str_tmp += "]";
		int output_node_id = kinetic_model.GetOutputList()[0];
		string output_name = kinetic_model.GetNodeInfo(output_node_id)->get_var_name();
		string output_str_tmp = "[" + Float2Str(exp_dataset.output_measurement[0][0].mean);
		for (unsigned int output_row = 1; output_row < exp_dataset.output_measurement.size(); output_row++)
			output_str_tmp += " " + Float2Str(exp_dataset.output_measurement[output_row][0].mean);
		output_str_tmp += "]";
		MATLAB_code += "function output = " + function_name + "_error(x)\n";
		MATLAB_code += "\toutput = repmat(" + function_name + "(x," + input_str_tmp + ") - " + output_str_tmp + ",1,5);\n";
		MATLAB_code += "end\n";
	}
	return MATLAB_code;
}

string ModelAndExpData::get_id() {
	return paper_ref.pub_id + "_" + experiment_id;
}

PAMLib::PAMLib(string directory_path) {
	import_from_a_folder(directory_path);
	dna_component_table.propagate_family_information(reaction_table.variant_to_family);
	//print(&std::cout);
	//dna_component_table.print(&std::cout);
}

void PAMLib::import_from_a_folder(string folder_path_str) {
	boost::filesystem::path folder_path(folder_path_str);
	if (!boost::filesystem::exists(folder_path))
		sbrome_error_print("There is no folder " + folder_path.string() + "\n");
	boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
	for (boost::filesystem::directory_iterator itr(folder_path); itr != end_itr; ++itr) {
		if (is_directory(itr->status()))
	      import_from_a_folder(itr->path().string());
	    else if (is_regular_file(itr->status()))
	    	import_from_a_file(itr->path().string());
	    else
	    	sbrome_error_print("Can not process the folder or the file: " + itr->path().filename().string());
	  }
}

void PAMLib::import_from_a_file(string filename) {
	if (filename[filename.length() - 1] == '~')
		return;
	//cout << filename << endl;
	boost::property_tree::ptree pt;
	boost::property_tree::read_json(filename, pt);
	Publication publication;
	publication.read(&pt);
	publication_list.push_back(publication);
	for (vector<DnaComponent>::iterator dna_component = publication.dna_component_list.begin(); dna_component != publication.dna_component_list.end(); dna_component++)
		dna_component_table.add_component(&(*dna_component));
	for (vector<Reaction>::iterator reaction = publication.reaction_list.begin(); reaction != publication.reaction_list.end(); reaction++)
		reaction_table.reaction_list.push_back(*reaction);
	for (vector<Strain>::iterator strain = publication.strain_list.begin(); strain != publication.strain_list.end(); strain++)
		strain_table[strain->name] = (*strain);
}

void PAMLib::print(ostream* f, int detail_level) {
	if (detail_level == ALL) {
		for (unsigned int i = 0; i < publication_list.size(); i++) {
			*f << "===============================" << endl;
			publication_list[i].print(f);
		}
	}
	else {
		int data_point_num = 0;
		*f << "# Publications: " << publication_list.size() << endl;
		int paper_with_exp_data_num = 0, exp_number = 0, exp_curve_num = 0, exp_point_num = 0;
		IdMap origin_freq, ar_freq, promoter_freq, RBS_freq, CDS_freq, terminator_freq;
		string experiment_name_list = "";
		for (vector<Publication>::iterator paper = publication_list.begin(); paper != publication_list.end(); paper++) {
			*f << "#";
			//paper->reference.print(f);
			if (!paper->experiment_list.empty())
				paper_with_exp_data_num++;
			for (vector<Experiment>::iterator experiment = paper->experiment_list.begin(); experiment != paper->experiment_list.end(); experiment++) {
				experiment_name_list += paper->reference.pub_id + "_" + experiment->reference + "\n";
				exp_number++;
				int exp_data_point_num = experiment->dataset.output_measurement.size();
				if (exp_data_point_num > 1)
					exp_curve_num++;
				else
					exp_point_num++;
				data_point_num += exp_data_point_num;
				for (StringList::iterator plasmid_name = experiment->genotype.plasmid_name_list.begin(); plasmid_name != experiment->genotype.plasmid_name_list.end(); plasmid_name++) {
					vector<DnaFeature> dna_feature_list;
					dna_component_table.get_full_feature(*plasmid_name, dna_feature_list);
					for (vector<DnaFeature>::iterator dna_feature = dna_feature_list.begin(); dna_feature != dna_feature_list.end(); dna_feature++) {
						DnaComponent tmp;
						dna_component_table.get_component(dna_feature->name, tmp);
						if (tmp.type.compare("origin") == STR_EQ)
							frequency_counting(tmp.name, origin_freq);
						else if (tmp.type.compare("Antibiotic Resistance") == STR_EQ)
							frequency_counting(tmp.name, ar_freq);
						else if (tmp.type.compare("promoter") == STR_EQ)
							frequency_counting(tmp.name, promoter_freq);
						else if (tmp.type.compare("RBS") == STR_EQ)
							frequency_counting(tmp.name, RBS_freq);
						else if (tmp.type.compare("CDS") == STR_EQ)
							frequency_counting(tmp.name, CDS_freq);
						else if (tmp.type.compare("terminator") == STR_EQ)
							frequency_counting(tmp.name, terminator_freq);
						else
							sbrome_error_print("Part type " + tmp.type + " of (" + dna_feature->name + ") not found!\n");
					}
				}
			}
		}
		*f << "# Paper with exp data: " << paper_with_exp_data_num << "\t# Experiment: " << exp_number << "\tCurve: " << exp_curve_num << "\tSingle point:" << exp_point_num << "\t#Data point: " << data_point_num << endl;
		cout << experiment_name_list << endl;
		// Origin
		for (IdMap::iterator it = origin_freq.begin(); it != origin_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
		for (IdMap::iterator it = ar_freq.begin(); it != ar_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
		for (IdMap::iterator it = promoter_freq.begin(); it != promoter_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
		for (IdMap::iterator it = RBS_freq.begin(); it != RBS_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
		for (IdMap::iterator it = CDS_freq.begin(); it != CDS_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
		for (IdMap::iterator it = terminator_freq.begin(); it != terminator_freq.end(); it++)
			*f << it->first << "\t" << it->second << endl;
	}
}

void PAMLib::get_strain_information(string strain_name, StringList& modification_info) {
	if (strain_name.compare("K-12") == STR_EQ)
		return;
	else if (strain_table.find(strain_name) == strain_table.end())
		sbrome_error_print("Strain " + strain_name + " not found!\n");
	else {
		StringList tmp;
		get_strain_information(strain_table[strain_name].parent, tmp);
		modification_info.clear();
		for (StringList::iterator it = tmp.begin(); it != tmp.end(); it++)
			modification_info.push_back((*it));
		for (StringList::iterator it = strain_table[strain_name].gene_modification_list.begin(); it != strain_table[strain_name].gene_modification_list.end(); it++)
			modification_info.push_back((*it));
	}
}

void PAMLib::convert_genotype_to_biological_network(Genotype genotype, BiologicalNetwork& logic_model) {
	// Add all regulations of the chromosome
	StringList strain_info;
	get_strain_information(genotype.host_strain_name, strain_info);
	bool is_lacZ_knocked_out = false, is_lacI_knocked_out = false, is_araC_knocked_out = false;
	map<string,bool> plasmid_in_chromosome_or_not;
	for (StringList::iterator it = strain_info.begin(); it != strain_info.end(); it++) {
		if (it->compare("delta lacZ") == STR_EQ)
			is_lacZ_knocked_out = true;
		if (it->compare("delta lacX74") == STR_EQ || it->compare("delta lacZYA-argF") == STR_EQ || it->compare("lacI22") == STR_EQ)
			is_lacI_knocked_out = true;
		if (it->compare("delta (ara, leu)7697") == STR_EQ || it->compare("delta (ara-leu)7697") == STR_EQ || it->compare("delta araC") == STR_EQ)
			is_araC_knocked_out = true;
		if (it->find("integrate ") != string::npos)
			plasmid_in_chromosome_or_not[it->substr(10, it->size() - 1)] = true;
	}
	if (!is_lacI_knocked_out) {	// lac operon
		int n1 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::mRNA, "placI", 1));
		int n2 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::PROTEIN, "lacI", 1, MolecularSpecies::UNKNOWN_MOLECULE, "RBS-placI", 1));
		logic_model.InsertEdge(n1, n2, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
		if (!is_lacZ_knocked_out) {
			int n3 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::mRNA, "pLAC", 1));
			int n4 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::PROTEIN, "lacZ", 1, MolecularSpecies::UNKNOWN_MOLECULE, "RBS-pLAC", 1));
			logic_model.InsertEdge(n3, n4, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
		}
	}
	if (!is_araC_knocked_out) {	// L-arabinose operon
		int n1 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::mRNA, "pC", 1));
		int n2 = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::PROTEIN, "araC", 1, MolecularSpecies::UNKNOWN_MOLECULE, "RBS-pC", 1));
		logic_model.InsertEdge(n1, n2, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));

	}
	for (StringList::iterator plasmid_name = genotype.plasmid_name_list.begin(); plasmid_name != genotype.plasmid_name_list.end(); plasmid_name++)
		plasmid_in_chromosome_or_not[(*plasmid_name)] = false;
	// Add more regulations from plasmids
	for (map<string,bool>::iterator it = plasmid_in_chromosome_or_not.begin(); it != plasmid_in_chromosome_or_not.end(); it++) {
		IdConverter part_id_to_node_id;
		vector<DnaFeature> dna_feature_list;
		dna_component_table.get_full_feature(it->first, dna_feature_list);
		double copy_number = (it->second)? 1:dna_component_table.get_copy_number(&dna_feature_list);
		// TEST HERE
		//for (unsigned int i = 0; i < dna_feature_list.size(); i++)
		//	cout << dna_feature_list[i].name << "\t";
		//cout << endl;
		// END TEST
		// Add nodes
		for (unsigned int i = 0; i < dna_feature_list.size(); i++) {
			string part_name = dna_feature_list[i].name;
			DnaComponent dna_comp;
			dna_component_table.get_component(part_name, dna_comp);
			int part_type = DnaComponent::Name2Id(dna_comp.type);
			if (part_type == DnaComponent::PROMOTER || part_type == DnaComponent::CDS)
				part_id_to_node_id[i] = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::PartType2MoleculeType(part_type), part_name, copy_number));
		}
		// Add edges
		for (unsigned int i = 0; i < dna_feature_list.size(); i++) {
			string part_name = dna_feature_list[i].name;
			DnaComponent dna_comp;
			dna_component_table.get_component(part_name, dna_comp);
			int part_type = DnaComponent::Name2Id(dna_comp.type);
			if (part_type == DnaComponent::PROMOTER) {
				int step = (dna_feature_list[i].direction.compare("backward") == STR_EQ)?-1:1;
				string current_RBS;
				for (unsigned int j = i + step; j >= 0 && j < dna_feature_list.size(); j += step) {
					string downstream_part_name = dna_feature_list[j].name;
					DnaComponent downstream_dna_comp;
					dna_component_table.get_component(downstream_part_name, downstream_dna_comp);
					int dest_type = DnaComponent::Name2Id(downstream_dna_comp.type);
					if (dest_type == DnaComponent::RBS)
						current_RBS = downstream_part_name;
					else if (dest_type == DnaComponent::CDS) {
						logic_model.GetNodeInfo(part_id_to_node_id[j])->aux_node_list.push_back(create_a_bio_net_node_info(MolecularSpecies::UNKNOWN_MOLECULE,current_RBS,1));
						logic_model.InsertEdge(part_id_to_node_id[i], part_id_to_node_id[j], create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
						current_RBS = "";
					}
					else
						break;
				}
			}
		}
	}
}

void PAMLib::add_IO_node_to_biological_network(Environment environment, Measurement measurement, BiologicalNetwork& logic_model) {
	// Add inputs
	logic_model.ClearInput();
	for (unsigned int i = 0; i < environment.supplemented_chemical_name_list.size(); i++) {
		string input_molecular_species = environment.supplemented_chemical_name_list[i].first;
		if (!logic_model.IsNodeNameAvailable(input_molecular_species)) {
			int node_id_tmp = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::LIGAND, input_molecular_species, 1));
			logic_model.InsertInput(node_id_tmp);
		}
		else
			logic_model.InsertInput(logic_model.GetFirstNodeId(input_molecular_species));
	}
	// Add outputs
	logic_model.ClearOutput();
	for (unsigned int i = 0; i < measurement.molecular_species_list.size(); i++) {
		string output_molecular_species = measurement.molecular_species_list[i].first;
		if (logic_model.IsNodeNameAvailable(output_molecular_species))
			logic_model.InsertOutput(logic_model.GetFirstNodeId(output_molecular_species));
		else
			sbrome_error_print("The output " + output_molecular_species + " is irrelavant with the model \n");
	}
}

void PAMLib::add_molecular_species_interaction_to_biological_network(BiologicalNetwork& logic_model) {
	// Add ligand or complex nodes which are produced from proteins
	/*
	bool is_new_species_added;
	do {
		is_new_species_added = false;
		for (map<string,int>::iterator it = available_species.begin(); it != available_species.end(); it++) {
			if (product_list_.find(it->first) != product_list_.end()) {
				StringList product_list = product_list_.at(it->first);
				for (int i = 0; i < product_list.size(); i++) {
					int product_type = species_name_to_type_.at(product_list[i]);
					if (available_species.find(product_list[i]) == available_species.end() &&
						(product_type == LIGAND || product_type == RNA_COMPLEX || product_type == PROTEIN_COMPLEX || product_type == LIGAND_PROTEIN_COMPLEX)) {
							StringList source_list = source_list_.at(product_list[i]);
							bool are_all_source_available = true;
							for (int j = 0; j < source_list.size(); j++) {
								if (available_species.find(source_list[j]) == available_species.end()) {
									are_all_source_available = false;
									break;
								}
							}
							if (are_all_source_available) {
								available_species[product_list[i]] = module_editor->InsertNode(new BioNetNode(species_name_to_type_.at(product_list[i]),product_list[i]));
								is_new_species_added = true;
							}
					}
				}
			}
		}
	}
	while (is_new_species_added);
	*/
	//logic_model.print(&std::cout);

	// Add edges
	//model_editor->print(&std::cout);
	map<StringPair,bool> reaction_tracking;
	bool is_new_reaction_added;
	do {
		is_new_reaction_added = false;
		IdList node_id_list = logic_model.GetNodeList();
		for (IdList::iterator src_node_id = node_id_list.begin(); src_node_id != node_id_list.end(); src_node_id++)
			for (IdList::iterator dest_node_id = src_node_id; dest_node_id != node_id_list.end(); dest_node_id++) {
				string src_name = logic_model.GetNodeInfo(*src_node_id)->name;
				string dest_name = logic_model.GetNodeInfo(*dest_node_id)->name;
				string interaction_type = reaction_table.find_interaction(src_name, dest_name);
				if (!interaction_type.empty() && reaction_tracking.find(make_pair(src_name, dest_name)) == reaction_tracking.end()) {
					reaction_tracking[make_pair(src_name, dest_name)] = true;
					is_new_reaction_added = true;
					if (interaction_type.compare("ligand-repressor") == STR_EQ) {
						if (logic_model.GetNodeInfo(*src_node_id)->type == MolecularSpecies::LIGAND)
							logic_model.InsertEdge(*src_node_id, *dest_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::INHIBITORY));
						else
							logic_model.InsertEdge(*dest_node_id, *src_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::INHIBITORY));
					}
					else if (interaction_type.compare("repressor-promoter") == STR_EQ) {
						int src_node_type = logic_model.GetNodeInfo(*src_node_id)->type;
						if (src_node_type == MolecularSpecies::PROTEIN || src_node_type == MolecularSpecies::LIGAND)
							logic_model.InsertEdge(*src_node_id, *dest_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::INHIBITORY));
						else
							logic_model.InsertEdge(*dest_node_id, *src_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::INHIBITORY));
					}
					else if (interaction_type.compare("activator-promoter") == STR_EQ) {
						if (logic_model.GetNodeInfo(*src_node_id)->type == MolecularSpecies::mRNA)
							logic_model.InsertEdge(*dest_node_id, *src_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
						else
							logic_model.InsertEdge(*src_node_id, *dest_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
					}
					else if (interaction_type.compare("ligand-activator") == STR_EQ) {
						string product = reaction_table.find_product(src_name, dest_name);
						int complex_node_id = logic_model.InsertNode(create_a_bio_net_node_info(MolecularSpecies::LIGAND_PROTEIN_COMPLEX, product, 1));
						logic_model.InsertEdge(*src_node_id, complex_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
						logic_model.InsertEdge(*dest_node_id, complex_node_id, create_a_bio_net_edge_info(BioNetEdgeInfo::ACTIVATORY));
					}
				}
			}
	}
	while (is_new_reaction_added);
}

void PAMLib::get_parameter_value_bound(map<string,int> parameter_name_map, StringList& parameter_name_list, ValueList& lb, ValueList& ub) {
	parameter_name_list.resize(parameter_name_map.size());
	lb.resize(parameter_name_map.size());
	ub.resize(parameter_name_map.size());
	for (map<string,int>::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++) {
		parameter_name_list[it->second] = it->first;
		StringList para_word_list;
		split_string(it->first, '_', para_word_list);
		/*if (it->first.compare("alpha_BBa_J23117") == STR_EQ) {
			lb[it->second] = 0.01;	// J23117
			ub[it->second] = 0.02;	// pL
		}
		else if (it->first.compare("alpha_BBa_J23114") == STR_EQ) {
			lb[it->second] = 0.02;	// J23117
			ub[it->second] = 0.03;	// pL
		}
		else if (it->first.compare("alpha_BBa_J23115") == STR_EQ) {
			lb[it->second] = 0.04;	// J23117
			ub[it->second] = 0.07;	// pL
		}
		else if (it->first.compare("alpha_BBa_J23105") == STR_EQ) {
			lb[it->second] = 0.09;	// J23117
			ub[it->second] = 0.12;	// pL
		}
		else if (it->first.compare("alpha_BBa_J23106") == STR_EQ) {
			lb[it->second] = 0.2;	// J23117
			ub[it->second] = 0.25;	// pL
		}
		else if (it->first.compare("alpha_BBa_J23101") == STR_EQ) {
			lb[it->second] = 0.95;	// J23117
			ub[it->second] = 1.05;	// pL
		}
		else if (it->first.compare("alpha_placIq") == STR_EQ) {
			lb[it->second] = 0.04;	// J23117
			ub[it->second] = 0.05;	// pL
		}
		else if (it->first.compare("ALPHA_BBa_B0030") == STR_EQ) {
			lb[it->second] = 0.95;
			ub[it->second] = 1.05;
		}
		else if (it->first.compare("ALPHA_BBa_B0032") == STR_EQ) {
			lb[it->second] = 0.45;
			ub[it->second] = 0.55;
		}
		else if (it->first.compare("ALPHA_BBa_B0033") == STR_EQ) {
			lb[it->second] = 0.015;
			ub[it->second] = 0.02;
		}
		else*/ if (para_word_list[0].compare("alpha") == STR_EQ) {
			lb[it->second] = 0.002;	// J23117
			ub[it->second] = 50;	// pL
		}
		else if (para_word_list[0].compare("beta") == STR_EQ) {
			lb[it->second] = 1e-4;
			ub[it->second] = 1;
		}
		else if (para_word_list[0].compare("ALPHA") == STR_EQ) {
			lb[it->second] = 0.03;	// B0030
			ub[it->second] = 5;		// B0034
		}
		else if (para_word_list[0].compare("K") == STR_EQ) {
			if (para_word_list[1].compare("p") == STR_EQ) {	// for protein binding affinity
				lb[it->second] = 1e-3;
				ub[it->second] = 1;
			}
			else {	// for ligand binding affinity
				lb[it->second] = 1e-3;
				ub[it->second] = 50;
			}
		}
		else if (para_word_list[0].compare("KE1") == STR_EQ) {
			lb[it->second] = 1e-3;
			ub[it->second] = 50;
		}
		else if (para_word_list[0].compare("KE2") == STR_EQ) {
			lb[it->second] = 1e-3;
			ub[it->second] = 50;
		}
		else if (para_word_list[0].compare("KM1") == STR_EQ) {
			lb[it->second] = 1e-3;
			ub[it->second] = 50;
		}
		else if (para_word_list[0].compare("KM2") == STR_EQ) {
			lb[it->second] = 1e-3;
			ub[it->second] = 50;
		}
		else if (para_word_list[0].compare("KM3") == STR_EQ) {
			lb[it->second] = 1e-3;
			ub[it->second] = 50;
		}
		else if (para_word_list[0].compare("n") == STR_EQ) {
			if (para_word_list[1].compare("p") == STR_EQ) {
				lb[it->second] = 1;
				ub[it->second] = 6;
			}
			else {
				lb[it->second] = 1;
				ub[it->second] = 4;
			}
		}
		else if (para_word_list[0].compare("scale") == STR_EQ) {
			lb[it->second] = 1e-2;
			ub[it->second] = 10;
		}

		else {
			lb[it->second] = EPSILON;
			ub[it->second] = 1;
			cout << "Error: There is no bound information for the parameter " << it->first << endl;
		}
	}
}

void PAMLib::convert_experiment_to_model_and_data(Experiment experiment, ModelAndExpData &model_and_exp_data) {
	KineticModel& kinetic_model = model_and_exp_data.kinetic_model;
	// Add product nodes from the chromosome and plasmids and the regulation from their order
	convert_genotype_to_biological_network(experiment.genotype, kinetic_model);
	// Add input nodes and output node
	add_IO_node_to_biological_network(experiment.environment, experiment.measurement, kinetic_model);
	// Add interactions between molecular species and extra ones
	add_molecular_species_interaction_to_biological_network(kinetic_model);
	// Add pool nodes to merge duplicated nodes
	kinetic_model.AddPoolNode();
	// Remove all irrelevant nodes
	kinetic_model.RemoveIrrelevantNode();
	// Add the experimental data
	model_and_exp_data.exp_dataset = experiment.dataset;
}

string PAMLib::get_math_model_type_name(int model_type) {
	switch (model_type) {
		case KineticModel::_MODEL_OPTION_ELLIS:
			return "Ellis";
		case KineticModel::_MODEL_OPTION_MOON:
			return "Moon";
		case KineticModel::_MODEL_OPTION_HILL:
			return "Hill";
		default:
			return "NA";
	}
}
string PAMLib::Export_MATLAB_code_for_simultaneous_fitting(vector<ModelAndExpData>* model_and_exp_data_list, IdMap* combined_parameter_name_map, int model_type) {
	string final_code_str;
	// call the solver
	StringList parameter_name_list;
	ValueList lb, ub;
	get_parameter_value_bound(*combined_parameter_name_map, parameter_name_list, lb, ub);
	final_code_str =  "function x_opt = simultaneous_fitting(ensemble_size, ite_num, func_ite_num)\n";
	// parameter name
	StringList parameter_name_list_tmp = parameter_name_list;
	for (unsigned int i = 0; i < parameter_name_list_tmp.size(); i++) {
		if (parameter_name_list_tmp[i].length() > 25)
			parameter_name_list_tmp[i] = parameter_name_list_tmp[i].substr(0, 24);
	}
	final_code_str += "\tparameter_name = {'" + parameter_name_list_tmp[0] + "'";
	for (unsigned int i = 1; i < lb.size(); i++)
		final_code_str += "; '" + parameter_name_list_tmp[i] + "'";
	final_code_str += "};\n";
	// lb
	final_code_str += "\tlb = [" + Float2Str(lb[0]);
	for (unsigned int i = 1; i < lb.size(); i++)
		final_code_str += " " + Float2Str(lb[i]);
	final_code_str += "];\n";
	// ub
	final_code_str += "\tub = [" + Float2Str(ub[0]);
	for (unsigned int i = 1; i < ub.size(); i++)
		final_code_str += " " + Float2Str(ub[i]);
	final_code_str += "];\n";
	// solver
	final_code_str += "\t[CI_lb, CI_ub, x_opt, solution_ensemble] = fit_a_model(ensemble_size, ite_num, func_ite_num, @simultaneous_fitting_error, lb, ub);\n";
	// For GNUPlot
	final_code_str += "\tparameter_fileID = fopen('Plot_data/sim_CI_" + get_math_model_type_name(model_type) + ".dat','w');\n";
	final_code_str += "\tfor i = 1:length(parameter_name)\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%s\\t', parameter_name{i});\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\t', lb(i));\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\t', ub(i));\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\t', CI_lb(i));\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\t', CI_ub(i));\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\n', x_opt(i));\n";
	final_code_str += "\tend\n";
	final_code_str += "\tfclose(parameter_fileID);\n";
	// For JSON
	final_code_str += "\tprint_parameter_to_json('JSON/inferred_parameter_value.json', lb, ub, CI_lb, CI_ub, x_opt, solution_ensemble, parameter_name);\n";
	// Sensitivity analysis
	final_code_str += "\tSA_fileID = fopen('Plot_data/sensitivity_analysis_" + get_math_model_type_name(model_type) + ".dat','w');\n";
	final_code_str += "\topt_error = norm(simultaneous_fitting_error(x_opt));\n";
	final_code_str += "\tpertubation = [0.01 0.05 0.1 0.5 2 10 50 100];\n";
	final_code_str += "\terror_diff = pertubation;\n";
	final_code_str += "\tfor i = 1:length(parameter_name)\n";
	final_code_str += "\t\tx = x_opt;\n";
	final_code_str += "\t\tfor j = 1:length(pertubation)\n";
	final_code_str += "\t\t\tx(i) = pertubation(j)*x_opt(i);\n";
	final_code_str += "\t\t\terror = norm(simultaneous_fitting_error(x));\n";
	final_code_str += "\t\t\terror_diff(j) = abs(error - opt_error)/opt_error;\n";
	final_code_str += "\t\tend\n";
	final_code_str += "\t\tfprintf(SA_fileID, '%s\t', parameter_name{i});\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\t', mean(error_diff));\n";
	final_code_str += "\t\tfprintf(parameter_fileID, '%f\\n', std(error_diff)/sqrt(length(error_diff)));\n";
	final_code_str += "\tend\n";
	final_code_str += "\tfclose(SA_fileID);\n";
	final_code_str += "end\n";
	// error function
	int data_point_number = 0;
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++)
		data_point_number += model_it->exp_dataset.input_supplement.size();
	final_code_str += "function error = simultaneous_fitting_error(x)\n";
	final_code_str += "\terror = zeros(" + Int2Str(data_point_number) + ",1);\n";
	final_code_str += "\tcurrent_base_index = 0;\n";
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		final_code_str += "\t%==============\n";
		// Write the input list
		double max_output_value = 1e-10;
		for (unsigned int i = 0; i < model_it->exp_dataset.output_measurement.size(); i++)
			for (unsigned int j = 0; j < model_it->exp_dataset.output_measurement[i].size(); j++)
				if (max_output_value < model_it->exp_dataset.output_measurement[i][j].mean)
					max_output_value = model_it->exp_dataset.output_measurement[i][j].mean;
		double weight = 1/max_output_value;
		string function_name = correct_to_a_var_name(model_it->get_id());
		if (model_it->kinetic_model.GetInputList().size() == 0) {
			final_code_str += "\terror(current_base_index + 1) = " + Float2Str(weight) + "*(simf_" + function_name + "(x) - " + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean) + ");\n";
			final_code_str += "\tcurrent_base_index = current_base_index + 1;\n";
		}
		else if (model_it->kinetic_model.GetInputList().size() == 1) {
			int input_node_id = model_it->kinetic_model.GetInputList()[0];
			final_code_str += "\t" + model_it->kinetic_model.GetNodeInfo(input_node_id)->get_var_name() + " = [" + Float2Str(model_it->exp_dataset.input_supplement[0][0]);
			for (unsigned int input_row = 1; input_row < model_it->exp_dataset.input_supplement.size(); input_row++)
				final_code_str += " " + Float2Str(model_it->exp_dataset.input_supplement[input_row][0]);
			final_code_str += "];\n";
			// Write the output list
			int output_node_id = model_it->kinetic_model.GetOutputList()[0];
			string output_name = model_it->kinetic_model.GetNodeInfo(output_node_id)->get_var_name();
			final_code_str += "\t" + output_name + " = [" + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean);
			for (unsigned int output_row = 1; output_row < model_it->exp_dataset.output_measurement.size(); output_row++)
				final_code_str += " " + Float2Str(model_it->exp_dataset.output_measurement[output_row][0].mean);
			final_code_str += "];\n";
			final_code_str += "\toutput = " + Float2Str(weight) + "*(simf_" + function_name + "(x, " + model_it->kinetic_model.GetNodeInfo(input_node_id)->get_var_name() + ") - " + output_name + ");\n";
			final_code_str += "\tfor i = 1:length(output)\n";
			final_code_str += "\t\terror(current_base_index + i) = output(i);\n";
			final_code_str += "\tend\n";
			final_code_str += "\tcurrent_base_index = current_base_index + length(output);\n";
		}
	}
	final_code_str += "end\n";
	return final_code_str;

}

string PAMLib::Export_MATLAB_code_for_simultaneous_fitting_LOOCV(vector<ModelAndExpData>* model_and_exp_data_list, IdMap* combined_parameter_name_map, int model_type) {
	string sim_LOOCV_str;
	IdMap parameter_appearance_count;
	int total_data_point_number = 0;
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		total_data_point_number += model_it->exp_dataset.output_measurement.size();
		IdMap parameter_name_map = model_it->extract_parameter_name_map(model_type);
		for (IdMap::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++)
			if (parameter_appearance_count.find(it->first) == parameter_appearance_count.end())
					parameter_appearance_count[it->first] = 1;
				else
					parameter_appearance_count[it->first]++;
	}

	//sbrome_print(&std::cout, parameter_appearance_count);

	StringList parameter_name_list;
	ValueList lb, ub;
	get_parameter_value_bound(*combined_parameter_name_map, parameter_name_list, lb, ub);
	sim_LOOCV_str =  "function sim_LOOCV(ensemble_size, ite_num, func_ite_num)\n";
	// lb
	sim_LOOCV_str += "\tlb = [" + Float2Str(lb[0]);
	for (unsigned int i = 1; i < lb.size(); i++)
		sim_LOOCV_str += " " + Float2Str(lb[i]);
	sim_LOOCV_str += "];\n";
	// ub
	sim_LOOCV_str += "\tub = [" + Float2Str(ub[0]);
	for (unsigned int i = 1; i < ub.size(); i++)
		sim_LOOCV_str += " " + Float2Str(ub[i]);
	sim_LOOCV_str += "];\n";
	sim_LOOCV_str += "\tLOOCV_sim_data_fileID = fopen('Plot_data/LOOCV_sim_data_" + get_math_model_type_name(model_type) + ".dat','w');\n";

	string sim_LOOCV_error_str = "";
	string sim_LOOCV_fitting_str = "";
	StringList LOO_model_name_list;
	int current_data_point_number = 0;
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		bool is_CV = true;
		IdMap parameter_name_map = model_it->extract_parameter_name_map(model_type);
		for (map<string,int>::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++)
			if (parameter_appearance_count[it->first] == 1) {
				is_CV = false;
				break;
			}
		if (is_CV) {
			string function_name = correct_to_a_var_name(model_it->get_id());
			LOO_model_name_list.push_back(function_name);
			// Fitting
			sim_LOOCV_fitting_str += "\t[x, y, x_opt, z] = fit_a_model(ensemble_size, ite_num, func_ite_num, @LOOCV_sim_error_" + function_name + ", lb, ub);\n";
			// Validation
			if (model_it->kinetic_model.GetInputList().size() == 0) {
				sim_LOOCV_fitting_str += "\tdesired_output = " + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean) + ";\n";
				sim_LOOCV_fitting_str += "\toutput = simf_" + function_name + "(x_opt);\n";
			}
			else if (model_it->kinetic_model.GetInputList().size() == 1) {
				int input_node_id = model_it->kinetic_model.GetInputList()[0];
				string signal_name = model_it->kinetic_model.GetNodeInfo(input_node_id)->get_var_name();
				sim_LOOCV_fitting_str += "\t" + signal_name + " = [" + Float2Str(model_it->exp_dataset.input_supplement[0][0]);
				for (unsigned int input_row = 1; input_row < model_it->exp_dataset.input_supplement.size(); input_row++)
					sim_LOOCV_fitting_str += " " + Float2Str(model_it->exp_dataset.input_supplement[input_row][0]);
				sim_LOOCV_fitting_str += "];\n";
				// Write the output list
				sim_LOOCV_fitting_str += "\tdesired_output = [" + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean);
				for (unsigned int output_row = 1; output_row < model_it->exp_dataset.output_measurement.size(); output_row++)
					sim_LOOCV_fitting_str += " " + Float2Str(model_it->exp_dataset.output_measurement[output_row][0].mean);
				sim_LOOCV_fitting_str += "];\n";
				sim_LOOCV_fitting_str += "\toutput = simf_" + function_name + "(x_opt, " + signal_name + ");\n";
			}
			sim_LOOCV_fitting_str += "\tfprintf(LOOCV_sim_data_fileID, '#" + function_name + "\\n');\n";
			sim_LOOCV_fitting_str += "\tfor i = 1:length(output)\n";
			sim_LOOCV_fitting_str += "\t\tfprintf(LOOCV_sim_data_fileID, '%f\\t%f\\n', desired_output(i), output(i));\n";
			sim_LOOCV_fitting_str += "\tend\n";
			// Generate error functions
			sim_LOOCV_error_str += "function error = LOOCV_sim_error_" + function_name + "(x)\n";
			sim_LOOCV_error_str += "\terror_all = simultaneous_fitting_error(x) ;\n";
			sim_LOOCV_error_str += "\terror = remove_a_sub_array(error_all, " + Int2Str(current_data_point_number) + ", " + Int2Str(current_data_point_number + model_it->exp_dataset.output_measurement.size() - 1) + ");\n";
			sim_LOOCV_error_str += "end\n";
			current_data_point_number += model_it->exp_dataset.output_measurement.size();
		}
	}
	cout << "# LOOCV model: " << LOO_model_name_list.size() << endl;
	for (unsigned int i = 0; i < LOO_model_name_list.size(); i++)
		cout << LOO_model_name_list[i] << endl;
	sim_LOOCV_str += "\tLOO_model_name_list = {'" + LOO_model_name_list[0] + "'";
	for (unsigned int i = 1; i < LOO_model_name_list.size(); i++)
		sim_LOOCV_str += ", '" + LOO_model_name_list[i] + "'";
	sim_LOOCV_str += "};\n";
	sim_LOOCV_str += sim_LOOCV_fitting_str;
	sim_LOOCV_str += "\tfclose(LOOCV_sim_data_fileID);\n";
	sim_LOOCV_str += "end\n";
	sim_LOOCV_str += sim_LOOCV_error_str;
	return sim_LOOCV_str;
}

string PAMLib::Export_MATLAB_code_for_printing_simultaneous_fitting_result(vector<ModelAndExpData>* model_and_exp_data_list, IdMap* combined_parameter_name_map, int model_type) {
	string print_data_function_str = "function print_plotting_sim_data(sim_x_opt)\n";
	print_data_function_str += "\tsim_pair_fileID = fopen('Plot_data/sim_pair_" + get_math_model_type_name(model_type) + ".dat','w');\n";
	print_data_function_str += "\tjson_fileID = fopen('JSON/simulation_data_" + get_math_model_type_name(model_type) + ".json','w');\n";
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		print_data_function_str += "\t%%%%%%%%%%%%%\n";
		// Write the input list
		string function_name = correct_to_a_var_name(model_it->get_id());
		if (model_it->kinetic_model.GetInputList().size() == 0) {
			// For GNUPlot
			print_data_function_str += "\toutput_sim = simf_" + function_name + "(sim_x_opt);\n";
			print_data_function_str += "\tpoint_pair = [output_sim \t" + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean) + "];\n";
			print_data_function_str += "\tfprintf(sim_pair_fileID, '%f \\t %f\\n', point_pair);\n";
			// For JSON
			print_data_function_str += "\tprint_simulated_data_to_json(json_fileID, '" + model_it->paper_ref.pub_id + "', '" + model_it->experiment_id + "', [], output_sim, " + model_it->write_simulation_code(model_type) + ");\n";
		}
		else if (model_it->kinetic_model.GetInputList().size() == 1) {
			int input_node_id = model_it->kinetic_model.GetInputList()[0];
			string signal_name = model_it->kinetic_model.GetNodeInfo(input_node_id)->get_var_name();
			print_data_function_str += "\t" + signal_name + " = [" + Float2Str(model_it->exp_dataset.input_supplement[0][0]);
			for (unsigned int input_row = 1; input_row < model_it->exp_dataset.input_supplement.size(); input_row++)
				print_data_function_str += " " + Float2Str(model_it->exp_dataset.input_supplement[input_row][0]);
			print_data_function_str += "];\n";
			// Write the output list
			int output_node_id = model_it->kinetic_model.GetOutputList()[0];
			string output_name = model_it->kinetic_model.GetNodeInfo(output_node_id)->get_var_name();
			print_data_function_str += "\t" + output_name + " = [" + Float2Str(model_it->exp_dataset.output_measurement[0][0].mean);
			for (unsigned int output_row = 1; output_row < model_it->exp_dataset.output_measurement.size(); output_row++)
				print_data_function_str += " " + Float2Str(model_it->exp_dataset.output_measurement[output_row][0].mean);
			print_data_function_str += "];\n";
			// For GNUPlot
			print_data_function_str += "\toutput_sim = simf_" + function_name + "(sim_x_opt," + signal_name + ");\n";
			print_data_function_str += "\tsim_pair_matrix = [output_sim; " + output_name + "];\n";
			print_data_function_str += "\tfprintf(sim_pair_fileID, '%f \\t %f\\n', sim_pair_matrix);\n";
			// For JSON
			print_data_function_str += "\tprint_simulated_data_to_json(json_fileID, '" + model_it->paper_ref.pub_id + "', '" + model_it->experiment_id + "', " + signal_name + ", output_sim, " + model_it->write_simulation_code(model_type) + ");\n";
		}
		else {
			sbrome_error_print("Model has more than input!\n");
		}
	}
	print_data_function_str += "\tfclose(sim_pair_fileID);\n";
	print_data_function_str += "\tfclose(json_fileID);\n";
	print_data_function_str += "end\n";
	return print_data_function_str;
}

string PAMLib::Export_MATLAB_code_for_sequential_fitting(vector<ModelAndExpData>* model_and_exp_data_list, IdMap* combined_parameter_name_map, int model_type) {
	string sequential_fitting_prob_str = "function parameter_opt = sequential_fitting(N1, N2)\n";
	sequential_fitting_prob_str += "\terror_opt = 1e6;\n";
	sequential_fitting_prob_str += "\tparameter_opt = zeros(1," + Int2Str(combined_parameter_name_map->size()) + ");\n";
	sequential_fitting_prob_str += "\tfor i = 1:N1\n";
	sequential_fitting_prob_str += "\t\tmodel_order = randperm(" + Int2Str(model_and_exp_data_list->size()) + ");\n";
	sequential_fitting_prob_str += "\t\tparameter = sequential_fitting_for_one_order(model_order, N2);\n";
	sequential_fitting_prob_str += "\t\terror = norm(simultaneous_fitting_error(parameter));\n";
	sequential_fitting_prob_str += "\t\tif (error < error_opt)\n";
	sequential_fitting_prob_str += "\t\t\terror_opt = error;\n";
	sequential_fitting_prob_str += "\t\t\tparameter_opt = parameter;\n";
	sequential_fitting_prob_str += "\t\tend\n";
	sequential_fitting_prob_str += "\tend\n";
	StringList parameter_name_list(combined_parameter_name_map->size());
	for (map<string,int>::iterator it = combined_parameter_name_map->begin(); it != combined_parameter_name_map->end(); it++)
		parameter_name_list[it->second] = it->first;
	sequential_fitting_prob_str += "\tparameter_name_list = {'" + parameter_name_list[0] + "'";
	for (unsigned int i = 1; i < parameter_name_list.size(); i++)
		sequential_fitting_prob_str += "; '" + parameter_name_list[i] + "'";
	sequential_fitting_prob_str += "};\n";
	sequential_fitting_prob_str += "\tparameter_fileID = fopen('Plot_data/sequential_parameter.dat','w');\n";
	sequential_fitting_prob_str += "\tfor i = 1:length(parameter_name_list)\n";
	sequential_fitting_prob_str += "\t\tfprintf(parameter_fileID, '%s\\t', parameter_name_list{i});\n";
	sequential_fitting_prob_str += "\t\tfprintf(parameter_fileID, '%f\\n', parameter_opt(i));\n";
	sequential_fitting_prob_str += "\tend\n";
	sequential_fitting_prob_str += "\tfclose(parameter_fileID);\n";
	sequential_fitting_prob_str += "end\n";

	sequential_fitting_prob_str += "function best_parameter = sequential_fitting_for_one_order(model_order, N)\n";
	sequential_fitting_prob_str += "\toptions = optimset('MaxIter', 1000, 'MaxFunEvals', 30000);\n";
	sequential_fitting_prob_str += "\tbest_parameter = zeros(1," + Int2Str(combined_parameter_name_map->size()) + ");\n";
	sequential_fitting_prob_str += "\tfor model_id = 1:length(model_order)\n";
	sequential_fitting_prob_str += "\t\tswitch (model_id)\n";
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		string case_fitting_str = "";
		case_fitting_str += "case model_order(model_id)\n";
		map<string,int> non_fixed_parameter_name_map = model_it->extract_parameter_name_map(model_type);
		StringList parameter_name_list;
		ValueList lb, ub;
		get_parameter_value_bound(non_fixed_parameter_name_map, parameter_name_list, lb, ub);
		// lb
		case_fitting_str += "\tlb = [" + Float2Str(lb[0]);
		for (unsigned int i = 1; i < lb.size(); i++)
			case_fitting_str += " " + Float2Str(lb[i]);
		case_fitting_str += "];\n";
		// ub
		case_fitting_str += "\tub = [" + Float2Str(ub[0]);
		for (unsigned int i = 1; i < ub.size(); i++)
			case_fitting_str += " " + Float2Str(ub[i]);
		case_fitting_str += "];\n";
		for (unsigned int para_id = 0; para_id < parameter_name_list.size(); para_id++) {
			int combined_para_id = combined_parameter_name_map->at(parameter_name_list[para_id]);
			case_fitting_str += "\tif (best_parameter("+ Int2Str(combined_para_id + 1) + ") > 1e-10)\n";
			case_fitting_str += "\t\tlb(" + Int2Str(para_id + 1) + ") = best_parameter("+ Int2Str(combined_para_id + 1) + ") - 1e-10;\n";
			case_fitting_str += "\t\tub(" + Int2Str(para_id + 1) + ") = best_parameter("+ Int2Str(combined_para_id + 1) + ") + 1e-10;\n";
			case_fitting_str += "\tend\n" ;
		}
		// solver
		case_fitting_str += "\tparameter_num = length(lb);\n";
		case_fitting_str += "\tx0 = zeros(1, parameter_num);\n";
		case_fitting_str += "\tx_opt = zeros(1, parameter_num);\n";
		case_fitting_str += "\terror_opt = 1e6;\n";
		case_fitting_str += "\tfor i = 1:N\n";
		case_fitting_str += "\t\tfor j = 1:parameter_num\n";
		case_fitting_str += "\t\t\tx0(j) = lb(j) + (ub(j) - lb(j))*rand;\n";
		case_fitting_str += "\t\tend\n";
		string function_name = correct_to_a_var_name(model_it->get_id());
		case_fitting_str += "\t\t[x_min, error_min] = lsqnonlin(@simf_" + function_name + "_error, x0, lb, ub, options);\n";
		case_fitting_str += "\t\tif (error_min < error_opt)\n";
		case_fitting_str += "\t\t\terror_opt = error_min;\n";
		case_fitting_str += "\t\t\tx_opt = x_min;\n";
		case_fitting_str += "\t\tend\n";
		case_fitting_str += "\tend\n";
		for (unsigned int para_id = 0; para_id < parameter_name_list.size(); para_id++) {
			int combined_para_id = combined_parameter_name_map->at(parameter_name_list[para_id]);
			case_fitting_str += "\tbest_parameter(" + Int2Str(combined_para_id + 1) + ") = x_opt(" + Int2Str(para_id + 1) + ");\n" ;
		}
		add_tab_to_each_sentence(case_fitting_str);
		add_tab_to_each_sentence(case_fitting_str);
		add_tab_to_each_sentence(case_fitting_str);
		sequential_fitting_prob_str += case_fitting_str;
	}
	sequential_fitting_prob_str += "\t\tend\n";
	sequential_fitting_prob_str += "\tend\n";
	sequential_fitting_prob_str += "end\n";
	return sequential_fitting_prob_str;
}

string PAMLib::Export_MATLAB_code_for_sequential_fitting_LOOCV(vector<ModelAndExpData>* model_and_exp_data_list, IdMap* combined_parameter_name_map, int model_type) {
	IdMap parameter_appearance_count;
	int total_data_point_number = 0;
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list->begin(); model_it != model_and_exp_data_list->end(); model_it++) {
		total_data_point_number += model_it->exp_dataset.output_measurement.size();
		IdMap parameter_name_map = model_it->extract_parameter_name_map(model_type);
		for (IdMap::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++)
			if (parameter_appearance_count.find(it->first) == parameter_appearance_count.end())
					parameter_appearance_count[it->first] = 1;
				else
					parameter_appearance_count[it->first]++;
	}

	StringList parameter_name_list;
	ValueList lb, ub;
	get_parameter_value_bound(*combined_parameter_name_map, parameter_name_list, lb, ub);
	string seq_LOOCV_str =  "function seq_LOOCV(N1, N2)\n";
	seq_LOOCV_str += "\tLOOCV_seq_data_fileID = fopen('Plot_data/LOOCV_sim_data_" + get_math_model_type_name(model_type) + ".dat','w');\n";
	//string sim_LOOCV_error_str = "";
	string seq_LOOCV_fitting_str = "";
	StringList LOO_model_name_list;

	for (unsigned int LOO_model_id = 0; LOO_model_id < model_and_exp_data_list->size(); LOO_model_id++) {
		bool is_CV = true;
		IdMap parameter_name_map = model_and_exp_data_list->at(LOO_model_id).extract_parameter_name_map(model_type);
		for (map<string,int>::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++)
			if (parameter_appearance_count[it->first] == 1) {
				is_CV = false;
				break;
			}
		if (is_CV) {
			string function_name = correct_to_a_var_name(model_and_exp_data_list->at(LOO_model_id).get_id());
			LOO_model_name_list.push_back(function_name);
			// Fitting
			seq_LOOCV_fitting_str += "\terror_opt = 1e6;\n";
			seq_LOOCV_fitting_str += "\tx_opt = zeros(1," + Int2Str(combined_parameter_name_map->size()) + ");\n";
			seq_LOOCV_fitting_str += "\tfor i = 1:N1\n";
			seq_LOOCV_fitting_str += "\t\tmodel_order = LOOCV_randperm(" + Int2Str(model_and_exp_data_list->size()) + "," + Int2Str(LOO_model_id) + ");\n";
			seq_LOOCV_fitting_str += "\t\tparameter = sequential_fitting_for_one_order(model_order, N2);\n";
			seq_LOOCV_fitting_str += "\t\terror = norm(simultaneous_fitting_error(parameter));\n";
			seq_LOOCV_fitting_str += "\t\tif (error < error_opt)\n";
			seq_LOOCV_fitting_str += "\t\t\terror_opt = error;\n";
			seq_LOOCV_fitting_str += "\t\t\tx_opt = parameter;\n";
			seq_LOOCV_fitting_str += "\t\tend\n";
			seq_LOOCV_fitting_str += "\tend\n";
			// Validation
			if (model_and_exp_data_list->at(LOO_model_id).kinetic_model.GetInputList().size() == 0) {
				seq_LOOCV_fitting_str += "\tdesired_output = " + Float2Str(model_and_exp_data_list->at(LOO_model_id).exp_dataset.output_measurement[0][0].mean) + ";\n";
				seq_LOOCV_fitting_str += "\toutput = simf_" + function_name + "(x_opt);\n";
			}
			else if (model_and_exp_data_list->at(LOO_model_id).kinetic_model.GetInputList().size() == 1) {
				int input_node_id = model_and_exp_data_list->at(LOO_model_id).kinetic_model.GetInputList()[0];
				string signal_name = model_and_exp_data_list->at(LOO_model_id).kinetic_model.GetNodeInfo(input_node_id)->get_var_name();
				seq_LOOCV_fitting_str += "\t" + signal_name + " = [" + Float2Str(model_and_exp_data_list->at(LOO_model_id).exp_dataset.input_supplement[0][0]);
				for (unsigned int input_row = 1; input_row < model_and_exp_data_list->at(LOO_model_id).exp_dataset.input_supplement.size(); input_row++)
					seq_LOOCV_fitting_str += " " + Float2Str(model_and_exp_data_list->at(LOO_model_id).exp_dataset.input_supplement[input_row][0]);
				seq_LOOCV_fitting_str += "];\n";
				// Write the output list
				seq_LOOCV_fitting_str += "\tdesired_output = [" + Float2Str(model_and_exp_data_list->at(LOO_model_id).exp_dataset.output_measurement[0][0].mean);
				for (unsigned int output_row = 1; output_row < model_and_exp_data_list->at(LOO_model_id).exp_dataset.output_measurement.size(); output_row++)
					seq_LOOCV_fitting_str += " " + Float2Str(model_and_exp_data_list->at(LOO_model_id).exp_dataset.output_measurement[output_row][0].mean);
				seq_LOOCV_fitting_str += "];\n";
				seq_LOOCV_fitting_str += "\toutput = simf_" + function_name + "(x_opt, " + signal_name + ");\n";
			}
			seq_LOOCV_fitting_str += "\tfprintf(LOOCV_seq_data_fileID, '#" + function_name + "\\n');\n";
			seq_LOOCV_fitting_str += "\tfor i = 1:length(output)\n";
			seq_LOOCV_fitting_str += "\t\tfprintf(LOOCV_seq_data_fileID, '%f\\t%f\\n', desired_output(i), output(i));\n";
			seq_LOOCV_fitting_str += "\tend\n";
		}
	}
	seq_LOOCV_str += "\tLOO_model_name_list = {'" + LOO_model_name_list[0] + "'";
	for (unsigned int i = 1; i < LOO_model_name_list.size(); i++)
		seq_LOOCV_str += ", '" + LOO_model_name_list[i] + "'";
	seq_LOOCV_str += "};\n";
	seq_LOOCV_str += "\tLOOCV_error = zeros(" + Int2Str(LOO_model_name_list.size()) + ",1);\n";
	seq_LOOCV_str += seq_LOOCV_fitting_str;
	seq_LOOCV_str += "\tLOOCV_fileID = fopen('Plot_data/LOOCV_seq.dat','w');\n";
	seq_LOOCV_str += "\tfor i = 1:length(LOO_model_name_list)\n";
	seq_LOOCV_str += "\t\tfprintf(LOOCV_fileID, '%s\\t', LOO_model_name_list{i});\n";
	seq_LOOCV_str += "\t\tfprintf(LOOCV_fileID, '%f\\n', LOOCV_error(i));\n";
	seq_LOOCV_str += "\tend\n";
	seq_LOOCV_str += "\tfclose(LOOCV_fileID);\n";
	seq_LOOCV_str += "end\n";
	return seq_LOOCV_str;
}


void PAMLib::Export_MATLAB_code(string filename, int model_type) {
	vector<ModelAndExpData> model_and_exp_data_list;
	for (vector<Publication>::iterator paper = publication_list.begin(); paper != publication_list.end(); paper++) {
		paper->reference.print(&std::cout);
		for (unsigned int experiment_id = 0; experiment_id < paper->experiment_list.size(); experiment_id++) {
			ModelAndExpData model_and_exp_data;
			model_and_exp_data.paper_ref = paper->reference;
			model_and_exp_data.experiment_id = paper->experiment_list[experiment_id].reference;
			convert_experiment_to_model_and_data(paper->experiment_list[experiment_id], model_and_exp_data);
			//model_and_exp_data.print(&std::cout);
			model_and_exp_data_list.push_back(model_and_exp_data);
		}
	}
	if (model_and_exp_data_list.empty()) {
		sbrome_error_print("Empty model set !!!\n");
		return;
	}
	// Combine all parameters
	IdMap combined_parameter_name_map;
	int parameter_num = 0;
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list.begin(); model_it != model_and_exp_data_list.end(); model_it++) {
		IdMap parameter_name_map = model_it->extract_parameter_name_map(model_type);
		for (IdMap::iterator it = parameter_name_map.begin(); it != parameter_name_map.end(); it++)
			if (combined_parameter_name_map.find(it->first) == combined_parameter_name_map.end()) {
				combined_parameter_name_map[it->first] = parameter_num;
				parameter_num++;
			}
	}
	// Write the file
	ofstream MATLAB_file(filename);
	MATLAB_file << "function real_benchmark_fitting" << endl;
	MATLAB_file << "\t%sim_x_opt = simultaneous_fitting(1,1000,50000);" << endl;
	MATLAB_file << "\t%print_plotting_sim_data(sim_x_opt);" << endl;
	MATLAB_file << "\t%sim_LOOCV(10,100,5000);" << endl;
	MATLAB_file << "\tseq_x_opt = sequential_fitting(10, 10);" << endl;
	MATLAB_file << "\tseq_LOOCV(10,10);" << endl;
	MATLAB_file << "end" << endl;

	// simulation functions
	for (vector<ModelAndExpData>::iterator model_it = model_and_exp_data_list.begin(); model_it != model_and_exp_data_list.end(); model_it++)
		MATLAB_file << model_it->export_MATLAB_simulation_code(combined_parameter_name_map, model_type) << endl;
	MATLAB_file << Export_MATLAB_code_for_simultaneous_fitting(&model_and_exp_data_list, &combined_parameter_name_map, model_type) << endl;
	MATLAB_file << Export_MATLAB_code_for_printing_simultaneous_fitting_result(&model_and_exp_data_list, &combined_parameter_name_map, model_type) << endl;
	MATLAB_file << Export_MATLAB_code_for_simultaneous_fitting_LOOCV(&model_and_exp_data_list, &combined_parameter_name_map, model_type) << endl;
	MATLAB_file << Export_MATLAB_code_for_sequential_fitting(&model_and_exp_data_list, &combined_parameter_name_map, model_type) << endl;
	MATLAB_file << Export_MATLAB_code_for_sequential_fitting_LOOCV(&model_and_exp_data_list, &combined_parameter_name_map, model_type) << endl;
	MATLAB_file.close();
}
