/*
 * biological_spec.cpp
 *
 *  Created on: May 2, 2012
 *      Author: linh
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include "biological_spec.h"

MolecularSpecies::MolecularSpecies(string s1, string s2, string s3) {
	name = s1;
	type = s2;
	decoded_dna_component = s3;
}

void MolecularSpecies::print(ostream* f) {
	*f << "Molecular species: " << name << "\t from DNA: " << decoded_dna_component << endl;
}

void MolecularSpecies::read (boost::property_tree::ptree* pt) {
	name = pt->get<string>("name");
	type = pt->get<string>("type","");
	decoded_dna_component = pt->get<string>("encoded dna","");
}

int MolecularSpecies::PartType2MoleculeType (int part_type) {
	switch (part_type) {
		case DnaComponent::PROMOTER:
			return mRNA;
		case DnaComponent::CDS:
			return PROTEIN;
		default:
			sbrome_error_print("There is no part type " + Int2Str(part_type) + "!\n");
			return UNKNOWN_MOLECULE;
	}
}

IdMap MolecularSpecies::initialize_name_map() {
	IdMap name_map_tmp;
	name_map_tmp ["Ligand"] 					= LIGAND;
	name_map_tmp ["Protein"] 					= PROTEIN;
	name_map_tmp ["mRNA"] 						= mRNA;
	name_map_tmp ["LigandProteinComplex"]		= LIGAND_PROTEIN_COMPLEX;
	return name_map_tmp;
}

const IdMap MolecularSpecies::name_map = MolecularSpecies::initialize_name_map();

int MolecularSpecies::Name2Id(string name) {
	return name_map.at(name);
}

string MolecularSpecies::Id2Name(int id) {
	for (IdMap::const_iterator it = name_map.begin(); it != name_map.end(); it++)
		if (it->second == id)
			return it->first;
	sbrome_error_print("There is no molecular species type Id = " + Int2Str(id) + "!\n");
	return "";
}

void Reaction::print(ostream* f) {
	*f << "Reaction with mechanism: " << mechanism_type << endl;
	*f << "Reactant list: ";
	sbrome_print(f, reactant_list);
	*f << "Product list: ";
	sbrome_print(f, product_list);
}

void Reaction::read(boost::property_tree::ptree* pt) {
	mechanism_type = pt->get<string>("mechanism");
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("reactant", empty_ptree)) {
		reactant_list.push_back(v.second.get_value<string>());
	}
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("product", empty_ptree)) {
		product_list.push_back(v.second.get_value<string>());
	}
}

BioNetEdgeInfo create_a_bio_net_edge_info(int edge_type) {
	BioNetEdgeInfo tmp;
	tmp.type = edge_type;
	return tmp;
}

void BioNetNodeInfo::print(ostream* f) const {
	*f << name << "\t" << MolecularSpecies::Id2Name(type) << " (" << copy_number << ")";
	if (!aux_node_list.empty()) {
		if (type == MolecularSpecies::PROTEIN)
			*f << " RBS: " << aux_node_list[0].name;
	}
}

string BioNetNodeInfo::get_var_name() {
	return correct_to_a_var_name(name);
}

BioNetNodeInfo create_a_bio_net_node_info(int type, string name, int n) {
	BioNetNodeInfo tmp;
	tmp.type = type;
	tmp.name = name;
	tmp.copy_number = n;
	return tmp;
}
BioNetNodeInfo create_a_bio_net_node_info(int type, string name, int n, int aux_type, string aux_name, int aux_n) {
	return create_a_bio_net_node_info(type, name, n, create_a_bio_net_node_info(aux_type, aux_name, aux_n));
}
BioNetNodeInfo create_a_bio_net_node_info(int type, string name, int n, BioNetNodeInfo aux_node_info) {
	BioNetNodeInfo tmp = create_a_bio_net_node_info(type, name, n);
	tmp.aux_node_list.push_back(aux_node_info);
	return tmp;
}

BiologicalNetwork::BiologicalNetwork() {
	available_node_id = 1;
}

void BiologicalNetwork::print(ostream* f) {
	for (map<int,BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++) {
		*f << "(" << it->first << ")\t";
		it->second.bio_net_node_info.print(f);
		*f << endl;
	}
	for (map<int,BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++)
		for (map<int,BioNetEdgeInfo>::iterator e_it = it->second.edge_list.begin(); e_it != it->second.edge_list.end(); e_it++)
			*f << it->first << "\t" << e_it->first << "\t" << ((e_it->second.type == BioNetEdgeInfo::ACTIVATORY)?"ACT":"REP") << endl;
	*f << "Inputs: ";
	for (unsigned int i = 0; i < input_list.size(); i++)
		*f << input_list[i] << " ";
	*f << endl << "Outputs: ";
	for (unsigned int i = 0; i < output_list.size(); i++)
		*f << output_list[i] << " ";
	*f << endl;
}

void BiologicalNetwork::print_debug(ostream* f) {
	print(f);
	for (map<string,IdList>::iterator it = name_to_node_id_list.begin(); it != name_to_node_id_list.end(); it++) {
		*f << it->first << "\t";
		sbrome_print(f, it->second);
	}
	*f << endl << "Node id list: ";
	sbrome_print(f, node_id_list);
	*f << endl;
}

int BiologicalNetwork::InsertNode(BioNetNodeInfo node_info) {
	BioNetNode tmp;
	tmp.bio_net_node_info = node_info;
	node_list[available_node_id] = tmp;
	if (name_to_node_id_list.find(node_info.name) == name_to_node_id_list.end()) {
		IdList list_tmp(1, available_node_id);
		name_to_node_id_list[node_info.name] = list_tmp;
	}
	else
		name_to_node_id_list[node_info.name].push_back(available_node_id);
	node_id_list.push_back(available_node_id);
	return available_node_id++;
}

void BiologicalNetwork::InsertInput(int node_id) {
	input_list.push_back(node_id);
}

void BiologicalNetwork::InsertOutput(int node_id) {
	output_list.push_back(node_id);
}

void BiologicalNetwork::InsertEdge(int src_id, int dest_id, BioNetEdgeInfo edge_info) {
	//cout << "Add edge " << src_id << "\t" << dest_id << endl;
	if (node_list[src_id].edge_list.count(dest_id) != 0)
		sbrome_error_print("Insert an edge that is already in the graph!\n");
	else
		node_list[src_id].edge_list[dest_id] = edge_info;
}

void BiologicalNetwork::InsertEdge(string src_name, string dest_name, BioNetEdgeInfo edge_info) {
	InsertEdge(GetFirstNodeId(src_name), GetFirstNodeId(dest_name), edge_info);
}

void BiologicalNetwork::DeleteNode(int node_id) {
	//cout << "================================" << endl;
	//print_debug(&std::cout);
	//cout << "Delete " << node_id << endl;
	//print_debug(&std::cout);
	//cout << "================================" << endl;
	string node_name = this->GetNodeInfo(node_id)->name;
	// Delete node
	node_list.erase(node_id);
	// Delete from node_id list, input and output
	delete_first_pos(node_id_list, node_id);
	delete_first_pos(input_list, node_id);
	delete_first_pos(output_list, node_id);
	// Delete all edges
	for (map<int, BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++)
		it->second.edge_list.erase(node_id);
	// Delete from name map
	delete_first_pos(name_to_node_id_list[node_name], node_id);
	if (name_to_node_id_list[node_name].empty())
		name_to_node_id_list.erase(node_name);
}

void BiologicalNetwork::DeleteEdge(int src, int dest) {
	if (node_list[src].edge_list.count(dest) == 0)
		sbrome_error_print("Delete an edge from " + Int2Str(src) + " to " + Int2Str(dest) + " that does not exist!\n");
	node_list[src].edge_list.erase(dest);
}

void BiologicalNetwork::ClearInput() {
	input_list.clear();
}

void BiologicalNetwork::ClearOutput() {
	output_list.clear();
}

bool BiologicalNetwork::IsNodeNameAvailable(string node_name) {
	return (name_to_node_id_list.find(node_name) != name_to_node_id_list.end());
}

int BiologicalNetwork::GetFirstNodeId(string node_name) {
	if (name_to_node_id_list.find(node_name) == name_to_node_id_list.end())
		sbrome_error_print("there is no node with name = " + node_name + " in the graph!\n");
	else {
		if (name_to_node_id_list[node_name].size() > 1)
			sbrome_error_print("there is more than one node with name = " + node_name + " in the graph!\n");
		else if (name_to_node_id_list[node_name].empty())
			sbrome_error_print("there is the name in the graph but the number of nodes is zero!\n");
		else
			return name_to_node_id_list[node_name][0];
	}
	return UNKNOWN;
}

BioNetNodeInfo* BiologicalNetwork::GetNodeInfo(int node_id) {
	return &(node_list[node_id].bio_net_node_info);
}

BioNetNodeInfo* BiologicalNetwork::GetNodeInfo(string node_name) {
	return &(node_list[GetFirstNodeId(node_name)].bio_net_node_info);
}

BioNetEdgeInfo* BiologicalNetwork::GetEdgeInfo(int src_id, int dest_id) {
	if (node_list.count(src_id) == 0)
		sbrome_error_print("Access an edge from an non-exist node " + Int2Str(src_id) + "!\n");
	if (node_list[src_id].edge_list.count(dest_id) == 0)
		sbrome_error_print("Access a non-exist edge from node " + Int2Str(src_id) + " to node " + Int2Str(dest_id) + "!\n");
	return &(node_list[src_id].edge_list[dest_id]);
}

BioNetEdgeInfo* BiologicalNetwork::GetEdgeInfo(string src_name, string dest_name) {
	return GetEdgeInfo(GetFirstNodeId(src_name), GetFirstNodeId(dest_name));
}

IdList BiologicalNetwork::GetNodeList() {
	return node_id_list;
}

IdList BiologicalNetwork::GetInNodeList(int node_id) {
	IdList tmp;
	for (map<int,BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++)
		if (it->second.edge_list.count(node_id) > 0)
			tmp.push_back(it->first);
	return tmp;
}

IdList BiologicalNetwork::GetOutNodeList(int node_id) {
	IdList tmp;
	for (map<int,BioNetEdgeInfo>::iterator it = node_list[node_id].edge_list.begin(); it != node_list[node_id].edge_list.end(); it++)
		tmp.push_back(it->first);
	return tmp;
}


IdList BiologicalNetwork::GetInputList() const {
	return input_list;
}

IdList BiologicalNetwork::GetOutputList() const {
	return output_list;
}

int BiologicalNetwork::NodeCount() const {
	return node_list.size();
}

int BiologicalNetwork::InEdgeCount(int node_id) {
	int n = 0;
	for (map<int,BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++)
		n += it->second.edge_list.count(node_id);
	return n;
}

int BiologicalNetwork::OutEdgeCount(int node_id) const {
	return node_list.at(node_id).edge_list.size();
}

void BiologicalNetwork::AddPoolNode() {

}

bool BiologicalNetwork::Topo_Sort(IdList* topo_order_list) {
	bool is_null_input = false;
	if (topo_order_list == NULL) {
		topo_order_list = new IdList;
		is_null_input = true;
	}
	bool is_cyclic;
	unsigned int number_of_nodes = node_id_list.size();
	map<int,bool> marked;
	map<int,int> in_deg;		// in-degree list of nodes
	topo_order_list->resize(number_of_nodes);
	for (unsigned int i = 0; i < number_of_nodes; i++) {
		in_deg[node_id_list[i]] = this->InEdgeCount(node_id_list[i]);
		marked[node_id_list[i]] = false;
	}
	int current_index = 0;
	do {
		is_cyclic = true;
		for (unsigned int i = 0; i < number_of_nodes; i++)
			if (!marked[node_id_list[i]] && in_deg[node_id_list[i]] == 0) {
				int out_num = this->OutEdgeCount(node_id_list[i]);
				for (int j = 0; j < out_num; j++)
					in_deg[this->GetOutNodeList(node_id_list[i])[j]]--;
				marked[node_id_list[i]] = true;
				topo_order_list->at(current_index) = node_id_list[i];
				current_index++;
				is_cyclic = false;
			}
		if (is_cyclic)
			break;
	}
	while (current_index < number_of_nodes);
	if (is_null_input)
		delete topo_order_list;
	//if (is_cyclic) {
	//	cout << "---------------------------" << endl;
	//	this->TestPrint(&std::cout);
	//	cout << "---------------------------" << endl;
	//}
	return (not is_cyclic);
}

void BiologicalNetwork::RemoveIrrelevantNode() {
	// Remove all abundant nodes
	bool continue_searching;
	map<int,bool> is_deleted_node;
	for (unsigned int i = 0; i < node_id_list.size(); i++)
		is_deleted_node[node_id_list[i]] = true;
	for (unsigned int i = 0; i < input_list.size(); i++)
		is_deleted_node[input_list[i]] = false;
	for (unsigned int i = 0; i < output_list.size(); i++)
		is_deleted_node[output_list[i]] = false;
	do {
		continue_searching = false;
		for (unsigned int i = 0; i < node_id_list.size(); i++) {
			if (!is_deleted_node[node_id_list[i]]) {
				//cout << "xxx " << node_id_list[i] << endl;
				IdList out_list = this->GetOutNodeList(node_id_list[i]);
				//std::cout << "Out "<< endl;
				//sbrome_print(&std::cout, out_list);
				//std::cout << endl;
				for (unsigned int j = 0; j < out_list.size(); j++)
					if (is_deleted_node[out_list[j]]) {
						is_deleted_node[out_list[j]] = false;
						continue_searching = true;
					}
				IdList in_list = this->GetInNodeList(node_id_list[i]);
				//std::cout << "In "<< endl;
				//sbrome_print(&std::cout, in_list);
				//std::cout << endl;
				for (unsigned int j = 0; j < in_list.size(); j++)
					if (is_deleted_node[in_list[j]]) {
						is_deleted_node[in_list[j]] = false;
						continue_searching = true;
					}
			}
		}
	}
	while (continue_searching);
	//for (unsigned int i = 0; i < node_id_list.size(); i++)
	//	cout << node_id_list[i] << "\t" << is_deleted_node[node_id_list[i]];
	//cout << endl;
	for (map<int,bool>::iterator it = is_deleted_node.begin(); it != is_deleted_node.end(); it++)
		if (it->second)
			DeleteNode(it->first);
}

void KineticModel::add_parameter(string parameter_name, IdMap& parameter_name_map) {
	if (parameter_name_map.find(parameter_name) == parameter_name_map.end())
		parameter_name_map[parameter_name] = parameter_name_map.size() - 1;
}

IdMap KineticModel::ExtractParameterMap(string experiment_name, int model_type) {
	IdMap parameter_name_map;
	for (map<int,BioNetNode>::iterator it = node_list.begin(); it != node_list.end(); it++) {
		string node_var_name = it->second.bio_net_node_info.get_var_name();
		switch (it->second.bio_net_node_info.type) {
			case MolecularSpecies::mRNA:
				if (InEdgeCount(it->first) == 0)
					add_parameter("alpha_" + node_var_name, parameter_name_map);
				else {
					bool is_hill_eq = false;
					if (model_type == KineticModel::_MODEL_OPTION_ELLIS) {
						if (node_var_name.compare("pLAC") == STR_EQ || node_var_name.compare("pLlacO_1") == STR_EQ || node_var_name.compare("placUV5") == STR_EQ || node_var_name.compare("pTAC") == STR_EQ
							|| node_var_name.compare("pTET") == STR_EQ || node_var_name.compare("pLtetO_1") == STR_EQ || node_var_name.compare("pTET_star_") == STR_EQ) {
							add_parameter("alpha_" + node_var_name, parameter_name_map);
							add_parameter("beta_" + node_var_name, parameter_name_map);
							add_parameter("KE1_" + node_var_name, parameter_name_map);
							add_parameter("KE2_" + node_var_name, parameter_name_map);
						}
						else if (node_var_name.compare("pBAD") == STR_EQ) {
							is_hill_eq = true;
						}
						else {
							sbrome_error_print("Ellis model does not include " + node_var_name);
						}
					}
					if (model_type == KineticModel::_MODEL_OPTION_MOON) {
						if (node_var_name.compare("pLAC") == STR_EQ || node_var_name.compare("pLlacO_1") == STR_EQ || node_var_name.compare("placUV5") == STR_EQ || node_var_name.compare("pTAC") == STR_EQ
							|| node_var_name.compare("pTET") == STR_EQ || node_var_name.compare("pLtetO_1") == STR_EQ || node_var_name.compare("pTET_star_") == STR_EQ) {
							add_parameter("alpha_" + node_var_name, parameter_name_map);
							add_parameter("KM1_" + node_var_name, parameter_name_map);
							add_parameter("KM2_" + node_var_name, parameter_name_map);
						}
						else if (node_var_name.compare("pBAD") == STR_EQ) {
							add_parameter("alpha_" + node_var_name, parameter_name_map);
							add_parameter("KM1_" + node_var_name, parameter_name_map);
							add_parameter("KM2_" + node_var_name, parameter_name_map);
							add_parameter("KM3_" + node_var_name, parameter_name_map);
						}
						else {
							sbrome_error_print("Moon model does not include " + node_var_name);
						}
					}
					if (model_type == KineticModel::_MODEL_OPTION_HILL || is_hill_eq) {
						add_parameter("alpha_" + node_var_name, parameter_name_map);
						add_parameter("beta_" + node_var_name, parameter_name_map);
						add_parameter("n_" + node_var_name, parameter_name_map);
						add_parameter("K_" + node_var_name, parameter_name_map);
					}
				}
				break;
			case MolecularSpecies::PROTEIN: {
				if (it->second.bio_net_node_info.aux_node_list.empty())
					sbrome_error_print("RBS in the uptream of CDS " + it->second.bio_net_node_info.name + " is not declared!\n");
				else {
					string RBS_name = it->second.bio_net_node_info.aux_node_list[0].get_var_name();
					add_parameter("ALPHA_" + RBS_name, parameter_name_map);
				}
				if (output_list[0] == it->first) {
					add_parameter("scale_" + node_var_name + "_" + correct_to_a_var_name(experiment_name), parameter_name_map);
				}
				break;
			}
			case MolecularSpecies::LIGAND:
				add_parameter("n_" + node_var_name, parameter_name_map);
				add_parameter("K_" + node_var_name, parameter_name_map);
				break;
		}
	}
	return parameter_name_map;
}
