/*
 * biological_spec.h
 *
 *  Created on: May 2, 2012
 *      Author: linh
 */

#ifndef BIOLOGICAL_SPEC_H_
#define BIOLOGICAL_SPEC_H_

#include "physical_spec.h"

class MolecularSpecies: public JSONEntityInferace {
public:
	static const int UNKNOWN_MOLECULE 			= -1;
	static const int LIGAND 					= 0;
	static const int RNA						= 1;
	static const int mRNA 						= 2;
	static const int tRNA						= 3;
	static const int PROTEIN 					= 4;
	static const int LIGAND_COMPLEX 			= 5;
	static const int RNA_COMPLEX 				= 6;
	static const int PROTEIN_COMPLEX 			= 7;
	static const int LIGAND_PROTEIN_COMPLEX 	= 8;
	static int PartType2MoleculeType(int);
	// Converter
	static int Name2Id(string);
	static string Id2Name(int);

	MolecularSpecies(string = "", string = "", string = "");
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string name;
	string decoded_dna_component;
	string type;
private:
	static IdMap initialize_name_map();
	static const IdMap name_map;
};

class Reaction: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	vector<string> reactant_list;
	vector<string> product_list;
	string mechanism_type;
};

struct BioNetNodeInfo {
	void print(ostream* f) const;
	string get_var_name();
	int type;
	string name;
	int copy_number;
	vector<BioNetNodeInfo> aux_node_list;	// For complex molecules, pool node or IO node
};

BioNetNodeInfo create_a_bio_net_node_info(int, string, int);
BioNetNodeInfo create_a_bio_net_node_info(int, string, int, int, string, int);
BioNetNodeInfo create_a_bio_net_node_info(int, string, int, BioNetNodeInfo);

struct BioNetEdgeInfo {
	static const int ACTIVATORY			= 1;
	static const int INHIBITORY			= -1;
	void print(ostream* f) const;
	int type;
};
BioNetEdgeInfo create_a_bio_net_edge_info(int);

struct BioNetNode {
	BioNetNodeInfo bio_net_node_info;
	map<int, BioNetEdgeInfo> edge_list;
};

class BiologicalNetwork {
public:
	BiologicalNetwork();
	void print(ostream* f);
	// Insert
	int InsertNode(BioNetNodeInfo);
	void InsertInput(int);
	void InsertOutput(int);
	void InsertEdge(int, int, BioNetEdgeInfo);
	void InsertEdge(string, string, BioNetEdgeInfo);
	// Delete
	void DeleteNode(int);
	void DeleteEdge(int, int);
	void ClearInput();
	void ClearOutput();
	// Check
	bool IsNodeNameAvailable(string);
	// Get
	int GetFirstNodeId(string);
	BioNetNodeInfo* GetNodeInfo(int);
	BioNetNodeInfo* GetNodeInfo(string);
	BioNetEdgeInfo* GetEdgeInfo(int, int);
	BioNetEdgeInfo* GetEdgeInfo(string, string);
	IdList GetNodeList();
	IdList GetInNodeList(int);
    IdList GetOutNodeList(int);
    IdList GetInputList() const;
    IdList GetOutputList() const;
    // Count
    int NodeCount() const;
    int InEdgeCount(int);
    int OutEdgeCount(int) const;

    // Functions on a network
	bool Topo_Sort(IdList* topo_order_list = NULL);		// return false if there is a loop
	void AddPoolNode();									// Add pool nodes to merge duplicated nodes
	void RemoveIrrelevantNode();						// Remove all irrelevant nodes that have no interaction with inputs and outputs
protected:
	void print_debug(ostream* f);
	int available_node_id;
	map<int,BioNetNode> node_list;
	map<string,IdList> name_to_node_id_list;
	IdList input_list, output_list, node_id_list;
};

class KineticModel: public BiologicalNetwork {
public:
	IdMap ExtractParameterMap(string experiment_name, int model_type);
	static const int _MODEL_OPTION_HILL = 1;
	static const int _MODEL_OPTION_MOON = 2;
	static const int _MODEL_OPTION_ELLIS = 3;
private:
	void add_parameter(string, IdMap&);
};

#endif /* BIOLOGICAL_SPEC_H_ */
