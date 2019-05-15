/*
 * functional_spec.h
 *
 *  Created on: Apr 20, 2012
 *      Author: linh, UC Davis
 */

#ifndef FUNCTIONAL_SPEC_H_
#define FUNCTIONAL_SPEC_H_

#include "utilities.h"
#include "physical_spec.h"

class Environment: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string medium;
	double temperature;
	vector<StringPair> supplemented_chemical_name_list;
};

class Measurement: public JSONEntityInferace {
public:
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	double growth_rate;
	vector<StringPair> molecular_species_list;
};

struct ExperimentalDataPoint {
	ExperimentalDataPoint(double = UNKNOWN);
	ExperimentalDataPoint(ValueList);
	void print(ostream* f);
	ValueList replicated_value_list;
	double mean;
	double std_dev;		// Standard deviation
	double std_error;	// Standard error
};
typedef vector<ExperimentalDataPoint> ExperimentalDataList;
typedef vector<ExperimentalDataList> ExperimentalDataMatrix;

class MeasurementData: public JSONEntityInferace {
public:
	//MeasurementData();
	//MeasurementData(vector<double>, ValueMatrix, ExperimentalDataMatrix);
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	vector<double> time_point_list;
	ValueMatrix input_supplement;				// Each row is a combination of supplemented concentrations
	ExperimentalDataMatrix output_measurement;	// Each row corresponds to a row of the input matrix, each tuple of n column (n = number of outputs) is for a time point
};

class Experiment: public JSONEntityInferace {
public:
	//Experiment();
	//Experiment(string, Genotype, Environment, Measurement, MeasurementData);
	void print(ostream*);
	void read(boost::property_tree::ptree*);
	string reference;
	Genotype genotype;
	Environment environment;
	Measurement measurement;
	MeasurementData dataset;
};

#endif /* FUNCTIONAL_SPEC_H_ */
