/*
 * functional_spec.cpp
 *
 *  Created on: Apr 20, 2012
 *      Author: linh, UC Davis
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include "functional_spec.h"

void Environment::print(ostream* f) {
	*f << "Environment: " << endl;
	*f << "Medium: " << medium << "\t T = " << temperature << endl;
	for (unsigned int i = 0; i < supplemented_chemical_name_list.size(); i++)
		*f << supplemented_chemical_name_list[i].first << " " << supplemented_chemical_name_list[i].second << "\t";
	*f << endl;
}

void Environment::read(boost::property_tree::ptree* pt) {
	medium = pt->get<string>("medium");
	temperature = pt->get<double>("temperature");
	supplemented_chemical_name_list.clear();
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("supplement", empty_ptree)) {
		boost::property_tree::ptree::iterator i1 = v.second.begin();
		boost::property_tree::ptree::iterator i2 = v.second.begin();
		i2++;
		supplemented_chemical_name_list.push_back(make_pair(i1->second.get_value<string>(), i2->second.get_value<string>()));
	}
}

void Measurement::print(ostream* f) {
	*f << "Phenotype: " << endl;
	for (unsigned int i = 0; i < molecular_species_list.size(); i++)
		*f << molecular_species_list[i].first << "\t" << molecular_species_list[i].second << endl;
}

void Measurement::read(boost::property_tree::ptree* pt) {
	molecular_species_list.clear();
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("molecular_species", empty_ptree)) {
		boost::property_tree::ptree::iterator i1 = v.second.begin();
		boost::property_tree::ptree::iterator i2 = v.second.begin();
		i2++;
		molecular_species_list.push_back(make_pair(i1->second.get_value<string>(), i2->second.get_value<string>()));
	}
}

ExperimentalDataPoint::ExperimentalDataPoint(double x) {
	if (x < 0)
		mean = std_dev= std_error = UNKNOWN;
	else {
		mean = x;
		std_dev = std_error = 1;
	}
}

ExperimentalDataPoint::ExperimentalDataPoint (ValueList value_list) {
	summarize(value_list, mean, std_dev, std_error);
}

void ExperimentalDataPoint::print(ostream* f) {
	if (std_dev >= 0)
		*f << mean << "  +/-  " << std_dev;
	else
		*f << mean;
}

void MeasurementData::print(ostream* f) {
	for (unsigned int i = 0; i < time_point_list.size(); i++)
		*f << time_point_list[i] << "\t";
	cout << endl;
	for (unsigned int i = 0; i < input_supplement.size(); i++) {
		for (unsigned int j = 0; j < input_supplement[i].size(); j++)
			*f << input_supplement[i][j] << "\t";
		for (unsigned int j = 0; j < output_measurement[i].size(); j++) {
			output_measurement[i][j].print(f);
			*f << "\t";
		}
		*f << endl;
	}
}

void MeasurementData::read(boost::property_tree::ptree* pt) {
	// Only for the case for one input + one output with no std
	time_point_list.clear();
	boost::property_tree::ptree empty_ptree;
	BOOST_FOREACH (boost::property_tree::ptree::value_type &v, pt->get_child("input_output", empty_ptree)) {
		if (v.second.size() == 1) {
			vector<double> x;
			input_supplement.push_back(x);
			ExperimentalDataPoint tmp(v.second.begin()->second.get_value<double>());
			vector<ExperimentalDataPoint> y(1, tmp);
			output_measurement.push_back(y);
		}
		else if (v.second.size() == 2) {
			boost::property_tree::ptree::iterator it = v.second.begin();
			vector<double> x(1, it->second.get_value<double>());
			input_supplement.push_back(x);
			ExperimentalDataPoint tmp((++it)->second.get_value<double>());
			vector<ExperimentalDataPoint> y(1, tmp);
			output_measurement.push_back(y);
		}
		else
			sbrome_error_print("Each data record has only 1 or 2 fields");
	}
}

void Experiment::print(ostream* f) {
	*f << "Reference: " << reference << endl;
	genotype.print(f);
	environment.print(f);
	measurement.print(f);
	*f << "Dataset: " << endl;
	dataset.print(f);
}

void Experiment::read(boost::property_tree::ptree* pt) {
	reference = pt->get<string>("reference", "");
	genotype.read(&(pt->get_child("genotype")));
	environment.read(&(pt->get_child("environment")));
	measurement.read(&(pt->get_child("phenotype")));
	dataset.read(&(pt->get_child("dataset")));
}


