/*
 * sbrome.cpp
 *
 *  Created on: Apr 8, 2012
 *      Author: linh, UC Davis
 */

#include <fstream>
#include "utilities.h"
#include "PAMLib.h"
#include "sbrome.h"

using namespace std;

void ACS_2016_Experiment() {
	PAMLib pam_lib("Database/ACS_2016");
	pam_lib.print(&std::cout, PAMLib::_PRINT_OPTION_SUMMARY);
	pam_lib.Export_MATLAB_code("MATLAB/ACS_2016/Hill_ACS_2016_benchmark_fitting.m", KineticModel::_MODEL_OPTION_HILL);
	//pam_lib.Export_MATLAB_code("MATLAB/ACS_2016/Ellis_ACS_2016_benchmark_fitting.m", KineticModel::_MODEL_OPTION_ELLIS);
	//pam_lib.Export_MATLAB_code("MATLAB/ACS_2016/Moon_ACS_2016_benchmark_fitting.m", KineticModel::_MODEL_OPTION_MOON);
}

void PAMLib_Experiment() {
	//combine_all_json_file("Database/CuratedLiterature", "MATLAB/JSON/PAMLib_curated_data.json");
	PAMLib pam_lib("Database/PAMLib");
	//ofstream debug_log_file;
	//debug_log_file.open("debug_log.txt");
	//pam_lib.print(&std::cout, ALL);
	pam_lib.Export_MATLAB_code("MATLAB/Hill_Real_benchmark_fitting.m", KineticModel::_MODEL_OPTION_HILL);
	cout << "FINISH!" << endl;
	//debug_log_file.close();
}

int main() {
	ACS_2016_Experiment();
	cout << "FINISH" << endl;
}
