/*
 * utilities.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: linh
 */

#include <fstream>
#include <boost/filesystem.hpp>
#include "utilities.h"

void sbrome_print(ostream* f, IdList id_list) {
	for (IdList::iterator it = id_list.begin(); it != id_list.end(); it++)
		*f << *it << "\t";
	*f << "\n";
}

void sbrome_print(ostream* f, StringList string_list) {
	for (StringList::iterator it = string_list.begin(); it != string_list.end(); it++)
		*f << *it << "\t";
	*f << "\n";
}
void sbrome_print(ostream* f, IdMap id_map) {
	for (IdMap::iterator it = id_map.begin(); it != id_map.end(); it++)
		*f << it->first << "\t" << it->second << endl;
}

void sbrome_error_print (string msg) {
	std::cerr << "SBROME ERROR: " << msg << endl;
}

string correct_to_a_var_name(string str) {
	string tmp = "";
	for (unsigned int i = 0; i < str.size(); i++) {
		switch (str[i]) {
		case '-':
			tmp += "_";
			break;
		case ':':
			tmp += "_";
			break;
		case '/':
			tmp += "_";
			break;
		case '.':
			tmp += "_";
			break;
		case ' ':
			tmp += "_";
			break;
		case '*':
			tmp += "_star_";
			break;
		case '+':
			tmp += "_plus_";
			break;
		default:
			tmp += str[i];
		}
	}
	if (tmp[0] >= '0' && tmp[0] <= '9')
		tmp = "C_" + tmp;
	return tmp;
}

void add_tab_to_each_sentence(string& paragraph) {
	string tmp = "\t";
	for (unsigned int i = 0; i < paragraph.size(); i++) {
		tmp += paragraph[i];
		if (paragraph[i] == '\n' && i != paragraph.size() - 1)
			tmp += "\t";
	}
	paragraph = tmp;
}

void split_string(string str, char c, vector<string> &word_list) {
	word_list.clear();
	int pre_index = 0;
	for (unsigned int i = 0; i < str.length(); i++)
		if (str[i] == c) {
			word_list.push_back(str.substr(pre_index, i - pre_index));
			pre_index = i + 1;
		}
	word_list.push_back(str.substr(pre_index, str.length() - pre_index));
}

string Int2Str(int num) {
	stringstream ss;
	ss << num;
	string tmp(ss.str());
	return tmp;
}

string Float2Str(double num) {
	num = round(num*1000)/1000;
	std::ostringstream strs;
	strs << num;
	return strs.str();
}

bool delete_first_pos(IdList& list, int x) {
	for (IdList::iterator it = list.begin(); it != list.end(); it++) {
		if (*it == x) {
			list.erase(it);
			return true;
		}
	}
	return false;
}

void frequency_counting(string name, IdMap& name_freq) {
	if (name_freq.find(name) == name_freq.end())
		name_freq[name] = 1;
	else
		name_freq[name]++;
}

void summarize(const ValueList& value_list, double &mean, double &std_dev, double &std_error) {
	if (value_list.empty())
		sbrome_error_print("Calculate mean & stdev for an empty dataset");
	if (value_list.size() > 1) {
		mean = 0;
		for (unsigned int i = 0; i < value_list.size(); i++)
			mean += value_list[i];
		mean = mean/value_list.size();
		std_dev = 0;
		for (unsigned int i = 0; i < value_list.size(); i++)
			std_dev += (value_list[i] - mean)*(value_list[i] - mean);
		std_dev = sqrt(std_dev/(value_list.size() - 1));
		std_error = std_dev/sqrt(value_list.size());
	}
	else {
		mean = value_list[0];
		std_dev = std_error = 0;
	}
}

void combine_all_json_file(string folder_path_str, string combined_filename) {
	std::ofstream combined_file(combined_filename);
	write_all_json_file(folder_path_str, &combined_file);
	combined_file.close();
}

void write_all_json_file(string folder_path_str, ostream* f) {
	boost::filesystem::path folder_path(folder_path_str);
	if (!boost::filesystem::exists(folder_path))
		sbrome_error_print("There is no folder " + folder_path.string() + "\n");
	boost::filesystem::directory_iterator end_itr;
	for (boost::filesystem::directory_iterator itr(folder_path); itr != end_itr; ++itr) {
		if (is_directory(itr->status()))
			write_all_json_file(itr->path().string(), f);
	    else if (is_regular_file(itr->status())) {
	    	string filename = itr->path().string();
	    	if (filename[filename.length() - 1] != '~') {
	    		std::ifstream file(filename);
	    		while (!file.eof()) {
	    			string line;
	    			file >> line;
	    			*f << " " << line;
	    		}
	    		*f << endl;
	    	}
	    }
	    else
	    	sbrome_error_print("Can not process the folder or the file: " + itr->path().filename().string());
	  }
}
