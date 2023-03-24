#include "include/check/chain_check.h"

#include <map>
#include <iostream>
#include <vector>
#include <fstream>
#include <exception>

#include "include/statistics/align_statistics.h"
#include "include/statistics/map_statistics.h"
#include "include/statistics/match_trigger_statistics.h"

namespace ribll {

const std::map<std::string, std::vector<std::string>> node_dependencies = {
	{"map-0", {}},
	{"map-1", {"map-0"}},
	{"map-2", {}},
	{"map-3", {"align"}},
	{"align", {"map-0"}},
	{"match-vt", {"map-0", "map-3"}},
	{"match-tof", {"map-0"}},
	{"match-vtof", {"map-3", "match-vt"}},
	{"match-xppac", {"map-0"}},
	{"match-vppac", {"map-3", "match-vt"}},
	{"match-t0d1", {"map-0", "map-2"}},
	{"match-t0d2", {"map-0", "map-1"}},
	{"match-t0d3", {"map-0", "map-1"}},
	{"match-t1d1", {"map-0", "map-1"}},
	{"match-t0s1", {"map-0"}},
	{"match-t0s2", {"map-0"}},
	{"match-t0s3", {"map-0"}},
	{"match-t1s1", {"map-0"}},
	{"match-taf0", {"map-3", "match-vt"}},
	{"match-taf1", {"map-3", "match-vt"}},
	{"match-taf2", {"map-0"}},
	{"match-taf3", {"map-0"}},
	{"match-taf4", {"map-0"}},
	{"match-taf5", {"map-0"}},
	{"match-tab0", {"map-3", "match-vt"}},
	{"match-tab1", {"map-3", "match-vt"}},
	{"match-tab2", {"map-3", "match-vt"}},
	{"match-tab3", {"map-3", "match-vt"}},
	{"match-tab4", {"map-3", "match-vt"}},
	{"match-tab5", {"map-3", "match-vt"}},
	{"match-t0csi", {"map-0"}},
	{"match-t1csi", {"map-0"}},
	{"match-tafcsi", {"map-0"}},
	{"match-tabcsi", {"map-0"}},
};


bool NodeCheck(const std::string &key) {
	auto search = node_dependencies.find(key);
	return search != node_dependencies.end();
}


struct NodeInfo {
	// key of node in dependency map
	std::string key;
	// unix timestap of this node
	time_t time;
	// wheter this node is up to date
	bool up_to_date;
};


template<typename Statistics>
time_t ReadEditTime(const std::string &file_name, const std::string &key) {
	// statistics to read
	Statistics statistics;
	// statistics file stream
	std::ifstream fin(file_name);
	if (!fin.good()) {
		std::cout << "Error: Open statistics file "
			<< file_name << " failed.\n";
		return -1;
	}
	// buffer to read first title line
	std::string buffer;
	// read first title line
	std::getline(fin, buffer);
	// read statistics file line by line
	while (fin.good()) {
		fin >> statistics;
		if (statistics.Key() >= key) break;
	}
	// close file
	fin.close();

	// statistics found
	if (statistics.Key() == key) return statistics.Time();

	// statistics not found
	return -1;
}



/// @brief recursion function to check time
/// @param[in] run run number
/// @param[in] key key of node to check
/// @param[in] checked_nodes checked nodes
/// @returns true if node is up to date, false if outdated
///
bool CheckTime(
	unsigned int run,
	const std::string &key,
	std::vector<NodeInfo> &checked_nodes
) {
	// check whether the key is dependency map
	auto search = node_dependencies.find(key);
	// node not found, fail
	if (search == node_dependencies.end()) {
		throw std::runtime_error("Invalid node key " + key);
	}

	// check whether the node is in checked nodes
	for (const auto &info : checked_nodes) {
// std::cout << key << " node checked.\n";
		if (info.key == key) return info.up_to_date;
	}

	// last stored time of this node read from statistics file
	time_t store_time = -1;
	// read edit time and check whether the process has run successfully
	// node type, map, align, match...
	std::string type = key.substr(0, key.find_first_of('-'));
	// name of statistics file that stores information of this node
	std::string file_name =
		std::string(kGenerateDataPath)+ "statistics/" + type + ".csv";
	if (type == "map") {
		unsigned int crate = key[key.find_first_of('-') + 1] -'0';
		MapStatistics statistics(run, crate, true);
		store_time =
			ReadEditTime<MapStatistics>(file_name, statistics.Key());
	} else if (type == "align") {
		double tmp[2];
		AlignStatistics statistics(run, 0, 0, tmp);
		store_time =
			ReadEditTime<AlignStatistics>(file_name, statistics.Key());
	} else if (type == "match") {
		std::string detector = key.substr(key.find_first_of('-') + 1);
		MatchTriggerStatistics statistics(run, detector, "", "", 0, 0);
		store_time =
			ReadEditTime<MatchTriggerStatistics>(file_name, statistics.Key());
	}
	if (store_time < 0) {
		// this process hasn't run yet
		checked_nodes.push_back({key, -1, false});
// std::cout << key << " store time not found.\n";
		return false;
	}

	// The node must be up to date if it doesn't have dependencies.
	if (search->second.empty()) {
		checked_nodes.push_back({key, store_time, true});
// std::cout << key << " standalone.\n";
		return true;
	}

	// check dependencies
	for (const auto &node : search->second) {
		if (!CheckTime(run, node, checked_nodes)) {
			checked_nodes.push_back({key, store_time, false});
// std::cout << key << " dependencies out of date.\n";
			return false;
		}
	}

	// dependencies is up to date, compare this node and its dependencies
	for (const auto &node : search->second) {
		for (const auto &info : checked_nodes) {
			if (info.key != node) continue;
			// found the dependency in checked_nodes
			if (info.time > store_time) {
				checked_nodes.push_back({key, store_time, false});
// std::cout << key << " later than dependence " << info.key << "\n";
				return false;
			}
		}
	}

	checked_nodes.push_back({key, store_time, true});
	return true;
}




bool ChainCheck(unsigned int run, const std::string &key) {
	// nodes to check
	std::vector<NodeInfo> checked_nodes;

	// check time
	if (CheckTime(run, key, checked_nodes)) {
		return true;
	}

	// Not all node up to date and show outdated nodes.
	for (const auto &info : checked_nodes) {
		if (info.up_to_date) continue;
		std::cout << "\033[1;31mOutdated\033[0m " << info.key << "\n";
	}

	return false;
}

}	// namespace ribll