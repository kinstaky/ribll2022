#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <ctime>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "include/defs.h"
#include "include/statistics/csv_line_reader.h"


namespace ribll {

class Statistics {
public:

	/// @brief default constructor
	///
	Statistics() = default;


	/// @brief constructor
	/// @param[in] run  run number
	///
	Statistics(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Statistics() = default;


	/// @brief write this statistics entry to file
	/// @tparam Entry type of statistics to write
	/// @param type type of statistics
	///
	template<typename Entry>
	void Write(const std::string &type);


	/// @brief write this statistics entry to file
	///
	virtual void Write() = 0;


	/// @brief print the statistics to stdout
	///
	virtual void Print() const = 0;


	/// @brief get the title of this type statistics entry
	/// @returns title in string format, separated by ','
	///
	virtual std::string Title() const = 0;


	/// @brief return run number and detector name pair as key in map
	/// @returns run number and detector name pair
	///
	virtual std::string Key() const;


	/// @brief get store time
	/// @returns store time
	///
	inline virtual time_t Time() const {
		return store_time_;
	}

protected:
	// run number
	unsigned int run_;
	// unix time
	time_t store_time_;
};


template<typename Entry>
void Statistics::Write(const std::string &type) {
	// update store time
	store_time_ = time(NULL);

	// statistics file name
	std::string file_name =
		std::string(kGenerateDataPath) + "statistics/" + type + ".csv";
	// statistics file stream
	std::ifstream fin(file_name);

	// map of entries
	std::map<std::string, Entry> entries;
	// single entry to read lines
	Entry entry;

	if (!fin.good()) {
		std::cout << "Open statistic file "
			<< file_name << " for reading failed.\n";
	} else {
		// read statistics data from file
		// buffer to read line
		char buffer[1024];
		// read the first title line
		fin.getline(buffer, sizeof(buffer));

		// loop to read entries
		while (fin.good()) {
			fin >> entry;
			entries.insert(std::make_pair(entry.Key(), entry));
		}
		// close file
		fin.close();
	}

	// use entry to store statistics
	entry = *(static_cast<Entry*>(this));
	// try to find the entry with the key(run number and detector name)
	auto search = entries.find(entry.Key());
	if (search != entries.end()) {
		search->second = entry;
	} else {
		entries.insert(std::make_pair(entry.Key(), entry));
	}

	// write statistics
	std::ofstream fout(file_name);
	if (!fout.good()) {
		std::cerr << "Error: Open statistic file "
			<< file_name << " for writing failed.\n";
		return;
	}
	fout << this->Title() << "\n";
	for (const auto &[key, value] : entries) {
		fout << value << "\n";
	}
	fout.close();
	return;
}


const std::string title_time = ",year,month,day,hour,minute,second,unix_time";


/// @brief write statistics store time at the end of line
/// @param[in] os output stream
/// @param[in] statistics_time store time
///
void WriteStatisticsTime(std::ostream &os, time_t statistics_time);



}

#endif 		// __STATISTICS_H__