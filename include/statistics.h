#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <ctime>

#include <string>
#include <fstream>
#include <iostream>
#include <map>

#include "include/defs.h"


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

protected:
	// run number
	unsigned int run_;
	// unix time
	time_t store_time_;
};



class MapStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	MapStatistics() = default;

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] crate crate number
	///
	MapStatistics(unsigned int run, unsigned int crate);


	/// @brief default destructor
	///
	virtual ~MapStatistics() = default;


	/// @brief write this statistics entry to file
	///
	virtual void Write() override;


	/// @brief print the statistics to stdout
	///
	virtual void Print() const override;


	/// @brief get the title of this type statistics entry
	/// @returns title in string format, separated by ','
	///
	virtual std::string Title() const override;


	/// @brief return run number and detector name pair as key in map
	/// @returns run number and detector name pair
	///
	virtual std::string Key() const override;


	/// @brief overloaded stream input function
	/// @param[in] is input stream
	/// @param[in] statistics statistics entry to read
	/// @returns input stream
	///
	friend std::istream& operator>>(
		std::istream &is,
		MapStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const MapStatistics &statistics
	);

private:
	unsigned int crate_;
};



class AlignStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	AlignStatistics() = default;

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] xia_events total number of XIA events
	/// @param[in] vme_events total number of VME events
	/// @param[in] calibration_parameters array of calibration parameters
	///
	AlignStatistics(
		unsigned int run,
		long long xia_events,
		long long vme_events,
		double *calibration_parameters
	);


	/// @brief default destructor
	///
	virtual ~AlignStatistics() = default;


	/// @brief write this statistics entry to file
	///
	virtual void Write() override;


	/// @brief print the statistics to stdout
	///
	virtual void Print() const override;


	/// @brief get the title of this type statistics entry
	/// @returns title in string format, separated by ','
	///
	virtual std::string Title() const override;


	/// @brief overloaded stream input function
	/// @param[in] is input stream
	/// @param[in] statistics statistics entry to read
	/// @returns input stream
	///
	friend std::istream& operator>>(
		std::istream &is,
		AlignStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream &operator<<(
		std::ostream &os,
		const AlignStatistics &statistics
	);


	// number of VME events can be aligned
	long long align_events;
	// number of VME events match more than one VME trigger in XIA
	long long oversize_events;

private:
	// detector name
	std::string detector_;
	// total number of XIA events
	long long xia_events_;
	// total number of VME events
	long long vme_events_;
	// pointer to calibration parameters
	double calibration_param_[2];
};


//-----------------------------------------------------------------------------
//							implementation of template functions
//-----------------------------------------------------------------------------

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


}

#endif 		// __STATISTICS_H__