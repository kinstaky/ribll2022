#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <ctime>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "include/defs.h"


namespace ribll {

class CsvLineReader {
public:

	/// @brief constructor
	/// @param[in] is istream to read
	/// 
	CsvLineReader(std::istream &is);

	/// @brief read store time and other time information
	/// @returns store time
	///
	time_t ReadTime();


	/// @brief read one value from istream
	/// @tparam T type of variable to read
	/// @param[in] t varaiable to read
	/// @returns reference to itself
	/// 
	template<typename T>
	CsvLineReader& operator>>(T &t) {
		std::string value;
		std::getline(line_, value, ',');
		std::stringstream ss;
		ss.str(value);
		ss >> t;
		return *this;
	}
private:
	std::stringstream line_;
};

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


	/// @brief return run number and crate id in single string
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
	friend std::ostream& operator<<(
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




class MatchTriggerStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	MatchTriggerStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] detector detector name
	/// @param[in] tag trigger tag, empty for origin trigger
	/// @param[in] extract_tag extract tag, empty for just matching
	/// @param[in] reference_events total number of reference trigger events
	/// @param[in] mapped_events total number of input mapped events
	///
	MatchTriggerStatistics(
		unsigned int run,
		const std::string &detector,
		const std::string &tag,
		const std::string &extract_tag,
		long long reference_events,
		long long mapeed_events
	);


	/// @brief default destructor
	///
	virtual ~MatchTriggerStatistics() = default;


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


	/// @brief return run number and detector name in single string
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
		MatchTriggerStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const MatchTriggerStatistics &statistics
	);


	// number of output fundamental events
	long long match_events;
	// number of used input mapped events
	long long used_events;
	// number of invalid input mapped events, because of noise or ...
	long long oversize_events;
	// number of conflict events, which can match multiple triggers
	long long conflict_events;

private:
	// detector name
	std::string detector_;
	// trigger tag
	std::string tag_;
	// extract
	std::string extract_tag_;
	// total number of main trigger
	long long reference_events_;
	// total number of input mapped events
	long long mapped_events_;
};


class XiaTriggerPeriodStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	XiaTriggerPeriodStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] period period time, in nanosecond
	///
	XiaTriggerPeriodStatistics(unsigned int run, double period);


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
		XiaTriggerPeriodStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const XiaTriggerPeriodStatistics &statistics
	);

private:
	double period_;
};



class BeamIdentifyStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	BeamIdentifyStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] total total number of entries
	///
	BeamIdentifyStatistics(unsigned int run, long long total);


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


	/// @brief return run number and beam type in single string
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
		BeamIdentifyStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const BeamIdentifyStatistics &statistics
	);

	// constant value of gaus fitting of 14C
	double const14c;
	// mean value of gaus fitting of 14C
	double mean14c;
	// sigma value of gaus fitting of 14C
	double sigma14c;
	// total number of 14C
	long long c14;

private:
	long long total_;
};



//-----------------------------------------------------------------------------
//					implementation of template functions
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