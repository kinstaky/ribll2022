#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <vector>
#include <map>

#include <TGraph.h>

#include "include/statistics/align_statistics.h"

namespace ribll {

class Alignment {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] group_num groups to divide
	/// @param[in] search_window width of search window, in nanoseconds
	/// @param[in] search_low_bound lower bound of search window, in nanoseconds
	/// @param[in] search_high_bound higher bound of search window, in nanoseconds
	///
	Alignment(
		unsigned int run,
		size_t group_num,
		double search_window,
		double search_low_bound,
		double search_high_bound
	);


	/// @brief destructor
	///
	virtual ~Alignment();


	/// @brief finding aligned time from VME trigger recored by XIA for vme events
	/// @returns 0 for successful, 1 otherwise
	///
	virtual int Align();


	/// @brief change verbose parameters for printing process
	/// @param[in] verbose value to change
	///
	inline void SetVerbose(bool verbose = true) {
		verbose_ = verbose;
	}


	/// @brief get statistics
	/// @returns align statistics
	inline AlignStatistics GetStatistics() const {
		return statistics_;
	}

private:

	/// @brief read VME trigger time recored by XIA
	/// @returns 0 for successful, -1 otherwise
	///
	int ReadXiaTime();


	/// @brief read VME time recored in v830 plugin (sdc in data)
	/// @returns 0 for successful, -1 otherwise
	///
	int ReadVmeTime();


	/// @brief align XIA and VME events in groups
	///
	void GroupAlignment();


	/// @brief align VME's second v1190
	/// @returns 0 for success, -1 otherwise
	///
	int AlignGdc();

	// run number
	unsigned int run_;

	// group number
	int group_num_;
	// search window size(step), low and high bound
	double search_window_;
	double search_low_bound_;
	double search_high_bound_;

	// xia data map
	std::vector<long long> xia_times_;
	// vme time list
	std::vector<long long> vme_times_;

	// align statistics
	AlignStatistics statistics_;

	// print processing information
	bool verbose_;
};

}


#endif			// __ALIGNMENT_H__