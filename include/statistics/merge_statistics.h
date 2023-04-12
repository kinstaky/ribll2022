#ifndef __MERGE_STATISTICS_H__
#define __MERGE_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class MergeStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	MergeStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] detector detector name
	/// @param[in] tag trigger tag, empty for origin trigger
	///
	MergeStatistics(
		unsigned int run,
		const std::string &detector,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~MergeStatistics() = default;


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
		MergeStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const MergeStatistics &statistics
	);


	// total number of valid events
	long long total;
	// number of merged events
	long long merged;
	// one hit merged events
	long long one_hit;
	// two hit merged events
	long long two_hit;
	// three hit merged events
	long long three_hit;

private:
	// detector name
	std::string detector_;
	// trigger tag
	std::string tag_;
};


}	// namespace ribll

#endif	// __MERGE_STATISTICS_H__