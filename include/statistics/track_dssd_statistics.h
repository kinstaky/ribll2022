#ifndef __TRACK_DSSD_STATISTICS_H__
#define __TRACK_DSSD_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class TrackDssdStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	TrackDssdStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] detector detector name
	/// @param[in] tag trigger tag, empty for origin trigger
	///
	TrackDssdStatistics(
		unsigned int run,
		const std::string &telescope,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~TrackDssdStatistics() = default;


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
		TrackDssdStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const TrackDssdStatistics &statistics
	);


	// total number of valid events
	long long total;
	// events number with flag 1
	long long flag1;
	// events number with flag 3
	long long flag3;
	// events number with flag 5
	long long flag5;
	// events number with flag 6
	long long flag6;
	// events number with flag 7
	long long flag7;

private:
	// detector name
	std::string telescope_;
	// trigger tag
	std::string tag_;
};

}		// namespace ribll

#endif	// __TRACK_DSSD_STATISTICS_H__