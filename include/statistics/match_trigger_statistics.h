#ifndef __MATCH_STATISTICS_H__
#define __MATCH_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

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

}		// namespace ribll


#endif		// __MATCH_STATISTICS_H__