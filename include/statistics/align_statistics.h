#ifndef __ALIGN_STATISTICS_H__
#define __ALIGN_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class AlignStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	AlignStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] xia_events total number of XIA events
	/// @param[in] vme_events total number of VME events
	/// @param[in] group group number
	///
	AlignStatistics(
		unsigned int run,
		long long xia_events,
		long long vme_events,
		unsigned int group
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


	/// @brief get total number of XIA events
	/// @returns XIA events count
	///
	inline long long XiaEvents() const {
		return xia_events_;
	}


	/// @brief get total number of VME events
	/// @returns VME events count
	inline long long VmeEvents() const {
		return vme_events_;
	}


	// number of events aligned in first window
	long long first_align_events;
	// number of events aligned in second window
	long long second_align_events;
	// number of events aligned in third window
	long long third_align_events;
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
	// group number
	unsigned int group_;
};

}		// namespace ribll

#endif		// __ALIGN_STATISTICS_H__