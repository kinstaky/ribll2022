#ifndef __XT_PERIOD_STATISTICS_H__
#define __XT_PERIOD_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

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

}		// namespace ribll

#endif		// __XT_PERIOD_STATISTICS_H__