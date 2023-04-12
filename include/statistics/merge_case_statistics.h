#ifndef __MERGE_CASE_STATISTICS_H__
#define __MERGE_CASE_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class MergeCaseStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	MergeCaseStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] detector detector name
	/// @param[in] tag trigger tag, empty for origin trigger
	/// @param[in] case_name name of this case, major case
	/// @param[in] minor_case minor case in number
	/// @param[in] tolerance energy difference tolerance
	///
	MergeCaseStatistics(
		unsigned int run,
		const std::string &detector,
		const std::string &tag,
		const std::string &case_name,
		unsigned int minor_case,
		double tolerance
	);


	/// @brief default destructor
	///
	virtual ~MergeCaseStatistics() = default;


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
		MergeCaseStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const MergeCaseStatistics &statistics
	);


	// total number of valid events
	long long total;
	// merged events
	long long merged;

private:
	// detector name
	std::string detector_;
	// trigger tag
	std::string tag_;
	// merge case
	std::string case_name_;
	// minor case index
	unsigned int minor_case_;
	// energy difference tolerance
	double tolerance_;
};

}	// namespace ribll

#endif // __MERGE_CASE_STATISTICS_H__