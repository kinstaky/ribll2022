#ifndef __NORMALIZE_STATISTICS_H__
#define __NORMALIZE_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class NormalizeStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	NormalizeStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] detector detector name
	/// @param[in] tag trigger tag, empty for origin trigger
	/// @param[in] end_run end run number
	/// @param[in] iteration iteration mode
	///
	NormalizeStatistics(
		unsigned int run,
		const std::string &detector,
		const std::string &tag,
		unsigned int end_run,
		int iteration
	);


	/// @brief default destructor
	///
	virtual ~NormalizeStatistics() = default;


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
		NormalizeStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const NormalizeStatistics &statistics
	);


private:
	// detector name
	std::string detector_;
	// trigger tag
	std::string tag_;
	// end run chained
	unsigned int end_run_;
	// iteration
	int iteration_;
};


}	// namespace ribll

#endif	// __NORMALIZE_STATISTICS_H__