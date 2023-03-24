#ifndef __MAP_STATISTICS_H__
#define __MAP_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class MapStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	MapStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] crate crate number
	/// @param[in] threshold cut at threshold ? 
	///
	MapStatistics(unsigned int run, unsigned int crate, bool threshold);


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
	bool threshold_;
};

}		// namespace ribll

#endif 	//	__MAP_STATISTICS_H__