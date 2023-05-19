#ifndef __CENTER_STATISTICS_H__
#define __CENTER_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class CenterStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	CenterStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] telescope telescope name
	/// @param[in] tag trigger tag, empty for origin trigger
	///
	CenterStatistics(
		unsigned int run,
		const std::string &telescope,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~CenterStatistics() = default;


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


	/// @brief return run number and telescope name in single string
	/// @returns run number and telescope name pair
	///
	virtual std::string Key() const override;


	/// @brief overloaded stream input function
	/// @param[in] is input stream
	/// @param[in] statistics statistics entry to read
	/// @returns input stream
	///
	friend std::istream& operator>>(
		std::istream &is,
		CenterStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const CenterStatistics &statistics
	);


	/// @brief get tag
	/// @returns tag
	///
	inline virtual std::string Tag() const {
		return tag_;
	}


	/// @brief get telescope name
	/// @returns telescope name
	///
	inline virtual std::string Telescope() const {
		return telescope_;
	}


	// x offset
	double x_offset[3];
	// y offset
	double y_offset[3];

private:
	// telescope name
	std::string telescope_;
	// tag
	std::string tag_;
};

}	// namespace ribll

#endif		// __CENTER_STATISTICS_H__