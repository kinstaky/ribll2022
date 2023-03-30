#ifndef __T0_HIT_STATISTICS_H__
#define __T0_HIT_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class T0HitStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	T0HitStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	/// @param[in] total total entries
	///
	T0HitStatistics(
		unsigned int run,
		const std::string &tag,
		long long total
	);


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
		T0HitStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const T0HitStatistics &statistics
	);


	long long d1_single_hit;
	long long d1_multi_hit;
	long long d2_single_hit;
	long long d2_multi_hit;
	long long d1d2_single_hit;
	long long d1d2_multi_hit;

private:
	std::string tag_;
	long long total_;
};

}		// namespace ribll

#endif		// __T0_HIT_STATISTICS_H__