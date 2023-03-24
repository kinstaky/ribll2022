#ifndef __BEAM_IDENTIFY_STATISTICS_H__
#define __BEAM_IDENTIFY_STATISTICS_H__

#include "include/statistics/statistics.h"

namespace ribll {

class BeamIdentifyStatistics : public Statistics {
public:

	/// @brief default constructor
	///
	BeamIdentifyStatistics() = default;


	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] total total number of entries
	///
	BeamIdentifyStatistics(unsigned int run, long long total);


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


	/// @brief return run number and beam type in single string
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
		BeamIdentifyStatistics &statistics
	);


	/// @brief overloaded stream output function
	/// @param[in] os output stream
	/// @param[in] statistics object to output
	/// @returns the output stream
	///
	friend std::ostream& operator<<(
		std::ostream &os,
		const BeamIdentifyStatistics &statistics
	);

	// constant value of gaus fitting of 14C
	double const14c;
	// mean value of gaus fitting of 14C
	double mean14c;
	// sigma value of gaus fitting of 14C
	double sigma14c;
	// total number of 14C
	long long c14;

private:
	long long total_;
};

}	// namespace ribll

#endif		// __BEAM_IDENTIFY_STATISTICS_H__