#ifndef __CSV_LINE_READER_H__
#define __CSV_LINE_READER_H__

#include <ctime>

#include <sstream>
#include <string>

namespace ribll {

class CsvLineReader {
public:

	/// @brief constructor
	/// @param[in] is istream to read
	/// 
	CsvLineReader(std::istream &is);

	/// @brief read store time and other time information
	/// @returns store time
	///
	time_t ReadTime();


	/// @brief read one value from istream
	/// @tparam T type of variable to read
	/// @param[in] t varaiable to read
	/// @returns reference to itself
	/// 
	template<typename T>
	CsvLineReader& operator>>(T &t) {
		std::string value;
		std::getline(line_, value, ',');
		std::stringstream ss;
		ss.str(value);
		ss >> t;
		return *this;
	}
private:
	std::stringstream line_;
};

} 		// namespace ribll

#endif 		// __CSV_LINE_READER_H__