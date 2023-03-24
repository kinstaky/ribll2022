#include "include/statistics/csv_line_reader.h"

namespace ribll {

CsvLineReader::CsvLineReader(std::istream &is) {
	std::string buffer;
	std::getline(is, buffer);
	line_.str(buffer);
}


time_t CsvLineReader::ReadTime() {
	// temporary variable to store time
	std::string tmp;
	for (int i = 0; i < 6; ++i) std::getline(line_, tmp, ',');
	// statistics time in unix time
	time_t result;
	line_ >> result;
	return result;
}

}		// namespace ribll