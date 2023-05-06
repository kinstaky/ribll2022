#ifndef __TRACE_MAPPER_H__
#define __TRACE_MAPPER_H__

namespace ribll {

class TraceMapper {
public:

	/// @brief constructor
	/// @param[in] run run number
/// @param[in] crate crate id, 0, 1 or 2
	///
	TraceMapper(unsigned int run, unsigned int crate);


	/// @brief default destructor
	///
	virtual ~TraceMapper() = default;


	/// @brief map the trace
	/// @returns 0 if sucess, -1 otherwise
	///
	virtual int Map();

protected:
	unsigned int run_;
	unsigned int crate_;
};

}	// namespace ribll

#endif // __TRA