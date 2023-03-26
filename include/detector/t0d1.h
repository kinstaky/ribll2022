#ifndef __T0D1_H__
#define __T0D1_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d1 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	T0d1(unsigned int run);


	/// @brief default destructor
	///
	virtual ~T0d1() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	inline size_t FrontStrip() const override {
		return 64;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const override {
		return 64;
	}

// private:
// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for back strip
// 	///
// 	/// @param[in] correlation corrlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns true if pass, or false otherwise
// 	///
// 	virtual bool NormalizeFrontEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;


// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for front strip
// 	///
// 	/// @param[in] correlation correlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns ture if pass, or false otherwise
// 	///
// 	virtual bool NormalizeBackEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;
};



}		// namespace ribll

#endif		// __T0D1_H__