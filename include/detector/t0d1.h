#ifndef __T0D1_H__
#define __T0D1_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d1 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d1(unsigned int run, const std::string &tag);


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

	// //-------------------------------------------------------------------------
	// //								merge
	// //-------------------------------------------------------------------------

	// /// @brief merge adjacent event in the same side and merge events of two sides
	// /// @param[in] energy_diff tolerant energy relateive difference
	// /// @returns 0 if success, -1 otherwise
	// ///
	// virtual int Merge(double energy_diff) override;

protected:

	//-------------------------------------------------------------------------
	//								geometry
	//-------------------------------------------------------------------------

	/// @brief calculate the position from strip index
	/// @param[in] front_strip front strip
	/// @param[in] back_strip front strip
	/// @returns vector point to the position
	///
	virtual ROOT::Math::XYZVector CalculatePosition(
		double front_strip,
		double back_strip
	) const override;


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(TChain *chain, int iteration) override;


	/// @brief check whether energy is suitable for fitting
	/// @param[in] side side to normalize
	/// @param[in] event fundamental event
	/// @returns true if pass check, false not pass
	///
	virtual bool NormEnergyCheck(
		size_t side,
		const DssdNormalizeEvent &event
	) const override;


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