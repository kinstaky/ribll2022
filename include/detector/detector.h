#ifndef __DETECTOR_H__
#define __DETECTOR_H__

#include <string>
#include <map>
#include <iostream>
#include <memory>

#include <TTree.h>
#include <TGraph.h>

namespace ribll {

class Detector {
public:

	Detector(
		int run,
		const std::string &name,
		unsigned int single_side_window,
		unsigned int double_sides_window
	);


	/// @brief destructor
	///
	virtual ~Detector();


	/// @brief correlate events in detector
	///
	/// @returns 0 if success, -1 otherwise
	/// 
	virtual int Correlate();


	/// @brief calculate normalize parameters
	///
	/// @param[in] length number of run to chain
	/// @param[in] ref_front reference front strip
	/// @param[in] ref_back reference back strip
	/// @param[in] iteration run in iteration mode
	/// @returns 0 for succes, -1 otherwise
	///
	virtual int Normalize(
		int length,
		unsigned short ref_front,
		unsigned short ref_back,
		bool iteration
	);


	/// @brief build normalize result
	virtual int NormalCorrelate();


	/// @brief merge events in adjacent strip
	///
	/// @param[in] energy_cut front back energy cut
	/// @returns 0 if success, -1 otherwise
	/// 
	virtual int Merge(double energy_cut = 0.02);


	/// @brief get front strip strip
	///
	/// @returns front strip number
	///
	virtual inline size_t FrontStrip() {
		return 0;
	}


	/// @brief get back strip number
	///
	/// @returns back strip number
	/// 
	virtual inline size_t BackStrip() {
		return 0;
	}


protected:
	// correlated event
	struct CorrelatedEvent {
		unsigned short index;
		unsigned short front_hit;
		unsigned short back_hit;
		long long timestamp;
		unsigned short front_strip[8];
		unsigned short back_strip[8];
		double front_time[8];
		double back_time[8];
		double front_energy[8];
		double back_energy[8];
	};


	// correlated and merged event
	struct MergedEvent {
		unsigned short hit;
		long long timestamp;
		double front_strip[4];
		double back_strip[4];
		double time[4];
		double energy[4];
	};


	// run number
	int run_;
	// detector name
	std::string name_;
	// residual correlated event
	CorrelatedEvent correlation_;
	// merged event
	MergedEvent merged_;
	// events read from corrleated file
	std::multimap<long long, CorrelatedEvent> correlated_events_;
	// residual correlated events file
	TFile *correlated_file_;
	// residual correlated events tree
	TTree *correlated_tree_;
	// correlated and merged events file
	TFile *merged_file_;
	// correlated and merged events tree
	TTree *merged_tree_;
	// residual correlated file
	TFile *residual_correlated_file_;
	// residual correlated tree
	TTree *residual_correlated_tree_;
	// normalize back parameters
	double normalize_back_param_[64][2];
	// normalize front parameters
	double normalize_front_param_[64][2];


	//-------------------------------------------------------------------------
	//						
	//-------------------------------------------------------------------------


	/// @brief setup residual single side output file and tree
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int SetupSingleSideOutput();


	/// @brief setup residual correlated output file and tree
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int SetupCorrelatedOutput();


	/// @brief setup residual correlated output file
	///
	/// @returns 0 for success, -1 otherwise
	/// 
	virtual int SetupResidualCorrelated();


	/// @brief abstrct function setup output tree of correlated events
	///
	/// @param[in] tree output correlated tree to setup
	/// 
	virtual void SetupCorrelatedTree(TTree *tree) = 0;


	/// @brief read correlated events from root file
	///
	/// @returns 0 for success, -1 otherwise
	/// 
	virtual void SetupCorrelatedInputTree(TTree *tree) = 0;


	/// @brief setup merged output file and tree
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int SetupMergedOutput();

	
	/// @brief calculate normalize parameters of first side
	///
	/// @param[in] chain input tree 
	/// @param[in] ref_front reference front strip
	/// @param[in] ref_back reference back strip
	/// @param[in] g_front_back_energy list of graph of front-back energy 
	/// @param[in] g_back_front_energy list of graph of back-front energy
	/// @param[in] iteration iteration mode
	/// 
	virtual void NormalizeFirstSide(
		TTree *chain,
		unsigned short ref_front,
		unsigned short ref_back,
		TGraph **g_front_back_energy,
		TGraph **g_back_front_energy,
		bool iteration
	);


	/// @brief calculate normalize parameters of first side
	///
	/// @param[in] chain input tree 
	/// @param[in] ref_front reference front strip
	/// @param[in] ref_back reference back strip
	/// @param[in] g_front_back_energy list of graph of front-back energy 
	/// @param[in] g_back_front_energy list of graph of back-front energy
	/// @param[in] iteration iteration mode
	/// 
	virtual void NormalizeSecondSide(
		TTree *chain,
		unsigned short ref_front,
		unsigned short ref_back,
		TGraph **g_front_back_energy,
		TGraph **g_back_front_energy,
		bool iteration
	);


	/// @brief check whether the energy of correlated event can be use in normalization for back strip
	///
	/// @param[in] correlation corrlated event 
	/// @param[in] iteration iteration mode 
	/// @returns true if pass, or false otherwise
	///
	virtual bool NormalizeFrontEnergyCheck(const CorrelatedEvent &correlation, bool iteration) = 0;


	/// @brief check whether the energy of correlated event can be use in normalization for front strip
	///
	/// @param[in] correlation correlated event 
	/// @param[in] iteration iteration mode 
	/// @returns ture if pass, or false otherwise
	/// 
	virtual bool NormalizeBackEnergyCheck(const CorrelatedEvent &correlation, bool iteration) = 0;


	/// @brief calculate normalized energy
	///
	/// @param[in] side 0 for front, 1 for back
	/// @param[in] strip strip number
	/// @param[in] energy energy before normalization
	/// @returns energy after normalization
	///
	inline double NormalizeEnergy(size_t side, unsigned short strip, double energy) {
		if (side == 0) {
			return normalize_front_param_[strip][0] + normalize_front_param_[strip][1] * energy;
		} else {
			return normalize_back_param_[strip][0] + normalize_back_param_[strip][1] * energy;
		}
	}


private:
	// single event before correlation
	struct Event {
		bool used;
		unsigned short index;
		unsigned short side;
		unsigned short strip;
		double time;
		double energy;
	};

	// single side event after single-side-correlation before double-sides-correlation
	struct SingleSideEvent {
		bool correlated;
		unsigned short index;
		unsigned short side;
		unsigned short hit;
		long long timestamp;
		unsigned short strip[8];
		double time[8];
		double energy[8];
	};

	
	// single side correlation window
	unsigned int single_side_window_;
	// double sides correlation window
	unsigned int double_sides_window_;
	// events read from mapping file
	std::multimap<long long, SingleSideEvent> single_side_events_;
	// residual single event
	SingleSideEvent single_side_;

	// output files and output trees
	// residual single side events file
	TFile *single_side_file_;
	// residual single side events tree
	TTree *single_side_tree_;

	/// @brief read detector events from mapped file and store in map
	///
	/// @returns 0 for success, -1 otherwise
	///
	int ReadEvents();


	/// @brief sort events in single-side-event by strip
	///
	/// @param[in] event pointer to single-side-event
	/// 
	void BubbleSort(SingleSideEvent *event);


	/// @brief read normalize parameters from file
	///
	/// @returns 0 for success, -1 otherwise
	///
	int ReadNormalizeParameters();


	/// @brief write normalize parameters to file
	///
	/// @returns 0 for success, -1 otherwise
	/// 
	int WriteNormalizeParameters();


	/// @brief store residual single-side events
	///
	void StoreResidualSingleSideEvents();
};


class DSSD : public Detector {
public:
	
	DSSD(
		int run,
		const std::string &name,
		unsigned int single_side_window = 200,
		unsigned int double_sides_window = 200
	);


	/// @brief default destructor
	///
	virtual ~DSSD() = default;


	/// @brief get front strip number 
	///
	/// @returns front strip number
	/// 
	inline size_t FrontStrip() override {
		return 32;
	}


	/// @brief get back strip number
	///
	/// @returns back strip number
	/// 
	virtual inline size_t BackStrip() {
		return 32;
	}

private:

	/// @brief setup output tree of correlated events
	///
	/// @param[in] tree output correlated tree to setup
	/// 
	virtual void SetupCorrelatedTree(TTree *tree) override;


	/// @brief read correlated events from root file
	///
	/// @returns 0 for success, -1 otherwise
	/// 
	virtual void SetupCorrelatedInputTree(TTree *tree) override;

};


class T0D1 : public DSSD {
public:
	
	T0D1(
		int run,
		const std::string &name,
		unsigned int single_side_window = 200,
		unsigned int double_sides_window = 200
	);


	/// @brief default destructor
	///
	virtual ~T0D1() = default;


	/// @brief get front strip number 
	///
	/// @returns front strip number
	/// 
	inline size_t FrontStrip() override {
		return 64;
	}


	/// @brief get back strip number
	///
	/// @returns back strip number
	/// 
	virtual inline size_t BackStrip() {
		return 64;
	}

private:
	/// @brief check whether the energy of correlated event can be use in normalization for back strip
	///
	/// @param[in] correlation corrlated event 
	/// @param[in] iteration iteration mode 
	/// @returns true if pass, or false otherwise
	///
	virtual bool NormalizeFrontEnergyCheck(const CorrelatedEvent &correlation, bool iteration) override;


	/// @brief check whether the energy of correlated event can be use in normalization for front strip
	///
	/// @param[in] correlation correlated event 
	/// @param[in] iteration iteration mode 
	/// @returns ture if pass, or false otherwise
	/// 
	virtual bool NormalizeBackEnergyCheck(const CorrelatedEvent &correlation, bool iteration) override;


	/// @brief calculate normalize parameters of first side
	///
	/// @param[in] chain input tree
	/// @param[in] ref_front reference front strip 
	/// @param[in] ref_back reference back strip
	/// @param[in] g_front_back_energy list of graph of front-back energy 
	/// @param[in] g_back_front_energy list of graph of back-front energy
	/// @param[in] iteration iteration mode
	/// 
	virtual void NormalizeFirstSide(
		TTree *chain,
		unsigned short ref_front,
		unsigned short ref_back,
		TGraph **g_front_back_energy,
		TGraph **g_back_front_energy,
		bool iteration
	) override;


	/// @brief calculate normalize parameters of first side
	///
	/// @param[in] chain input tree 
	/// @param[in] ref_front reference front strip
	/// @param[in] ref_back reference back strip
	/// @param[in] g_front_back_energy list of graph of front-back energy 
	/// @param[in] g_back_front_energy list of graph of back-front energy
	/// @param[in] iteration iteration mode
	/// 
	virtual void NormalizeSecondSide(
		TTree *chain,
		unsigned short ref_front,
		unsigned short ref_back,
		TGraph **g_front_back_energy,
		TGraph **g_back_front_energy,
		bool iteration
	) override;

};


class T0D3 : public DSSD {
public:
	
	T0D3(
		int run,
		const std::string &name,
		unsigned int single_side_window = 200,
		unsigned int double_sides_window = 200
	);


	/// @brief default destructor
	///
	virtual ~T0D3() = default;


	/// @brief get front strip number 
	///
	/// @returns front strip number
	/// 
	inline size_t FrontStrip() override {
		return 32;
	}


	/// @brief get back strip number
	///
	/// @returns back strip number
	/// 
	virtual inline size_t BackStrip() {
		return 32;
	}
};


}

#endif			// __DETECTOR_H__