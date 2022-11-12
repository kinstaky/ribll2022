#ifndef __DETECTOR_H__
#define __DETECTOR_H__

#include <string>
#include <map>
#include <iostream>

#include <TTree.h>

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
	int Correlate();



protected:
	// correlated event
	struct CorrelatedEvent {
		long long timestamp;
		unsigned short index;
		unsigned short x_hit;
		unsigned short y_hit;
		unsigned short x_strip[8];
		unsigned short y_strip[8];
		double x_time[8];
		double y_time[8];
		double x_energy[8];
		double y_energy[8];
	};


	// correlated and merged event
	struct MergedEvent {
		long long timestamp;
		unsigned short index;
		unsigned short hit;
		unsigned short x_strip[4];
		unsigned short y_strip[4];
		double time;
		double energy;
	};


	// residual correlated event
	CorrelatedEvent correlation_;
	// residual correlated events file
	TFile *correlated_file_;
	// residual correlated events tree
	TTree *correlated_tree_;


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


	/// @brief abstrct function setup output tree of correlated events
	///
	/// @param[in] tree output correlated tree to setup
	/// 
	virtual void SetupCorrelatedTree() = 0;

	
	/// @brief read detector events from mapped file and store in map
	///
	/// @returns 0 for success, -1 otherwise
	///
	int ReadEvents();


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
		bool used;
		unsigned short index;
		unsigned short side;
		unsigned short hit;
		unsigned short strip[8];
		double time[8];
		double energy[8];
	};

	// run number
	int run_;
	// detector name
	std::string name_;
	// single side correlation window
	unsigned int single_side_window_;
	// double sides correlation window
	unsigned int double_sides_window_;
	// events read from mapping file
	std::multimap<long long, SingleSideEvent> events_;
	// residual single event
	SingleSideEvent single_side_;

	// output files and output trees
	// residual single side events file
	TFile *single_side_file_;
	// residual single side events tree
	TTree *single_side_tree_;
	// correlated and merged events file
	TFile *merged_file_;
	// correlated and merged events tree
	TTree *merged_tree_;
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

private:

	/// @brief setup output tree of correlated events
	///
	/// @param[in] tree output correlated tree to setup
	/// 
	virtual void SetupCorrelatedTree() override;

};


}

#endif			// __DETECTOR_H__