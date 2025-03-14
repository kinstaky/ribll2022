#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <fstream>
#include <chrono>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TCutG.h>
#include <TH2F.h>
#include <TRandom3.h>

#include "include/event/dssd_event.h"
#include "include/event/ssd_event.h"
#include "include/event/t0_event.h"

using namespace ribll;

constexpr bool print_debug = false;
constexpr double d1_range = 500.0;
constexpr double d2_range = 1000.0;
constexpr double d3_range = 200.0;

// straight parameters
constexpr double sa12 = 0.47;
constexpr double sb12 = -0.042;
constexpr double sa23 = 0.59;
constexpr double sb23 = -0.042;


/// @brief get fixed energy after straight
/// @param[in] de energy(channel) in the first layer
/// @param[in] e energy(channel) in the second layer 
/// @param[in] a straight parameter a
/// @param[in] b straight parameter b
/// @returns fixed energy after straight
/// 
inline double Straight(
	const double de,
	const double e,
	const double a, 
	const double b
) {
	return sqrt(de*e + a*de*de) + b*e;
}

// const double cy_param[6][2] = {
// 	{-0.451425,	0.00679791},
// 	{0.217737, 0.00697266},
// 	{0.121048, 0.00565964},
// 	{-0.578595, 0.00231362},
// 	{1.03268, 0.00235314},
// 	{-2.97632, 0.00217422}
// };

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run               Set run number.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag
) {
	// initialize
	help = false;
	trigger_tag.clear();
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else if (argv[result][1] == 't') {
			// option of trigger tag
			// get tag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			trigger_tag = argv[result];
		} else {
			return -result;
		}
	}
	return result;
}



struct MergeEvent {
	// used fundamental event flag
	unsigned int flag;
	// tag, 0-f1b1, 1-f1b2, 2-f2b1, 3-f2b2
	unsigned int tag;
	// is it hole?
	bool hole;
	// position
	double x, y;
	// energy
	double energy;
};


struct MergeEventNode {
	// index in merge event vector
	int index;
	// parent's index in tree
	int prev;
	// total number of merge events, also means layer
	int num;
	// used event flag
	unsigned int flag;
	// is it leaf?
	bool leaf;
	// is it hole?
	bool hole;
};


class HoleInfo {
public:
	bool simulate;
	bool hole_flag[32][32];
	double hole_possibility[32][32];
	TRandom3 *generator;

	bool IsHole(int fs, int bs) const {
		if (simulate) {
			double random = generator->Rndm();
			return random < hole_possibility[fs][bs];
		} else {
			return hole_flag[fs][bs];
		}
		return false;
	}
};


/// @brief get all possible merged events from fundamental events
/// @param[in] event DSSD fundamental event
/// @param[in] range front-back correlation, energy difference in range
/// @param[in] positionFromStrip function to calculate postion from strip
/// @param[in] hole T0D2 hole flag
/// @param[out] result DSSD merge events
/// @param[out] hole wheter this merged event includes hole pixel
///
void GetMergeEvents(
	const DssdFundamentalEvent &event,
	double range,
	void (*positionFromStrip)(double, double, double&, double&),
	const HoleInfo *hole_info,
	std::vector<DssdMergeEvent> &result,
	std::vector<std::vector<bool>> &hole
) {
	// initialize
	result.clear();

	// for convenient
	const int fhit = event.front_hit;
	const int bhit = event.back_hit;
	const unsigned short *fs = event.front_strip;
	const unsigned short *bs = event.back_strip;
	const double *fe = event.front_energy;
	const double *be = event.back_energy;

	// record front and back event index
	// if it's single strip, set the second paramters as -1
	std::vector<std::pair<int, int>> front_events;
	std::vector<std::pair<int, int>> back_events;
	// loop front events
	// single strip events
	for (int i = 0; i < fhit; ++i) {
		if (fe[i] < 200.0) break;
		front_events.push_back(std::make_pair(i, -1));
	}
	// adjacent strip events
	for (int i = 0; i < fhit; ++i) {
		if (fe[i] < 200.0) break;
		for (int j = i+1; j < fhit; ++j) {
			if (fe[j] < 200.0) break;
			// possible adjacent strip
			if (abs(fs[i]-fs[j]) == 1) {
				front_events.push_back(std::make_pair(i, j));
			}
		}
	}
	// loop back events
	// single strip events
	for (int i = 0; i < bhit; ++i) {
		if (be[i] < 200.0) break;
		back_events.push_back(std::make_pair(i, -1));
	}
	// adjacent strip events
	for (int i = 0; i < bhit; ++i) {
		if (be[i] < 200.0) break;
		for (int j = i+1; j < bhit; ++j) {
			if (be[j] < 200.0) break;
			// possible adjacent strip
			if (abs(bs[i]-bs[j]) == 1) {
				back_events.push_back(std::make_pair(i, j));
			}
		}
	}

	// all possible merged events
	std::vector<MergeEvent> merge_events;
	// front-back correlation
	for (const auto &front : front_events) {
		// front energy sum
		double fe_sum = fe[front.first];
		// front flag
		unsigned int fflag = 0x1 << front.first;
		// front tag
		unsigned int ftag = 0;
		// front strip
		double fstrip = fs[front.first];
		// adjacent strip
		if (front.second != -1) {
			fe_sum += fe[front.second];
			fflag |= 0x1 << front.second;
			ftag = 2;
			fstrip += fs[front.first] > fs[front.second] ? -0.5 : 0.5;
		}

		for (const auto &back : back_events) {
			// back energy sum
			double be_sum = be[back.first];
			// back flag
			unsigned int bflag = 0x100 << back.first;
			// back tag
			unsigned int btag = 0;
			// back strip
			double bstrip = bs[back.first];
			// adjacent strip
			if (back.second != -1) {
				be_sum += be[back.second];
				bflag |= 0x100 << back.second;
				btag = 1;
				bstrip += bs[back.first] > bs[back.second] ? -0.5 : 0.5;
			}
			if (fabs(fe_sum - be_sum) < range) {
				MergeEvent merge;
				merge.flag = fflag | bflag;
				merge.tag = ftag + btag;
				merge.energy = fe_sum;
				positionFromStrip(fstrip, bstrip, merge.x, merge.y);
				// check hole
				merge.hole = false;
				if (hole_info) {
					bool is_hole = false;
					int pixels[4][2] = {
						{fs[front.first], bs[back.first]},
						{fs[front.first], bs[back.second]},
						{fs[front.second], bs[back.first]},
						{fs[front.second], bs[back.second]}
					};
					// check pixels
					for (int i = 0; i < 1; ++i) {
						int fstrip = pixels[i][0];
						int bstrip = pixels[i][1];
						if (fstrip < 0 || bstrip < 0) continue;
						if (hole_info->IsHole(fstrip, bstrip)) is_hole = true;
					}
					// assign hole flag
					merge.hole = is_hole;
				}
				merge_events.push_back(merge);
			}
		}
	}

	// merge event nodes, store in array but it's in tree structure
	std::vector<MergeEventNode> nodes;
	// generate forest (serveral trees) from merge events
	for (size_t i = 0; i < merge_events.size(); ++i) {
		// make each event as root of tree
		MergeEventNode root;
		root.index = i;
		root.prev = -1;
		root.num = 1;
		root.flag = merge_events[i].flag;
		root.leaf = true;
		nodes.push_back(root);

		// new nodes will be add to tree later
		std::vector<MergeEventNode> new_nodes;
		// check all nodes and append this event as new child if possible
		for (size_t j = 0; j < nodes.size(); ++j) {
			// it's possible if all fundamental event is not used
			if ((merge_events[i].flag & nodes[j].flag) == 0) {
				MergeEventNode node;
				node.index = i;
				node.prev = j;
				node.num = nodes[j].num + 1;
				node.flag = nodes[j].flag | merge_events[i].flag;
				node.leaf = true;
				new_nodes.push_back(node);
				// set it's parent as NOT leaf
				nodes[j].leaf = false;
			}
		}

		// move new nodes into existed forest
		for (const auto &node : new_nodes) {
			nodes.push_back(node);
		}
	}

	std::vector<std::vector<MergeEvent>> branches;
	for (size_t i = 0; i < nodes.size(); ++i) {
		// jump if it's not leaf
		if (!nodes[i].leaf) continue;

		// branch with serveral nodes, the same as DSSDMergeEvent
		std::vector<MergeEvent> branch;
		// loop until root and fill to the branch
		for (int j = i; j >= 0; j = nodes[j].prev) {
			branch.push_back(merge_events[nodes[j].index]);
		}
		branches.push_back(branch);
	}

	// check and delete subset
	// branch is subset of other branch ?
	std::vector<bool> subset;
	for (size_t i = 0; i < branches.size(); ++i) {
		bool found_superset = false;
		for (size_t j = i+1; j < branches.size(); ++j) {
			// the superset must have more elements than the subset
			if (branches[j].size() <= branches[i].size()) continue;
			// is subset
			bool is_subset = true;
			// check nodes in branch
			for (size_t k = 0; k < branches[i].size(); ++k) {
				// search for the same node in possible superset
				// found
				bool found = false;
				for (size_t l = 0; l < branches[j].size(); ++l) {
					if (branches[i][k].flag == branches[j][l].flag) {
						found = true;
						break;
					}
				}
				// There is at least one element out of this set.
				// It can't be the superset.
				if (!found) {
					is_subset = false;
					break;
				}
			}
			// All element can be found in this set. We found superset.
			if (is_subset) {
				found_superset = true;
				break;
			}
		}
		subset.push_back(found_superset);
	}


	for (size_t i = 0; i < branches.size(); ++i) {
		if (subset[i]) continue;
		// fill to result
		DssdMergeEvent dssd_merge;
		dssd_merge.hit = 0;
		int &num = dssd_merge.hit;
		std::vector<bool> dssd_hole;
		for (int j = branches[i].size()-1; j >= 0; --j) {
			if (num >= 8) break;
			dssd_merge.flag[num] = branches[i][j].flag;
			dssd_merge.merge_tag[num] = branches[i][j].tag;
			dssd_merge.x[num] = branches[i][j].x;
			dssd_merge.y[num] = branches[i][j].y;
			dssd_merge.z[num] = 0.0;
			dssd_merge.energy[num] = branches[i][j].energy;
			// ignore time at current stage
			if (hole_info) {
				dssd_hole.push_back(branches[i][j].hole);
			}
			++num;
		}
		result.push_back(dssd_merge);
		hole.push_back(dssd_hole);
	}

	// put an empty event if nothing is found
	if (result.empty()) {
		DssdMergeEvent empty_event;
		empty_event.hit = 0;
		result.push_back(empty_event);
		std::vector<bool> empty_hole;
		hole.push_back(empty_hole);
	}


	if (print_debug) {
		std::cout << "---------- Possible merge events ----------\n"
			<< "flag, tag, x, y, energy\n";
		for (const auto &merge : merge_events) {
			std::cout << std::hex << merge.flag << ", ";
			if (merge.tag == 0) std::cout << "f1b1, ";
			else if (merge.tag == 1) std::cout << "f1b2, ";
			else if (merge.tag == 2) std::cout << "f2b1, ";
			else if (merge.tag == 3) std::cout << "f2b2, ";
			else std::cout << ", ";
			std::cout << merge.x << ", "
				<< merge.y << ", "
				<< merge.energy << "\n";
		}

		std::cout << "---------- Nodes ----------\n"
			<< "index, prev, num, flag, leaf\n";
		for (const auto &node : nodes) {
			std::cout << std::dec << node.index << ", "
				<< node.prev << ", "
				<< node.num << ", "
				<< std::hex << node.flag << ", "
				<< (node.leaf ? "Y" : "N") << "\n";
		}
	}
}


/// @brief calculate T0D1 position from strip
/// @param[in] front_strip front strip number
/// @param[in] back_strip back strip number
/// @param[out] x position x
/// @param[out] y position y
///
void T0D1PositionFromStrip(
	double front_strip,
	double back_strip,
	double &x,
	double &y
) {
	x = back_strip - 31.5;
	y = front_strip - 31.5;
}


/// @brief calculate T0D2 position from strip
/// @param[in] front_strip front strip number
/// @param[in] back_strip back strip number
/// @param[out] x position x
/// @param[out] y position y
///
void T0D2PositionFromStrip(
	double front_strip,
	double back_strip,
	double &x,
	double &y
) {
	x = 2.0 * front_strip - 31.0 - 1.03;
	y = 2.0 * back_strip - 31.0 - 0.86;
}


/// @brief calculate T0D3 position from strip
/// @param[in] front_strip front strip number
/// @param[in] back_strip back strip number
/// @param[out] x position x
/// @param[out] y position y
///
void T0D3PositionFromStrip(
	double front_strip,
	double back_strip,
	double &x,
	double &y
) {
	x = 2.0 * front_strip - 31.0 - 0.95;
	y = 2.0 * back_strip - 31.0 - 0.8;
}


/// @brief search for T0D2 f2b1 binding events
/// @param[in] d2_event T0D2 fundamental event
/// @param[in] range T0D2 front-back correlation energy range
/// @param[in] hole T0D2 hole flag
/// @param[out] d2_merge T0D2 merge events
/// @returns true if T0D2 binding events is found
///
bool SearchD2Front2Back1(
	const DssdFundamentalEvent &d2_event,
	const double range,
	const HoleInfo *hole,
	std::vector<DssdMergeEvent> &d2_merge
) {
	// result
	bool found = false;
	// variables for convenient
	const int &fhit = d2_event.front_hit;
	const int &bhit = d2_event.back_hit;
	const unsigned short *fs = d2_event.front_strip;
	const unsigned short *bs = d2_event.back_strip;
	const double *fe = d2_event.front_energy;
	const double *be = d2_event.back_energy;

	// search for f2b1 event
	for (int i = 0; i < bhit; ++i) {
		for (int j = 0; j < fhit; ++j) {
			if (hole->IsHole(fs[j], bs[i])) continue;
			for (int k = j+1; k < fhit; ++k) {
				if (hole->IsHole(fs[k], bs[i])) continue;
				if (abs(fs[j]-fs[k]) == 1) continue;
				// energy difference
				double de = be[i] - fe[j] - fe[k];
				if (fabs(de) >= range) continue;
				DssdMergeEvent merge;
				merge.hit = 2;
				merge.flag[0] = (0x100<<i) | (0x1<<j);
				merge.merge_tag[0] = 5;
				merge.energy[0] = fe[j];
				T0D2PositionFromStrip(fs[j], bs[i], merge.x[0], merge.y[0]);
				merge.flag[1] = (0x100<<i) | (0x1<<k);
				merge.merge_tag[1] = 5;
				merge.energy[1] = fe[k];
				T0D2PositionFromStrip(fs[k], bs[i], merge.x[1], merge.y[1]);
				found = true;
				d2_merge.push_back(merge);
			}
		}
	}

	return found;
}


/// @brief search for T0D2 f1b2 binding events
/// @param[in] d2_event T0D2 fundamental event
/// @param[in] range T0D2 front-back correlation energy range
/// @param[in] hole T0D2 hole flag
/// @param[out] d2_merge T0D2 merge events
/// @returns true if T0D2 binding events is found
///
bool SearchD2Front1Back2(
	const DssdFundamentalEvent &d2_event,
	const double range,
	const HoleInfo *hole,
	std::vector<DssdMergeEvent> &d2_merge
) {
	// result
	bool found = false;
	// variables for convenient
	const int &fhit = d2_event.front_hit;
	const int &bhit = d2_event.back_hit;
	const unsigned short *fs = d2_event.front_strip;
	const unsigned short *bs = d2_event.back_strip;
	const double *fe = d2_event.front_energy;
	const double *be = d2_event.back_energy;

	// search for f1b2 event
	for (int i = 0; i < fhit; ++i) {
		for (int j = 0; j < bhit; ++j) {
			if (hole->IsHole(fs[i], bs[j])) continue;
			for (int k = j+1; k < bhit; ++k) {
				if (hole->IsHole(fs[i], bs[k])) continue;
				if (abs(bs[j]-bs[k]) == 1) continue;
				// energy difference
				double de = fe[i] - be[j] - be[k];
				if (fabs(de) >= range) continue;
				DssdMergeEvent merge;
				merge.hit = 2;
				merge.flag[0] = (0x1<<i) | (0x100<<j);
				merge.merge_tag[0] = 4;
				merge.energy[0] = fe[i] * be[j] / (be[j] + be[k]);
				T0D2PositionFromStrip(fs[i], bs[j], merge.x[0], merge.y[0]);
				merge.flag[1] = (0x1<<i) | (0x100<<k);
				merge.merge_tag[1] = 4;
				merge.energy[1] = fe[i] * be[k] / (be[j] + be[k]);
				T0D2PositionFromStrip(fs[i], bs[k], merge.x[1], merge.y[1]);
				found = true;
				d2_merge.push_back(merge);
			}
		}
	}

	return found;
}


/// @brief get binding events from T0D1 and T0D2 fundamental events
/// @param[in] d1_event T0D1 fundamental event
/// @param[in] d2_event T0D2 fundamental event
/// @param[in] range1 T0D1 front-back energy correlation range
/// @param[in] range2 T0D2 front-back energy correlation range
/// @param[in] hole T0D2 hole flag
/// @param[out] d1_merge1 T0D1 f1b2 binding or f2ab2 merge events
/// @param[out] d1_merge2 T0D2 f2b1 binding merge events
/// @param[out] d2_merge1 T0D1 f2b1 binding or f2b2a merge events
/// @param[out] d2_merge2 T0D2 f1b2 binding merge events
///
void GetBindEvents(
	const DssdFundamentalEvent &d1_event,
	const DssdFundamentalEvent &d2_event,
	const double range1,
	const double range2,
	const HoleInfo *hole,
	std::vector<DssdMergeEvent> &d1_merge1,
	std::vector<DssdMergeEvent> &d1_merge2,
	std::vector<DssdMergeEvent> &d2_merge1,
	std::vector<DssdMergeEvent> &d2_merge2
) {
	// for convenient
	int d1_fhit = d1_event.front_hit;
	int d1_bhit = d1_event.back_hit;
	const double *d1fe = d1_event.front_energy;
	const double *d1be = d1_event.back_energy;
	const unsigned short *d1fs = d1_event.front_strip;
	const unsigned short *d1bs = d1_event.back_strip;

	// search for D1 f1b2 event
	bool found_d1_f1b2 = false;
	for (int i = 0; i < d1_fhit; ++i) {
		if (d1fe[i] < 200.0) break;
		if (d1fe[i] < 15000.0) continue;
		// search for back twin events in energy range
		for (int j = 0; j < d1_bhit; ++j) {
			if (d1be[j] < 10000.0) continue;
			for (int k = j+1; k < d1_bhit; ++k) {
				if (d1be[k] > 10000.0) continue;
				// energy difference
				double de = d1fe[i] - d1be[j] - d1be[k];
				if (fabs(de) >= range1) continue;
				// found
				DssdMergeEvent merge;
				merge.case_tag = 0;
				merge.hit = 2;
				merge.flag[0] = (0x1<<i) | (0x100<<j);
				merge.merge_tag[0] = 4;
				merge.energy[0] = d1fe[i] * d1be[j] / (d1be[j]+d1be[k]);
				T0D1PositionFromStrip(
					d1fs[i], d1bs[j], merge.x[0], merge.y[0]
				);
				merge.flag[1] = (0x1<<i) | (0x100<<k);
				merge.merge_tag[1] = 4;
				merge.energy[1] = d1fe[i] * d1be[k] / (d1be[j]+d1be[k]);
				T0D1PositionFromStrip(
					d1fs[i], d1bs[k], merge.x[1], merge.y[1]
				);
				merge.z[0] = merge.z[1] = 100.0;
				merge.time[0] = d1_event.front_time[i];
				merge.time[1] = d1_event.front_time[i];
				merge.time_flag[0] = merge.time_flag[1] = 0;
				d1_merge1.push_back(merge);
				found_d1_f1b2 = true;
			}
		}
	}

	// search for D1 f2ab2 event
	bool found_d1_f2ab2 = false;
	for (int i = 0; i < d1_fhit; ++i) {
		if (d1fe[i] < 10000.0) break;
		for (int j = i+1; j < d1_fhit; ++j) {
			if (d1fe[j] > 10000.0) continue;
			if (abs(d1fs[i]-d1fs[j]) != 1) continue;
			// search for back events
			for (int k = 0; k < d1_bhit; ++k) {
				double de1 = d1fe[i] - d1be[k];
				if (fabs(de1) >= range1) continue;
				for (int l = k+1; l < d1_bhit; ++l) {
					double de2 = d1fe[j] - d1be[l];
					if (fabs(de2) >= range1) continue;
					// found
					DssdMergeEvent merge;
					merge.hit = 2;
					merge.flag[0] = (0x1<<i) | (0x100<<k);
					merge.merge_tag[0] = 0;
					merge.energy[0] = d1fe[i];
					T0D1PositionFromStrip(
						d1fs[i], d1bs[k], merge.x[0], merge.y[0]
					);
					merge.flag[1] = (0x1<<j) | (0x100<<l);
					merge.merge_tag[1] = 0;
					merge.energy[1] = d1fe[j];
					T0D1PositionFromStrip(
						d1fs[j], d1bs[l], merge.x[1], merge.y[1]
					);
					merge.z[0] = merge.z[1] = 100.0;
					merge.time[0] = d1_event.front_time[i];
					merge.time[1] = d1_event.front_time[j];
					merge.time_flag[0] = merge.time_flag[1] = 0;
					d1_merge1.push_back(merge);
					found_d1_f2ab2 = true;
				}
			}
		}
	}

	if (found_d1_f1b2 || found_d1_f2ab2) {
		bool found_d2_f2b1 =
			SearchD2Front2Back1(d2_event, range2, hole, d2_merge1);
		if (!found_d2_f2b1) {
			d1_merge1.clear();
			d2_merge1.clear();
		}
	}


	// search for D1 f2b1 event
	bool found_d1_f2b1 = false;
	for (int i = 0; i < d1_bhit; ++i) {
		if (d1be[i] < 15000.0) continue;
		// search for front twin events in energy range
		for (int j = 0; j < d1_fhit; ++j) {
			if (d1fe[j] < 10000.0) continue;
			for (int k = j+1; k < d1_fhit; ++k) {
				if (d1fe[k] > 10000.0) continue;
				// energy difference
				double de = d1be[i] - d1fe[j] - d1fe[k];
				if (fabs(de) >= range1) continue;
				// found
				DssdMergeEvent merge;
				merge.hit = 2;
				merge.flag[0] = (0x100<<i) | (0x1<<j);
				merge.merge_tag[0] = 5;
				merge.energy[0] = d1fe[j];
				T0D1PositionFromStrip(
					d1fs[j], d1bs[i], merge.x[0], merge.y[0]
				);
				merge.flag[1] = (0x100<<i) | (0x1<<k);
				merge.merge_tag[1] = 5;
				merge.energy[1] = d1fe[k];
				T0D1PositionFromStrip(
					d1fs[k], d1bs[i], merge.x[1], merge.y[1]
				);
				merge.z[0] = merge.z[1] = 100.0;
				merge.time[0] = d1_event.front_time[j];
				merge.time[1] = d1_event.front_time[k];
				merge.time_flag[0] = merge.time_flag[1] = 0;
				d1_merge2.push_back(merge);
				found_d1_f2b1 = true;
			}
		}
	}

	// search for D1 f2b2a event
	bool found_d1_f2b2a = false;
	for (int i = 0; i < d1_bhit; ++i) {
		if (d1be[i] < 10000.0) break;
		for (int j = i+1; j < d1_bhit; ++j) {
			if (d1be[j] > 10000.0) continue;
			if (abs(d1bs[i]-d1bs[j]) != 1) continue;
			// search for front events
			for (int k = 0; k < d1_fhit; ++k) {
				double de1 = d1be[i] - d1fe[k];
				if (fabs(de1) >= range1) continue;
				for (int l = k+1; l < d1_fhit; ++l) {
					double de2 = d1be[j] - d1fe[l];
					if (fabs(de2) >= range1) continue;
					// found
					DssdMergeEvent merge;
					merge.hit = 2;
					merge.flag[0] = (0x100<<i) | (0x1<<k);
					merge.merge_tag[0] = 0;
					merge.energy[0] = d1fe[k];
					T0D1PositionFromStrip(
						d1fs[k], d1bs[i], merge.x[0], merge.y[0]
					);
					merge.flag[1] = (0x100<<j) | (0x1<<l);
					merge.merge_tag[1] = 0;
					merge.energy[1] = d1fe[l];
					T0D1PositionFromStrip(
						d1fs[l], d1bs[j], merge.x[1], merge.y[1]
					);
					merge.z[0] = merge.z[1] = 100.0;
					merge.time[0] = d1_event.front_time[k];
					merge.time[1] = d1_event.front_time[l];
					merge.time_flag[0] = merge.time_flag[1] = 0;
					d1_merge2.push_back(merge);
					found_d1_f2b2a = true;
				}
			}
		}
	}

	if (found_d1_f2b1 || found_d1_f2b2a) {
		bool found_d2_f1b2 =
			SearchD2Front1Back2(d2_event, range2, hole, d2_merge2);
		if (!found_d2_f1b2) {
			d1_merge2.clear();
			d2_merge2.clear();
		}
	}
}


struct ParticleCut {
	unsigned short charge;
	unsigned short mass;
	std::unique_ptr<TCutG> cut;
};


std::unique_ptr<TCutG> ReadCut(const char *name) {
	// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name
	);
	// open cut file to read points
	std::ifstream fin(cut_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< cut_file_name << " failed.\n";
		return nullptr;
	}
	// result
	std::unique_ptr<TCutG> result = std::make_unique<TCutG>();
	// point index
	int point;
	// point positions
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		result->SetPoint(point, x, y);
	}
	// close file
	fin.close();
	return result;
}


struct T0Cut {
	// T0D1-D2 cuts
	std::vector<ParticleCut> d1d2_cuts;
	// T0D1-D2 tail cuts
	std::vector<ParticleCut> d1d2_tail_cuts;
	// T0D2-D3 cuts
	std::vector<ParticleCut> d2d3_cuts;
	// T0D2-D3 tail cuts
	std::vector<ParticleCut> d2d3_tail_cuts;
	// T0D3-S1 cuts
	std::vector<ParticleCut> d3s1_cuts;
	// T0D3-S1 tail cuts
	std::vector<ParticleCut> d3s1_tail_cuts;
	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_cuts;
	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_tail_cuts;
	// T0S2-S3 cuts
	std::vector<ParticleCut> s2s3_cuts;
	// T0S2-S3 pass cuts
	std::vector<ParticleCut> s2s3_tail_cuts;
	// T0D1-D2 hole cuts
	std::vector<ParticleCut> d1d2_hole_cuts;
	// T0D1-D2 hole tail cuts
	std::vector<ParticleCut> d1d2_hole_tail_cuts;
	// T0D2-D3 hole cuts
	std::vector<ParticleCut> d2d3_hole_cuts;
	// T0D2-D3 hole tail cuts
	std::vector<ParticleCut> d2d3_hole_tail_cuts;
};


class Slice {
public:
	unsigned short charge;
	unsigned short mass;
	bool tail;
	bool straight;
	int index[2];

	Slice(
		unsigned short c,
		unsigned short m,
		bool t,
		bool s,
		int i1,
		int i2
	)
	: charge(c)
	, mass(m)
	, tail(t)
	, straight(s)
	{
		index[0] = i1;
		index[1] = i2;
	}

	inline bool Connectable(const Slice &next) {
		return (
			charge == next.charge
			&& (mass == next.mass || mass == 0)
			&& index[1] == next.index[0]
			&& tail
		);
	}
};


struct StraightInfo {
	unsigned short charge;
	unsigned short mass;
	double min;
	double max;
};


/// @brief build DSSD-DSSD slice with hole, e.g. T0D1-D2, T0D2-D3
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[in] hole_cuts PID cut for hole area and particle stopped
/// @param[in] hole_tail_cuts PID cut for hole area and particle passed
/// @param[in] hole1 true the first is hole
/// @param[in] hole2 true the second is hole
/// @param[in] offset durable position offset for two layers
/// @param[out] slices built slices
///
void BuildDssdSliceHole(
	const DssdMergeEvent &layer1,
	const DssdMergeEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	const std::vector<ParticleCut> &hole_cuts,
	const std::vector<ParticleCut> &hole_tail_cuts,
	const bool *hole1,
	const bool *hole2,
	const double sa,
	const double sb,
	const std::map<int, StraightInfo> &straight_info,
	double offset,
	std::vector<Slice> &slices
) {
	for (int i = 0; i < layer1.hit; ++i) {
		for (int j = 0; j < layer2.hit; ++j) {
			double xoffset = layer1.x[i] - layer2.x[j];
			double yoffset = layer1.y[i] - layer2.y[j];
			// check x, y offset
			if (fabs(xoffset) < offset && fabs(yoffset) < offset) {
				const std::vector<ParticleCut> *used_cuts = nullptr;
				const std::vector<ParticleCut> *used_tail_cuts = nullptr;
				if ((hole1 && hole1[i]) || (hole2 && hole2[j])) {
					used_cuts = &hole_cuts;
					used_tail_cuts = &hole_tail_cuts;
				} else if ((hole1 && !hole1[i]) || (hole2 && !hole2[j])) {
					used_cuts = &cuts;
					used_tail_cuts = &tail_cuts;
				} else {
					std::cerr << "Error: hole1 and hole2 invalid.\n";
					return;
				}
				// check stopped PID curv
				for (const auto &cut : *(used_cuts)) {
					if (cut.cut->IsInside(
						layer2.energy[j], layer1.energy[i]
					)) {
						double ef = Straight(
							layer1.energy[i], layer2.energy[j], sa, sb
						);
						auto search = straight_info.find(cut.charge*100+cut.mass);
						bool is_straight = false;
						if (
							search != straight_info.end()
							&& ef > search->second.min
							&& ef < search->second.max
						) {
							is_straight = true;
						}

						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							false,
							is_straight,
							i, j
						);
					}
				}
				// check pass PID curv
				for (const auto &cut : *(used_tail_cuts)) {
					if (cut.cut->IsInside(
						layer2.energy[j], layer1.energy[i]
					)) {
						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							true, false,
							i, j
						);
					}
				}
			}
		}
	}
}


/// @brief build DSSD-SSD slice, e.g. T0D3-S1
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[out] slices built slices
///
void BuildDssdSsdSlice(
	const DssdMergeEvent &layer1,
	const SsdEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	std::vector<Slice> &slices
) {
	// check valid, although energy is 0 if is invalid
	if (layer2.time < -9e4) return;
	for (int i = 0; i < layer1.hit; ++i) {
		for (const auto &cut : cuts) {
			if (cut.cut->IsInside(layer2.energy, layer1.energy[i])) {
				// build slice
				slices.emplace_back(
					cut.charge,
					cut.mass,
					false, false,
					i, 0
				);
			}
		}
		for (const auto &cut : tail_cuts) {
			if (cut.cut->IsInside(layer2.energy, layer1.energy[i])) {
				// build slice
				slices.emplace_back(
					cut.charge,
					cut.mass,
					true, false,
					i, 0
				);
			}
		}
	}
	return;
}


/// @brief build SSD-SSD slice, e.g. T0S1-S2, T0S2-S3
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[out] slices built slices
///
void BuildSsdSlice(
	const SsdEvent &layer1,
	const SsdEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	std::vector<Slice> &slices
) {
	// check valid, although energy is 0 if is invalid
	if (layer1.time < -9e4 || layer2.time < -9e4) return;
	for (const auto &cut : cuts) {
		if (cut.cut->IsInside(layer2.energy, layer1.energy)) {
			slices.emplace_back(
				cut.charge, cut.mass,
				false, false,
				0, 0
			);
		}
	}
	for (const auto &cut : tail_cuts) {
		if (cut.cut->IsInside(layer2.energy, layer1.energy)) {
			slices.emplace_back(
				cut.charge, cut.mass,
				true, false,
				0, 0
			);
		}
	}
}


/// @brief build slices for all layers
/// @param[in] d1 T0D1 merged event
/// @param[in] d2 T0D2 merged event
/// @param[in] d3 T0D3 merged event
/// @param[in] s1 T0S1 event
/// @param[in] s2 T0S2 event
/// @param[in] s3 T0S3 event
/// @param[in] t0_cut PID cut of T0 all layers
/// @param[out] slices built slices
///
void BuildSlice(
	const DssdMergeEvent &d1,
	const DssdMergeEvent &d2,
	const DssdMergeEvent &d3,
	const SsdEvent &s1,
	const SsdEvent &s2,
	const SsdEvent &s3,
	const std::vector<bool> &hole,
	const T0Cut &t0_cut,
	std::vector<Slice> *slices
) {
	bool d2_hole[8];
	for (size_t i = 0; i < hole.size(); ++i) d2_hole[i] = hole[i];
	std::map<int, StraightInfo> d1d2_straight_info;
	d1d2_straight_info.insert(
		std::make_pair(410, StraightInfo{4, 10, 20200, 21200})
	);
	d1d2_straight_info.insert(
		std::make_pair(204, StraightInfo{2, 4, 6100, 6600})
	);
	BuildDssdSliceHole(
		d1, d2,
		t0_cut.d1d2_cuts, t0_cut.d1d2_tail_cuts,
		t0_cut.d1d2_hole_cuts, t0_cut.d1d2_hole_tail_cuts,
		nullptr, d2_hole,
		sa12, sb12,
		d1d2_straight_info,
		4.0,
		slices[0]
	);

	std::map<int, StraightInfo> d2d3_straight_info;
	d2d3_straight_info.insert(
		std::make_pair(410, StraightInfo{4, 10, 24000, 25000})
	);
	d2d3_straight_info.insert(
		std::make_pair(204, StraightInfo{2, 4, 7000, 8000})
	);
	BuildDssdSliceHole(
		d2, d3,
		t0_cut.d2d3_cuts, t0_cut.d2d3_tail_cuts,
		t0_cut.d2d3_hole_cuts, t0_cut.d2d3_hole_tail_cuts,
		d2_hole, nullptr,
		sa23, sb23,
		d2d3_straight_info,
		4.2,
		slices[1]
	);
	BuildDssdSsdSlice(
		d3, s1, t0_cut.d3s1_cuts, t0_cut.d3s1_tail_cuts, slices[2]
	);
	BuildSsdSlice(
		s1, s2, t0_cut.s1s2_cuts, t0_cut.s1s2_tail_cuts, slices[3]
	);
	BuildSsdSlice(
		s2, s3, t0_cut.s2s3_cuts, t0_cut.s2s3_tail_cuts, slices[4]
	);
}


struct SliceNode {
	// layer 0-T0D2, 1-T0D3, 2-T0S1, 3-T0S2, 4-T0S3, 5-CsI
	int layer;
	// index in slices
	int index;
	// parent node, 0 for root
	int parent;
	// flag of merge event, includes all layers' information
	// prevent to use events repeatedly
	unsigned int flag;
	// particle type get from cut, mass*10+charge, e.g. 4He is 42, 10Be is 104
	int type;
	// true if it's last layer is not in tail
	bool leaf;
};


/// @brief convert slices into tree structue
/// @param[in] slices slice nodes in two dimensional array
/// @param[in] nodes slice nodes store in vector, but connected by pointer like tree
///
void BuildSliceTree(
	std::vector<Slice> *slices,
	std::vector<SliceNode> &nodes
) {
	nodes.clear();
	// fill root node
	nodes.push_back(SliceNode{-1, -1, 0, 0, 0, false});
	// start index of each layer in nodes
	int start[6] = {
		1, 1, 1, 1, 1, 1
	};
	// offset of each layer in flag
	unsigned int offsets[6] = {
		0, 8, 16, 24, 25, 26
	};
	// fill first layer
	for (int i = 0; i < int(slices[0].size()); ++i) {
		// T0D1 flag, with offset 0
		unsigned int flag = 0x1 << slices[0][i].index[0];
		// T0D2 flag, with offset 8
		flag |= 0x100 << slices[0][i].index[1];
		// particle type
		int type = slices[0][i].mass * 10 + slices[0][i].charge;
		// insert new node
		nodes.push_back(SliceNode{
			0, i, 0, flag, type, !(slices[0][i].tail)
		});
		// increase start index
		++start[1];
	}


	for (int layer = 1; layer < 5; ++layer) {
		// intialize the next layer start value
		start[layer+1] = start[layer];
		// loop slice in this layer
		for (int i = 0; i < int(slices[layer].size()); ++i) {
			Slice &current_layer_slice = slices[layer][i];
			// search for connectable slice in nodes
			for (int j = start[layer-1]; j < start[layer]; ++j) {
				if (print_debug) {
					std::cout << "layer " << layer
						<< ", slice index " << i
						<< ", node index " << j << "\n";
				}
				const SliceNode &node = nodes[j];
				Slice &last_layer_slice = slices[layer-1][node.index];
				if (last_layer_slice.Connectable(current_layer_slice)) {
					// flag of current layer
					unsigned int flag = 0x1 << (
						current_layer_slice.index[1] + offsets[layer+1]
					);
					// or with the previous layers
					flag |= node.flag;
					// particle type
					int type = current_layer_slice.mass*10
						+ current_layer_slice.charge;
					if (layer == 4 && current_layer_slice.tail) {
						// insert node that stop in CsI
						nodes.push_back(SliceNode{
							5, i, j, flag, type, true
						});
					} else {
						// insert new node
						nodes.push_back(SliceNode{
							layer, i, j, flag, type, !(current_layer_slice.tail)
						});
					}
					// increase start index
					++start[layer+1];
				}
			}
		}
	}

	if (print_debug) {
		std::cout << "Slice:\n";
		std::cout << "Charge, Mass, Penetrate, Index1, Index2\n";
		for (size_t i = 0; i < 5; ++i) {
			std::cout << "Layer " << i << "\n";
			for (size_t j = 0; j < slices[i].size(); ++j) {
				std::cout << slices[i][j].charge << ", "
					<< slices[i][j].mass << ", "
					<< slices[i][j].tail << ", "
					<< slices[i][j].index[0] << ", "
					<< slices[i][j].index[1] << "\n";
			}
		}
		std::cout << "Node:\n";
		std::cout << ", Layer, Index, Parent, Flag, Leaf\n";
		for (size_t i = 0; i < nodes.size(); ++i) {
			std::cout << i << " " << nodes[i].layer << ", "
				<< nodes[i].index << ", " << nodes[i].parent
				<< ", " << std::hex << nodes[i].flag << ", "
				<< std::dec << nodes[i].leaf << "\n";
		}
	}

	return;
}


struct LeaveGroup {
	int prev;
	int index;
	int count;
};


/// @brief check whether flag is contraditory to the existing group
/// @param[in] leaves all leave nodes
/// @param[in] groups all current groups
/// @param[in] group_index index of group to check
/// @param[in] leave_flag flag of leaf to checks
/// @returns true if flag is valid, false if it's contraditory
///
bool CheckSliceFlag(
	const std::vector<SliceNode> &leaves,
	const std::vector<LeaveGroup> &groups,
	int group_index,
	unsigned int leaf_flag
) {
	for (int i = group_index; i != 0; i = groups[i].prev) {
		if ((leaves[groups[i].index].flag & leaf_flag) != 0) {
			return false;
		}
	}
	return true;
}


/// @brief Choose a group of leaves(particles)
/// @param[in] nodes all slices
/// @param[out] group selected group
/// @param[out] behe_found found 10Be and 4He in this group
/// @param[out] used_slices number of slices used in this group
///
void PickSliceGroup(
	const std::vector<SliceNode> &nodes,
	std::vector<SliceNode> &group,
	bool &behe_found,
	int &used_slices
) {
	// leaves
	std::vector<SliceNode> leaves;
	// get leaves
	for (const auto &node : nodes) {
		if (node.leaf) leaves.push_back(node);
	}

	// groups
	std::vector<LeaveGroup> groups;
	// input the root slot
	groups.push_back(LeaveGroup{-1, -1, 0});

	// insert the first leave in groups
	for (int i = 0; i < int(leaves.size()); ++i) {
		groups.push_back(LeaveGroup{0, i, 1});
	}

	// start index to add more leaves
	int start = 1;
	// tail index to add more leaves, exclusive
	int tail = groups.size();
	// combine the leaves into groups
	// if start is larger than or equal to tail, no new node has been add
	// to the groups, so, terminate the loop
	while (start < tail) {
		for (int i = start; i < tail; ++i) {
			for (int j = groups[i].index+1; j < int(leaves.size()); ++j) {
				if (CheckSliceFlag(leaves, groups, i, leaves[j].flag)) {
					// add new node to group
					groups.push_back(LeaveGroup{i, j, groups[i].count+1});
				}
			}
		}
		start = tail;
		tail = groups.size();
	}

	// // search for the group with maximum count
	// // max count
	// int max_count = 0;
	// // index with maximum count
	// int found_index = 0;
	// // loop to search
	// for (int i = 0; i < int(groups.size()); ++i) {
	// 	if (groups[i].count > max_count) {
	// 		max_count = groups[i].count;
	// 		found_index = i;
	// 	}
	// }

	behe_found = false;
	used_slices = 0;
	int max_count = 0;
	int found_index = 0;
	for (int i = 1; i < int(groups.size()); ++i) {
		int behe_flag = 0;
		int slice_num = 0;
		for (auto g = groups[i]; g.index != -1; g = groups[g.prev]) {
			const auto &leave = leaves[g.index];
			if (leave.type == 104) behe_flag |= 1;
			else if (leave.type == 42) behe_flag |= 2;
			slice_num += leave.layer+1;
			if (print_debug) {
				std::cout << "-----------------------------\n"
					<< "group " << i << ": " << g.index
					<< ", prev " << g.prev << "\n"
					<< "type " << leave.type
					<< ", flag " << behe_flag << "\n"
					<< "slice number " << slice_num << "\n";
			}
		}
		if (behe_flag == 3) {
			if (!behe_found) {
				behe_found = true;
				used_slices = slice_num;
				found_index = i;
			} else {
				if (groups[i].count > max_count) {
					max_count = groups[i].count;
					used_slices = slice_num;
					found_index = i;
				} else if (
					groups[i].count == max_count
					&& slice_num > used_slices
				) {
					used_slices = slice_num;
					found_index = i;
				}
			}
		} else {
			if (
				(
					!behe_found
					&& groups[i].count > max_count
				) || (
					!behe_found
					&& groups[i].count == max_count
					&& slice_num > used_slices
				)
			) {
				max_count = groups[i].count;
				used_slices = slice_num;
				found_index = i;
			}
		}
		if (print_debug) {
			std::cout << "found index " << found_index
				<< ", behe " << (behe_found ? "Y" : "N")
				<< ", used slice " << used_slices << "\n";
		}
	}

	group.clear();
	for (int i = found_index; i != 0; i = groups[i].prev) {
		group.push_back(leaves[groups[i].index]);
	}
}


struct ExtraInfo {
	int run;
	long long entry;
};


/// @brief track T0 events in slice method
/// @param[in] d1_events possible T0D1 events
/// @param[in] d2_events possible T0D2 events
/// @param[in] d3_events possible T0D3 events
/// @param[in] s1_event T0S1 event
/// @param[in] s2_event T0S2 event
/// @param[in] s3_event T0S3 event
/// @param[in] d2_hole T0D2 hole flags
/// @param[in] t0_cut T0 cuts
/// @param[in] info extra information
/// @param[out] t0_event tracked T0 event
/// @param[out] merge_event selecteed merged event
/// @param[out] found_behe_times times Be+He is found
/// @returns 0 if success, -1 otherwise
///
int SliceTrack(
	const std::vector<DssdMergeEvent> &d1_events,
	const std::vector<DssdMergeEvent> &d2_events,
	const std::vector<DssdMergeEvent> &d3_events,
	const SsdEvent &s1_event,
	const SsdEvent &s2_event,
	const SsdEvent &s3_event,
	const std::vector<std::vector<bool>> &d2_hole,
	const std::vector<DssdMergeEvent> &d1_bind_events1,
	const std::vector<DssdMergeEvent> &d2_bind_events1,
	const std::vector<DssdMergeEvent> &d1_bind_events2,
	const std::vector<DssdMergeEvent> &d2_bind_events2,
	const T0Cut &t0_cut,
	const ExtraInfo &info,
	T0Event &t0_event,
	DssdMergeEvent *merge_event,
	int &found_behe_times
) {
	// picked slices
	std::vector<Slice> slices[5];
	// picked nodes
	std::vector<SliceNode> nodes;
	// picked nodess group
	std::vector<SliceNode> group;
	// found Be and He times
	found_behe_times = 0;
	// found particles of selected group
	int found_particles = 0;
	// used slices of selected group
	int used_slices = 0;
	// selected events
	// const DssdMergeEvent *selected_d1_event = nullptr;
	// const DssdMergeEvent *selected_d2_event = nullptr;
	// const DssdMergeEvent *selected_d3_event = nullptr;
	const std::vector<bool> *selected_d2_hole = nullptr;
	// loop possible events
	for (size_t i = 0; i < d1_events.size(); ++i) {
		for (size_t j = 0; j < d2_events.size(); ++j) {
			for (size_t k = 0; k < d3_events.size(); ++k) {
				// used events
				const DssdMergeEvent *used_d1_event = &(d1_events[i]);
				const DssdMergeEvent *used_d2_event = &(d2_events[j]);
				const DssdMergeEvent *used_d3_event = &(d3_events[k]);

				// slices without supplementary
				std::vector<Slice> normal_slices[5];
				// nodes without supplementary
				std::vector<SliceNode> normal_nodes;
				// temporary picked nodes group
				std::vector<SliceNode> normal_group;
				// found 10Be and 4He in normal group ?
				bool normal_found_behe;
				// number of used slices in normal group
				int normal_used_slices;
				// build slice from DSSD and SSD events
				BuildSlice(
					*used_d1_event, *used_d2_event, *used_d3_event,
					s1_event, s2_event, s3_event,
					d2_hole[j],
					t0_cut,
					normal_slices
				);
				// convert to tree structure
				BuildSliceTree(normal_slices, normal_nodes);
				// pick one group
				PickSliceGroup(
					normal_nodes, normal_group,
					normal_found_behe, normal_used_slices
				);

				// this group is selected?
				bool selected = false;
				// select it if it's the best group
				if (found_behe_times == 0 && normal_found_behe) {
					selected = true;
				} else if (found_behe_times == 0 || normal_found_behe) {
					if (int(normal_group.size()) > found_particles) selected = true;
					else if (int(normal_group.size()) == found_particles) {
						if (normal_used_slices > used_slices) selected = true;
					}
				}
				// increase found Be He times
				if (normal_found_behe) ++found_behe_times;
				// record the selected slices
				if (selected) {
					found_particles = normal_group.size();
					used_slices = normal_used_slices;
					for (int l = 0; l < 5; ++l) slices[l] = normal_slices[l];
					nodes = normal_nodes;
					group = normal_group;
					merge_event[0] = d1_events[i];
					merge_event[1] = d2_events[j];
					merge_event[2] = d3_events[k];
					selected_d2_hole = &(d2_hole[j]);
				}
				// std::cout << "Select: " << (selected ? "Y" : "N")
				// 	<< ", BeHe " << (normal_found_behe ? "Y" : "N")
				// 	<< ", number " << normal_group.size()
				// 	<< ", slices " << normal_used_slices << "\n";
			}
		}
	}

	// make T0D2 binding events 1 hole
	std::vector<std::vector<bool>> d2_bind_hole1;
	for (size_t i = 0; i < d2_bind_events1.size(); ++i) {
		std::vector<bool> temp;
		for (int j = 0; j < d2_bind_events1[i].hit; ++j) {
			temp.push_back(false);
		}
		d2_bind_hole1.push_back(temp);
	}
	// loop possible events
	for (size_t i = 0; i < d1_bind_events1.size(); ++i) {
		if (found_behe_times > 0) break;
		for (size_t j = 0; j < d2_bind_events1.size(); ++j) {
			for (size_t k = 0; k < d3_events.size(); ++k) {
				// used events
				const DssdMergeEvent *used_d1_event = &(d1_bind_events1[i]);
				const DssdMergeEvent *used_d2_event = &(d2_bind_events1[j]);
				const DssdMergeEvent *used_d3_event = &(d3_events[k]);

				// slices
				std::vector<Slice> bind_slices[5];
				// nodes
				std::vector<SliceNode> bind_nodes;
				// picked nodes group
				std::vector<SliceNode> bind_group;
				// found 10Be and 4He in group ?
				bool bind_found_behe;
				// number of used slices
				int bind_used_slices;
				// build slice from DSSD and SSD events
				BuildSlice(
					*used_d1_event, *used_d2_event, *used_d3_event,
					s1_event, s2_event, s3_event,
					d2_bind_hole1[j],
					t0_cut,
					bind_slices
				);
				// convert to tree structure
				BuildSliceTree(bind_slices, bind_nodes);
				// pick one group
				PickSliceGroup(
					bind_nodes, bind_group,
					bind_found_behe, bind_used_slices
				);

				// this group is selected?
				bool selected = false;
				// select it if it's the best group
				if (found_behe_times == 0 && bind_found_behe) {
					selected = true;
				}
				if (bind_found_behe) ++found_behe_times;
				// record the selected slices
				if (selected) {
					found_particles = bind_group.size();
					used_slices = bind_used_slices;
					for (int l = 0; l < 5; ++l) slices[l] = bind_slices[l];
					nodes = bind_nodes;
					group = bind_group;
					merge_event[0] = d1_bind_events1[i];
					merge_event[1] = d2_bind_events1[j];
					merge_event[2] = d3_events[k];
					selected_d2_hole = &(d2_bind_hole1[j]);
				}
			}
		}
	}

	// make T0D2 binding events 1 hole
	std::vector<std::vector<bool>> d2_bind_hole2;
	for (size_t i = 0; i < d2_bind_events2.size(); ++i) {
		std::vector<bool> temp;
		for (int j = 0; j < d2_bind_events2[i].hit; ++j) {
			temp.push_back(false);
		}
		d2_bind_hole2.push_back(temp);
	}
	// loop possible events
	for (size_t i = 0; i < d1_bind_events2.size(); ++i) {
		if (found_behe_times > 0) break;
		for (size_t j = 0; j < d2_bind_events2.size(); ++j) {
			for (size_t k = 0; k < d3_events.size(); ++k) {
				// used events
				const DssdMergeEvent *used_d1_event = &(d1_bind_events2[i]);
				const DssdMergeEvent *used_d2_event = &(d2_bind_events2[j]);
				const DssdMergeEvent *used_d3_event = &(d3_events[k]);

				// slices
				std::vector<Slice> bind_slices[5];
				// nodes
				std::vector<SliceNode> bind_nodes;
				// picked nodes group
				std::vector<SliceNode> bind_group;
				// found 10Be and 4He in group ?
				bool bind_found_behe;
				// number of used slices
				int bind_used_slices;
				// build slice from DSSD and SSD events
				BuildSlice(
					*used_d1_event, *used_d2_event, *used_d3_event,
					s1_event, s2_event, s3_event,
					d2_bind_hole2[j],
					t0_cut,
					bind_slices
				);
				// convert to tree structure
				BuildSliceTree(bind_slices, bind_nodes);
				// pick one group
				PickSliceGroup(
					bind_nodes, bind_group,
					bind_found_behe, bind_used_slices
				);

				// this group is selected?
				bool selected = false;
				// select it if it's the best group
				if (found_behe_times == 0 && bind_found_behe) {
					selected = true;
				}
				if (bind_found_behe) ++found_behe_times;
				// record the selected slices
				if (selected) {
					found_particles = bind_group.size();
					used_slices = bind_used_slices;
					for (int l = 0; l < 5; ++l) slices[l] = bind_slices[l];
					nodes = bind_nodes;
					group = bind_group;
					merge_event[0] = d1_bind_events2[i];
					merge_event[1] = d2_bind_events2[j];
					merge_event[2] = d3_events[k];
					selected_d2_hole = &(d2_bind_hole2[j]);
				}
			}
		}
	}

	if (print_debug) {
		std::cout << "---------- Slice Track ----------\n";
		std::cout << "Slices:\n";
		for (int i = 0; i < 5; ++i) {
			for (const auto &slice : slices[i]) {
				std::cout << "Slice layer " << i
					<< ": Z = " << slice.charge << ", A = " << slice.mass
					<< ", tail = " << slice.tail << ", index = " << slice.index[0]
					<< ", " << slice.index[1] << "\n";
			}
		}
		std::cout << "Slice tree:\n";
		for (const auto &node : nodes) {
			std::cout << "layer " << node.layer << ", index = " << node.index
				<< ", parent " << node.parent << ", flag = " << std::hex << node.flag
				<< std::dec << ", type = " << node.type << ", leaf = "
				<< node.leaf << "\n";
		}
		std::cout << "Slice group:\n";
		for (const auto &node : group) {
			std::cout << "layer " << node.layer << ", index = " << node.index
				<< ", parent " << node.parent << ", flag = " << std::hex << node.flag
				<< std::dec << ", type = " << node.type << ", leaf = "
				<< node.leaf << "\n";
		}
	}


	// fill output event
	t0_event.num = int(group.size());
	// if (t0_event.num > 8) t0_event.num = 8;
	for (int i = 0; i < t0_event.num; ++i) {
		int detect_layer = group[i].layer == 5 ? 4 : group[i].layer;
		if (group[i].index >= int(slices[detect_layer].size())) {
			std::cerr
				<< "Error: group[i] index over range in run " << info.run
				<< ", entry " << info.entry
				<< ", layer " << detect_layer
				<< ", index is " << group[i].index
				<< ", range " << slices[detect_layer].size()
				<< " .\n";
			return -1;
		}
		Slice *slice = &(slices[detect_layer][group[i].index]);

		t0_event.layer[i] = group[i].layer + 1;
		t0_event.flag[i] = group[i].layer == 0 ? 0x3 : 0x7;
		t0_event.ssd_flag = group[i].flag >> 24;
		t0_event.charge[i] = slice->charge;
		t0_event.mass[i] = slice->mass;
		t0_event.straight[i] = true;

		SliceNode *node = &(group[i]);
		for (int layer = detect_layer; layer >= 0; --layer) {
			if (node->index >= int(slices[layer].size())) {
				std::cerr
					<< "Error: node index over range in run " << info.run
					<< ", entry " << info.entry
					<< ", layer " << layer
					<< ", index is " << node->index
					<< ", range " << slices[layer].size()
					<< " .\n";
				return -1;
			}
			slice = &(slices[layer][node->index]);
			if (layer == 0) {
				// fill T0D1 information
				int d1_index = slice->index[0];
				if (d1_index >= merge_event[0].hit) {
					std::cerr << "Error: D1 index over hit in run "
						<< info.run << ", entry "
						<< info.entry << ", index " << d1_index
						<< ", hit " << merge_event[0].hit << "\n";
					return -1;
				}
				t0_event.energy[i][0] = merge_event[0].energy[d1_index];
				t0_event.time[i][0] = merge_event[0].time[d1_index];
				t0_event.x[i][0] = merge_event[0].x[d1_index];
				t0_event.y[i][0] = merge_event[0].y[d1_index];
				t0_event.z[i][0] = 100.0;
				t0_event.dssd_flag[i][0] = merge_event[0].flag[d1_index];
				// fill T0D2 information
				int d2_index = slice->index[1];
				if (d2_index >= merge_event[1].hit) {
					std::cerr << "Error: D2 index over hit in run "
						<< info.run << ", entry "
						<< info.entry << ", index " << d2_index
						<< ", hit " << merge_event[1].hit << "\n";
					return -1;
				}
				t0_event.energy[i][1] = merge_event[1].energy[d2_index];
				t0_event.time[i][1] = merge_event[1].time[d2_index];
				t0_event.x[i][1] = merge_event[1].x[d2_index];
				t0_event.y[i][1] = merge_event[1].y[d2_index];
				t0_event.z[i][1] = 111.76;
				t0_event.dssd_flag[i][1] = merge_event[1].flag[d2_index];
				t0_event.hole[i] = selected_d2_hole->at(d2_index);
				if (detect_layer == 0 && !slice->straight) {
					t0_event.straight[i] = false;
				} 
			} else if (layer == 1) {
				// fill T0D3 information
				int d3_index = slice->index[1];
				if (d3_index >= merge_event[2].hit) {
					std::cerr << "Error: D1 index over hit in run "
						<< info.run << ", entry "
						<< info.entry << ", index " << d3_index
						<< ", hit " << merge_event[2].hit << " .\n";
					return -1;
				}
				t0_event.energy[i][2] = merge_event[2].energy[d3_index];
				t0_event.time[i][2] = merge_event[2].time[d3_index];
				t0_event.x[i][2] = merge_event[2].x[d3_index];
				t0_event.y[i][2] = merge_event[2].y[d3_index];
				t0_event.z[i][2] = 123.52;
				t0_event.dssd_flag[i][2] = merge_event[2].flag[d3_index];
				if (detect_layer == 1 && !slice->straight) {
					t0_event.straight[i] = false;
				}
			} else if (layer == 2) {
				// fill T0S1 information
				t0_event.ssd_energy[0] = s1_event.energy;
				t0_event.ssd_time[0] = s1_event.time;
			} else if (layer == 3) {
				t0_event.ssd_energy[1] = s2_event.energy;
				t0_event.ssd_time[1] = s2_event.time;
			} else if (layer == 4) {
				t0_event.ssd_energy[2] = s3_event.energy;
				t0_event.ssd_time[2] = s3_event.time;
			}
			node = &(nodes[node->parent]);
		}
	}

	return 0;
}



int main(int argc, char **argv) {
	if (argc < 2) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	if (pos_start >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Miss run argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[pos_start]);

	// T0D1 normalize result event file name
	TString d1_file_name = TString::Format(
		"%s%st0d1-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// T0D1 normalize result file
	TFile d1_file(d1_file_name, "read");
	// tree
	TTree *ipt = (TTree*)d1_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< d1_file_name << " failed.\n";
		return -1;
	}
	// T0D1 noramlzie result event
	DssdFundamentalEvent d1;
	// setup input branches
	d1.SetupInput(ipt);

	// add T0D2 data
	// T0D2 file name
	TString d2_file_name = TString::Format(
		"%s%st0d2-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("d2=tree", d2_file_name);
	// T0D2 normalize result event
	DssdFundamentalEvent d2;
	// setup input branches
	d2.SetupInput(ipt, "d2.");

	// add T0D3 data
	// T0D3 file name
	TString d3_file_name = TString::Format(
		"%s%st0d3-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("d3=tree", d3_file_name);
	// T0D3 normalize result event
	DssdFundamentalEvent d3;
	// setup input branches
	d3.SetupInput(ipt, "d3.");

	// add T0S1 data
	// T0S1 file name
	TString s1_file_name = TString::Format(
		"%s%st0s1-fundamental-%s%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("s1=tree", s1_file_name);
	// T0S1 event
	SsdEvent s1;
	// setup input branches
	s1.SetupInput(ipt, "s1.");

	// add T0S2 data
	// T0S2 file name
	TString s2_file_name = TString::Format(
		"%s%st0s2-fundamental-%s%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("s2=tree", s2_file_name);
	// T0S2 event
	SsdEvent s2;
	// setup input branches
	s2.SetupInput(ipt, "s2.");

	// add T0S3 data
	// T0S3 file name
	TString s3_file_name = TString::Format(
		"%s%st0s3-fundamental-%s%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("s3=tree", s3_file_name);
	// T0S3 event
	SsdEvent s3;
	// setup input branches
	s3.SetupInput(ipt, "s3.");

	HoleInfo hole_info;
	if (run == 2) {
		hole_info.simulate = true;
		hole_info.generator = new TRandom3(201);

		// T0D2 average hole possibility file name
		TString hole_possibility_file_name = TString::Format(
			"%s%saverage-hole.root", kGenerateDataPath, kHoleDir
		);
		// T0D2 average hole file
		TFile hole_possibility_file(hole_possibility_file_name, "read");
		// average hole histogram
		TH2F *hahp = (TH2F*)hole_possibility_file.Get("hahp");
		if (!hahp) {
			std::cerr << "Error: Get hole average from "
				<< hole_possibility_file_name << " failed.\n";
			return -1;
		}
		for (int fs = 0; fs < 32; ++fs) {
			for (int bs = 0; bs < 32; ++bs) {
				hole_info.hole_possibility[fs][bs] =
					hahp->GetBinContent(fs+1, bs+1);
			}
		}
		// close file
		hole_possibility_file.Close();
	} else {
		hole_info.simulate = false;

		// T0D2 pixel resolution file name
		TString hole_flag_file_name = TString::Format(
			"%s%shole-flag.root", kGenerateDataPath, kHoleDir
		);
		// pixel resolution file
		TFile hole_flag_file(hole_flag_file_name, "read");
		// resolution histogram
		TH2F *hhsr = (TH2F*)hole_flag_file.Get("hhsr");
		if (!hhsr) {
			std::cerr << "Error: Get hole start run from "
				<< hole_flag_file_name << " failed.\n";
			return -1;
		}
		for (int fs = 0; fs < 32; ++fs) {
			for (int bs = 0; bs < 32; ++bs) {
				hole_info.hole_flag[fs][bs] = false;
				int start_run = hhsr->GetBinContent(fs+1, bs+1);
				if (start_run == 0) continue;
				if (run >= start_run) hole_info.hole_flag[fs][bs] = true;
			}
		}
		hole_flag_file.Close();
	}


	// output file name
	TString output_file_name = TString::Format(
		"%s%st0-telescope-%sv2-%04d.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "merge and slice track T0 events");
	// output event
	T0Event t0;
	// setup output branches
	t0.SetupOutput(&opt);

	// output merge file
	TString merge_file_name[3];
	for (int i = 0; i < 3; ++i) {
		merge_file_name[i].Form(
			"%s%st0d%d-mt-merge-%s%04d.root",
			kGenerateDataPath,
			kMergeDir,
			i+1,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
	}
	// output merge file
	TFile *merge_file[3];
	// output merge tree
	TTree *merge_tree[3];
	// output merge event
	DssdMergeEvent merge_event[3];
	for (int i = 0; i < 3; ++i) {
		merge_file[i] = new TFile(merge_file_name[i], "recreate");
		merge_tree[i] = new TTree("tree", "merge-track middle result");
		merge_event[i].SetupOutput(merge_tree[i]);
	}

	// T0 cuts
	T0Cut t0_cut;

	// T0D1-D2 cuts
	t0_cut.d1d2_cuts.push_back({2, 3, ReadCut("t0-d1d2-p1i1-3He")});
	t0_cut.d1d2_cuts.push_back({2, 4, ReadCut("t0-d1d2-p1i1-4He")});
	t0_cut.d1d2_cuts.push_back({2, 6, ReadCut("t0-d1d2-p1i1-6He")});
	t0_cut.d1d2_cuts.push_back({3, 6, ReadCut("t0-d1d2-p1i1-6Li")});
	t0_cut.d1d2_cuts.push_back({3, 7, ReadCut("t0-d1d2-p1i1-7Li")});
	t0_cut.d1d2_cuts.push_back({3, 8, ReadCut("t0-d1d2-p1i1-8Li")});
	t0_cut.d1d2_cuts.push_back({3, 9, ReadCut("t0-d1d2-p1i1-9Li")});
	t0_cut.d1d2_cuts.push_back({4, 7, ReadCut("t0-d1d2-p1i1-7Be")});
	t0_cut.d1d2_cuts.push_back({4, 9, ReadCut("t0-d1d2-p1i1-9Be")});
	t0_cut.d1d2_cuts.push_back({4, 10, ReadCut("t0-d1d2-p1i1-10Be")});
	t0_cut.d1d2_cuts.push_back({5, 10, ReadCut("t0-d1d2-p1i1-10B")});
	t0_cut.d1d2_cuts.push_back({5, 11, ReadCut("t0-d1d2-p1i1-11B")});
	t0_cut.d1d2_cuts.push_back({5, 12, ReadCut("t0-d1d2-p1i1-12B")});
	t0_cut.d1d2_cuts.push_back({5, 13, ReadCut("t0-d1d2-p1i1-13B")});
	t0_cut.d1d2_cuts.push_back({6, 12, ReadCut("t0-d1d2-p1i1-12C")});
	t0_cut.d1d2_cuts.push_back({6, 13, ReadCut("t0-d1d2-p1i1-13C")});
	t0_cut.d1d2_cuts.push_back({6, 14, ReadCut("t0-d1d2-p1i1-14C")});
	// t0_cut.d1d2_cuts.push_back({6, 15, ReadCut("t0-d1d2-p1i1-15C")});
	// ChenYing's cuts for 16C run
	// t0_cut.d1d2_cuts.push_back({2, 4, ReadCut("t0-d1d2-cy-4He")});
	// t0_cut.d1d2_cuts.push_back({2, 6, ReadCut("t0-d1d2-cy-6He")});
	// t0_cut.d1d2_cuts.push_back({4, 10, ReadCut("t0-d1d2-cy-10Be")});
	// t0_cut.d1d2_cuts.push_back({4, 12, ReadCut("t0-d1d2-cy-12Be")});

	// T0D1-D2 tail cuts
	t0_cut.d1d2_tail_cuts.push_back({2, 0, ReadCut("t0-d1d2-tail-p1i1-He")});
	t0_cut.d1d2_tail_cuts.push_back({3, 0, ReadCut("t0-d1d2-tail-p1i1-Li")});
	t0_cut.d1d2_tail_cuts.push_back({4, 0, ReadCut("t0-d1d2-tail-p1i1-Be")});
	// ChenYing's cuts for 16C run
	// t0_cut.d1d2_tail_cuts.push_back({2, 0, ReadCut("t0-d1d2-tail-cy-He")});
	// t0_cut.d1d2_tail_cuts.push_back({4, 0, ReadCut("t0-d1d2-tail-cy-Be")});

	// T0D2-D3 cuts
	t0_cut.d2d3_cuts.push_back({2, 3, ReadCut("t0-d2d3-p1i1-3He")});
	t0_cut.d2d3_cuts.push_back({2, 4, ReadCut("t0-d2d3-p1i1-4He")});
	t0_cut.d2d3_cuts.push_back({2, 6, ReadCut("t0-d2d3-p1i1-6He")});
	t0_cut.d2d3_cuts.push_back({3, 6, ReadCut("t0-d2d3-p1i1-6Li")});
	t0_cut.d2d3_cuts.push_back({3, 7, ReadCut("t0-d2d3-p1i1-7Li")});
	t0_cut.d2d3_cuts.push_back({3, 8, ReadCut("t0-d2d3-p1i1-8Li")});
	t0_cut.d2d3_cuts.push_back({3, 9, ReadCut("t0-d2d3-p1i1-9Li")});
	t0_cut.d2d3_cuts.push_back({4, 9, ReadCut("t0-d2d3-p1i1-9Be")});
	t0_cut.d2d3_cuts.push_back({4, 10, ReadCut("t0-d2d3-p1i1-10Be")});
	// ChenYing's cuts for 16C run
	// t0_cut.d2d3_cuts.push_back({2, 4, ReadCut("t0-d2d3-cy-4He")});
	// t0_cut.d2d3_cuts.push_back({2, 6, ReadCut("t0-d2d3-cy-6He")});
	// t0_cut.d2d3_cuts.push_back({4, 10, ReadCut("t0-d2d3-cy-10Be")});
	// t0_cut.d2d3_cuts.push_back({4, 12, ReadCut("t0-d2d3-cy-12Be")});

	// T0D2-D3 tail cuts
	t0_cut.d2d3_tail_cuts.push_back({2, 0, ReadCut("t0-d2d3-tail-p1i1-He")});
	t0_cut.d2d3_tail_cuts.push_back({3, 0, ReadCut("t0-d2d3-tail-p1i1-Li")});
	// ChenYing's cuts for 16C run
	// t0_cut.d2d3_tail_cuts.push_back({2, 0, ReadCut("t0-d2d3-tail-cy-He")});


	// T0D3-S1 cuts
	t0_cut.d3s1_cuts.push_back({2, 3, ReadCut("t0-d3s1-p1i1-3He")});
	t0_cut.d3s1_cuts.push_back({2, 4, ReadCut("t0-d3s1-p1i1-4He")});
	t0_cut.d3s1_cuts.push_back({2, 6, ReadCut("t0-d3s1-p1i1-6He")});
	t0_cut.d3s1_cuts.push_back({3, 6, ReadCut("t0-d3s1-p1i1-6Li")});
	t0_cut.d3s1_cuts.push_back({3, 7, ReadCut("t0-d3s1-p1i1-7Li")});
	t0_cut.d3s1_cuts.push_back({3, 8, ReadCut("t0-d3s1-p1i1-8Li")});
	t0_cut.d3s1_cuts.push_back({3, 9, ReadCut("t0-d3s1-p1i1-9Li")});
	// ChenYing's cuts for 16C run
	// t0_cut.d3s1_cuts.push_back({2, 4, ReadCut("t0-d3s1-cy-4He")});
	// t0_cut.d3s1_cuts.push_back({2, 6, ReadCut("t0-d3s1-cy-6He")});

	// T0D3-S1 tail cuts
	t0_cut.d3s1_tail_cuts.push_back({2, 4, ReadCut("t0-d3s1-tail-p1i1-He")});
	// ChenYing's cuts for 16C run
	// t0_cut.d3s1_tail_cuts.push_back({2, 0, ReadCut("t0-d3s1-tail-cy-He")});

	// T0S1-S2 cuts
	t0_cut.s1s2_cuts.push_back({2, 3, ReadCut("t0-s1s2-p1i1-3He")});
	t0_cut.s1s2_cuts.push_back({2, 4, ReadCut("t0-s1s2-p1i1-4He")});
	t0_cut.s1s2_cuts.push_back({2, 6, ReadCut("t0-s1s2-p1i1-6He")});
	// ChenYing's cuts for 16C run
	// t0_cut.s1s2_cuts.push_back({2, 4, ReadCut("t0-s1s2-cy-4He")});
	// t0_cut.s1s2_cuts.push_back({2, 6, ReadCut("t0-s1s2-cy-6He")});

	// T0S1-S2 tail cuts
	t0_cut.s1s2_tail_cuts.push_back({2, 4, ReadCut("t0-s1s2-tail-p1i1-He")});
	t0_cut.s1s2_tail_cuts.push_back({2, 6, ReadCut("t0-s1s2-tail-p1i1-6He")});
	// ChenYing's cuts for 16C run
	// t0_cut.s1s2_tail_cuts.push_back({2, 0, ReadCut("t0-s1s2-tail-cy-He")});

	// T0S2-S3 cuts
	t0_cut.s2s3_cuts.push_back({2, 3, ReadCut("t0-s2s3-p1i1-3He")});
	t0_cut.s2s3_cuts.push_back({2, 4, ReadCut("t0-s2s3-p1i1-4He")});
	t0_cut.s2s3_cuts.push_back({2, 6, ReadCut("t0-s2s3-p1i1-6He")});
	// ChenYing's cuts for 16C run
	// t0_cut.s2s3_cuts.push_back({2, 4, ReadCut("t0-s2s3-cy-4He")});
	// t0_cut.s2s3_cuts.push_back({2, 6, ReadCut("t0-s2s3-cy-6He")});

	// T0S2-S3 pass cuts
	t0_cut.s2s3_tail_cuts.push_back({2, 4, ReadCut("t0-s2s3-tail-p1i1-4He")});
	t0_cut.s2s3_tail_cuts.push_back({2, 6, ReadCut("t0-s2s3-tail-p1i1-6He")});
	// ChenYing's cuts for 16C run
	// t0_cut.s2s3_tail_cuts.push_back({2, 4, ReadCut("t0-s2s3-tail-cy-4He")});
	// t0_cut.s2s3_tail_cuts.push_back({2, 6, ReadCut("t0-s2s3-tail-cy-6He")});

	// T0D1-D2 hole cuts
	t0_cut.d1d2_hole_cuts.push_back(
		{2, 4, ReadCut("t0-d1d2-hole-p1i1-4He")}
	);
	t0_cut.d1d2_hole_cuts.push_back(
		{4, 10, ReadCut("t0-d1d2-hole-p1i1-10Be")}
	);
	// T0D1-D2 hole tail cuts
	t0_cut.d1d2_hole_tail_cuts.push_back(
		{2, 0, ReadCut("t0-d1d2-hole-tail-p1i1-He")}
	);
	t0_cut.d1d2_hole_tail_cuts.push_back(
		{4, 0, ReadCut("t0-d1d2-hole-tail-p1i1-Be")}
	);

	// T0D2-D3 hole cuts
	t0_cut.d2d3_hole_cuts.push_back(
		{2, 4, ReadCut("t0-d2d3-hole-p1i1-4He")});
	t0_cut.d2d3_hole_cuts.push_back(
		{4, 10, ReadCut("t0-d2d3-hole-p1i1-10Be")}
	);
	// T0D2-D3 hole tail cuts
	t0_cut.d2d3_hole_tail_cuts.push_back(
		{2, 0, ReadCut("t0-d2d3-hole-tail-p1i1-He")}
	);


	// statistics
	int found_behe_counts[4];
	for (int i = 0; i < 4; ++i) found_behe_counts[i] = 0;

	// // time statistics
	// std::chrono::duration<double, std::milli> get_data_time =
	// 	std::chrono::duration<double, std::milli>::zero();
	// std::chrono::duration<double, std::milli> d1_merge_time = get_data_time;
	// std::chrono::duration<double, std::milli> d2_merge_time = get_data_time;
	// std::chrono::duration<double, std::milli> d3_merge_time = get_data_time;
	// std::chrono::duration<double, std::milli> track_time = get_data_time;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Merge and slice track   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// const auto start_time_point =
		// 	std::chrono::high_resolution_clock::now();

		// get data
		ipt->GetEntry(entry);

		// const auto get_data_time_point =
		// 	std::chrono::high_resolution_clock::now();
		// get_data_time += get_data_time_point - start_time_point;

		if (print_debug) {
			std::cout << "========== entry " << entry << " ==========\n";
		}


		std::vector<DssdMergeEvent> d1_merge;
		std::vector<std::vector<bool>> d1_hole;
		GetMergeEvents(
			d1, d1_range, T0D1PositionFromStrip, nullptr,
			d1_merge, d1_hole
		);
		if (print_debug) {
			std::cout << "---------- D1 events ----------\n"
				<< "Possible cases " << d1_merge.size() << "\n";
			for (size_t i = 0; i < d1_merge.size(); ++i) {
				std::cout << "Case " << i << ", hit " << d1_merge[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d1_merge[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d1_merge[i].flag[j] << ", "
						<< std::dec << d1_merge[i].merge_tag[j] << ", "
						<< d1_merge[i].x[j] << ", "
						<< d1_merge[i].y[j] << ", "
						<< d1_merge[i].energy[j] << "\n";
				}
			}
		}

		// const auto d1_merge_time_point =
		// 	std::chrono::high_resolution_clock::now();
		// d1_merge_time += d1_merge_time_point - get_data_time_point;

		std::vector<DssdMergeEvent> d2_merge;
		std::vector<std::vector<bool>> d2_hole;
		GetMergeEvents(
			d2, d2_range, T0D2PositionFromStrip, &hole_info,
			d2_merge, d2_hole
		);
		if (print_debug) {
			std::cout << "---------- D2 events ----------\n"
				<< "Possible cases " << d2_merge.size() << "\n";
			for (size_t i = 0; i < d2_merge.size(); ++i) {
				std::cout << "Case " << i << ", hit " << d2_merge[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d2_merge[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d2_merge[i].flag[j] << ", "
						<< std::dec << d2_merge[i].merge_tag[j] << ", "
						<< d2_merge[i].x[j] << ", "
						<< d2_merge[i].y[j] << ", "
						<< d2_merge[i].energy[j] << "\n";
				}
			}
		}

		// const auto d2_merge_time_point =
		// 	std::chrono::high_resolution_clock::now();
		// d2_merge_time += d2_merge_time_point - d1_merge_time_point;

		std::vector<DssdMergeEvent> d3_merge;
		std::vector<std::vector<bool>> d3_hole;
		GetMergeEvents(
			d3, d3_range, T0D3PositionFromStrip, nullptr,
			d3_merge, d3_hole
		);
		if (print_debug) {
			std::cout << "---------- D3 events ----------\n"
				<< "Possible cases " << d3_merge.size() << "\n";
			for (size_t i = 0; i < d3_merge.size(); ++i) {
				std::cout << "Case " << i << ", hit " << d3_merge[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d3_merge[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d3_merge[i].flag[j] << ", "
						<< std::dec << d3_merge[i].merge_tag[j] << ", "
						<< d3_merge[i].x[j] << ", "
						<< d3_merge[i].y[j] << ", "
						<< d3_merge[i].energy[j] << "\n";
				}
			}
		}

		// const auto d3_merge_time_point =
		// 	std::chrono::high_resolution_clock::now();
		// d3_merge_time += d3_merge_time_point - d2_merge_time_point;

		std::vector<DssdMergeEvent> d1_bind_merge1;
		std::vector<DssdMergeEvent> d1_bind_merge2;
		std::vector<DssdMergeEvent> d2_bind_merge1;
		std::vector<DssdMergeEvent> d2_bind_merge2;
		GetBindEvents(
			d1, d2, d1_range, d2_range, &hole_info,
			d1_bind_merge1, d1_bind_merge2,
			d2_bind_merge1, d2_bind_merge2
		);
		if (print_debug) {
			std::cout << "---------- D1 bind events 1 ----------\n"
				<< "Possible cases " << d1_bind_merge1.size() << "\n";
			for (size_t i = 0; i < d1_bind_merge1.size(); ++i) {
				std::cout << "Case " << i
					<< ", hit " << d1_bind_merge1[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d1_bind_merge1[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d1_bind_merge1[i].flag[j] << ", "
						<< std::dec << d1_bind_merge1[i].merge_tag[j] << ", "
						<< d1_bind_merge1[i].x[j] << ", "
						<< d1_bind_merge1[i].y[j] << ", "
						<< d1_bind_merge1[i].energy[j] << "\n";
				}
			}
			std::cout << "---------- D1 bind events 2 ----------\n"
				<< "Possible cases " << d1_bind_merge2.size() << "\n";
			for (size_t i = 0; i < d1_bind_merge2.size(); ++i) {
				std::cout << "Case " << i
					<< ", hit " << d1_bind_merge2[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d1_bind_merge2[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d1_bind_merge2[i].flag[j] << ", "
						<< std::dec << d1_bind_merge2[i].merge_tag[j] << ", "
						<< d1_bind_merge2[i].x[j] << ", "
						<< d1_bind_merge2[i].y[j] << ", "
						<< d1_bind_merge2[i].energy[j] << "\n";
				}
			}
			std::cout << "---------- D2 bind events 1 ----------\n"
				<< "Possible cases " << d2_bind_merge1.size() << "\n";
			for (size_t i = 0; i < d2_bind_merge1.size(); ++i) {
				std::cout << "Case " << i
					<< ", hit " << d2_bind_merge1[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d2_bind_merge1[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d2_bind_merge1[i].flag[j] << ", "
						<< std::dec << d2_bind_merge1[i].merge_tag[j] << ", "
						<< d2_bind_merge1[i].x[j] << ", "
						<< d2_bind_merge1[i].y[j] << ", "
						<< d2_bind_merge1[i].energy[j] << "\n";
				}
			}
			std::cout << "---------- D2 bind events 2 ----------\n"
				<< "Possible cases " << d2_bind_merge2.size() << "\n";
			for (size_t i = 0; i < d2_bind_merge2.size(); ++i) {
				std::cout << "Case " << i
					<< ", hit " << d2_bind_merge2[i].hit << "\n"
					<< "index, flag, tag, x, y, energy\n";
				for (int j = 0; j < d2_bind_merge2[i].hit; ++j) {
					std::cout << std::dec << j << ", "
						<< std::hex << d2_bind_merge2[i].flag[j] << ", "
						<< std::dec << d2_bind_merge2[i].merge_tag[j] << ", "
						<< d2_bind_merge2[i].x[j] << ", "
						<< d2_bind_merge2[i].y[j] << ", "
						<< d2_bind_merge2[i].energy[j] << "\n";
				}
			}
		}

		// for (auto &m : d1_merge) {
		// 	for (int i = 0; i < m.hit; ++i) {
		// 		m.energy[i] = cy_param[0][0] + cy_param[0][1] * m.energy[i];
		// 	}
		// }
		// for (auto &m : d2_merge) {
		// 	for (int i = 0; i < m.hit; ++i) {
		// 		m.energy[i] = cy_param[1][0] + cy_param[1][1] * m.energy[i];
		// 	}
		// }
		// for (auto &m : d3_merge) {
		// 	for (int i = 0; i < m.hit; ++i) {
		// 		m.energy[i] = cy_param[2][0] + cy_param[2][1] * m.energy[i];
		// 	}
		// }
		// s1.energy = cy_param[3][0] + cy_param[3][1] * s1.energy;
		// s2.energy = cy_param[4][0] + cy_param[4][1] * s2.energy;
		// s3.energy = cy_param[5][0] + cy_param[5][1] * s3.energy;

		// extra information, inclueds run and entry
		ExtraInfo extract_info;
		extract_info.run = run;
		extract_info.entry = entry;
		// times we found Be+He in this event
		int found_behe_times;
		if (SliceTrack(
			d1_merge, d2_merge, d3_merge,
			s1, s2, s3,
			d2_hole,
			d1_bind_merge1, d2_bind_merge1,
			d1_bind_merge2, d2_bind_merge2,
			t0_cut, extract_info,
			t0, merge_event, found_behe_times
		)) {
			std::cerr << "Error: Slice track in run " << run
				<< ", entry " << entry << " failed.\n";
			return -1;
		}

		if (found_behe_times > 0 && found_behe_times <= 3) {
			++found_behe_counts[found_behe_times-1];
		} else if (found_behe_times > 3) {
			++found_behe_counts[3];
		}
		// if (found_behe_times > 0) std::cout << entry << "\n";

		// const auto track_time_point =
		// 	std::chrono::high_resolution_clock::now();
		// track_time += track_time_point - d3_merge_time_point;

		opf.cd();
		opt.Fill();
		for (int i = 0; i < 3; ++i) {
			merge_file[i]->cd();
			merge_tree[i]->Fill();
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	std::vector<bool> identify;
	identify.push_back(true);

	// save tree
	opf.cd();
	opt.Write();
	opf.WriteObject(&identify, "identify");
	opf.Close();

	for (int i = 0; i < 3; ++i) {
		merge_file[i]->cd();
		merge_tree[i]->Write();
		merge_file[i]->Close();
	}

	d1_file.Close();


	// print statistics
	std::cout << "Found Be+He times:\n";
	for (int i = 0; i < 3; ++i) {
		std::cout << i+1 << ": " << found_behe_counts[i] << "\n";
	}
	std::cout << ">3: " << found_behe_counts[3] << "\n";

	// // print time statistics
	// std::cout << "Get Data " << get_data_time.count() << " ms\n"
	// 	<< "D1 merge " << d1_merge_time.count() << " ms\n"
	// 	<< "D2 merge " << d2_merge_time.count() << " ms\n"
	// 	<< "D3 merge " << d3_merge_time.count() << " ms\n"
	// 	<< "Track " << track_time.count() << " ms\n";

	return 0;
}