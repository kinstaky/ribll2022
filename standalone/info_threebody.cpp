#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/channel_event.h"
#include "include/event/dssd_event.h"
#include "include/event/particle_event.h"
#include "include/event/particle_type_event.h"
#include "include/event/ta_event.h"
#include "include/event/t0_event.h"
#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;

int FillResult(
	unsigned short flag,
	int result_hit,
	double *result_channel,
	double *result_time,
	unsigned short *result_strip,
	int &hit,
	double *channel,
	double *time,
	unsigned int *strip
) {
	hit = 0;
	for (int i = 0; i < 8; ++i) {
		if ((flag & (1 << i)) == 0) continue;
		if (i >= result_hit) {
			std::cerr << "Error: Found flag bit over result hit.\n";
			return -1;
		}
		channel[hit] = result_channel[i];
		time[hit] = result_time[i];
		strip[hit] = result_strip[i];
		++hit;
	}
	return 0;
}


/// @brief rebuild threebody reaction process
/// @param[inout] event input event
/// @param[in] vppac use VPPAC to calculate
/// @returns Q value
///
double ThreeBodyProcess(ThreeBodyInfoEvent &event, bool vppac = false) {
	double tx = vppac ? event.vptx : event.xptx;
	double ty = vppac ? event.vpty : event.xpty;

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - tx,
		event.be_y[0] - ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - tx,
		event.he_y[0] - ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - tx,
		event.d_y - ty,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinematic energy
	event.c14_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = event.t0_energy[0] + event.t0_energy[1]
		+ event.taf_energy - event.c14_kinetic;

	return q;
}


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options]\n"
		"  run               Set run number.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -r recoil_mass    Set recoil particle mass number.\n"
		"  -s                Use simulated data.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] recoil_mass recoil particle mass number
/// @param[out] sim use simulated data
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	int &recoil_mass,
	bool &sim
) {
	// initialize
	help = false;
	recoil_mass = 2;
	sim = false;
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
		} else if (argv[result][1] == 'r') {
			// option of recoil mass number
			// get mass number in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			recoil_mass = atoi(argv[result]);
		} else if (argv[result][1] == 's') {
			sim = true;
		} else {
			return -result;
		}
	}
	return result;
}

int main(int argc, char **argv) {
	// help flag
	bool help = false;
	// recoil particle mass number
	int recoil_mass = 2;
	// use simulated data
	bool sim = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, recoil_mass, sim);

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

	if (recoil_mass <= 0 || recoil_mass > 3) {
		std::cerr << "Error: Invalid recoil particle mass "
			<< recoil_mass << "\n";
		return -1;
	}

	std::string tag = sim ? "sim-ta" : "ta";

	// output file name
	TString output_file_name = TString::Format(
		"%sinfo/threebody%s%s.root",
		kGenerateDataPath,
		argc == 1 ? "" : TString::Format("-%dH", recoil_mass).Data(),
		sim ? "-sim" : ""
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "information of three body");
	// output data
	ThreeBodyInfoEvent event;
	int vppac_num, xppac_num;
	// setup output branches
	event.SetupOutput(&opt);
	opt.Branch("vnum", &vppac_num, "vnum/I");
	opt.Branch("xnum", &xppac_num, "xnum/I");

	// possible 2H stopped in ADSSD
	int possible_2H_num[4];
	for (int i = 0; i < 4; ++i) {
		possible_2H_num[i] = 0;
	}
	// total number of possible 2H events
	int total_possible_2H = 0;

	elc::CsiEnergyCalculator calculator("2H");

	// statistics for calculating efficiency
	size_t t0_valid_count = 0;
	size_t taf_valid_count = 0;
	size_t ppac_valid_count = 0;
	size_t valid_count = 0;
	size_t total_count = 0;


	for (
		unsigned int run = sim ? 0 : 618;
		run <= (sim ? 0 : 716);
		++run
	) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		std::cout << "Processing run " << run << "\n";
		// input file name
		TString channel_file_name = TString::Format(
			"%s%sC14-10Be-4He-%dH-%s%04d.root",
			kGenerateDataPath,
			kChannelDir,
			recoil_mass,
			sim ? "sim-" : "",
			run
		);
		// input file
		TFile channel_file(channel_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)channel_file.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< channel_file_name << " failed.\n";
			return -1;
		}
		// input channel event
		ChannelEvent channel;
		// t0 index
		int t0_index[8];
		// PPAC flag
		int ppac_flag;
		// setup input branches
		channel.SetupInput(ipt);
		ipt->SetBranchAddress("t0_index", t0_index);
		ipt->SetBranchAddress("ppac_flag", &ppac_flag);

		std::vector<long long> valid_entries;
		std::vector<int> taf_indexes;
		std::vector<int> be10_indexes;
		std::vector<int> he4_indexes;
		std::vector<int> ppac_flags;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			be10_indexes.push_back(t0_index[0]);
			he4_indexes.push_back(t0_index[1]);
			ppac_flags.push_back(ppac_flag);
		}
		// close file
		channel_file.Close();

		// open telescope file to get more information
		TString t0_tele_file_name = TString::Format(
			"%s%st0-telescope-%s-%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			tag.c_str(),
			run
		);
		// t0 file
		TFile t0_tele_file(t0_tele_file_name, "read");
		// t0 tree
		TTree *t0_tree = (TTree*)t0_tele_file.Get("tree");
		if (!t0_tree) {
			std::cerr << "Error: Get tree from "
				<< t0_tele_file_name << " failed.\n";
			return -1;
		}
		// identify flag
		std::vector<bool> *identify =
			(std::vector<bool>*)t0_tele_file.Get("identify");
		if (!identify || !(identify->at(0))) {
			// add T0 particle identify file
			t0_tree->AddFriend("pid=tree", TString::Format(
				"%s%st0-particle-type-%s-%04u.root",
				kGenerateDataPath,
				kParticleIdentifyDir,
				tag.c_str(),
				run
			));
		}
		// add T0D1 normalize result friend
		t0_tree->AddFriend("t0d1=tree", TString::Format(
			"%s%st0d1-result-%s-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag.c_str(),
			run
		));
		// add T0D2 normalize result friend
		t0_tree->AddFriend("t0d2=tree", TString::Format(
			"%s%st0d2-result-%s-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag.c_str(),
			run
		));
		// add T0D3 normalize result friend
		t0_tree->AddFriend("t0d3=tree", TString::Format(
			"%s%st0d3-result-%s-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag.c_str(),
			run
		));
		for (int i = 0; i < 6; ++i) {
			// add TAF telescope tree as friend
			t0_tree->AddFriend(
				TString::Format("taf%d=tree", i),
				TString::Format(
					"%s%staf%d-telescope-%s-%04u.root",
					kGenerateDataPath,
					kTelescopeDir,
					i,
					tag.c_str(),
					run
				)
			);
			// add TAF ADSSD fundamental tree as friend
			if (sim) continue;
			t0_tree->AddFriend(
				TString::Format("tafd%d=tree", i),
				TString::Format(
					"%s%stafd%d-fundamental-ta-%04u.root",
					kGenerateDataPath,
					kFundamentalDir,
					i,
					run
				)
			);
		}
		// add xppac as friend
		t0_tree->AddFriend(
			"xppac=tree",
			TString::Format(
				"%s%sxppac-particle-%s-%04u.root",
				kGenerateDataPath,
				kParticleDir,
				tag.c_str(),
				run
			)
		);
		// add vppac as friend
		t0_tree->AddFriend(
			"vppac=tree",
			TString::Format(
				"%s%svppac-particle-%s-%04u.root",
				kGenerateDataPath,
				kParticleDir,
				tag.c_str(),
				run
			)
		);
		// input events
		// T0 telescope event
		T0Event t0;
		ParticleTypeEvent t0_type;
		// T0Dx normalize result event
		DssdFundamentalEvent dssd_result[3];
		// TAF telescope events
		TaEvent taf[6];
		// TAF ADSSD fundamental events
		DssdFundamentalEvent tafd[6];
		// XPPAC events
		ParticleEvent xppac;
		// VPPAC events
		ParticleEvent vppac;

		// setup input branches
		t0.SetupInput(t0_tree);
		if (!identify || !(identify->at(0))) {
			t0_type.SetupInput(t0_tree, "pid.");
		} else {
			t0_tree->SetBranchAddress("mass", t0_type.mass);
			t0_tree->SetBranchAddress("charge", t0_type.charge);
			t0_tree->SetBranchAddress("layer", t0_type.layer);
		}
		dssd_result[0].SetupInput(t0_tree, "t0d1.");
		dssd_result[1].SetupInput(t0_tree, "t0d2.");
		dssd_result[2].SetupInput(t0_tree, "t0d3.");
		for (int i = 0; i < 6; ++i) {
			taf[i].SetupInput(t0_tree, TString::Format("taf%d.", i).Data());
			if (sim) continue;
			tafd[i].SetupInput(t0_tree, TString::Format("tafd%d.", i).Data());
		}
		xppac.SetupInput(t0_tree, "xppac.");
		t0_tree->SetBranchAddress("xppac.xflag", &event.xppac_xflag);
		t0_tree->SetBranchAddress("xppac.yflag", &event.xppac_yflag);
		vppac.SetupInput(t0_tree, "vppac.");
		t0_tree->SetBranchAddress("vppac.xflag", &event.vppac_xflag);
		t0_tree->SetBranchAddress("vppac.yflag", &event.vppac_yflag);


		// total number of entries in this run
		size_t entries = valid_entries.size();
		total_count += entries;
		// 1/100 of entries, for showing process
		size_t entry100 = entries / 100 + 1;
		// show start with simulated data
		if (sim) {
			printf("Collecting information   0%%");
			fflush(stdout);
		}
		for (size_t i = 0; i < valid_entries.size(); ++i) {
			// show process
			if (sim && i % entry100 == 0) {
				printf("\b\b\b\b%3ld%%", i / entry100);
				fflush(stdout);
			}

			t0_tree->GetEntry(valid_entries[i]);

			if (event.xppac_xflag != 0 && event.xppac_yflag != 0) {
				++ppac_valid_count;
			}
			if (taf_indexes[i] == -1) {
				event.taf_flag = 3;
				opt.Fill();
				continue;
			} else if (taf_indexes[i] == -2) {
				event.taf_flag = 4;
				++taf_valid_count;
				opt.Fill();
				continue;
			} else if (taf_indexes[i] == -4) {
				event.taf_flag = 4;
				++t0_valid_count;
				opt.Fill();
				continue;
			} else if (taf_indexes[i] == -6) {
				event.taf_flag = 4;
				opt.Fill();
				continue;
			}
			++t0_valid_count;
			++taf_valid_count;
			if (event.xppac_xflag != 0 && event.xppac_yflag != 0) {
				++valid_count;
			}

			// T0
			// T0 layer
			event.layer[0] = t0_type.layer[be10_indexes[i]];
			event.layer[1] = t0_type.layer[he4_indexes[i]];
			// T0 channels
			for (int j = 0; j < 3; ++j) {
				event.be_channel[j] = t0.energy[be10_indexes[i]][j];
				event.he_channel[j] = t0.energy[he4_indexes[i]][j];
			}
			for (int j = 0; j < 3; ++j) {
				event.ssd_channel[j] = t0.ssd_energy[j];
			}

			// calculate energy
			// 10Be kinetic energy
			// T0D1
			event.t0_energy[0] =
				t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
			// T0D2
			event.t0_energy[0] +=
				t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
			// T0D3
			if (event.layer[0] > 1) {
				event.t0_energy[0] +=
					t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
			}
			// T0S1
			if (event.layer[0] > 2) {
				event.t0_energy[0] +=
					t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
			}
			if (event.layer[0] > 3) {
				event.t0_energy[0] +=
					t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
			}
			if (event.layer[0] > 4) {
				event.t0_energy[0] +=
					t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
			}

			// 4He kinetic energy
			// T0D1
			event.t0_energy[1] =
				t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
			// T0D2
			event.t0_energy[1] +=
				t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
			// T0D3
			if (event.layer[1] > 1) {
				event.t0_energy[1] +=
					t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
			}
			if (event.layer[1] > 2) {
				event.t0_energy[1] +=
					t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
			}
			if (event.layer[1] > 3) {
				event.t0_energy[1] +=
					t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
			}
			if (event.layer[1] > 4) {
				event.t0_energy[1] +=
					t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
			}

			// position
			for (int j = 0; j < 3; ++j) {
				event.be_x[j] = t0.x[be10_indexes[i]][j];
				event.be_y[j] = t0.y[be10_indexes[i]][j];
				event.he_x[j] = t0.x[he4_indexes[i]][j];
				event.he_y[j] = t0.y[he4_indexes[i]][j];
			}

			// fill T0D1 10Be result
			if (FillResult(
				t0.dssd_flag[be10_indexes[i]][0] >> 8,
				dssd_result[0].back_hit,
				dssd_result[0].back_energy,
				dssd_result[0].back_time,
				dssd_result[0].back_strip,
				event.be_x_hit[0],
				event.be_x_channel[0],
				event.be_x_time[0],
				event.be_x_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 10Be T0D1 back. \n";
				return -1;
			}
			if (FillResult(
				t0.dssd_flag[be10_indexes[i]][0],
				dssd_result[0].front_hit,
				dssd_result[0].front_energy,
				dssd_result[0].front_time,
				dssd_result[0].front_strip,
				event.be_y_hit[0],
				event.be_y_channel[0],
				event.be_y_time[0],
				event.be_y_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 10Be T0D1 front. \n";
				return -1;
			}

			// fill T0D1 4He result
			if (FillResult(
				t0.dssd_flag[he4_indexes[i]][0] >> 8,
				dssd_result[0].back_hit,
				dssd_result[0].back_energy,
				dssd_result[0].back_time,
				dssd_result[0].back_strip,
				event.he_x_hit[0],
				event.he_x_channel[0],
				event.he_x_time[0],
				event.he_x_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 4He T0D1 back. \n";
				return -1;
			}
			if(FillResult(
				t0.dssd_flag[he4_indexes[i]][0],
				dssd_result[0].front_hit,
				dssd_result[0].front_energy,
				dssd_result[0].front_time,
				dssd_result[0].front_strip,
				event.he_y_hit[0],
				event.he_y_channel[0],
				event.he_y_time[0],
				event.he_y_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 4He T0D1 front. \n";
				return -1;
			}

			// fill T0D2 and T0D3 10Be result
			for (int j = 1; j < event.layer[0]+1; ++j) {
				if (j >= 3) break;
				if (FillResult(
					t0.dssd_flag[be10_indexes[i]][j],
					dssd_result[j].front_hit,
					dssd_result[j].front_energy,
					dssd_result[j].front_time,
					dssd_result[j].front_strip,
					event.be_x_hit[j],
					event.be_x_channel[j],
					event.be_x_time[j],
					event.be_x_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 10Be T0D" << j+1 << " front. \n";
					return -1;
				}
				if (FillResult(
					t0.dssd_flag[be10_indexes[i]][j] >> 8,
					dssd_result[j].back_hit,
					dssd_result[j].back_energy,
					dssd_result[j].back_time,
					dssd_result[j].back_strip,
					event.be_y_hit[j],
					event.be_y_channel[j],
					event.be_y_time[j],
					event.be_y_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 10Be T0D" << j+1 << " back. \n";
					return -1;
				}
			}

			// fill T0D2 and T0D3 4He result
			for (int j = 1; j < event.layer[1]+1; ++j) {
				if (j >= 3) break;
				// 4He result
				if (FillResult(
					t0.dssd_flag[he4_indexes[i]][j],
					dssd_result[j].front_hit,
					dssd_result[j].front_energy,
					dssd_result[j].front_time,
					dssd_result[j].front_strip,
					event.he_x_hit[j],
					event.he_x_channel[j],
					event.he_x_time[j],
					event.he_x_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 4He T0D" << j+1 << " front. \n";
					return -1;
				}
				if (FillResult(
					t0.dssd_flag[he4_indexes[i]][j] >> 8,
					dssd_result[j].back_hit,
					dssd_result[j].back_energy,
					dssd_result[j].back_time,
					dssd_result[j].back_strip,
					event.he_y_hit[j],
					event.he_y_channel[j],
					event.he_y_time[j],
					event.he_y_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 4He T0D" << j+1 << " back. \n";
					return -1;
				}
			}


			// TAF
			// TAF possible 2H number
			int possible_2H = 0;
			if (taf_indexes[i] == -1) {
				event.taf_flag = 1;
				event.csi_index = -1;
				// loop to search for possible 2H
				for (int j = 0; j < 6; ++j) {
					if (taf[j].num == 1 && taf[j].flag[0] == 0x1) {
						++possible_2H;
						taf_indexes[i] = j;
					}
				}
				if (possible_2H > 0 && possible_2H < 4) {
					++possible_2H_num[possible_2H-1];
				} else if (possible_2H >= 4) {
					++possible_2H_num[3];
				}
				++total_possible_2H;
				if (possible_2H != 1) continue;
				// TAFD energy
				event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
				// try to calculate CsI energy
				event.csi_energy = calculator.Energy(
					taf[taf_indexes[i]].theta[0], event.tafd_energy, 150.0
				);
				if (event.csi_energy < 0) continue;
				if (event.csi_energy < 6.0) {
					event.taf_energy = event.tafd_energy + event.csi_energy;
					event.taf_flag = 2;
				} else {
					event.taf_energy = taf[taf_indexes[i]].energy[0][0];
				}
			} else {
				// TAF flag
				event.taf_flag = 0;
				// CsI index
				if (taf[taf_indexes[i]].flag[0] == 0x3) {
					event.csi_index = taf_indexes[i]*2;
				} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
					event.csi_index = taf_indexes[i]*2 + 1;
				} else {
					event.csi_index = -1;
				}
				// CsI channel
				event.csi_channel = taf[taf_indexes[i]].energy[0][1];
				// TAFD energy
				event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
				// TAF-CsI energy
				// calibrated parameters
				double a0 = csi_param[event.csi_index][0];
				double a1 = csi_param[event.csi_index][1];
				double a2 = csi_param[event.csi_index][2];
				// CsI energy
				event.csi_energy = pow(
					(event.csi_channel - a2) / a0,
					1.0 / a1
				);
				// recoil particle energy
				event.taf_energy = event.tafd_energy + event.csi_energy;
			}

			// position
			ROOT::Math::Polar3DVector recoil_position(
				taf[taf_indexes[i]].radius[0],
				taf[taf_indexes[i]].theta[0],
				taf[taf_indexes[i]].phi[0]
			);
			event.d_x = recoil_position.X();
			event.d_y = recoil_position.Y();

			// channel
			event.d_x_channel = tafd[taf_indexes[i]].front_energy[0];
			event.d_y_channel = tafd[taf_indexes[i]].back_energy[0];
			// time
			event.d_x_time = tafd[taf_indexes[i]].front_time[0];
			event.d_y_time = tafd[taf_indexes[i]].back_time[0];
			// strip
			event.d_x_strip = tafd[taf_indexes[i]].front_strip[0];
			event.d_y_strip = tafd[taf_indexes[i]].back_strip[0];

			// PPAC
			vppac_num = vppac.num;
			xppac_num = xppac.num;
			event.ppac_flag = 0;
			if (event.xppac_xflag != 0 && event.xppac_yflag != 0) {
				event.ppac_flag |= 1;
				for (int j = 0; j < 3; ++j) {
					event.xppac_x[j] = xppac.x[j];
					event.xppac_y[j] = xppac.y[j];
					if (!sim) {
						event.xppac_x[j] -= ppac_correct[0][j];
						event.xppac_y[j] -= ppac_correct[1][j];
					}
				}
				// PPAC x track
				event.xppac_track[0] = TrackPpac(
					event.xppac_xflag, ppac_xz, event.xppac_x,
					event.t0_energy[0], event.t0_energy[1], event.taf_energy,
					event.be_x[0], event.he_x[0], event.d_x, event.d_y,
					event.xptx
				);
				// PPAC y track
				event.xppac_track[1] = TrackPpac(
					event.xppac_yflag, ppac_yz, event.xppac_y,
					event.t0_energy[0], event.t0_energy[1], event.taf_energy,
					event.be_y[0], event.he_y[0], event.d_y, event.d_x,
					event.xpty
				);

			}
			if (event.vppac_xflag != 0 && event.vppac_yflag != 0) {
				event.ppac_flag |= 2;
				for (int j = 0; j < 3; ++j) {
					event.vppac_x[j] = vppac.x[j];
					event.vppac_y[j] = vppac.y[j];
					if (!sim) {
						event.vppac_x[j] -= ppac_correct[0][j];
						event.vppac_y[j] -= ppac_correct[1][j];
					}
				}
				// PPAC x track
				event.vppac_track[0] = TrackPpac(
					event.vppac_xflag, ppac_xz, event.vppac_x,
					event.t0_energy[0], event.t0_energy[1], event.taf_energy,
					event.be_x[0], event.he_x[0], event.d_x, event.d_y,
					event.vptx
				);
				// PPAC y track
				event.vppac_track[1] = TrackPpac(
					event.vppac_yflag, ppac_yz, event.vppac_y,
					event.t0_energy[0], event.t0_energy[1], event.taf_energy,
					event.be_y[0], event.he_y[0], event.d_y, event.d_x,
					event.vpty
				);
			}

			// Q
			event.q = ThreeBodyProcess(event);
			event.vq = ThreeBodyProcess(event, true);

			// hole
			event.hole[0] = t0.hole[be10_indexes[i]];
			event.hole[1] = t0.hole[he4_indexes[i]];

			// extra information
			// run number
			event.run = run;
			// entry
			event.entry = valid_entries[i];

			opt.Fill();
		}
		t0_tele_file.Close();
	}
	// show finish
	if (sim) printf("\b\b\b\b100%%\n");
	if (sim) {
		// show efficiency
		std::cout << "Total efficiency: " << valid_count << " / " << total_count
			<< "  " << double(valid_count) / double(total_count) << "\n"
			<< "T0 efficency: " << t0_valid_count << " / " << total_count
			<< "  " << double(t0_valid_count) / double(total_count) << "\n"
			<< "TAF efficiency " << taf_valid_count << " / " << total_count
			<< "  " << double(taf_valid_count) / double(total_count) << "\n"
			<< "XPPAC efficiency: " << ppac_valid_count << " / " << total_count
			<< "  " << double(ppac_valid_count) / double(total_count) << "\n";
	} else {
		// show possible 2H statistics
		std::cout << "Possible 2H events total " << total_possible_2H << "\n";
		for (int i = 0; i < 4; ++i) {
			std::cout << "With TAF " << i+1 << ": " << possible_2H_num[i] << "\n";
		}
	}


	opf.cd();
	// save trees
	opt.Write();
	// close files
	opf.Close();
	return 0;
}