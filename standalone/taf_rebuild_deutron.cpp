#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/ta_event.h"
#include "include/event/particle_event.h"

using namespace ribll;


constexpr double taf_straight_strip_pars[12][16][2] = {
	{
		{278, -0.00170},
		{279, -0.00169},
		{265, -0.00177},
		{264, -0.00172},
		{253, -0.00178},
		{248, -0.00184},
		{280, -0.00169},
		{240, -0.00181},
		{223, -0.00192},
		{232, -0.00189},
		{244, -0.00182},
		{199, -0.00204},
		{216, -0.00197},
		{205, -0.00227},
		{365, -0.00187},
		{319, -0.00163}
	},
	{
		{258, -0.00160},
		{253, -0.00175},
		{254, -0.00172},
		{252, -0.00164},
		{249, -0.00164},
		{234, -0.00180},
		{251, -0.00170},
		{248, -0.00170},
		{222, -0.00183},
		{238, -0.00187},
		{225, -0.00190},
		{197, -0.00208},
		{203, -0.00207},
		{259, -0.00208},
		{316, -0.00160},
		{308, -0.00161}
	},
	{
		{223, -0.00157},
		{251, -0.00164},
		{290, -0.00148},
		{260, -0.00162},
		{257, -0.00152},
		{251, -0.00175},
		{266, -0.00157},
		{258, -0.00162},
		{228, -0.00181},
		{229, -0.00180},
		{230, -0.00180},
		{212, -0.00185},
		{212, -0.00191},
		{292, -0.00159},
		{292, -0.00159},
		{292, -0.00159}
	},
	{
		{248, -0.00142},
		{279, -0.00153},
		{298, -0.00138},
		{307, -0.00137},
		{207, -0.00147},
		{265, -0.00156},
		{285, -0.00138},
		{293, -0.00142},
		{267, -0.00156},
		{267, -0.00156},
		{230, -0.00180},
		{247, -0.00177},
		{227, -0.00186},
		{193, -0.00229},
		{179, -0.00201},
		{333, -0.00141}
	},
	{
		{330, -0.00136},
		{328, -0.00132},
		{314, -0.00146},
		{331, -0.00134},
		{327, -0.00139},
		{310, -0.00147},
		{307, -0.00149},
		{284, -0.00155},
		{305, -0.00150},
		{296, -0.00151},
		{274, -0.00171},
		{278, -0.00174},
		{286, -0.00166},
		{264, -0.00175},
		{393, -0.00132},
		{393, -0.00132}
	},
	{
		{323, -0.00145},
		{341, -0.00135},
		{324, -0.00140},
		{311, -0.00143},
		{311, -0.00138},
		{259, -0.00152},
		{268, -0.00155},
		{286, -0.00139},
		{297, -0.00140},
		{246, -0.00149},
		{246, -0.00148},
		{256, -0.00158},
		{257, -0.00185},
		{257, -0.00185},
		{351, -0.00136},
		{351, -0.00136}
	},
	{
		{325, -0.00149},
		{330, -0.00138},
		{327, -0.00148},
		{320, -0.00149},
		{313, -0.00150},
		{298, -0.00160},
		{290, -0.00155},
		{279, -0.00147},
		{274, -0.00160},
		{276, -0.00163},
		{274, -0.00160},
		{249, -0.00171},
		{253, -0.00165},
		{242, -0.00180},
		{376, -0.00139},
		{375, -0.00139}
	},
	{
		{307, -0.00146},
		{291, -0.00144},
		{285, -0.00161},
		{295, -0.00154},
		{291, -0.00159},
		{282, -0.00154},
		{278, -0.00154},
		{276, -0.00157},
		{265, -0.00169},
		{241, -0.00177},
		{237, -0.00186},
		{227, -0.00203},
		{233, -0.00165},
		{357, -0.00147},
		{357, -0.00147}
	},
	{
		{363, -0.00119},
		{352, -0.00120},
		{357, -0.00132},
		{350, -0.00130},
		{352, -0.00130},
		{339, -0.00129},
		{331, -0.00142},
		{317, -0.00148},
		{272, -0.00181},
		{304, -0.00154},
		{293, -0.00154},
		{304, -0.00163},
		{280, -0.00170},
		{277, -0.00171},
		{427, -0.00123},
		{427, -0.00123},
	},
	{
		{310, -0.00129},
		{293, -0.00122},
		{327, -0.00128},
		{344, -0.00127},
		{349, -0.00127},
		{319, -0.00126},
		{311, -0.00137},
		{306, -0.00139},
		{300, -0.00143},
		{290, -0.00143},
		{283, -0.00147},
		{253, -0.00166},
		{260, -0.00174},
		{186, -0.00210},
		{391, -0.00122},
		{391, -0.00122}
	},
	{
		{439, -0.00122},
		{439, -0.00122},
		{439, -0.00122},
		{439, -0.00122},
		{372, -0.00119},
		{369, -0.00124},
		{369, -0.00129},
		{345, -0.00140},
		{313, -0.00143},
		{313, -0.00150},
		{309, -0.00149},
		{273, -0.00163},
		{264, -0.00170},
		{232, -0.00174},
		{439, -0.00122},
		{439, -0.00122}
	},
	{
		{435, -0.00129},
		{435, -0.00129},
		{435, -0.00129},
		{435, -0.00129},
		{373, -0.00119},
		{362, -0.00136},
		{338, -0.00151},
		{330, -0.00155},
		{322, -0.00145},
		{321, -0.00151},
		{316, -0.00152},
		{301, -0.00153},
		{297, -0.00155},
		{256, -0.00181},
		{435, -0.00129},
		{435, -0.00129}
	}
};

constexpr double taf_strip_fit_d[12][16][2] = {
	{
		{70.9563, 2.91151},
		{71.5869, 2.61046},
		{71.6433, 2.86328},
		{72.1522, 2.91233},
		{72.409, 2.85624},
		{73.4555, 3.03226},
		{75.7065, 2.97385},
		{74.7924, 2.72869},
		{74.83, 2.72539},
		{76.2501, 2.56183},
		{77.5893, 2.65883},
		{76.156, 2.71645},
		{77.682, 2.96011},
		{76.0293, 3.23551},
		{81.8738, 1},		// bad
		{63.7151, 9.11853},	// bad
	},
	{
		{70.1927, 2.61913},
		{69.9992, 2.51606},
		{71.008, 2.52115},
		{71.804, 2.59668},
		{72.3324, 2.88918},
		{72.3028, 2.57013},
		{73.9351, 2.79208},
		{73.8755, 2.96614},
		{73.5435, 2.62488},
		{74.642, 2.41157},
		{74.9252, 2.79383},
		{73.7568, 2.94455},
		{74.3591, 2.59287},
		{75.4891, 3.07303},
		{71.1859, 9.99773},	// bad
		{73.1789, 10},		// bad
	},
	{
		{70.0229, 3.26877},
		{72.2865, 2.90135},
		{75.3183, 3.09839},
		{75.0111, 3.22725},
		{70.1663, 4.52677},
		{75.6148, 3.08326},
		{77.9126, 3.40648},
		{78.8249, 3.0249},
		{78.5111, 2.88591},
		{72.1873, 4.80619},
		{80.4705, 3.07033},
		{80.1554, 3.00339},
		{80.2595, 3.42284},
		{83.0955, 3.42624},	// bad
		{73.8772, 9.9998},	// bad
		{86.8104, 10},		// bad
	},
	{
		{72.821, 3.11485},
		{74.6757, 2.08471},
		{77.9893, 3.18651},
		{79.5343, 2.84952},
		{72.797, 4.55495},
		{79.8537, 3.22941},
		{81.8852, 3.15576},
		{82.7574, 2.95464},
		{82.1238, 3.04673},
		{75.4212, 4.53953},
		{81.72, 2.95845},
		{82.9792, 2.85219},
		{82.0828, 3.01708},
		{79.5779, 3.05013},
		{89.9999, 1},		// bad
		{68.8591, 9.84076},	// bad
	},
	{
		{80.278, 3.02623},
		{81.4021, 2.81276},
		{81.3156, 2.78278},
		{82.8441, 3.31716},
		{83.5145, 3.28173},
		{83.9544, 3.23711},
		{85.0607, 2.84483},
		{85.0486, 3.06979},
		{87.0133, 2.99593},
		{87.5536, 3.01965},
		{86.6726, 3.00914},
		{87.7781, 3.15065},
		{86.9635, 3.86769},
		{85.8217, 3.20234},
		{73.5797, 10},		// bad
		{89.3968, 1.93405},	// bad
	},
	{
		{77.7485, 3.28248},
		{80.2737, 3.13674},
		{80.7508, 3.01884},
		{80.932, 3.42156},
		{82.0975, 3.72343},
		{81.6138, 3.51401},
		{83.3395, 3.24172},
		{86.0342, 3.11961},
		{87.3417, 3.2062},
		{86.0449, 3.04998},
		{86.5602, 3.04414},
		{87.511, 2.89273},
		{83.3365, 3.99953},
		{82.0115, 4.02462},
		{64.2344, 9.73732},	// bad
		{80.7314, 10},		// bad
	},
	{
		{79.7099, 3.43727},
		{80.3111, 4.32788},
		{82.0424, 3.03992},
		{83.1203, 2.90103},
		{83.4035, 3.17767},
		{84.0336, 2.95874},
		{84.261, 2.91611},
		{83.7614, 2.94033},
		{84.3689, 2.88425},
		{85.8647, 2.72095},
		{86.5072, 2.78788},
		{86.2993, 2.77627},
		{87.315, 2.8253},
		{85.6197, 3.11107},
		{76.9014, 10},		// bad
		{86.981, 1.51767},	// bad
	},
	{
		{78.6776, 3.12589},
		{78.4785, 4.5388},
		{79.2403, 2.91489},
		{80.6575, 2.99326},
		{80.3583, 3.14053},
		{81.7984, 2.75277},
		{81.6708, 2.9812},
		{80.5219, 2.99007},
		{81.1745, 2.6778},
		{81.5756, 2.65063},
		{81.0937, 2.59421},
		{79.9447, 2.55721},
		{82.2385, 2.6701},
		{85.057, 3.54485},
		{79.5521, 1},		// bad
		{80.3022, 10},		// bad
	},
	{
		{81.5867, 3.82439},
		{81.9336, 3.59812},
		{83.5966, 3.48835},
		{84.1434, 3.41045},
		{84.9714, 3.21931},
		{84.8308, 3.30384},
		{86.0818, 3.0832},
		{85.9711, 2.93532},
		{84.0966, 2.76215},
		{86.8203, 2.88312},
		{88.28, 2.70538},
		{88.6277, 2.82232},
		{88.3752, 3.09797},
		{87.8112, 3.27379},
		{83.6053, 1},		// bad
		{76.4577, 10},		// bad
	},
	{
		{75.7813, 2.92513},
		{77.0914, 3.07556},
		{79.8701, 2.95605},
		{80.8612, 3.32668},
		{82.1127, 3.34133},
		{81.6932, 3.57405},
		{82.2699, 2.8817},
		{82.8297, 2.80227},
		{83.643, 2.97491},
		{84.0323, 2.69008},
		{84.8532, 2.72008},
		{83.6818, 2.98885},
		{83.9002, 2.70891},
		{79.6825, 2.03753},
		{68.7654, 9.99711},	// bad
		{90, 10},			// bad
	},
	{
		{85, 3},			// bad
		{85, 3},			// bad
		{85, 3},			// bad
		{85, 3},			// bad
		{83.2193, 3.61056},
		{83.8147, 3.64171},
		{84.6001, 3.30656},
		{84.8716, 3.29717},
		{83.6797, 3.29011},
		{84.6914, 2.97043},
		{85.611, 2.71875},
		{84.7213, 3.40686},
		{84.4221, 2.99822},
		{84.3533, 3.50011},
		{71.0885, 10},		// bad
		{85, 3},			// bad
	},
	{
		{85, 3},			// bad
		{85, 3},			// bad
		{85, 3},			// bad
		{85, 3},			// bad
		{85.5228, 3.41961},
		{85.5448, 3.44281},
		{85.7701, 3.60119},
		{86.6629, 3.30912},
		{86.9838, 3.1692},
		{88.4769, 3.03806},
		{89.443, 3.01299},
		{90.4083, 3.52338},
		{90.2928, 3.03731},
		{89.4199, 2.98574},
		{72.9787, 10},		// bad
		{85, 3},			// bad
	}
};


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run taf_index\n"
		"  run               Set run number.\n"
		"  taf_index         Set TAF index, 0-5.\n"
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

	if (pos_start + 1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss taf_index argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	int run = atoi(argv[pos_start]);
	// taf index
	int taf_index = atoi(argv[pos_start+1]);

	// input file name
	TString input_file_name = TString::Format(
		"%s%staf%d-telescope-%s%04d.root",
		kGenerateDataPath,
		kTelescopeDir,
		taf_index,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	TaEvent taf_event;
	// setup input branches
	taf_event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%staf%d-particle-%sv2-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		taf_index,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of d-sigma
	TH1F hist_d_sigma("hds", "d sigma in straight PID", 1000, -10, 10);
	// output tree
	TTree opt("tree", "rebuilt taf deutron particle");
	// output event
	ParticleEvent particle;
	// setup output branches
	particle.SetupOutput(&opt);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		// initialize
		particle.num = 0;
		particle.px[0] = 0.0;
		particle.py[0] = 0.0;
		particle.pz[0] = 0.0;
		particle.status[0] = 0;

		// check number
		if (taf_event.num != 1) {
			opt.Fill();
			continue;
		}
		// check flag
		if (taf_event.flag[0] != 0x3 && taf_event.flag[0] != 0x5) {
			opt.Fill();
			continue;
		}
		// check strips
		if (
			taf_event.front_strip[0] >= 14
			|| (taf_index == 5 && taf_event.front_strip[0] <= 3)
		) {
			opt.Fill();
			continue;
		}

		int csi_index = taf_event.flag[0] == 0x3 ? 0 : 1;
		int fs = taf_event.front_strip[0];
		double de = taf_event.energy[0][0];
		double e = taf_event.energy[0][1];
		double a = taf_straight_strip_pars[taf_index*2+csi_index][fs][0];
		double b = taf_straight_strip_pars[taf_index*2+csi_index][fs][1];
		double ef = sqrt(de*e + a*de*de) + b*e;
		double mean = taf_strip_fit_d[taf_index*2+csi_index][fs][0];
		double sigma = taf_strip_fit_d[taf_index*2+csi_index][fs][1];
		if (fabs((ef - mean) / sigma) < 2.0) {
			particle.num = 1;
			// charge and mass
			particle.charge[0] = 1;
			particle.mass[0] = 2;
			// energy
			double a0 = power_csi_param[taf_index*2+csi_index][0];
			double a1 = power_csi_param[taf_index*2+csi_index][1];
			double a2 = power_csi_param[taf_index*2+csi_index][2];
			double csi_energy = pow((e - a2 ) / a0, 1.0 / a1);
			particle.energy[0] = de + csi_energy;
			// time
			particle.time[0] = taf_event.time[0][0];
			// position
			ROOT::Math::Polar3DVector position(
				taf_event.radius[0], taf_event.theta[0], taf_event.phi[0]
			);
			particle.x[0] = position.X();
			particle.y[0] = position.Y();
			particle.z[0] = position.Z();
			// index
			particle.index[0] = csi_index;
		}

		opt.Fill();

		hist_d_sigma.Fill((ef-mean)/sigma);
	}
	
	// save tree
	hist_d_sigma.Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}