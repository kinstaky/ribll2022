#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;

constexpr double taf_straight_pars[12][2] = {
	{319, -0.00163},
	{308, -0.00161},
	{292, -0.00159},
	{333, -0.00141},
	{393, -0.00132},
	{351, -0.00136},
	{376, -0.00139},
	{357, -0.00147},
	{427, -0.00123},
	{391, -0.00122},
	{439, -0.00122},
	{435, -0.00129}
};


constexpr double taf_fit_d[12][2] = {
	{65.4205, 2.70253},
	{64.2327, 2.87975},
	{66.7161, 3.42465},
	{69.8217, 3.48572},
	{74.7281, 3.0812},
	{73.1064, 3.49608},
	{73.9381, 3.25488},
	{71.3212, 3.72456},
	{75.6524, 3.13069},
	{72.1416, 3.05611},
	{72.9958, 3.11807},
	{75.4693, 2.98265}
};


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


inline double VelocityFromMomentum(double momentum, double mass) {
	return momentum / sqrt(pow(momentum, 2.0) + pow(mass, 2.0));
}

inline double GammaFromMomentum(double momentum, double mass) {
	return sqrt(1.0 + pow(momentum/mass, 2.0));
}


void CenterMomentum(
	ROOT::Math::XYZVector momentum1,
	ROOT::Math::XYZVector momentum2,
	double mass1,
	double mass2,
	ROOT::Math::XYZVector &momentum1_center,
	ROOT::Math::XYZVector &momentum2_center
) {
	// gamma
	double gamma1 = GammaFromMomentum(momentum1.R(), mass1);
	double gamma2 = GammaFromMomentum(momentum2.R(), mass2);
	// total energy
	double total_energy1 = gamma1 * mass1;
	double total_energy2 = gamma2 * mass2;


	// mass of center of mass
	double effect_center_mass = gamma1 * mass1 + gamma2 * mass2;
	// center of mass velocity
	ROOT::Math::XYZVector center_velocity =
		(momentum1 + momentum2) / effect_center_mass;
	// gamma of center of mass
	double center_gamma = 1.0 / sqrt(1.0 - pow(center_velocity.R(), 2.0));

	// momentum1 project to center of mass momentum
	ROOT::Math::XYZVector momentum1_parallel =
		center_velocity.Unit().Dot(momentum1) * center_velocity.Unit();
	// momentum1 orthometric to center of mass momentum
	ROOT::Math::XYZVector momentum1_ortho = momentum1 - momentum1_parallel;

	// momentum2 project to center of mass momentum
	ROOT::Math::XYZVector momentum2_parallel =
		center_velocity.Unit().Dot(momentum2) * center_velocity.Unit();
	// momentum2 orthometric to center of mass momentum
	ROOT::Math::XYZVector momentum2_ortho = momentum2 - momentum2_parallel;

	// transmission to center of mass frame
	ROOT::Math::XYZVector momentum1_center_parallel =
		(
			center_gamma * momentum1_parallel.R()
			- sqrt(pow(center_gamma, 2.0) - 1.0) * total_energy1
		) * momentum1_parallel.Unit();
	ROOT::Math::XYZVector momentum2_center_parallel =
		(
			center_gamma * momentum2_parallel.R()
			- sqrt(pow(center_gamma, 2.0) - 1.0) * total_energy2
		) * momentum2_parallel.Unit();

	// reconstruct momentum in center of mass
	momentum1_center = momentum1_center_parallel + momentum1_ortho;
	momentum2_center = momentum2_center_parallel + momentum2_ortho;
}


int main(int argc, char **argv) {
	if (argc > 2) {
		std::cout << "Usage: " << argv[0] << " [tag]\n"
			<< "  tag    9Be or 10Be\n";
		return -1;
	}
	std::string tag = "10Be";
	if (argc == 2) tag = std::string(argv[1]);

	// input file name
	TString info_file_name = TString::Format(
		"%s%sthreebody-%s.root",
		kGenerateDataPath, kInformationDir, tag.c_str()
	);
	// input file
	TFile ipf(info_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< info_file_name << " failed.\n";
		return -1;
	}
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// spectrum V2 file name
	TString spectrum_v2_file_name = TString::Format(
		"%s%sthreebody-%s-2.root", kGenerateDataPath, kSpectrumDir, tag.c_str()
	);
	// add friend
	ipt->AddFriend("s=tree", spectrum_v2_file_name);
	// input data
	double spectrum_c_kinetic[4], spectrum_c_momentum[4];
	double spectrum_q[4];
	int spectrum_be_state[4];
	double spectrum_excited_energy[4];
	// setup input branches
	ipt->SetBranchAddress("s.c_kinetic", spectrum_c_kinetic);
	ipt->SetBranchAddress("s.c_momentum", spectrum_c_momentum);
	ipt->SetBranchAddress("s.q", spectrum_q);
	ipt->SetBranchAddress("s.be_state", spectrum_be_state);
	ipt->SetBranchAddress("s.excited_energy_target", spectrum_excited_energy);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-extra-%s.root",
		kGenerateDataPath, kInformationDir, tag.c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "extra threebody information");
	// output data
	bool valid;
	int hole;
	int layer[2];
	double cd_angle;
	double cd_velocity;
	double cd_center_velocity;
	bool in_target;
	bool in_center;
	double c14_kinetic_sigma;
	double c14_momentum_sigma;
	double q, c14_excited_energy;
	int be_state;
	int csi_index;
	double d_ef, d_sigma;
	double d_strip_ef, d_strip_sigma;
	bool narrow_correlation;
	// setup output branches
	opt.Branch("valid", &valid, "valid/O");
	opt.Branch("hole", &hole, "hole/I");
	opt.Branch("layer", layer, "layer[2]/I");
	opt.Branch("cd_angle", &cd_angle, "cda/D");
	opt.Branch("cd_velocity", &cd_velocity, "cdv/D");
	opt.Branch("cd_center_velocity", &cd_center_velocity, "cdcv/D");
	opt.Branch("in_target", &in_target, "it/O");
	opt.Branch("in_center", &in_center, "ic/O");
	opt.Branch("c14_kinetic_sigma", &c14_kinetic_sigma, "c14ks/D");
	opt.Branch("c14_momentum_sigma", &c14_momentum_sigma, "c14ps/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("c14_excited_energy", &c14_excited_energy, "c14ex/D");
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("d_ef", &d_ef, "def/D");
	opt.Branch("d_sigma", &d_sigma, "ds/D");
	opt.Branch("d_strip_ef", &d_strip_ef, "dsef/D");
	opt.Branch("d_strip_sigma", &d_strip_sigma, "dss/D");
	opt.Branch("narrow_correlation", &narrow_correlation, "nc/O");

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		// initialize
		valid = true;
		hole = 0;
		hole |= event.hole[0] ? 1 : 0;
		hole |= event.hole[1] ? 2 : 0;
		layer[0] = event.layer[0];
		layer[1] = event.layer[1];
		in_target = false;
		in_center = false;
		c14_kinetic_sigma = 10.0;
		c14_momentum_sigma = 10.0;

		// check PPAC
		if (
			(event.ppac_flag & 1) == 0
			|| event.xppac_track[0] == 0
			|| event.xppac_track[1] == 0
		) valid = false;
		// check TAF flag
		if (event.taf_flag != 0) valid = false;
		// check binding events
		if (event.bind != 0) valid = false;

		if (!valid) {
			opt.Fill();
			continue;
		}

		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
		// 10Be momentum vector direction
		ROOT::Math::XYZVector d_be(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		d_be = d_be.Unit();
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be = d_be * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
		// 4He momentum vector direction
		ROOT::Math::XYZVector d_he(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		d_he = d_he.Unit();
		// 4He momentum vector
		ROOT::Math::XYZVector p_he = d_he * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
		// 2H momentum vector
		ROOT::Math::XYZVector d_d(
			event.d_x - event.xptx,
			event.d_y - event.xpty,
			135.0
		);
		d_d = d_d.Unit();
		// 2H momoentum vector
		ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// excited 14C momentum vector
		ROOT::Math::XYZVector p_excited_c = p_be + p_he;
		// excited 14C momentum
		double excited_c_momentum = p_excited_c.R();

		// 14C and 2H angle in lab frame
		cd_angle = acos(p_excited_c.Unit().Dot(p_d.Unit()));

		// excited 14C velocity value

		double velocity_c = VelocityFromMomentum(excited_c_momentum, mass_14c);
		// excited 14C velocity vector
		ROOT::Math::XYZVector v_c = p_excited_c.Unit() * velocity_c;

		// 2H velocity value
		double velocity_d = VelocityFromMomentum(d_momentum, mass_2h);
		// 2H velocity vector
		ROOT::Math::XYZVector v_d = p_d.Unit() * velocity_d;

		// relative velocity vector
		ROOT::Math::XYZVector cd_velocity_vec = v_c - v_d;
		cd_velocity = cd_velocity_vec.R();

		// transfrom to center of mass coordinate
		ROOT::Math::XYZVector pcc, pdc;
		CenterMomentum(p_excited_c, p_d, mass_14c, mass_2h, pcc, pdc);

		// get velocity
		ROOT::Math::XYZVector vcc = pcc.Unit() * VelocityFromMomentum(pcc.R(), mass_14c);
		ROOT::Math::XYZVector vdc = pdc.Unit() * VelocityFromMomentum(pdc.R(), mass_2h);

		// relative velocity
		cd_center_velocity = (vcc - vdc).R();

		// check reaction position
		in_target = (pow(event.xptx-2.0, 2.0) + pow(event.xpty+1.0, 2.0)) < 200.0;
		in_center =
			event.xptx > -10 && event.xptx < -5
			&& event.xpty > -3 && event.xpty < 3;

		// get kinetic energy and momentum
		double c14_kinetic = spectrum_c_kinetic[3];
		double c14_momentum = spectrum_c_momentum[3];
		if (in_center) {
			c14_kinetic_sigma = (c14_kinetic - 376.0) / 7.4;
			c14_momentum_sigma = (c14_momentum - 3147.0) / 37.3;
		} else {
			c14_kinetic_sigma = (c14_kinetic - 384.0) / 4.6;
			c14_momentum_sigma = (c14_momentum - 3184.0) / 23.7;
		}
		q = spectrum_q[3];
		be_state = spectrum_be_state[3];
		c14_excited_energy = spectrum_excited_energy[3];

		csi_index = event.csi_index;
		// CsI direction from origin point
		ROOT::Math::XYZVector d_cd(event.d_x, event.d_y, 135.0);
		// thickness corrected TAFD energy
		double taf_cde = event.tafd_energy * cos(d_cd.Theta());
		// CsI energy
		double taf_e = event.csi_channel;
		// straight parameters
		double taf_a = taf_straight_pars[csi_index][0];
		double taf_b = taf_straight_pars[csi_index][1];
		// fixed energy
		d_ef =
			sqrt(taf_cde*taf_e + taf_a*taf_cde*taf_cde) + taf_b*taf_e;
		// sigma
		d_sigma = (d_ef - taf_fit_d[csi_index][0]) / taf_fit_d[csi_index][1];

		// single strip
		double taf_de = event.tafd_energy;
		int taf_fs = event.d_x_strip;
		double taf_as = taf_straight_strip_pars[csi_index][taf_fs][0];
		double taf_bs = taf_straight_strip_pars[csi_index][taf_fs][1];
		d_strip_ef =
			sqrt(taf_de*taf_e + taf_as*taf_de*taf_de) + taf_bs*taf_e;
		d_strip_sigma = (d_strip_ef - taf_strip_fit_d[csi_index][taf_fs][0])
			/ taf_strip_fit_d[csi_index][taf_fs][1];

		narrow_correlation = true;
		if (event.bind != 0) narrow_correlation = false;
		if (event.hole[0] || event.hole[1]) narrow_correlation = false;
		// 10Be
		if (event.be_x_hit[0] != 1 || event.be_y_hit[0] != 1) {
			narrow_correlation = false;
		} else if (
			event.be_x_channel[0][0] - event.be_y_channel[0][0] < -220
			|| event.be_x_channel[0][0] - event.be_y_channel[0][0] > 150
		) {
			narrow_correlation = false;
		}
		if (event.layer[0] == 2) narrow_correlation = false;
		// 4He T0D1
		if (event.he_x_hit[0] !=1 || event.he_y_hit[0] != 1) {
			narrow_correlation = false;
		} else if (
			event.he_x_channel[0][0] - event.he_y_channel[0][0] < -274
			&& event.he_x_channel[0][0] - event.he_y_channel[0][0] > 170
		) {
			narrow_correlation = false;
		}
		// 4He T0D2
		if (event.he_x_hit[1] != 1 || event.he_y_hit[1] != 1) {
			narrow_correlation = false;
		} else if (
			event.he_x_channel[1][0] - event.he_y_channel[1][0] < -236
			|| event.he_x_channel[1][0] - event.he_y_channel[1][0] > 256
		) {
			narrow_correlation = false;
		}
		// 4He T0D3
		if (event.layer[1] > 1) {
			if (event.he_x_hit[2] != 1 || event.he_y_hit[2] != 1) {
				narrow_correlation = false;
			} else if (
				event.he_x_channel[2][0] - event.he_y_channel[2][0] < -64
				|| event.he_x_channel[2][0] - event.he_y_channel[2][0] > 120
			) {
				narrow_correlation = false;
			}
		}

		// fill
		opt.Fill();
	}

	// save
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}