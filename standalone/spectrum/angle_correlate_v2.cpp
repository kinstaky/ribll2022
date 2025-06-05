#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <Math/Vector3D.h>

#include "include/defs.h"

using namespace ribll;

constexpr double beam_mass = mass_14c;
constexpr double recoil_mass = mass_2h;
constexpr double parent_mass = mass_14c;
constexpr double fragment_mass[2] = {mass_10be, mass_4he};

constexpr int psi_bins = 6;

void RelativeVelocity(
	const ROOT::Math::XYZVector &p1,
	const ROOT::Math::XYZVector &p2,
	const double mass1,
	const double mass2,
	ROOT::Math::XYZVector &v1,
	ROOT::Math::XYZVector &v2
) {
	ROOT::Math::XYZVector total_p = p1 + p2;
	double energy1 = sqrt(p1.Mag2() + mass1 * mass1);
	double energy2 = sqrt(p2.Mag2() + mass2 * mass2);
	double total_energy = energy1 + energy2;
	ROOT::Math::XYZVector velocity_cm = total_p / total_energy;
	double beta = velocity_cm.R();
	double gamma = 1.0 / sqrt(1.0 - beta * beta);

	double p1_parallel = p1.Dot(velocity_cm.Unit());
	ROOT::Math::XYZVector p1_vertical = p1 - p1_parallel * velocity_cm.Unit();
	double p2_parallel = p2.Dot(velocity_cm.Unit());
	ROOT::Math::XYZVector p2_vertical = p2 - p2_parallel * velocity_cm.Unit();

	double p1_parallel_cm = -gamma*beta*energy1 + gamma*p1_parallel;
	ROOT::Math::XYZVector p1_cm =
		p1_parallel_cm * velocity_cm.Unit() + p1_vertical;
	double p2_parallel_cm = -gamma*beta*energy2 + gamma*p2_parallel;
	ROOT::Math::XYZVector p2_cm =
	    p2_parallel_cm * velocity_cm.Unit() + p2_vertical;

	double energy1_cm = gamma*energy1 - gamma*beta*p1_parallel;
	double energy2_cm = gamma*energy2 - gamma*beta*p2_parallel;

	v1 = p1_cm / energy1_cm;
	v2 = p2_cm / energy2_cm;
	return;
}


int main(int argc, char **argv) {
	// check arguments
	if (argc != 1) {
		std::cout << "Usage: " << argv[0] << "\n";
		return 0;
	}

	// input threebody info file name
	TString input_file_name = TString::Format(
        "%s%sC14-10Be-4He-2H-v3.root",
        kGenerateDataPath,
		kSpectrumDir
    );
    // open input file
    TFile ipf(input_file_name, "read");
    // get input tree
    TTree *ipt = (TTree*)ipf.Get("tree");
    if (!ipt) {
        std::cerr << "Error: Failed to find input tree in file: "
			<< input_file_name << "\n";
        return -1;
	}
    // input event
	int valid;
	double be_kinetic, he_kinetic, d_kinetic;
	double c_momentum, c_kinetic;
	// particle direction
	double be_dx, be_dy, be_dz;
	double he_dx, he_dy, he_dz;
	double d_dx, d_dy, d_dz;
	// Q value
	double q;
	// 10Be state, under different T0D2 energy method
	int be_state;
	// excited energy with target energy lost
	double excited_energy;

	// setup input branches
    ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("be_kinetic", &be_kinetic);
	ipt->SetBranchAddress("he_kinetic", &he_kinetic);
	ipt->SetBranchAddress("d_kinetic", &d_kinetic);
	ipt->SetBranchAddress("c_kinetic", &c_kinetic);
	ipt->SetBranchAddress("c_momentum", &c_momentum);
	ipt->SetBranchAddress("q", &q);
	ipt->SetBranchAddress("be_state", &be_state);
	ipt->SetBranchAddress("excited_energy", &excited_energy);
	ipt->SetBranchAddress("be_dx", &be_dx);
	ipt->SetBranchAddress("be_dy", &be_dy);
	ipt->SetBranchAddress("be_dz", &be_dz);
	ipt->SetBranchAddress("he_dx", &he_dx);
	ipt->SetBranchAddress("he_dy", &he_dy);
	ipt->SetBranchAddress("he_dz", &he_dz);
	ipt->SetBranchAddress("d_dx", &d_dx);
	ipt->SetBranchAddress("d_dy", &d_dy);
	ipt->SetBranchAddress("d_dz", &d_dz);

	// output file name
	TString output_file_name = TString::Format(
        "%s%sangle-correlation-v2.root",
        kGenerateDataPath,
		kSpectrumDir
    );
    // open output file
    TFile opf(output_file_name, "recreate");
	// output tree
    TTree opt("tree", "angle correlation");
	// output data
	double recoil_vx, recoil_vy, recoil_vz;
	double parent_vx, parent_vy, parent_vz;
	double fragment_vx[2], fragment_vy[2], fragment_vz[2];
	double parent_recoil_angle, fragments_angle;
	double angle_theta_star;
	double angle_psi;
	double angle_chi;
	// setup output branches
	opt.Branch("valid", &valid, "valid/I");
	opt.Branch("q", &q, "q/D");
	opt.Branch("excited_energy", &excited_energy, "ex/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("recoil_vx", &recoil_vx, "rvx/D");
	opt.Branch("recoil_vy", &recoil_vy, "rvy/D");
	opt.Branch("recoil_vz", &recoil_vz, "rvz/D");
	opt.Branch("parent_vx", &parent_vx, "pvx/D");
	opt.Branch("parent_vy", &parent_vy, "pvy/D");
	opt.Branch("parent_vz", &parent_vz, "pvz/D");
	opt.Branch("fragment_vx", fragment_vx, "fvx[2]/D");
	opt.Branch("fragment_vy", fragment_vy, "fvy[2]/D");
	opt.Branch("fragment_vz", fragment_vz, "fvz[2]/D");
	opt.Branch("parent_recoil_angle", &parent_recoil_angle, "pra/D");
	opt.Branch("fragments_angle", &fragments_angle, "ffa/D");
	opt.Branch("angle_theta_star", &angle_theta_star, "athetas/D");
	opt.Branch("angle_psi", &angle_psi, "apsi/D");
	opt.Branch("angle_chi", &angle_chi, "achi/D");

	// total number of entries
    long long entries = ipt->GetEntries();
    // 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Calculating angle   0%%");
	fflush(stdout);
	// loop events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
            fflush(stdout);
		}
		// get data
	    ipt->GetEntry(entry);
		if (valid != 0) {
			opt.Fill();
			continue;
		}

		// recoil particle direction
		ROOT::Math::XYZVector d_direction(d_dx, d_dy, d_dz);
        // fragment1 particle directions
		ROOT::Math::XYZVector be_direction(be_dx, be_dy, be_dz);
        // fragment2 particle directions
		ROOT::Math::XYZVector he_direction(he_dx, he_dy, he_dz);
		// // recoil particle
		// double d_energy = mass_2h + d_kinetic;
		// // fragment particle 1
        // double be_energy = mass_10be + be_kinetic;
        // // fragment particle 2
        // double he_energy = mass_4he + he_kinetic;

		// get momentum
		// recoil particle
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
		ROOT::Math::XYZVector dp = d_direction * d_momentum;
		// fragment particle 1
		double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
		ROOT::Math::XYZVector bep = be_direction * be_momentum;
		// fragment particle 2
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
		ROOT::Math::XYZVector hep = he_direction * he_momentum;
		// beam particle
		ROOT::Math::XYZVector cp = dp + bep + hep;
		// parent particle
		ROOT::Math::XYZVector xcp = bep + hep;

		// get beam direction
		ROOT::Math::XYZVector c_direction = cp.Unit();

		// // get beam energy
		// double c_energy = sqrt(pow(cp.R(), 2.0) + pow(mass_14c, 2.0));

		// velocity
		ROOT::Math::XYZVector xcv, dv;
		RelativeVelocity(
			xcp, dp,
			mass_14c, mass_2h,
			xcv, dv
		);
		ROOT::Math::XYZVector bev, hev;
		RelativeVelocity(
			bep, hep,
            mass_10be, mass_4he,
            bev, hev
		);

		// relative velocity
		ROOT::Math::XYZVector parent_recoil_velocity = xcv - dv;
		ROOT::Math::XYZVector fragments_velocity = bev - hev;

		// angle between parent and recoil velocity
		parent_recoil_angle = parent_recoil_velocity.Theta();
		// angle between fragments' velocity
		fragments_angle = fragments_velocity.Theta();

		// velocity angle between relative velocity and beam
		// 母核、反冲核的相对速度矢量与入射束流方向矢量的夹角，韩家兴论文里的 θ*
		angle_theta_star =
			acos(parent_recoil_velocity.Unit().Dot(c_direction));
		// angle between fragments relative velocity and beam direction
		// 衰变子核的相对速度矢量和入射束流方向的夹角，韩家兴论文里的 Psi 角
		angle_psi =
            acos(fragments_velocity.Unit().Dot(c_direction));
		// 反应平面和衰变平面的夹角，即两个平面的法向量的夹角
		angle_chi = fabs(acos(
			parent_recoil_velocity.Cross(c_direction).Unit().Dot(
				fragments_velocity.Cross(c_direction).Unit()
			)
		));
		if (angle_chi > pi) angle_chi = 2*pi - angle_chi;

		// record velocity
		recoil_vx = dv.X();
		recoil_vy = dv.Y();
		recoil_vz = dv.Z();
		parent_vx = xcv.X();
		parent_vy = xcv.Y();
		parent_vz = xcv.Z();
		fragment_vx[0] = bev.X();
		fragment_vy[0] = bev.Y();
		fragment_vz[0] = bev.Z();
		fragment_vx[1] = hev.X();
		fragment_vy[1] = hev.Y();
		fragment_vz[1] = hev.Z();

		// fill to tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// // output pdf file name
	// TString pdf_file_name = TString::Format(
	// 	"%s%sground-state-energy-psi.pdf",
	// 	kGenerateDataPath,
	// 	kImageDir
	// );
	// TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	// c1->cd();
	// c1->Print(pdf_file_name+"[");
	// c1->Print(pdf_file_name+"]");

	// save
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}