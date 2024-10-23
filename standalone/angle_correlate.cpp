#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"

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
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " tag\n";
		return 0;
	}
	std::string tag(argv[1]);

	// input threebody info file name
	TString input_file_name = TString::Format(
        "%s%sthreebody-%s.root",
        kGenerateDataPath,
		kInformationDir,
		tag.c_str()
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
    ThreeBodyInfoEvent event;
    // setup input branches
    event.SetupInput(ipt);

	// spectrum V2 tree
	TString spectrum_v2_file_name = TString::Format(
        "%s%sthreebody-%s-2.root",
        kGenerateDataPath,
        kSpectrumDir,
        tag.c_str()
    );
	// add friend
	ipt->AddFriend("spectrum=tree", spectrum_v2_file_name);
	// input spectrum V2 data
	// straight flag
	int straight[2];
	// Q value
	double input_q[4];
	// excited energy ignore 10Be state
	double stateless_excited_energy[3][4];
    // setup output branches
    ipt->SetBranchAddress("spectrum.straight", straight);
    ipt->SetBranchAddress("spectrum.q", input_q);
	ipt->SetBranchAddress(
		"spectrum.stateless_excited_energy",
		stateless_excited_energy
	);

	// output file name
	TString output_file_name = TString::Format(
        "%s%sangle-correlation-%s.root",
        kGenerateDataPath,
		kSpectrumDir,
		tag.c_str()
    );
    // open output file
    TFile opf(output_file_name, "recreate");
	// histogram of ground state, divided by cos\psi
	TH1F hist_ex_psi[psi_bins];
	for (int i = 0; i < psi_bins; ++i) {
		hist_ex_psi[i] = TH1F(
			TString::Format("hep%d", i),
			TString::Format(
				"ground state excited energy for |cos(Psi)| from %d/%d to %d/%d ",
				i, psi_bins, i+1, psi_bins
			),
			50, 12.1, 27.1
		);
	}

    // create output tree
    TTree opt("tree", "angle correlation");
	// output data
	bool valid;
	int be_state, c_state;
	double q, excited_energy;
	double beam_theta, beam_phi;
	double recoil_theta, recoil_phi;
	double parent_theta, parent_phi;
	double fragment_theta[2], fragment_phi[2];
	double beam_energy, recoil_energy, parent_energy, fragment_energy[2];
	double recoil_vx, recoil_vy, recoil_vz;
	double parent_vx, parent_vy, parent_vz;
	double fragment_vx[2], fragment_vy[2], fragment_vz[2];
	double parent_recoil_angle, fragments_angle;
	double angle_theta_star;
	double angle_psi;
	double angle_chi;
	// setup output branches
	opt.Branch("valid", &valid, "valid/O");
	opt.Branch("q", &q, "q/D");
	opt.Branch("excited_energy", &excited_energy, "ex/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("c_state", &c_state, "cs/I");
	opt.Branch("beam_theta", &beam_theta, "btheta/D");
	opt.Branch("beam_phi", &beam_phi, "bphi/D");
	opt.Branch("recoil_theta", &recoil_theta, "rtheta/D");
	opt.Branch("recoil_phi", &recoil_phi, "rphi/D");
	opt.Branch("parent_theta", &parent_theta, "ptheta/D");
	opt.Branch("parent_phi", &parent_phi, "pphi/D");
	opt.Branch("fragment_theta", fragment_theta, "ftheta[2]/D");
	opt.Branch("fragment_phi", fragment_phi, "fphi[2]/D");
	opt.Branch("beam_energy", &beam_energy, "be/D");
	opt.Branch("recoil_energy", &recoil_energy, "re/D");
	opt.Branch("parent_energy", &parent_energy, "pe/D");
	opt.Branch("fragment_energy", fragment_energy, "fe[2]/D");
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

		valid = true;
		if ((event.ppac_flag & 1) == 0) valid = false;
		if (event.taf_flag != 0) valid = false;
		if (event.target_flag != 1) valid = false;
		if (event.bind != 0) valid = false;
		if (event.hole[0] || event.hole[1]) valid = false;
		if (straight[0] != 3) valid = false;

		q = input_q[0];
		be_state = -1;
		c_state = -1;
		excited_energy = 0.0;
		if (q < -11 && q > -13) {
			be_state = 0;
			excited_energy = stateless_excited_energy[0][0];
			if (excited_energy > 18.0 && excited_energy < 18.6) {
				c_state = 183;
			}
		} else if (q < -14.5 && q > -16) {
			be_state = 1;
			excited_energy = stateless_excited_energy[1][0];
		} else if (q < -17 && q > -20) {
			be_state = 2;
			excited_energy = stateless_excited_energy[2][0];
		}

		// recoil particle direction
		ROOT::Math::XYZVector recoil_direction(
            event.d_x - event.xptx,
			event.d_y - event.xpty,
			135.0
        );
        recoil_direction = recoil_direction.Unit();
        // recoil_theta = recoil_direction.Theta();
        // recoil_phi = recoil_direction.Phi();

        // fragment1 particle directions
		ROOT::Math::XYZVector fragment1_direction(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		fragment1_direction = fragment1_direction.Unit();
		// fragment_theta[0] = fragment1_direction.Theta();
		// fragment_phi[0] = fragment1_direction.Phi();

        // fragment2 particle directions
		ROOT::Math::XYZVector fragment2_direction(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		fragment2_direction = fragment2_direction.Unit();
		// fragment_theta[1] = fragment2_direction.Theta();
		// fragment_phi[1] = fragment2_direction.Phi();

		// calculate kinetic energy and total energy
		// recoil particle
		double recoil_kinetic = event.taf_energy;
		recoil_energy = recoil_mass + recoil_kinetic;
		// fragment particle 1
		double fragment1_kinetic = event.t0_energy[0];
        fragment_energy[0] = fragment_mass[0] + fragment1_kinetic;
        // fragment particle 2
        double fragment2_kinetic = event.t0_energy[1];
        fragment_energy[1] = fragment_mass[1] + fragment2_kinetic;

		// get momentum
		// recoil particle
		double recoil_momentum =
			MomentumFromKinetic(recoil_mass, recoil_kinetic);
		ROOT::Math::XYZVector recoil_p =
			recoil_direction * recoil_momentum;
		// fragment particle 1
		double fragment1_momentum =
			MomentumFromKinetic(fragment_mass[0], fragment1_kinetic);
		ROOT::Math::XYZVector fragment1_p =
			fragment1_direction * fragment1_momentum;
		// fragment particle 2
		double fragment2_momentum =
			MomentumFromKinetic(fragment_mass[1], fragment2_kinetic);
		ROOT::Math::XYZVector fragment2_p =
			fragment2_direction * fragment2_momentum;
		// beam particle
		ROOT::Math::XYZVector beam_p = recoil_p + fragment1_p + fragment2_p;
		// parent particle
		ROOT::Math::XYZVector parent_p = fragment1_p + fragment2_p;

		// get beam direction
		ROOT::Math::XYZVector beam_direction = beam_p.Unit();

		// get beam energy
		beam_energy = sqrt(pow(beam_p.R(), 2.0) + pow(beam_mass, 2.0));

		// velocity
		ROOT::Math::XYZVector parent_v, recoil_v;
		RelativeVelocity(
			parent_p, recoil_p,
			parent_mass, recoil_mass,
			parent_v, recoil_v
		);
		ROOT::Math::XYZVector fragment1_v, fragment2_v;
		RelativeVelocity(
			fragment1_p, fragment2_p,
            fragment_mass[0], fragment_mass[1],
            fragment1_v, fragment2_v
		);

		// relative velocity
		ROOT::Math::XYZVector parent_recoil_velocity = parent_v - recoil_v;
		ROOT::Math::XYZVector fragments_velocity = fragment1_v - fragment2_v;

		// angle between parent and recoil velocity
		parent_recoil_angle = parent_recoil_velocity.Theta();
		// angle between fragments' velocity
		fragments_angle = fragments_velocity.Theta();

		// velocity angle between relative velocity and beam
		// 母核、反冲核的相对速度矢量与入射束流方向矢量的夹角，韩家兴论文里的 θ*
		angle_theta_star =
			acos(parent_recoil_velocity.Unit().Dot(beam_direction));
		// angle between fragments relative velocity and beam direction
		// 衰变子核的相对速度矢量和入射束流方向的夹角，韩家兴论文里的 Psi 角
		angle_psi =
            acos(fragments_velocity.Unit().Dot(beam_direction));
		// 反应平面和衰变平面的夹角，即两个平面的法向量的夹角
		angle_chi = fabs(acos(
			parent_recoil_velocity.Cross(beam_direction).Unit().Dot(
				fragments_velocity.Cross(beam_direction).Unit()
			)
		));
		if (angle_chi > pi) angle_chi = 2*pi - angle_chi;

		// record velocity
		recoil_vx = recoil_v.X();
		recoil_vy = recoil_v.Y();
		recoil_vz = recoil_v.Z();
		parent_vx = parent_v.X();
		parent_vy = parent_v.Y();
		parent_vz = parent_v.Z();
		fragment_vx[0] = fragment1_v.X();
		fragment_vy[0] = fragment1_v.Y();
		fragment_vz[0] = fragment1_v.Z();
		fragment_vx[1] = fragment2_v.X();
		fragment_vy[1] = fragment2_v.Y();
		fragment_vz[1] = fragment2_v.Z();

		// fill to tree
		opt.Fill();

		// fill to histograms
		int index = int(fabs(cos(angle_psi)) * psi_bins);
		if (index >= 0 && index < psi_bins && valid && be_state == 0) {
			hist_ex_psi[index].Fill(excited_energy);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// output pdf file name
	TString pdf_file_name = TString::Format(
		"%s%sground-state-energy-psi.pdf",
		kGenerateDataPath,
		kImageDir
	);
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();
	c1->Print(pdf_file_name+"[");
	for (int i = 0; i < psi_bins; ++i) {
		hist_ex_psi[i].Draw();
		c1->Print(pdf_file_name);
	}
	c1->Print(pdf_file_name+"]");

	// save
	opf.cd();
	for (int i = 0; i < psi_bins; ++i) hist_ex_psi[i].Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}