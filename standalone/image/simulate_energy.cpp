/*
 * simulate_energy.cpp 用来画模拟之后的重建的能量和设定的差别
 */

// #define BE
// #define HE
// #define H
#define EH

void simulate_energy() {
	TString channel_file_name = "/mnt/d/data/ribll2022/channel/C14-10Be-4He-2H-v2-sim.root";
	TFile *channel_file = new TFile(channel_file_name, "read");
	TTree *channel_tree = (TTree*)channel_file->Get("tree");
	int valid;
	bool tafd_edge;
	long long entry;
	double fragment_kinetic[2];
	double recoil_kinetic;
	channel_tree->SetBranchAddress("valid", &valid);
	channel_tree->SetBranchAddress("tafd_edge", &tafd_edge);
	channel_tree->SetBranchAddress("entry", &entry);
	channel_tree->SetBranchAddress("fragment_kinetic", fragment_kinetic);
	channel_tree->SetBranchAddress("recoil_kinetic", &recoil_kinetic);
	// Open the ROOT file
	TString generate_file_name = "/mnt/d/data/ribll2022/simulate/generate-0002.root";
    TFile *generate_file = new TFile(generate_file_name, "read");
	TTree *generate_tree = (TTree*)generate_file->Get("tree");
	double gxfkit, gfkit[2], grkit, gfkat[2], grkat;
	generate_tree->SetBranchAddress("excited_fragment0_kinetic_in_target", &gxfkit);
	generate_tree->SetBranchAddress("fragment_kinetic_in_target", gfkit);
	generate_tree->SetBranchAddress("recoil_kinetic_in_target", &grkit);
	generate_tree->SetBranchAddress("fragment_kinetic_after_target", gfkat);
	generate_tree->SetBranchAddress("recoil_kinetic_after_target", &grkat);

	TH1F *hbek = new TH1F("hbek", "energy difference of 10Be", 100, -10, 10);
	TH1F *hhek = new TH1F("hhek", "energy difference of 4He", 100, -10, 10);
	TH1F *hdk = new TH1F("hdk", "energy difference of 2H", 100, -10, 10);
	TH1F *hedk = new TH1F("hedk", "energy difference of edge 2H", 100, -10, 10);
	for (long long e = 0; e < channel_tree->GetEntriesFast(); ++e) {
		channel_tree->GetEntry(e);
		if (valid != 0) continue;
		generate_tree->GetEntry(entry);
		hbek->Fill(fragment_kinetic[0] - gfkit[0]);
		hhek->Fill(fragment_kinetic[1] - gfkit[1]);
		hdk->Fill(recoil_kinetic - grkit);
		if (tafd_edge) hedk->Fill(recoil_kinetic - grkit);
	}

#ifdef BE
	hbek->SetStats(0);
	hbek->Draw();
#endif
#ifdef HE
	hhek->SetStats(0);
	hhek->Draw();
#endif
#ifdef H
	hdk->SetStats(0);
	hdk->Draw();
#endif
#ifdef EH
	hedk->SetStats(0);
	hedk->Draw();
#endif


}
