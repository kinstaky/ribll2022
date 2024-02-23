#include "include/event/threebody_info_event.h"

namespace ribll {

void ThreeBodyInfoEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"csi_index").c_str(), &csi_index);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"taf_flag").c_str(), &taf_flag);
	tree->SetBranchAddress((prefix+"csi_channel").c_str(), &csi_channel);
	tree->SetBranchAddress((prefix+"be_channel").c_str(), be_channel);
	tree->SetBranchAddress((prefix+"he_channel").c_str(), he_channel);
	tree->SetBranchAddress((prefix+"ssd_channel").c_str(), ssd_channel);
	tree->SetBranchAddress((prefix+"t0_energy").c_str(), t0_energy);
	tree->SetBranchAddress((prefix+"tafd_energy").c_str(), &tafd_energy);
	tree->SetBranchAddress((prefix+"csi_energy").c_str(), &csi_energy);
	tree->SetBranchAddress((prefix+"taf_energy").c_str(), &taf_energy);
	tree->SetBranchAddress((prefix+"be_x").c_str(), be_x);
	tree->SetBranchAddress((prefix+"be_y").c_str(), be_y);
	tree->SetBranchAddress((prefix+"he_x").c_str(), he_x);
	tree->SetBranchAddress((prefix+"he_y").c_str(), he_y);
	tree->SetBranchAddress((prefix+"d_x").c_str(), &d_x);
	tree->SetBranchAddress((prefix+"d_y").c_str(), &d_y);
	tree->SetBranchAddress((prefix+"xptx").c_str(), &xptx);
	tree->SetBranchAddress((prefix+"xpty").c_str(), &xpty);
	tree->SetBranchAddress((prefix+"vptx").c_str(), &vptx);
	tree->SetBranchAddress((prefix+"vpty").c_str(), &vpty);
	tree->SetBranchAddress((prefix+"ppac_flag").c_str(), &ppac_flag);
	tree->SetBranchAddress((prefix+"xppac_xflag").c_str(), &xppac_xflag);
	tree->SetBranchAddress((prefix+"xppac_yflag").c_str(), &xppac_yflag);
	tree->SetBranchAddress((prefix+"vppac_xflag").c_str(), &vppac_xflag);
	tree->SetBranchAddress((prefix+"vppac_yflag").c_str(), &vppac_yflag);
	tree->SetBranchAddress((prefix+"xppac_x").c_str(), xppac_x);
	tree->SetBranchAddress((prefix+"xppac_y").c_str(), xppac_y);
	tree->SetBranchAddress((prefix+"vppac_x").c_str(), vppac_x);
	tree->SetBranchAddress((prefix+"vppac_y").c_str(), vppac_y);
	tree->SetBranchAddress((prefix+"xppac_track").c_str(), xppac_track);
	tree->SetBranchAddress((prefix+"vppac_track").c_str(), vppac_track);
	tree->SetBranchAddress((prefix+"be_x_hit").c_str(), be_x_hit);
	tree->SetBranchAddress((prefix+"be_y_hit").c_str(), be_y_hit);
	tree->SetBranchAddress((prefix+"he_x_hit").c_str(), he_x_hit);
	tree->SetBranchAddress((prefix+"he_y_hit").c_str(), he_y_hit);
	tree->SetBranchAddress((prefix+"be_x_channel").c_str(), be_x_channel);
	tree->SetBranchAddress((prefix+"be_y_channel").c_str(), be_y_channel);
	tree->SetBranchAddress((prefix+"he_x_channel").c_str(), he_x_channel);
	tree->SetBranchAddress((prefix+"he_y_channel").c_str(), he_y_channel);
	tree->SetBranchAddress((prefix+"d_x_channel").c_str(), &d_x_channel);
	tree->SetBranchAddress((prefix+"d_y_channel").c_str(), &d_y_channel);
	tree->SetBranchAddress((prefix+"be_x_time").c_str(), be_x_time);
	tree->SetBranchAddress((prefix+"be_y_time").c_str(), be_y_time);
	tree->SetBranchAddress((prefix+"he_x_time").c_str(), he_x_time);
	tree->SetBranchAddress((prefix+"he_y_time").c_str(), he_y_time);
	tree->SetBranchAddress((prefix+"d_x_time").c_str(), &d_x_time);
	tree->SetBranchAddress((prefix+"d_y_time").c_str(), &d_y_time);
	tree->SetBranchAddress((prefix+"be_x_strip").c_str(), be_x_strip);
	tree->SetBranchAddress((prefix+"be_y_strip").c_str(), be_y_strip);
	tree->SetBranchAddress((prefix+"he_x_strip").c_str(), he_x_strip);
	tree->SetBranchAddress((prefix+"he_y_strip").c_str(), he_y_strip);
	tree->SetBranchAddress((prefix+"d_x_strip").c_str(), &d_x_strip);
	tree->SetBranchAddress((prefix+"d_y_strip").c_str(), &d_y_strip);
	tree->SetBranchAddress((prefix+"c14_kinetic").c_str(), &c14_kinetic);
	tree->SetBranchAddress((prefix+"q").c_str(), &q);
	tree->SetBranchAddress((prefix+"vq").c_str(), &vq);
	tree->SetBranchAddress((prefix+"hole").c_str(), hole);
	tree->SetBranchAddress((prefix+"run").c_str(), &run);
	tree->SetBranchAddress((prefix+"entry").c_str(), &entry);
}


void ThreeBodyInfoEvent::SetupOutput(TTree *tree) {
	tree->Branch("csi_index", &csi_index, "ci/I");
	tree->Branch("layer", layer, "layer[2]/I");
	tree->Branch("taf_flag", &taf_flag, "tafflag/I");
	tree->Branch("csi_channel", &csi_channel, "csic/D");
	tree->Branch("be_channel", be_channel, "bec[3]/D");
	tree->Branch("he_channel", he_channel, "hec[3]/D");
	tree->Branch("ssd_channel", ssd_channel, "ssdc[3]/D");
	tree->Branch("t0_energy", t0_energy, "t0e[2]/D");
	tree->Branch("tafd_energy", &tafd_energy, "tafde/D");
	tree->Branch("csi_energy", &csi_energy, "csie/D");
	tree->Branch("taf_energy", &taf_energy, "tafe/D");
	tree->Branch("be_x", be_x, "bex[3]/D");
	tree->Branch("be_y", be_y, "bey[3]/D");
	tree->Branch("he_x", he_x, "hex[3]/D");
	tree->Branch("he_y", he_y, "hey[3]/D");
	tree->Branch("d_x", &d_x, "dx/D");
	tree->Branch("d_y", &d_y, "dy/D");
	tree->Branch("xptx", &xptx, "xptx/D");
	tree->Branch("xpty", &xpty, "xpty/D");
	tree->Branch("vptx", &vptx, "vptx/D");
	tree->Branch("vpty", &vpty, "vpty/D");
	tree->Branch("ppac_flag", &ppac_flag, "pflag/I");
	tree->Branch("xppac_xflag", &xppac_xflag, "xpxflag/s");
	tree->Branch("xppac_yflag", &xppac_yflag, "xpyflag/s");
	tree->Branch("vppac_xflag", &vppac_xflag, "vpxflag/s");
	tree->Branch("vppac_yflag", &vppac_yflag, "vpyflag/s");
	tree->Branch("xppac_x", xppac_x, "xppacx[3]/D");
	tree->Branch("xppac_y", xppac_y, "xppacy[3]/D");
	tree->Branch("vppac_x", vppac_x, "vppacx[3]/D");
	tree->Branch("vppac_y", vppac_y, "vppacy[3]/D");
	tree->Branch("xppac_track", xppac_track, "xptrack[2]/I");
	tree->Branch("vppac_track", vppac_track, "vptrack[2]/I");
	tree->Branch("be_x_hit", be_x_hit, "bexhit[3]/I");
	tree->Branch("be_y_hit", be_y_hit, "beyhit[3]/I");
	tree->Branch("he_x_hit", he_x_hit, "hexhit[3]/I");
	tree->Branch("he_y_hit", he_y_hit, "heyhit[3]/I");
	tree->Branch("be_x_channel", be_x_channel, "be_x_channel[3][2]/D");
	tree->Branch("be_y_channel", be_y_channel, "be_y_channel[3][2]/D");
	tree->Branch("he_x_channel", he_x_channel, "he_x_channel[3][2]/D");
	tree->Branch("he_y_channel", he_y_channel, "he_y_channel[3][2]/D");
	tree->Branch("d_x_channel", &d_x_channel, "d_x_channel/D");
	tree->Branch("d_y_channel", &d_y_channel, "d_y_channel/D");
	tree->Branch("be_x_time", be_x_time, "be_x_time[3][2]/D");
	tree->Branch("be_y_time", be_y_time, "be_y_time[3][2]/D");
	tree->Branch("he_x_time", he_x_time, "he_x_time[3][2]/D");
	tree->Branch("he_y_time", he_y_time, "he_y_time[3][2]/D");
	tree->Branch("d_x_time", &d_x_time, "d_x_time/D");
	tree->Branch("d_y_time", &d_y_time, "d_y_time/D");
	tree->Branch("be_x_strip", be_x_strip, "be_x_strip[3][2]/i");
	tree->Branch("be_y_strip", be_y_strip, "be_y_strip[3][2]/i");
	tree->Branch("he_x_strip", he_x_strip, "he_x_strip[3][2]/i");
	tree->Branch("he_y_strip", he_y_strip, "he_y_strip[3][2]/i");
	tree->Branch("d_x_strip", &d_x_strip, "d_x_strip/i");
	tree->Branch("d_y_strip", &d_y_strip, "d_y_strip/i");
	tree->Branch("c14_kinetic", &c14_kinetic, "c14k/D");
	tree->Branch("q", &q, "q/D");
	tree->Branch("vq", &vq, "vq/D");
	tree->Branch("hole", hole, "hole[2]/O");
	tree->Branch("run", &run, "run/I");
	tree->Branch("entry", &entry, "entry/L");
}

}	// ribll
