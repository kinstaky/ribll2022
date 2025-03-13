void csi_error_fit() {
	TGraph *g = new TGraph;
	g->AddPoint(10, 0.5);
	g->AddPoint(15, 0.75);
	g->AddPoint(20, 1.04);
	g->AddPoint(25, 1.42);
	g->AddPoint(30, 1.76);
	g->AddPoint(35, 2);
	// fit
	TF1 *f = new TF1("f", "pol1", 5, 40);
	g->Fit(f, "R+");
	// change style
	g->SetMarkerStyle(33);
	g->SetMarkerSize(3);
	g->Draw("AP");
}