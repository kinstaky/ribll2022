#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLine.h>

void rotational_band() {
	const double marker_size = 6;

    // Define x values
    double theory_x_values[] = {0, 2, 4, 6, 8};
    double observe_x_values[] = {0, 2, 4};

    // Define y values
    double theory_pi_values[] = {14.64, 15.73, 17.98, 21.80, 27.25};
    double theory_sigma_values[] = {22.16, 22.93, 24.30, 26.45, 29.39};
    double observe_pi_values[] = {13.9, 14.9, 17.3};
    // double observe_sigma_values[] = {22.3, 23.2, 24.6};
	double observe_sigma_values[] = {22.08, 22.78, 24.41};
    // double observe_sigma_values[] = {22.2, 22.9, 24.4};

    // // Create a canvas
    // TCanvas *c1 = new TCanvas("c1", "Rotational Band", 1200, 800);

    // Create graphs
    TGraph *g1 = new TGraph(5, theory_x_values, theory_pi_values);
    TGraph *g2 = new TGraph(5, theory_x_values, theory_sigma_values);
    TGraph *g3 = new TGraph(3, observe_x_values, observe_pi_values);
    TGraph *g4 = new TGraph(3, observe_x_values, observe_sigma_values);

    // Set marker styles, sizes, and colors
    g1->SetMarkerStyle(89);
    g1->SetMarkerSize(marker_size);
    g1->SetMarkerColor(kGreen);
    g1->SetLineColor(kGreen);
    g1->SetLineStyle(2);
	g1->SetLineWidth(4);

    g2->SetMarkerStyle(90);
    g2->SetMarkerSize(marker_size);
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
	g2->SetLineWidth(4);

    g3->SetMarkerStyle(43);
    g3->SetMarkerSize(marker_size);
    g3->SetMarkerColor(kOrange);

    g4->SetMarkerStyle(33);
    g4->SetMarkerSize(marker_size);
    g4->SetMarkerColor(kRed);

	gStyle->SetOptTitle(0);
    // Draw graphs
    g1->Draw("APL");
    g2->Draw("PL SAME");
    g3->Draw("P SAME");
    g4->Draw("P SAME");

    // // Add a dashed line
    // TLine *line0 = new TLine(-0.5, 12.0125, 1.5, 12.0125);
    // line0->SetLineStyle(2); // Dashed line
    // line0->Draw();

	// TLine *line1 = new TLine(-0.5, 12.0125+3.368, 1.5, 12.0125+3.368);
    // line1->SetLineStyle(2); // Dashed line
    // line1->Draw();

	// TLine *line2 = new TLine(-0.5, 12.0125+6.179, 1.5, 12.0125+6.179);
    // line2->SetLineStyle(2); // Dashed line
    // line2->Draw();


    // Set axis labels and limits
    g1->GetXaxis()->SetLimits(-0.5, 8.5);
    g1->GetYaxis()->SetRangeUser(8, 35);
    // g1->GetXaxis()->SetTitle("");
    // g1->GetYaxis()->SetTitle("");
    // g1->GetXaxis()->SetTitleSize(0.05);
    // g1->GetYaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetLabelSize(0.08);
    g1->GetYaxis()->SetLabelSize(0.08);

    // Set primary axis dimension to 5
    g1->GetXaxis()->SetNdivisions(5);

    // Create a legend
    // Remove legend edge
	gStyle->SetLegendBorderSize(0);
	TLegend *legend_theory = new TLegend(0.65, 0.15, 0.9, 0.3);
    legend_theory->AddEntry(g1, "AMD #pi linear chain", "pl");
    legend_theory->AddEntry(g2, "AMD #sigma linear chain", "pl");
    legend_theory->SetTextSize(0.04);
    legend_theory->Draw();


	TLegend *legend_observed = new TLegend(0.15, 0.7, 0.32, 0.85);
    legend_observed->AddEntry(g4, "this work", "p");
    legend_observed->AddEntry(g3, "previous work", "p");
    legend_observed->SetTextSize(0.04);
    legend_observed->Draw();

}
