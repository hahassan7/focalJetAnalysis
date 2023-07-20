#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

void chatGPT() {
    const int nBins = 10;
    const double xMin = 0.0;
    const double xMax = 100.0;

    // Arrays to hold histograms for different types and R values
    TH1F* hist_pT_R[3]; // hist_pT_R[0], hist_pT_R[1], hist_pT_R[2]
    TH1F* hist_E_R[3]; // hist_E_R[0], hist_E_R[1], hist_E_R[2]
    TH1F* hist_con_R[3]; // hist_con_R[0], hist_con_R[1], hist_con_R[2]

    TH1F* hist_pT_R_e[3][3]; // hist_pT_R_e[0][0], hist_pT_R_e[0][1], ..., hist_pT_R_e[2][2]
    TH1F* hist_E_R_e[3][3]; // hist_E_R_e[0][0], hist_E_R_e[0][1], ..., hist_E_R_e[2][2]
    TH1F* hist_con_R_e[3][3]; // hist_con_R_e[0][0], hist_con_R_e[0][1], ..., hist_con_R_e[2][2]

    // Read histograms from the ROOT file
    TFile* file = TFile::Open("path_to_your_root_file.root", "READ");

    for (int R = 0; R < 3; ++R) {
        // Read hist_pT_R[3] histograms from the file using Clone method
        hist_pT_R[R] = dynamic_cast<TH1F*>(file->Get(Form("histo_pT_R%d", R))->Clone());
        hist_pT_R[R]->SetTitle(Form("pT Distribution for R = %d; pT; Entries", R));

        // Read hist_E_R[3] histograms from the file using Clone method
        hist_E_R[R] = dynamic_cast<TH1F*>(file->Get(Form("histo_E_R%d", R))->Clone());
        hist_E_R[R]->SetTitle(Form("E Distribution for R = %d; E; Entries", R));

        // Read hist_con_R[3] histograms from the file using Clone method
        hist_con_R[R] = dynamic_cast<TH1F*>(file->Get(Form("histo_con_R%d", R))->Clone());
        hist_con_R[R]->SetTitle(Form("Constituent Distribution for R = %d; Constituent; Entries", R));

        for (int e = 0; e < 3; ++e) {
            // Read hist_pT_R_e[3][3] histograms from the file using Clone method
            hist_pT_R_e[R][e] = dynamic_cast<TH1F*>(file->Get(Form("histo_pT_R%d_e%d", R, e))->Clone());
            hist_pT_R_e[R][e]->SetTitle(Form("pT Distribution for R = %d, e = %d; pT; Entries", R, e));

            // Read hist_E_R_e[3][3] histograms from the file using Clone method
            hist_E_R_e[R][e] = dynamic_cast<TH1F*>(file->Get(Form("histo_E_R%d_e%d", R, e))->Clone());
            hist_E_R_e[R][e]->SetTitle(Form("E Distribution for R = %d, e = %d; E; Entries", R, e));

            // Read hist_con_R_e[3][3] histograms from the file using Clone method
            hist_con_R_e[R][e] = dynamic_cast<TH1F*>(file->Get(Form("histo_con_R%d_e%d", R, e))->Clone());
            hist_con_R_e[R][e]->SetTitle(Form("Constituent Distribution for R = %d, e = %d; Constituent; Entries", R, e));
        }
    }

    // Close the ROOT file
    file->Close();

    // Draw histograms for hist_pT_R[3]
    TCanvas* canvas_pT_R = new TCanvas("canvas_pT_R", "Histogram Canvas (pT_R)", 800, 600);
    canvas_pT_R->SetWindowSize(800, 600);
    canvas_pT_R->SetCanvasSize(800, 600);

    canvas_pT_R->cd();
    hist_pT_R[0]->Draw();
    for (int R = 1; R < 3; ++R) {
        hist_pT_R[R]->Draw("SAME");
    }

    // Draw histograms for hist_E_R[3]
    TCanvas* canvas_E_R = new TCanvas("canvas_E_R", "Histogram Canvas (E_R)", 800, 600);
    canvas_E_R->SetWindowSize(800, 600);
    canvas_E_R->SetCanvasSize(800, 600);

    canvas_E_R->cd();
    hist_E_R[0]->Draw();
    for (int R = 1; R < 3; ++R) {
        hist_E_R[R]->Draw("SAME");
    }

    // Draw histograms for hist_con_R[3]
    TCanvas* canvas_con_R = new TCanvas("canvas_con_R", "Histogram Canvas (con_R)", 800, 600);
    canvas_con_R->SetWindowSize(800, 600);
    canvas_con_R->SetCanvasSize(800, 600);

    canvas_con_R->cd();
    hist_con_R[0]->Draw();
    for (int R = 1; R < 3; ++R) {
        hist_con_R[R]->Draw("SAME");
    }

    // Draw the hist_pT_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_pT_R_e = new TCanvas(Form("canvas_pT_R%d_e", R), Form("Histogram Canvas (pT_R = %d, e)", R), 800, 600);
        canvas_pT_R_e->SetWindowSize(800, 600);
        canvas_pT_R_e->SetCanvasSize(800, 600);

        canvas_pT_R_e->cd();
        hist_pT_R_e[R][0]->Draw();
        for (int e = 1; e < 3; ++e) {
            hist_pT_R_e[R][e]->Draw("SAME");
        }
    }

    // Draw the hist_E_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_E_R_e = new TCanvas(Form("canvas_E_R%d_e", R), Form("Histogram Canvas (E_R = %d, e)", R), 800, 600);
        canvas_E_R_e->SetWindowSize(800, 600);
        canvas_E_R_e->SetCanvasSize(800, 600);

        canvas_E_R_e->cd();
        hist_E_R_e[R][0]->Draw();
        for (int e = 1; e < 3; ++e) {
            hist_E_R_e[R][e]->Draw("SAME");
        }
    }

    // Draw the hist_con_R_e histograms for each R value in separate canvases and one plot each
    for (int R = 0; R < 3; ++R) {
        TCanvas* canvas_con_R_e = new TCanvas(Form("canvas_con_R%d_e", R), Form("Histogram Canvas (con_R = %d, e)", R), 800, 600);
        canvas_con_R_e->SetWindowSize(800, 600);
        canvas_con_R_e->SetCanvasSize(800, 600);

        canvas_con_R_e->cd();
        hist_con_R_e[R][0]->Draw();
        for (int e = 1; e < 3; ++e) {
            hist_con_R_e[R][e]->Draw("SAME");
        }
    }
}
