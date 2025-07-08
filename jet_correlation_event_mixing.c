// jet_correlation_event_mixing.C
// This ROOT macro conceptually demonstrates jet-like correlation analysis
// using simulated data and introduces the idea of event mixing for background subtraction.
// It focuses on azimuthal correlations (Delta Phi) between a "trigger" particle (e.g., a high-pT hadron)
// and "associated" particles.
//
// To run this macro in ROOT:
// .L jet_correlation_event_mixing.C
// jet_correlation_event_mixing();

#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <TLegend.h> // Include TLegend for the legend

// --- Particle Struct (Simpler for this example) ---
// Represents a generic particle with its relevant kinematic properties.
struct Track
{
    double pt;  // Transverse momentum
    double phi; // Azimuthal angle
    double eta; // Pseudorapidity
    int id;     // Simple ID (e.g., 0 for associated, 1 for trigger)
};

// --- Main Analysis Function ---
void jet_correlation_event_mixing()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(1111);

    // 1. Define Histograms
    // Signal + Background (Same Event)
    TH1F *hDeltaPhi_SE = new TH1F("hDeltaPhi_SE", "Same-Event #Delta#phi (Trigger-Associated);#Delta#phi (rad);Counts", 100, -TMath::Pi() / 2, 3 * TMath::Pi() / 2);
    TH2F *hDeltaPhiEta_SE = new TH2F("hDeltaPhiEta_SE", "Same-Event #Delta#phi vs #Delta#eta; #Delta#phi (rad); #Delta#eta", 100, -TMath::Pi() / 2, 3 * TMath::Pi() / 2, 100, -2.0, 2.0);

    // Background (Mixed Event) - will be scaled later
    TH1F *hDeltaPhi_ME = new TH1F("hDeltaPhi_ME", "Mixed-Event #Delta#phi (Trigger-Associated);#Delta#phi (rad);Counts", 100, -TMath::Pi() / 2, 3 * TMath::Pi() / 2);
    TH2F *hDeltaPhiEta_ME = new TH2F("hDeltaPhiEta_ME", "Mixed-Event #Delta#phi vs #Delta#eta; #Delta#phi (rad); #Delta#eta", 100, -TMath::Pi() / 2, 3 * TMath::Pi() / 2, 100, -2.0, 2.0);

    // Corrected Correlation Function (Signal)
    TH1F *hDeltaPhi_Corr = new TH1F("hDeltaPhi_Corr", "Corrected #Delta#phi (Jet-like Correlation);#Delta#phi (rad);Counts", 100, -TMath::Pi() / 2, 3 * TMath::Pi() / 2);

    // 2. Simulation Parameters
    TRandom3 *rand = new TRandom3();
    int numEvents = 50000;
    int maxParticlesPerEvent = 50; // Average multiplicity
    double triggerPtMin = 3.0;     // Min pT for trigger particle
    double associatedPtMin = 0.5;  // Min pT for associated particle
    double maxEtaRange = 0.8;      // Pseudorapidity acceptance (corrected variable name)

    // Event buffer for event mixing (stores events to mix with)
    // In a real analysis, you'd store full event information (multiplicity, event plane, etc.)
    // For this conceptual example, we just store a vector of tracks.
    std::vector<std::vector<Track>> eventBuffer;
    int eventBufferDepth = 10; // Number of previous events to mix with

    // 3. Event Loop
    for (int i = 0; i < numEvents; ++i)
    {
        std::vector<Track> currentEventTracks;
        std::vector<Track> triggerParticles;

        // Generate particles for the current event
        int numParticles = rand->Poisson(maxParticlesPerEvent); // Poisson distribution for multiplicity
        for (int j = 0; j < numParticles; ++j)
        {
            Track t;
            t.pt = rand->Uniform(0.2, 10.0); // pT from 0.2 to 10 GeV/c
            t.phi = rand->Uniform(0, 2 * TMath::Pi());
            t.eta = rand->Uniform(-maxEtaRange, maxEtaRange); // Corrected variable name
            t.id = 0;                                         // Default to associated particle

            if (t.pt > triggerPtMin)
            {
                t.id = 1; // Mark as trigger candidate
                triggerParticles.push_back(t);
            }
            currentEventTracks.push_back(t);
        }

        // --- Same-Event Analysis (Signal + Background) ---
        // Loop over trigger particles in the current event
        for (const auto &trigger : triggerParticles)
        {
            // Loop over all associated particles in the same event
            for (const auto &associated : currentEventTracks)
            {
                if (&trigger == &associated)
                    continue; // Skip self-correlation

                if (associated.pt > associatedPtMin && TMath::Abs(associated.eta) < maxEtaRange)
                { // Corrected variable name
                    double deltaPhi = trigger.phi - associated.phi;
                    double deltaEta = trigger.eta - associated.eta;

                    // Normalize DeltaPhi to be within [-pi/2, 3pi/2]
                    while (deltaPhi < -TMath::Pi() / 2)
                        deltaPhi += 2 * TMath::Pi();
                    while (deltaPhi >= 3 * TMath::Pi() / 2)
                        deltaPhi -= 2 * TMath::Pi();

                    hDeltaPhi_SE->Fill(deltaPhi);
                    hDeltaPhiEta_SE->Fill(deltaPhi, deltaEta);
                }
            }
        }

        // --- Mixed-Event Analysis (Background Estimation) ---
        // Mix current event triggers with associated particles from buffered events
        for (const auto &trigger : triggerParticles)
        {
            for (const auto &mixedEvent : eventBuffer)
            {
                for (const auto &associated : mixedEvent)
                {
                    if (associated.pt > associatedPtMin && TMath::Abs(associated.eta) < maxEtaRange)
                    { // Corrected variable name
                        double deltaPhi = trigger.phi - associated.phi;
                        double deltaEta = trigger.eta - associated.eta;

                        // Normalize DeltaPhi
                        while (deltaPhi < -TMath::Pi() / 2)
                            deltaPhi += 2 * TMath::Pi();
                        while (deltaPhi >= 3 * TMath::Pi() / 2)
                            deltaPhi -= 2 * TMath::Pi();

                        hDeltaPhi_ME->Fill(deltaPhi);
                        hDeltaPhiEta_ME->Fill(deltaPhi, deltaEta);
                    }
                }
            }
        }

        // Update event buffer
        eventBuffer.push_back(currentEventTracks);
        if (eventBuffer.size() > eventBufferDepth)
        {
            eventBuffer.erase(eventBuffer.begin()); // Remove oldest event
        }
    }

    // 4. Normalize and Subtract Background
    // The mixed-event background needs to be normalized to the same-event yield.
    // This is often done by scaling the mixed-event histogram by the ratio of
    // the total number of trigger-associated pairs in same-event vs mixed-event.
    // For simplicity here, we'll just normalize by the integral.
    if (hDeltaPhi_ME->Integral() > 0)
    {
        double scaleFactor = hDeltaPhi_SE->Integral() / hDeltaPhi_ME->Integral();
        hDeltaPhi_ME->Scale(scaleFactor);
    }
    else
    {
        std::cerr << "Warning: Mixed-event histogram has zero integral. Cannot normalize." << std::endl;
    }

    // Subtract background to get the corrected signal
    hDeltaPhi_Corr->Add(hDeltaPhi_SE);
    hDeltaPhi_Corr->Add(hDeltaPhi_ME, -1); // Subtract mixed-event background

    // 5. Draw Results
    TCanvas *c1 = new TCanvas("c1", "Delta Phi Correlations", 1200, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    hDeltaPhi_SE->SetTitle("Same-Event #Delta#phi");
    hDeltaPhi_SE->SetLineColor(kBlue);
    hDeltaPhi_SE->Draw();
    hDeltaPhi_ME->SetLineColor(kRed);
    hDeltaPhi_ME->Draw("SAME");
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.85);
    leg1->AddEntry(hDeltaPhi_SE, "Same Event (Signal+Bkg)", "l");
    leg1->AddEntry(hDeltaPhi_ME, "Mixed Event (Bkg)", "l");
    leg1->Draw();

    c1->cd(2);
    hDeltaPhi_Corr->SetTitle("Corrected #Delta#phi (Jet-like Signal)");
    hDeltaPhi_Corr->SetLineColor(kGreen + 2);
    hDeltaPhi_Corr->Draw();
    c1->SaveAs("jet_delta_phi_correlation.png");

    TCanvas *c2 = new TCanvas("c2", "Delta Phi vs Delta Eta (Same Event)", 800, 600);
    hDeltaPhiEta_SE->Draw("COLZ");
    c2->SaveAs("jet_delta_phi_eta_SE.png");

    TCanvas *c3 = new TCanvas("c3", "Delta Phi vs Delta Eta (Mixed Event)", 800, 600);
    hDeltaPhiEta_ME->Draw("COLZ");
    c3->SaveAs("jet_delta_phi_eta_ME.png");

    // Clean up
    // delete hDeltaPhi_SE; delete hDeltaPhiEta_SE;
    // delete hDeltaPhi_ME; delete hDeltaPhiEta_ME;
    // delete hDeltaPhi_Corr;
    // delete rand;
    // delete c1; delete c2; delete c3;
}
