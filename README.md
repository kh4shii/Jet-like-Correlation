# Jet-like Correlation Analysis with Event Mixing

ROOT/C++ implementation of two-particle angular correlations using event mixing for background estimation. Developed as preparation for charm-hadron correlation studies in the **BEACON project**.

## 🎯 Physics Goal
Study jet-like structures via Δη–Δφ correlation functions:
- Same-event pairs: signal + combinatorial background
- Mixed-event pairs: background-only estimate (multiplicity-matched)
- Correlation function: C(Δη,Δφ) = Same / Mixed

## 🔧 How to Run
```bash
root -l -b -q jet_correlation_event_mixing.c
