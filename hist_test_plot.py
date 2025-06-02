import ROOT
import sys
from array import array

file = ROOT.TFile("rootoutputs/output.root", "READ")

histograms = [
    "h_muon_pt", "h_muon_eta", "h_muon_phi", "h_muon_energy", 
    "h_muon_mass", "h_muon_leading", "h_muon_subleading", 
    "h_Z_pt", "h_Z_eta", "h_Z_phi", "h_Z_energy", "h_Z_mass", "h_Z_mass_eq", "h_Z_mass_fine"
]

ROOT.gStyle.SetOptStat(0)

for hist_name in histograms:
    hist = file.Get(hist_name)

    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlue)
    canvas = ROOT.TCanvas(f"c_{hist_name}", f"{hist_name} Canvas", 800, 800)
    hist.SetYTitle("Number of events")
    canvas.SetLogy()
    hist.SetTitle("")
    if hist_name != "h_Z_mass_fine":
        hist.SetMinimum(1)

    titles = {
        "h_muon_pt": ("p_{T#mu} (GeV)"),
        "h_muon_eta": ("#eta_{#mu}"),
        "h_muon_phi": ("#phi_{#mu}"),
        "h_muon_energy": ("E_{#mu} (GeV)"),
        "h_muon_mass": ("M_{#mu} (GeV)"),
        "h_muon_leading": ("p_{T#mu} (GeV)"),
        "h_muon_subleading": ("p_{T#mu} (GeV)"),
        "h_Z_pt": ("p_{T#mu#mu} (GeV)"),
        "h_Z_eta": ("#eta_{#mu#mu}"),
        "h_Z_phi": ("#phi_{#mu#mu}"),
        "h_Z_energy": ("E_{#mu#mu} (GeV)"),
        "h_Z_mass": ("M_{#mu#mu} (GeV)"),
        "h_Z_mass_eq": ("M_{#mu#mu} (GeV)"),
        "h_Z_mass_fine": ("M_{#mu#mu} (GeV)"),
    }

    xtitle = titles.get(hist_name, (""))
    hist.SetXTitle(xtitle)

    if hist_name == "h_Z_mass":
        canvas.SetLogx()

    elif hist_name == "h_Z_mass_eq":
        canvas.SetLogx()

    elif hist_name == "h_Z_mass_fine":
        combined_func = ROOT.TF1("combined_func", "[0]/((x*x - [1]*[1])^2 + ([1]*[2])^2) + [3] + [4]*x + [5]*x*x", 80, 100)
        combined_func.SetParameters(1e10, 91.2, 2.5, 1e3, -1.0, 0.1)
        combined_func.SetParNames("Norm", "M (Mass)", "Gamma (Width)", "D", "E", "F")
        hist.GetXaxis().SetRangeUser(80, 100)
        fit_result = hist.Fit(combined_func, "S, R", "", 85, 95)
        combined_func.Draw("SAME")
        canvas.Modified()
        canvas.Update()
        
    
    hist.Draw()

    if hist_name != "h_Z_mass_fine":
        legend = ROOT.TLegend(0.75, 0.85, 0.9, 0.9)
        legend.AddEntry(hist, "Data", "l")
    if hist_name == "h_Z_mass_fine":
        legend = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
        legend.AddEntry(hist, "Data", "l")
        legend.AddEntry(combined_func, "Approximation", "l")
    legend.Draw()

    canvas.SaveAs(f"hist_real_data_hist/{hist_name}.png")
    canvas.Close()

file.Close()

