import ROOT
import sys
from array import array


lumi = 16494
xsecs = {
    "sim1": 6422,
    "sim2": 20480,
    "tt": 687.1,
    "ww": 75.88,
    "wz": 27.56,
    "zz": 12.08,
    "twtop": 32.45,
    "twantitop": 32.51,
    "tchantop": 119.7,
    "tchanantitop": 71.74,
}
weightsums = {
    "sim1": 4.96939e+07,
    "sim2": 3.69428e+07,
    "tt": 4.32772e+07,
    "ww": 1.5821e+07,
    "wz": 7.584e+06,
    "zz": 1.151e+06,
    "twtop": 3.36824e+06,
    "twantitop": 3.65433e+06,
    "tchantop": 5.54612e+07,
    "tchanantitop": 2.92337e+07,
}


file_real = ROOT.TFile.Open("rootoutput/output.root", "READ")
file_sim1 = ROOT.TFile.Open("rootoutput/simoutput1.root", "READ")
file_sim2 = ROOT.TFile.Open("rootoutput/simoutput2.root", "READ")
file_tt = ROOT.TFile.Open("rootoutput/simoutputtt.root", "READ")
file_twtop = ROOT.TFile.Open("rootoutput/simoutputtwtop.root", "READ")
file_twantitop = ROOT.TFile.Open("rootoutput/simoutputtwantitop.root", "READ")
file_tchantop = ROOT.TFile.Open("rootoutput/simoutputtchantop.root", "READ")
file_tchanantitop = ROOT.TFile.Open("rootoutput/simoutputtchanantitop.root", "READ")
file_ww = ROOT.TFile.Open("rootoutput/simoutputww.root", "READ")
file_wz = ROOT.TFile.Open("rootoutput/simoutputwz.root", "READ")
file_zz = ROOT.TFile.Open("rootoutput/simoutputzz.root", "READ")

hist_pairs = {
    "h_muon_pt": "simh_muon_pt",
    "h_muon_eta": "simh_muon_eta",
    "h_muon_phi": "simh_muon_phi",
    "h_muon_energy": "simh_muon_energy",
    "h_muon_mass": "simh_muon_mass",
    "h_muon_leading": "simh_muon_leading",
    "h_muon_subleading": "simh_muon_subleading",
    "h_Z_pt": "simh_Z_pt",
    "h_Z_eta": "simh_Z_eta",
    "h_Z_phi": "simh_Z_phi",
    "h_Z_energy": "simh_Z_energy",
    "h_Z_mass": "simh_Z_mass",
    "h_Z_mass_eq": "simh_Z_mass_eq",
    "h_Z_mass_fine": "simh_Z_mass_fine",
}

for real_name, sim_name in hist_pairs.items():
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("canvas_"+real_name, "Real vs Simulated", 800, 800)

    pad1 = ROOT.TPad("pad1", "Main Plot", 0, 0.3, 1, 1)
    pad2 = ROOT.TPad("pad2", "Ratio Plot", 0, 0, 1, 0.3)
    pad1.SetBottomMargin(0)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad1.SetGrid()
    pad2.SetGrid()
    pad1.Draw()
    pad2.Draw()
    pad1.cd()


    hist_real = file_real.Get(real_name)
    

    hist_sim1 = file_sim1.Get(sim_name)
    hist_sim1.Scale(xsecs["sim1"]*lumi/weightsums["sim1"])
    hist_sim2 = file_sim2.Get(sim_name)
    hist_sim2.Scale(xsecs["sim2"]*lumi/weightsums["sim2"])
    
    sim_name_DYtau = sim_name.replace("simh_", "simh_DYtau_")
    hist_sim1_DYtau = file_sim1.Get(sim_name_DYtau)
    hist_sim1_DYtau.Scale(xsecs["sim1"]*lumi/weightsums["sim1"])
    hist_sim2_DYtau = file_sim2.Get(sim_name_DYtau)
    hist_sim2_DYtau.Scale(xsecs["sim2"]*lumi/weightsums["sim2"])
    
    hist_sim_combined = hist_sim1.Clone("hist_sim_combined")
    hist_sim_combined.Add(hist_sim2)
    
    hist_sim_DYtau_combined = hist_sim1_DYtau.Clone("hist_sim_DYtau_combined")
    hist_sim_DYtau_combined.Add(hist_sim2_DYtau)


    hist_tt = file_tt.Get(sim_name)
    hist_tt.Scale(xsecs["tt"]*lumi/weightsums["tt"]/2)


    hist_twtop = file_twtop.Get(sim_name)
    hist_twtop.Scale(xsecs["twtop"]*lumi/weightsums["twtop"])
    hist_twantitop = file_twantitop.Get(sim_name)
    hist_twantitop.Scale(xsecs["twantitop"]*lumi/weightsums["twantitop"])
    hist_tchantop = file_tchantop.Get(sim_name)
    hist_tchantop.Scale(xsecs["tchantop"]*lumi/weightsums["tchantop"])
    hist_tchanantitop = file_tchanantitop.Get(sim_name)
    hist_tchanantitop.Scale(xsecs["tchanantitop"]*lumi/weightsums["tchanantitop"])

    hist_st_combined = hist_twtop.Clone("hist_st_combined")
    hist_st_combined.Add(hist_twantitop)
    hist_st_combined.Add(hist_tchantop)
    hist_st_combined.Add(hist_tchanantitop)


    hist_ww = file_ww.Get(sim_name)
    hist_ww.Scale(xsecs["ww"]*lumi/weightsums["ww"])
    hist_wz = file_wz.Get(sim_name)
    hist_wz.Scale(xsecs["wz"]*lumi/weightsums["wz"])
    hist_zz = file_zz.Get(sim_name)
    hist_zz.Scale(xsecs["zz"]*lumi/weightsums["zz"])

    hist_ew_combined = hist_ww.Clone("hist_ew_combined")
    hist_ew_combined.Add(hist_wz)
    hist_ew_combined.Add(hist_zz)


    pad1.SetLogy()
    hist_real.SetMarkerColor(ROOT.kBlack)
    hist_real.SetLineColor(ROOT.kBlack)
    hist_real.SetMarkerSize(0.8)
    hist_real.SetMarkerStyle(20)

    hist_sim_combined.SetLineColor(ROOT.kRed+2)
    hist_sim_combined.SetLineWidth(2)
    hist_sim_combined.SetFillColor(ROOT.kRed+2)
    hist_sim_combined.SetFillStyle(1001)

    hist_sim_DYtau_combined.SetLineColor(ROOT.kOrange+2)
    hist_sim_DYtau_combined.SetLineWidth(2)
    hist_sim_DYtau_combined.SetFillColor(ROOT.kOrange+2)
    hist_sim_DYtau_combined.SetFillStyle(1001)
    
    hist_tt.SetLineColor(ROOT.kBlue+2)
    hist_tt.SetLineWidth(2)
    hist_tt.SetFillColor(ROOT.kBlue+2)
    hist_tt.SetFillStyle(1001)

    hist_st_combined.SetLineColor(ROOT.kMagenta+2)
    hist_st_combined.SetLineWidth(2)
    hist_st_combined.SetFillColor(ROOT.kMagenta+2)
    hist_st_combined.SetFillStyle(1001)

    hist_ew_combined.SetLineColor(ROOT.kGreen+2)
    hist_ew_combined.SetLineWidth(2)
    hist_ew_combined.SetFillColor(ROOT.kGreen+2)
    hist_ew_combined.SetFillStyle(1001)

    stack = ROOT.THStack("stack", "")
    stack.Add(hist_ew_combined)
    stack.Add(hist_st_combined)
    stack.Add(hist_tt)
    stack.Add(hist_sim_DYtau_combined)
    stack.Add(hist_sim_combined)
    stack.SetMinimum(1)

    if real_name in ["h_Z_mass", "h_Z_mass_eq"]:
        pad1.SetLogx()

    hist_real.Draw("PE")
    stack.Draw("HIST")
    hist_real.Draw("PE SAME")
    hist_real.Draw("SAMEAXIS")
    
    stack.GetYaxis().SetTitle("Number of Events")
    stack.GetYaxis().SetTitleSize(20)
    stack.GetYaxis().SetTitleFont(43)
    stack.GetYaxis().SetTitleOffset(1.5)
    stack.GetYaxis().SetLabelSize(20)
    stack.GetYaxis().SetLabelFont(43)
    if real_name in ["h_Z_mass", "h_Z_mass_eq"]:
        stack.GetXaxis().SetMoreLogLabels()
        stack.GetXaxis().SetNoExponent()
        hist_real.GetXaxis().SetMoreLogLabels()
        hist_real.GetXaxis().SetNoExponent()

    legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    legend.AddEntry(hist_real, "Data", "pl")
    legend.AddEntry(hist_sim_combined, "DY #rightarrow #mu#mu", "f")
    legend.AddEntry(hist_sim_DYtau_combined, "DY #rightarrow #tau#tau", "f")
    legend.AddEntry(hist_tt, "t#bar{t}", "f")
    legend.AddEntry(hist_st_combined, "Single top", "f")
    legend.AddEntry(hist_ew_combined, "EW", "f")
    legend.Draw()

    pad2.cd()
    pad2.Clear()

    hist_mc_total = hist_ew_combined.Clone("hist_mc_total")
    hist_mc_total.Add(hist_tt)
    hist_mc_total.Add(hist_st_combined)
    hist_mc_total.Add(hist_sim_DYtau_combined)
    hist_mc_total.Add(hist_sim_combined)

    hist_ratio = hist_real.Clone("hist_ratio")
    hist_ratio.Divide(hist_mc_total)
    hist_ratio.GetYaxis().SetRangeUser(0.5, 1.5)

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

    xtitle = titles.get(real_name, (""))
    hist_ratio.SetXTitle(xtitle)

    if real_name == "h_Z_mass_eq":
        pad2.SetLogx()
    if real_name == "h_Z_mass":
        pad2.SetLogx()
    if real_name == "h_Z_mass_fine":
        hist_ratio.GetXaxis().SetRangeUser(80, 100)

    hist_ratio.SetLineColor(ROOT.kBlack)
    hist_ratio.SetMarkerColor(ROOT.kBlack)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetMarkerSize(0.65)
    hist_ratio.SetTitle("")

    hist_ratio.GetYaxis().SetTitle("Data / MC")
    hist_ratio.GetYaxis().CenterTitle()
    hist_ratio.GetYaxis().SetNdivisions(505)
    hist_ratio.GetYaxis().SetTitleSize(20)
    hist_ratio.GetYaxis().SetTitleFont(43)
    hist_ratio.GetYaxis().SetTitleOffset(1.5)
    hist_ratio.GetYaxis().SetLabelFont(43)
    hist_ratio.GetYaxis().SetLabelSize(20)
    hist_ratio.GetXaxis().SetTitleSize(25)
    hist_ratio.GetXaxis().SetTitleFont(43)
    hist_ratio.GetXaxis().SetTitleOffset(1)
    hist_ratio.GetXaxis().SetLabelFont(43)
    hist_ratio.GetXaxis().SetLabelSize(20)

    hist_ratio.Draw("PE")

    if real_name == "h_Z_mass":
        line = ROOT.TLine(0, 1, 1000, 1)
    elif real_name == "h_Z_mass_eq":
        line = ROOT.TLine(30, 1, 1000, 1)
    elif real_name == "h_Z_mass_fine":
        line = ROOT.TLine(80, 1, 100, 1)
    else:
        line = ROOT.TLine(hist_ratio.GetXaxis().GetXmin(), 1, hist_ratio.GetXaxis().GetXmax(), 1)

    line.SetLineColor(ROOT.kBlack)
    line.Draw()

    canvas.SaveAs(f"hist_both_data_hist/{real_name}.png")