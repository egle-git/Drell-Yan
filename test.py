import ROOT
import sys
import os

lumi = 16494
xsecs = {
    "sim1": 6422,
    "sim2": 20480,
    "tt": 687.1,
    # "ttnew": 756.2,
    "ww": 75.88,
    "wz": 27.56,
    "zz": 12.08,
    "twtop": 32.45,
    "twantitop": 32.51,
    "tchantop": 119.7,
    "tchanantitop": 71.74,
}
# weightsums = {
#     "sim1": 4.96939e+07,
#     "sim2": 3.69428e+07,
#     "tt": 4.32772e+07,
#     # "ttnew": ?,
#     "ww": 1.5821e+07,
#     "wz": 7.584e+06,
#     "zz": 1.151e+06,
#     "twtop": 3.36824e+06,
#     "twantitop": 3.65433e+06,
#     "tchantop": 5.54612e+07,
#     "tchanantitop": 2.92337e+07,
# }

weightsums = {
    "sim1": 544907,
    "sim2":494561,
    "tt": 373938,
    # "ttnew": 169658,
    "ww": 788000,
    "wz": 331000,
    "zz": 483000,
    "twtop": 333250,
    "twantitop": 673726,
    "tchantop": 539780,
    "tchanantitop": 835376,
}

varbins = [40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000]
nbins = 36

filesandtrees = {
    "real": {"file": "rootoutputs/treeout.root", "tree": "Events"},
    "sim1": {"file": "rootoutputs/treeout1.root", "tree": "Muon"},
    "sim2": {"file": "rootoutputs/treeout2.root", "tree": "Muon"},
    "sim1_DYtau": {"file": "rootoutputs/treeout1.root", "tree": "Tau"},
    "sim2_DYtau": {"file": "rootoutputs/treeout2.root", "tree": "Tau"},
    "tt": {"file": "rootoutputs/treeouttt.root", "tree": "Events"},
    # "ttnew": {"file": "rootoutputs/treeoutttnew.root", "tree": "Events"},
    # "twtop": {"file": "rootoutputs/treeouttwtop.root", "tree": "Events"},
    "twantitop": {"file": "rootoutputs/treeouttwantitop.root", "tree": "Events"},
    "tchantop": {"file": "rootoutputs/treeouttchantop.root", "tree": "Events"},
    "tchanantitop": {"file": "rootoutputs/treeouttchanantitop.root", "tree": "Events"},
    # "ww": {"file": "rootoutputs/treeoutww.root", "tree": "Events"},
    # "wz": {"file": "rootoutputs/treeoutwz.root", "tree": "Events"},
    # "zz": {"file": "rootoutputs/treeoutzz.root", "tree": "Events"},
}

histograms = {
    "muon_pt": {"branches": {"muon_pt1", "muon_pt2"}, "Bins": (100, 0, 150)},
    "muon_eta": {"branches": {"muon_eta1", "muon_eta2"}, "Bins": (100, -2.5, 2.5)},
    "muon_phi": {"branches": {"muon_phi1", "muon_phi2"}, "Bins": (100, -3.14, 3.14)},
    "muon_energy": {"branches": {"muon_energy1", "muon_energy2"}, "Bins": (100, 0, 150)},
    "muon_mass": {"branches": {"muon_mass1", "muon_mass2"}, "Bins": (100, 0.1055, 0.1059)},
    "muon_p": {"branches": {"muon_p1", "muon_p2"}, "Bins": (100, 0, 150)},
    "muon_px": {"branches": {"muon_px1", "muon_px2"}, "Bins": (100, 0, 150)},
    "muon_py": {"branches": {"muon_py1", "muon_py2"}, "Bins": (100, 0, 150)},
    "muon_pz": {"branches": {"muon_pz1", "muon_pz2"}, "Bins": (100, 0, 150)},
    "muon_charge": {"branches": {"muon_charge1", "muon_charge2"}, "Bins": (3, -1.5, 1.5)},
    "muon_leading": {"branches": {"muon_pt1"}, "Bins": (100, 0, 150)},
    "muon_subleading": {"branches": {"muon_pt2"}, "Bins": (100, 0, 150)},

    "Z_pt": {"branches": {"Z_pt"}, "Bins": (100, 0, 300)},
    "Z_eta": {"branches": {"Z_eta"}, "Bins": (100, -2.5, 2.5)},
    "Z_phi": {"branches": {"Z_phi"}, "Bins": (100, -3.14, 3.14)},
    "Z_energy": {"branches": {"Z_energy"}, "Bins": (100, 0, 500)},
    "Z_mass": {"branches": {"Z_mass"}, "Bins": (1000, 30, 1000)},
    "Z_mass_eq": {"branches": {"Z_mass"}, "Bins": (1000, 30, 1000)},
    "Z_mass_fine": {"branches": {"Z_mass"}, "Bins": (150, 70, 110)},
    "Z_px": {"branches": {"Z_px"}, "Bins": (100, 0, 300)},
    "Z_py": {"branches": {"Z_py"}, "Bins": (100, 0, 300)},
    "Z_pz": {"branches": {"Z_pz"}, "Bins": (100, 0, 300)},
}

open_files = {}
for dataname, datainfo in filesandtrees.items():
    file = ROOT.TFile.Open(datainfo["file"])
    tree = file.Get(datainfo["tree"])
    open_files[dataname] = {"file": file, "tree": tree}

def fill_hist(tree, branch, bins, xmin, xmax, weight_branch=None):
    hist = ROOT.TH1F("hist", "", bins, xmin, xmax)
    draw_cmd = f"{branch}>>hist"
    weight_opt = weight_branch if weight_branch else "1"
    tree.Draw(draw_cmd, weight_opt, "goff")
    return hist

for histname, histinfo in histograms.items():
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("canvas_" + histname, "Real vs Simulated", 800, 800)
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

    branches = list(histinfo["branches"])
    bins, xmin, xmax = histinfo["Bins"]

    def fill_multi_branch_hist(tree, branches, bins, xmin, xmax, weight_branch=None):
        hist = ROOT.TH1F("hist", "", bins, xmin, xmax)
        for br in branches:
            draw_cmd = f"{br}>>+hist"
            weight_opt = weight_branch if weight_branch else "1"
            tree.Draw(draw_cmd, weight_opt, "goff")
        return hist

    hist_real = fill_multi_branch_hist(open_files["real"]["tree"], branches, bins, xmin, xmax)
    hist_real.SetMarkerColor(ROOT.kBlack)
    hist_real.SetLineColor(ROOT.kBlack)
    hist_real.SetMarkerSize(0.8)
    hist_real.SetMarkerStyle(20)

    def scale_hist(name, branch):
        t = open_files[name]["tree"]
        weight = xsecs.get(name.replace("_DYtau", ""), 1) * lumi / weightsums.get(name.replace("_DYtau", ""), 1)
        weight_expr = f"muon_norm_weight*{weight}"
        return fill_hist(t, branch, bins, xmin, xmax, weight_branch=weight_expr)

    branch = branches[0]

    tree.Draw("muon_norm_weight>>hw(100, 0, 5)", "", "goff")
    hw = ROOT.gDirectory.Get("hw")
    print("Mean:", hw.GetMean(), "Entries:", hw.GetEntries())


    sim1 = scale_hist("sim1", branch)
    sim2 = scale_hist("sim2", branch) if "sim2" in open_files else ROOT.TH1F("empty2", "", bins, xmin, xmax)
    sim1_DYtau = scale_hist("sim1_DYtau", branch)
    sim2_DYtau = scale_hist("sim2_DYtau", branch) if "sim2_DYtau" in open_files else ROOT.TH1F("empty2tau", "", bins, xmin, xmax)

    DY = sim1.Clone("DY")
    DY.Add(sim2)
    DYtau = sim1_DYtau.Clone("DYtau")
    DYtau.Add(sim2_DYtau)

    tt = scale_hist("tt", branch)

    st = ROOT.TH1F("st", "", bins, xmin, xmax)
    added_any = False
    for st_name in ["twtop", "twantitop", "tchantop", "tchanantitop"]:
        if st_name in open_files:
            if not added_any:
                st = scale_hist(st_name, branch)
                added_any = True
            else:
                st.Add(scale_hist(st_name, branch))

    ew = ROOT.TH1F("ew", "", bins, xmin, xmax)
    added_any = False
    for ew_name in ["ww", "wz", "zz"]:
        if ew_name in open_files:
            if not added_any:
                ew = scale_hist(ew_name, branch)
                added_any = True
            else:
                ew.Add(scale_hist(ew_name, branch))

    for h, color in zip([DY, DYtau, tt, st, ew],
                        [ROOT.kRed + 2, ROOT.kOrange + 2, ROOT.kBlue + 2, ROOT.kMagenta + 2, ROOT.kGreen + 2]):
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetFillColor(color)
        h.SetFillStyle(1001)

    stack = ROOT.THStack("stack", "")
    for h in [ew, st, tt, DYtau, DY]:
        stack.Add(h)
    stack.SetMinimum(1)

    pad1.SetLogy()
    if histname in ["Z_mass", "Z_mass_eq"]:
        pad1.SetLogx()

    # hist_real.Draw("PE")
    # stack.Draw("HIST SAME")
    # hist_real.Draw("PE SAME")
    # hist_real.Draw("AXIS SAME")

    stack.Draw("HIST")
    hist_real.Draw("PE SAME")
    hist_real.Draw("AXIS SAME")

    legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    legend.AddEntry(hist_real, "Data", "pl")
    legend.AddEntry(DY, "DY #rightarrow #mu#mu", "f")
    legend.AddEntry(DYtau, "DY #rightarrow #tau#tau", "f")
    legend.AddEntry(tt, "t#bar{t}", "f")
    legend.AddEntry(st, "Single top", "f")
    legend.AddEntry(ew, "EW", "f")
    legend.Draw()

    pad2.cd()
    total_mc = ew.Clone("total_mc")
    total_mc.Add(st)
    total_mc.Add(tt)
    total_mc.Add(DYtau)
    total_mc.Add(DY)

    ratio = hist_real.Clone("ratio")
    ratio.Divide(total_mc)
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)

    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.65)

    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(20)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.5)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(20)
    ratio.GetXaxis().SetTitleSize(25)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(1)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(20)
    ratio.SetTitle("")
    ratio.GetXaxis().SetTitle(branch)

    if histname in ["Z_mass", "Z_mass_eq"]:
        pad2.SetLogx()
    if histname == "Z_mass_fine":
        ratio.GetXaxis().SetRangeUser(80, 100)

    ratio.Draw("PE")
    line = ROOT.TLine(xmin, 1, xmax, 1)
    line.SetLineColor(ROOT.kBlack)
    line.Draw()

    canvas.SaveAs(f"hist_tree_both_data_hist/{histname}.png")