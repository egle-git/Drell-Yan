import ROOT
import array

file = ROOT.TFile("rootoutputs/treeout.root", "READ")
tree = file.Get("Events")

varbins = [40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000]
nbins = 36

muon_vars = ["pt", "eta", "phi", "energy", "mass", "p", "px", "py", "pz", "charge"]
Z_vars = ["pt", "eta", "phi", "energy", "mass", "px", "py", "pz"]
histograms = {
    "muon_pt": (100, 0, 150),
    "muon_eta": (100, -2.5, 2.5),
    "muon_phi": (100, -3.14, 3.14),
    "muon_energy": (100, 0, 150),
    "muon_mass": (100, 0.1055, 0.1059),
    "muon_p": (100, 0, 150),
    "muon_px": (100, 0, 150),
    "muon_py": (100, 0, 150),
    "muon_pz": (100, 0, 150),
    "muon_charge": (3, -1.5, 1.5),
    "muon_leading": (100, 0, 150),
    "muon_subleading": (100, 0, 150),

    "Z_pt": (100, 0, 300),
    "Z_eta": (100, -2.5, 2.5),
    "Z_phi": (100, -3.14, 3.14),
    "Z_energy": (100, 0, 500),
    "Z_mass": (1000, 30, 1000),
    "Z_mass_eq": (1000, 30, 1000),
    "Z_mass_fine": (150, 70, 110),
    "Z_px": (100, 0, 300),
    "Z_py": (100, 0, 300),
    "Z_pz": (100, 0, 300)
}

# ROOT.gStyle.SetOptStat(0)


##############################################################################

index = []
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    if tree.muon_tight1 and tree.muon_tight2:
        index.append(i)

print(index)
print("\n")
print(len(index))
print("\n")
print(tree.GetEntries())
print("\n")

for hist_name, (bins, xmin, xmax) in histograms.items():
    canvas = ROOT.TCanvas(f"canvas_{hist_name}", hist_name, 800, 800)
    if hist_name == "Z_mass":
        hist = ROOT.TH1F(hist_name, hist_name, nbins, array.array("d", varbins))
        for i in index:
            tree.GetEntry(i)
            hist.Fill(tree.Z_mass)
    else:
        hist = ROOT.TH1F(hist_name, hist_name, bins, xmin, xmax)
        for i in index:
            tree.GetEntry(i)
            if hist_name.startswith("muon_"):
                if hist_name == "muon_leading":
                    hist.Fill(tree.muon_pt1)
                elif hist_name == "muon_subleading":
                    hist.Fill(tree.muon_pt2)
                else:
                    hist.Fill(getattr(tree, f"{hist_name}1"))
                    hist.Fill(getattr(tree, f"{hist_name}2"))

            else:
                if (hist_name == "Z_mass_eq") or (hist_name == "Z_mass_fine"):
                    hist.Fill(tree.Z_mass)
                else:
                    hist.Fill(getattr(tree, f"{hist_name}"))
    
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlue)
    hist.SetYTitle("Number of events")
    canvas.SetLogy()
    hist.SetTitle("")

    if hist_name != "Z_mass_fine":
        hist.SetMinimum(1)

    titles = {
        "muon_pt": ("p_{T#mu} (GeV)"),
        "muon_eta": ("#eta_{#mu}"),
        "muon_phi": ("#phi_{#mu}"),
        "muon_energy": ("E_{#mu} (GeV)"),
        "muon_mass": ("M_{#mu} (GeV)"),
        "muon_p": ("p_{#mu} (GeV)"),
        "muon_px": ("p_{x#mu} (GeV)"),
        "muon_py": ("p_{y#mu} (GeV)"),
        "muon_pz": ("p_{z#mu} (GeV)"),
        "muon_charge": ("charge"),
        "muon_leading": ("p_{T#mu} (GeV)"),
        "muon_subleading": ("p_{T#mu} (GeV)"),
        "Z_pt": ("p_{T#mu#mu} (GeV)"),
        "Z_eta": ("#eta_{#mu#mu}"),
        "Z_phi": ("#phi_{#mu#mu}"),
        "Z_energy": ("E_{#mu#mu} (GeV)"),
        "Z_mass": ("M_{#mu#mu} (GeV)"),
        "Z_mass_eq": ("M_{#mu#mu} (GeV)"),
        "Z_mass_fine": ("M_{#mu#mu} (GeV)"),
        "Z_px": ("p_{x#mu#mu} (GeV)"),
        "Z_py": ("p_{y#mu#mu} (GeV)"),
        "Z_pz": ("p_{z#mu#mu} (GeV)"),
    }

    xtitle = titles.get(hist_name, (""))
    hist.SetXTitle(xtitle)

    if hist_name == "Z_mass":
        canvas.SetLogx()

    elif hist_name == "Z_mass_eq":
        canvas.SetLogx()

    elif hist_name == "Z_mass_fine":
        combined_func = ROOT.TF1("combined_func", "[0]/((x*x - [1]*[1])^2 + ([1]*[2])^2) + [3] + [4]*x + [5]*x*x", 80, 100)
        combined_func.SetParameters(1e10, 91.2, 2.5, 1e3, -1.0, 0.1)
        combined_func.SetParNames("Norm", "M (Mass)", "Gamma (Width)", "D", "E", "F")
        hist.GetXaxis().SetRangeUser(80, 100)
        fit_result = hist.Fit(combined_func, "S, R", "", 85, 95)
        combined_func.Draw("SAME")
        canvas.Modified()
        canvas.Update()
        

    hist.Draw()

    if hist_name != "Z_mass_fine":
        legend = ROOT.TLegend(0.75, 0.85, 0.9, 0.9)
        legend.AddEntry(hist, "Data", "l")
    if hist_name == "Z_mass_fine":
        legend = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
        legend.AddEntry(hist, "Data", "l")
        legend.AddEntry(combined_func, "Approximation", "l")
    # legend.Draw()


    canvas.SaveAs(f"hist_tree_real_data_tightmuon_hist/new_{hist_name}.png")
    canvas.Close()

file.Close()



################################################################################


# for hist_name, (bins, xmin, xmax) in histograms.items():
#     canvas = ROOT.TCanvas(f"canvas_{hist_name}", hist_name, 800, 800)
#     if hist_name == "Z_mass":
#         hist = ROOT.TH1F(hist_name, hist_name, nbins, array.array("d", varbins))
#         draw_expr = f"{hist_name} >> {hist_name}"
#         tree.Draw(draw_expr, "", "goff")
#     else:
#         hist = ROOT.TH1F(hist_name, hist_name, bins, xmin, xmax)

#         if hist_name.startswith("muon_"):
#             if hist_name == "muon_leading":
#                 draw_expr1 = f"muon_pt1 >>+{hist_name}"
#                 tree.Draw(draw_expr1, "", "goff")
#             elif hist_name == "muon_subleading":
#                 draw_expr2 = f"muon_pt2 >>+{hist_name}"
#                 tree.Draw(draw_expr2, "", "goff")
#             else:
#                 base = hist_name.replace("muon_", "")
#                 draw_expr1 = f"muon_{base}1 >>+{hist_name}"
#                 tree.Draw(draw_expr1, "", "goff")
#                 draw_expr2 = f"muon_{base}2 >>+{hist_name}"
#                 tree.Draw(draw_expr2, "", "goff")
        
#         else:
#             if hist_name == "Z_mass_eq":
#                 draw_expr = f"Z_mass >> {hist_name}"
#                 tree.Draw(draw_expr, "", "goff")
#             elif hist_name == "Z_mass_fine":
#                 draw_expr = f"Z_mass >> {hist_name}"
#                 tree.Draw(draw_expr, "", "goff")
#             else:
#                 draw_expr = f"{hist_name} >> {hist_name}"
#                 tree.Draw(draw_expr, "", "goff")
        
#     hist.SetLineWidth(2)
#     hist.SetLineColor(ROOT.kBlue)
#     hist.SetYTitle("Number of events")
#     canvas.SetLogy()
#     hist.SetTitle("")

#     if hist_name != "Z_mass_fine":
#         hist.SetMinimum(1)

#     titles = {
#         "muon_pt": ("p_{T#mu} (GeV)"),
#         "muon_eta": ("#eta_{#mu}"),
#         "muon_phi": ("#phi_{#mu}"),
#         "muon_energy": ("E_{#mu} (GeV)"),
#         "muon_mass": ("M_{#mu} (GeV)"),
#         "muon_p": ("p_{#mu} (GeV)"),
#         "muon_px": ("p_{x#mu} (GeV)"),
#         "muon_py": ("p_{y#mu} (GeV)"),
#         "muon_pz": ("p_{z#mu} (GeV)"),
#         "muon_charge": ("charge"),
#         "muon_leading": ("p_{T#mu} (GeV)"),
#         "muon_subleading": ("p_{T#mu} (GeV)"),
#         "Z_pt": ("p_{T#mu#mu} (GeV)"),
#         "Z_eta": ("#eta_{#mu#mu}"),
#         "Z_phi": ("#phi_{#mu#mu}"),
#         "Z_energy": ("E_{#mu#mu} (GeV)"),
#         "Z_mass": ("M_{#mu#mu} (GeV)"),
#         "Z_mass_eq": ("M_{#mu#mu} (GeV)"),
#         "Z_mass_fine": ("M_{#mu#mu} (GeV)"),
#         "Z_px": ("p_{x#mu#mu} (GeV)"),
#         "Z_py": ("p_{y#mu#mu} (GeV)"),
#         "Z_pz": ("p_{z#mu#mu} (GeV)"),
#     }

#     xtitle = titles.get(hist_name, (""))
#     hist.SetXTitle(xtitle)

#     if hist_name == "Z_mass":
#         canvas.SetLogx()

#     elif hist_name == "Z_mass_eq":
#         canvas.SetLogx()

#     elif hist_name == "Z_mass_fine":
#         combined_func = ROOT.TF1("combined_func", "[0]/((x*x - [1]*[1])^2 + ([1]*[2])^2) + [3] + [4]*x + [5]*x*x", 80, 100)
#         combined_func.SetParameters(1e10, 91.2, 2.5, 1e3, -1.0, 0.1)
#         combined_func.SetParNames("Norm", "M (Mass)", "Gamma (Width)", "D", "E", "F")
#         hist.GetXaxis().SetRangeUser(80, 100)
#         fit_result = hist.Fit(combined_func, "S, R", "", 85, 95)
#         combined_func.Draw("SAME")
#         canvas.Modified()
#         canvas.Update()
        

#     hist.Draw()

#     if hist_name != "Z_mass_fine":
#         legend = ROOT.TLegend(0.75, 0.85, 0.9, 0.9)
#         legend.AddEntry(hist, "Data", "l")
#     if hist_name == "Z_mass_fine":
#         legend = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
#         legend.AddEntry(hist, "Data", "l")
#         legend.AddEntry(combined_func, "Approximation", "l")
#     # legend.Draw()

#     canvas.SaveAs(f"hist_tree_real_data_hist/{hist_name}.png")
#     canvas.Close()

# file.Close()