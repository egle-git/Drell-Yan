import ROOT

file_sim1 = ROOT.TFile.Open("rootoutput/simoutputtrue1.root", "READ")
file_sim2 = ROOT.TFile.Open("rootoutput/simoutputtrue2.root", "READ")


hist_names = {"simh_particle_pt", "simh_particle_eta", "simh_particle_phi",
              "simh_particle_energy", "simh_particle_mass", "simh_particle_leading",
              "simh_particle_subleading", "simh_Z_pt", "simh_Z_eta", "simh_Z_phi",
              "simh_Z_energy", "simh_Z_mass", "simh_Z_mass_fine",
}

for sim_name in hist_names:
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("canvas", "True", 800, 800)

    hist_sim1 = file_sim1.Get(sim_name)
    hist_sim2 = file_sim2.Get(sim_name)
    hist_sim1.Scale(6422*16494/4.96939e+07)
    hist_sim2.Scale(20480*16494/3.69428e+07)

    hist_sim1.SetLineColor(ROOT.kRed)
    hist_sim1.SetLineWidth(2)
    hist_sim2.SetLineColor(ROOT.kBlue)
    hist_sim2.SetLineWidth(2)
    hist_sim2.SetYTitle("Number of events")
    canvas.SetLogy()
    hist_sim2.SetMinimum(1)
    hist_sim2.SetTitle("")

    titles = {
        "simh_particle_pt": ("p_{T#mu} (GeV)"),
        "simh_particle_eta": ("#eta_{#mu}"),
        "simh_particle_phi": ("#phi_{#mu}"),
        "simh_particle_energy": ("E_{#mu} (GeV)"),
        "simh_particle_mass": ("M_{#mu} (GeV)"),
        "simh_particle_leading": ("p_{T#mu} (GeV)"),
        "simh_particle_subleading": ("p_{T#mu} (GeV)"),
        "simh_Z_pt": ("p_{T#mu#mu} (GeV)"),
        "simh_Z_eta": ("#eta_{#mu#mu}"),
        "simh_Z_phi": ("#phi_{#mu#mu}"),
        "simh_Z_energy": ("E_{#mu#mu} (GeV)"),
        "simh_Z_mass": ("M_{#mu#mu} (GeV)"),
        "simh_Z_mass_eq": ("M_{#mu#mu} (GeV)"),
        "simh_Z_mass_fine": ("M_{#mu#mu} (GeV)"),
    }

    xtitle = titles.get(sim_name, (""))
    hist_sim2.SetXTitle(xtitle)

    if sim_name == "simh_Z_mass":
        canvas.SetLogx()
        hist_sim2.GetXaxis().SetMoreLogLabels()
        hist_sim2.GetXaxis().SetNoExponent()

    elif sim_name == "simh_Z_mass_fine":
        hist_sim2.GetXaxis().SetRangeUser(80, 100)

    hist_sim2.Draw("HIST")
    hist_sim1.Draw("HIST SAME")

    legend = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
    legend.AddEntry(hist_sim1, "M: 50+ GeV", "l")
    legend.AddEntry(hist_sim2, "M: 10-50 GeV", "l")
    legend.Draw()

    canvas.SaveAs(f"hist_true_data_hist/{sim_name}.png")

file_sim1.Close()
file_sim2.Close()
