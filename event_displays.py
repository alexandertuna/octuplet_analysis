import argparse
import copy
import os
import sys
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

def main():

    # setup
    ops = options()
    if not ops.i:
        fatal("Please provide an input ROOT file with -i")
    if not os.path.isfile(ops.i):
        fatal("Input ROOT file does not exist: %i" % (ops.i))

    # i/o
    fi = ROOT.TFile.Open(ops.i)
    pdf = "event_displays.pdf"

    # get events
    events = []
    tdir = fi.Get("event_displays")
    for key in tdir.GetListOfKeys():
        name = key.GetName()
        ev = name.split("_")[1]
        events.append(ev)
    events = sorted(list(set(events)), key=lambda ev: int(ev))

    # draw
    for ev in events:
        first = ev==events[0]
        print "Drawing %s :: %i / %i" % (ev, events.index(ev)+1, len(events))
        draw(fi, ev, pdf, first)

    close_pdf(pdf)

def draw(fi, ev, pdf, first):

    rootlogon()
    black = ROOT.kGray+3
    mmall = copy.deepcopy(fi.Get("event_displays/track2D_%s_MMall" % (ev)))
    mmfit = copy.deepcopy(fi.Get("event_displays/track2D_%s_MMfit" % (ev)))
    tpfit = copy.deepcopy(fi.Get("event_displays/track2D_%s_TPfit" % (ev)))

    gr_tp = None
    gr_mm = None
    gr_mf = None
    li_mm = None
    x_avg = None

    # draw the MM event
    for prim in mmall.GetListOfPrimitives():
        if isinstance(prim, ROOT.TGraph):
            prim.SetMarkerColor(black)
            prim.SetLineColor(black)
            gr_mm = copy.copy(prim)
        if isinstance(prim, ROOT.TMultiGraph):
            prim.GetXaxis().SetLabelColor(black)
            prim.GetXaxis().SetTitleColor(black)
            prim.GetYaxis().SetLabelColor(black)
            prim.GetYaxis().SetTitleColor(black)
            prim.GetXaxis().SetAxisColor(black)
            prim.GetYaxis().SetAxisColor(black)
            for gr in prim.GetListOfGraphs():
                if gr.GetN()==2 and abs(gr.GetY()[0]-gr.GetY()[1]) < 0.1:
                    gr.SetLineColor(black)
                else:
                    gr.SetLineWidth(1)
                    gr.SetLineColor(black)
                    li_mm = copy.copy(gr)
                    x_avg = average(gr.GetX(), gr.GetN())

    mmall.SetFillStyle(4000)
    mmall.Draw()

    # draw the TP
    for prim in tpfit.GetListOfPrimitives():
        if isinstance(prim, ROOT.TGraph):
            prim.SetMarkerColor(ROOT.kRed)
            prim.SetMarkerStyle(24)
            prim.SetMarkerSize(2)
            prim.Draw("P")
            gr_tp = copy.copy(prim)
            break

    # draw the MM fit
    for prim in mmfit.GetListOfPrimitives():
        if isinstance(prim, ROOT.TGraph):
            prim.SetMarkerColor(ROOT.kBlue)
            prim.SetMarkerStyle(5)
            prim.SetMarkerSize(2)
            #prim.Draw("P")
            gr_mf = copy.copy(prim)
            break

    # keep the MM graph on top
    gr_mm.Draw("P")

    # legend
    x1 = 0.63 if x_avg < 100 else 0.28
    leg = ROOT.TLegend(x1, 0.45, x1+0.15, 0.62)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextColor(black)
    leg.SetTextFont(132)
    leg.SetTextSize(0.04)
    leg.AddEntry(gr_mm, " MM cluster", "p")
    leg.AddEntry(li_mm, " MM track",   "l")
    leg.AddEntry(gr_tp, " TP hit",     "p")
    leg.Draw()

    # annotate
    descr = ROOT.TLatex(0.2, 0.96, "Event %s" % (ev.lstrip("0")))
    descr.SetNDC()
    descr.SetTextColor(black)
    descr.SetTextSize(0.04)
    descr.SetTextFont(132)
    descr.Draw("same")

    # save
    if first:
        mmall.Print(pdf+"(", "pdf")
    else:
        mmall.Print(pdf, "pdf")

def average(shit, n):
    return sum([shit[it] for it in xrange(n)]) / float(n)

def rootlogon():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetFillColor(10)

def close_pdf(pdf):
    rootlogon()
    canv = ROOT.TCanvas("close", "close", 800, 800)
    canv.Draw()
    canv.Print(pdf+")", "pdf")

def options():
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input ROOT file")
    return parser.parse_args()

if __name__ == "__main__":
    main()

