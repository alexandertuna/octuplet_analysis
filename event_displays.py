import argparse
import copy
import os
import sys
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

zoom = False

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
    storage = []

    # draw
    first = True
    for ev in events:
        #if not is_overlap(ev):
        #    continue
        print "Drawing %s :: %i / %i" % (ev, events.index(ev)+1, len(events))
        storage.append( draw(fi, ev, pdf, first) )
        first = False

    close_pdf(pdf)

def draw(fi, ev, pdf, first):


    rootlogon()
    black = ROOT.kGray+3
    mmall = fi.Get("event_displays/track2D_%s_MMall" % (ev))
    mmfit = fi.Get("event_displays/track2D_%s_MMfit" % (ev))
    tpfit = fi.Get("event_displays/track2D_%s_TPfit" % (ev))

    gr_tp = None
    gr_mm = None
    gr_mf = None
    li_mm = None
    x_avg = None
    xrang = (None, None)
    yrang = (None, None)

    # draw the MM event
    for prim in mmall.GetListOfPrimitives():
        if isinstance(prim, ROOT.TGraph):
            prim.SetMarkerColor(black)
            prim.SetLineColor(black)
            # gr_mm = copy.copy(prim)
            gr_mm = prim
        if isinstance(prim, ROOT.TMultiGraph):
            prim.GetXaxis().SetLabelColor(black)
            prim.GetXaxis().SetTitleColor(black)
            prim.GetYaxis().SetLabelColor(black)
            prim.GetYaxis().SetTitleColor(black)
            prim.GetXaxis().SetAxisColor(black)
            prim.GetYaxis().SetAxisColor(black)
            xrang = (prim.GetXaxis().GetXmin(), prim.GetXaxis().GetXmax())
            yrang = (prim.GetYaxis().GetXmin(), prim.GetYaxis().GetXmax())
            for gr in prim.GetListOfGraphs():
                # horizontal lines
                if gr.GetN()==2 and abs(gr.GetY()[0]-gr.GetY()[1]) < 0.1:
                    gr.SetLineColor(black)
                    gr.SetLineColor(black)
                # the fit
                else:
                    gr.SetLineWidth(2)
                    gr.SetLineColor(black)
                    # li_mm = copy.copy(gr)
                    li_mm = gr
                    x_avg = average(gr.GetX(), gr.GetN())

    mmall.SetFillStyle(4000)
    mmall.Draw()

    # draw the TP
    if tpfit:
        for prim in tpfit.GetListOfPrimitives():
            if isinstance(prim, ROOT.TGraph):
                prim.SetMarkerColor(ROOT.kRed)
                prim.SetMarkerStyle(24)
                prim.SetMarkerSize(2)
                prim.Draw("P")
                # gr_tp = copy.copy(prim)
                gr_tp = prim
                break

    # draw the MM fit
    pxs, pys = [], []
    for prim in mmfit.GetListOfPrimitives():
        if isinstance(prim, ROOT.TGraph):
            prim.SetMarkerColor(ROOT.kBlue)
            prim.SetMarkerStyle(5)
            prim.SetMarkerSize(2)
            # prim.Draw("P")
            # gr_mf = copy.copy(prim)
            gr_mf = prim
            for i in xrange(gr_mf.GetN()):
                px = ROOT.Double(0)
                py = ROOT.Double(0)
                _ = gr_mf.GetPoint(i, px, py)
                pxs.append(float(px))
                pys.append(float(py))
            break

    # draw roads?
    box = []
    roadsize = 8
    strip2mm = 0.4
    mmfeavgx = average2( filter(lambda px: isX( pys[pxs.index(px)] ), pxs) )
    xmin, xmax = 9999, -9999
    xdn, xup = xrang
    ydn, yup = yrang
    for ro in xrange(1000):
        roaddn = roadsize*strip2mm*ro
        roadup = roadsize*strip2mm*(ro+1)
        if mmfeavgx > roaddn and mmfeavgx < roadup:
            box.append( ROOT.TBox(roaddn, ydn, roadup, yup) )
            box[-1].SetFillColorAlpha(ROOT.kWhite, 0.1)
            for neighb in xrange(-3, 4):
                if neighb == 0:
                    continue
                delta   = roadsize*strip2mm*neighb
                overlap = roadsize*strip2mm/2
                # road
                box.append( ROOT.TBox(roaddn+delta, ydn, roadup+delta, yup) )
                box[-1].SetFillColorAlpha(ROOT.kWhite if neighb%2==0 else ROOT.kGray, 0.3)
                # overlap
                box.append( ROOT.TBox(roaddn+delta+overlap, ydn, roadup+delta+overlap, yup) )
                box[-1].SetFillColorAlpha(ROOT.kWhite, 0.0)
                box[-1].SetLineStyle(2)
                xmin = min(xmin, roaddn+delta)
                xmax = max(xmax, roadup+delta)
            break
    else:
        fatal("Couldnt find a road which contained mmfeavgx = %s" % (mmfeavgx))
        
    # keep the MM graph on top
    gr_mm.Draw("P")

    # draw order
    for prim in mmall.GetListOfPrimitives():
        if zoom and isinstance(prim, ROOT.TMultiGraph):
            padding = 3*roadsize*strip2mm
            prim.GetXaxis().SetLimits(xmin-padding, xmax+padding)
    mmall.Draw()
    if zoom:
        for bx in box:
            bx.Draw("l")
    if tpfit:
        gr_tp.Draw("P")
    # gr_mf.Draw("P")
    gr_mm.Draw("P")
    mmall.Draw()
    for prim in mmall.GetListOfPrimitives():
        if isinstance(prim, ROOT.TMultiGraph):
            prim.Draw()

    # legend
    if zoom:
        x1 = 0.18
    else:
        x1 = 0.63 if x_avg < 100 else 0.28
    leg = ROOT.TLegend(x1, 0.45, x1+0.27, 0.62)
    leg.SetBorderSize(1 if zoom else 0)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetFillStyle(1001 if zoom else 0)
    leg.SetTextColor(black)
    leg.SetTextFont(132)
    leg.SetTextSize(0.04)
    leg.AddEntry(gr_mm, " MM cluster", "p")
    leg.AddEntry(li_mm, " MM track",   "l")
    if tpfit:
        leg.AddEntry(gr_tp, " TP hit", "p")
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

    return (mmall, mmfit, tpfit)

def isX(zpos):
    return zpos < 20 or zpos > 130

def fatal(msg):
    sys.exit("Fatal: %s" % (msg))

def average2(shit):
    return sum(shit) / len(shit)

def average(shit, n):
    # for TGraphs
    return sum([shit[it] for it in xrange(n)]) / float(n)

def is_overlap(ev):
    return int(ev) in [
        3179,
        6177,
        6877,
        17544,
        22154,
        23454,
        49852,
        52259,
        53356,
        53725,
        56479,
        62070,
        63145,
        64304,
        68644,
        77273,
        78001,
        78061,
        80915,
        82149,
        84274,
        86822,
        87946,
        95852,
        96749,
        105393,
        106616,
        106944,
        112277,
        114175,
        118417,
        126944,
        127974,
        138513,
        139776,
        140299,
        141283,
        ]

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

