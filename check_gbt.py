"""
check_gbt.py :: a script to help visualize the GBT data.

Run like:
> python check_gbt.py
"""

import argparse
import math
import os
import ROOT
ROOT.gROOT.SetBatch()

def main():

    rootlogon()
    ops = options()
    if not ops.i:
        fatal("Need input GBT ROOT file (-i) to run over")
    if not os.path.isfile(ops.i):
        fatal("Input file (-i) does not exist: %s" % (ops.i))
    if not ops.o or not ops.o.endswith(".pdf"):
        fatal("Need output pdf name (-o)")
    if not ops.l and not ops.n:
        fatal("Please give either a number of events (-n) or a list of event numbers (-l)")
    #if ops.l and ops.n:
    #    fatal("Please give either a number of events (-n) or a list of event numbers (-l), but not both")
    if ops.l and not os.path.isfile(ops.l):
        fatal("Please give list of event numbers (-l) which exists")

    # loop parameters
    events    = int(ops.n) if ops.n else 9e9
    events_ok = 0
    min_angle = 10*math.pi/180

    # i/o
    fname = ops.i # "/Volumes/TunaPassport/research/teststand/3522/intermediate/Run3522_GBT_decoded.root"
    froot = ROOT.TFile(ops.i)
    fout  = ROOT.TFile(ops.o.replace("pdf", "root"), "recreate")
    tr    = froot.Get("GBT_data")
    ents  = tr.GetEntries()
    
    # event list?
    list_of_events = []
    if ops.l:
        flist = open(ops.l)
        for line in flist:
            list_of_events.append( int(line.strip()) )
        events = min([events, len(list_of_events)])

    # histogram parameters
    hname  = "gbtmap_%06i"
    htitle = ";BCID;Board number;Hits in road"
    xbins, xlo, xhi = 14, -1.5, 12.5
    ybins, ylo, yhi = 8, -0.5, 7.5
    gbtmap = {}

    # announce
    print
    print "Input           : %s" % (ops.i)
    print "Entries         : %s" % (ents)
    print "Events to scan  : %s" % (events)
    print "Min. angle, deg : %s" % (min_angle/math.pi*180)
    print "Output          : %s" % (ops.o)
    print
    if ops.n and events > ents:
        fatal("N(events) > N(entries): %s vs %s" % (events, ents))

    # loop
    for ent in xrange(ents):

        _ = tr.GetEntry(ent)
        if events_ok >= events:
            break

        # i/o
        gbt_evt = tr.EventNum
        if ops.l and not gbt_evt in list_of_events:
            continue
        gbt_vmm = list(tr.gbt_VMM)
        gbt_ch  = list(tr.gbt_CH)
        gbt_fe  = list(tr.gbt_MMFE8)
        gbt_bc  = list(tr.gbt_BCID)
        gbt_n   = xrange(len(gbt_vmm))

        # make triggers
        trig = Trigger()
        for hit in zip(gbt_fe, gbt_vmm, gbt_ch, gbt_bc):
            trig.add(*hit)

        # pick your favorite
        roads = trig.roads()
        if not roads:
            continue
        most_hits = max([len(trig.filter(road)) for road in roads])
        for road in roads:
            if len(trig.filter(road)) == most_hits:
                hits  = trig.filter(road)
                hits  = sorted(hits, key=lambda tup: tup[3])
                if not ops.l:
                    if chi2_ndof(hits) > 10:
                        continue
                    if angle(hits) < min_angle:
                        continue
                    if trig.nx() < 2 or trig.nu()+trig.nv() < 2:
                        continue
                    if not trig.coin(road):
                        continue

                texs = []
                bcids = [bc for (_, _, _, bc) in hits]
                gbtmap[gbt_evt] = ROOT.TH2F(hname % gbt_evt, htitle, xbins, xlo, xhi, ybins, ylo, yhi)
                for (bo, vmm, ch, bc) in hits:
                    bc_adjusted = bc - min(bcids)
                    gbtmap[gbt_evt].Fill(bc_adjusted, bo)
                    texs.append(ROOT.TLatex(bc_adjusted,      bo+0.22, "%i" % (vmm)))
                    texs.append(ROOT.TLatex(bc_adjusted+0.01, bo+0.02, "#minus"))
                    texs.append(ROOT.TLatex(bc_adjusted,      bo-0.22, "%i" % (ch)))
                texs.append(ROOT.TLatex(8, 7.8, "GBT event %i" % (gbt_evt)))

                # canvas!
                canvname = "canv_%06i" % (gbt_evt)
                canv = ROOT.TCanvas(canvname, canvname, 800, 800)
                canv.Draw()
                style(gbtmap[gbt_evt])
                gbtmap[gbt_evt].Draw("colzsame")
                for (it, tex) in enumerate(texs):
                    tex.SetTextAlign(22)
                    tex.SetTextSize(0.04)
                    if (it % 3 == 1):
                        tex.SetTextSize(0.08)
                        tex.SetTextFont(42)
                    tex.Draw()

                ROOT.gPad.RedrawAxis()
                if   events_ok == 0:        canv.Print(ops.o+"(", "pdf")
                elif events_ok == events-1: canv.Print(ops.o+")", "pdf")
                else:                       canv.Print(ops.o,     "pdf")

                events_ok += 1
                break

    fout.Write()
    fout.Close()


class Trigger(object):

    def __init__(self):
        self.hits      = []
        self.nroads    = 16
        self.roadsize  = 64
        self.neighbors = [-1, 0, 1]

    def add(self, mmfe, vmm, ch, bcid):
        self.hits.append([mmfe2board(mmfe), vmm, ch, bcid])

    def has(self, board):
        for [bo, vmm, ch, bcid] in self.hits:
            if bo==board:
                return 1
        return 0

    def nx(self): 
        return self.has(0) + self.has(1) + self.has(6) + self.has(7)
    def nu(self): 
        return self.has(2) + self.has(4)
    def nv(self): 
        return self.has(3) + self.has(5)

    def roads(self):
        roads_satisfied = []
        for road in xrange(self.nroads):
            boards = []
            for (bo, vmm, ch, bc) in self.hits:
                road_hit = tp_strip(bo, vmm, ch) / self.roadsize
                if road_hit - road in self.neighbors:
                    boards.append(bo)
            boards = list(set(boards))
            if coincidence(boards):
                roads_satisfied.append(road)
        return sorted(roads_satisfied)

    def filter(self, road):
        these_hits = []
        for (bo, vmm, ch, bc) in self.hits:
            road_hit = tp_strip(bo, vmm, ch) / self.roadsize
            if road_hit - road in self.neighbors:
                these_hits.append( [bo, vmm, ch, bc] )
        return these_hits

    def coin(self, road):
        boards = []
        for (bo, vmm, ch, bc) in self.filter(road):
            boards.append(bo)
        return coincidence(boards)
        
def tp_strip(board, vmm, ch):
    strip = 64*vmm + ch - 1
    if board in [0, 3, 5, 6]:
        strip = 511 - strip
    if   board in [0, 1, 6, 7]: strip += 64
    elif board in [2, 4]:       strip += 58
    elif board in [3, 5]:       strip += 71
    return int(strip)
    
def coincidence(boards):
    horiz_ok  = (0 in boards or 1 in boards) and (6 in boards or 7 in boards)
    stereo_ok = any([(2 in boards or 4 in boards) and (3 in boards or 5 in boards),
                     (2 in boards and 3 in boards),
                     (4 in boards and 5 in boards)])
    return horiz_ok and stereo_ok

def chi2_ndof(hits):
    return 1.0

def angle(hits):
    return 1.0

def mmfe2board(mmfe):
    if   mmfe == 118: return 0
    elif mmfe == 116: return 1
    elif mmfe == 102: return 2
    elif mmfe == 119: return 3
    elif mmfe == 106: return 4
    elif mmfe == 107: return 5
    elif mmfe == 117: return 6
    elif mmfe == 105: return 7
    else: fatal("Cant convert mmfe2board on %s" % (mmfe))

def rootlogon():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetFillColor(10)
    import array
    ncontours = 200
    stops = array.array("d", [0.00, 0.50, 1.00])
    red   = array.array("d", [1.00, 1.00, 0.00])
    blue  = array.array("d", [1.00, 0.00, 0.00])
    green = array.array("d", [1.00, 1.00, 1.00])
    ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def style(hist):
    size = 0.045
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.6)
    hist.GetZaxis().SetTitleOffset(1.0)
    hist.GetZaxis().SetTitleSize(size)
    hist.GetZaxis().SetLabelSize(size)
    hist.GetYaxis().SetLabelOffset(0.01)
    hist.GetZaxis().SetNdivisions(502)

def fatal(msg):
    import sys
    sys.exit("Fatal error: %s" % (msg))
    
def options():
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input GBT decoded ROOT file")
    parser.add_argument("-o", help="Output ROOT or PDF file")
    parser.add_argument("-l", help="List of GBT EventNumbers to run over")
    parser.add_argument("-n", help="Number of GBT events to consider")
    return parser.parse_args()

if __name__ == "__main__":
    main()
