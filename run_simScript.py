#!/usr/bin/env python
import os
import sys
import ROOT

import shipunit as u
import shipRoot_conf
import rootUtils as ut
from ShipGeoConfig import ConfigRegistry
from argparse import ArgumentParser

DownScaleDiMuon = False

mcEngine = "TGeant4"
simEngine = "Pythia8"  # "Genie" # Ntuple

inclusive = "c"  # True = all processes if "c" only ccbar -> HNL, if "b" only bbar -> HNL, if "bc" only Bc+/Bc- -> HNL, and for darkphotons: if meson = production through meson decays, pbrem = proton bremstrahlung, qcd = ffbar -> DP.

MCTracksWithHitsOnly = False  # copy particles which produced a hit and their history
MCTracksWithEnergyCutOnly = True  # copy particles above a certain kin energy cut
MCTracksWithHitsOrEnergyCut = False  # or of above, factor 2 file size increase compared to MCTracksWithEnergyCutOnly

inputFile = "/eos/experiment/ship/data/Charm/Cascade-parp16-MSTP82-1-MSEL4-978Bpot.root"
defaultInputFile = True

parser = ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument(
    "--Pythia8", dest="pythia8", help="Use Pythia8", required=False, action="store_true"
)
parser.add_argument(
    "--PG", dest="pg", help="Use Particle Gun", required=False, action="store_true"
)
parser.add_argument(
    "--pID",
    dest="pID",
    help="id of particle used by the gun (default=22)",
    required=False,
    default=22,
    type=int,
)
parser.add_argument(
    "--Estart",
    help="start of energy range of particle gun (default=10 GeV)",
    default=10,
    type=float,
)
parser.add_argument(
    "--Eend",
    help="end of energy range of particle gun (default=10 GeV)",
    default=10,
    type=float,
)
parser.add_argument(
    "--FixedTarget", help="Use FixedTarget generator", action="store_true"
)
parser.add_argument(
    "--Ntuple",
    dest="ntuple",
    help="Use ntuple as input",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--MuonBack",
    dest="muonback",
    help="Generate events from muon background file, --Cosmics=0 for cosmic generator data",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--FollowMuon",
    dest="followMuon",
    help="Make muonshield active to follow muons",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--FastMuon",
    dest="fastMuon",
    help="Only transport muons for a fast muon only background estimate",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--phiRandom",
    dest="phiRandom",
    help="only relevant for muon background generator, random phi",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--Cosmics",
    dest="cosmics",
    help="Use cosmic generator, argument switch for cosmic generator 0 or 1",
    required=False,
    default=None,
)
parser.add_argument(
    "-n",
    "--nEvents",
    dest="nEvents",
    help="Number of events to generate",
    required=False,
    default=100,
    type=int,
)
parser.add_argument(
    "-i",
    "--firstEvent",
    dest="firstEvent",
    help="First event of input file to use",
    required=False,
    default=0,
    type=int,
)
parser.add_argument(
    "-s",
    "--seed",
    dest="theSeed",
    help="Seed for random number. Only for experts, see TRrandom::SetSeed documentation",
    required=False,
    default=0,
    type=int,
)
parser.add_argument(
    "-S",
    "--sameSeed",
    dest="sameSeed",
    help="can be set to an integer for the muonBackground simulation with specific seed for each muon, only for experts!",
    required=False,
    default=False,
    type=int,
)
group.add_argument(
    "-f",
    dest="inputFile",
    help="Input file if not default file",
    required=False,
    default=False,
)
parser.add_argument(
    "-g",
    dest="geofile",
    help="geofile for muon shield geometry, for experts only",
    required=False,
    default=None,
)
parser.add_argument(
    "-o",
    "--output",
    dest="outputDir",
    help="Output directory",
    required=False,
    default=".",
)
parser.add_argument(
    "-F",
    dest="deepCopy",
    help="default = False: copy only stable particles to stack, except for HNL events",
    required=False,
    action="store_true",
)
parser.add_argument(
    "-t",
    "--test",
    dest="testFlag",
    help="quick test",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--dry-run",
    dest="dryrun",
    help="stop after initialize",
    required=False,
    action="store_true",
)
parser.add_argument(
    "-D",
    "--display",
    dest="eventDisplay",
    help="store trajectories",
    required=False,
    action="store_true",
)
parser.add_argument(
    "--debug",
    help="1: print weights and field 2: make overlap check",
    required=False,
    default=0,
    type=int,
    choices=range(0, 3),
)

options = parser.parse_args()

if options.pythia8:
    simEngine = "Pythia8"
if options.pg:
    simEngine = "PG"
if options.FixedTarget:
    simEngine = "FixedTarget"
if options.ntuple:
    simEngine = "Ntuple"
if options.muonback:
    simEngine = "MuonBack"
if options.cosmics:
    simEngine = "Cosmics"
    Opt_high = int(options.cosmics)
if options.inputFile:
    if options.inputFile == "none":
        options.inputFile = None
    inputFile = options.inputFile
    defaultInputFile = False
if options.testFlag:
    inputFile = "$FAIRSHIP/files/Cascade-parp16-MSTP82-1-MSEL4-76Mpot_1_5000.root"


print("FairShip setup for", simEngine, "to produce", options.nEvents, "events")
if (simEngine == "Ntuple" or simEngine == "MuonBack") and defaultInputFile:
    print("input file required if simEngine = Ntuple or MuonBack")
    print(
        " for example -f /eos/experiment/ship/data/Mbias/pythia8_Geant4-withCharm_onlyMuons_4magTarget.root"
    )
    sys.exit()
ROOT.gRandom.SetSeed(
    options.theSeed
)  # this should be propagated via ROOT to Pythia8 and Geant4VMC
shipRoot_conf.configure(0)  # load basic libraries, prepare atexit for python
ship_geo = ConfigRegistry.loadpy(
    "$FAIRSHIP../ruffiano/charm-geometry_config.py",
)

# Output file name, add dy to be able to setup geometry with ambiguities.
if simEngine == "PG":
    tag = simEngine + "_" + str(options.pID) + "-" + mcEngine
else:
    tag = simEngine + "-" + mcEngine
if options.eventDisplay:
    tag = tag + "_D"
if not os.path.exists(options.outputDir):
    os.makedirs(options.outputDir)
outFile = f"{options.outputDir}/ship.{tag}.root"

# rm older files !!!
for x in os.listdir(options.outputDir):
    if not x.find(tag) < 0:
        os.system(f"rm {options.outputDir}/{x}")
# Parameter file name
parFile = f"{options.outputDir}/ship.params.{tag}.root"

# In general, the following parts need not be touched
# ========================================================================

# -----Timer--------------------------------------------------------
timer = ROOT.TStopwatch()
timer.Start()
# ------------------------------------------------------------------------
# -----Create simulation run----------------------------------------
run = ROOT.FairRunSim()
run.SetName(mcEngine)  # Transport engine
run.SetSink(ROOT.FairRootFileSink(outFile))  # Output file
run.SetUserConfig("g4Config.C")  # user configuration file default g4Config.C
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
# import shipMuShield_only as shipDet_conf # special use case for an attempt to convert active shielding geometry for use with FLUKA
# import shipTarget_only as shipDet_conf
import charmDet_conf as shipDet_conf

modules = shipDet_conf.configure(run, ship_geo)


def addScoringPlane(
    anindex=0, xpos=0.0, ypos=0.0, zpos=0.0, xhalfw=500.0, yhalfh=500.0
):
    izstring = "ScoringPlane" + str(anindex)
    # put the 3rd arg as True if you want to stop particles being tracked at this plane
    scoringplane = ROOT.ScoringPlane(
        izstring, ROOT.kTRUE, ROOT.kFALSE, xhalfw, yhalfh, 0.1
    )
    scoringplane.SetVetoPointName("sco" + str(anindex) + "_")
    scoringplane.SetXYZposition(xpos, ypos, zpos)
    print(
        "    defined "
        + izstring
        + " at x,y,z = "
        + str(xpos)
        + " , "
        + str(ypos)
        + " , "
        + str(zpos)
        + " cm (halfW/halfH = "
        + str(xhalfw)
        + " , "
        + str(yhalfh)
        + ")"
    )
    return scoringplane


scoring_planes = [
    addScoringPlane(1, zpos=0.5 * u.m, xhalfw=1 * u.m, yhalfh=1 * u.m),
    addScoringPlane(2, zpos=1 * u.m, xhalfw=1 * u.m, yhalfh=1 * u.m),
    addScoringPlane(3, zpos=6 * u.m, xhalfw=1 * u.m, yhalfh=1 * u.m),
    addScoringPlane(4, zpos=6.5 * u.m, xhalfw=1 * u.m, yhalfh=1 * u.m),
]

for scoring_plane in scoring_planes:
    run.AddModule(scoring_plane)

# -----Create PrimaryGenerator--------------------------------------
primGen = ROOT.FairPrimaryGenerator()
if simEngine == "Pythia8":
    primGen.SetTarget(ship_geo.target.z0, 0.0)
    # -----Pythia8--------------------------------------
    primGen.SetTarget(0.0, 0.0)  # vertex is set in pythia8Generator
    ut.checkFileExists(inputFile)
    if ship_geo.Box.gausbeam:
        primGen.SetBeam(
            0.0, 0.0, 0.5, 0.5
        )  # more central beam, for hits in downstream detectors
        primGen.SmearGausVertexXY(True)  # sigma = x
    else:
        primGen.SetBeam(
            0.0, 0.0, ship_geo.Box.TX - 1.0, ship_geo.Box.TY - 1.0
        )  # Uniform distribution in x/y on the target (0.5 cm of margin at both sides)
        primGen.SmearVertexXY(True)
    P8gen = ROOT.Pythia8Generator()
    P8gen.UseExternalFile(inputFile, options.firstEvent)
    P8gen.SetTarget(
        "target_1", 0.0, 0.0
    )  # will distribute PV inside target, beam offset x=y=0.
    primGen.AddGenerator(P8gen)
if simEngine == "FixedTarget":
    P8gen = ROOT.FixedTargetGenerator()
    P8gen.SetTarget("target_1", 0.0, 0.0)
    P8gen.SetMom(400.0 * u.GeV)
    P8gen.SetEnergyCut(0.0)
    P8gen.SetHeartBeat(100000)
    P8gen.SetG4only()
    primGen.AddGenerator(P8gen)

# -----Particle Gun-----------------------
if simEngine == "PG":
    myPgun = ROOT.FairBoxGenerator(options.pID, 1)
    myPgun.SetPRange(options.Estart, options.Eend)
    myPgun.SetPhiRange(0, 360)  # // Azimuth angle range [degree]
    myPgun.SetXYZ(0.0 * u.cm, 0.0 * u.cm, 0.0 * u.cm)
    myPgun.SetThetaRange(0, 0)  # // Polar angle in lab system range [degree]
    primGen.AddGenerator(myPgun)
if simEngine == "Ntuple":
    # reading previously processed muon events, [-50m - 50m]
    ut.checkFileExists(inputFile)
    primGen.SetTarget(ship_geo.target.z0 + 50 * u.m, 0.0)
    Ntuplegen = ROOT.NtupleGenerator()
    Ntuplegen.Init(inputFile, options.firstEvent)
    primGen.AddGenerator(Ntuplegen)
    options.nEvents = min(options.nEvents, Ntuplegen.GetNevents())
    print("Process ", options.nEvents, " from input file")
#
if simEngine == "MuonBack":
    # reading muon tracks from previous Pythia8/Geant4 simulation with charm replaced by cascade production
    fileType = ut.checkFileExists(inputFile)
    if fileType == "tree":
        # 2018 background production
        primGen.SetTarget(ship_geo.target.z0 + 70.845 * u.m, 0.0)
    else:
        primGen.SetTarget(ship_geo.target.z0 + 50 * u.m, 0.0)
    #
    MuonBackgen = ROOT.MuonBackGenerator()
    # MuonBackgen.FollowAllParticles() # will follow all particles after hadron absorber, not only muons
    MuonBackgen.Init(inputFile, options.firstEvent, options.phiRandom)
    MuonBackgen.SetSmearBeam(5 * u.cm)  # radius of ring, thickness 8mm
    if DownScaleDiMuon:
        testf = ROOT.TFile.Open(inputFile)
        if not testf.FileHeader.GetTitle().find("diMu100.0") < 0:
            MuonBackgen.SetDownScaleDiMuon()  # avoid interference with boosted channels
            print("MuonBackgenerator: set downscale for dimuon on")
        testf.Close()
    if options.sameSeed:
        MuonBackgen.SetSameSeed(options.sameSeed)
    primGen.AddGenerator(MuonBackgen)
    options.nEvents = min(options.nEvents, MuonBackgen.GetNevents())
    MCTracksWithHitsOnly = True  # otherwise, output file becomes too big
    print(
        "Process ",
        options.nEvents,
        " from input file, with Phi random=",
        options.phiRandom,
        " with MCTracksWithHitsOnly",
        MCTracksWithHitsOnly,
    )
    if options.followMuon:
        options.fastMuon = True
        modules["Veto"].SetFollowMuon()
    if options.fastMuon:
        modules["Veto"].SetFastMuon()

    # optional, boost gamma2muon conversion
    # ROOT.kShipMuonsCrossSectionFactor = 100.
#
run.SetGenerator(primGen)
# ------------------------------------------------------------------------

# ---Store the visualiztion info of the tracks, this make the output file very large!!
# --- Use it only to display but not for production!
if options.eventDisplay:
    run.SetStoreTraj(ROOT.kTRUE)
else:
    run.SetStoreTraj(ROOT.kFALSE)

# -----Configure external decayer globally------------------------------------
run.SetPythiaDecayer('DecayConfigTEvtGen.C')
print('Using TEvtGenDecayer for J/psi and quarkonium decays with EvtGen')

# -----Initialize simulation run------------------------------------
run.Init()
if options.dryrun:  # Early stop after setting up Pythia 8
    sys.exit(0)
gMC = ROOT.TVirtualMC.GetMC()
fStack = gMC.GetStack()
EnergyCut = 100.0 * u.MeV

if MCTracksWithHitsOnly:
    fStack.SetMinPoints(1)
    fStack.SetEnergyCut(-100.0 * u.MeV)
elif MCTracksWithEnergyCutOnly:
    fStack.SetMinPoints(-1)
    fStack.SetEnergyCut(EnergyCut)
elif MCTracksWithHitsOrEnergyCut:
    fStack.SetMinPoints(1)
    fStack.SetEnergyCut(EnergyCut)
elif options.deepCopy:
    fStack.SetMinPoints(0)
    fStack.SetEnergyCut(0.0 * u.MeV)

if options.eventDisplay:
    # Set cuts for storing the trajectories, can only be done after initialization of run (?!)
    trajFilter = ROOT.FairTrajFilter.Instance()
    trajFilter.SetStepSizeCut(1 * u.mm)
    trajFilter.SetVertexCut(
        -20 * u.m,
        -20 * u.m,
        ship_geo.target.z0 - 1 * u.m,
        20 * u.m,
        20 * u.m,
        200.0 * u.m,
    )
    trajFilter.SetMomentumCutP(0.1 * u.GeV)
    trajFilter.SetEnergyCut(0.0, 400.0 * u.GeV)
    trajFilter.SetStorePrimaries(ROOT.kTRUE)
    trajFilter.SetStoreSecondaries(ROOT.kTRUE)

# The VMC sets the fields using the "/mcDet/setIsLocalMagField true" option in "gconfig/g4config.in"
import geomGeant4
# geomGeant4.setMagnetField() # replaced by VMC, only has effect if /mcDet/setIsLocalMagField  false

# Print VMC fields and associated geometry objects
if options.debug == 1:
    geomGeant4.printVMCFields()
    geomGeant4.printWeightsandFields(
        onlyWithField=True,
        exclude=[
            "DecayVolume",
            "Tr1",
            "Tr2",
            "Tr3",
            "Tr4",
            "Veto",
            "Ecal",
            "Hcal",
            "MuonDetector",
            "SplitCal",
        ],
    )
# Plot the field example
# fieldMaker.plotField(1, ROOT.TVector3(-9000.0, 6000.0, 50.0), ROOT.TVector3(-300.0, 300.0, 6.0), 'Bzx.png')
# fieldMaker.plotField(2, ROOT.TVector3(-9000.0, 6000.0, 50.0), ROOT.TVector3(-400.0, 400.0, 6.0), 'Bzy.png')

# -----Start run----------------------------------------------------
run.Run(options.nEvents)
# -----Runtime database---------------------------------------------
kParameterMerged = ROOT.kTRUE
parOut = ROOT.FairParRootFileIo(kParameterMerged)
parOut.open(parFile)
rtdb.setOutput(parOut)
rtdb.saveOutput()
rtdb.printParamContexts()
getattr(rtdb, "print")()
# ------------------------------------------------------------------------
run.CreateGeometryFile(f"{options.outputDir}/geofile_full.{tag}.root")
# save ShipGeo dictionary in geofile
import saveBasicParameters

saveBasicParameters.execute(f"{options.outputDir}/geofile_full.{tag}.root", ship_geo)

# checking for overlaps
if options.debug == 2:
    fGeo = ROOT.gGeoManager
    fGeo.SetNmeshPoints(10000)
    fGeo.CheckOverlaps(0.1)  # 1 micron takes 5minutes
    fGeo.PrintOverlaps()
    # check subsystems in more detail
    for x in fGeo.GetTopNode().GetNodes():
        x.CheckOverlaps(0.0001)
        fGeo.PrintOverlaps()
# -----Finish-------------------------------------------------------
timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print(" ")
print("Macro finished successfully.")
print("Output file is ", outFile)
print("Parameter file is ", parFile)
print("Real time ", rtime, " s, CPU time ", ctime, "s")

# remove empty events
if simEngine == "MuonBack":
    tmpFile = outFile + "tmp"
    xxx = outFile.split("/")
    check = xxx[len(xxx) - 1]
    fin = False
    for ff in ROOT.gROOT.GetListOfFiles():
        nm = ff.GetName().split("/")
        if nm[len(nm) - 1] == check:
            fin = ff
    if not fin:
        fin = ROOT.TFile.Open(outFile)
    t = fin.Get("cbmsim")
    fout = ROOT.TFile(tmpFile, "recreate")
    fSink = ROOT.FairRootFileSink(fout)

    sTree = t.CloneTree(0)
    nEvents = 0
    pointContainers = []
    for x in sTree.GetListOfBranches():
        name = x.GetName()
        if not name.find("Point") < 0:
            pointContainers.append(
                "sTree." + name + ".GetEntries()"
            )  # makes use of convention that all sensitive detectors fill XXXPoint containers
    for n in range(t.GetEntries()):
        rc = t.GetEvent(n)
        empty = True
        for x in pointContainers:
            if eval(x) > 0:
                empty = False
        if not empty:
            rc = sTree.Fill()
            nEvents += 1

    branches = ROOT.TList()
    branches.SetName("BranchList")
    branches.Add(ROOT.TObjString("MCTrack"))
    branches.Add(ROOT.TObjString("vetoPoint"))
    branches.Add(ROOT.TObjString("ShipRpcPoint"))
    branches.Add(ROOT.TObjString("TargetPoint"))
    branches.Add(ROOT.TObjString("TTPoint"))
    branches.Add(ROOT.TObjString("ScoringPoint"))
    branches.Add(ROOT.TObjString("strawtubesPoint"))
    branches.Add(ROOT.TObjString("EcalPoint"))
    branches.Add(ROOT.TObjString("sEcalPointLite"))
    branches.Add(ROOT.TObjString("smuonPoint"))
    branches.Add(ROOT.TObjString("TimeDetPoint"))
    branches.Add(ROOT.TObjString("MCEventHeader"))
    branches.Add(ROOT.TObjString("sGeoTracks"))

    sTree.AutoSave()
    fSink.WriteObject(branches, "BranchList", ROOT.TObject.kSingleKey)
    fSink.SetOutTree(sTree)

    fout.Close()
    print("removed empty events, left with:", nEvents)
    rc1 = os.system("rm  " + outFile)
    rc2 = os.system("mv " + tmpFile + " " + outFile)
    fin.SetWritable(False)  # bpyass flush error

# ------------------------------------------------------------------------
import checkMagFields


def visualizeMagFields():
    checkMagFields.run()


def checkOverlapsWithGeant4():
    # after /run/initialize, but prints warning messages, problems with TGeo volume
    mygMC = ROOT.TGeant4.GetMC()
    mygMC.ProcessGeantCommand("/geometry/test/recursion_start 0")
    mygMC.ProcessGeantCommand("/geometry/test/recursion_depth 2")
    mygMC.ProcessGeantCommand("/geometry/test/run")
