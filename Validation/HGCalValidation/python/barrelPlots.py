from __future__ import print_function
import os
import sys
import copy
import collections

import ROOT
from ROOT import TFile, TString
from ROOT import gDirectory
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from Validation.RecoTrack.plotting.plotting import Plot, PlotGroup, PlotFolder, Plotter, PlotOnSideGroup
from Validation.RecoTrack.plotting.html import PlotPurpose
import Validation.RecoTrack.plotting.plotting as plotting
import Validation.RecoTrack.plotting.validation as validation
import Validation.RecoTrack.plotting.html as html

from Validation.HGCalValidation.BarrelValidator_cfi import barrelValidator
from Validation.HGCalValidation.PostProcessorBarrel_cfi import lcToCP_linking#, simDict, variables

barrelVal_dqm = "DQMData/Run 1/CLUE-barrel/Run summary/BarrelValidator/"
geometryscenario = 88

layerscheme = { 'lastLayerEB' : 1,
		'lastLayerHB' : 5,
		'lastLayerHO' : 6,
		'maxLayer': 6}

lastLayerEB = layerscheme['lastLayerEB']
lastLayerHB = layerscheme['lastLayerHB']
lastLayerHO = layerscheme['lastLayerHO']
maxLayer = layerscheme['maxLayer']

_common = {"stat":True, "drawStyle": "hist", "staty": 0.65 }
_legend_common = {"legendDx": -0.3,
		  "legendDy": -0.05,
		  "legendDw": 0.1}

_SelectedCaloParticles = PlotGroup("SelectedCaloParticles", [
	Plot("num_caloparticle_eta", xtitle="", **_common),
	Plot("caloparticle_energy", xtitle="", **_common),
	Plot("caloparticle_pt", xtitle="", **_common),
	Plot("caloparticle_phi", xtitle="", **_common),
	Plot("Eta vs Zorigin", xtitle="", **_common),
       ])

_common = {"stat": True, "drawStyle": "hist", "statx": 0.38, "staty": 0.68 }
_num_reco_cluster_eta = PlotGroup("num_reco_cluster_eta", [
  Plot("num_reco_cluster_eta", xtitle="", **_common),
  ], ncols=1)


_common = {"stat": False, "drawStyle": "hist", "staty": 0.65 }

_totclusternum_layer_EB = PlotGroup("totclusternum_layer_EB", [
  Plot("totclusternum_layer_00", xtitle="", **_common)
  ], ncols=1)

_totclusternum_layer_HB = PlotGroup("totclusternum_layer_HB", [
  Plot("totclusternum_layer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerEB, lastLayerHB)
  ], ncols = 4)

_totclusternum_layer_HO = PlotGroup("totclusternum_layer_HO", [
  Plot("totclusternum_layer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerHB, lastLayerHO)
  ], ncols=1)

_energyclustered_perlayer_EB = PlotGroup("energyclustered_perlayer_EB", [
  Plot("energyclustered_perlayer00", xtitle="", **_common)], ncols=1)

_energyclustered_perlayer_HB = PlotGroup("energyclustered_perlayer_HB", [
  Plot("energyclustered_perlayer{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerEB, lastLayerHB)
], ncols=4)

_energyclustered_perlayer_HO = PlotGroup("energyclustered_perlayer_HO", [
  Plot("energyclustered_perlayer{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerHB, lastLayerHO)
], ncols=1)


_common_distance = {}
_common_distance.update(_common)
_common_distance.update(_legend_common)
_common_distance["xmax"] = 150
_common_distance["stat"] = False
_common_distance["ymin"] = 1e-3
_common_distance["ymax"] = 10000
_common_distance["ylog"] = True

_distancetomaxcell_perlayer_EB = PlotGroup("distancetomaxcell_perlayer_EB", [
  Plot("distancetomaxcell_perlayer_00", xtitle="", **_common)], ncols=1)

_distancetomaxcell_perlayer_HB = PlotGroup("distancetomaxcell_perlayer_HB", [
  Plot("distancetomaxcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerEB, lastLayerHB)
], ncols=4)

_distancetomaxcell_perlayer_HO = PlotGroup("distancetomaxcell_perlayer_HO", [
  Plot("distancetomaxcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerHB, lastLayerHO)
], ncols=1)

_distancebetseedandmaxcell_perlayer_EB = PlotGroup("distancebetseedandmaxcell_perlayer_EB", [
  Plot("distancebetseedandmaxcell_perlayer_00", xtitle="", **_common)], ncols=1)

_distancebetseedandmaxcell_perlayer_HB = PlotGroup("distancebetseedandmaxcell_perlayer_HB", [
  Plot("distancebetseedandmaxcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerEB, lastLayerHB)
], ncols=4)

_distancebetseedandmaxcell_perlayer_HO = PlotGroup("distancebetseedandmaxcell_perlayer_HO", [
  Plot("distancebetseedandmaxcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerHB, lastLayerHO)
], ncols=1)

_distancetoseedcell_perlayer_EB = PlotGroup("distancetoseedcell_perlayer_EB", [
  Plot("distancetoseedcell_perlayer_00", xtitle="", **_common)], ncols=1)

_distancetoseedcell_perlayer_HB = PlotGroup("distancetoseedcell_perlayer_HB", [
  Plot("distancetoseedcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerEB, lastLayerHB)
], ncols=4)

_distancetoseedcell_perlayer_HO = PlotGroup("distancetoseedcell_perlayer_HO", [
  Plot("distancetoseedcell_perlayer_{:02d}".format(i), xtitle="", **_common) for i in range(lastLayerHB, lastLayerHO)
], ncols=1)

#TODO insert energy weighted plots

_common = {"stat": True, "drawStyle": "hist", "staty": 0.65}

_common_score = {"title": "Score CaloParticle to LayerClusters", 
		 "stat": False,
		 "ymin": 0.1,
		 "ymax": 1000,
		 "xmin": 0,
		 "xmax": 1,
		 "drawStyle": "hist",
		 "lineWidth": 1,
		 "ylog": True
		}
_common_score.update(_legend_common)
_score_caloparticle_to_layercluster = PlotGroup("score_caloparticle_to_layercluster", [
  Plot("Score_caloparticle2layercl_perlayer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common) for i in range(lastLayerHO)], ncols=6)

_common_score = {"title": "Score LayerCluster to CaloParticles", 
		 "stat": False,
		 "ymin": 0.1,
		 "ymax": 1000,
		 "xmin": 0, 
		 "xmax": 1,
		 "drawStyle": "hist",
		 "lineWidth": 1,
		 "ylog": True
		}
_common_score.update(_legend_common)
_score_layercluster_to_caloparticles = PlotGroup("score_layercluster_to_caloparticle", [
  Plot("Score_layercl2caloparticle_perlayer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common) for i in range(lastLayerHO)], ncols=6)

_common_shared = {"title": "Shared Energy CaloParticle To Layer Cluster",
		  "stat": False,
		  "legend": False,
		  }
_common_shared.update(_legend_common)
_shared_plots = [Plot("SharedEnergy_caloparticle2layercl_perlayuer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common) for i in range(lastLayerHO)]
_shared_plots = [Plot("SharedEnergy_caloparticle2layercl_vs_eta_perlayer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common) for i in range(lastLayerHO)]
_shared_plots = [Plot("SharedEnergy_caloparticle2layercl_vs_phi_perlayer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common) for i in range(lastLayerHO)]
_sharedEnergy_caloparticle_to_layercluster = PlotGroup("sharedEnergy_caloparticle_to_layercluster", _shared_plots, ncols=6)

_common_assoc = {"stat": False,
		 "legend": False,
		 "xbinlables": ["", "TN(pur)", "FN(ineff.)", "FP(fake)", "TP(eff)"],
		 "xbinlabeloption": "h",
		 "drawStyle": "hist",
		 "ymin": 0.1,
		 "ymax": 10000,
		 "ylog": True}
_common_assoc.update(_legend_common)
_cell_association_table = PlotGroup("cellAssociaiton_table", [
	Plot("cellAssociation_perlayer{:02d}".format(i), xtitle="Layer {:02d}".format(i), **_common_assoc) for i in range(lastLayerHO)], ncols=6)

_bin_count = 0
_xbinlabels = [ "{:02d}".format(i) for i in range(lastLayerHO)]
_xtitle = "Layer Number"
_common_eff = {"stat":False, "legend": False, "ymin": 0.0, "ymax": 1.1, "xbinlabeloption": "d"}
_effplots_eta = [Plot("effic_eta_layer{:02d}".format(i), xtitle="", **_common_eff) for i in range(lastLayerHO)]
_effplots_phi = [Plot("effic_phi_layer{:02d}".format(i), xtitle="", **_common_eff) for i in range(lastLayerHO)]

_common_eff = {"stat": False, "legend": False, "xbinlabels": _xbinlabels, "xbinlablesize": 12, "xbinlabeloption": "v", "ymin": 0.0, "ymax": 1.1}
_common_eff["xmin"] = _bin_count
_common_eff["xmax"] = maxLayer

_effplots = [Plot("globalEfficiencies", xtitle=_xtitle, ytitle="Efficiency", **_common_eff)]
_efficiencies_eta = PlotGroup("Efficiencies_vs_eta", _effplots_eta, ncols=6)
_efficiencies_phi = PlotGroup("Efficiencies_vs_phi", _effplots_phi, ncols=6)
_efficiencies = PlotGroup("Efficiencies_vs_layer", _effplots, ncols=1)

_common_dup = {"stat": False, "legend": False, "ymin": 0.0, "ymax": 1.1}
_dupplots_eta = [Plot("duplicate_eta_layer{:02d}".format(i), xtitle="", **_common_dup) for i in range(lastLayerHO)]
_dupplots_phi = [Plot("duplicate_phi_layer{:02d}".format(i), xtitle="", **_common_dup) for i in range(lastLayerHO)]

_common_dup = {"stat": False, "legend": False, "title": "Global Duplicates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloption": "v", "ymin": 0.0, "ymax": 1.1}
_common_dup["xmin"] = _bin_count
_common_dup["xmax"] = maxLayer
_dupplots = [Plot("globaleEfficiencies", xtitle=_xtitle, ytitle="Duplicates", **_common_dup)]
_duplicates_eta = PlotGroup("Duplicates_vs_eta", _dupplots_eta, ncols=6)
_duplicates_phi = PlotGroup("Duplicates_vs_phi", _dupplots_phi, ncols=6)
_duplicates = PlotGroup("Duplicates_vs_layer", _dupplots, ncols=1)

_common_fake = {"stat": False, "legend": False, "ymin":0.0, "ymax":1.1}
_fakeplots_eta = [Plot("fake_eta_layer{:02d}".format(i), xtitle="", **_common_fake) for i in range(lastLayerHO)]
_fakeplots_phi = [Plot("fake_phi_layer{:02d}".format(i), xtitle="", **_common_fake) for i in range(lastLayerHO)]

_common_fake = {"stat":True, "legend": False, "title": "Global Fake Rates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloption": "v", "ymin": 0.0, "ymax": 1.1}
_common_fake["xmin"] = _bin_count
_common_fake["xmax"] = maxLayer
_common_fake["xbinlabelsize"] = 10.
_fakeplots = [Plot("gloablEfficiencies", xtitle=_xtitle, ytitle="Fake rate", **_common_fake)]
_fakes_eta = PlotGroup("FakeRate_vs_eta", _fakeplots_eta, ncols=6)
_fakes_phi = PlotGroup("FakeRate_vs_phi", _fakeplots_phi, ncols=6)
_fakes = PlotGroup("FakeRate_vs_layuer", _fakeplots, ncols=1)

_common_merge = {"stat": False, "legend": False, "ymin": 0.0, "ymax": 1.1}
_mergeplots_eta = [Plot("merge_eta_layer{:02d}".format(i), xtitle="", **_common_merge) for i in range(lastLayerHO)]
_mergeplots_phi = [Plot("merge_phi_layer{:02d}".format(i), xtitle="", **_common_merge) for i in range(lastLayerHO)]

_common_merge = {"stat": False, "legend": False, "title": "Global Merge Rates", "xbinlabels": _xbinlabels, "xbinlabelsize": 12, "xbinlabeloption": "v", "ymin": 0.0, "ymax": 1.1}
_common_merge["xmin"] = _bin_count
_common_merge["xmax"] = maxLayer
_common_merge["xbinlabelsize"] = 10.
_mergeplots = [Plot("globalEfficiencies", xtitle=_xtitle, ytitle="Merge Rate", **_common_merge)]
_merges_eta = PlotGroup("MergeRate_vs_eta", _mergeplots_eta, ncols=6)
_merges_phi = PlotGroup("MergeRate_vs_phi", _mergeplots_phi, ncols=6)
_merges = PlotGroup("MergeRate_vs_layer", _mergeplots, ncols=1)


_common_energy_score = dict(removeEmptyBins=False, xbinlablesize=10, legend=False,
    stat=True,
    xbinlableoption="d",
    ncols=1,
    xmin=0.001,
    xmax=1.,
    ymin=0.01,
    ymax=1.)

_energyscore_cp2lc = PlotGroup("Energy_vs_Score_CP2LC", [Plot("Energy_vs_Score_caloparticle2layer_perlayer{:02d}".format(i), title="Energy_vs_Score_CP2LC", xtitle="Layer {}".format(i), drawStyle="COLZ", adjustMarginRight=0.1, **_common_energy_score) for i in range(lastLayerHO)], ncols=6)

_common_energy_score["xmin"] = -0.1
_energyscore_lc2cp = PlotGroup("Energy_vs_Score_LC2CP", [Plot("Energy_vs_Score_layer2caloparticle_perlayer{:02d}".format(i), title="Energy_vs_Score_LC2CP", xtitle="Layer {}".format(i), drawStyle="COLZ", adjustMarginRight=0.1, **_common_energy_score) for i in range(lastLayerHO)], ncols=6)

#======================================================================
barrelLayerClustersPlotter = Plotter()
layerClustersLabel = "Layer Clusters"

lc_general_clusterlevel = [
  _totclusternum_layer_EB,
  _totclusternum_layer_HB,
  _totclusternum_layer_HO,
  _num_reco_cluster_eta,
  _energyclustered_perlayer_EB,
  _energyclustered_perlayer_HB,
  _energyclustered_perlayer_HO
]

lc_clusterlevel = [
  _totclusternum_layer_EB,
  _totclusternum_layer_HB,
  _totclusternum_layer_HO,
  _energyclustered_perlayer_EB,
  _energyclustered_perlayer_HB,
  _energyclustered_perlayer_HO
]

lc_cp_association = [
  _efficiencies,
  _efficiencies_eta,
  _efficiencies_phi,
  _duplicates,
  _duplicates_eta,
  _duplicates_phi,
  _fakes,
  _fakes_eta,
  _fakes_phi,
  _merges,
  _merges_eta,
  _merges_phi,
  _score_caloparticle_to_layercluster,
  _score_layercluster_to_caloparticles,
  #_sharedEnergy_layercluster_to_caloparticle,
  _sharedEnergy_caloparticle_to_layercluster,
  _energyscore_cp2lc,
  _energyscore_lc2cp
]

lc_extended = [
  _distancetomaxcell_perlayer_EB,
  _distancetomaxcell_perlayer_HB,
  _distancetomaxcell_perlayer_HO,
  _distancetoseedcell_perlayer_EB,
  _distancetoseedcell_perlayer_HB,
  _distancetoseedcell_perlayer_HO,
  _distancebetseedandmaxcell_perlayer_EB,
  _distancebetseedandmaxcell_perlayer_HB,
  _distancebetseedandmaxcell_perlayer_HO
]

def append_barrelLayerClustersPlots(collection = barrelValidator.label_layerClusterPlots._InputTag__moduleLabel, name_collection = layerClustersLabel, extended = False):
  print('extended : ', extended)
  regions_clusterLevel = ['General: Cluster Level']
  regions_LCtoCP_association = ['General: LC_CP association']

  plots_lc_general_clusterlevel = lc_general_clusterlevel
  plots_lc_clusterlevel_perlayer = lc_clusterlevel
  plots_lc_cp_association = lc_cp_association
  
  if extended:
    plots_lc_extended = lc_extended

  setPlots_ClusterLevel = [plots_lc_general_clusterlevel, plots_lc_clusterlevel_perlayer]
  setPlots_LCtoCPAssociation = [plots_lc_cp_association]

  for reg, setPlot in zip(regions_clusterLevel, setPlots_ClusterLevel):
    barrelLayerClustersPlotter.append(collection+"_", [
		      _barrelFolders(collection+"/ClusterLevel")
		      ], PlotFolder(
		      *setPlot,
		      loopSubFolders = False,
		      purpose=PlotPurpose.Timing, page=layerClustersLabel, section=reg))
  for reg, setPlot in zip(regions_LCtoCP_association, setPlots_LCtoCPAssociation):
    barrelLayerClustersPlotter.append(collection+"_", [
		      _barrelFolders(collection + "/" + lcToCP_linking)
		      ], PlotFolder(
		      *setPlot,
		      loopSubFolders=False,
		      purpose=PlotPurpose.Timing, page=layerClustersLabel, section=reg))

def _barrelFolders(lastDirName="barrelLayerClusters"):
  return barrelVal_dqm+lastDirName
