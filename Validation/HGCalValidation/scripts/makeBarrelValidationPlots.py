#!/usr/bin/env python3

from __future__ import print_function
import os
import argparse
from time import time

#from Validation.BarrelValidation.PostProcessorBarrel_cfi import tracksterLabels as trackstersIters

from Validation.RecoTrack.plotting.validation import SeparateValidation, SimpleValidation, SimpleSample
from Validation.HGCalValidation.BarrelValidator_cfi import barrelValidator
import Validation.HGCalValidation.barrelPlots as barrelPlots
import Validation.RecoTrack.plotting.plotting as plotting

simClustersIters = [
    barrelValidator.label_SimClustersLevel._InputTag__moduleLabel, "ticlSimTracksters"]

hitCalLabel = 'hitCalibration'
hitValLabel = 'hitValidation'
layerClustersLabel = 'barrelLayerClusters'
trackstersLabel = 'tracksters'
trackstersWithEdgesLabel = 'trackstersWithEdges'
simLabel = 'simulation'
allLabel = 'all'

collection_choices = [allLabel]
collection_choices.extend(#[hitCalLabel]+[hitValLabel]+
			 [layerClustersLabel]
                         #+[trackstersLabel]+[trackstersWithEdges]+[simLabel]
			 )



def main(opts):
    drawArgs = {}
    extendedFlag = False
    if opts.no_ratio:
        drawArgs['ratio'] = False
    if opts.separate:
        drawArgs['separate'] = True
    if opts.png:
        drawArgs['saveFormat'] = '.png'
    if opts.extended:
        extendedFlag = True
    if opts.verbose:
        plotting.verbose = True

    filenames = [(f, f.replace('.root', '')) for fi in opts.files]
    sample = SimpleSample(opts.subdirprefix[0], opts.html_sample, filenames)

    val = SimpleValidation([sample], opts.outputDir[0])
    if opts.separate:
        val = SeparateValidation([sample], opts.outputDir[0])
    htmlReport = val.createHtmlReport(
        validationName=opts.html_validation_name[0])

    # layerClusters
    def plot_LC():
        barrellayclus = [barrelPlots.barrelLayerClustersPlotter]
        barrelPlots.append_barrelLayerClustersPlots(
            barrelValidator.label_layerClusterPlots._InputTag__moduleLabel, "Layer Clusters", extendedFlag)
        val.doPlots(barrellayclus, plotterDrawArgs=drawArgs)

    plotDict = {layerClustersLabel: [plot_LC]}

    if (opts.collection != allLabel):
        for task in plotDict[opts.collection]:
            task()
    else:
        for label in plotDict:
            if (label == trackstersLabel):
                continue
            for task in plotDict[label]:
                task()

    if opts.no_html:
        print("Plots created into directory '%s'." % opts.outputDir)
    else:
        print("Plots adn HTML report created into directory '%s'. You can just move it to some www area and access the pages via web browser" % (
            ','.join(opts.outputDir)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create set of TICL-barrel validation plots from one or more DQM files.")
    parser.add_argument("files", metavar="file", type=str, nargs="+",
                        default="DQM_V0001_R000000001__Global_CMSSW_X_Y_Z__RECO.root",
                        help="DQM file to plot the validation plots from")
    parser.add_argument("-o", "--outputDir", type=str, default=["plots1", "plots2"], nargs="+",
                        help="Plot output directories (default: 'plots1'")
    parser.add_argument("--subdirprefix", type=str, default=["plots1", "plots2"], nargs="+",
                        help="Prefix for subdirectories inside outputDir (default:'plots1')")
    parser.add_argument("--no-ratio", action="store_true", default=False,
                        help="Disable ratio pads")
    parser.add_argument("--separate", action="store_true", default=False,
                        help="Save all plots separately instead of grouping them")
    parser.add_argument("--png", action="store_true",
                        help="Save plots in PNG instead of PDF")
    parser.add_argument("--no-html", action="store_true", default=False,
                        help="Disable HTML page generation")
    parser.add_argument("--html-sample", default="Sample",
                        help="Sample name for HTML page generation (deafult 'Sample')")
    parser.add_argument("--html-validation-name", type=str, default=["", ""], nargs="+",
                        help="Validation name for HTML page generation (enters to <title> element) (default '')")
    parser.add_argument("--collection", choices=collection_choices, default=layerClustersLabel,
                        help="Choose output plots collections among possible choices")
    parser.add_argument("--extended", action="store_true", default=False,
                        help="Include extended set of plots (e.g. bunch of distributions; default off)")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Be verbose")

    opts = parser.parse_args()

    for f in opts.files:
        if not os.path.exists(f):
            parser.error("DQM file %s does not exists" % f)

    main(opts)
