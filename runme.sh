#!/bin/sh -e

#converter/convert.py
## for the combine tool!
#export PATH=/home/walten/git/HiggsAnalysis/CombinedLimit/build/bin/:/home/walten/git/HiggsAnalysis/CombinedLimit/test:$PATH
#export LD_LIBRARY_PATH=/home/walten/git/HiggsAnalysis/CombinedLimit/build/lib/:$LD_LIBRARY_PATH
# export PATH=/home/walten/git/HiggsAnalysis/CombinedLimit/scripts/:$PATH
#export PYTHONPATH=/home/walten/git/HiggsAnalysis/CombinedLimit/build/lib/python:$PYTHONPATH

# python3 ./HiggsAnalysis/CombinedLimit/test/datacardConvert.py /home/walten/git/combine2pyhf//validation/cards/combine/combine2pyhf/sus20004/2DTChiHH600_LSP150.txt --bbl --normshape --prune --out /home/walten/git/combine2pyhf//validation/cards/combine/combine2pyhf/sus20004/2DTChiHH600_LSP150
# ./HiggsAnalysis/CombinedLimit/test/datacardConvert.py ~/Downloads/20_004/datacards/datacards/TChiHH/2DTChiHH800_LSP300.txt  --out ./sus20004

# datacardConvert.py ~/Downloads/20_004/datacards/datacards/TChiHH/2DTChiHH800_LSP300.txt  --out ./sus20004
datacardConvert.py ./cards/combine/one-bin/one-bin-sys-histosys-corr.txt --out ./onebin
