#!/bin/sh

#converter/convert.py
PYTHONPATH=/home/walten/git/HiggsAnalysis/CombinedLimit/build/lib/python:$PYTHONPATH

python3 ./HiggsAnalysis/CombinedLimit/test/datacardConvert.py /home/walten/git/combine2pyhf//validation/cards/combine/combine2pyhf/sus20004/2DTChiHH600_LSP150.txt --bbl --normshape --prune --out /home/walten/git/combine2pyhf//validation/cards/combine/combine2pyhf/sus20004/2DTChiHH600_LSP150
