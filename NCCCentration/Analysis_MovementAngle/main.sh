#!/usr/bin/env bash

if [ "$1" = "clean" ]; then
  pushd src/SimulationCode
  make clean
  popd
  exit 0
fi

pushd src/SimulationCode
make &&
./mtsim -m 0 &&
cp out_* ../../Result/Simulation/MTFixed/ &&
./mtsim -m 1 &&
cp out_* ../../Result/Simulation/MTVariable/
popd
pwd

pushd src/AnalysisCode/
python main.py &&
Rscript wilcoxonTest.R
popd
