#!/usr/bin/env bash

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

if [ "$1" = "clean" ]; then
  pushd src/SimulationCode
  make clean
  popd
  exit 0
fi

pushd src/SimulationCode
make &&
./mtsim -m 0 &&
cp out_* ../../Result/Simulation/led/MTFixed/ &&
./mtsim -m 0 -s&&
cp out_* ../../Result/Simulation/sled/MTFixed/ &&
./mtsim -m 1 &&
cp out_* ../../Result/Simulation/MTVariable/
popd
pwd

pushd src/AnalysisCode/
python main.py &&
python main_sq.py &&
Rscript wilcoxonTest.R
popd
