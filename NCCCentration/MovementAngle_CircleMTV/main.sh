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
./mtsim -m 1 &&
cp out_* ../../Result/Simulation/ &&
popd
pwd

pushd src/AnalysisCode/
python main.py &&
popd
