#!/usr/bin/env bash

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

if [ "$1" = "clean" ]; then
  pushd ./src/SimulationCode/celegans/
  make clean
  popd
  pushd ./src/SimulationCode/sea-urchin/
  make clean
  popd
  exit 0
fi

pushd ./src/SimulationCode/celegans/
./because_of_graph.sh -m 0 -l 0 &&
mv result.csv ../../../Result/Simulation/Cel_MTFixed.csv &&
./because_of_graph.sh -m 1 -l 0 &&
mv result.csv ../../../Result/Simulation/Cel_MTVariable.csv
./because_of_graph_sq.sh -m 1 -l 0 -s &&
mv result.csv ../../../Result/Simulation/Cel_led_MTVariable.csv
popd

pushd ./src/SimulationCode/sea-urchin/
LF=$(printf '\\\012_')
LF=${LF%_}
make &&
./mtsim -m 0 < enter.txt > result.csv &&

sed -i -e "/NN/d;/aspect/d;/to/d;/^3.[1-4]/d" result.csv &&
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv &&
mv result.csv ../../../Result/Simulation/Su_MTFixed.csv &&

./mtsim -m 1 < enter.txt > result.csv &&

sed -i -e "/NN/d;/aspect/d;/to/d;/^3.[1-4]/d" result.csv &&
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv &&
mv result.csv ../../../Result/Simulation/Su_MTVariable.csv
make clean
popd

pushd src/AnalysisCode
python analysis_cel.py -s -e &&
python analysis_cel_sq.py -s &&
python analysis_su.py -e
popd
