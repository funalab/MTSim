#!/usr/bin/env sh
cd src/SimulationCode
make
./mtsim -m 0
cp out_* ../../Result/Simulation/MTFixed/
./mtsim -m 1
cp out_* ../../Result/Simulation/MTVariable/
cd ../..

cd src/AnalysisCode/
python main.py
Rscript wilcoxonTest.R
cd ../../
