cd src/SimulationCode
make
./mtsim -m 0
cp out_* ../../Result/Simulation/MTFixed/
./mtsim -m 1
cp out_* ../../Result/Simulation/MTVariable/
cd ../..