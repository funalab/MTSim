cd ./src/SimulationCode
make CRADS="15*pow\(10,-6\)"
./mtsim
cp vectordata/centering_vector_strain8_rad25_rads15.dat ../../Result/Simulation
make clean

make CRADS="25*pow\(10,-6\)"
./mtsim
cp vectordata/centering_vector_strain8_rad25_rads25.dat ../../Result/Simulation
make clean

make CRADS="10*pow\(10,-6\)"
./mtsim
cp vectordata/centering_vector_strain8_rad25_rads10.dat ../../Result/Simulation
make clean

cd ../../
