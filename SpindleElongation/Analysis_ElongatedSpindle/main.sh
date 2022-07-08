cd ./src/SimulationCode/celegans/
./because_of_graph.sh -m 0 -l 0
mv result.csv ../../../Result/Simulation/Cel_MTFixed.csv
./because_of_graph.sh -m 1 -l 0
mv result.csv ../../../Result/Simulation/Cel_MTVariable.csv

cd ../sea-urchin/
LF=$(printf '\\\012_')
LF=${LF%_}
make
./mtsim -m 0 < enter.txt > result.csv

sed -i -e "/NN/d;/aspect/d;/to/d;/^3.[1-4]/d" result.csv
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv
mv result.csv ../../../Result/Simulation/Su_MTFixed.csv
rm result.csv-e

./mtsim -m 1 < enter.txt > result.csv

sed -i -e "/NN/d;/aspect/d;/to/d;/^3.[1-4]/d" result.csv
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv
mv result.csv ../../../Result/Simulation/Su_MTVariable.csv
rm result.csv-e

make clean
cd ../../../

cd src/AnalysisCode
python analysis_cel.py -s -e
python analysis_su.py -e
cd ../../
