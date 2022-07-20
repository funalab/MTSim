## How to build
```sh
% docker build -t funasoul/mtsim .
```

## How to run
### If you are using macOS
1. Start XQuartz
2. Open Preferences
3. Go to Security Settings and ensure that ["Allow connections from network clients" is on](https://gist.github.com/sorny/969fe55d85c9b0035b0109a31cbcb088)

```sh
% cd ..
% docker run --rm -it -e DISPLAY="$(hostname):0" -v ~/.Xauthority:/root/.Xauthority -v $PWD:/root/mtsim funasoul/mtsim

# Compares the steady state between MT Fixed and MT Variable for the translation and rotation of the central body
docker$ cd /root/mtsim/NCCCentration/MovementAngle/
docker$ ./main.sh

# Perform steady state analysis for a wide range of initial positions and rotation angles for the centrosome
docker$ cd /root/mtsim/NCCCentration/SteadyState/
docker$ ./main.sh

# Simulates the movement of the nuclus-centrosome complex in spindle elongation based on MT Fixed and MT Variable and compares their steady states
docker$ cd /root/mtsim/SpindleElongation/
docker$ ./main.sh
```
