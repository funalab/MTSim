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

docker$ cd /root/mtsim
docker$ cd NCCCentration/MovementAngle/
docker$ ./main.sh
```
