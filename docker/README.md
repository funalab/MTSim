## How to build
```sh
% docker build -t funasoul/mtsim .
```

## How to run
```sh
% cd ..
% docker run --rm -it -e DISPLAY="$(hostname):0" -v ~/.Xauthority:/root/.Xauthority -v $PWD:/root/mtsim funasoul/mtsim
```
