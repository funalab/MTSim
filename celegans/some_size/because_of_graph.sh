#!/bin/zsh

# set variable
local -A opthash
zparseopts -D -A opthash -- h m l:
model=0
len=0
LF=$(printf '\\\012_')
LF=${LF%_}

if [[ -n "${opthash[(i)-h]}" ]]; then
    echo "Usage : ./mtsim [option]"
    echo " -m # : specify model (ex. -m 1 )"
    echo "     0: MT angle fixed"
    echo "     1: MT angle variable"
    echo " -l # : metaspindle length (ex. -l 1 )"
    echo "     0: Rad vs metaspindle length"
    echo "     1: aspect ratio vs metaspindle length / Rad"
    return
fi

if [[ -n "${opthash[(i)-m]}" ]]; then
    if [[ ${opthash[-m]} -eq 1 ]]; then
        model=1
    else
        model=0
    fi
fi

if [[ -n "${opthash[(i)-l]}" ]]; then
    if [[ ${opthash[-l]} -eq 1 ]]; then
        len=1
    else
        len=0
    fi
fi

# show simulation environment
echo "1.0 < aspect ratio < 3.0"
echo "step = 0.1"
echo ""

if [[ ${model} -eq 0 ]]; then
    echo "MT_fixed_model"
elif [[ ${model} -eq 1 ]]; then
    echo "MT_variable_model"
else
    echo "model error"
    return
fi

if [[ ${len} -eq 0 ]]; then
    echo "Rad vs metaspindle length"
elif [[ ${len} -eq 1 ]]; then
    echo "aspect ratio vs metaspindle length / Rad"
else
    echo "length error"
    return
fi

# exec simulation
sed -i -e "s/ar_step = 0.5/ar_step = 0.1/g" main.c

make -s
./mtsim -m ${model} -l 0 < enter.txt > result.csv

sed -i -e "s/ar_init = 1.0/ar_init = 1.5/g" main.c
make -s
./mtsim -m ${model} -l 0 < enter.txt >> result.csv

sed -i -e "s/ar_init = 1.5/ar_init = 2.0/g" main.c
make -s
./mtsim -m ${model} -l 0 < enter.txt >> result.csv

sed -i -e "s/ar_init = 2.0/ar_init = 2.5/g" main.c
make -s
./mtsim -m ${model} -l 0 < enter.txt >> result.csv

sed -i -e "s/ar_init = 2.5/ar_init = 3.0/g" main.c
make -s
./mtsim -m ${model} -l 0 < enter.txt >> result.csv

# arrange result
sed -i -e "/NN/d;/aspect/d;/to/d;/3.[1-4]/d" result.csv
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

# clear phase
echo ""
rm main.c-e result.csv-e
make clean
git reset --hard master
