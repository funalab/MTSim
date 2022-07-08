#!/bin/zsh

# set variable
local -A opthash
zparseopts -D -A opthash -- h m: l: p: a:
model=0
len=0
para=4
ar=1.0
LF=$(printf '\\\012_')
LF=${LF%_}

if [[ -n "${opthash[(i)-h]}" ]]; then
    echo "Usage : ./mtsim [option]"
    echo " -m # : specify model (ex. -m 1 )"
    echo "     0: MT angle fixed"
    echo "     1: MT angle variable"
    return
fi

if [[ -n "${opthash[(i)-m]}" ]]; then
    if [[ ${opthash[-m]} -eq 1 ]]; then
        model=1
    else
        model=0
    fi
fi

# simulation
step=2
step_num=5
bef_angle=130
aft_angle=$((${bef_angle}+${step}))
last_angle=$((${bef_angle}+${step}*${step_num}))
echo "aspect ratio = ${ar}"
echo "${bef_angle} < MT init angle < ${last_angle}"
echo ""

if [[ ${model} -eq 0 ]]; then
    echo "MT_fixed_model"
    angle_def=137
elif [[ ${model} -eq 1 ]]; then
    echo "MT_variable_model"
    angle_def=120
else
    echo "model error"
    return
fi

# exec simulation
sed -i -e "s/ar_init = 1.0/ar_init = ${ar}/g;s/\(p<\)6/\11/g" main.c

sed -i -e "s/MTInitAngle_degree = ${angle_def}/MTInitAngle_degree = ${bef_angle}/g" main.c
make -s
./mtsim -m ${model} < enter.txt > result.csv

for i in {1..5}; do
    sed -i -e "s/MTInitAngle_degree = ${bef_angle}/MTInitAngle_degree = ${aft_angle}/g" main.c
    make -s
    ./mtsim -m ${model} < enter.txt >> result.csv

    bef_angle=$((${bef_angle}+${step}))
    aft_angle=$((${aft_angle}+${step}))
done

# arrange result
sed -i -e "/NN/d;/aspect/d;/to/d" result.csv
sed -i '' -ne "/^${ar}/p" result.csv
sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

# clear phase
echo ""
rm *-e
make clean
git reset --hard master
