#!/bin/zsh

rm out_f_MTP.csv
# set variable
local -A opthash
zparseopts -D -A opthash -- h m: l: p: a:
model=0
len=0
para=4
ar=1.6
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
    echo " -p # : effect of parameter fluctuations (ex. -p 1 )"
    echo "     0: volume vary"
    echo "     1: MT init angle vary"
    echo "     2: MT density vary"
    echo "     3: attraction coefficient vary"
    echo "     4: default"
    echo " -a # : aspect ratio (ex. -a 1 )"
    echo "     0: aspect ratio = 1.6"
    echo "     1: aspect ratio = 2.6"
    return
fi

if [[ -n "${opthash[(i)-p]}" ]]; then
    para=${opthash[-p]}
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

if [[ -n "${opthash[(i)-a]}" ]]; then
    if [[ ${opthash[-a]} -eq 1 ]]; then
        ar=2.6
    else
        ar=1.6
    fi
fi

# commentout XQuarts
sed -i -e "s/draw_graphs(/\/\/ draw_graphs(/g" main.c
sed -i -e "s/XFlush/\/\/ XFlush/g" main.c

# volume vary
if [[ ${para} -eq 0 ]]; then
    # show simulation environment
    echo "aspect ratio = ${ar}"
    echo "0.5 < volume < 4.0"
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
    
    # ipython lambda func for calculating Cir_Rad
    # refunc = lambda v : math.pow((3 * v / 4 / math.pi) , (1/3))
    # exec simulation
    sed -i -e "s/ar_init = 1.0/ar_init = ${ar}/g" main.c

    sed -i -e "s/19.487570/10.60784417947055/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt > result.csv
    
    sed -i -e "s/10.60784417947055/13.365046175719757/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/13.365046175719757/15.299158709729346/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/15.299158709729346/16.838903009606295/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/16.838903009606295/18.139158392989046/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv

    sed -i -e "s/18.139158392989046/19.27573210407049/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/19.27573210407049/20.292075912899804/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len}< enter.txt >> result.csv
    
    sed -i -e "s/20.292075912899804/21.2156883589411/g" mtsim.h
    make -s
    ./mtsim -m ${model} -l ${len}< enter.txt >> result.csv
    
    # arrange result
    sed -i -e "/NN/d;/aspect/d;/to/d" result.csv
    sed -i '' -ne "/^${ar}/p" result.csv
    sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

elif [[ ${para} -eq 1 ]]; then
    # show simulation environment
    step=10
    bef_angle=100
    aft_angle=${bef_angle}+${step}
    echo "aspect ratio = ${ar}"
    echo "${bef_angle} < MT init angle < 180"
    echo ""
    
    if [[ ${model} -eq 0 ]]; then
        echo "MT_fixed_model"
        angle_def=126
    elif [[ ${model} -eq 1 ]]; then
        echo "MT_variable_model"
        angle_def=97
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
    sed -i -e "s/ar_init = 1.0/ar_init = ${ar}/g" main.c

    sed -i -e "s/MTInitAngle_degree = ${angle_def}/MTInitAngle_degree = ${bef_angle}/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt > result.csv
    
    for i in {1..8}; do
        sed -i -e "s/MTInitAngle_degree = ${bef_angle}/MTInitAngle_degree = ${aft_angle}/g" main.c
        make -s
        ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv

        bef_angle=${bef_angle}+${step}
        aft_angle=${aft_angle}+${step}
    done
    
    # arrange result
    sed -i -e "/NN/d;/aspect/d;/to/d" result.csv
    sed -i '' -ne "/^${ar}/p" result.csv
    sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

elif [[ ${para} -eq 2 ]]; then
    # show simulation environment
    step=3
    n=9
    bef_dens=1
    aft_dens=${bef_dens}+${step}
    last=$((${bef_dens}+${step}*${n}))
    echo "aspect ratio = ${ar}"
    echo "${bef_dens} < MT density < ${last}"
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
    sed -i -e "s/ar_init = 1.0/ar_init = ${ar}/g" main.c

    sed -i -e "s/MTDiv90 = 6/MTDiv90 = ${bef_dens}/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt > result.csv
    
    for i in {1..${n}}; do
        sed -i -e "s/MTDiv90 = ${bef_dens}/MTDiv90 = ${aft_dens}/g" main.c
        make -s
        ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv

        bef_angle=${bef_dens}+${step}
        aft_angle=${aft_dens}+${step}
    done
    
    # arrange result
    sed -i -e "/NN/d;/aspect/d;/to/d" result.csv
    sed -i '' -ne "/^${ar}/p" result.csv
    sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

elif [[ ${para} -eq 3 ]]; then
    # show simulation environment
    step=1
    n=4
    bef_atCoe=-4
    aft_atCoe=${bef_atCoe}+${step}
    last=$((${bef_atCoe}+${step}*${n}))
    echo "aspect ratio = ${ar}"
    echo "10^${bef_atCoe} < attraction coefficient order < 10^${last}"
    echo ""
    
    if [[ ${model} -eq 0 ]]; then
        echo "MT_fixed_model"
        atCoe_def=0.25
    elif [[ ${model} -eq 1 ]]; then
        echo "MT_variable_model"
        atCoe_def=0.09
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
    sed -i -e "s/ar_init = 1.0/ar_init = ${ar}/g" main.c

    sed -i -e "s/ForceCoef2 = ${atCoe_def} * pow(10,-2)/ForceCoef2 = ${atCoe_def} * pow(10,${bef_atCoe})/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt > result.csv
    
    for i in {1..${n}}; do
        sed -i -e "s/ForceCoef2 = ${atCoe_def} * pow(10,${bef_atCoe})/ForceCoef2 = ${atCoe_def} * pow(10,${aft_atCoe})/g" main.c
        make -s
        ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv

        bef_atCoe=${bef_atCoe}+${step}
        aft_atCoe=${aft_atCoe}+${step}
    done
    
    # arrange result
    sed -i -e "/NN/d;/aspect/d;/to/d" result.csv
    sed -i '' -ne "/^${ar}/p" result.csv
    sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv

else
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
    ./mtsim -m ${model} -l ${len} < enter.txt > result.csv
    
    sed -i -e "s/ar_init = 1.0/ar_init = 1.5/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/ar_init = 1.5/ar_init = 2.0/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len} < enter.txt >> result.csv
    
    sed -i -e "s/ar_init = 2.0/ar_init = 2.5/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len}< enter.txt >> result.csv
    
    sed -i -e "s/ar_init = 2.5/ar_init = 3.0/g" main.c
    make -s
    ./mtsim -m ${model} -l ${len}< enter.txt >> result.csv
    
    # arrange result
    sed -i -e "/NN/d;/aspect/d;/to/d;/^3.[1-4]/d" result.csv
    sed -i -e '1s/^/#aspect_ratio,Rad,RadS,MetaSpindle_L,spindle_length,time[sec],spindle_length\/(Rad*2)'"$LF"'/' result.csv
fi

# clear phase
echo ""
rm *-e
make clean
git reset --hard master
