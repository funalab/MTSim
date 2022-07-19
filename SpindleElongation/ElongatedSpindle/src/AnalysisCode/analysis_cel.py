#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# FileName: 	main
# CreatedDate:  2017-11-20 14:05:59
#

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd


def main():
    # setting option
    parser = argparse.ArgumentParser(
        description="draw graph for simulation vs. experiment"
    )
    parser.add_argument(
        "-s",
        "--simulation",
        action="store_true",
        help="draw simulation data",
        default=False,
    )
    parser.add_argument(
        "-e",
        "--experiment",
        action="store_true",
        help="draw experiment data",
        default=False,
    )
    args = parser.parse_args()

    # path
    data_path = os.path.join("../../../")
    fixed_path = os.path.join(
        data_path, "ElongatedSpindle/Result/Simulation/Cel_MTFixed.csv"
    )
    variable_path = os.path.join(
        data_path, "ElongatedSpindle/Result/Simulation/Cel_MTVariable.csv"
    )
    exp_path = os.path.join(data_path, "Experiment_ElongatedSpindle/exp_cel.csv")

    # read data
    fixed = pd.read_csv(fixed_path)
    variable = pd.read_csv(variable_path)
    experiment = pd.read_csv(exp_path)

    # draw graph
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")
    plt.rc("xtick", direction="in")
    plt.rc("ytick", direction="in")
    if args.simulation:
        plt.plot(
            fixed["#aspect_ratio"],
            fixed["spindle_length/(Rad*2)"],
            lw=5,
            zorder=1,
            c="blue",
            label="MT Fixed Model",
        )
        plt.plot(
            variable["#aspect_ratio"],
            variable["spindle_length/(Rad*2)"],
            zorder=2,
            lw=5,
            c="green",
            label="MT Variable Model",
        )
    if args.experiment:
        exp_data1 = experiment[
            experiment["name"].str.contains("TH27_")
            | experiment["name"].str.contains("TH32_17")
        ]
        exp_data2 = experiment[
            experiment["name"].str.contains("dpy-11")
            | experiment["name"].str.contains("DPY11")
        ]
        exp_data3 = experiment[experiment["name"].str.contains("C27D9.1")]
        exp_data4 = experiment[
            experiment["name"].str.contains("CAL1628")
            & ~experiment["name"].str.contains("C27D9.1")
        ]
        plt.scatter(
            exp_data1["aspect ratio"],
            exp_data1["distance from centrosome/long axis"],
            zorder=6,
            s=80,
            facecolor="red",
            edgecolors="black",
            label="TH27,TH32",
        )
        plt.scatter(
            exp_data2["aspect ratio"],
            exp_data2["distance from centrosome/long axis"],
            zorder=4,
            marker="^",
            s=80,
            facecolor="magenta",
            edgecolors="black",
            label=r"TH27,TH32;\textit{dpy-11}(RNAi)",
        )
        plt.scatter(
            exp_data4["aspect ratio"],
            exp_data4["distance from centrosome/long axis"],
            zorder=3,
            marker="*",
            s=120,
            facecolor="yellow",
            edgecolors="black",
            label="CAL1628",
        )
        plt.scatter(
            exp_data3["aspect ratio"],
            exp_data3["distance from centrosome/long axis"],
            zorder=5,
            marker=",",
            s=80,
            facecolor="cyan",
            edgecolors="black",
            label=r"CAL1628;\textit{C27D9.1}(RNAi)",
        )

    # settig graph
    plt.xlabel(r"Aspect ratio", fontsize=18)
    plt.xlim(1, 3)
    plt.ylabel("Elongated pole-to-pole distance/Cell size", fontsize=18)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", ncol=2, scatterpoints=1, fontsize=11)
    plt.tick_params(labelsize=16, pad=8)
    plt.tight_layout()
    plt.savefig("../../Result/Analysis/Cel_sim_vs_exp.pdf")


if __name__ == "__main__":
    main()
