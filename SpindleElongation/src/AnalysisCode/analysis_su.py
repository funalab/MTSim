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
import numpy as np
import pandas as pd


def main():
    # setting option
    parser = argparse.ArgumentParser(
        description="draw graph for simulation vs. experiment"
    )
    parser.add_argument(
        "-e",
        "--experiment",
        action="store_true",
        help="draw experimental data",
        default=False,
    )
    parser.add_argument(
        "-b",
        "--before-graph",
        action="store_true",
        help="draw before graph",
        default=False,
    )
    parser.add_argument(
        "-el",
        "--ellipsoid",
        action="store_true",
        help="draw ellipsoid data",
        default=False,
    )
    args = parser.parse_args()

    # check option
    data_type = "_cy_fit"
    if args.ellipsoid:
        data_type = "_el"
    x_name = "#aspect_ratio"
    y_name = "spindle_length/(Rad*2)"
    y_lim_size = 1
    y_lim_name = "Elongated pole-to-pole distance/Cell size"
    if args.before_graph:
        y_name = "spindle_length"
        y_lim_size = 200
        y_lim_name = r"Elongated pole-to-pole distance[${\rm \mu m}$]"

    # path
    data_path = os.path.join(os.path.abspath("../../"))
    fixed_path = os.path.join(data_path, "Result", "Simulation", "Su_MTFixed.csv")
    variable_path = os.path.join(data_path, "Result", "Simulation", "Su_MTVariable.csv")
    exp_path = os.path.join(data_path, "Experiment", "exp_su.csv")

    # read data
    fixed = pd.read_csv(fixed_path)
    variable = pd.read_csv(variable_path)
    experiment = pd.read_csv(exp_path)

    # draw graph
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")
    plt.rc("xtick", direction="in")
    plt.rc("ytick", direction="in")
    plt.plot(
        fixed[x_name], fixed[y_name], lw=5, zorder=1, c="blue", label="The MT-Fixed Model"
    )
    plt.plot(
        variable[x_name],
        variable[y_name],
        zorder=2,
        lw=5,
        c="green",
        label="The MT-Variable Model",
    )
    if args.experiment:
        plt.scatter(
            experiment[x_name],
            experiment[y_name],
            zorder=6,
            s=80,
            facecolor="red",
            edgecolors="black",
            label="Experimental data",
        )

    # settig graph
    n = len(experiment)
    plt.xlabel(r"Aspect ratio", fontsize=18)
    plt.xlim(1.0, 4.0)
    plt.ylabel(y_lim_name, fontsize=18)
    plt.ylim(0, y_lim_size)
    if args.before_graph:
        plt.yticks(np.arange(0, 250, 50))
    plt.legend(loc="upper left", ncol=2, scatterpoints=1, fontsize=12)
    plt.tick_params(labelsize=16, pad=8)
    plt.tight_layout()
    plt.savefig("../../Result/Analysis/Su_sim_vs_exp.pdf")


if __name__ == "__main__":
    main()
