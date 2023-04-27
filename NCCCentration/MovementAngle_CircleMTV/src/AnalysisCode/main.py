#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# FileName: 	main
# CreatedDate:  2019-06-14 19:20:26 +0900
# LastModified: 2019-10-05 16:42:26 +0900
#

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from barplot_annotate_brackets import barplot_annotate_brackets


class draw_graph_ncc_centering:
    def __init__(self):
        # variable
        base_dir = Path("../../Result/",
                        "Simulation")
        variable_path = base_dir / 'out_terminal.csv'

        self.output_pos = Path.cwd() / '../../Result/Analysis/position.pdf'
        self.output_ang = Path.cwd() / '../../Result/Analysis/angle.pdf'

        self.pos_str = 'Nuc_distance[per]'
        self.ang_str = 'angel[degree]'

        self.pos_ylab = 'Migration length of\n nucleus-centrosome complex / cell size [\%]'
        self.ang_ylab = 'Angle of\n nucleus-centrosome complex rotation $\phi$ [deg]'

        self.label = ['The MT-Variable\nmodel with Circle Cell']
        self.width = 0.4

        # read
        self.variable_model = pd.read_csv(variable_path)

        self.pos_val = stats.ttest_ind(self.variable_model[self.pos_str])
        self.ang_val = stats.ttest_ind(self.variable_model[self.ang_str])

        # setting matplotlib
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.rc('font', size=18)
        plt.rc('xtick', direction='in')
        plt.rc('ytick', direction='in')

    def draw_graph(self, mode: str):
        if mode == 'pos':
            name = self.output_pos
            str = self.pos_str
            val = self.pos_val
            lab = self.pos_ylab
        else:
            name = self.output_ang
            str = self.ang_str
            val = self.ang_val
            lab = self.ang_ylab

        heights = [self.variable_model[str].max()]
        pos = np.arange(0, len(heights)) / 2

        plt.clf()
        plt.boxplot([self.variable_model[str]],
                    widths=self.width,
                    positions=pos,
                    labels=self.label)
        plt.xlim(pos[0] - self.width * 0.7, pos[1] + self.width * 0.7)
        plt.ylabel(lab, fontsize = 14)
        # plt.ylim(0, )

        plt.tight_layout()
        plt.savefig(name)


if __name__ == "__main__":
    c = draw_graph_ncc_centering()
    c.draw_graph('pos')
    c.draw_graph('ang')
