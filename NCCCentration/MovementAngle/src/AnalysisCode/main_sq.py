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
        led_path = base_dir / 'MTFixed/led/out_terminal.csv'
        sled_path = base_dir / 'MTFixed/sled/out_terminal.csv'

        self.output_pos = Path.cwd() / '../../Result/Analysis/PositionAngle/position_sled.pdf'
        self.output_ang = Path.cwd() / '../../Result/Analysis/PositionAngle/angle_sled.pdf'

        self.pos_str = 'Nuc_distance[per]'
        self.ang_str = 'angel[degree]'

        self.pos_ylab = 'Migration length of\n nucleus-centrosome complex / cell size [\%]'
        self.ang_ylab = 'Angle of\n nucleus-centrosome complex rotation $\phi$ [deg]'

        self.label = ['The MT-Fixed\n length', 'The MT-Fixed\n squared length']
        self.width = 0.4

        # read
        self.led_model = pd.read_csv(led_path)
        self.sled_model = pd.read_csv(sled_path)

        self.pos_val = stats.ttest_ind(self.led_model[self.pos_str],
                                       self.sled_model[self.pos_str])
        self.ang_val = stats.ttest_ind(self.led_model[self.ang_str],
                                       self.sled_model[self.ang_str])

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

        heights = [self.led_model[str].max(), self.sled_model[str].max()]
        pos = np.arange(0, len(heights)) / 2

        plt.clf()
        plt.boxplot([self.led_model[str], self.sled_model[str]],
                    widths=self.width,
                    positions=pos,
                    labels=self.label)
        plt.xlim(pos[0] - self.width * 0.7, pos[1] + self.width * 0.7)
        plt.ylabel(lab, fontsize = 14)
        # plt.ylim(0, )
        # barplot_annotate_brackets(0, 1, val[1], pos, heights)

        plt.tight_layout()
        plt.savefig(name)


if __name__ == "__main__":
    c = draw_graph_ncc_centering()
    c.draw_graph('pos')
    c.draw_graph('ang')
