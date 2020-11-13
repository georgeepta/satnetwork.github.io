import numpy as np
import matplotlib.pyplot as plt

def parse_stretch_hop_metrics(multi_motif_metrics_file, exclude_values):
    stretch_list = []
    hop_list = []
    with open(multi_motif_metrics_file) as fp:
        for line in fp:
            line = line.strip().split(",")
            if float(line[1]) not in exclude_values:
                stretch_list.append(float(line[1]))
            if float(line[2]) not in exclude_values:
                hop_list.append(float(line[2]))
    return stretch_list, hop_list

def parse_grid(grid_metrics_file):
    stretch_list = []
    hop_list = []
    with open(grid_metrics_file) as fp:
        for line in fp:
            line = line.strip().split(",")
            stretch_list.append(float(line[0]))
            hop_list.append(float(line[1]))
    return stretch_list, hop_list

def reproduce_figure_5(rrg_stretch_list, rrg_hop_list, grid_stretch_list, grid_hop_list, save_plot_path):

    grid_stretch_data_sorted = np.sort(grid_stretch_list)
    grid_stretch_p = 1. * np.arange(len(grid_stretch_list)) / (len(grid_stretch_list) - 1)
    plt.plot(grid_stretch_data_sorted, grid_stretch_p, label = "+Grid stretch")

    grid_hop_data_sorted = np.sort(grid_hop_list)
    grid_hop_p = 1. * np.arange(len(grid_hop_list)) / (len(grid_hop_list) - 1)
    plt.plot(grid_hop_data_sorted, grid_hop_p, label = "+Grid hop count")

    rrg_stretch_data_sorted = np.sort(rrg_stretch_list)
    rrg_stretch_p = 1. * np.arange(len(rrg_stretch_list)) / (len(rrg_stretch_list) - 1)
    plt.plot(rrg_stretch_data_sorted, rrg_stretch_p, label = "RRG stretch")

    rrg_hop_data_sorted = np.sort(rrg_hop_list)
    rrg_hop_p = 1. * np.arange(len(rrg_hop_list)) / (len(rrg_hop_list) - 1)
    plt.plot(rrg_hop_data_sorted, rrg_hop_p, label = "RRG hop count")


    plt.legend()
    plt.xlabel("City-city stretch or hop count")
    plt.ylabel("CDF across city pairs")
    plt.savefig(save_plot_path)

if __name__ == '__main__':

    motif_level0_40_40_53deg_5014_stretch_list, motif_level0_40_40_53deg_5014_hop_list = \
        parse_stretch_hop_metrics(
            "../output_data/multi_motif_40_40_53deg_5014/level_0_motif_metrics.txt",
            [99999.0]
        )

    grid_40_40_53deg_5014_stretch_list, grid_40_40_53deg_5014_hop_list = \
        parse_grid("results/Grid_metrics.txt")

    reproduce_figure_5(motif_level0_40_40_53deg_5014_stretch_list,
                       motif_level0_40_40_53deg_5014_hop_list,
                       grid_40_40_53deg_5014_stretch_list,
                       grid_40_40_53deg_5014_hop_list,
                       "plots/figure5.png"
                       )