import numpy as np
import matplotlib.pyplot as plt

def parse_stretch_hop_metrics(multi_motif_metrics_file, exclude_values):
    stretch_list = []
    hop_list = []
    m1_list = []
    with open(multi_motif_metrics_file) as fp:
        for line in fp:
            line = line.strip().split(",")
            if float(line[1]) not in exclude_values:
                stretch_list.append(float(line[1]))
            if float(line[2]) not in exclude_values:
                hop_list.append(float(line[2]))
            if float(line[3]) not in exclude_values:
                m1_list.append(float(line[3]))
    return stretch_list, hop_list, m1_list

def parse_grid(grid_metrics_file):
    stretch_list = []
    hop_list = []
    m1_list = []
    with open(grid_metrics_file) as fp:
        for line in fp:
            line = line.strip().split(",")
            stretch_list.append(float(line[0]))
            hop_list.append(float(line[1]))
            m1_list.append(float(line[2]))
    return stretch_list, hop_list, m1_list

def reproduce_figure_5(rrg_stretch_list, rrg_hop_list, grid_stretch_list, grid_hop_list, save_plot_path):

    grid_stretch_data_sorted = np.sort(grid_stretch_list)
    grid_stretch_p = 1. * np.arange(len(grid_stretch_list)) / (len(grid_stretch_list) - 1)
    plt.figure(5)
    plt.plot(grid_stretch_data_sorted, grid_stretch_p, label = "+Grid stretch")

    grid_hop_data_sorted = np.sort(grid_hop_list)
    grid_hop_p = 1. * np.arange(len(grid_hop_list)) / (len(grid_hop_list) - 1)
    plt.figure(5)
    plt.plot(grid_hop_data_sorted, grid_hop_p, label = "+Grid hop count")

    rrg_stretch_data_sorted = np.sort(rrg_stretch_list)
    rrg_stretch_p = 1. * np.arange(len(rrg_stretch_list)) / (len(rrg_stretch_list) - 1)
    plt.figure(5)
    plt.plot(rrg_stretch_data_sorted, rrg_stretch_p, label = "RRG stretch")

    rrg_hop_data_sorted = np.sort(rrg_hop_list)
    rrg_hop_p = 1. * np.arange(len(rrg_hop_list)) / (len(rrg_hop_list) - 1)
    plt.figure(5)
    plt.plot(rrg_hop_data_sorted, rrg_hop_p, label = "RRG hop count")


    plt.legend()
    plt.xlabel("City-city stretch or hop count")
    plt.ylabel("CDF across city pairs")
    plt.savefig(save_plot_path)

def reproduce_figure_10(grid_m1_list, best_m1_list, save_plot_path):

    grid_m1_data_sorted = np.sort(grid_m1_list)
    grid_m1_p = 1. * np.arange(len(grid_m1_list)) / (len(grid_m1_list) - 1)
    plt.figure(10)
    plt.plot(grid_m1_data_sorted, grid_m1_p, label="+Grid 53")

    best_m1_data_sorted = np.sort(best_m1_list)
    best_m1_p = 1. * np.arange(len(best_m1_list)) / (len(best_m1_list) - 1)
    plt.figure(10)
    plt.plot(best_m1_data_sorted, best_m1_p, label="Best motif 53")

    plt.legend()
    plt.xlabel("City-city M1")
    plt.ylabel("CDF across city pairs")
    plt.savefig(save_plot_path)



if __name__ == '__main__':

    motif_level0_40_40_53deg_5014_stretch_list, motif_level0_40_40_53deg_5014_hop_list, motif_level0_40_40_53deg_5014_m1_list = \
        parse_stretch_hop_metrics(
            "../output_data/multi_motif_40_40_53deg_5014/level_0_motif_metrics.txt",
            [99999.0]
        )

    motif_level1_40_40_53deg_5014_stretch_list, motif_level1_40_40_53deg_5014_hop_list, motif_level1_40_40_53deg_5014_m1_list = \
        parse_stretch_hop_metrics(
            "../output_data/multi_motif_40_40_53deg_5014/level_1_motif_metrics.txt",
            [99999.0]
        )

    motif_level2_40_40_53deg_5014_stretch_list, motif_level2_40_40_53deg_5014_hop_list, motif_level2_40_40_53deg_5014_m1_list = \
        parse_stretch_hop_metrics(
            "../output_data/multi_motif_40_40_53deg_5014/level_2_motif_metrics.txt",
            [99999.0]
        )

    motif_40_40_53deg_5014_stretch_list = motif_level0_40_40_53deg_5014_stretch_list + \
                                          motif_level1_40_40_53deg_5014_stretch_list + \
                                          motif_level2_40_40_53deg_5014_stretch_list

    motif_40_40_53deg_5014_hop_list = motif_level0_40_40_53deg_5014_hop_list + \
                                      motif_level1_40_40_53deg_5014_hop_list + \
                                      motif_level2_40_40_53deg_5014_hop_list

    motif_40_40_53deg_5014_m1_list = motif_level0_40_40_53deg_5014_m1_list + \
                                     motif_level1_40_40_53deg_5014_m1_list + \
                                     motif_level2_40_40_53deg_5014_m1_list

    grid_40_40_53deg_1467_stretch_list, \
    grid_40_40_53deg_1467_hop_list,  \
    grid_40_40_53deg_1467_m1_list = parse_grid("results/Grid_metrics.txt")

    reproduce_figure_5(motif_40_40_53deg_5014_stretch_list,
                       motif_40_40_53deg_5014_hop_list,
                       grid_40_40_53deg_1467_stretch_list,
                       grid_40_40_53deg_1467_hop_list,
                       "plots/figure5.png"
                       )

    reproduce_figure_10(grid_40_40_53deg_1467_m1_list,
                        motif_40_40_53deg_5014_m1_list,
                        "plots/figure10.png"
                        )