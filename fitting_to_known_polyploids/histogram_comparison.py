import math
import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import config
import ks_parsers
from results_viewer import batch_histogrammer, curve_fitting

class MyTestCase(unittest.TestCase):

    MBE_dpi=350
    def test_coffee_histogram(self):

        coffee_num=20#coffee_num=17
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim42_coffee/Allo_Coffea{0}".format(coffee_num)
        specks_csv_file = "Allo_Coffea{0}_ML_rep0_Ks_by_GeneTree.csv".format(coffee_num)
        ksrates_csv_file="coffea.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        species_for_plot_title= 'Coffea arabica'
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        #bin_size=0.001
        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.015
        density = 1

        make_both_histograms(bin_size, color, specks_out_folder,
                             wgd_ks, max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path)

    def test_poplar_histogram(self):

        pop_num=23
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim42_Poplar/Allo_Poplar{0}".format(pop_num)
        specks_csv_file = "Allo_Poplar{0}_ML_rep0_Ks_by_GeneTree.csv".format(pop_num)
        ksrates_csv_file = "poplar.ks.tsv"

        splat = specks_csv_file.split("_")
        species_run_name = splat[0] + splat[1]
        species_for_plot_title= 'Populus trichocarpa'
        specks_full_path = os.path.join(specks_out_folder, specks_csv_file)
        real_full_path = os.path.join(ksrates_out_folder, ksrates_csv_file)

        bin_size = 0.002
        max_Ks = 0.5
        color = 'blue'
        wgd_ks=0.21
        density = 1

        make_both_histograms(bin_size, color, specks_out_folder, wgd_ks,
                             max_Ks, density, real_full_path,
                             species_run_name, species_for_plot_title, specks_full_path)
    def test_maize_histogram(self):

        maize_num="32"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim42_Maize/Allo{0}_Maize".format(maize_num)
        specks_csv_file = "Allo{0}_Maize_ML_rep0_Ks_by_GeneTree.csv".format(maize_num)
        ksrates_csv_file="mays.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        species_for_plot_title= 'Zea mays'
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.3
        color='blue'
        wgd_ks=0.13
        density = 1
        make_both_histograms(bin_size, color,  specks_out_folder, wgd_ks,
                         max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path)

        hist_by_type_of_paralog(bin_size, color, density, max_Ks, species_run_name, specks_full_path,
                                     specks_out_folder, wgd_ks)

    def test_Single_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d25y2024_h17m47s14/Allo_Maize/8_final_results"

        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d26y2024_h12m30s41/Auto_Poplar/8_final_results"
        #specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_maize"

        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d26y2024_h13m45s40/Auto_Olive/8_final_results"

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        #specks_csv_file="Allo_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        #specks_csv_file = "Auto_Poplar_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        specks_csv_file = "Auto_Poplar_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="mays.ks.tsv"
        ksrates_csv_file="poplar.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)

        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)
        real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.3
        WGD_ks=0.001
        color=config.allo_color
        density = None

        make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks,max_Ks, density,real_ks_results,
                                  species_run_name, specks_ks_results)

def hist_by_type_of_paralog(bin_size, color, density, max_Ks, species_run_name, specks_full_path,
                            specks_out_folder, wgd_ks):
    ks_by_paralog_type = ks_parsers.read_Ks_csv_by_dup_type(specks_full_path)
    out_png4 = os.path.join(specks_out_folder, "specks_" + species_run_name + "_layered.png")
    fig = plt.figure(figsize=(10, 10), dpi=MyTestCase.MBE_dpi)
    [n_WGD, bins] = make_simple_histogram(ks_by_paralog_type['WGD'], species_run_name, bin_size,
                                          color, wgd_ks,
                                          max_Ks, density, out_png4)
    n_WGD_SSD, bins, patches = plt.hist(ks_by_paralog_type['WGD'], bins=bins,
                                        facecolor='blue', alpha=0.5,
                                        label="WGD-WGD")
    n_WGD_SSD, bins, patches = plt.hist(ks_by_paralog_type['WGD-SSD'], bins=bins,
                                        facecolor='green', alpha=0.5,
                                        label="WGD-SSD")
    n_SSD_SSD, bins, patches = plt.hist(ks_by_paralog_type['SSD-SSD'], bins=bins,
                                        facecolor='yellow', alpha=0.5,
                                        label="SSD-SSD")

    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.".format(species_run_name))
    plt.savefig(out_png4)
    plt.clf()
    plt.close()
def make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks, max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path):

    specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)
    real_ks_results = ks_parsers.parse_external_ksfile(real_full_path)

    if not os.path.exists(hist_comparison_out_folder):
        os.makedirs(hist_comparison_out_folder)
    out_png1 = os.path.join(hist_comparison_out_folder, "specks_" + species_run_name + "_out.png")
    out_png2 = os.path.join(hist_comparison_out_folder, "specks_" + species_run_name + "_fit.png")
    specks_hist_data =make_simple_histogram(specks_ks_results, species_run_name, bin_size, color, WGD_ks,
                          max_Ks, density, out_png1)
    #fit_fxns_to_Ks(specks_ks_results, species_run_name,5000,400,
    #               max_Ks, out_png2)

    out_png1 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_out.png")
    out_png2 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_fit.png")
    real_hist_data = make_simple_histogram(real_ks_results, species_run_name, bin_size, color, WGD_ks,
                          max_Ks, density, out_png1)

    out_png3 = os.path.join(hist_comparison_out_folder, "overlay_" + species_run_name + "_fit.png")
    overlay_histogram(species_run_name, species_for_plot_title,
                      [specks_hist_data,real_hist_data], WGD_ks, out_png3)

    out_png4 = os.path.join(hist_comparison_out_folder, "differences_" + species_run_name + "_fit.png")
    overlay_differences_in_curves(species_for_plot_title, [specks_hist_data,real_hist_data],
                                  WGD_ks, out_png4)

def make_simple_histogram(Ks_results, species_name, bin_size, color,WGD_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    fig = plt.figure(figsize=(10, 10), dpi=MyTestCase.MBE_dpi)
    x = Ks_results
    # print(PAML_hist_out_file)
    label="hist for " + os.path.basename(out_png).replace("_out.png","")
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])


    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs=sum(n)
    num_after_wgd=sum([n[i] for i in range(0,len(n)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1} pairs of genes. ~{2} retained from WGD.".format(
        species_name,num_pairs,num_after_wgd))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [n,bins]


def overlay_histogram(species_name, species_for_plot_title, list_of_hist_data, WGD_ks,out_png):

    colors = ['blue','green']
    labels = ['SpecKS','truth']

    fig = plt.figure(dpi=MyTestCase.MBE_dpi)
    ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    ax2 = plt.subplot2grid((4, 1), (3, 0), rowspan=1)

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])/2
        xs=[b + i*width for b in bins[0:len(bins)-1]]
        ax1.bar(xs,n,width=width,
                color=colors[i], alpha=1, label=labels[i])
    #ax[0].axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    diffs = [(list_of_hist_data[0][0][j]-list_of_hist_data[1][0][j]) for j in range(0,len(list_of_hist_data[1][0]))]
    rmse=math.sqrt(sum([d*d for d in diffs])/len(diffs))
    ax2.bar(xs, diffs , width=width,
            color='orange', alpha=1, label="error\nrmse={0}".format(round(rmse,3)))


    print(species_for_plot_title)
    num_pairs = sum(n)
    num_after_wgd = sum([n[i] for i in range(0, len(n)) if bins[i] > WGD_ks])
    ax1.legend()
    ax2.legend()
    plt.xlabel("Ks")
    ax2.set(ylabel="Error")
    ax1.set(ylabel="Density")
    ax1.set(xlim=[0,0.1])
    ax2.set(xlim=[0,0.1])
    #ax1.set(title ='Ks histogram for ' +  species_for_plot_title)
    #ax1.set(title ='$\Gamma + \mathit{\Gamma}$')
    #\n{1} pairs of genes. ~{2} retained from WGD.".format(
    #    species_name, num_pairs, num_after_wgd))
    fig.suptitle(species_for_plot_title, style='italic')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return n, bins

def overlay_differences_in_curves(species_for_plot_title, list_of_hist_data, WGD_ks,out_png):

    colors = [config.allo_color,config.color_blind_friendly_color_cycle_analogs['green'],
              config.color_blind_friendly_color_cycle_analogs['brown']]
    labels = ['SpecKS','truth']


    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])
        bar_plot_xs=[b + i*width for b in bins[0:len(bins)-1]] #to match bar-plot axes
        sub_bins=bins[0:len(bins) - 1]
        label=labels[i]
        #if label=='truth':

            #x_value_of_ymax = get_Ks_at_max(WGD_ks, bar_plot_xs, n, sub_bins)
            #label = label +" (peak at Ks=" +str(WGD_ks) + ")"
            #plt.axvline(WGD_ks, color=config.color_blind_friendly_color_cycle_analogs['green'],
            #            linestyle=':')
        plt.plot(bar_plot_xs,n,c=colors[i], alpha=1, label=label)
        #plt.plot(sub_bins, n, c=colors[i], alpha=1, label=label)

    diffs = [(list_of_hist_data[0][0][j]-list_of_hist_data[1][0][j]) for j in range(0,len(list_of_hist_data[1][0]))]
    rmse=math.sqrt(sum([d*d for d in diffs])/len(diffs))
    hist_data0 = list_of_hist_data[1]
    ys0=hist_data0[0]
    xs0=hist_data0[1]

    for i in range(0,len(ys0)):
        if i==0:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15, label="rmse: "+ str(round(rmse,4))))
        else:
            ax.add_patch(Rectangle((xs0[i], ys0[i]), width, diffs[i], color=colors[2],
                                alpha=0.15))

    plt.legend()
    plt.xlabel("Ks")
    ax.set(ylabel="Density")
    plt.title(species_for_plot_title, style='italic')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return n, bins


def get_Ks_at_max(WGD_ks, bar_plot_xs, n, sub_bins):
    true_peak_xs = [];
    true_peak_ys = []
    true_peak_bar_xs = []
    padding = 0.02
    for bin_i in range(0, len(sub_bins)):
        xi = sub_bins[bin_i]
        yi = n[bin_i]
        bi = bar_plot_xs[bin_i]
        if (xi > WGD_ks - padding) and (xi < WGD_ks + padding):
            true_peak_xs.append(xi)
            true_peak_ys.append(yi)
            true_peak_bar_xs.append(bi)
    center_of_mass, x_value_of_ymax = curve_fitting.get_mode_and_cm(
        true_peak_ys, true_peak_xs)
    return x_value_of_ymax


if __name__ == '__main__':
    unittest.main()