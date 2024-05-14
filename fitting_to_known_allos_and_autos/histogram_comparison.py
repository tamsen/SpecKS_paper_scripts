import math
import os
import unittest
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import config
from results_viewer import batch_histogrammer

class MyTestCase(unittest.TestCase):

    def test_coffee_histogram(self):

        coffee_num=13
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_coffee/Allo_Coffea{0}".format(coffee_num)
        specks_csv_file = "Allo_Coffea{0}_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv".format(coffee_num)
        ksrates_csv_file="coffea.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        species_for_plot_title= 'Coffea arabica'
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.002
        density = 1

        make_both_histograms(bin_size, color, specks_out_folder,
                             wgd_ks, max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path)

    def test_poplar_histogram(self):

        pop_num=14
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim41_Poplar/Allo_Poplar{0}".format(pop_num)
        specks_csv_file = "Allo_Poplar{0}_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv".format(pop_num)
        ksrates_csv_file = "poplar.ks.tsv"

        splat = specks_csv_file.split("_")
        species_run_name = splat[0] + splat[1]
        species_for_plot_title= 'Populus trichocarpa'
        specks_full_path = os.path.join(specks_out_folder, specks_csv_file)
        real_full_path = os.path.join(ksrates_out_folder, ksrates_csv_file)

        bin_size = 0.002
        max_Ks = 0.5
        color = 'blue'
        wgd_ks=0.18
        density = 1

        make_both_histograms(bin_size, color, specks_out_folder, wgd_ks,
                             max_Ks, density, real_full_path,
                             species_run_name, species_for_plot_title, specks_full_path)
    def test_maize_histogram(self):

        maize_num="12"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_Maize/Allo{0}_Maize".format(maize_num)
        specks_csv_file = "Allo{0}_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv".format(maize_num)
        ksrates_csv_file="mays.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        species_for_plot_title= 'Zea mays'
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.3
        color='blue'
        wgd_ks=0.125
        density = 1

        make_both_histograms(bin_size, color,  specks_out_folder, wgd_ks,
                         max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path)
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
        real_ks_results = parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.3
        WGD_ks=0.001
        color=config.allo_color
        density = None

        make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks,max_Ks, density,real_ks_results,
                                  species_run_name, specks_ks_results)

def make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks, max_Ks, density, real_full_path,
                         species_run_name, species_for_plot_title, specks_full_path):

    specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)
    real_ks_results = parse_external_ksfile(real_full_path)

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

    fig = plt.figure(figsize=(10, 10), dpi=100)
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

    fig = plt.figure()
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
    rmse=math.sinh(sum([d*d for d in diffs])/len(diffs))
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
        xs=[b + i*width for b in bins[0:len(bins)-1]]
        plt.plot(xs,n,c=colors[i], alpha=1, label=labels[i])

    diffs = [(list_of_hist_data[0][0][j]-list_of_hist_data[1][0][j]) for j in range(0,len(list_of_hist_data[1][0]))]
    rmse=math.sinh(sum([d*d for d in diffs])/len(diffs))
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

def parse_external_ksfile(ks_file):

    if ".fa" in ks_file: #1KP file
        olea_ks_df = pd.read_csv(ks_file, sep='\t', header=0)
        olea_ks_array = olea_ks_df.loc[:, "Node Ks"]
        ks_data = [k for k in olea_ks_array.tolist() if k <= 2]
    else:    #KS rates input file
        ks_data = parse_ks_rates(ks_file)
    return ks_data

def parse_ks_rates(input_file):
    ks_df=pd.read_csv(input_file, sep='\t', header=0)
    ks_array = ks_df.loc[:,"Ks"]
    ks_list = [k for k in ks_array.tolist() if k <= 2]
    #print(ks_list[1:10])
    return ks_list
if __name__ == '__main__':
    unittest.main()