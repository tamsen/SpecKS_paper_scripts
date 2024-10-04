import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
import batch_histogrammer
import two_d_colors, config
from results_viewer.batch_analyzer import smooth_data


class TimeSeriesHistogrammer(unittest.TestCase):

    def test_make_time_series_histograms_for_off_diagonals(self):

        batch_name="sim53_offdiags" #
        input_folder=os.path.join( "/home/tamsen/Data/Specks_outout_from_mesx",batch_name)
        output_folder=input_folder#os.path.join(input_folder,"analysis")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = batch_histogrammer.get_ks_data_from_folders(input_folder)
        params_by_polyploid = batch_histogrammer.get_truth_from_name_dict(csvfiles_by_polyploid_by_species_rep_by_algorithm)
        alg="ML"

        print("Making plots..")
        bin_size = 0.001
        max_Ks = 0.85
        max_Ys= [100]

        #off-diagonals

        filter = "Allo"
        colors=two_d_colors.two_d_colors.high_Ne_low_dT
        plot_title=("High Ne, low ΔT")
        for spec in ['polyploid']:#species:
            get_time_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     "rep0", alg, params_by_polyploid,
                                                      max_Ys, max_Ks, bin_size,filter,colors)

        filter = "Auto"
        colors=two_d_colors.two_d_colors.low_Ne_high_dT
        plot_title=("Low Ne, high ΔT")
        for spec in ['polyploid']:#species:
            get_time_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     "rep0", alg, params_by_polyploid,
                                                      max_Ys, max_Ks, bin_size,filter,colors)


    def test_make_time_series_histograms_for_diagonals(self):

        batch_name="sim52_diagonals" #
        input_folder=os.path.join( "/home/tamsen/Data/Specks_outout_from_mesx",batch_name)
        output_folder=input_folder#os.path.join(input_folder,"analysis")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = batch_histogrammer.get_ks_data_from_folders(input_folder)
        params_by_polyploid = batch_histogrammer.get_truth_from_name_dict(csvfiles_by_polyploid_by_species_rep_by_algorithm)
        alg="ML"

        print("Making plots..")
        bin_size = 0.001
        max_Ks = 0.85
        max_Ys= [100]

        #off-diagonals

        filter = "Allo"
        colors=two_d_colors.two_d_colors.high_Ne_high_dT
        plot_title=("High Ne, high ΔT")
        for spec in ['polyploid']:#species:
            get_time_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     "rep0", alg, params_by_polyploid,
                                                      max_Ys, max_Ks, bin_size,filter,colors)

        filter = "Auto"
        colors=two_d_colors.two_d_colors.low_Ne_low_dT
        plot_title=("Low Ne, low ΔT")
        for spec in ['polyploid']:#species:
            get_time_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     "rep0", alg, params_by_polyploid,
                                                      max_Ys, max_Ks, bin_size,filter,colors)


def get_time_series_histograms_for_runs_in_batch(out_folder, sample_name, csvfiles_by_polyploid_by_rep_by_algorthim,
                                     spec, replicate, alg, params_by_polyploid,max_Ys, max_Ks_for_x_axis,
                                     bin_size,filter, base_color):

    result_names=(list(csvfiles_by_polyploid_by_rep_by_algorthim.keys()))
    result_names.sort()
    ordered_results=[n for n in result_names if filter in n]
    metrics_by_result_names= {}

    # making subplots
    spec_times= [params_by_polyploid[name].DIV_time_MYA for name in ordered_results]
    num_spec_times=len(spec_times)
    num_modes_of_speciation = 1#auto and allow

    fig, ax = plt.subplots(1, num_modes_of_speciation,figsize=(5,5))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")
    #ax.set_title(filter+ "\n",fontsize=15)
    red_color_interval = 1.4 / num_spec_times
    blue_color_interval = 2 / num_spec_times
    for sim_idx in range(0, num_spec_times):

        for spec_mode_idx in range(0,num_modes_of_speciation):

            print(spec_mode_idx)
            allo_result_name=ordered_results[sim_idx]

            print(allo_result_name)

            if ((base_color ==two_d_colors.two_d_colors.low_Ne_high_dT)
                    or (base_color == two_d_colors.two_d_colors.low_Ne_low_dT)) :
                amount = 0.4 + sim_idx * red_color_interval
            else:
                amount = 0.2 + sim_idx *blue_color_interval

            new_color = lighten_color(base_color, amount=amount)

            alpha = 1
            csvs_for_allo_result= csvfiles_by_polyploid_by_rep_by_algorthim[allo_result_name]
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]
            params=params_by_polyploid[allo_result_name]
            hist_data=make_time_series_histogram_subplot(ax, alpha, ks_for_allo_result,
                                                         bin_size, params,
                                                         max_Ks_for_x_axis,
                                                         max_Ys[spec_mode_idx],new_color )

            out_file_name = os.path.join(out_folder, allo_result_name + ".hist.csv")
            save_hist_to_csv(hist_data, out_file_name)

            #TODO
            ax.set(xlabel="~MYA ({0}*Ks)".format(1.0/config.SpecKS_config.Ks_per_Myr))

    ax.set(ylabel=config.histogram_y_label)
    plt.tight_layout()
    out_file_name=os.path.join(out_folder, "histogram" + "_plot_" + spec +
                               "_" + filter + "_" + str(max_Ks_for_x_axis) + ".png")
    plt.savefig(out_file_name)
    plt.close()
    return metrics_by_result_names


def make_time_series_histogram_subplot(this_ax, alpha, Ks_results, bin_size, params,
                                       max_Ks, maxY, plot_color):

    x = Ks_results
    kernel_size=4
    bins = np.arange(0, max_Ks + 0.1, bin_size)
    n, bins, patches = this_ax.hist(x, bins=bins, facecolor='none')#, alpha=alpha,
    #                                label='SPC/WGD time: {0}/{1} MYA'.format(
    #                                    params.SPC_time_MYA, params.WGD_time_MYA) )
    trimmed_bins=bins[0:len(n)]
    smoothed_ys = smooth_data(kernel_size, n)
    this_ax.plot(trimmed_bins, smoothed_ys, color=plot_color,
                 label='DIV/WGD time: {0}/{1} MYA'.format(params.DIV_time_MYA, params.WGD_time_MYA))
    hist_result=[n, bins]
    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6

    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    xticks=this_ax.get_xticks()
    as_my=[x/config.SpecKS_config.Ks_per_Myr for x in xticks]
    tick_labels=[str(round(my,3)) for my in as_my]
    #this_ax.set_xticks(ticks)
    this_ax.set_xticklabels(tick_labels)

    this_ax.legend()
    this_ax.set(ylim=[0, yaxis_limit])
    #this_ax.set(ylim=[0, 150])

    this_ax.set(xlim=[0, max_Ks])


    return hist_result


def save_hist_to_csv(hist_result, out_file_name):

    [n, bins]=hist_result
    with open(out_file_name, 'w') as f:
        data_headers= ['bin','# in bin']
        f.writelines(",".join(data_headers) +"\n")

        for i in range(0,len(n)):
            data_list=[str(bins[i]),str(n[i])]
            f.writelines(",".join(data_list) +"\n")

#from  https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
#from, https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    result= colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    positive_r=[max([0,r])for r in result]
    print(result)
    return positive_r