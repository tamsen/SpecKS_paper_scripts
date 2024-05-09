import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
import batch_histogrammer,time_series_histogrammer
import config

class MyTestAcrossNe(unittest.TestCase):
    def test_across_ne(self):

        batch_names = ["sim37_N0p1","sim37_N1","sim37_N5","sim36_N10","sim37_N20"]
        plot_title=("The effect of ancestral Ne on Ks histogram shape")

        main_dir="/home/tamsen/Data/Specks_outout_from_mesx"
        output_folder=os.path.join(main_dir)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        print("Reading csv files..")
        SPC_filter=10
        WGD_filter=5
        auto_csv_files_by_batch={}
        allo_csv_files_by_batch={}
        for batch in batch_names:
            input_folder = os.path.join(main_dir, batch)
            csvfiles_by_polyploid_by_species_rep_by_algorithm = batch_histogrammer.get_ks_data_from_folders(input_folder)
            params_by_polyploid = batch_histogrammer.get_truth_from_name_dict(csvfiles_by_polyploid_by_species_rep_by_algorithm)
            alg="ML"
            runs=params_by_polyploid.keys()
            for run in runs:
                if (params_by_polyploid[run].SPC_time_MYA == SPC_filter):
                    if params_by_polyploid[run].WGD_time_MYA == SPC_filter:
                        #its the autopolyploid we want
                        auto_csv_files_by_batch[batch]=[csvfiles_by_polyploid_by_species_rep_by_algorithm[run],
                                                    params_by_polyploid[run]]
                    elif  params_by_polyploid[run].WGD_time_MYA == WGD_filter:
                        allo_csv_files_by_batch[batch]=[csvfiles_by_polyploid_by_species_rep_by_algorithm[run],
                                                    params_by_polyploid[run]]

        print("Making plots..")
        bin_size = 0.001
        max_Ks = 0.6
        max_Y= 150

        #return
        get_Ne_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     auto_csv_files_by_batch,
                                                    allo_csv_files_by_batch,
                                                         'polyploid',
                                                     "rep0", alg,
                                                      max_Y, max_Ks, bin_size)


def get_Ne_series_histograms_for_runs_in_batch(out_folder, sample_name,
                                                 auto_csv_files_by_batch,
                                                 allo_csv_files_by_batch,
                                     spec, replicate, alg,max_Y, max_Ks_for_x_axis,
                                     bin_size):


    batches=list(allo_csv_files_by_batch.keys())
    ne_string = [s.split('_')[-1].replace('N','').replace('p','.') for s in batches]
    num_batches=len(batches)
    fig, ax = plt.subplots(2, 5,figsize=(20,10))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")

    for batch_idx in range(0,num_batches):

            row=0;col=batch_idx
            ax[row,col].set_title("Ne*Gt="+ne_string[batch_idx], fontsize=20)
            [csvs_for_allo_result,params]= allo_csv_files_by_batch[batches[batch_idx]]
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]
            hist_data=make_Ne_series_histogram_subplot(ax[row,col], spec, ks_for_allo_result,
                                                       bin_size, params,
                                                       max_Ks_for_x_axis, max_Y,config.allo_color)
            out_file_name = os.path.join(out_folder, "allo" + batches[batch_idx] + ".hist.csv")
            time_series_histogrammer.save_hist_to_csv(hist_data, out_file_name)

            row=1;col=batch_idx
            [csvs_for_auto_result,params]= auto_csv_files_by_batch[batches[batch_idx]]
            ks_for_auto_result= csvs_for_auto_result[spec][replicate][alg]
            hist_data=make_Ne_series_histogram_subplot(ax[row,4-col], spec, ks_for_auto_result,
                                             bin_size, params,
                                              max_Ks_for_x_axis, max_Y,config.auto_color)

            ax[1, batch_idx].set(xlabel="Ks")

            out_file_name = os.path.join(out_folder, "auto" + batches[batch_idx] + ".hist.csv")
            time_series_histogrammer.save_hist_to_csv(hist_data, out_file_name)



    ax[1,0].set(ylabel=config.histogram_y_label)
    ax[0,0].set(ylabel=config.histogram_y_label)
    plt.tight_layout()
    out_file_name=os.path.join(out_folder, "histogram_across_ne" + "_plot_" + spec +
                               "_" + replicate + "_" + str(max_Ks_for_x_axis) + ".png")
    plt.savefig(out_file_name)
    plt.close()
    return


def make_Ne_series_histogram_subplot(this_ax, spec, Ks_results, bin_size, params,
                                     max_Ks, maxY, plot_color):

    x = Ks_results
    bins = np.arange(0, max_Ks + 0.1, bin_size)
    n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color,
                                    alpha=1)
    hist_result=[n, bins]
    yaxis_limit=maxY
    this_ax.set(ylim=[0, yaxis_limit])
    this_ax.set(xlim=[0, max_Ks])


    return hist_result


if __name__ == '__main__':
    unittest.main()
