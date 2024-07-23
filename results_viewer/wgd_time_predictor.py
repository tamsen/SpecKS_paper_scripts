import math
import os
import random
import statistics
import unittest
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from sklearn.metrics import mean_squared_error
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

import config
import two_d_colors
from results_viewer import curve_fitting
from results_viewer.batch_histogrammer import get_truth_from_name_list


class WGDPredictor(unittest.TestCase):

    def test_wgd2_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1 = "genes_remaining_vs_wgd_time.csv"
        file2="mode_vs_spec_time.csv"
        file3="mode_vs_wgd_time.csv"
        full_path1 = os.path.join(out_folder, file1)
        full_path2 = os.path.join(out_folder, file2)
        #full_path3 = os.path.join(out_folder, file3)

        wgd_xs,genes_remaining,wgd_sims = read_xs_ys_csv(full_path1)
        genes_remaining_dict = dict(zip(wgd_sims,genes_remaining))
        wgd_popt=[ 43.89486241,574.88694572,143.86691067]
        print("wgd_popt:\t" + str(wgd_popt))

        spec_xs,mode,spec_sims = read_xs_ys_csv(full_path2)
        mode_dict = dict(zip(spec_sims, mode))
        spec_popt=[97.05019317,-0.18999441]
        print("spec_popt:\t" + str(spec_popt))

        truth_by_sim_name=get_truth_from_name_list(wgd_sims)
        allo_vs_auto_truth_by_sim={}
        allo_vs_auto_prediction_by_sim={}
        print("starting truth assemby")
        sims_names=truth_by_sim_name.keys()
        predicted_indexes_for_true_autos=[]
        predicted_indexes_for_true_allos=[]
        for sim_name in sims_names:

            polyploid_params=truth_by_sim_name[sim_name]
            how_auto_t=polyploid_params.SPC_time_MYA - polyploid_params.WGD_time_MYA
            truth=[polyploid_params.SPC_time_MYA,polyploid_params.WGD_time_MYA,how_auto_t]
            allo_vs_auto_truth_by_sim[sim_name]=truth
            print("true spec {0}\ttrue wgd {1}".format(polyploid_params.SPC_time_MYA,polyploid_params.WGD_time_MYA))

            genes_remaing= genes_remaining_dict[sim_name]
            wgd_prediction=curve_fitting.logfit(genes_remaing, *wgd_popt)
            mode= mode_dict[sim_name]
            spec_prediction=curve_fitting.linear(mode, *spec_popt)
            how_auto_predicition=spec_prediction - wgd_prediction
            prediction = [spec_prediction, wgd_prediction,how_auto_predicition]
            allo_vs_auto_prediction_by_sim[sim_name]=prediction
            #print("predicted index:\t" + str(how_auto_predicition))
            print("predicted spec {0}\tpredicted wgd {1}".format(spec_prediction, wgd_prediction))
            print("\n")

            if (polyploid_params.WGD_time_MYA==polyploid_params.SPC_time_MYA):
                predicted_indexes_for_true_autos.append(how_auto_predicition)
            else:
                predicted_indexes_for_true_allos.append(how_auto_predicition)

        highest_predicted_index_for_true_autos=max(predicted_indexes_for_true_autos)
        lowest_predicted_index_for_true_allos=min(predicted_indexes_for_true_allos)
        print("highest_predicted_index_for_true_autos:\t" + str(highest_predicted_index_for_true_autos))
        print("lowest_predicted_index_for_true_allos:\t" + str(lowest_predicted_index_for_true_allos))

        discrim_criteria_midpoint= 0.5 *(highest_predicted_index_for_true_autos + lowest_predicted_index_for_true_allos)
        print("midpoint:\n" + str(discrim_criteria_midpoint))

        tests=["SPEC","WGD","degree of allopolyploidy (true ΔT)"]

        # convert the list of key-value pairs to a dictionary
        errors_by_test = {test: [] for test in tests}
        sims_names_list=list(sims_names)
        num_data_points=len(sims_names)
        for test_i in range(0,len(tests)):
            ava_truth=[allo_vs_auto_truth_by_sim[s][test_i] for s in sims_names_list]
            ava_predictions=[allo_vs_auto_prediction_by_sim[s][test_i] for s in sims_names_list]
            errors_for_test = [abs(ava_predictions[j]-ava_truth[j]) for j in range(0,len(ava_truth))]
            errors_by_test[tests[test_i]]=errors_for_test
            plot_data = [sims_names_list,ava_truth,ava_predictions]

            plot_file = plot_data_and_CI(ava_predictions, ava_truth, discrim_criteria_midpoint, num_data_points,
                                              out_folder, sims_names_list, test_i, tests)

            data_file = plot_file.replace("png", "csv")
            save_metrics_to_csv(plot_data,data_file)

        for test_i in range(0,len(tests)):
            ava_truth=[allo_vs_auto_truth_by_sim[s][test_i] for s in sims_names_list]
            ava_predictions=[allo_vs_auto_prediction_by_sim[s][test_i] for s in sims_names_list]
            errors_for_test = errors_by_test[tests[test_i]]
            plot_data_and_CI(ava_predictions, ava_truth, discrim_criteria_midpoint, num_data_points, out_folder,
                             sims_names_list, test_i, tests)

            metric_name = "degree of allopolyploidy (true ΔT)" #allopolyploid_index"
            metric = [allo_vs_auto_truth_by_sim[s][2] for s in sims_names_list]
            plot_error_vs_metric(errors_for_test, metric, metric_name,
                                 sims_names_list, test_i, tests, out_folder)

        mode_predictions = [100 * mode_dict[s] for s in sims_names_list]
        wgd_truths = [allo_vs_auto_truth_by_sim[s][1] for s in sims_names_list]
        mode_predictions_auto_only = [100 * mode_dict[s] for s in sims_names_list if "Auto" in s]
        wgd_truths_auto_only = [allo_vs_auto_truth_by_sim[s][1] for s in sims_names_list if "Auto" in s]
        mode_predictions_allo_only = [100 * mode_dict[s] for s in sims_names_list if "Allo" in s]
        wgd_truths_allo_only = [allo_vs_auto_truth_by_sim[s][1] for s in sims_names_list if "Allo" in s]

        tests=["mode vs WGD time"]
        errors=[mode_predictions[j]-wgd_truths[j] for j in range(0,len(mode_predictions))]

        plot_data_and_CI_for_mode_prediction(mode_predictions_allo_only,wgd_truths_allo_only,
                                             mode_predictions_auto_only, wgd_truths_auto_only,
                                 out_folder, sims_names_list, 0, tests)

        metric = [allo_vs_auto_truth_by_sim[s][2] for s in sims_names_list]
        plot_error_vs_metric(errors, metric, "degree of allopolyploidy (true ΔT)",
                             sims_names_list, 0, tests, out_folder)

    def test_highN_vs_lowN_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="metric3.csv"
        full_path1 = os.path.join(out_folder, file1)

        spc_xs,metric,wgd_sims = read_xs_ys_csv(full_path1)
        low_N_sims= ["Auto","sim37_N0p1", "sim37_N1"]
        med_N_sims=["sim37_N5"]# "sim35_log"]
        high_N_sims=["sim36_N10","sim37_N20"]

        truth_by_sim={}
        metrics_by_category = {"Low":[],"Medium":[],"High":[]}
        spec_by_category = {"Low":[],"Medium":[],"High":[]}
        cutoff=60
        for i in range(0,len(wgd_sims)):
            sim = wgd_sims[i]
            spc_time=spc_xs[i]
            #true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
            #truth_by_sim[sim]=true_category
            if spc_time >= cutoff:
                continue
            if metric[i] <= 0:
                continue
            true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
            #truth_by_sim[sim]=true_category
            #print(metric[i])
            #plus_1=metric[i]+1
            truth_by_sim[sim] = true_category
            metrics_by_category[true_category].append(metric[i])
            spec_by_category[true_category].append(spc_xs[i])


        max_low_N_result=max(metrics_by_category["Low"])
        min_med_N_result=min(metrics_by_category["Medium"])
        max_med_N_result=max(metrics_by_category["Medium"])
        min_high_N_result=min(metrics_by_category["High"])

        print(statistics.stdev(metrics_by_category["Low"]))
        low_n_mean=(statistics.mean(metrics_by_category["Low"]))
        med_n_mean=(statistics.mean(metrics_by_category["Medium"]))
        high_n_mean=(statistics.mean(metrics_by_category["High"]))
        low_vs_medium_discrimination = 0.5 * sum([max_low_N_result, min_med_N_result ]) #0.01
        medium_vs_high_discrimination = 0.5 * sum([max_med_N_result, min_high_N_result]) #0.05

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors_by_category= {"Low":"gray","Medium":"cyan","High":"blue"}

        accuracy_by_sim={}
        for i in range(0,len(wgd_sims)):
            sim = wgd_sims[i]
            true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
            plt.scatter(spc_xs[i],metric[i], alpha=0.25,
                        c=colors_by_category[true_category])

            predicted_category=predict_category(
                low_vs_medium_discrimination, medium_vs_high_discrimination,metric[i])
            accuracy_by_sim[sim]=[true_category,predicted_category]

        ax.axhline(y=low_vs_medium_discrimination, color='cyan', linestyle='--', label="lvm disc. criteria"
                + " ({0})".format(round(low_vs_medium_discrimination, 2)))

        ax.axhline(y=medium_vs_high_discrimination, color='purple', linestyle='--', label="mvh disc. criteria"
                + " ({0})".format(round(medium_vs_high_discrimination, 2)))

        ax.set(yscale='log')
        plt.legend()
        plot_file = os.path.join(out_folder, "highN_vs_lowN_predictor.png")
        plt.savefig(plot_file)
        plt.close()

        plot_confusion(accuracy_by_sim, colors_by_category, high_n_mean, low_n_mean, low_vs_medium_discrimination,
                            med_n_mean, medium_vs_high_discrimination, out_folder)

        make_box_plots(colors_by_category, metrics_by_category, out_folder)

        #data_file = os.path.join(out_folder,"highN_vs_lowN_truth_and_predictions.csv")
        #save_metrics_to_csv(plot_data, data_file)


def plot_error_vs_metric(error, metric, metric_name,
                     sims_names_list, test_i, tests, out_folder):
    print(sims_names_list)
    colors = [two_d_colors.get_color_from_dT(s) for s in sims_names_list]
    ci_percent=99.99
    alphas = [0.25 if "Auto" in s else 0.25 for s in sims_names_list]
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    #TODO - add legend..?
    plt.scatter( metric,error, alpha=alphas,
                c=colors)

    foo = pd.DataFrame({'error': error,
                        'metric': metric})

    sns.regplot(data=foo, x='metric', y='error', color='k',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(all data)".format(ci_percent),
                line_kws={"linewidth":1})

    #else:
    ax.set(xlabel=metric_name)
    ax.set(ylabel="prediction error")
    ax.set(ylim=[-10,65])
    plt.legend()
    plot_file = os.path.join(out_folder, tests[test_i] + "_error_vs_" + metric_name +".png")
    plt.savefig(plot_file)
    plt.close()
    return plot_file

def plot_data_and_CI_for_mode_prediction(ava_predictions_allo, ava_truth_allo,
                                         ava_predictions_auto, ava_truth_auto,out_folder,
                                        sims_names_list, test_i, tests):

    ci_percent=99.99#60#99.99
    ci_shading_auto = ["CI at {0}% (auto only)".format(ci_percent) for s in sims_names_list if "Auto" in s]
    ci_shading_allo = ["CI at {0}% (allo only)".format(ci_percent) for s in sims_names_list if "Allo" in s]

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    #TODO: get color based on delta T
    dt_colors_auto = [two_d_colors.get_color_from_dT(s) for s in sims_names_list if "Auto" in s]
    dt_colors_allo = [two_d_colors.get_color_from_dT(s) for s in sims_names_list if "Allo" in s]
    have_tan_label=False
    have_orange_label=False
    plt.scatter(ava_truth_auto, ava_predictions_auto, alpha=1,
                c=dt_colors_auto, label="low ΔT")

    for i in range(0,len(ava_truth_allo)):

        color=dt_colors_allo[i]

        if color=='tan':
            if have_tan_label:
                plt.scatter(ava_truth_allo[i], ava_predictions_allo[i], alpha=0.5,
                            c=dt_colors_allo[i])
            else:
                plt.scatter(ava_truth_allo[i], ava_predictions_allo[i], alpha=0.5,
                    c=dt_colors_allo[i], label="high ΔT")
                have_tan_label=True

        if color=='orange':
            if have_orange_label:
                plt.scatter(ava_truth_allo[i], ava_predictions_allo[i], alpha=0.5,
                    c=dt_colors_allo[i])
            else:
                plt.scatter(ava_truth_allo[i], ava_predictions_allo[i], alpha=0.5,
                    c=dt_colors_allo[i], label="medium ΔT")
                have_orange_label=True
    foo_auto = pd.DataFrame({'truth': ava_truth_auto,
                        'prediction': ava_predictions_auto,
                        'CI': ci_shading_auto})

    foo_allo = pd.DataFrame({'truth': ava_truth_allo,
                        'prediction': ava_predictions_allo,
                        'CI': ci_shading_allo})

    #sns.lineplot(data=foo, x='truth', y='prediction', hue='CI', palette=['k'],
    #             errorbar=('ci', ci_percent))
    sns.regplot(data=foo_auto, x='truth', y='prediction', color='darkred',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(low ΔT only)".format(ci_percent),
                line_kws={"linewidth":1})
    '''
    sns.regplot(data=foo_allo, x='truth', y='prediction', color='darkblue',ci=ci_percent,
                marker=None,scatter=False,
                label="CI at {0}%\n(allos only)".format(ci_percent),
                line_kws={"linewidth":1})

    sns.regplot(data=foo_allo, x='truth', y='prediction', color='darkblue',ci=60,
                marker=None,scatter=False,
                label="CI at {0}%\n(allos only)".format(60),
                line_kws={"linewidth":1})
    '''
    ax.set(xlabel="truth (MY)")
    ax.set(ylabel="prediction (MY)")
    plt.legend()
    plot_file = os.path.join(out_folder, tests[test_i] + "_truth_vs_predictions.png")
    plt.savefig(plot_file)
    plt.close()
    return plot_file

def plot_data_and_CI(ava_predictions, ava_truth, discrim_criteria_midpoint, num_data_points, out_folder,
                     sims_names_list, test_i, tests):
    #colors = ["red" if "Auto" in s else "blue" for s in sims_names_list]
    ci_percent=99.99



    #colors = [config.low_Ne_low_dT_color if "Auto" in s else config.allo_color for s in sims_names_list]
    dt_colors = [two_d_colors.get_color_from_dT(s) for s in sims_names_list]

    dTs=[two_d_colors.get_dT_from_name(s) for s in sims_names_list]
    print("dts:" + str(dTs))
    dNs=[two_d_colors.get_Ne_from_name(s) for s in sims_names_list]
    print("dNs:" + str(dNs))
    ci_shading = ["CI at {0}% (auto and allo)".format(ci_percent) for s in sims_names_list]
    #alphas = [1 if "Auto" in s else 0.5 for s in sims_names_list]

    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    have_low_dT_label=False; have_tan_label=False; have_orange_label=False

    #TODO, fix labels

    for i in range(0,len(ava_truth)):

        color=dt_colors[i]

        #if not have_low_dT_label:
        #    plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
        #                    c=dt_colors[i], s=100, label="ΔT")
        #    dT_labeled = True
        #else:
        #    plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
        #                    c=dt_colors[i], s=50)

        if color == two_d_colors.two_d_colors.low_dT:
            if have_low_dT_label:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i])
            else:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i], label="low ΔT")
                have_low_dT_label = True

        elif color == 'tan':
            if have_tan_label:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i])
            else:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i], label="high ΔT")
                have_tan_label = True
        elif color == 'orange':
            if have_orange_label:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i])
            else:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
                        c=dt_colors[i], label="medium ΔT")
                have_orange_label = True

    #for i in range(0,len(ava_truth)):

    #    if not Ne_labeled:
    #        plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
    #                c=ne_colors[i], s=10, marker="X", label=''r'$\pi$')
    #        Ne_labeled = True
    #    else:
    #        plt.scatter(ava_truth[i], ava_predictions[i], alpha=0.5,
    #            c=ne_colors[i], s=10, marker="X")

    foo = pd.DataFrame({'truth': ava_truth, 'prediction': ava_predictions, 'CI': ci_shading})
    #sns.lineplot(data=foo, x='truth', y='prediction', hue='CI', palette=['k'],errorbar=('ci', ci_percent) )
    sns.regplot(data=foo, x='truth', y='prediction', color='k',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(all data)".format(ci_percent),
                line_kws={"linewidth":1})

    if test_i == 2:
        ax.set(xlabel="<-- true ΔT -->")
        ax.set(ylabel="<-- predicted ΔT (MY) -->")
        ax.set(title="Predictions of polyploid index (ΔT, in MY)")

        ax.axhline(y=discrim_criteria_midpoint, color='gray', linestyle='--', label="disc. criteria"
                                                                                 + " (y={0})".format(
            round(discrim_criteria_midpoint, 2)))

        ax.axhline(y=30, color='gray', linestyle='--', label="disc. criteria"
                                                                                 + " (y=TODO)")
    else:
        ax.set(xlabel="truth (MY)")
        ax.set(ylabel="prediction (MY)")
    plt.legend()
    plot_file = os.path.join(out_folder, tests[test_i] + "_truth_vs_predictions.png")
    plt.savefig(plot_file)
    plt.close()
    return plot_file

    #https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
def plot_confusion(accuracy_by_sim, colors_by_category, high_n_mean, low_n_mean, low_vs_medium_discrimination,
                   med_n_mean, medium_vs_high_discrimination, out_folder):
    means_by_category = {"Low": low_n_mean, "Medium": med_n_mean, "High": high_n_mean}
    stutter = 0.5
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.add_patch(Rectangle((0, 0), low_vs_medium_discrimination, 1, color='gray', alpha=0.15))
    ax.add_patch(Rectangle((0, 0), 1, low_vs_medium_discrimination, color='gray', alpha=0.15))
    ax.add_patch(Rectangle((low_vs_medium_discrimination, 0),
                           medium_vs_high_discrimination - low_vs_medium_discrimination, 1,
                           color='cyan', alpha=0.15))
    ax.add_patch(Rectangle((0, low_vs_medium_discrimination),
                           1, medium_vs_high_discrimination - low_vs_medium_discrimination,
                           color='cyan', alpha=0.15))
    ax.add_patch(Rectangle((medium_vs_high_discrimination, 0),
                           1 - medium_vs_high_discrimination, 1,
                           color='blue', alpha=0.15))
    ax.add_patch(Rectangle((0, medium_vs_high_discrimination),
                           1, 1 - medium_vs_high_discrimination,
                           color='blue', alpha=0.15))
    for sim, results in accuracy_by_sim.items():
        [true_category, predicted_category] = accuracy_by_sim[sim]
        x_value = means_by_category[true_category] * (1 + stutter * random.random())
        y_value = means_by_category[predicted_category] * (1 + stutter * random.random())

        plt.scatter(x_value, y_value, alpha=0.25,
                    c=colors_by_category[true_category])
    ax.axvline(x=low_vs_medium_discrimination, color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                                   + " ({0})".format(
        round(low_vs_medium_discrimination, 2)))
    # ax.add_patch(Rectangle((0, 0), low_vs_medium_discrimination, 1), color='cyan')
    ax.axvline(x=medium_vs_high_discrimination, color='purple', linestyle='--', label="mvh disc. criteria"
                                                                                      + " ({0})".format(
        round(medium_vs_high_discrimination, 2)))
    ax.set(xlabel="<-- truth -->")
    ax.set(ylabel="<-- prediction -->")
    ax.set(xscale='log')
    ax.set(yscale='log')
    plt.legend()
    plot_file = os.path.join(out_folder, "highN_vs_lowN_accuracy.png")
    plt.savefig(plot_file)
    plt.close()


def predict_category(low_vs_medium_discrimination, medium_vs_high_discrimination, metric):
    if metric < low_vs_medium_discrimination:
        return "Low"
    elif metric < medium_vs_high_discrimination:
        return "Medium"
    return "High"

def make_box_plots(colors_by_category, metrics_by_category, out_folder):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    # https://matplotlib.org/stable/gallery/statistics/boxplot_color.html#sphx-glr-gallery-statistics-boxplot-color-py
    boxplot_data = [metrics_by_category["Low"], metrics_by_category["Medium"], metrics_by_category["High"]]
    labels = ["Low", "Medium", "High"]
    bplot = ax.boxplot(boxplot_data, labels=labels, patch_artist=True)
    for patch, label in zip(bplot['boxes'], labels):
        patch.set_facecolor(colors_by_category[label])
    plot_file = os.path.join(out_folder, "highN_vs_lowN_boxplots.png")
    plt.savefig(plot_file)
    plt.close()


def categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim):
    for n in low_N_sims:
        if n in sim:
            return "Low"
    for n in med_N_sims:
        if n in sim:
            return "Medium"
    for n in high_N_sims:
        if n in sim:
            return "High"
    else:
        return False


if __name__ == '__main__':
    unittest.main()

def save_metrics_to_csv(plot_data_for_sims, out_file_name):

    with open(out_file_name, 'w') as f:
        data_headers= ['sims','truth','prediction']
        f.writelines(",".join(data_headers) +"\n")
        [sims_names, ava_truth, prediction] = plot_data_for_sims
        for i in range(0,len(sims_names)):
            data_list=[sims_names[i], ava_truth[i], prediction[i]]
            data_list_string=[str(d) for d in data_list]
            f.writelines(",".join(data_list_string) +"\n")


def read_xs_ys_csv(csv_file):

    batches=[]
    sim_names=[]
    xs=[]
    ys=[]
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "batch" in line:
                continue
            if len(line)==0:
                break
            if "NA" in line:
                continue
            data = line.strip().split(",")
            #batches.append(data[0])
            sim_names.append(data[1]+ "_" + data[0])
            xs.append(float(data[2]))
            ys.append(float(data[3]))

    return xs,ys,sim_names