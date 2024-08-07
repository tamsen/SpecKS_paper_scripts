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
from results_viewer import curve_fitting
from results_viewer.batch_histogrammer import get_truth_from_name_list


class AlloAutoPredictor(unittest.TestCase):
    def test_spec_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="mode_vs_spec_time.csv"
        full_path=os.path.join(out_folder,file1)
        spec_xs,mode_ys, sims = read_xs_ys_csv(full_path)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors=['b' if "Allo" in s else 'r' for s in sims]
        colors=['k' for s in sims]
        plt.scatter(mode_ys, spec_xs, c=colors,label="empirical data",alpha=0.25)

        fit_spec_times, xs, goodness_of_fit=curve_fitting.fit_curve_to_xs_and_ys(
            mode_ys, spec_xs, curve_fitting.linear)

        linear_fit_popt=goodness_of_fit.popt
        ax.set(xlabel="mode (Ks space)")
        ax.set(ylabel="input SPC (MYA)")
        #rms2 = mean_squared_error(xs,fit_spec_times, squared=False)
        parameter_string=[str(round(p,2)) for p in linear_fit_popt]
        r = scipy.stats.pearsonr(xs,fit_spec_times)
        print(r)
        plt.plot(xs,fit_spec_times, c='k',label="line fit (model)"
                                                              +"\nparameters: "+str(parameter_string)
                                                              + "\nresidual: " + str(round(r[0],2)))


        #example_modes=np.arange(0.2, 1.0, 0.2)
        #predicted_spec_times= [curve_fitting.linear(x, *linear_fit_popt) for x in example_modes]
        #plt.bar(example_modes,predicted_spec_times, width=0.002, color='pink',label="model input (mode)")
        #plt.scatter(example_modes,predicted_spec_times, c='r',label="model output (est. spec time)")
        #print("linear_fit_popt:\t" + str(linear_fit_popt))

        plot_file=full_path.replace(".csv",".new.png")
        plt.legend()
        plt.savefig(plot_file)
        plt.close()

    def test_wgd_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="genes_remaining_vs_wgd_time.csv"
        full_path=os.path.join(out_folder,file1)

        wgd_xs,genes_remaining,sims = read_xs_ys_csv(full_path)
        popt, pcov = curve_fit(curve_fitting.logfit, genes_remaining, wgd_xs)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(genes_remaining, wgd_xs, c='k',label="empirical data",alpha=0.25)
        genes_remaining.sort()
        parameter_string=[str(round(p,2)) for p in popt]
        fit_wgd_times = [curve_fitting.logfit(x, *popt) for x in genes_remaining]
        plt.plot(genes_remaining, fit_wgd_times, c='k', label="log fit (model)\n"
                                                              +"parameters: " +str(parameter_string))
        print("logfit_fit_popt:\t" + str(popt))

        #example_genes_remaining=np.arange(500, 3000, 500)
        #predicted_wgd_times= [curve_fitting.logfit(x, *popt) for x in example_genes_remaining]
        #plt.bar(example_genes_remaining,predicted_wgd_times, width=10, color='pink',label="model input (# genes)")
        #plt.scatter(example_genes_remaining,predicted_wgd_times, c='r',label="model output (est. wgd time)")

        ax.set(xlabel="# genes remaining")
        ax.set(ylabel="input WGD (MYA)")
        plt.legend()
        plot_file=full_path.replace(".csv",".new.png").replace("shed","remaining")

        plt.savefig(plot_file)
        plt.close()

    def test_log_wgd_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1 = "genes_remaining_vs_wgd_time.csv"
        full_path = os.path.join(out_folder, file1)

        wgd_xs, genes_remaining, sims = read_xs_ys_csv(full_path)
        log_genes_remaining=[math.log(g) for g in genes_remaining]
        popt, pcov = curve_fit(curve_fitting.linear,log_genes_remaining, wgd_xs)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(log_genes_remaining, wgd_xs, c='k', label="empirical data",alpha=0.25)
        genes_remaining.sort()
        parameter_string = [str(round(p, 2)) for p in popt]
        fit_wgd_times = [curve_fitting.linear(x, *popt) for x in log_genes_remaining]
        r = scipy.stats.pearsonr(wgd_xs,fit_wgd_times)
        plt.plot(log_genes_remaining, fit_wgd_times, c='k', label="log fit (model)\n"
                                                              + "parameters: " + str(parameter_string)
                                                                  + "\nresidual: " + str(round(r[0], 2)))
        print("logfit_fit_popt:\t" + str(popt))

        # example_genes_remaining=np.arange(500, 3000, 500)
        # predicted_wgd_times= [curve_fitting.logfit(x, *popt) for x in example_genes_remaining]
        # plt.bar(example_genes_remaining,predicted_wgd_times, width=10, color='pink',label="model input (# genes)")
        # plt.scatter(example_genes_remaining,predicted_wgd_times, c='r',label="model output (est. wgd time)")

        ax.set(xlabel="ln( # genes remaining)")
        ax.set(ylabel="input WGD (MYA)")
        plt.legend()
        plot_file = full_path.replace(".csv", ".new.log.png").replace("shed", "remaining")

        plt.savefig(plot_file)
        plt.close()

    def test_auto_vs_allo_predictor(self):

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
            how_auto_t= polyploid_params.DIV_time_MYA - polyploid_params.WGD_time_MYA
            truth=[polyploid_params.DIV_time_MYA, polyploid_params.WGD_time_MYA, how_auto_t]
            allo_vs_auto_truth_by_sim[sim_name]=truth
            print("true spec {0}\ttrue wgd {1}".format(polyploid_params.DIV_time_MYA, polyploid_params.WGD_time_MYA))

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

            if (polyploid_params.WGD_time_MYA==polyploid_params.DIV_time_MYA):
                predicted_indexes_for_true_autos.append(how_auto_predicition)
            else:
                predicted_indexes_for_true_allos.append(how_auto_predicition)

        highest_predicted_index_for_true_autos=max(predicted_indexes_for_true_autos)
        lowest_predicted_index_for_true_allos=min(predicted_indexes_for_true_allos)
        print("highest_predicted_index_for_true_autos:\t" + str(highest_predicted_index_for_true_autos))
        print("lowest_predicted_index_for_true_allos:\t" + str(lowest_predicted_index_for_true_allos))

        discrim_criteria_midpoint= 0.5 *(highest_predicted_index_for_true_autos + lowest_predicted_index_for_true_allos)
        print("midpoint:\n" + str(discrim_criteria_midpoint))

        tests=["SPEC","WGD","allopolyploid_index"]

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

            metric_name = "allopolyploid_index"
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
        plot_error_vs_metric(errors, metric, "allopolyploid_index",
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

        make_box_plots(colors_by_category, metrics_by_category, out_folder)

        #data_file = os.path.join(out_folder,"highN_vs_lowN_truth_and_predictions.csv")
        #save_metrics_to_csv(plot_data, data_file)


def plot_error_vs_metric(error, metric, metric_name,
                     sims_names_list, test_i, tests, out_folder):
    colors = [config.auto_color if "Auto" in s else config.allo_color for s in sims_names_list]
    ci_percent=99.99
    alphas = [0.25 if "Auto" in s else 0.25 for s in sims_names_list]
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.scatter( metric,error, alpha=alphas,
                c=colors)

    foo = pd.DataFrame({'error': error,
                        'metric': metric})

    sns.regplot(data=foo, x='metric', y='error', color='k',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(all data)".format(ci_percent),
                line_kws={"linewidth":1})

    #if test_i == 2:
    #    ax.set(xlabel="<--auto    truth (MY)   allo-->")
    #    ax.set(ylabel="<--auto  prediction (MY) allo-->")
    #    ax.set(title="Predictions of polyploid index (SPEC - WGD time, in MY)")

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
    plt.scatter(ava_truth_auto, ava_predictions_auto, alpha=1,
                c=config.low_Ne_low_dT_color, label="auto")
    plt.scatter(ava_truth_allo, ava_predictions_allo, alpha=0.5,
                c=config.allo_color, label="allo")

    foo_auto = pd.DataFrame({'truth': ava_truth_auto,
                        'prediction': ava_predictions_auto,
                        'CI': ci_shading_auto})

    foo_allo = pd.DataFrame({'truth': ava_truth_allo,
                        'prediction': ava_predictions_allo,
                        'CI': ci_shading_allo})

    #sns.lineplot(data=foo, x='truth', y='prediction', hue='CI', palette=['k'],
    #             errorbar=('ci', ci_percent))
    sns.regplot(data=foo_auto, x='truth', y='prediction', color='darkred',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(autos only)".format(ci_percent),
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
    colors = [config.auto_color if "Auto" in s else config.allo_color for s in sims_names_list]
    ci_shading = ["CI at {0}% (auto and allo)".format(ci_percent) for s in sims_names_list]
    alphas = [1 if "Auto" in s else 0.5 for s in sims_names_list]
    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    auto_labeled=False; allo_labeled=False

    for i in range(0,len(ava_truth)):
        color=colors[i]
        if color==config.auto_color and not auto_labeled:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=alphas[i],
                            c=color, label="auto")
                auto_labeled=True
        elif color==config.allo_color and not allo_labeled:
                plt.scatter(ava_truth[i], ava_predictions[i], alpha=alphas[i],
                            c=color, label="allo")
                allo_labeled=True
        else:
            plt.scatter(ava_truth[i], ava_predictions[i], alpha=alphas[i],
                c=color)

    foo = pd.DataFrame({'truth': ava_truth, 'prediction': ava_predictions, 'CI': ci_shading})
    #sns.lineplot(data=foo, x='truth', y='prediction', hue='CI', palette=['k'],errorbar=('ci', ci_percent) )
    sns.regplot(data=foo, x='truth', y='prediction', color='k',ci=ci_percent,marker=None,scatter=False,
                label="CI at {0}%\n(all data)".format(ci_percent),
                line_kws={"linewidth":1})

    if test_i == 2:
        ax.set(xlabel="<--auto    truth (MY)   allo-->")
        ax.set(ylabel="<--auto  prediction (MY) allo-->")
        ax.set(title="Predictions of polyploid index (ΔT), in MY)")

        ax.axhline(y=discrim_criteria_midpoint, color='r', linestyle='--', label="discimination criteria"
                                                                                 + " (y={0})".format(
            round(discrim_criteria_midpoint, 2)))
    else:
        ax.set(xlabel="truth (MY)")
        ax.set(ylabel="prediction (MY)")
    plt.legend()
    plot_file = os.path.join(out_folder, tests[test_i] + "_truth_vs_predictions.png")
    plt.savefig(plot_file)
    plt.close()
    return plot_file

    #https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python

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