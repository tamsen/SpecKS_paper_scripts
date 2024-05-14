import os
import unittest

import config
import results_viewer.batch_analyzer
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from results_viewer import curve_fitting, batch_analyzer

from results_viewer.run_metrics import run_metrics, plot_and_save_metrics
def analyze_histogram(bins, n, WGD_time_MYA, SPC_time_MYA,
                      kernel_size, maxY, right_most_ssd_peak, out_file_name):
    linewidth = 4
    polyploid_name=os.path.basename(out_file_name).replace(".png","")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.suptitle('Histogram analysis for ' + polyploid_name)
    WGD_as_Ks = WGD_time_MYA * config.SpecKS_config.Ks_per_Myr
    SPEC_as_Ks =SPC_time_MYA * config.SpecKS_config.Ks_per_Myr
    #plt.bar(bins, n, width=0.001, color=hist_plot_color, label="hist",alpha=0.5)

    hist_maximum=max(n)

    #kernel_size=50
    smoothed_color=config.color_blind_friendly_color_cycle_analogs['purple']
    smoothed_ys = batch_analyzer.smooth_data(kernel_size, n)
    plt.plot(bins, smoothed_ys,color=smoothed_color,
             linewidth=linewidth/2.0,
             label="smoothed data")


    smoothed_minima = batch_analyzer.find_local_minima(bins, smoothed_ys)
    smoothed_maxima = batch_analyzer.find_global_maxima(bins, smoothed_ys, right_most_ssd_peak)
    #plt.scatter([m[0] for m in smoothed_minima], [m[1] for m in smoothed_minima], color="red", label="minima", marker='*')

    ssd_end, next_min =batch_analyzer.smallest_min(smoothed_minima,smoothed_maxima)
    ssds_xs_to_subtract, ssds_ys_to_subtract = batch_analyzer.estimate_overlap(bins, ssd_end, next_min)

    plt.scatter(ssd_end[0],3*ssd_end[1],color="black", label="wgd_start",marker='*',s=50)
    ssd_xs, ssd_ys, wgd_xs, wgd_ys = batch_analyzer.sort_ssds_and_wgds(bins, n, ssd_end, ssds_ys_to_subtract)

    wgd_maxima = batch_analyzer.find_global_maxima(wgd_xs, wgd_ys, 0.075)
    ssd_xs, ssd_ys, wgd_xs, wgd_ys = batch_analyzer.sort_ssds_and_wgds(bins, n, ssd_end, ssds_ys_to_subtract)

    wgd_maxima = batch_analyzer.find_global_maxima(wgd_xs, wgd_ys, 0.075)
    #plt.scatter(wgd_maxima[0], -0.05 * maxY,
    #            color='blue', marker='^', s=80, label="wgd max")

    wgd_max_d=batch_analyzer.get_loc_of_max_derivative(kernel_size, wgd_ys, wgd_xs)
    raw_cm, raw_x_value_of_ymax = curve_fitting.get_mode_and_cm(wgd_ys, wgd_xs)
    ax.set(xlabel="Ks")
    ax.set(ylabel="# paralogs")
    width=wgd_xs[2]-wgd_xs[1]
    plt.bar(ssd_xs, ssd_ys, width=width, color="lightgray", label="ssd")
    plt.bar(wgd_xs, wgd_ys, width=width,
            color=config.color_blind_friendly_color_cycle_analogs['gray'], label="wgd")
    #plt.bar(ssds_xs_to_subtract, ssds_ys_to_subtract, width=0.001, color="lightgray", label="ssd under wgd")

    gaussian_fit_curve_ys1, xs_for_wgd, gaussian_goodness_of_fit = \
            curve_fitting.fit_curve_to_xs_and_ys(wgd_xs, wgd_ys, curve_fitting.wgd_gaussian)

    lognorm_fit_curve_ys1, xs_for_wgd,lognorm_goodness_of_fit =  \
            curve_fitting.fit_curve_to_xs_and_ys(wgd_xs, wgd_ys, curve_fitting.wgd_lognorm )


    if lognorm_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(lognorm_goodness_of_fit.RMSE,4))
            pser_str= str(round(lognorm_goodness_of_fit.pearsons_corr_result[2],4))
            plt.plot(xs_for_wgd,lognorm_fit_curve_ys1,
                 color=config.color_blind_friendly_color_cycle_analogs['blue'], linestyle=':',
                     linewidth=linewidth,
                     label="Lognorm fit \n(RMSE={0})".format(rmse_str,pser_str))

    if gaussian_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(gaussian_goodness_of_fit.RMSE,4))
            pser_str= str(round(gaussian_goodness_of_fit.pearsons_corr_result[2],4))
            plt.plot(xs_for_wgd,gaussian_fit_curve_ys1,
                 color=config.color_blind_friendly_color_cycle_analogs['red'], linestyle=':',
                     linewidth=linewidth,
                     label="Gaussian fit \n(RMSE={0})".format(rmse_str,pser_str))

    '''
    if lognorm_fit_curve_ys1 and (hist_maximum>0):
            plt.scatter(lognorm_goodness_of_fit.cm,-0.05*maxY,
                 color='darkgreen', marker='o', s=100, label="cm",)

            plt.scatter(lognorm_goodness_of_fit.mode,-0.05*maxY,
                 color='cyan', marker='^', s=80, label="mode")

    plt.scatter(wgd_maxima[0], -0.05 * maxY,
                color='green', marker='^', s=80, label="wgd abs max")

    plt.scatter(wgd_max_d[0], -0.05 * maxY,
                color='k', marker='^', s=80, label="wgd max derivative")
    '''
    #plt.bar(ssd_xs, ssd_ys, width=0.001, color="lightgray", label="ssd")
    #plt.bar(wgd_xs, wgd_ys, width=0.001, color="gray", label="wgd")
    #plt.bar(ssds_xs_to_subtract, ssds_ys_to_subtract, width=0.001, color="lightgray", label="ssd under wgd")

    plt.legend()
    plt.savefig(out_file_name)
    plt.close()

    fit_data = run_metrics(polyploid_name[0:4], polyploid_name,
                           SPC_time_MYA, WGD_time_MYA,
                           raw_cm, raw_x_value_of_ymax,
                           wgd_maxima[0],wgd_max_d[0],
                           lognorm_goodness_of_fit, gaussian_goodness_of_fit)

    return fit_data

