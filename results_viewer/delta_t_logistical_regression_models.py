
import os
import unittest
import numpy as np
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import config, two_d_colors



#https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
# https://stackoverflow.com/questions/28716241/controlling-the-threshold-in-logistic-regression-in-scikit-learn

class DeltaT_TestLinearRegression(unittest.TestCase):

    def test_delta_t_regression(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"

        #input file
        input_file_basename= "allopolyploid_index_truth_vs_predictions.csv"

        # Note, the input file is created by "test_auto_vs_allo_predictor,"
        # found in "wgd_time_predictor.py"
        # "genes_remaining_vs_wgd_time.csv" and "mode_vs_spec_time.csv"
        # which are made by "test_parse_agg_results" in the batch_aggregator.py

        csv_file=os.path.join(out_folder,input_file_basename)
        obs_data1, truth1, target1 = load_dt_data(csv_file, two_d_colors.low_dt_limit)
        colors = [two_d_colors.rgb_colors.nice_low_dt if m == 0 else
                  two_d_colors.rgb_colors.nice_tan for m in target1]
        linear_regression_threshold_color = 'black'
        likely_range_for_threshold = np.arange(0.8, 2, 0.001)
        custom_prob_thresholds =  np.arange(0, 1, 0.01)
        ROC_plot_base_name = "ROC for low vs medium & high ΔT classification"
        threshold_plot_title = "Lower Threshold applied to ΔT estimation"
        threshold_plot_labels_by_case={0:"true low ΔT",1:"true high&medium ΔT"}
        case0_color=config.auto_color
        xlabel = "true ΔT"
        ylabel = "inferred ΔT"
        linear_regression_threshold1,clf= do_LG_and_ROC_plots(ROC_plot_base_name, colors, obs_data1, likely_range_for_threshold,
                            custom_prob_thresholds,
                            linear_regression_threshold_color, case0_color,
                            out_folder, truth1, xlabel, ylabel, target1,
                            threshold_plot_title, threshold_plot_labels_by_case)




        obs_data2, truth2, target2 = load_dt_data(csv_file, two_d_colors.medium_dt_limit)
        colors = [two_d_colors.rgb_colors.nice_dkorange if m == 0
                  else two_d_colors.rgb_colors.nice_tan  for m in target2]
        linear_regression_threshold_color = 'black'
        likely_range_for_threshold = np.arange(30, 40, 0.001)
        custom_prob_thresholds = np.arange(0, 1, 0.01)
        ROC_plot_base_name = "ROC for low & medium vs high ΔT classification"
        threshold_plot_title = "Medium Threshold applied to ΔT estimation"
        threshold_plot_labels_by_case = {0: "true medium ΔT", 1: "true high&medium ΔT"}
        case0_color = config.auto_color
        xlabel = "true ΔT"
        ylabel = "inferred ΔT"
        linear_regression_threshold2,clf = do_LG_and_ROC_plots(ROC_plot_base_name, colors, obs_data2, likely_range_for_threshold,
                            custom_prob_thresholds,
                            linear_regression_threshold_color, case0_color,
                            out_folder, truth2, xlabel, ylabel, target2,
                            threshold_plot_title, threshold_plot_labels_by_case)

        #plot the two thresholds together

        target3=[target1[i]+target2 [i] for i in range(0,len(target2))]
        linear_regression_threshold_color1='black'
        linear_regression_threshold_color2='black'
        threshold_plot_title3 = "Low, medium and high ΔT classification"

        colors3 = [two_d_colors.rgb_colors.nice_tan if m == 2 else
                   two_d_colors.rgb_colors.nice_dkorange if m == 1 else
                   two_d_colors.rgb_colors.nice_low_dt for m in target3]

        plot_2_thresholds_against_data(colors3, obs_data2, linear_regression_threshold1,
                                            linear_regression_threshold2,
                                            linear_regression_threshold_color1,
                                            linear_regression_threshold_color2,
                                            out_folder, truth2, threshold_plot_title3)


def make_basic_ROC_plot(file_basename, clf, out_folder,
                        X_test, y_test, accuracy_string):

    y_pred_proba = clf.predict_proba(X_test)[::, 1]
    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_proba, drop_intermediate=False)
    auc = metrics.roc_auc_score(y_test, y_pred_proba)
    auc_string=str(round(auc, 4))
    print("thresholds:\t" + str(thresholds))
    title = file_basename.replace(".csv", " ROC plot")
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.plot(fpr, tpr, label="auc={0},accuracy={1}".format(auc_string,accuracy_string),
             color='k', marker='o')
    ax.set(xlabel="false positive rate")
    ax.set(ylabel="true positive rate")
    plt.legend(loc=4)
    ax.set(title=title)
    plt.tight_layout()
    plot_file = os.path.join(out_folder, title.replace(" ", "_") + ".png")
    plt.savefig(plot_file)
    plt.close()

def do_LG_and_ROC_plots(ROC_plot_base_name, colors, observed_data, likely_range_for_threshold,
                        custom_prob_thresholds,
                        linear_regression_threshold_color, case0_color, out_folder,
                        truth, xlabel, ylabel, target, threshold_plot_title, plot_labels_by_case):
    target_as_array = np.array(target)
    data_as_array = np.array([observed_data]).transpose()
    X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33,
                                                        random_state=44)
    clf = LogisticRegression(penalty='l2', C=0.1)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    accuracy = metrics.accuracy_score(y_test, y_pred)
    accuracy_string = str(round(accuracy, 4))
    make_basic_ROC_plot(ROC_plot_base_name, clf, out_folder, X_test, y_test, accuracy_string)


    linear_regression_threshold = get_linear_regession_metric_threshold(clf, likely_range_for_threshold)
    plot_threshold_against_data(colors, observed_data, linear_regression_threshold,
                                linear_regression_threshold_color, case0_color,
                                out_folder, threshold_plot_title,
                                xlabel, ylabel,
                                truth, plot_labels_by_case)
    print(linear_regression_threshold)
    return linear_regression_threshold,clf

def get_linear_regession_metric_threshold(clf, likely_range_for_threshold):
    model_predictions = clf.predict(np.array([likely_range_for_threshold]).transpose())
    model_prediction_change = [model_predictions[i + 1] - model_predictions[i] for i in
                               range(0, len(likely_range_for_threshold) - 1)]
    i_of_threshold = [i for i in range(0, len(model_prediction_change)) if model_prediction_change[i] == 1]
    linear_regression_threshold = 0.5 * (
                likely_range_for_threshold[i_of_threshold[0]] + likely_range_for_threshold[i_of_threshold[0] + 1])
    return linear_regression_threshold

def plot_threshold_against_data(colors, data, linear_regression_threshold,
                                linear_regression_threshold_color, case0_color,
                                out_folder, plot_title, xlabel, ylabel,
                                specs, plot_labels_by_case):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    labeled_0=False
    labeled_1=False
    alpha=0.2
    num_data_points=len(specs)
    for i in range(0,num_data_points):

        if colors[i]==case0_color and labeled_0==0:
            plt.scatter(specs[i], data[i], label=plot_labels_by_case[0], marker='o',
                        c=colors[i],alpha=alpha)
            labeled_0=True
        elif colors[i]!=case0_color and labeled_1==0:
            plt.scatter(specs[i], data[i], label=plot_labels_by_case[1], marker='o',
                        c=colors[i],alpha=alpha)
            labeled_1=True
        else:
            plt.scatter(specs[i], data[i], marker='o', c=colors[i],alpha=alpha)


    ax.axhline(y=linear_regression_threshold, color=linear_regression_threshold_color, linestyle='--',
               label="thresh"
                     + " ({0})".format(round(linear_regression_threshold, 4)))
    ax.set(xlabel=xlabel)
    ax.set(ylabel=ylabel)
    #plt.legend(loc=4)
    plt.legend()
    ax.set(title=plot_title + "\nn="+str(num_data_points))
    plot_file = os.path.join(out_folder, plot_title + ".png")
    plt.savefig(plot_file)
    plt.close()


#cvs data looks like this:
#sims	truth	prediction
#Allo6_S030W010.hist_Ks_hist_fit2.0_sim37_N0p1	20	19.4126670765897
def load_dt_data(csv_file, threshold):
    with open(csv_file, "r") as f:
        lines = f.readlines()

    #"low" is <
    estimate = []
    target = []
    truth = []
    for l in lines:
        if "truth" in l:
            continue
        splat=l.split(',')
        truth_i=float(splat[1])
        estimate_i=float(splat[2])
        estimate.append(estimate_i)
        truth.append(truth_i)

        if truth_i < threshold:
            target.append(0)
        else:
            target.append(1)

    return estimate,truth, target


def plot_2_thresholds_against_data(colors3, data, linear_regression_threshold1, linear_regression_threshold2,
                                   linear_regression_threshold_color1, linear_regression_threshold_color2,
                                   out_folder, specs, threshold_plot_title3):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    alpha = 0.2
    n=len(specs)
    for i in range(0, len(specs)):
        plt.scatter(specs[i], data[i], marker='o', c=colors3[i], alpha=alpha)
    ax.axhline(y=linear_regression_threshold1, color=linear_regression_threshold_color1, linestyle='--',
               label="lvm thresh"
                     + " ({0})".format(round(linear_regression_threshold1, 4)))
    ax.axhline(y=linear_regression_threshold2, color=linear_regression_threshold_color2, linestyle='--',
               label="mvh thresh"
                     + " ({0})".format(round(linear_regression_threshold2, 4)))
    ax.set(xlabel=" true ΔT ")
    ax.set(ylabel=" est. ΔT ")
    plt.legend()#loc=1)
    plt.tight_layout(pad=3)
    ax.set(title=threshold_plot_title3 + "\nn={}".format(n))
    plot_file = os.path.join(out_folder, threshold_plot_title3 + ".png")
    plt.savefig(plot_file)
    plt.close()

if __name__ == '__main__':
    unittest.main()
