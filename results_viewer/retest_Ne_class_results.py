import math
import os
import unittest
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np


class Ne_TestLogisticalRegression(unittest.TestCase):

    def test_Ne_regression_on_high_vs_low_N_data(self):
        # color scheme:
        # High NE - allo color
        # Medium NE - light blue
        # low NE - gray
        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        xlabel = "SPC time"
        ylabel = "ln(metric3)"
        medium_color = 'blue'
        # Note, the input file is created by "metric3.csv"
        # created by "test_parse_agg_results" in the batch_aggregator.py

        cutoff = 100  # MY (change to 60 or 100 as desired)
        sims, specs, true_category, data = get_highN_vs_lowN_truth(cutoff)

        all_highs= [data[i] for i in range(0,len(data)) if true_category[i]=="High" ]
        min_highs=min(all_highs)
        print("min of high set",min_highs )

        all_meds= [data[i] for i in range(0,len(data)) if true_category[i]=="Medium" ]
        max_meds=max(all_meds)
        print("max of med set", max_meds)

        target = [1 if m == "High" else 0 for m in true_category]
        target_as_array = np.array(target)
        data_as_array = np.array([data]).transpose()
        X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33,
                                                            random_state=44)
        clf = LogisticRegression(penalty='l2', C=0.1)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy1 = metrics.accuracy_score(y_test, y_pred)
        accuracy_string1 = str(round(accuracy1, 4))

        y_pred2 = clf.predict(data_as_array)
        accuracy2 = metrics.accuracy_score(target_as_array, y_pred2)
        accuracy_string2 = str(round(accuracy2, 4))

        ROC_plot_base_name = "Retest Ne classification"
        #make_basic_ROC_plot(ROC_plot_base_name, clf, out_folder, X_test, y_test, accuracy_string1)
        make_basic_ROC_plot(ROC_plot_base_name, clf, out_folder, data_as_array, target_as_array, accuracy_string2)


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
    plt.savefig(plot_file, dpi=350)
    plt.close()

def get_highN_vs_lowN_truth(cutoff):

    out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
    file1="metric3.csv"
    full_path1 = os.path.join(out_folder, file1)

    spc_xs,metric,wgd_sims = read_xs_ys_csv(full_path1)
    low_N_sims= ["Auto","sim37_N0p1", "sim37_N1"]
    med_N_sims=["sim37_N5"]#, "sim35_log"]
    high_N_sims=["sim36_N10","sim37_N20"]

    sims=[]
    specs=[]
    category=[]
    data=[]

    #cutoff=60
    for i in range(0,len(wgd_sims)):
        sim = wgd_sims[i]
        spc_time=spc_xs[i]

        if spc_time >= cutoff:
            continue
        if metric[i] <= 0:
            #either its machine error, or a really bad fit.
            #either way, we shouldnt throw out the data point..
            metric[i]=abs(metric[i])

        true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
        sims.append(sim)
        specs.append(spc_time)
        category.append(true_category)
        data.append(math.log(metric[i]))

    return sims,specs,category,data


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