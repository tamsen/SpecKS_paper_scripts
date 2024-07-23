import pandas as pd


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

def read_Ks_csv_by_dup_type(csv_file):

    ks_results = []
    ks_results_by_type={"WGD":[],"WGD-SSD":[],"SSD-SSD":[]}
    with open(csv_file, "r") as f:

        reading_header=True
        while True:
            line = f.readline()
            if "ersion" in line:
                continue
            if "Git" in line:
                continue
            if "leaf names" in line:
                continue
            if not line:
                break
            if len(line)==0:
                break
            if reading_header:
                reading_header=False
                continue
            data = line.split(",")
            #print(data)
            gene_names=data[1].split()
            count_SSDs=sum([1 for n in gene_names if 'duplicate' in n])
            ks_value=float(data[2])
            ks_results.append(ks_value)
            if count_SSDs ==0:
                ks_results_by_type['WGD'].append(ks_value)
            elif count_SSDs == 1:
                ks_results_by_type['WGD-SSD'].append(ks_value)
            else:
                ks_results_by_type['SSD-SSD'].append(ks_value)

    return ks_results_by_type