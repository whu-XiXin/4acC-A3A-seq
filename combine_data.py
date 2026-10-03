import pandas as pd
import numpy as np
import os

file_exp1 = "/home/A3A.Ecoil.rep1/methylation_data/A3A.Ecoil.rep1.deduplicated.CX_report.dss.txt"
file_exp2 = "/home/A3A.Ecoil.rep2/methylation_data/A3A.Ecoil.rep2.deduplicated.CX_report.dss.txt"
file_con1 = "/home/Base.A3A.Ecoil.rep1/methylation_data/Base.A3A.Ecoil.rep2.deduplicated.CX_report.dss.txt"
file_con2 = "/home/Base.A3A.Ecoil.rep2/methylation_data/Base.A3A.Ecoil.rep2.deduplicated.CX_report.dss.txtss.txt"

temp_file_1 = "/home/bismark/temp_combined_1_rmNA.tsv"
temp_file_2 = "/home/bismark/temp_combined_2_rmNA.tsv"
final_output = "/home/bismark/4acC_A3A_seq_combine.txt"

TEST_MODE = False
nrows = 10000 if TEST_MODE else None


def add_prefix(df, prefix):
    protected = {"chr", "pos", "strand"}
    new_cols = []
    for c in df.columns:
        if c in protected:
            new_cols.append(c)
        else:
            new_cols.append(f"{prefix}_{c}")
    df.columns = new_cols
    return df


def read_with_methylation_rate(file_path):
    df = pd.read_csv(file_path, sep="\t", header=0, nrows=nrows)
    df["Methylation_rate"] = df["X"] / df["N"]
    return df


if not os.path.exists(temp_file_1):
    print("Merging exp1 and con1...")
    data_exp1 = read_with_methylation_rate(file_exp1)
    data_con1 = read_with_methylation_rate(file_con1)

    combined_1 = pd.merge(data_exp1, data_con1, on=["chr", "pos"], how="outer")
    combined_1 = add_prefix(combined_1, "CtoT-rep1")
    combined_1 = add_prefix(combined_1, "4acCtoT-rep1")

    combined_1.to_csv(temp_file_1, sep="\t", index=False)
    del data_exp1, data_con1, combined_1
    print(f"Saved to {temp_file_1}")
else:
    print(f"Temporary file {temp_file_1} already exists. Skipping merging step 1.")


if not os.path.exists(temp_file_2):
    print("Merging exp2 and con2...")
    data_exp2 = read_with_methylation_rate(file_exp2)
    data_con2 = read_with_methylation_rate(file_con2)

    combined_2 = pd.merge(data_exp2, data_con2, on=["chr", "pos"], how="outer")
    combined_2 = add_prefix(combined_2, "CtoT-rep2")
    combined_2 = add_prefix(combined_2, "4acCtoT-rep2")

    combined_2.to_csv(temp_file_2, sep="\t", index=False)
    del data_exp2, data_con2, combined_2
    print(f"Saved to {temp_file_2}")
else:
    print(f"Temporary file {temp_file_2} already exists. Skipping merging step 2.")


print("Final merging of two temporary files...")
combined_1 = pd.read_csv(temp_file_1, sep="\t")
combined_2 = pd.read_csv(temp_file_2, sep="\t")
data_combine = pd.merge(combined_1, combined_2, on=["chr", "pos"], how="outer")
del combined_1, combined_2

data_combine = data_combine.dropna()

print("Saving final result...")
data_combine.to_csv(final_output, sep="\t", index=False)
print(f"Final file saved to {final_output}")