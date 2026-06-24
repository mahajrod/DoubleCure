#!/usr/bin/env python
__author__ = 'mahajrod'
import argparse
from copy import deepcopy
from pathlib import Path
from functools import partial

import pandas as pd

from RouToolPa.GeneralRoutines import FileRoutines

AGP_MAIN_COLUMNS_LIST = ["scaffold", "start", "end", "part_number", "part_type",
                             "part_id/gap_length", "part_start/gap_type", "part_end/linkage",
                             "orientation/evidence"]

AGP_GAP_TAG_LIST = ["U"]

RECOGNIZED_SEX_CHR_TAGS = ["W", "Z", "X", "Y"]

for sex_chr in ["W", "Z", "X", "Y"]:
    RECOGNIZED_SEX_CHR_TAGS += [f"{sex_chr}{i}" for i in range(1,51)]

RECOGNIZED_CURATION_TAGS = ["HAP1", "HAP2", "Haplotig", "Painted", "W", "Unloc", "Z", "X", "Y"]


def copy_tag_lambdda (s, tag):
    #print(s)
    return  tag in s

def parse_pretext_agp(pretext_agp_file):
    AGP_MAIN_COLUMNS_LIST = ["scaffold", "start", "end", "part_number", "part_type",
                             "part_id/gap_length", "part_start/gap_type", "part_end/linkage",
                             "orientation/evidence"]
    main_columns_list = []
    tag_set = set()
    with FileRoutines.metaopen(str(pretext_agp_file), "r") as agp_fd:
        for line in agp_fd:
            if line[0] == "#":
                continue
            line_list = line.strip().split("\t")
            main_columns_list.append(line_list[:9])
            if len(line_list) > 9:
                main_columns_list[-1].append(line_list[9:])
                tag_set = tag_set | set(line_list[9:])
            else:
                main_columns_list[-1].append([])

    agp_df = pd.DataFrame.from_records(main_columns_list, columns=AGP_MAIN_COLUMNS_LIST + ["TAG_TMP"])
    tag_columns = sorted(list(tag_set))
    for tag in tag_columns:
        agp_df[f"CT_{tag}"] = agp_df["TAG_TMP"].apply(partial(copy_tag_lambdda, tag=tag))
        agp_df[f"CT_{tag}"] = agp_df[f"CT_{tag}"].astype("boolean")

    if "CT_Haplotig" not in agp_df.columns:
        agp_df["CT_Haplotig"] = False
        agp_df["CT_Haplotig"] = agp_df["CT_Haplotig"].astype("boolean")
        tag_columns.append("Haplotig")
    if "CT_Unloc" not in agp_df.columns:
        agp_df["CT_Unloc"] = False
        agp_df["CT_Unloc"] = agp_df["CT_Unloc"].astype("boolean")
        tag_columns.append("Unloc")

    tag_columns = sorted(tag_columns)
    tag_columns = list(map(lambda s: f"CT_{s}", tag_columns ))
    return agp_df[AGP_MAIN_COLUMNS_LIST + tag_columns].set_index("scaffold")

def add_original_length_to_agp_df(row, len_df): # Row level function
    if row.iloc[3] in AGP_GAP_TAG_LIST:
        return pd.NA
    else:
        return len_df.loc[row.iloc[4], "length"]


def check_if_unsplit(row, len_df): # Row level function, checks if original scaffold was uncut. Returns 'False' for gap
    if row.iloc[3] in AGP_GAP_TAG_LIST: #part_type, check for gaps
        return pd.NA
    if int(row.iloc[5]) == 1:  #"part_start/gap_type"

        if len_df.loc[row.iloc[4], "length"] == int(row.iloc[6]):
            return True
        else:
            return False
    else:
        return False

def add_singleton_tag(agp_df): # DF level function

    part_number_df = agp_df[["start"]].groupby(by="scaffold").count()
    part_number_df.columns = pd.Index(["parts"])

    agp_df["PT_singleton"] = False
    agp_df["PT_singleton"] = agp_df["PT_singleton"].astype("boolean")
    agp_df.loc[agp_df.index.isin(part_number_df[part_number_df["parts"] == 1].index), "PT_singleton"] = True
    agp_df.loc[agp_df["part_type"].isin(AGP_GAP_TAG_LIST), "PT_singleton"] = pd.NA
    #return agp_df


def check_if_unprocessed(agp_df, curation_tag_set):
    unprocessed_bool_sr = agp_df["PT_unsplit"] & agp_df["PT_singleton"]

    for curation_tag in curation_tag_set:
        unprocessed_bool_sr &= ~agp_df[curation_tag]

    return unprocessed_bool_sr

def check_if_part_of_sex_chromosome(agp_df): # DF level function
    present_sex_chr_tag_list = []
    agp_df["PT_sex_chr"] = pd.NA
    agp_df["PT_sex_chr"] = agp_df["PT_sex_chr"].astype("object")
    for sex_chr in RECOGNIZED_SEX_CHR_TAGS:
        sex_chr_tag = f"CT_{sex_chr}"
        if sex_chr_tag in agp_df.columns:
            present_sex_chr_tag_list.append(sex_chr_tag)
            agp_df.loc[ agp_df[sex_chr_tag], "PT_sex_chr"] = sex_chr

    if not present_sex_chr_tag_list:
        agp_df["PT_part_of_sex_chr"] = False
        agp_df["PT_part_of_sex_chr"] = agp_df["PT_part_of_sex_chr"].astype("object")
    else:
        agp_df["PT_part_of_sex_chr"] = agp_df[present_sex_chr_tag_list[0]]
        agp_df["PT_part_of_sex_chr"] = agp_df["PT_part_of_sex_chr"].astype("object")
        if len(present_sex_chr_tag_list) > 1:
            for i in range(1, len(present_sex_chr_tag_list)):
                agp_df["PT_part_of_sex_chr"] |= agp_df[present_sex_chr_tag_list[i]]


def detect_haplotypes(agp_df):
    return set(list(map(lambda s: s.split(HAPLOTYPE_LABEL_SEPARATOR)[0],
             set(agp_df["part_id/gap_length"].loc[~agp_df["part_type"].isin(AGP_GAP_TAG_LIST)]))))


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--agp_folder", action="store", dest="agp_folder", required=True,
                    help="Folder with .agp files. Required.")
parser.add_argument("-l", "--len_file", action="store", dest="len_file", required=True,
                    help="Two-column tab-separated file with lengths of scaffolds. "
                         "The first column must contain scaffold ids, the second - its lengths, respectively. Required.")
parser.add_argument("-y", "--haplotype_prefix", action="store", dest="haplotype_prefix", default="hap",
                    help="Prefix of scaffold ids used to label haplotypes."
                         "Scaffold ids must follow this template: <haplotype_prefix><haplotype_number><separator><scaffold_subid>."
                         "For example, hap1.scaf256. "
                         "Default: 'hap'")
parser.add_argument("-s", "--haplotype_separator", action="store", dest="haplotype_separator", default=".",
                    help="Separator used to divide haplotype and scaffold subid within scaffold ids."
                         "Scaffold ids must follow this template: <haplotype_prefix><haplotype_number><separator><scaffold_subid>. "
                         "For example, hap1.scaf256. "
                         "Default: '.'")
parser.add_argument("-l", "--len_file", action="store", dest="len_file", required=True,
                    help="Two-column tab-separated file with lengths of scaffolds. "
                         "The first column must contain scaffold ids, the second - its lengths, respectively. Required.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files. Required.")

args = parser.parse_args()

agp_folder = args.agp_folder
len_file = args.len_file
output_prefix = args.output_prefix

HAPLOTYPE_PREFIX = args.haplotype_prefix
HAPLOTYPE_LABEL_SEPARATOR = args.haplotype_separator

#---- Parse agp files and add tags ----
agp_folder_path = Path(agp_folder)
agp_filelist = list(agp_folder_path.glob("*.agp"))

agp_filelist_dict = {f"AGP_{i+1}": agp_filelist[i] for i in range(0, len(agp_filelist))}

agp_df_dict = { agp_id: parse_pretext_agp(agp_filelist_dict[agp_id]) for agp_id in agp_filelist_dict}

agp_id_list = list(agp_df_dict.keys())

len_df = pd.read_csv(len_file, sep="\t", header=None, index_col=0, names=["scaffold", "length"])

all_curation_tags_set = set()

for agp_id in agp_id_list: # get curation tags from all curation files
    all_curation_tags_set |= set(agp_df_dict[agp_id].columns)

all_curation_tags_set -= set(AGP_MAIN_COLUMNS_LIST)
all_curation_tags_list = list(all_curation_tags_set)

number_of_tags = len(all_curation_tags_set)


for agp_id in agp_id_list: # addd missing tags to all files
    for curation_tag in all_curation_tags_set:
        if curation_tag not in list(agp_df_dict[agp_id].columns):
            agp_df_dict[agp_id][curation_tag] = False

for agp_id in agp_id_list:
    agp_df_dict[agp_id]["PT_original_length"] = agp_df_dict[agp_id].apply(partial(add_original_length_to_agp_df, len_df=len_df), axis=1)
    agp_df_dict[agp_id]["PT_unsplit"] = agp_df_dict[agp_id].apply(partial(check_if_unsplit, len_df=len_df), axis=1).astype("boolean")
    add_singleton_tag(agp_df_dict[agp_id])
    agp_df_dict[agp_id]["PT_unprocessed"] = check_if_unprocessed(agp_df_dict[agp_id], all_curation_tags_set)
    check_if_part_of_sex_chromosome(agp_df_dict[agp_id])
    agp_df_dict[agp_id]["PT_cut_segment"] = (~agp_df_dict[agp_id]["PT_unsplit"]) & (~agp_df_dict[agp_id]["CT_Painted"])

haplotype_list = detect_haplotypes(agp_df_dict[agp_id_list[0]])
#----

#---- extract rows for processed_scaffolds ----
agp_df_processed_dict = {}
agp_df_unprocessed_dict = {}

for agp_index in range(0, len(agp_id_list)):
    agp_id = agp_id_list[agp_index]
    # extract processed (modified or merged) scaffolds
    agp_df_processed_dict[agp_id] = agp_df_dict[agp_id].loc[~agp_df_dict[agp_id]["PT_unprocessed"]].reset_index(drop=False, inplace=False)
    # label individual dataframes
    agp_df_processed_dict[agp_id]["agp_id"] = agp_id
    # make 'agp_id' the first column
    agp_df_processed_dict[agp_id] = agp_df_processed_dict[agp_id][["agp_id"] + list(agp_df_processed_dict[agp_id].columns[:-1])]
    # extract processed (modified or merged) scaffolds
    agp_df_unprocessed_dict[agp_id] = agp_df_dict[agp_id].loc[agp_df_dict[agp_id]["PT_unprocessed"]].reset_index(drop=False, inplace=False)
    # label individual dataframes
    agp_df_unprocessed_dict[agp_id]["agp_id"] = agp_id
    # make 'agp_id' the first column
    agp_df_unprocessed_dict[agp_id] = agp_df_unprocessed_dict[agp_id][["agp_id"] + list(agp_df_unprocessed_dict[agp_id].columns[:-1])]
    #extract ids of scaffolds that were not processed in any of dataframes
    if agp_index == 0:
        unprocessed_scaffolds_set = set(agp_df_unprocessed_dict[agp_id]["part_id/gap_length"])
    else:
        unprocessed_scaffolds_set &= set(agp_df_unprocessed_dict[agp_id]["part_id/gap_length"])
    print(f"{agp_id}\t{len(agp_df_processed_dict[agp_id])}\t{len(agp_df_unprocessed_dict[agp_id])}\t{agp_filelist_dict[agp_id]}")

merged_processed_df = pd.concat(agp_df_processed_dict.values(), axis=0).reset_index(drop=True)
unprocessed_scaffolds_df = agp_df_unprocessed_dict[agp_id_list[0]][agp_df_unprocessed_dict[agp_id_list[0]]["part_id/gap_length"].isin(unprocessed_scaffolds_set)]
#----

#---- check for segments that were cut from contigs, but not included in scaffolds ----
cut_segment_df = merged_processed_df[(merged_processed_df["PT_cut_segment"]) & (merged_processed_df["part_type"] != "U")]
if not cut_segment_df.empty:
    cut_segment_df["filename"] = cut_segment_df["agp_id"].apply(lambda agp_id: agp_filelist_dict[agp_id].name)
    cut_segment_df = cut_segment_df.set_index("scaffold")
    cut_segment_df[["part_id/gap_length", "agp_id", "filename"]].to_csv(f"{output_prefix}.cut_segments.info", sep="\t", header=True, index=False)
cut_segment_df.to_csv(f"{output_prefix}.cut_segments.extended_info", sep="\t" , header=True, index=True)
#----

#---- check for scaffolds curated multiple times and report them separately ----

no_gaps_merged_processed_df = merged_processed_df[~merged_processed_df["part_type"].isin(AGP_GAP_TAG_LIST)]
tmp_no_gaps_merged_processed_df = deepcopy(no_gaps_merged_processed_df)
curated_times_df = no_gaps_merged_processed_df[["agp_id", "part_id/gap_length"]].sort_values(by="part_id/gap_length").drop_duplicates().groupby(by="part_id/gap_length").count()
multiple_curated_df = curated_times_df[curated_times_df["agp_id"] > 1]
multiple_curated_records_df = []

if len(multiple_curated_df) > 0:
    print("WARNING!!! Some scaffolds were processed twice or more times")
    multiply_curated_main_info_df = no_gaps_merged_processed_df[["agp_id", "part_id/gap_length"]].sort_values(by="part_id/gap_length").drop_duplicates().set_index("part_id/gap_length").loc[multiple_curated_df.index]
    multiply_curated_main_info_df["filename"] = multiply_curated_main_info_df["agp_id"].apply(lambda s: agp_filelist_dict[s] )
    multiply_curated_main_info_df.to_csv(f"{output_prefix}.main_info.tsv", sep="\t", header=True, index=True)

    multiple_curated_scaffolds_ids = set(multiple_curated_df.index)
    multiple_curated_scaffolds_id_df = no_gaps_merged_processed_df[no_gaps_merged_processed_df["part_id/gap_length"].isin(multiple_curated_scaffolds_ids)][["agp_id", "scaffold"]].drop_duplicates().set_index("agp_id")
    multiple_curated_scaffolds_id_df.to_csv(f"{output_prefix}.double_curated.scaffolds.ids.tsv", sep="\t", header=True, index=False)

    for agp_id in multiple_curated_scaffolds_id_df.index:
        for sup_scaffold_id in multiple_curated_scaffolds_id_df["scaffold"].loc[[agp_id]]:
            print(agp_id, sup_scaffold_id)
            double_curated_selection_series = (merged_processed_df["agp_id"] == agp_id) & (merged_processed_df["scaffold"] == sup_scaffold_id)
            multiple_curated_records_df.append(merged_processed_df[double_curated_selection_series])
            merged_processed_df = merged_processed_df[~double_curated_selection_series]
            double_curated_selection_series = (no_gaps_merged_processed_df["agp_id"] == agp_id) & (no_gaps_merged_processed_df["scaffold"] == sup_scaffold_id)
            no_gaps_merged_processed_df = no_gaps_merged_processed_df[~double_curated_selection_series]

    multiple_curated_records_df = pd.concat(multiple_curated_records_df, axis=0)
    multiple_curated_records_df.to_csv(f"{output_prefix}.double_curated.scaffolds.tsv", sep="\t", header=True, index=False)

# remove haplotigs
haplotig_df = no_gaps_merged_processed_df[no_gaps_merged_processed_df["CT_Haplotig"]]
painted_haplotig_df = haplotig_df[haplotig_df["CT_Painted"]]
if len(painted_haplotig_df) > 0:
    print("WARNING!!! Some of haplotig scafffolds are labeled as 'Painted'")
painted_haplotig_df.to_csv(f"{output_prefix}.haplotigs.painted.tsv", sep="\t", header=True, index=False)
haplotig_df.to_csv(f"{output_prefix}.haplotigs.tsv", sep="\t", header=True, index=False)

no_gaps_merged_processed_df = no_gaps_merged_processed_df[~no_gaps_merged_processed_df["CT_Haplotig"]]

if not multiple_curated_records_df.empty:
    contig_id_to_file_df = multiple_curated_records_df[multiple_curated_records_df["part_type"] != "U"][["agp_id","part_id/gap_length"]].drop_duplicates().set_index("part_id/gap_length")
    tmpfg_df = contig_id_to_file_df.groupby("part_id/gap_length").count()
    double_curated_contigs_df = tmpfg_df[tmpfg_df["agp_id"] > 1]
    double_curated_contigs_info_df = contig_id_to_file_df.loc[double_curated_contigs_df.index]
    double_curated_contigs_info_df["filename"] = double_curated_contigs_info_df["agp_id"].apply(lambda agp_id: agp_filelist_dict[agp_id].name)
    double_curated_contigs_info_df.to_csv(f"{output_prefix}.double_curated.contigs.info", sep="\t", header=True, index=True)

#----

#TODO: rewite to maintain polyploid genomes

def split_scaffolds(df): # function to handle df provided by groupby(by=["agp_id", "scaffold"])

    for sex_chr in RECOGNIZED_SEX_CHR_TAGS:
        sex_chr_tag = f"CT_{sex_chr}"
        if sex_chr_tag in df.columns:
            df.loc[df[sex_chr_tag], "scaffold"] = "chr" + sex_chr
    number_of_unlocated_scaffolds = len(df.loc[df["CT_Unloc"]])
    if number_of_unlocated_scaffolds > 0:
        #print(df.loc[df["CT_Unloc"], "scaffold"] + [f"_unloc{i}" for i in range(1, number_of_unlocated_scaffolds + 1)])
        df.loc[df["CT_Unloc"], "scaffold"] = df.loc[df["CT_Unloc"], "scaffold"] + [f"_unloc{i}" for i in range(1, number_of_unlocated_scaffolds + 1)]
    df.loc[df["CT_HAP1"], "scaffold"] = "hap1." + df.loc[df["CT_HAP1"], "scaffold"]
    df.loc[df["CT_HAP2"], "scaffold"] = "hap2." + df.loc[df["CT_HAP2"], "scaffold"]

    return df

def correct_coordinates_and_insert(df): # function to handle df provided by groupby(by=["agp_id", "scaffold"]). Must run after application of split_scaffolds_df
    if len(df) == 1:
        df["start"] = 1
        df["end"] = df["part_end/linkage"].astype(int) - df["part_start/gap_type"].astype(int) + 1
        df["part_number"] = 1
        return df
    else:
        row_list = []
        row_number = len(df)
        last_row_index = row_number - 1
        for row_index in range(0, row_number):
            row_df = df.iloc[[row_index]]
            if row_index == 0:
                row_df["start"] = 1
                row_df["part_number"] = 1
            else:
                row_df["start"] = row_list[-1]["end"].iloc[0]  + 1
                row_df["part_number"] = row_list[-1]["part_number"].iloc[0] + 1
            row_df["PT_original_length"] = row_df["PT_original_length"].astype('Int64')
            row_df["end"] = row_df["start"].astype(int) + row_df["part_end/linkage"].astype(int) - row_df["part_start/gap_type"].astype(int)
            row_list.append(row_df)

            if row_index < last_row_index: # add gap
                gap_row_df = row_df.copy(deep=True)
                gap_row_df["start"] = row_df["end"] + 1
                gap_row_df["end"] = gap_row_df["start"] + 99
                gap_row_df["part_number"] = row_df["part_number"] + 1
                gap_row_df["part_type"] = "U"
                gap_row_df["part_id/gap_length"] = "100"
                gap_row_df["part_start/gap_type"] = "scaffold"
                gap_row_df["part_end/linkage"] = "yes"
                gap_row_df["orientation/evidence"] = "proximity_ligation"
                gap_row_df[all_curation_tags_list] = False
                gap_row_df["PT_original_length"] = pd.NA
                gap_row_df[["PT_unsplit", "PT_singleton", "PT_unprocessed"]] = pd.NA

                row_list.append(gap_row_df)
        return pd.concat(row_list, axis=0)

no_gaps_merged_processed_split_scaffolds_df = no_gaps_merged_processed_df.groupby(by=["agp_id", "scaffold"]).apply(split_scaffolds).reset_index(inplace=False, drop=True)
#print(no_gaps_merged_processed_split_scaffolds_df[~no_gaps_merged_processed_split_scaffolds_df["scaffold"].apply(lambda s: (s[:4] == "hap1") or (s[:4] == "hap2"))]) # TODO add check for missing tags
prefinal_unsorted_df = no_gaps_merged_processed_split_scaffolds_df.groupby(by=["agp_id", "scaffold"]).apply(correct_coordinates_and_insert).reset_index(inplace=False, drop=True)

#create syn dict for renaming autosomes according the length
#TODO: add handling a case if processed scaffold from hap2 doesnt have a homolog
# select HAP1 processed scaffolds without unlocated contigs
hap1_bool_series = prefinal_unsorted_df["CT_HAP1"] & (~prefinal_unsorted_df["CT_Unloc"])
unlock_bool_series = prefinal_unsorted_df["CT_Unloc"].copy(deep=True)
hap2_bool_series = prefinal_unsorted_df["CT_HAP2"] & (~prefinal_unsorted_df["CT_Unloc"])
# remove_sex_chromosomes

for sex_chr in RECOGNIZED_SEX_CHR_TAGS:
    sex_chr_tag = f"CT_{sex_chr}"
    if sex_chr_tag in prefinal_unsorted_df.columns:
        hap1_bool_series &= ~prefinal_unsorted_df[sex_chr_tag]
        hap2_bool_series &= ~prefinal_unsorted_df[sex_chr_tag]
        unlock_bool_series &= ~prefinal_unsorted_df[sex_chr_tag]

hap1_syn_df = prefinal_unsorted_df[hap1_bool_series][["agp_id", "scaffold", "end"]].groupby(by=["agp_id", "scaffold"]).max().sort_values(by="end", ascending=False) #.reset_index(inplace=False, drop=False)
hap1_syn_df["scaffold_syn"] = [f"hap1.aut{i}" for i in range(1, len(hap1_syn_df) + 1) ]
hap2_syn_df = prefinal_unsorted_df[hap2_bool_series][["agp_id", "scaffold", "end"]].groupby(by=["agp_id", "scaffold"]).max().reset_index(drop=False)
hap2_syn_df["hap1_homolog"] = hap2_syn_df["scaffold"].apply(lambda s: s.replace("hap2", "hap1"))
hap2_syn_df.set_index(["agp_id", "hap1_homolog"], inplace=True)
hap2_syn_df["hap1_homolog_syn"] = hap1_syn_df["scaffold_syn"].loc[hap2_syn_df.index]
hap2_syn_df["scaffold_syn"] = hap2_syn_df["hap1_homolog_syn"].apply(lambda s: s.replace("hap1", "hap2"))
hap2_syn_df = hap2_syn_df.reset_index(drop=False)[["agp_id", "scaffold", "end", "scaffold_syn"]].set_index(["agp_id", "scaffold"])

syn_df = pd.concat([hap1_syn_df, hap2_syn_df], axis=0)

unlocated_renaming_df = prefinal_unsorted_df[unlock_bool_series][["agp_id", "scaffold", "end"]]
unlocated_renaming_df["prefix"] = unlocated_renaming_df["scaffold"].apply(lambda s: s.split("_unloc")[0])
unlocated_renaming_df["suffix"] = unlocated_renaming_df["scaffold"].apply(lambda s: "_" + s.split("_")[-1])
unlocated_renaming_df.set_index(["agp_id", "prefix"], inplace=True)

unlocated_renaming_df["prefix_syn"] = syn_df["scaffold_syn"].loc[unlocated_renaming_df.index]
unlocated_renaming_df["scaffold_syn"] = unlocated_renaming_df["prefix_syn"] + unlocated_renaming_df["suffix"]
unlocated_renaming_df.reset_index(drop=False, inplace=True)
unlocated_renaming_df = unlocated_renaming_df[["agp_id", "scaffold", "end", "scaffold_syn" ]].set_index(["agp_id", "scaffold"])

syn_df = pd.concat([syn_df, unlocated_renaming_df], axis=0)
syn_df.reset_index(drop=False, inplace=True)
syn_df["full_scaffold_id"] = syn_df["agp_id"] + "@" +  syn_df["scaffold"]
rename_dict = syn_df[["full_scaffold_id", "scaffold_syn" ]].set_index("full_scaffold_id")["scaffold_syn"].to_dict()

# rename processed scaffolds by syn dict
prefinal_unsorted_df["full_scaffold_id"] = prefinal_unsorted_df["agp_id"] + "@" +  prefinal_unsorted_df["scaffold"]

prefinal_unsorted_df.loc[prefinal_unsorted_df["PT_part_of_sex_chr"].fillna(False), "full_scaffold_id"] =   prefinal_unsorted_df["scaffold"]

prefinal_unsorted_df["scaffold_syn"] = prefinal_unsorted_df["full_scaffold_id"].copy(deep=True)
prefinal_unsorted_df["scaffold_syn"] = prefinal_unsorted_df["scaffold_syn"].replace(rename_dict)
#prefinal_unsorted_df[list(prefinal_unsorted_df.columns[-2:]) + list(prefinal_unsorted_df.columns[:-2])]
prefinal_unsorted_df["scaffold"] = prefinal_unsorted_df["scaffold_syn"]

#rename cut segments
prefinal_unsorted_cut_segments_df = prefinal_unsorted_df[prefinal_unsorted_df["PT_cut_segment"]==True][["scaffold", "part_id/gap_length"]]
#safety check, all cut segments should have only one contig.
contig_counts_per_cut_segment_df = prefinal_unsorted_cut_segments_df.groupby("scaffold").count()
cut_segments_multiple_contigs_df = contig_counts_per_cut_segment_df[contig_counts_per_cut_segment_df["part_id/gap_length"] > 1]

if not cut_segments_multiple_contigs_df.empty:
    print("WARNING!!! Some cut segments contain multiple contigs which should be impossible. Something have gone wrong!")

correct_cut_segments_index = contig_counts_per_cut_segment_df[contig_counts_per_cut_segment_df["part_id/gap_length"] == 1].index

correct_cut_segments_df = prefinal_unsorted_cut_segments_df.set_index("scaffold").loc[correct_cut_segments_index]
correct_cut_segments_df["prefix"]  = correct_cut_segments_df["part_id/gap_length"].apply(lambda s: s.split(HAPLOTYPE_LABEL_SEPARATOR)[0])
correct_cut_segments_df["suffix"] = [f"cutseg{i}" for i in range(1, len(correct_cut_segments_df) + 1) ]

correct_cut_segments_df["scaffold_syn"] = correct_cut_segments_df["prefix"] + HAPLOTYPE_LABEL_SEPARATOR + correct_cut_segments_df["suffix"]

cut_segments_rename_dict = correct_cut_segments_df["scaffold_syn"].to_dict()
prefinal_unsorted_df["scaffold"] = prefinal_unsorted_df["scaffold"].replace(cut_segments_rename_dict )

prefinal_unsorted_df = prefinal_unsorted_df[AGP_MAIN_COLUMNS_LIST + all_curation_tags_list + ["PT_original_length","PT_unsplit","PT_singleton", "PT_unprocessed", "PT_sex_chr", "PT_part_of_sex_chr"]]

#rename unprocessed scaffolds by length
unprocessed_scaffolds_df = unprocessed_scaffolds_df.sort_values(by=["end"], ascending=False)
unprocessed_scaffolds_df["suffix"] = [f"unloc{i}" for i in range(1, len(unprocessed_scaffolds_df) + 1) ]
unprocessed_scaffolds_df["prefix"] = unprocessed_scaffolds_df["part_id/gap_length"].apply(lambda s: s.split(HAPLOTYPE_LABEL_SEPARATOR)[0])
unprocessed_scaffolds_df["scaffold"] = unprocessed_scaffolds_df["prefix"] + HAPLOTYPE_LABEL_SEPARATOR + unprocessed_scaffolds_df["suffix"]
unprocessed_scaffolds_df["CT_HAP1"] = unprocessed_scaffolds_df["prefix"] == "hap1"
unprocessed_scaffolds_df["CT_HAP2"] = unprocessed_scaffolds_df["prefix"] == "hap2"
unprocessed_scaffolds_df = unprocessed_scaffolds_df[AGP_MAIN_COLUMNS_LIST + all_curation_tags_list + ["PT_original_length","PT_unsplit","PT_singleton", "PT_unprocessed", "PT_sex_chr", "PT_part_of_sex_chr"]]

final_df = pd.concat([prefinal_unsorted_df, unprocessed_scaffolds_df])
final_df.to_csv(f"{output_prefix}.final.tsv", sep="\t", header=True, index=False)
final_df[AGP_MAIN_COLUMNS_LIST].to_csv(f"{output_prefix}.final.agp", sep="\t", header=False, index=False)

