#!/usr/bin/env python
__author__ = 'mahajrod'

from copy import deepcopy

import argparse

import pandas as pd

from RouToolPa.Parsers.AGP import CollectionAGP
from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fasta", action="store", dest="fasta", required=True,
                    help="Fasta file with original sequences. Required.")
parser.add_argument("-a", "--agp", action="store", dest="agp", required=True,
                    help="AGP file to correct. Required.")
parser.add_argument("-x", "--texel_size", action="store", dest="texel_size", required=True, type=int,
                    help="Size of the texel in the pretext map used to generate agp file. If you agp file is a unification of several files, use median or mean value. Required.")
parser.add_argument("-m", "--max_texel_dist", action="store", dest="max_texel_dist", default=2, type=int,
                    help="Maximal distance (in texels) from a breakpoint to the nearest gap. If there is a gap within such a distance, "
                         "coordinates of the breakpoint will be changed to the gap."
                         "Default: 2")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files. Required.")

args = parser.parse_args()

fasta_col = CollectionSequence(in_file=args.fasta, parsing_mode="parse")
fasta_col.get_stats_and_features()

agp_col = CollectionAGP(in_file=args.agp)

texel_size = args.texel_size
max_texel_dist = args.max_texel_dist
max_gap_dist = texel_size * max_texel_dist

output_prefix = args.output_prefix

#---- Add ids to gaps ----
fasta_col.gaps_bed.records["gap_id"] = list(map(lambda s: f"GAP{s}",
                                           range(0, len(fasta_col.gaps_bed.records))))
fasta_col.gaps_bed.records.to_csv(f"{output_prefix}.original_gaps.tsv", sep="\t", index=True, header=True)
#----

#---- Get processed parts ----
agp_parts_df = agp_col.get_parts_df(sort=True)

# Added columns:
# "full_sequence_length" - length of the contig from which the segment have originated
# "five_prime_is_start"  - five prime end (orientation of the segment within scaffold is not taken into account)
#                          of the segment corresponds to the five prime end of the original contig
# "three_prime_is_end"  - three prime end (orientation of the segment within scaffold is not taken into account)
#                          of the segment corresponds to the three prime end of the original contig
# "full_len_part"       - the segment is an unprocessed original contig

agp_parts_df["full_sequence_length"] = agp_parts_df["part_id/gap_length"].apply(lambda s: fasta_col.seq_lengths.loc[s,"length"])
agp_parts_df["five_prime_is_start"] = agp_parts_df["part_start/gap_type"] == 0
agp_parts_df["three_prime_is_end"] = agp_parts_df["part_end/linkage"] == agp_parts_df["full_sequence_length"]
agp_parts_df["full_len_part"] = agp_parts_df["five_prime_is_start"] & agp_parts_df["three_prime_is_end"]

agp_parts_df.to_csv(f"{output_prefix}.parts.tsv", sep="\t", index=True, header=True)

agp_processed_parts_df = agp_parts_df[~agp_parts_df["full_len_part"]]
agp_processed_parts_df.to_csv(f"{output_prefix}.processed_parts.tsv", sep="\t", index=True, header=True)
#----

#---- Assign uniq ids to the breakpoints ----

breakpoint_start_df = agp_processed_parts_df[~agp_processed_parts_df["five_prime_is_start"]][["part_id/gap_length", "part_start/gap_type"]]
breakpoint_start_df.columns = pd.Index(["scaffold", "position"])
breakpoint_start_df["BP_start"] = True
breakpoint_start_df["BP_end"] = False
breakpoint_end_df = agp_processed_parts_df[~agp_processed_parts_df["three_prime_is_end"]][["part_id/gap_length", "part_end/linkage"]]
breakpoint_end_df.columns = pd.Index(["scaffold", "position"])
breakpoint_end_df["BP_start"] = False
breakpoint_end_df["BP_end"] = True

breakpoint_df = pd.concat([breakpoint_start_df, breakpoint_end_df]).sort_values(by=["scaffold", "position"])
breakpoint_uniq_df = breakpoint_df[["scaffold", "position"]].drop_duplicates()
breakpoint_uniq_df["BP_id"] = list(map(lambda s: f"BP{s}",
                                           range(0, len(breakpoint_uniq_df))))
breakpoint_uniq_df.set_index(["scaffold", "position"], inplace=True)

breakpoint_df["BP_id"] = breakpoint_df[["scaffold", "position"]].apply(lambda s: breakpoint_uniq_df.loc[(s.iloc[0], s.iloc[1]), "BP_id"], axis=1)
agp_processed_parts_df.loc[:, "BP_start"] = breakpoint_df.loc[breakpoint_df["BP_start"], "BP_id"]
agp_processed_parts_df.loc[:, "BP_end"] = breakpoint_df.loc[breakpoint_df["BP_end"], "BP_id"]
agp_processed_parts_df.to_csv(f"{output_prefix}.processed_parts.extended.tsv", sep="\t", index=True, header=True)

breakpoint_uniq_df.reset_index(drop=False, inplace=True)
breakpoint_uniq_df.set_index("BP_id", inplace=True)
#----
#---- Find surrounding gaps ----

def find_surrounding_gaps(brk_df):
    scaffold_id = brk_df["scaffold"].iloc[0]
    bp_df = deepcopy(brk_df)
    bp_df[["nested_gap", "nested_start", "nested_end", "nested_distance",
                 "left_gap", "left_start", "left_end", "left_distance",
                 "right_gap", "right_start", "right_end", "right_distance"]] = pd.NA
    if scaffold_id in fasta_col.gaps_bed.records.index:
        gap_df = fasta_col.gaps_bed.records.loc[[scaffold_id]]

        for BP_id in bp_df.index:
            position = bp_df.loc[BP_id, "position"]
            gap_df["start_distance"] = position - gap_df["start"]
            gap_df["end_distance"] = position - gap_df["end"] + 1
            gap_df["nested"] = (gap_df["start_distance"] >= 0) & (gap_df["end_distance"] <= 0)
            gap_df["left"] = gap_df["end_distance"] > 0
            gap_df["right"] = gap_df["start_distance"] < 0

            if sum(gap_df["nested"]) > 1:
                raise ValueError(f"ERROR!!! More than one nested gap were detected for breakpoint {BP_id} ({position}) !")
            elif sum(gap_df["nested"]) == 1:
                tmp = gap_df[gap_df["nested"]]
                bp_df.loc[BP_id, "nested_gap"] = tmp["gap_id"].iloc[0]
                bp_df.loc[BP_id, "nested_start"] = tmp["start"].iloc[0]
                bp_df.loc[BP_id, "nested_end"] = tmp["end"].iloc[0]
                bp_df.loc[BP_id, "nested_distance"] = 0
            if not gap_df[gap_df["left"]].empty:
                tmp = gap_df[gap_df["left"]]
                min_dist = min(tmp["end_distance"])
                tmp = tmp[tmp["end_distance"] == min_dist]
                if len(tmp) > 1:
                    raise ValueError(f"ERROR!!! More than one closest left gap were detected for breakpoint {BP_id} ({position})! ")
                bp_df.loc[BP_id, "left_gap"] = tmp["gap_id"].iloc[0]
                bp_df.loc[BP_id, "left_start"] = tmp["start"].iloc[0]
                bp_df.loc[BP_id, "left_end"] = tmp["end"].iloc[0]
                bp_df.loc[BP_id, "left_distance"] = tmp["end_distance"].iloc[0]
            if not gap_df[gap_df["right"]].empty:
                tmp = gap_df[gap_df["right"]]
                min_dist = max(tmp["end_distance"])
                tmp = tmp[tmp["end_distance"] == min_dist]

                if len(tmp) > 1:
                    raise ValueError(f"ERROR!!! More than one closest right gap were detected for breakpoint {BP_id} ({position})! ")
                bp_df.loc[BP_id, "right_gap"] = tmp["gap_id"].iloc[0]
                bp_df.loc[BP_id, "right_start"] = tmp["start"].iloc[0]
                bp_df.loc[BP_id, "right_end"] = tmp["end"].iloc[0]
                bp_df.loc[BP_id, "right_distance"] = -tmp["start_distance"].iloc[0]

    else:
        print(f"WARNING!!! Scaffold {scaffold_id} has no gaps, but was processed!!!")
    for column in ["nested_start", "nested_end", "nested_distance",
                   "left_start", "left_end", "left_distance",
                   "right_start", "right_end", "right_distance"]:
        bp_df[column] = bp_df[column].astype("Int64")
    return bp_df

breakpoint_uniq_df = breakpoint_uniq_df.groupby(["scaffold"]).apply(find_surrounding_gaps).reset_index(level=0, drop=True) #
#----

#---- Select_closest gap ----
def select_gap(row_series):
    if not pd.isna(row_series.loc["nested_gap"]):
        series = row_series.loc[["nested_gap", "nested_start", "nested_end", "nested_distance"]]
    elif (not pd.isna(row_series.loc["left_gap"])) and (not pd.isna(row_series.loc["right_gap"])):
        if row_series.loc["left_distance"] <= row_series.loc["right_distance"]:
            series = row_series.loc[["left_gap", "left_start", "left_end", "left_distance"]]
        else:
            series = row_series.loc[["right_gap", "right_start", "right_end", "right_distance"]]
    elif not pd.isna(row_series.loc["left_gap"]):
        series = row_series.loc[["left_gap", "left_start", "left_end", "left_distance"]]
    elif not pd.isna(row_series.loc["right_gap"]):
        series = row_series.loc[["right_gap", "right_start", "right_end", "right_distance"]]
    else:
        series = pd.Series([pd.NA] * 4)

    series.index = pd.Index(["closest_gap", "closest_start", "closest_end", "closest_distance"])
    return pd.concat([row_series, series])

breakpoint_uniq_df = breakpoint_uniq_df.apply(select_gap, axis=1)
breakpoint_uniq_df["closest_distance/texel_size"] = breakpoint_uniq_df["closest_distance"] / texel_size
breakpoint_uniq_df["closest_distance_OK"] = breakpoint_uniq_df["closest_distance"] <= max_gap_dist

# Verify the selected gaps. Check following:
#
# 1. selected gap can't have coordinates smaller than previous breakpoints
# 2. selected gap can't have coordinates bigger than following breakpoints
# 3. same gap can't be selected for several breakpoints
# 4. selected gap can't precede any of gaps selected for previous breakpoints
# 5. selected gap must be precede all gaps selected for all following breakpoints
#
# If any of these requirements is not meet, prohibit the replacement for all breakpoints which are involved in the issue.

def verify_breakpoints(per_scaff_brk_uniq_df):
    tmp_brk_uniq_df = deepcopy(per_scaff_brk_uniq_df)
    tmp_brk_uniq_df_len = len(tmp_brk_uniq_df)
    # get index from gap id of the closest gap
    tmp_brk_uniq_df["internal_contig_breakpoint_index"] = list(range(0, tmp_brk_uniq_df_len))
    tmp_brk_uniq_df["closest_gap_index"] = \
                    tmp_brk_uniq_df["closest_gap"].apply(lambda s: pd.NA if pd.isnull(s) else int(s[3:])).astype("Int64")

    # replace index and coordinates for declined gaps, i.e gaps with "closest_distance" > max_gap_dist

    # initiate all criteria columns with false
    tmp_brk_uniq_df["closest_gap_coords_smaller_prev_breakpoints"] = False  # criterion 1
    tmp_brk_uniq_df["closest_gap_coords_bigger_follow_breakpoints"] = False # criterion 2
    tmp_brk_uniq_df["closest_gap_multiuse"] = False # criterion 3
    tmp_brk_uniq_df["closest_gap_preceeds_prev_gaps"] = False  # criterion 4
    tmp_brk_uniq_df["closest_gap_follows_follow_gaps"] = False  # criterion 5

    # case of single breakpoint contig
    if tmp_brk_uniq_df_len == 1:
        tmp_brk_uniq_df["use_closest"] = tmp_brk_uniq_df["closest_distance_OK"]
        return tmp_brk_uniq_df

    for internal_brk_index in tmp_brk_uniq_df["internal_contig_breakpoint_index"][tmp_brk_uniq_df["closest_distance_OK"]].iloc[1:]:
        # check criterion 1
        distance_ok_sr = tmp_brk_uniq_df["closest_distance_OK"].iloc[0:internal_brk_index]
        prev_brk_position_sr = tmp_brk_uniq_df["position"].iloc[0:internal_brk_index][distance_ok_sr]
        # count fails of criterion 1
        failed_cr1_sum = sum(prev_brk_position_sr >= tmp_brk_uniq_df["closest_start"].iloc[internal_brk_index])

        tmp_brk_uniq_df.loc[tmp_brk_uniq_df.index[internal_brk_index],
                            "closest_gap_coords_smaller_prev_breakpoints"] = True if failed_cr1_sum > 0 else False
        # check criterion 4
        prev_gap_indexes_sr = tmp_brk_uniq_df["closest_gap_index"].iloc[0:internal_brk_index][distance_ok_sr]
        failed_cr4_sum = sum(prev_gap_indexes_sr >= tmp_brk_uniq_df["closest_gap_index"].iloc[internal_brk_index])
        tmp_brk_uniq_df.loc[tmp_brk_uniq_df.index[internal_brk_index],
                            "closest_gap_preceeds_prev_gaps"] = True if failed_cr4_sum > 0 else False

    for internal_brk_index in tmp_brk_uniq_df["internal_contig_breakpoint_index"][tmp_brk_uniq_df["closest_distance_OK"]].iloc[:-1]:
        # check criterion 2
        distance_ok_sr = tmp_brk_uniq_df["closest_distance_OK"].iloc[internal_brk_index + 1:]
        follow_brk_position_sr = tmp_brk_uniq_df["position"].iloc[internal_brk_index + 1:][distance_ok_sr]
        # count fails of criterion 2
        failed_cr2_sum = sum(follow_brk_position_sr < tmp_brk_uniq_df["closest_end"].iloc[internal_brk_index])

        tmp_brk_uniq_df.loc[tmp_brk_uniq_df.index[internal_brk_index],
                            "closest_gap_coords_smaller_prev_breakpoints"] = True if failed_cr2_sum > 0 else False

        # check criterion 5
        follow_gap_indexes_sr = tmp_brk_uniq_df["closest_gap_index"].iloc[internal_brk_index + 1:][distance_ok_sr]
        failed_cr5_sum = sum(follow_gap_indexes_sr <= tmp_brk_uniq_df["closest_gap_index"].iloc[internal_brk_index])
        tmp_brk_uniq_df.loc[tmp_brk_uniq_df.index[internal_brk_index],
                            "closest_gap_follows_follow_gaps"] = True if failed_cr5_sum > 0 else False


    for internal_brk_index in tmp_brk_uniq_df["internal_contig_breakpoint_index"][tmp_brk_uniq_df["closest_distance_OK"]]:
        # check criterion 3
        distance_ok_closest_gap_sr = tmp_brk_uniq_df["closest_gap"][tmp_brk_uniq_df["closest_distance_OK"]]
        same_gap_id_counts = sum(distance_ok_closest_gap_sr == tmp_brk_uniq_df["closest_gap"].iloc[internal_brk_index])
        tmp_brk_uniq_df.loc[tmp_brk_uniq_df.index[internal_brk_index],
                            "closest_gap_multiuse"] = True if same_gap_id_counts > 1 else False


    tmp_brk_uniq_df["use_closest"] = ~(tmp_brk_uniq_df["closest_gap_coords_smaller_prev_breakpoints"] |
                                       tmp_brk_uniq_df["closest_gap_coords_bigger_follow_breakpoints"] |
                                       tmp_brk_uniq_df["closest_gap_multiuse"] |
                                       tmp_brk_uniq_df["closest_gap_preceeds_prev_gaps"] |
                                       tmp_brk_uniq_df["closest_gap_follows_follow_gaps"] ) & tmp_brk_uniq_df["closest_distance_OK"]

    return tmp_brk_uniq_df

breakpoint_uniq_df = breakpoint_uniq_df.groupby("scaffold").apply(verify_breakpoints).reset_index(level=0, drop=True)
breakpoint_uniq_df.to_csv(f"{output_prefix}.breakpoint.uniq.tsv", sep="\t", index=True, header=True)
#----
#---- Add corrected_coordinates to agp_processed_parts_df and write it ----
for column in ["corrected_start", "corrected_end"]:
    agp_processed_parts_df.loc[:, column] = pd.Series(pd.NA, index=agp_processed_parts_df.index, dtype='Int64')

# copy
agp_processed_parts_df.loc[agp_processed_parts_df["five_prime_is_start"], "corrected_start"] = agp_processed_parts_df.loc[agp_processed_parts_df["five_prime_is_start"], "part_start/gap_type"]
agp_processed_parts_df.loc[agp_processed_parts_df["three_prime_is_end"], "corrected_end"] = agp_processed_parts_df.loc[agp_processed_parts_df["three_prime_is_end"], "part_end/linkage"]

def get_breakpoint_start_coordinates(breakpoint_id):
    if breakpoint_uniq_df.loc[breakpoint_id, "use_closest"]:
        return breakpoint_uniq_df.loc[breakpoint_id, "closest_end"]
    else:
        # return original coordinate without correction
        return breakpoint_uniq_df.loc[breakpoint_id, "position"]

def get_breakpoint_end_coordinates(breakpoint_id):
    if breakpoint_uniq_df.loc[breakpoint_id, "use_closest"]:
        return breakpoint_uniq_df.loc[breakpoint_id, "closest_start"]

    else:
        # return original coordinate without correction
        return breakpoint_uniq_df.loc[breakpoint_id, "position"]

agp_processed_parts_df.loc[~agp_processed_parts_df["five_prime_is_start"], "corrected_start"] = \
        agp_processed_parts_df.loc[~agp_processed_parts_df["five_prime_is_start"],
                                   "BP_start"].apply(get_breakpoint_start_coordinates)

agp_processed_parts_df.loc[~agp_processed_parts_df["three_prime_is_end"], "corrected_end"] = \
        agp_processed_parts_df.loc[~agp_processed_parts_df["three_prime_is_end"],
                                   "BP_end"].apply(get_breakpoint_end_coordinates)

for column in ["corrected_start", "corrected_end"]:
    agp_processed_parts_df.loc[:, column] = agp_processed_parts_df[column].astype(int)

agp_processed_parts_df.to_csv(f"{output_prefix}.corrected_processed_segments.tsv", sep="\t", header=True, index=True)
#----

#---- Correct original agp ----
tmp = deepcopy(agp_processed_parts_df[["corrected_start", "corrected_end"]])
tmp.columns = pd.Index(["part_start/gap_type", "part_end/linkage"])
agp_col.records.loc[agp_processed_parts_df.index, ["part_start/gap_type",
                                                   "part_end/linkage"]] = tmp

tmp_index = 1

def correct_scaffold_coordinates_from_segments(scaf_df):
    scaf_df_len = len(scaf_df)
    tmp = deepcopy(scaf_df)
    tmp_columns = list(tmp.columns)
    tmp["seg_len"] = 0
    U_type_sr = tmp["part_type"] == "U"
    W_type_sr = tmp["part_type"] == "W"

    tmp.loc[U_type_sr, "seg_len"] = tmp.loc[U_type_sr, "part_id/gap_length"].astype(int)
    tmp.loc[W_type_sr, "seg_len"] = \
        tmp.loc[W_type_sr, "part_end/linkage"].astype(int) - tmp.loc[W_type_sr, "part_start/gap_type"].astype(int)

    tmp.loc[tmp.index[0], "start"] = 0
    tmp.loc[tmp.index[0], "end"] = tmp.loc[tmp.index[0], "seg_len"]

    if scaf_df_len > 1:
        for i in range(1, scaf_df_len):
            tmp.loc[tmp.index[i], "start"] = tmp.loc[tmp.index[i-1], "end"]
            tmp.loc[tmp.index[i], "end"] = tmp.loc[tmp.index[i], "start"] + tmp.loc[tmp.index[i], "seg_len"]

    return tmp[tmp_columns]

agp_col.records = agp_col.records.groupby(["scaffold"]).apply(correct_scaffold_coordinates_from_segments).reset_index(level=0, drop=True)
agp_col.write(f"{output_prefix}.agp")
#----

#---- Run tests ----
def check_agp_scaffolds_df(df):
    results_df = pd.DataFrame([[False] * 5], columns=pd.Index(["scaffold_part_number_test",
                                                               "scaffold_start_test",
                                                               "scaffold_fragment_length_test",
                                                               "scaffold_segment_end_test",
                                                               "scaffold_missing_bases_test"]))
    # check part number

    results_df.loc[results_df.index[0], "scaffold_part_number_test"] = True if df.loc[df.index[-1], "part_number"] == len(df) else False
    #print(part_number_test_flag)

    # check for scaffold start
    results_df.loc[results_df.index[0], "scaffold_start_test"] = True if df.loc[df.index[0], "start"] == 0 else False

    # check for length of segments
    results_df.loc[results_df.index[0], "scaffold_fragment_length_test"] = \
         True if sum((df.loc[agp_col.records["part_type"] != "U",
                             "end"] - df.loc[agp_col.records["part_type"] != "U",
                                             "start"]) == (df.loc[agp_col.records["part_type"] != "U",
                                                                  "part_end/linkage"] - df.loc[agp_col.records["part_type"] != "U",
                                                                                               "part_start/gap_type"])) ==  len(df.loc[agp_col.records["part_type"] != "U"]) else False

    # check scaffold segment ends (should be bigger then start)
    results_df.loc[results_df.index[0], "scaffold_segment_end_test"] = True if sum(df["end"] - df["start"] > 0) == len(df) else False

    # check for missing bases in scaffold
    if len(df) > 1:
        results_df.loc[results_df.index[0], "scaffold_missing_bases_test"] = True if sum((df["start"] - df["end"].shift(periods=1)).iloc[1:] == 0)  == len(df) - 1 else False
    else:
        results_df["scaffold_missing_bases_test"] = True

    return results_df

def check_agp_fragments_df(fragment_df):
    results_df = pd.DataFrame([[False, False, 0, 0, 0, 0, 0, 0]], columns=pd.Index(["overlap_test",
                                                                              "gap_test",
                                                                              "breakpoint_number",
                                                                              "gapless_breakpoint_number",
                                                                              "gap_number",
                                                                              "total_gap_length",
                                                                              "mean_gap_length",
                                                                              "median_gap_length"]))
    gap_length_series = (fragment_df["part_start/gap_type"] - fragment_df["part_end/linkage"].shift(periods=1)).iloc[1:]
    if len(fragment_df) > 1:
        results_df["overlap_test"] = True if sum(gap_length_series >= 0) == (len(fragment_df) - 1)  else False
        results_df["gap_presence_test"] = True if sum(gap_length_series > 0) == 0 else False
        results_df["breakpoint_number"] = len(gap_length_series[gap_length_series >= 0])
        results_df["gap_number"] = len(gap_length_series[gap_length_series > 0])
        results_df["gapless_breakpoint_number"] = results_df["breakpoint_number"] - results_df["gap_number"]
        results_df["total_gap_length"] = sum(gap_length_series)
        results_df["mean_gap_length"] = gap_length_series[gap_length_series > 0].mean()
        results_df["median_gap_length"] = gap_length_series[gap_length_series > 0].median()
    else:
        results_df["overlap_test"] = True
        results_df["gap_test"] = True
        results_df["breakpoint_number"] = 0
        results_df["gap_number"] = 0
        results_df["gapless_breakpoint_number"] = 0
        results_df["total_gap_length"] = 0
        results_df["mean_gap_length"] = pd.NA
        results_df["median_gap_length"] = pd.NA

    return results_df

def verify_agp(agp_records):
    fragment_test_df = agp_col.records[agp_col.records["part_type"] != "U"][["part_id/gap_length",
                                                          "part_start/gap_type",
                                                          "part_end/linkage"]].set_index("part_id/gap_length").sort_values(by=["part_id/gap_length",
                                                                                                                               "part_start/gap_type",
                                                                                                                               "part_end/linkage"])
    fragment_test_df = fragment_test_df.groupby("part_id/gap_length").apply(check_agp_fragments_df).reset_index(level=1, drop=True).sort_values("part_id/gap_length")
    scaffold_test_df = agp_col.records.groupby("scaffold").apply(check_agp_scaffolds_df).reset_index(level=1, drop=True)
    scaffold_test_df["all_tests_passed"] = scaffold_test_df["scaffold_part_number_test"] & \
                                           scaffold_test_df["scaffold_start_test"] & \
                                           scaffold_test_df["scaffold_fragment_length_test"] & \
                                           scaffold_test_df["scaffold_segment_end_test"] & \
                                           scaffold_test_df["scaffold_missing_bases_test"]
    return scaffold_test_df, fragment_test_df

test_results = verify_agp(agp_col.records)
test_results[0].to_csv(f"{output_prefix}.checks.scaffold.tsv", sep="\t", index=True, header=True)
test_results[1].to_csv(f"{output_prefix}.checks.fragments.tsv", sep="\t", index=True, header=True)
#----
