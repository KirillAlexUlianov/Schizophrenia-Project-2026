#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm
import bioframe as bf
import os
from typing import Set, Dict, Any, Union
from pathlib import Path
import argparse


def bedpe_to_bed(df, ann_col='loop_num', slop=False):
    """
    Extract unique bed regions from bedpe dataframe. Index is ignored.
    """
    df_reg1 = df[['chrom1', 'start1', 'end1', ann_col]].copy()
    df_reg1['anchor'] = 'left_anh'
    df_reg2 = df[['chrom2', 'start2', 'end2', ann_col]].copy()
    df_reg2['anchor'] = 'right_anh'
    df_reg1.columns = df_reg1.columns.str.rstrip('1')
    df_reg2.columns = df_reg2.columns.str.rstrip('2')

    df_bed = pd.concat([df_reg1, df_reg2], ignore_index=True)\
               .drop_duplicates(ignore_index=True)    
    if slop:
        df_bed.start = df_bed.start - slop
        df_bed.start = df_bed.start.astype(int)
        df_bed.end = df_bed.end + slop
        df_bed.end = df_bed.end.astype(int)

    return df_bed 

def bedpe_to_bed_multiple_cols(df, ann_cols=['num'], slope=0):
    """
    Extract unique bed regions from bedpe dataframe. Index is ignored.
    """
    df_reg1 = df[['chrom1', 'start1', 'end1'] + ann_cols].copy()
    df_reg1['anchor'] = 'left_anh'
    df_reg2 = df[['chrom2', 'start2', 'end2'] + ann_cols].copy()
    df_reg2['anchor'] = 'right_anh'
    df_reg1.columns = df_reg1.columns.str.rstrip('1')
    df_reg2.columns = df_reg2.columns.str.rstrip('2')

    df_bed = pd.concat([df_reg1, df_reg2], ignore_index=True)\
               .drop_duplicates(ignore_index=True)
    df_bed.start = df_bed.start - slope
    df_bed.end = df_bed.end + slope
    return df_bed 

def intersection_with_h3k27me3(merged_all_clusters, h3k27me3_post_mortem):
    merged_all_clusters_bed = bedpe_to_bed_multiple_cols(merged_all_clusters, ann_cols=['cluster'], slope=0)
    
    merged_all_clusters_stats_h3k27me3_kz = get_overlap_with_h3k27me3(merged_all_clusters_bed, h3k27me3_post_mortem, 'has_h3k27me3_peak_kozlenkov', factor_slop = 2)
    
    assert merged_all_clusters_stats_h3k27me3_kz.cluster.nunique() == merged_all_clusters.cluster.nunique()

    init_cluster_number = merged_all_clusters.cluster.nunique()
    assert  'has_h3k27me3_peak_kozlenkov' in merged_all_clusters_stats_h3k27me3_kz.columns.tolist()

    merged_all_clusters = merged_all_clusters.merge(merged_all_clusters_stats_h3k27me3_kz[['cluster', 'has_h3k27me3_peak_kozlenkov_anchor', 'has_h3k27me3_peak_kozlenkov']].drop_duplicates(), on=['cluster'], how='left')
    assert  'has_h3k27me3_peak_kozlenkov' in merged_all_clusters.columns.tolist()

    assert  merged_all_clusters.cluster.nunique() == init_cluster_number

    merged_all_clusters['h3k27me3_peak_status_kozlenkov'] = merged_all_clusters.has_h3k27me3_peak_kozlenkov.apply(lambda x: "has h3k27me3 peak" if x else 'do not have h3k27me3 peak')
    
    return merged_all_clusters

def get_peak_status(df, column_new_name):
    anchor_status_all = {}
    anchor_status_general = {}
    grouped = df.groupby(['cluster', 'anchor'])[column_new_name].sum().unstack(fill_value=0)
    print('Getting peak status...')

    for cluster, row in tqdm(grouped.iterrows()):
        left_anh = row.get('left_anh', 0)
        right_anh = row.get('right_anh', 0)

        if right_anh > 0 and left_anh > 0:
            anchor_status_all[cluster] = 'both'
            anchor_status_general[cluster] = True
        elif right_anh > 0 or left_anh > 0:
            anchor_status_all[cluster] = 'one_anh'
            anchor_status_general[cluster] = True
        else:
            anchor_status_all[cluster] = 'none'
            anchor_status_general[cluster] = False

    return anchor_status_all, anchor_status_general

def intersection_with_gene_annotation(merged_all_clusters, genes, factor_slop=2):
    clusters_allStats_bed = bedpe_to_bed(merged_all_clusters, ann_col='cluster', slop=20000*factor_slop)
    genes_in_clusters = bf.overlap(genes, clusters_allStats_bed, how='left', suffixes=('','_anno'))
    genes_in_clusters = genes_in_clusters[['gene_name', 'gene_type', 'cluster_anno']].rename(columns={'cluster_anno': 'cluster'})
    genes_in_clusters.drop_duplicates(inplace=True)
    merged_all_clusters = merged_all_clusters.merge(genes_in_clusters, how='left', on='cluster')
    merged_all_clusters.rename(columns={'gene_name': 'gene_name_anno', 'gene_type': 'gene_type_anno'}, inplace=True)
    return merged_all_clusters

def get_overlap_with_h3k27me3(df1, h3k27me3_encode, column_new_name, factor_slop=2):
    df1 = df1.drop_duplicates()
    df1.start = df1.start - 20000 * factor_slop
    df1.end = df1.end + 20000 * factor_slop
    df1.start = df1.start.astype(int)
    df1.end = df1.end.astype(int)

    df1_h3k27me3_encode = bf.overlap(df1, h3k27me3_encode, how='left', suffixes=('','_anno'))
    df1_h3k27me3_encode[f'{column_new_name}_temp'] = df1_h3k27me3_encode['start_anno'].apply(lambda x: True if isinstance(x, int) else False)
    anchor_status_all, anchor_status_general = get_peak_status(df1_h3k27me3_encode, f'{column_new_name}_temp')
    df1_h3k27me3_encode[f'{column_new_name}_anchor'] = df1_h3k27me3_encode['cluster'].map(anchor_status_all)
    df1_h3k27me3_encode[column_new_name] = df1_h3k27me3_encode['cluster'].map(anchor_status_general)
    df1_h3k27me3_encode.drop(columns=f'{column_new_name}_temp', inplace=True)

    df1_h3k27me3_encode.start = df1_h3k27me3_encode.start + 20000 * factor_slop
    df1_h3k27me3_encode.end = df1_h3k27me3_encode.end - 20000 * factor_slop
    df1_h3k27me3_encode.start = df1_h3k27me3_encode.start.astype(int)
    df1_h3k27me3_encode.end = df1_h3k27me3_encode.end.astype(int)

    return df1_h3k27me3_encode

def intersection_with_human_tfs(merged_all_clusters, tfs):
    cluster_to_tfs = {}
    df_grouped = merged_all_clusters.groupby("cluster")
    for cluster, group in tqdm(df_grouped):
        if tfs.intersection(group.gene_name_anno):
            cluster_to_tfs[cluster] = True
        else:
            cluster_to_tfs[cluster] = False
    merged_all_clusters['tf_in_cluster'] = merged_all_clusters.cluster.map(cluster_to_tfs)
    merged_all_clusters['tf_in_cluster_status'] = merged_all_clusters.tf_in_cluster.apply(lambda x: "contains TF" if x else "no TFs")
    return merged_all_clusters


def intersection_with_promoters(merged_all_clusters_stats, promoters, slop_factor=2):
    clusters_allStats_bed = bedpe_to_bed(merged_all_clusters_stats, ann_col='cluster', slop=20000*slop_factor)
    promoters_in_clusters = bf.overlap(promoters, clusters_allStats_bed, how='left', suffixes=('', '_promoters'))
    promoters_in_clusters = promoters_in_clusters[['gene_name', 'gene_type', 'cluster_promoters']].rename(columns={'gene_name': 'gene_name_promoter', 'gene_type': 'gene_type_promoter', 'cluster_promoters': 'cluster'})
    promoters_in_clusters.drop_duplicates(inplace=True)
    merged_all_clusters_stats = merged_all_clusters_stats.merge(promoters_in_clusters, how='left', on='cluster')
    return merged_all_clusters_stats

def intersection_with_polycomb_dots(merged_all_clusters_stats_trueLoops, pair_sites_bed, slop_factor=2):
    merged_all_clusters_stats_trueLoops_bed = bedpe_to_bed(merged_all_clusters_stats_trueLoops, ann_col='cluster', slop=20000*slop_factor)
    polycombDots_in_clusters = bf.overlap(merged_all_clusters_stats_trueLoops_bed, pair_sites_bed, how='left', suffixes=('', '_polycombDots'))

    groupped_pol = polycombDots_in_clusters[["cluster", 'chrom_polycombDots']].groupby("cluster")
    updated = {}
    for row in groupped_pol:
        all_status = row[1].chrom_polycombDots.tolist()
        if any(e != None for e in all_status):
            updated[row[0]] = True
        else:
            updated[row[0]] = False
    polycombDots_in_clusters_status, _ = get_peak_polycomb_dots_status(polycombDots_in_clusters)

    merged_all_clusters_stats_trueLoops['has_polycombDot_anchor'] = merged_all_clusters_stats_trueLoops.cluster.map(updated)
    merged_all_clusters_stats_trueLoops['has_polycombDot_anchor_status_detailed'] = merged_all_clusters_stats_trueLoops.cluster.map(polycombDots_in_clusters_status)

    return merged_all_clusters_stats_trueLoops, polycombDots_in_clusters


def get_peak_polycomb_dots_status(polycombDots_in_clusters):
    anchor_status_all = {}
    anchor_status_general = {}
    grouped = polycombDots_in_clusters.groupby(['cluster', 'anchor'])["chrom_polycombDots"].sum().unstack(fill_value=0)
    print('Getting peak status...')

    for cluster, row in tqdm(grouped.iterrows()):
        left_anh = row.get('left_anh', 0)
        right_anh = row.get('right_anh', 0)
        if isinstance(right_anh, str) and isinstance(left_anh, str):
            anchor_status_all[cluster] = 'both'
            anchor_status_general[cluster] = True
        elif isinstance(right_anh, str) or isinstance(left_anh, str):
            anchor_status_all[cluster] = 'one_anh'
            anchor_status_general[cluster] = True
        else:
            anchor_status_all[cluster] = 'none'
            anchor_status_general[cluster] = False

    assert len(anchor_status_all) == polycombDots_in_clusters['cluster'].nunique()
    assert len(anchor_status_general) == polycombDots_in_clusters['cluster'].nunique()
    print('Done')
    return anchor_status_all, anchor_status_general


def main():
    parser = argparse.ArgumentParser(description='Perform genomic intersection analysis')
    parser.add_argument('--loops', required=True, help='Path to loops file (BEDPE format)')
    
    base_path = Path('./0.additional_data')
    
    # Load input data
    print("Loading input data...")
    args = parser.parse_args()

    print("Reading loops file...")
    merged_all_clusters = pd.read_pickle(args.loops)
    
    # Load chromatin annotation
    print("Loading chromatin annotation...")
    chromatin_annotation = pd.read_pickle(base_path / 'chromatin_annotation.feather')
    promoters = chromatin_annotation[chromatin_annotation.type == "promoters"].reset_index(drop=True)
    genes = chromatin_annotation[chromatin_annotation.type == "gene"].reset_index(drop=True)
    
    # Load H3K27me3 data
    print("Loading H3K27me3 data...")
    h3k27me3_post_mortem = pd.read_csv(base_path / 'kozlenkov2018_merged_h3k27me3.bed',
                                      sep="\t", header=None,
                                      names=['chrom', 'start', 'end'])
    
    # Load TFs
    print("Loading transcription factors...")
    tfs = pd.read_csv(base_path / 'DatabaseExtract_v_1.01.csv')
    tfs = set(tfs[tfs['Is TF?'] == "Yes"]['HGNC symbol'].tolist())
    
    # Polycomb
    pair_sites = pd.read_csv(base_path / "pair_sites_polycomb_dot_anchors.csv")
    pair_sites['num'] = [i for i in range(pair_sites.shape[0])]
    pair_sites_bed = bedpe_to_bed(pair_sites, ann_col='num', slop=0)

    print("\nPerforming intersections...")
    merged_all_clusters = intersection_with_gene_annotation(merged_all_clusters, genes)
    merged_all_clusters = intersection_with_h3k27me3(merged_all_clusters, h3k27me3_post_mortem)
    merged_all_clusters = intersection_with_human_tfs(merged_all_clusters, tfs)
    merged_all_clusters, polycombDots_in_clusters = intersection_with_polycomb_dots(merged_all_clusters, pair_sites_bed)
    merged_all_clusters = intersection_with_promoters(merged_all_clusters, promoters)

    # Save results
    print("\nSaving results...")
    name = args.loops.split('/')[-1].split('.')[0]

    output_file = f'./loops_cooltools_data_noSexChromosomes/intersection_results/{name}_intersection.pickle'
    merged_all_clusters.to_pickle(output_file)
    print(f"Results saved to {output_file}")

if __name__ == '__main__':
    main()
