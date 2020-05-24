#!/usr/bin/env python

import collections
import pandas as pd

import tools.intervals
import tools.misc
import tools.mathOps
import tools.fileOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
import tools.procOps
from cat.consensus import *
from argparse import ArgumentParser

def get_df(db_path, ref_db_path, chrXOnly=False):
    tm_eval_df = load_transmap_evals(db_path)
    ref_df = tools.sqlInterface.load_annotation(ref_db_path)
    tx_modes = ['transMap']
    mrna_metrics_df = pd.concat([load_metrics_from_db(db_path, tx_mode, 'mRNA') for tx_mode in tx_modes])
    cds_metrics_df = pd.concat([load_metrics_from_db(db_path, tx_mode, 'CDS') for tx_mode in tx_modes])
    eval_df = pd.concat([load_evaluations_from_db(db_path, tx_mode) for tx_mode in tx_modes]).reset_index()
    hgm_df = pd.concat([load_hgm_vectors(db_path, tx_mode) for tx_mode in tx_modes])


    # if chrXOnly:
    #     cmd = [['grep', 'chrX', 'gencode.v30.annotation.gp'], ['cut', '-f', '1']]
    #     chrx_txs = set(tools.procOps.call_proc_lines(cmd))
    #     ref_df = ref_df[ref_df.TranscriptId.isin(chrx_txs)]
    # else:
    #     # remove chrY
    #     cmd = [['grep', 'chrY', 'gencode.v30.annotation.gp'], ['cut', '-f', '1']]
    #     chry_txs = set(tools.procOps.call_proc_lines(cmd))
    #     ref_df = ref_df[~ref_df.TranscriptId.isin(chry_txs)]


    num_txs = len(set(ref_df[ref_df.TranscriptBiotype == 'protein_coding'].TranscriptId))
    num_genes = len(set(ref_df[ref_df.TranscriptBiotype == 'protein_coding'].GeneId))


    ## code below is from the consensus module. I ripped out from the combine_and_filter_dfs method
    ## because it needs the genePred, but the info is also present elsewhere.

    #add the reference information to gain biotype information
    hgm_ref_df = pd.merge(hgm_df, ref_df, on=['GeneId', 'TranscriptId'])
    # combine in homGeneMapping results
    hgm_ref_tm_df = pd.merge(hgm_ref_df, tm_eval_df, on=['GeneId', 'TranscriptId'])
    # split merged_df into coding and noncoding
    coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype == 'protein_coding']
    non_coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype != 'protein_coding']
    # add metrics information to coding df
    metrics_df = pd.merge(mrna_metrics_df, cds_metrics_df, on='AlignmentId', suffixes=['_mRNA', '_CDS'])
    coding_df = pd.merge(coding_df, metrics_df, on='AlignmentId')
    # add evaluation information to coding df, where possible. This adds information on frame shifts.
    coding_df = pd.merge(coding_df, eval_df, on='AlignmentId', how='left')
    # fill the original intron values to 100 so we don't filter them out -- means a no-intron gene
    coding_df['OriginalIntronsPercent_mRNA'] = coding_df.OriginalIntronsPercent_mRNA.fillna(100)
    coding_df['OriginalIntronsPercent_CDS'] = coding_df.OriginalIntronsPercent_CDS.fillna(100)
    non_coding_df['TransMapOriginalIntronsPercent'] = non_coding_df.TransMapOriginalIntronsPercent.fillna(100)


    # rawest counts. homGeneMapping was ran on the unfiltered dataset, so use that
    # do this only on coding transcripts for now
    unique_genes = 0
    unique_txs = 0
    tmp = hgm_ref_df[hgm_ref_df.TranscriptBiotype == 'protein_coding']
    num_coding_genes = len(set(tmp.GeneId))
    num_coding_txs = len(set(tmp.TranscriptId))
    for gene_id, d in tmp.groupby('GeneId'):
        paralogy = collections.Counter(x.split('-')[0] for x in d.AlignmentId)
        if sum(paralogy.values()) == len(paralogy):
            unique_genes += 1
        for tx_id, dd in d.groupby('TranscriptId'):
            if len(dd) == 1:
                unique_txs += 1

    data = {}
    data['GenesFound'] = num_coding_genes
    data['GenesFoundPercent'] = 100.0 * num_coding_genes / num_genes
    data['GenesMultiplyMapping'] = num_genes - unique_genes
    data['GenesMultiplyMappingPercent'] = 100.0 * (num_genes - unique_genes) / num_genes
    data['TranscriptsFound'] = num_coding_txs
    data['TranscriptsFoundPercent'] = 100.0 * num_coding_txs / num_txs
    data['TranscriptsMultiplyMapping'] = num_txs - unique_txs
    data['TranscriptsMultiplyMappingPercent'] = 100.0 * (num_txs - unique_txs) / num_txs

    # full coverage
    full_cov_mrna = len(coding_df[coding_df.AlnCoverage_mRNA == 100])
    full_cov_cds = len(coding_df[coding_df.AlnCoverage_CDS == 100])
    data['FullmRNACoverage'] = full_cov_mrna
    data['FullmRNACoveragePercent'] = 100.0 * full_cov_mrna / num_txs
    data['FullCDSCoverage'] = full_cov_cds
    data['FullCDSCoveragePercent'] = 100.0 * full_cov_cds / num_txs

    # construct a stringent filter that requires the following:
    # 1) Has all original introns
    # 2) Full CDS Coverage
    # 3) No Frame-shift
    frameshift = len(coding_df[coding_df.Frameshift == True])
    original_introns = len(coding_df[coding_df.OriginalIntronsPercent_mRNA == 100])
    cov = len(coding_df[coding_df.AlnCoverage_CDS == 100])
    cov_frameshift = len(coding_df[(coding_df.AlnCoverage_CDS == 100) &
                                              (coding_df.Frameshift != True)])
    cov_frameshift_original_introns = len(coding_df[(coding_df.AlnCoverage_CDS == 100) &
                                                    (coding_df.Frameshift != True) &
                                                    (coding_df.OriginalIntronsPercent_mRNA == 100)])
    data['TranscriptsWithFrameshift'] = frameshift
    data['TranscriptsWithFrameshiftPercent'] = 100.0 * frameshift / num_txs
    data['TranscriptsWithOriginalIntrons'] = original_introns
    data['TranscriptsWithOriginalIntronsPercent'] = 100.0 * original_introns / num_txs
    data['TranscriptsWithFullCDSCoverage'] = cov
    data['TranscriptsWithFullCDSCoveragePercent'] = 100.0 * cov / num_txs
    data['TranscriptsWithFullCDSCoverageAndNoFrameshifts'] = cov_frameshift
    data['TranscriptsWithFullCDSCoverageAndNoFrameshiftsPercent'] = 100.0 * cov_frameshift / num_txs
    data['TranscriptsWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntrons'] = cov_frameshift_original_introns
    data['TranscriptsWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntronsPercent'] = 100.0 * cov_frameshift_original_introns / num_txs

    # naive gene level
    frameshift = len(set(coding_df[coding_df.Frameshift == True].GeneId))
    original_introns = len(set(coding_df[coding_df.OriginalIntronsPercent_mRNA == 100].GeneId))
    cov = len(set(coding_df[(coding_df.ProperOrf == True) & (coding_df.AlnCoverage_CDS == 100)].GeneId))
    cov_frameshift = len(set(coding_df[(coding_df.AlnCoverage_CDS == 100) &
                                       (coding_df.Frameshift != True)].GeneId))
    cov_frameshift_original_introns = len(set(coding_df[(coding_df.AlnCoverage_CDS == 100) &
                                                        (coding_df.Frameshift != True) &
                                                        (coding_df.OriginalIntronsPercent_mRNA == 100)].GeneId))
    data['GenesWithFrameshift'] = frameshift
    data['GenesWithFrameshiftPercent'] = 100.0 * frameshift / num_genes
    num_genes_all_shifted = 0
    for gene_id, d in coding_df.groupby('GeneId'):
        if len(d[d.Frameshift == True]) == len(d):
            num_genes_all_shifted += 1
    data['GenesWithFrameshiftAllIsoforms'] = num_genes_all_shifted
    data['GenesWithFrameshiftAllIsoformsPercent'] = 100.0 * num_genes_all_shifted / num_genes
    data['GenesWithOriginalIntrons'] = original_introns
    data['GenesWithOriginalIntronsPercent'] = 100.0 * original_introns / num_genes
    data['GenesWithFullCDSCoverage'] = cov
    data['GenesWithFullCDSCoveragePercent'] = 100.0 * cov / num_genes
    data['GenesWithFullCDSCoverageAndNoFrameshifts'] = cov_frameshift
    data['GenesWithFullCDSCoverageAndNoFrameshiftsPercent'] = 100.0 * cov_frameshift / num_genes
    data['GenesWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntrons'] = cov_frameshift_original_introns
    data['GenesWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntronsPercent'] = 100.0 * cov_frameshift_original_introns / num_genes

    missing = set(ref_df[ref_df.TranscriptBiotype == 'protein_coding'].GeneId) - set(tmp.GeneId)

    data['MissingGenes'] = len(missing)
    data['MissingGenesPercent'] = (100.0 * len(missing)) / num_genes

    data['Name'] = db_path.replace('.db', '')

    return data

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('ref_db')
    parser.add_argument('target_db', nargs='+')
    parser.add_argument('--chrXOnly', action='store_true')
    opts = parser.parse_args()
    columns = ['Name', 'GenesFound', 'TranscriptsFound', 'GenesMultiplyMapping', 'TranscriptsMultiplyMapping', 'FullmRNACoverage', 'FullCDSCoverage', 'TranscriptsWithFrameshift', 'TranscriptsWithOriginalIntrons', 'TranscriptsWithFullCDSCoverage', 'TranscriptsWithFullCDSCoverageAndNoFrameshifts', 'TranscriptsWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntrons', 'GenesWithFrameshift', 'GenesWithOriginalIntrons', 'GenesWithFullCDSCoverage', 'GenesWithFullCDSCoverageAndNoFrameshifts', 'GenesWithFullCDSCoverageAndNoFrameshiftsAndOriginalIntrons', 'MissingGenes']
    columns_w_percent = []
    for column in columns:
        if column != 'Name':
            columns_w_percent.extend([column, column + 'Percent'])
        else:
            columns_w_percent.append(column)
    df = pd.DataFrame(columns=columns_w_percent)
    for target_db in opts.target_db:
        df = df.append(get_df(target_db, opts.ref_db, opts.chrXOnly), ignore_index=True)
    print(df.to_csv())