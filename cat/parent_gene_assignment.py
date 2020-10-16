"""
A set of functions to perform parental gene assignment in the AugustusPB/AugustusCGP modules
"""
import pandas as pd
import itertools
import collections
# from tools.defaultOrderedDict import DefaultOrderedDict
import tools.procOps
import tools.fileOps
import tools.mathOps
import tools.transcripts
import tools.intervals
import tools.nameConversions


def assign_parents(filtered_tm_gp, unfiltered_tm_gp, denovo_gp, min_distance=0.4,
                   tm_jaccard_distance=0.25, stranded=True):
    """
    Main function for assigning parental genes. Parental gene assignment methodology:
    A) clusterGenes is used to cluster filtered transMap transcripts.
    B) If a denovo transcript is assigned to more than one gene, then this is attempted to be resolved.
    Resolution occurs by looking first at the transMap themselves. 
    If any transMap projections overlap each other with a Jaccard metric > tm_jaccard_distance,
    then we call this as a badAnnotOrTm.
    These will be discarded unless all splices are supported.
    C) Next, we look at the asymmetric distance between this prediction and the gene intervals. 
    If this difference in these distances is over min_distance for all comparisons,
    we call this rescued and it can be incorporated.
    Otherwise, this transcript is tagged ambiguousOrFusion.

    Additionally, we look at all of the transMap projections that were filtered out and apply those
    gene names to the AlternativeGeneIds column. This is a marker of possible paralogy.
    """

    def assign_type(s):
        if tools.nameConversions.aln_id_is_denovo(s.gene):
            return True
        return False

    filtered_transmap_dict = tools.transcripts.get_gene_pred_dict(filtered_tm_gp, stranded)
    unfiltered_transmap_dict = tools.transcripts.get_gene_pred_dict(unfiltered_tm_gp, stranded)
    denovo_dict = tools.transcripts.get_gene_pred_dict(denovo_gp, stranded)

    with tools.fileOps.TemporaryFilePath() as tmp:
        # -ignoreBases = 10, to deal with potential overlap in exons
        # -conflicted flag to track exon conflicts
        cmd = ['clusterGenes', '-ignoreBases=10', '-conflicted', tmp, 'no', unfiltered_tm_gp, denovo_gp]
        print ("cmd = ", cmd)
        if not stranded:
            cmd.append(['-ignoreStrand'])
        tools.procOps.run_proc(cmd)
        cluster_df = pd.read_csv(tmp, sep='\t')

    cluster_df['is_denovo'] = cluster_df.apply(assign_type, axis=1)

    r = []
    for _, d in cluster_df.groupby('#cluster'):
        if not any(d.is_denovo):
            continue
        unfiltered_overlapping_tm_txs = set()
        filtered_overlapping_tm_txs = set()
        denovo_txs = set()
        for tx_id, is_denovo in zip(d.gene, d.is_denovo):
            if is_denovo:
                denovo_txs.add(denovo_dict[tx_id])
            elif tx_id in filtered_transmap_dict:
                filtered_overlapping_tm_txs.add(filtered_transmap_dict[tx_id])
            else:
                unfiltered_overlapping_tm_txs.add(unfiltered_transmap_dict[tx_id])

        # extract only gene names for the filtered set
        filtered_gene_ids = {tx.name2 for tx in filtered_overlapping_tm_txs}
        filtered_tx_ids = {tx.name for tx in filtered_overlapping_tm_txs}
        tx_to_gene_ids = {}
        gene_to_tx_ids = {}
        for tx in filtered_overlapping_tm_txs:
            tx_to_gene_ids[tx.name] = tx.name2
            if tx.name2 in gene_to_tx_ids:
                gene_to_tx_ids[tx.name2].append(tx.name)
            else:
                gene_to_tx_ids[tx.name2] = [tx.name]

        # print("\ndenovo_txs: ", denovo_txs)
        
        for denovo_tx in denovo_txs:
            # denovo txs may not actually overlap despite being in the same cluster 
            # (such as a cluster involving a readthrough transcript)
            # use exon conflicts to resolve cases involving readthrough transcripts
   
            # Get a list of all the conflicting exons for the current denovo transcript 
            exon_conflict_df = cluster_df.loc[cluster_df['gene']==denovo_tx.name, ['exonConflicts']]
            exon_conflict_df = exon_conflict_df.dropna()
            if not exon_conflict_df.empty:
                # Convert into a numpy array
                exon_conflicts = exon_conflict_df.to_numpy()[0][0].split(',')[:-1]
                # only add to denovo exon conflicts if it is a de novo transcript
                denovo_exon_conflicts = set()
                for tx in exon_conflicts:
                    # print("current tx: ", tx)
                    if "augPB" in tx or "augCGP" in tx: 
                        denovo_exon_conflicts.add(tx)
       
                nonoverlapping_genes = {}
                nonoverlapping_gene_ids = set()
                nonoverlapping_tx_ids = set()
                refiltered_nonoverlapping_tm_txs = set()
                # Need to deal with case where different transcripts from same gene can be both
                # overlapping and nonoverlapping 
                for conflict in exon_conflicts:
                    tx_conflict = conflict.split(':')[1]
                    if tx_conflict in tx_to_gene_ids: #otherwise this is a denovo augPB transcript
                        if tx_to_gene_ids[tx_conflict] in nonoverlapping_genes:
                            nonoverlapping_genes[tx_to_gene_ids[tx_conflict]].append(tx_conflict)
                        else:
                            nonoverlapping_genes[tx_to_gene_ids[tx_conflict]] = [tx_conflict]
                        nonoverlapping_tx_ids.add(tx_conflict)
                        refiltered_nonoverlapping_tm_txs.add(filtered_transmap_dict[tx_conflict])
                for gene in gene_to_tx_ids:
                    if gene in nonoverlapping_genes and len(nonoverlapping_genes[gene]) == len(gene_to_tx_ids[gene]):
                        nonoverlapping_gene_ids.add(gene)

                # print("nonoverlapping_gene_ids: ", nonoverlapping_gene_ids)
                overlapping_gene_ids = filtered_gene_ids - nonoverlapping_gene_ids
                # print("overlapping_gene_ids: ", overlapping_gene_ids)
                refiltered_overlapping_tm_txs = filtered_overlapping_tm_txs - refiltered_nonoverlapping_tm_txs
                # print("\tfiltered_overlapping_tm_txs: ", filtered_overlapping_tm_txs)
                # print("\trefiltered_nonoverlapping_tm_txs: ", refiltered_nonoverlapping_tm_txs)
                # print("\trefiltered_overlapping_tm_txs: ", refiltered_overlapping_tm_txs)
                # print("\toverlapping_gene_ids: ", overlapping_gene_ids)
                # print("\tnonoverlapping_gene_ids: ", nonoverlapping_gene_ids)
                # print("\tfiltered_gene_ids: ", filtered_gene_ids)

                # If there are any exon conflicts (denovo transcripts which do not share any exons), 
                # do not consider those when resolving the transcript.
                if len(overlapping_gene_ids) == 0:
                    resolved_name = resolution_method = None  # we have no matches, which means putative novel
                elif len(nonoverlapping_gene_ids) > 1: 
                    resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, refiltered_overlapping_tm_txs,
                                                                              min_distance, tm_jaccard_distance)
                elif len(filtered_gene_ids) > 1:  # we have more than one match, so resolve it
                    # print("\nfiltered gene ids 1: ", filtered_gene_ids)
                    resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, filtered_overlapping_tm_txs,
                                                                              min_distance, tm_jaccard_distance)
                elif len(filtered_gene_ids) == 1:  # yay, we have exactly one match
                    resolved_name = list(filtered_gene_ids)[0]
                    resolution_method = None
                else:
                    resolved_name = resolution_method = None  # we have no matches, which means putative novel
                # find only genes for the unfiltered set that are not present in the filtered set
                alternative_gene_ids = {tx.name2 for tx in refiltered_overlapping_tm_txs} - {resolved_name}
                # print("alternative_gene_ids: ", alternative_gene_ids)
                # print("unfiltered_overlapping_tm_txs: ", unfiltered_overlapping_tm_txs)
                # print("resolved_name: ", resolved_name)
                # print("resolution_method: ", resolution_method)
                alternative_gene_ids = ','.join(sorted(alternative_gene_ids)) if len(alternative_gene_ids) > 0 else None
                r.append([denovo_tx.name, resolved_name, alternative_gene_ids, resolution_method])
            else:
                # print("filtered_overlapping_tm_txs 2: ", filtered_overlapping_tm_txs)
                # if the are no exon conflicts, resolve like normal 
                if len(filtered_gene_ids) > 1:  # we have more than one match, so resolve it
                    # print("\nfiltered gene ids 2: ", filtered_gene_ids)
                    resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, filtered_overlapping_tm_txs,
                                                                              min_distance, tm_jaccard_distance)
                elif len(filtered_gene_ids) == 1:  # yay, we have exactly one match
                    resolved_name = list(filtered_gene_ids)[0]
                    resolution_method = None
                else:
                    resolved_name = resolution_method = None  # we have no matches, which means putative novel
                # find only genes for the unfiltered set that are not present in the filtered set
                alternative_gene_ids = {tx.name2 for tx in unfiltered_overlapping_tm_txs} - {resolved_name}
                # print("alternative_gene_ids: ", alternative_gene_ids)
                # print("unfiltered_overlapping_tm_txs: ", unfiltered_overlapping_tm_txs)
                # print("resolved_name: ", resolved_name)
                # print("resolution_method: ", resolution_method)
                alternative_gene_ids = ','.join(sorted(alternative_gene_ids)) if len(alternative_gene_ids) > 0 else None
                r.append([denovo_tx.name, resolved_name, alternative_gene_ids, resolution_method])


    combined_alternatives = pd.DataFrame(r, columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds',
                                                     'ResolutionMethod'])
    combined_alternatives = combined_alternatives.set_index('TranscriptId')
    return combined_alternatives

def resolve_multiple_genes(denovo_tx, overlapping_tm_txs, min_distance, tm_jaccard_distance):
    """
    Resolve multiple assignments based on the following rules:
    """
    # use Jaccard metric to determine if the problem lies with transMap or annotation
    tm_txs_by_gene = tools.transcripts.group_transcripts_by_name2(overlapping_tm_txs)
    tm_jaccards = [find_highest_gene_jaccard(x, y) for x, y in itertools.combinations(list(tm_txs_by_gene.values()), 2)]
    # print("overlapping_tm_txs: ", overlapping_tm_txs)
    # print("jaccards: ", tm_jaccards)
    # print("jaccard distance: ", tm_jaccard_distance)
    # print("min distance: ", min_distance)
    if any(x > tm_jaccard_distance for x in tm_jaccards):
        return None, 'badAnnotOrTm'
    # calculate asymmetric difference for this prediction
    scores = collections.defaultdict(list)
    for tx in overlapping_tm_txs:
        scores[tx.name2].append(tools.intervals.calculate_bed12_asymmetric_jaccard(denovo_tx.exon_intervals,
         tx.exon_intervals))
    # print("scores: ", scores)
    best_scores = {gene_id: max(scores[gene_id]) for gene_id in scores}
    # print("best_scores: ", best_scores)
    high_score = max(best_scores.values())
    high_gene = max(best_scores, key=lambda key: best_scores[key])
    # print("best_scores: ", best_scores)
    # print("high_score: ", high_score)
    # print("high_gene: ", high_gene)
    # This currently will rescue transcripts that are ambiguous if two share the same best score
    if all(high_score - x >= min_distance for x in best_scores.values() if (x != high_score)):
        best = sorted(iter(best_scores.items()), key=lambda gene_id_score: gene_id_score[1])[-1][0]
        # print("best = ", best)
        if best != high_gene:
            # print("ambiguousOrFusion!!!")
            return None, 'ambiguousOrFusion'
        else:
            # print ("rescued: ", best)
            return best, 'rescued'
    else:
        # print("ambiguousOrFusion")
        return None, 'ambiguousOrFusion'


def find_highest_gene_jaccard(gene_list_a, gene_list_b):
    """
    Calculates the overall distance between two sets of transcripts by finding their distinct exonic intervals and then
    measuring the Jaccard distance.
    """
    def find_interval(gene_list):
        gene_intervals = set()
        for tx in gene_list:
            gene_intervals.update(tx.exon_intervals)
        gene_intervals = tools.intervals.gap_merge_intervals(gene_intervals, 0)
        return gene_intervals

    a_interval = find_interval(gene_list_a)
    b_interval = find_interval(gene_list_b)
    return tools.intervals.calculate_bed12_jaccard(a_interval, b_interval)
