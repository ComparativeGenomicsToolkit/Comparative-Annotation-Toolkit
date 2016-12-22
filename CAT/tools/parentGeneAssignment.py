"""
A set of functions to perform parental gene assignment in the AugustusPB/AugustusCGP modules
"""
import pandas as pd
import itertools
import collections
import procOps
import fileOps
import dataOps
import sqlInterface
import transcripts


def create_chrom_dict(tx_dict):
    """
    Split up a dictionary of Transcript objects by chromosome
    """
    chrom_dict = collections.defaultdict(dict)
    for tx_id, tx in tx_dict.iteritems():
        chrom_dict[tx.chromosome][tx_id] = tx
    return chrom_dict


def find_tm_overlaps(cgp_tx, tm_tx_dict):
    """Find overlap with transMap transcripts first on a genomic scale then an exonic scale"""
    r = []
    for tx in tm_tx_dict.itervalues():
        if tx.interval.intersection(cgp_tx.interval) is not None:
            # make sure that we have exon overlap
            if ensure_exon_overlap(tx, cgp_tx) is True:
                r.append(tx)
    return r


def ensure_exon_overlap(tx, cgp_tx):
    """Do these two transcripts have at least 1 exonic base of overlap?"""
    for tm_exon in tx.exon_intervals:
        for cgp_exon in cgp_tx.exon_intervals:
            if tm_exon.overlap(cgp_exon) is True:
                return True
    return False


def resolve_multiple_genes(cgp_tx, overlapping_tm_txs, gene_biotype_map, min_distance=0.2):
    """
    Resolve multiple assignments based on the following rules:
    If a CGP has multiple assignments, see if only 1 is protein_coding. If so, assign it.
    If the difference in Jaccard scores is >=min_distance, assign it to the higher score.
    If neither of these are satisfiable, discard the transcript.
    """
    filtered = [tm_tx for tm_tx in overlapping_tm_txs if gene_biotype_map[tm_tx.name2] == 'protein_coding']
    if len(filtered) == 0:
        return None
    elif len(filtered) == 1:
        return filtered[0].name2
    # attempt the Jaccard metric
    scores = calculate_jaccard(cgp_tx, filtered)
    high_score = max(scores.itervalues())
    if all(high_score - x >= min_distance for x in scores.itervalues() if x != high_score):
        best = sorted(scores.iteritems(), key=lambda (gene_id, score): score)[-1]
        results = best[0]
    else:
        results = None
    return results


def calculate_jaccard(cgp_tx, filtered_overlapping_tm_txs):
    """calculates the Jaccard distance metric."""
    results = collections.defaultdict(float)
    with fileOps.TemporaryFilePath() as cgp:
        with open(cgp, 'w') as outf:
            fileOps.print_row(outf, cgp_tx.get_bed())
        with fileOps.TemporaryFilePath() as tm:
            for tm_tx in filtered_overlapping_tm_txs:
                with open(tm, 'w') as outf:
                    fileOps.print_row(outf, tm_tx.get_bed())
                cmd = ['bedtools', 'jaccard', '-s', '-a', cgp, '-b', tm]
                r = procOps.call_proc_lines(cmd)
                j = float(r[-1].split()[-2])
                results[tm_tx.name2] = max(results[tm_tx.name2], j)
    return results


###
# Parent gene assignment Toil section
###


def assign_parents(job, ref_db_file_id, filtered_tm_gp_file_id, unfiltered_tm_gp_file_id, joined_gff_file_id):
    """
    Main function for assigning parental genes. Parental gene assignment methodology:
    A) Each CGP transcript is evaluated for overlapping any transMap transcripts. Overlap is defined as having at least
    1 exonic base shared.
    B) If a CGP transcript is assigned to more than one gene, then this is attempted to be resolved by first:
        1) Determine if only one of the possible genes is coding. If so, make the assignment.
        2) Use the Jaccard distance metric to attempt to find one best hit.
        If none of the above are satisfied, the CGP transcript is discarded as unresolvable.

    TODO: this process operates on genePreds, and so I needlessly transform to genePred then back. This could introduce
    errors. The better way would be a proper GFF parser. Mark's Gff3 parser does not work.
    """
    ref_genome_db = job.fileStore.readGlobalFile(ref_db_file_id)
    gene_biotype_map = sqlInterface.get_gene_biotype_map(ref_genome_db)
    filtered_tm_gp_file = job.fileStore.readGlobalFile(filtered_tm_gp_file_id)
    unfiltered_tm_gp_file = job.fileStore.readGlobalFile(unfiltered_tm_gp_file_id)
    filtered_transmap_dict = transcripts.get_gene_pred_dict(filtered_tm_gp_file)
    unfiltered_transmap_dict = transcripts.get_gene_pred_dict(unfiltered_tm_gp_file)
    filtered_ids = unfiltered_transmap_dict.viewkeys() - filtered_transmap_dict.viewkeys()

    # convert GFF to genePred
    cgp_gp = fileOps.get_tmp_toil_file()
    joined_gff = job.fileStore.readGlobalFile(joined_gff_file_id)
    cmd = ['gtfToGenePred', '-genePredExt', joined_gff, cgp_gp]
    procOps.run_proc(cmd)
    cgp_dict = transcripts.get_gene_pred_dict(cgp_gp)
    tm_chrom_dict = create_chrom_dict(unfiltered_transmap_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)

    final_gps = []
    for chrom, tm_tx_by_chromosome in tm_chrom_dict.iteritems():
        for cgp_chunk in dataOps.grouper(cgp_chrom_dict[chrom].iteritems(), 70):
            j = job.addChildJobFn(assign_parent_chunk, tm_tx_by_chromosome, cgp_chunk, gene_biotype_map, filtered_ids)
            final_gps.append(j.rv())
    return job.addFollowOnJobFn(merge_parent_assignment_chunks, final_gps).rv()


def assign_parent_chunk(job, tm_tx_by_chromosome, cgp_chunk, gene_biotype_map, filtered_ids):
    """
    Runs a chunk of CGP transcripts on the same chromosome as all transMap transcripts in tm_tx_by_chromosome

    :param tm_tx_by_chromosome: dictionary of GenePredTranscript objects all on the same chromosome
    :param cgp_chunk: Iterable of (cgp_tx_id, cgp_tx) tuples to be analyzed
    :param gene_biotype_map: dictionary mapping gene IDs to biotype
    :param filtered_ids: Set of transcript IDs which were filtered out
    :return: list to be turned into a dataframe. See merge_parent_assignment_chunks()
    """
    r = []
    for cgp_tx_id, cgp_tx in cgp_chunk:
        # find the names of both filtered and unfiltered transMap transcript IDs that overlap
        unfiltered_overlapping_tm_txs = find_tm_overlaps(cgp_tx, tm_tx_by_chromosome)
        filtered_overlapping_tm_txs = {tx for tx in unfiltered_overlapping_tm_txs if tx.name not in filtered_ids}
        # extract only gene names for the filtered set
        filtered_gene_ids = {tx.name2 for tx in filtered_overlapping_tm_txs}
        if len(filtered_gene_ids) > 0:  # we have at least one match
            if len(filtered_gene_ids) > 1:  # we have more than one match, so resolve it
                resolved_name = resolve_multiple_genes(cgp_tx, filtered_overlapping_tm_txs, gene_biotype_map)
            elif len(filtered_gene_ids) == 1:  # yay, we have exactly one match
                resolved_name = list(filtered_gene_ids)[0]
            if resolved_name is None:
                continue  # don't save transcripts that can't be resolved
        else:
            resolved_name = None  # we have no matches, which means putative novel

        # find only genes for the unfiltered set that are not present in the filtered set
        alternative_gene_ids = {tx.name2 for tx in unfiltered_overlapping_tm_txs} - {resolved_name}
        alternative_gene_ids = ','.join(alternative_gene_ids) if len(alternative_gene_ids) > 0 else None
        r.append([cgp_tx, resolved_name, alternative_gene_ids])
    return r


def merge_parent_assignment_chunks(job, final_gps):
    """
    Merge the chunks of transcripts produced by assign_parent_chunk, converting back to GFF
    :param final_gps: paired list of lists of GenePredTranscript objects and alternative tx dicts
    :return: fileStore ID to a output GFF and a pandas DataFrame to be written to sql
    """
    out_gp = fileOps.get_tmp_toil_file()
    out_gff = fileOps.get_tmp_toil_file()
    combined_alternatives = []
    with open(out_gp, 'w') as outf:
        for tx_rec, assigned_id, alt_ids in itertools.chain.from_iterable(final_gps):
            fileOps.print_row(outf, tx_rec.get_gene_pred())
            combined_alternatives.append([tx_rec.name, assigned_id, alt_ids])
    cmd = ['genePredToGtf', '-utr', '-honorCdsStat', '-source=AugustusCGP', 'file', out_gp, out_gff]
    procOps.run_proc(cmd)
    combined_alternatives = pd.DataFrame(combined_alternatives)
    combined_alternatives.columns = ['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds']
    combined_alternatives = combined_alternatives.set_index('TranscriptId')
    return job.fileStore.writeGlobalFile(out_gff), combined_alternatives
