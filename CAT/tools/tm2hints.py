"""
Converts transMap genePred entries into Augustus hints
"""
import psl
import tools.procOps

cmd = ['perl', '--ep_cutoff=0', '--ep_margin=12', '--min_intron_len=50', '--start_stop_radius=5',
       '--tss_tts_radius=10', '--utrend_cutoff=10', '--in=/dev/stdin', '--out=/dev/stdout']


def tm_to_hints(tm_tx, tm_psl, ref_psl, tm_to_hints_script):
    """
    Converts a genePred transcript to hints parseable by Augustus.

    Note for anyone reading this code: GFF coordinates are 1-based, genePred coordinates are 0-based
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_psl: PslRow object for the relationship between the source transcript and genome as made by
    GenePredToFakePsl
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param tm_to_hints_script: Path to transMap2hints.pl
    :return: GFF formatted string.
    """
    cmd.insert(1, tm_to_hints_script)
    ref_starts = psl.fix_ref_q_starts(ref_psl)
    intron_vector = ['1' if psl.is_fuzzy_intron(i, tm_psl, ref_starts) else '0' for i in tm_tx.intron_intervals]
    tm_gp = '\t'.join(tm_tx.get_gene_pred())
    tm_rec = ''.join([tm_gp, '\t', ','.join(intron_vector), '\n'])
    return tools.procOps.popen_catch(cmd, tm_rec)
