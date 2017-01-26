"""
Generate all plots for the pipeline. For biotype specific plots, all plots are generated as a multi page PDF. There
is a plot for each biotype on its own, and one for the combined results.
"""
import json
import matplotlib
import logging
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.use('Agg')
import itertools
import warnings
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import pandas as pd
import tools.psl
import tools.sqlInterface
import tools.nameConversions

logger = logging.getLogger(__name__)

# suppress all warnings to make logging cleaner. The only warnings should be the chained assignment warning from pandas
# as well as the bottom == top when plots have no data.
warnings.filterwarnings('ignore')
bar_width = 0.45
paralogy_bins = [0, 1, 2, 3, 4, float('inf')]
boxplot_saturation = 0.7


def generate_plots(args):
    """
    Generates the plots.
    :param args:
    :return:
    """
    tm_data = OrderedDict([[genome, json.load(open(tgt))] for genome, tgt in args.tm_jsons.iteritems()])
    consensus_data = OrderedDict([[genome, json.load(open(tgt))] for genome, tgt in args.metrics_jsons.iteritems()])
    tm_metrics = load_tm_metrics(args.dbs)
    para_data = load_para_data(args.dbs)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.annotation_db)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.annotation_db)
    biotypes = sorted(tools.sqlInterface.get_transcript_biotypes(args.annotation_db))
    args.ordered_genomes = list(args.ordered_genomes)  # weird bug in pandas

    # hack to bring coding to the top
    try:
        biotypes.insert(0, biotypes.pop(biotypes.index('protein_coding')))
    except ValueError:
        pass

    fail_rate_plot(consensus_data, args.ordered_genomes, biotypes, args.gene_failure, args.transcript_failure)
    tx_modes_plot(consensus_data, args.ordered_genomes, args.tx_modes)
    tm_filter_plots(tm_data, args.ordered_genomes, args.transmap_filtering, biotypes)
    tm_metrics_plot(tm_metrics, args.ordered_genomes, biotypes, transcript_biotype_map, args.tm_coverage,
                    args.tm_identity)
    tm_para_plot(para_data, args.ordered_genomes, biotypes, transcript_biotype_map, args.paralogy)
    consensus_metrics_plot(consensus_data, args.ordered_genomes, biotypes, args.coverage, args.identity)
    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Annotation Support', 'Exon Annotation Support', 'Original Introns'],
                           title='Reference annotation support',
                           tgt=args.consensus_annot_support)
    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Support', 'Exon Support'],
                           title='Single species extrinsic support' if args.in_species_rna_support_only else 'All species extrinsic support',
                           tgt=args.consensus_extrinsic_support)
    completeness_plot(consensus_data, args.ordered_genomes, biotypes, args.completeness, gene_biotype_map,
                      transcript_biotype_map)
    indel_plot(consensus_data, args.ordered_genomes, args.indel)
    if 'denovo' in args:
        denovo_plot(consensus_data, args.ordered_genomes, args.denovo)
    if 'split_genes' in args:
        split_genes_plot(tm_data, args.ordered_genomes, args.split_genes)
    if 'pb_support' in args:
        pb_support_plot(consensus_data, args.ordered_genomes, args.pb_genomes, args.pb_support)
    if 'improvement' in args:
        improvement_plot(consensus_data, args.ordered_genomes, args.improvement)


###
# Load metrics from transMap PSLs
###


def load_tm_metrics(dbs):
    """Loads transMap data from PSLs"""
    tm_metrics = {'transMap Coverage': OrderedDict(), 'transMap Identity': OrderedDict()}
    tm_name_map = {'TransMapCoverage': 'transMap Coverage', 'TransMapIdentity': 'transMap Identity'}
    for genome, db_path in dbs.iteritems():
        session = tools.sqlInterface.start_session(db_path)
        table = tools.sqlInterface.TmEval
        for classifier in ['TransMapCoverage', 'TransMapIdentity']:
            query = session.query(table.AlignmentId, table.value).filter(table.classifier == classifier)
            tm_metrics[tm_name_map[classifier]][genome] = dict(query.all())
    return tm_metrics


def load_para_data(dbs):
    para_data = OrderedDict()
    for genome, db_path in dbs.iteritems():
        session = tools.sqlInterface.start_session(db_path)
        table = tools.sqlInterface.TmEval
        query = session.query(table.TranscriptId, table.value).filter(table.classifier == 'Paralogy')
        para_data[genome] = dict(query.all())
    return para_data


###
# Plots
###


def tm_filter_plots(tm_data, ordered_genomes, tgt, biotypes):
    """Plots for the transMap filtering process. The data munging is a huge mess, sorry."""
    data = OrderedDict()
    for genome in ordered_genomes:
        for biotype in tm_data[genome]['Paralogy']:
            data[(genome, biotype)] = OrderedDict(sorted(tm_data[genome]['Paralogy'][biotype].items()))
    df = pd.DataFrame.from_dict(data)
    df = df.transpose().reset_index()
    df.columns = ['genome', 'biotype'] + list(df.columns[2:])
    df['genome'] = pd.Categorical(df.genome, ordered_genomes, ordered=True)
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        combined_df = df.groupby('genome').apply(sum).transpose()
        combined_df.columns = ordered_genomes
        combined_df = combined_df.iloc[3:][::-1]
        combined_df = combined_df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        title_string = 'transMap paralogous alignment resolution'
        if len(ordered_genomes) > 1:
            generic_stacked_barplot(combined_df, pdf, title_string, combined_df.index, 'Number of transcripts',
                                    combined_df.columns, 'method', bbox_to_anchor=(1.2, 0.7))
        else:
            generic_barplot(combined_df, pdf, ordered_genomes[0], 'Number of transcripts', title_string,
                            x=combined_df.index, y=ordered_genomes[0])
        for biotype in biotypes:
            b_df = biotype_filter(df, biotype)
            if b_df is None:
                continue
            b_df = b_df.transpose()
            try:  # we don't have all genomes. probably a tiny biotype
                b_df.columns = ordered_genomes
            except ValueError:
                b_df = b_df.transpose()
                b_df.genome = pd.Categorical(b_df.genome, ordered_genomes, ordered=True)
                missing = pd.DataFrame([[genome, biotype, 0, 0, 0, 0] for genome in b_df.genome.cat.categories if genome
                                        not in list(b_df.genome)])
                missing.columns = b_df.columns
                b_df = b_df.append(missing).sort_values('genome').transpose()
            b_df = b_df.iloc[3:][::-1]
            b_df = b_df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
            if b_df.sum(numeric_only=True).sum() == 0:
                continue
            title_string = 'transMap paralogous alignment resolution\nfor biotype {}'.format(biotype)
            if len(ordered_genomes) > 1:
                generic_stacked_barplot(b_df, pdf, title_string, b_df.index, 'Number of transcripts',
                                        b_df.columns, 'method', bbox_to_anchor=(1.2, 0.7))
            else:
                generic_barplot(b_df, pdf, ordered_genomes[0], 'Number of transcripts', title_string,
                                x=b_df.index, y=ordered_genomes[0])


def tm_metrics_plot(tm_metrics, ordered_genomes, biotypes, transcript_biotype_map, tm_coverage_tgt, tm_identity_tgt):
    """plots for transMap coverage, identity"""
    tm_iter = zip(*[['transMap Coverage', 'transMap Identity'],
                    [tm_coverage_tgt, tm_identity_tgt]])
    for mode, tgt in tm_iter:
        df = dict_to_df_with_biotype(tm_metrics[mode], transcript_biotype_map)
        df = pd.melt(df, id_vars='biotype', value_vars=ordered_genomes)
        df.columns = ['biotype', 'genome', mode]
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome')


def consensus_metrics_plot(consensus_data, ordered_genomes, biotypes, coverage_tgt, identity_tgt):
    """plots for consensus coverage, identity, score"""
    cons_iter = zip(*[['Coverage', 'Identity'],
                      [coverage_tgt, identity_tgt]])
    for mode, tgt in cons_iter:
        df = json_to_df_with_biotype(consensus_data, mode)
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome')


def consensus_support_plot(consensus_data, ordered_genomes, biotypes, modes, title, tgt):
    """grouped violin plots of original intron / intron annotation / exon annotation support"""
    def adjust_plot(g, this_title):
        g.set_xticklabels(rotation=90)
        g.fig.suptitle(this_title)
        g.fig.subplots_adjust(top=0.9)
        for ax in g.axes.flat:
            ax.set_ylabel('Percent supported')
            ax.set_ylim(-1, 101)

    dfs = []
    for i, mode in enumerate(modes):
        df = json_to_df_with_biotype(consensus_data, mode)
        if i > 0:
            df = df[mode]
        dfs.append(df)
    df = pd.concat(dfs, axis=1)
    df = pd.melt(df, value_vars=modes, id_vars=['genome', 'biotype'])
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        if len(ordered_genomes) > 1:
            g = sns.factorplot(data=df, y='value', x='genome', col='variable', col_wrap=2, kind='violin', sharex=True,
                               sharey=True, row_order=ordered_genomes, cut=0)
        else:
            g = sns.factorplot(data=df, y='value', x='variable', kind='violin', sharex=True,
                               sharey=True, row_order=ordered_genomes, cut=0)
        adjust_plot(g, title)
        multipage_close(pdf, tight_layout=False)
        title += ' for {}'
        for biotype in biotypes:
            this_title = title.format(biotype)
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                if len(ordered_genomes) > 1:
                    g = sns.factorplot(data=biotype_df, y='value', x='genome', col='variable', col_wrap=2,
                                       kind='violin', sharex=True, sharey=True, row_order=ordered_genomes, cut=0)
                else:
                    g = sns.factorplot(data=df, y='value', x='variable', kind='violin', sharex=True,
                                       sharey=True, row_order=ordered_genomes, cut=0)
                adjust_plot(g, this_title)
                multipage_close(pdf, tight_layout=False)


def tm_para_plot(para_data, ordered_genomes, biotypes, transcript_biotype_map, para_tgt):
    """transMap paralogy plots"""
    def generate_hists(ordered_genomes, df):
        hists = OrderedDict([genome, np.roll(np.histogram(df[genome].fillna(0), paralogy_bins)[0], -1)]
                             for genome in ordered_genomes)
        hists_df = pd.DataFrame.from_dict(hists)
        return hists_df

    df = dict_to_df_with_biotype(para_data, transcript_biotype_map)
    legend_labels = ['= {}'.format(x) for x in paralogy_bins[1:-2]] + [u'\u2265 {}'.format(paralogy_bins[-2])] + \
                    ['= {}'.format(paralogy_bins[0])]
    title_string = 'Proportion of transcripts that have multiple alignments'
    biotype_title_string = 'Proportion of {} transcripts that have multiple alignments'
    plot_fn = generic_unstacked_barplot if len(df.columns) <= 5 else generic_stacked_barplot
    box_label = 'Number of\nalignments'
    with para_tgt.open('w') as outf, PdfPages(outf) as pdf:
        hists_df = generate_hists(ordered_genomes, df)
        plot_fn(hists_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes, box_label)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                hists_df = generate_hists(ordered_genomes, biotype_df)
                title_string = biotype_title_string.format(biotype)
                plot_fn(hists_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes, box_label)


def fail_rate_plot(consensus_data, ordered_genomes, biotypes, gene_fail_plot_tgt, transcript_fail_plot_tgt):
    for mode, tgt in zip(*[['Gene', 'Transcript'], [gene_fail_plot_tgt, transcript_fail_plot_tgt]]):
        base_title = '{} outcomes in consensus gene set'.format(mode)
        ylabel = 'Number of {}s'.format(mode.lower())
        fail_df = json_biotype_counter_to_df(consensus_data, '{} Failed'.format(mode))
        missing_df = json_biotype_counter_to_df(consensus_data, '{} Missing'.format(mode))
        try:
            fail_df.columns = ['biotype', '{} Failed'.format(mode), 'genome']
        except ValueError:  # we have nothing here
            fail_df = pd.DataFrame([[biotype, 0, genome] for biotype, genome in itertools.product(biotypes,
                                                                                                  ordered_genomes)])
            fail_df.columns = ['biotype', '{} Failed'.format(mode), 'genome']
        try:
            missing_df.columns = ['biotype', '{} Missing'.format(mode), 'genome']
        except ValueError:  # we have nothing here
            missing_df = pd.DataFrame([[biotype, 0, genome] for biotype, genome in itertools.product(biotypes,
                                                                                                     ordered_genomes)])
            missing_df.columns = ['biotype', '{} Missing'.format(mode), 'genome']
        df = pd.merge(fail_df, missing_df, on=['genome', 'biotype'], how='outer')
        if len(df) == 0:
            continue
        df = pd.melt(df, id_vars=['biotype', 'genome'], value_vars=['{} Failed'.format(mode),
                                                                    '{} Missing'.format(mode)])
        df.columns = ['biotype', 'Genome', 'Outcome', 'value']
        with tgt.open('w') as outf, PdfPages(outf) as pdf:
            tot_df = pd.melt(df.groupby(by=['Genome', 'Outcome']).aggregate(np.sum).transpose())
            tot_df.columns = ['Genome', 'Outcome', 'value']
            if len(ordered_genomes) > 1:
                generic_barplot(tot_df, pdf, 'Genome', ylabel, base_title, x='Genome', y='value', hue='Outcome',
                                row_order=ordered_genomes)
            else:
                generic_barplot(tot_df, pdf, 'Outcome', ylabel, base_title, x='Outcome', y='value')
            for biotype in biotypes:
                biotype_df = biotype_filter(df, biotype)
                if biotype_df is None:
                    continue
                title = base_title + '\nfor biotype {}'.format(biotype)
                if len(ordered_genomes) > 1:
                    generic_barplot(biotype_df, pdf, 'Genome', ylabel, title, x='Genome', y='value', hue='Outcome',
                                    row_order=ordered_genomes)
                else:
                    generic_barplot(biotype_df, pdf, 'Outcome', ylabel, title, x='Outcome', y='value')


def tx_modes_plot(consensus_data, ordered_genomes, tx_mode_plot_tgt):
    ordered_groups = ['transMap', 'transMap+TM', 'transMap+TMR', 'transMap+TM+TMR', 'TM', 'TMR', 'TM+TMR', 'CGP', 'PB',
                      'Other']
    ordered_groups = OrderedDict([[frozenset(x.split('+')), x] for x in ordered_groups])

    def split_fn(s):
        return ordered_groups.get(frozenset(s['Transcript Modes'].replace('aug', '').split(',')), 'Other')

    modes_df = json_biotype_counter_to_df(consensus_data, 'Transcript Modes')
    df = modes_df.pivot(index='genome', columns='Transcript Modes').transpose().reset_index()
    df['Modes'] = df.apply(split_fn, axis=1)
    df = df[['Modes'] + ordered_genomes]
    ordered_values = [x for x in ordered_groups.itervalues() if x in set(df['Modes'])]
    with tx_mode_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        title_string = 'Transcript modes in protein coding consensus gene set'
        ylabel = 'Number of transcripts'
        if len(ordered_genomes) > 1:
            df['Ordered Modes'] = pd.Categorical(df['Modes'], ordered_values, ordered=True)
            df = df.sort_values('Ordered Modes')
            df = df[['Ordered Modes'] + ordered_genomes].set_index('Ordered Modes')
            generic_stacked_barplot(df, pdf, title_string, df.index, ylabel, ordered_genomes, 'Transcript mode(s)',
                                    bbox_to_anchor=(1.25, 0.7))

        else:
            generic_barplot(pd.melt(df, id_vars='Modes'), pdf, 'Transcript mode(s)', ylabel, title_string, x='Modes',
                            y='value', order=ordered_values)


def denovo_plot(consensus_data, ordered_genomes, denovo_tgt):
    with denovo_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_nested_counter_to_df(consensus_data, 'denovo')
        # fix column names because json_biotype_nested_counter_to_df makes assumptions
        df.columns = ['Result', 'Number of transcripts', 'Augustus mode', 'genome']
        has_pb = len(set(df['Augustus mode'])) == 2
        if len(ordered_genomes) > 1:
            if has_pb is True:
                ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result',
                                    hue='Augustus mode', col_wrap=2, row_order=ordered_genomes, sharex=True,
                                    sharey=False)
            else:
                ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result',
                                    col_wrap=2, row_order=ordered_genomes, sharex=True, sharey=False)
        else:
            if has_pb is True:
                ax = sns.factorplot(data=df, x='Result', y='Number of transcripts', kind='bar', hue='Augustus mode')
            else:
                ax = sns.factorplot(data=df, x='Result', y='Number of transcripts', kind='bar')
        ax.set_xticklabels(rotation=90)
        ax.fig.suptitle('Incorporation of de-novo predictions')
        ax.fig.subplots_adjust(top=0.9)
        multipage_close(pdf, tight_layout=False)


def split_genes_plot(tm_data, ordered_genomes, split_plot_tgt):
    with split_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_counter_to_df(tm_data, 'Split Genes')
        df.columns = ['category', 'count', 'genome']
        title = 'Split genes'
        if len(ordered_genomes) > 1:
            g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', col='category', xlabel='', col_wrap=2,
                                sharey=False, ylabel='Number of transcripts or genes', row_order=ordered_genomes,
                                title=title)
        else:
            g = generic_barplot(pdf=pdf, data=df, x='category', y='count', ylabel='Number of transcripts or genes',
                                title=title, xlabel='Category')


def pb_support_plot(consensus_data, ordered_genomes, pb_genomes, pb_support_tgt):
    with pb_support_tgt.open('w') as outf, PdfPages(outf) as pdf:
        pb_genomes = [x for x in ordered_genomes if x in pb_genomes]  # fix order
        df = json_biotype_counter_to_df(consensus_data, 'IsoSeq Transcript Validation')
        df.columns = ['IsoSeq Transcript Validation', 'Number of transcripts', 'genome']
        ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', hue='IsoSeq Transcript Validation',
                            kind='bar', row_order=pb_genomes)
        ax.set_xticklabels(rotation=90)
        ax.fig.suptitle('Isoforms validated by at least one IsoSeq read')
        multipage_close(pdf, tight_layout=False)


def completeness_plot(consensus_data, ordered_genomes, biotypes, completeness_plot_tgt, gene_biotype_map,
                      transcript_biotype_map):
    def adjust_plot(g, gene_count, tx_count):
        for ax, c in zip(*[g.axes[0], [gene_count, tx_count]]):
            _ = ax.set_ylim(0, c)
            ax.spines['top'].set_edgecolor('#e74c3c')
            ax.spines['top'].set_linewidth(2)
            ax.spines['top'].set_visible(True)
            ax.spines['top'].set_linestyle('dashed')

    df = json_grouped_biotype_nested_counter_to_df(consensus_data, 'Completeness')
    with completeness_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        tot_df = df.groupby(by=['genome', 'category']).aggregate(np.sum).reset_index()
        tot_df = sort_long_df(tot_df, ordered_genomes)
        title = 'Number of comparative genes/transcripts present'
        g = generic_barplot(pdf=pdf, data=tot_df, x='genome', y='count', col='category', xlabel='',
                            sharey=False, ylabel='Number of genes/transcripts', title=title,
                            col_order=['Gene', 'Transcript'], close=False, palette=choose_palette(ordered_genomes))
        adjust_plot(g, len(gene_biotype_map), len(transcript_biotype_map))
        multipage_close(pdf, tight_layout=False)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                biotype_df = sort_long_df(biotype_df, ordered_genomes)
                gene_biotype_count = len({i for i, b in gene_biotype_map.iteritems() if b == biotype})
                tx_biotype_count = len({i for i, b in transcript_biotype_map.iteritems() if b == biotype})
                title = 'Number of comparative genes/transcripts present for biotype {}'.format(biotype)
                g = generic_barplot(pdf=pdf, data=biotype_df, x='genome', y='count', col='category', xlabel='',
                                    sharey=False, ylabel='Number of genes/transcripts',
                                    title=title, col_order=['Gene', 'Transcript'], close=False,
                                    palette=choose_palette(ordered_genomes))
                adjust_plot(g, gene_biotype_count, tx_biotype_count)
                multipage_close(pdf, tight_layout=False)


def improvement_plot(consensus_data, ordered_genomes, improvement_tgt):
    def do_kdeplot(x, y, ax, n_levels=None, bw='scott'):
        try:
            sns.kdeplot(x, y, ax=ax, cut=0, cmap='Purples_d', shade=True, shade_lowest=False, n_levels=n_levels, bw=bw)
        except ValueError:
            logger.warning('Unable to do a KDE fit to AUGUSTUS improvement.')
            pass

    with improvement_tgt.open('w') as outf, PdfPages(outf) as pdf, sns.axes_style("whitegrid"):
        for genome in ordered_genomes:
            data = pd.DataFrame(consensus_data[genome]['Evaluation Improvement']['changes'])
            unchanged = consensus_data[genome]['Evaluation Improvement']['unchanged']
            if len(data) == 0:
                continue
            data.columns = ['transMap original introns',
                            'transMap intron annotation support',
                            'transMap intron RNA support',
                            'Original introns',
                            'Intron annotation support',
                            'Intron RNA support',
                            'transMap alignment goodness',
                            'Alignment goodness']
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2)
            for ax in [ax1, ax2, ax3]:  # goodness plots are allowed to auto-set scale
                ax.set_xlim(0, 100)
                ax.set_ylim(0, 100)
            goodness_min = min(data['Alignment goodness'])
            ax4.set_xlim(goodness_min, 100)
            ax4.set_ylim(goodness_min, 100)
            do_kdeplot(data['transMap original introns'], data['Original introns'], ax1, n_levels=25, bw=2)
            sns.regplot(x=data['transMap original introns'], y=data['Original introns'], ax=ax1,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7}, fit_reg=False)
            do_kdeplot(data['transMap intron annotation support'], data['Intron annotation support'], ax2,
                       n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron annotation support'], y=data['Intron annotation support'], ax=ax2,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7}, fit_reg=False)
            do_kdeplot(data['transMap intron RNA support'], data['Intron RNA support'], ax3, n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron RNA support'], y=data['Intron RNA support'], ax=ax3,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7}, fit_reg=False)
            do_kdeplot(data['transMap alignment goodness'], data['Alignment goodness'], ax4, n_levels=20, bw=1)
            sns.regplot(x=data['transMap alignment goodness'], y=data['Alignment goodness'], ax=ax4,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7}, fit_reg=False)
            fig.suptitle('AUGUSTUS metric improvements for {:,} transcripts in {}.\n'
                         '{:,} transMap transcripts were chosen.'.format(len(data), genome, unchanged))
            for ax in [ax1, ax2, ax3, ax4]:
                ax.set(adjustable='box-forced', aspect='equal')
            fig.subplots_adjust(hspace=0.3)
            multipage_close(pdf, tight_layout=False)


def indel_plot(consensus_data, ordered_genomes, indel_plot_tgt):
    with indel_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        tm_df = pd.concat([pd.DataFrame.from_dict(consensus_data[genome]['transMap Indels'], orient='index').T
                           for genome in ordered_genomes])
        tm_df['genome'] = ordered_genomes
        tm_df['transcript set'] = ['transMap'] * len(tm_df)
        consensus_df = pd.concat([pd.DataFrame.from_dict(consensus_data[genome]['Consensus Indels'], orient='index').T
                                  for genome in ordered_genomes])
        consensus_df['genome'] = ordered_genomes
        consensus_df['transcript set'] = ['Consensus'] * len(consensus_df)
        df = pd.concat([consensus_df, tm_df])
        df = pd.melt(df, id_vars=['genome', 'transcript set'],
                     value_vars=['CodingDeletion', 'CodingInsertion', 'CodingMult3Indel'])
        df.columns = ['Genome', 'Transcript set', 'Type', 'Percent of transcripts']
        g = sns.factorplot(data=df, x='Genome', y='Percent of transcripts', col='Transcript set',
                           hue='Type', kind='bar', row_order=ordered_genomes,
                           col_order=['transMap', 'Consensus'])
        g.set_xticklabels(rotation=90)
        g.fig.subplots_adjust(top=.8)
        g.fig.suptitle('Coding indels')
        multipage_close(pdf, tight_layout=False)


###
# shared plotting functions
###


def cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=None, y=None, xlabel=None):
    """violin plots for coverage and identity."""
    if xlabel is None:
        xlabel = 'Percent {}'.format(mode)
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        title = 'Overall {}'.format(mode)
        xmin = int(min(df[mode]))
        horizontal_violin_plot(df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=(xmin, 100))
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = '{} for biotype {}'.format(mode, biotype)
                xmin = int(min(df[mode]))
                horizontal_violin_plot(biotype_df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=(xmin, 100))


###
# generic plotting functions
###


def generic_barplot(data, pdf, xlabel, ylabel, title, row_order=None, x=None, y=None, hue=None, hue_order=None,
                    order=None, col=None, col_wrap=None, sharex=True, sharey=True, col_order=None, palette=None,
                    close=True):
    g = sns.factorplot(data=data, x=x, y=y, hue=hue, ci=None, kind='bar', hue_order=hue_order, row_order=row_order,
                       col=col, col_wrap=col_wrap, sharex=sharex, sharey=sharey, col_order=col_order, palette=palette,
                       order=order)
    g.set_xticklabels(rotation=90)
    g.fig.suptitle(title)
    g.fig.subplots_adjust(top=.8)
    g.set_axis_labels(xlabel, ylabel)
    try:  # depending on columns, axes could be flat or not
        axes = list(itertools.chain.from_iterable(g.axes))
    except TypeError:
        axes = g.axes
    for ax in axes:
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=10, steps=[1, 2, 5, 10], integer=True))
        ax.margins(y=0.15)
        ax.autoscale(enable=True, axis='y', tight=False)
        ax.set_ylim(0, ax.get_ylim()[1])
    if close is True:
        multipage_close(pdf, tight_layout=False)
    return g


def horizontal_violin_plot(data, ordered_genomes, title, xlabel, pdf, hue=None, x=None, y=None, xlim=None):
    """not so generic function that specifically produces a paired boxplot/violinplot"""
    fig, ax = plt.subplots()
    sns.violinplot(data=data, x=x, y=y, hue=hue, order=ordered_genomes, palette=choose_palette(ordered_genomes),
                   saturation=boxplot_saturation, orient='h', cut=0, scale='count', ax=ax)
    fig.suptitle(title)
    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    multipage_close(pdf, tight_layout=False)


def _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor):
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=bbox_to_anchor, frameon=True,
               title=box_label)
    ax.set_title(title_string)
    ax.set_ylabel(ylabel)
    set_ticks(names, ax)
    ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
    sns.despine(top=True, right=True)
    multipage_close(pdf)


def generic_unstacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label,
                              bbox_to_anchor=(1.12, 0.7)):
    fig, ax = plt.subplots()
    bars = []
    shorter_bar_width = bar_width / len(df)
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)) + shorter_bar_width * i, d, shorter_bar_width,
                           color=sns.color_palette()[i], linewidth=0.0))
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor)


def generic_stacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label, bbox_to_anchor=(1.12, 0.7)):
    fig, ax = plt.subplots()
    bars = []
    cumulative = np.zeros(len(df.columns))
    color_palette = choose_palette(legend_labels)
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)), d, bar_width, bottom=cumulative,
                           color=color_palette[i], linewidth=0.0))
        cumulative += d
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor)


###
# Shared functions
###


def json_flat_to_df(consensus_data, key):
    """converts cases where we have exactly genome:value pairs"""
    r = []
    for genome, d in consensus_data.iteritems():
        r.append([genome, d[key]])
    return pd.DataFrame(r)


def json_to_df_with_biotype(consensus_data, key):
    """converts JSON entries with many transcripts, such as those for coverage/identity"""
    dfs = []
    for genome, d in consensus_data.iteritems():
        for biotype, vals in d[key].iteritems():
            df = pd.DataFrame(vals)
            if len(df) > 0:
                df.columns = [key]
                df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
            dfs.append(df)
    return pd.concat(dfs)


def json_biotype_nested_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the first level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.iteritems():
        for biotype, vals in d[key].iteritems():
            df = pd.DataFrame(vals.items())
            if len(df) > 0:
                df.columns = [key, 'count']
                df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
                dfs.append(df)
    return pd.concat(dfs)


def json_grouped_biotype_nested_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the second level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.iteritems():
        for group, vals in d[key].iteritems():
            df = pd.DataFrame(vals.items())
            if len(df) > 0:
                df.columns = ['biotype', 'count']
                df = df.assign(category=[group] * len(df), genome=[genome] * len(df))
            dfs.append(df)
    return pd.concat(dfs)


def json_biotype_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the first level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.iteritems():
        vals = consensus_data[genome][key]
        df = pd.DataFrame(vals.items())
        if len(df) > 0:
            df.columns = [key, 'count']
            df = df.assign(genome=[genome] * len(df))
        dfs.append(df)
    return pd.concat(dfs)


def dict_to_df_with_biotype(data, transcript_biotype_map):
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.iteritems()]))
    try:
        df['biotype'] = [transcript_biotype_map[tx] for tx in df.index]
    except KeyError:
        # try removing names
        df['biotype'] = [transcript_biotype_map[tools.nameConversions.strip_alignment_numbers(tx)] for tx in df.index]
    return df


def biotype_filter(df, biotype):
    df = df[df.biotype == biotype]
    return df if len(df) > 0 else None


def multipage_close(pdf, tight_layout=True):
    """convenience function for closing up a pdf page"""
    if tight_layout:
        plt.tight_layout()
    pdf.savefig(bbox_inches='tight')
    plt.close('all')


def choose_palette(ordered_genomes):
    """choose palette in cases where genomes get different colors"""
    if len(ordered_genomes) <= 6:
        return sns.color_palette()
    else:
        return sns.color_palette("Set2", len(ordered_genomes))


def set_ticks(names, ax, nbins=10.0):
    ax.margins(y=0.15)
    ax.autoscale(enable=True, axis='y', tight=False)
    ax.set_ylim(0, plt.ylim()[1])
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nbins, steps=[1, 2, 5, 10], integer=True))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(len(names)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.xaxis.set_ticklabels(names, rotation=90)


def sort_long_df(df, ordered_genomes):
    """sorts a long form dataframe by ordered genomes"""
    ordered_index = dict(zip(ordered_genomes, range(len(ordered_genomes))))
    df['order'] = df['genome'].map(ordered_index)
    df = df.sort_values('order')
    return df.drop('order', axis=1)
