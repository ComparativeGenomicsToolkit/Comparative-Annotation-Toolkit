"""
Generate all plots for the pipeline. For biotype specific plots, all plots are generated as a multi page PDF. There
is a plot for each biotype on its own, and one for the combined results.
"""
import json
import matplotlib
import itertools
import warnings
matplotlib.use('Agg')
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
    fail_rate_plot(consensus_data, args.ordered_genomes, biotypes, args.gene_failure, args.transcript_failure)
    completeness_plot(consensus_data, args.ordered_genomes, biotypes, args.completeness, gene_biotype_map,
                      transcript_biotype_map)
    tx_modes_plot(consensus_data, args.ordered_genomes, args.tx_modes)
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
        generic_stacked_barplot(combined_df, pdf, title_string, combined_df.index, 'Number of transcripts',
                                combined_df.columns, 'method', bbox_to_anchor=(1.2, 0.7))
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
            title_string = 'transMap paralogous alignment resolution for biotype {}'.format(biotype)
            generic_stacked_barplot(b_df, pdf, title_string, b_df.index, 'Number of transcripts',
                                    b_df.columns, 'method', bbox_to_anchor=(1.2, 0.7))


def tm_metrics_plot(tm_metrics, ordered_genomes, biotypes, transcript_biotype_map, tm_coverage_tgt, tm_identity_tgt):
    """plots for transMap coverage, identity"""
    tm_iter = zip(*[['transMap Coverage', 'transMap Identity'],
                    [tm_coverage_tgt, tm_identity_tgt]])
    for mode, tgt in tm_iter:
        df = dict_to_df_with_biotype(tm_metrics[mode], transcript_biotype_map)
        for genome in ordered_genomes:
            df[genome] = df[genome]
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, xlim=(-1, 101))


def consensus_metrics_plot(consensus_data, ordered_genomes, biotypes, coverage_tgt, identity_tgt):
    """plots for consensus coverage, identity, score"""
    cons_iter = zip(*[['Coverage', 'Identity'],
                      [coverage_tgt, identity_tgt]])
    for mode, tgt in cons_iter:
        df = json_to_df_with_biotype(consensus_data, mode)
        df[mode] = df[mode]
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome', xlim=(-1, 101))


def consensus_support_plot(consensus_data, ordered_genomes, biotypes, modes, title, tgt):
    """grouped violin plots of original intron / intron annotation / exon annotation support"""
    def adjust_plot(g, this_title):
        g.set_xticklabels(rotation=60)
        g.fig.suptitle(this_title)
        g.fig.subplots_adjust(top=0.9)
        for ax in g.axes:
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
        g = sns.factorplot(data=df, y='value', x='genome', col='variable', col_wrap=2, kind='violin', sharex=True,
                           sharey=True, row_order=ordered_genomes, cut=0)
        adjust_plot(g, title)
        multipage_close(pdf, tight_layout=False)
        title += ' for {}'
        for biotype in biotypes:
            this_title = title.format(biotype)
            biotype_df = biotype_filter(df, biotype)
            g = sns.factorplot(data=biotype_df, y='value', x='genome', col='variable', col_wrap=2, kind='violin',
                               sharex=True, sharey=True, row_order=ordered_genomes, cut=0)
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
        base_title = '{} missing/failure rate in consensus gene set'.format(mode)
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
        df = pd.merge(fail_df, missing_df, on=['genome', 'biotype'])
        if len(df) == 0:
            continue
        df.genome = pd.Categorical(df.genome, ordered_genomes, ordered=True)
        with tgt.open('w') as outf, PdfPages(outf) as pdf:
            tot_df = df.groupby(by=['genome']).aggregate(np.sum).transpose()
            generic_unstacked_barplot(tot_df, pdf, base_title, tot_df.index, ylabel, tot_df.columns, 'category',
                                      bbox_to_anchor=(1.15, 0.7))
            cols = ['genome', '{} Failed'.format(mode), '{} Missing'.format(mode)]
            for biotype in biotypes:
                biotype_df = biotype_filter(df, biotype)
                if biotype_df is None:
                    continue
                biotype_df = biotype_df[cols]
                # re-add missing data
                missing = pd.DataFrame([[genome, 0, 0] for genome in biotype_df.genome.cat.categories if genome
                                        not in list(biotype_df.genome)])
                if len(missing) > 0:
                    missing.columns = [cols]
                    biotype_df = biotype_df.append(missing).sort_values('genome')
                title = base_title + ' for biotype {}'.format(biotype)
                biotype_df = biotype_df.set_index('genome').transpose()
                generic_unstacked_barplot(biotype_df, pdf, title, biotype_df.index, ylabel, biotype_df.columns,
                                          'category', bbox_to_anchor=(1.15, 0.7))


def tx_modes_plot(consensus_data, ordered_genomes, tx_mode_plot_tgt):
    modes_df = json_biotype_counter_to_df(consensus_data, 'Transcript Modes')
    df = modes_df.pivot(index='genome', columns='Transcript Modes').transpose()
    df.index = ['+'.join([y.replace('aug', '') for y in x.split(',')]) for x in df.index.get_level_values(1)]
    row_sum = df.sum(axis=1)
    total_sum = row_sum.sum()
    other = df.loc[row_sum / total_sum < 0.01].sum()
    df = df.loc[row_sum / total_sum >= 0.01]
    df.ix['Other'] = other
    df = df.sort_index(ascending=False)
    with tx_mode_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        title_string = 'Transcript modes in protein coding consensus gene set\n' \
                       '(combinations representing less than 1% of total are grouped in "Other")'
        ylabel = 'Number of transcripts'
        box_label = 'Transcript mode'
        generic_stacked_barplot(df, pdf, title_string, df.index, ylabel, ordered_genomes, box_label,
                                bbox_to_anchor=(1.15, 0.7))


def denovo_plot(consensus_data, ordered_genomes, denovo_tgt):
    with denovo_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_nested_counter_to_df(consensus_data, 'denovo')
        # fix column names because json_biotype_nested_counter_to_df makes assumptions
        df.columns = ['Result', 'Number of transcripts', 'Augustus mode', 'genome']
        has_pb = len(set(df['Augustus mode'])) == 2
        if has_pb is True:
            ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result',
                                hue='Augustus mode', col_wrap=2, row_order=ordered_genomes, sharex=True)
        else:
            ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result', col_wrap=2,
                                row_order=ordered_genomes, sharex=True)
        ax.set_xticklabels(rotation=60)
        ax.fig.suptitle('Incorporation of de-novo predictions')
        ax.fig.subplots_adjust(top=0.9)
        multipage_close(pdf, tight_layout=False)


def split_genes_plot(tm_data, ordered_genomes, split_plot_tgt):
    with split_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_counter_to_df(tm_data, 'Split Genes')
        df.columns = ['category', 'count', 'genome']
        title = 'Resolution methods of genes split across sequences'
        g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', col='category', xlabel='', col_wrap=2,
                            sharey=False, ylabel='Number of transcripts or genes', row_order=ordered_genomes,
                            title=title)


def pb_support_plot(consensus_data, ordered_genomes, pb_genomes, pb_support_tgt):
    with pb_support_tgt.open('w') as outf, PdfPages(outf) as pdf:
        pb_genomes = [x for x in ordered_genomes if x in pb_genomes]  # fix order
        df = json_biotype_counter_to_df(consensus_data, 'IsoSeq Transcript Valdiation')
        df.columns = ['IsoSeq Transcript Valdiation', 'Number of transcripts', 'genome']
        ax = sns.factorplot(data=df, x='genome', y='Number of transcripts', hue='IsoSeq Transcript Valdiation',
                            kind='bar', row_order=pb_genomes)
        ax.set_xticklabels(rotation=60)
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
    with improvement_tgt.open('w') as outf, PdfPages(outf) as pdf:
        for genome in ordered_genomes:
            data = pd.DataFrame(consensus_data[genome]['Evaluation Improvement'])
            data.columns = ['TransMapOriginalIntrons', 'TransMapIntronAnnotSupport', 'TransMapIntronRnaSupport',
                            'OriginalIntrons', 'IntronAnnotSupport', 'IntronRnaSupport']
            fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
            sns.regplot(x=data['TransMapOriginalIntrons'], y=data['OriginalIntrons'], ax=ax1)
            ax1.title('Original intron percent')
            sns.regplot(x=data['TransMapIntronAnnotSupport'], y=data['IntronAnnotSupport'], ax=ax2)
            ax2.title('Intron annotation support percent')
            sns.regplot(x=data['TransMapIntronRnaSupport'], y=data['IntronRnaSupport'], ax=ax3)
            ax3.title('Intron RNA support percent')
            fig.suptitle('AUGUSTUS metric improvement for {}'.format(genome))
            multipage_close(pdf, tight_layout=True)


###
# shared plotting functions
###


def cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=None, y=None, xlim=None, xlabel=None):
    """plots for coverage and identity. paired boxplot/violinplot format"""
    if xlabel is None:
        xlabel = 'Percent {}'.format(mode)
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        title = 'Overall {}'.format(mode)
        generate_boxplot_violin_pair(df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=xlim)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = '{} for biotype {}'.format(mode, biotype)
                generate_boxplot_violin_pair(biotype_df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=xlim)


###
# generic plotting functions
###


def generic_barplot(data, pdf, xlabel, ylabel, title, row_order=None, x=None, y=None, hue=None, hue_order=None,
                    col=None, col_wrap=None, sharex=True, sharey=True, col_order=None, palette=None, close=True):
    g = sns.factorplot(data=data, x=x, y=y, hue=hue, ci=None, kind='bar', hue_order=hue_order, row_order=row_order,
                       col=col, col_wrap=col_wrap, sharex=sharex, sharey=sharey, col_order=col_order, palette=palette)
    g.set_xticklabels(rotation=60)
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


def generate_boxplot_violin_pair(data, ordered_genomes, title, xlabel, pdf, hue=None, x=None, y=None, close=True,
                                 xlim=(0, 100)):
    """not so generic function that specifically produces a paired boxplot/violinplot"""
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=True)
    g1 = sns.boxplot(data=data, x=x, y=y, hue=hue, order=ordered_genomes, palette=choose_palette(ordered_genomes),
                     ax=ax1, saturation=boxplot_saturation, orient='h', width=bar_width, color=sns.color_palette()[0])
    g2 = sns.violinplot(data=data, x=x, y=y, hue=hue, order=ordered_genomes, palette=choose_palette(ordered_genomes),
                        ax=ax2, saturation=boxplot_saturation, orient='h', cut=0, scale='count',
                        color=sns.color_palette()[1])
    ax1.figure.suptitle(title)
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)
    sns.despine(trim=True, right=True, top=True, ax=ax1)
    sns.despine(trim=True, left=True, right=True, top=True, ax=ax2)
    ax2.set_ylabel('')
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    if close is True:
        multipage_close(pdf, tight_layout=False)
    return g1, g2


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
    ax.xaxis.set_ticklabels(names, rotation=60)


def sort_long_df(df, ordered_genomes):
    """sorts a long form dataframe by ordered genomes"""
    ordered_index = dict(zip(ordered_genomes, range(len(ordered_genomes))))
    df['order'] = df['genome'].map(ordered_index)
    df = df.sort_values('order')
    return df.drop('order', axis=1)
