"""
Generate all plots for the pipeline. For biotype specific plots, all plots are generated as a multi page PDF. There
is a plot for each biotype on its own, and one for the combined results.
"""
import json
import matplotlib
import itertools
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

    tm_splice_plot(consensus_data, args.ordered_genomes, biotypes, args.transmap_splice_support)
    tm_filter_plots(tm_data, args.ordered_genomes, args.transmap_filtering)
    tm_metrics_plot(tm_metrics, args.ordered_genomes, biotypes, transcript_biotype_map, args.tm_coverage,
                    args.tm_identity, args.tm_badness)
    tm_para_plot(para_data, args.ordered_genomes, biotypes, transcript_biotype_map, args.paralogy)
    consensus_metrics_plot(consensus_data, args.ordered_genomes, biotypes, args.coverage, args.identity, args.badness,
                           args.consensus_score)
    consensus_splice_plot(consensus_data, args.ordered_genomes, biotypes, args.consensus_splice_support)
    gene_fail_rescue_plot(consensus_data, args.ordered_genomes, biotypes, args.gene_failure)
    category_plot(consensus_data, args.ordered_genomes, biotypes, args.categories)
    completeness_plot(consensus_data, args.ordered_genomes, biotypes, args.completeness, gene_biotype_map,
                      transcript_biotype_map)
    if 'tx_modes' in args:
        tx_modes_plot(consensus_data, args.ordered_genomes, args.tx_modes)
    if 'novel' in args:
        novel_genes_plot(consensus_data, args.ordered_genomes, args.novel)
    if 'split_genes' in args:
        split_genes_plot(tm_data, args.ordered_genomes, args.split_genes)


###
# Load metrics from transMap PSLs
###


def load_tm_metrics(dbs):
    """Loads transMap data from PSLs"""
    tm_metrics = {'transMap Coverage': OrderedDict(), 'transMap Identity': OrderedDict(),
                  'transMap Badness': OrderedDict()}
    tm_name_map = {'TransMapCoverage': 'transMap Coverage', 'TransMapIdentity': 'transMap Identity',
                   'TransMapBadness': 'transMap Badness'}
    for genome, db_path in dbs.iteritems():
        session = tools.sqlInterface.start_session(db_path)
        table = tools.sqlInterface.TmEval
        for classifier in ['TransMapCoverage', 'TransMapIdentity', 'TransMapBadness']:
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


def tm_filter_plots(tm_data, ordered_genomes, tgt):
    """Plots for the transMap filtering process"""
    # massage the plots into something seaborn can handle. generate a background plot
    data = OrderedDict([[genome, tm_data[genome]['Paralogy']] for genome in ordered_genomes])
    df = pd.DataFrame.from_dict(data)
    df = df.transpose()
    df['combined_resolved'] = df['Transcripts resolved due to one non-zero'] + df['Transcripts resolved due to synteny delta']
    color_palette = sns.color_palette()
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        fig, ax = plt.subplots()
        xpos = np.arange(len(ordered_genomes))
        one_non_zero_bar = ax.bar(xpos, df['Transcripts resolved due to one non-zero'], bar_width,
                                  color=color_palette[1])[0]
        synteny_delta_bar = ax.bar(xpos, df['Transcripts resolved due to synteny delta'], bar_width,
                                   bottom=df['Transcripts resolved due to one non-zero'], color=color_palette[2])[0]
        # plot the single bar
        xpos = xpos + bar_width
        not_resolved_bar = ax.bar(xpos, df['Transcripts not resolved'], bar_width, color=color_palette[0])[0]
        set_ticks(ordered_genomes, ax)
        ax.set_xticks(xpos)
        ax.set_title('transMap paralogous alignment resolution by synteny score')
        legend_labels = ['Resolved due\nto one non-zero', 'Resolved due\nto synteny delta', 'Not resolved']
        fig.legend([one_non_zero_bar, synteny_delta_bar, not_resolved_bar], legend_labels,
                   bbox_to_anchor=(1.15, 0.6), frameon=True)
        sns.despine()
        multipage_close(pdf)
        ax = sns.barplot(data=pd.DataFrame([df['Alignments discarded']]))
        ax.set_title('Number of alignments discarded by synteny filtering')
        set_ticks(ordered_genomes, ax)
        ax.set_xticks(range(len(ordered_genomes)))
        sns.despine()
        multipage_close(pdf)


def tm_metrics_plot(tm_metrics, ordered_genomes, biotypes, transcript_biotype_map, tm_coverage_tgt, tm_identity_tgt,
                    tm_badness_tgt):
    """plots for transMap coverage, identity, badness"""
    tm_iter = zip(*[['transMap Coverage', 'transMap Identity', 'transMap Badness'],
                    [tm_coverage_tgt, tm_identity_tgt, tm_badness_tgt]])
    for mode, tgt in tm_iter:
        df = dict_to_df_with_biotype(tm_metrics[mode], transcript_biotype_map)
        for genome in ordered_genomes:
            df[genome] = 100 * df[genome]
        xlim = min(df.min(numeric_only=True, skipna=True)), max(df.max(numeric_only=True, skipna=True))
        cov_ident_badness_plot(biotypes, ordered_genomes, mode, tgt, df, xlim=xlim)


def consensus_metrics_plot(consensus_data, ordered_genomes, biotypes, coverage_tgt, identity_tgt, badness_tgt,
                           score_tgt):
    """plots for consensus coverage, identity, badness"""
    cons_iter = zip(*[['Coverage', 'Identity', 'Badness', 'Consensus Score'],
                      [coverage_tgt, identity_tgt, badness_tgt, score_tgt]])
    for mode, tgt in cons_iter:
        df = json_to_df_with_biotype(consensus_data, mode)
        if mode != 'Consensus Score':
            df[mode] = 100 * df[mode]
            xlabel = None
        else:
            xlabel = 'Consensus score'
        xlim = df[mode].min(skipna=True), df[mode].max(skipna=True)
        cov_ident_badness_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome', xlim=xlim, xlabel=xlabel)


def tm_para_plot(para_data, ordered_genomes, biotypes, transcript_biotype_map, para_tgt):
    """transMap paralogy plots"""
    def generate_hists(ordered_genomes, df):
        hists = {genome: np.roll(np.histogram(df[genome].fillna(0), paralogy_bins)[0], -1)
                 for genome in ordered_genomes}
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


def tm_splice_plot(consensus_data, ordered_genomes, biotypes, splice_tgt):
    """plots for splice junction suppport. paired boxplot/violinplot format"""
    def transform_biotype(d, biotype, genome):
        """transforms input data from JSON format to flat dataframe"""
        df = pd.DataFrame(d)
        df.columns = ['Percent of introns supported', 'Paralogy']
        df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
        df['Percent of introns supported'] = 100 * df['Percent of introns supported']
        return df

    def munge_legend():
        plt.gcf().axes[0].legend(bbox_to_anchor=(-0.25, 1), frameon=True, title='Paralogous\nalignment?')
        plt.gcf().axes[1].legend_.remove()

    dfs = []
    for genome in ordered_genomes:
        for biotype, d in consensus_data[genome]['transMap Splice Support'].iteritems():
            dfs.append(transform_biotype(d, biotype, genome))
    df = pd.concat(dfs)

    with splice_tgt.open('w') as outf, PdfPages(outf) as pdf:
        title = 'transMap splice junction support'
        xlabel = 'Percent splice junctions supported per transcript'
        generate_boxplot_violin_pair(df, ordered_genomes, title, xlabel, hue='Paralogy', y='genome',
                                     x='Percent of introns supported')
        munge_legend()
        multipage_close(pdf)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = 'transMap splice junction support for biotype {}'.format(biotype)
                generate_boxplot_violin_pair(biotype_df, ordered_genomes, title, xlabel, hue='Paralogy',
                                             y='genome', x='Percent of introns supported')
                munge_legend()
                multipage_close(pdf)


def consensus_splice_plot(consensus_data, ordered_genomes, biotypes, consensus_splice_tgt):
    """plots for splice junctions in consensus."""
    def transform_biotype(d, biotype, genome):
        """transforms input data from JSON format to flat dataframe"""
        df = pd.DataFrame(d)
        df.columns = ['Percent of introns supported']
        df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
        df['Percent of introns supported'] = 100 * df['Percent of introns supported']
        return df

    dfs = []
    for genome in ordered_genomes:
        for biotype, d in consensus_data[genome]['Splice Support'].iteritems():
            dfs.append(transform_biotype(d, biotype, genome))
    df = pd.concat(dfs)

    with consensus_splice_tgt.open('w') as outf, PdfPages(outf) as pdf:
        title = 'Consensus splice junction support'
        xlabel = 'Percent splice junctions supported per transcript'
        generate_boxplot_violin_pair(df, ordered_genomes, title, xlabel, y='genome', x='Percent of introns supported')
        multipage_close(pdf)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = 'Consensus splice junction support for biotype {}'.format(biotype)
                generate_boxplot_violin_pair(biotype_df, ordered_genomes, title, xlabel,
                                             y='genome', x='Percent of introns supported')
                multipage_close(pdf)


def category_plot(consensus_data, ordered_genomes, biotypes, category_plot_tgt):
    df = json_biotype_nested_counter_to_df(consensus_data, 'Transcript Categories')
    title = 'Consensus transcript categories'
    with category_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', hue='Transcript Categories', xlabel='',
                            ylabel='Number of transcripts', hue_order=['Excellent', 'Passing', 'Fail'],
                            row_order=ordered_genomes, title=title)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = 'Consensus transcript categories for biotype {}'.format(biotype)
                g = generic_barplot(pdf=pdf, data=biotype_df, x='genome', y='count', hue='Transcript Categories',
                                    xlabel='', ylabel='Number of transcripts',
                                    hue_order=['Excellent', 'Passing', 'Fail'],
                                    row_order=ordered_genomes, title=title)


def gene_fail_rescue_plot(consensus_data, ordered_genomes, biotypes, gene_fail_plot_tgt):
    fail_df = json_biotype_counter_to_df(consensus_data, 'Gene Failed')
    fail_df.columns = ['biotype', 'count', 'genome']
    fail_df['Category'] = ['Gene Failure'] * len(fail_df)
    rescue_df = json_biotype_counter_to_df(consensus_data, 'Gene Rescue')
    rescue_df.columns = ['biotype', 'count', 'genome']
    rescue_df['Category'] = ['Gene Rescue'] * len(rescue_df)
    df = pd.concat([fail_df, rescue_df])
    title = 'Rate of gene failure and gene rescue'
    with gene_fail_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        tot_df = df.groupby(by=['genome', 'Category']).aggregate(np.sum).reset_index()
        sort_long_df(tot_df, ordered_genomes)
        g = generic_barplot(pdf=pdf, data=tot_df, x='genome', y='count', hue='Category', xlabel='',
                            ylabel='Number of genes', hue_order=None, title=title)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                sort_long_df(tot_df, ordered_genomes)
                title = 'Rate of gene failure and gene rescue for biotype {}'.format(biotype)
                g = generic_barplot(pdf=pdf, data=biotype_df, x='genome', y='count', hue='Category', xlabel='',
                                    ylabel='Number of genes', hue_order=None, title=title)


def tx_modes_plot(consensus_data, ordered_genomes, tx_mode_plot_tgt):
    modes_df = json_biotype_counter_to_df(consensus_data, 'Transcript Modes')
    df = modes_df.pivot(index='genome', columns='Transcript Modes').transpose()
    with tx_mode_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        title_string = 'Transcript modes in consensus gene set'
        legend_labels = df.index.get_level_values('Transcript Modes')
        ylabel = 'Number of transcripts'
        box_label = 'Transcript\nmode'
        generic_stacked_barplot(df, pdf, title_string, legend_labels, ylabel, ordered_genomes, box_label)


def novel_genes_plot(consensus_data, ordered_genomes, novel_plot_tgt):
    with novel_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_counter_to_df(consensus_data, 'CGP')
        title = 'Novel genes/transcripts predicted by Augustus CGP'
        g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', hue='CGP', xlabel='',
                            ylabel='Number of transcripts/genes', row_order=ordered_genomes, title=title)
        df = json_flat_to_df(consensus_data, 'Novel isoforms')
        df.columns = ['genome', 'count']
        title = 'Transcripts predicted by Augustus CGP and assigned to comparative genes as novel isoforms'
        g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', hue=None, xlabel='',
                            ylabel='Number of transcripts',  row_order=ordered_genomes, title=title)


def split_genes_plot(tm_data, ordered_genomes, split_plot_tgt):
    with split_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        df = json_biotype_counter_to_df(tm_data, 'Split Genes')
        df.columns = ['category', 'count', 'genome']
        title = 'Resolution methods of genes split across sequences'
        g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', col='category', xlabel='', col_wrap=2,
                            sharey=False, ylabel='Number of transcripts or genes', row_order=ordered_genomes,
                            title=title)


def completeness_plot(consensus_data, ordered_genomes, biotypes, completeness_plot_tgt, gene_biotype_map,
                      transcript_biotype_map):
    def adjust_plot(g, gene_count, tx_count):
        for ax, c in zip(*[g.axes[0], [gene_count, tx_count]]):
            _ = ax.set_ylim(0, c)

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


###
# shared plotting functions
###


def cov_ident_badness_plot(biotypes, ordered_genomes, mode, tgt, df, x=None, y=None, xlim=None, xlabel=None):
    """plots for coverage, identity and badness. paired boxplot/violinplot format"""
    if xlabel is None:
        xlabel = 'Percent {}'.format(mode) if 'Badness' not in mode else 'Badness metric'
    with tgt.open('w') as outf, PdfPages(outf) as pdf:
        title = 'Overall {}'.format(mode)
        generate_boxplot_violin_pair(df, ordered_genomes, title, xlabel, x=x, y=y, xlim=xlim)
        multipage_close(pdf)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = '{} for biotype {}'.format(mode, biotype)
                generate_boxplot_violin_pair(biotype_df, ordered_genomes, title, xlabel, x=x, y=y, xlim=xlim)


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


def generate_boxplot_violin_pair(data, ordered_genomes, title, xlabel, hue=None, x=None, y=None, xlim=(0, 100)):
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


def _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label):
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=(1.12, 0.7), frameon=True,
               title=box_label)
    ax.set_title(title_string)
    ax.set_ylabel(ylabel)
    set_ticks(names, ax)
    ax.xaxis.set_ticks(np.arange(0, len(names)) + bar_width / 2.0)
    sns.despine(top=True, right=True)
    multipage_close(pdf)


def generic_unstacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label):
    fig, ax = plt.subplots()
    bars = []
    shorter_bar_width = bar_width / len(df)
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)) + shorter_bar_width * i, d, shorter_bar_width,
                           color=sns.color_palette()[i], linewidth=0.0))
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label)


def generic_stacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label):
    fig, ax = plt.subplots()
    bars = []
    cumulative = np.zeros(len(df.columns))
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)), d, bar_width, bottom=cumulative,
                           color=sns.color_palette()[i], linewidth=0.0))
        cumulative += d
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label)


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
    """converts JSON entries with many transcripts, such as those for coverage/badness/identity"""
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
    ordered_index = dict(zip(ordered_genomes ,range(len(ordered_genomes))))
    df['order'] = df['genome'].map(ordered_index)
    df = df.sort_values('order')
    return df.drop('order', axis=1)
