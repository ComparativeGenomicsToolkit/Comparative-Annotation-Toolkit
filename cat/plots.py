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
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.annotation_db)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.annotation_db)
    biotypes = sorted(tools.sqlInterface.get_transcript_biotypes(args.annotation_db))
    args.ordered_genomes = list(args.ordered_genomes)  # weird bug in pandas

    # hack to bring coding to the top
    try:
        biotypes.insert(0, biotypes.pop(biotypes.index('protein_coding')))
    except ValueError:
        pass

    tx_modes_plot(consensus_data, args.ordered_genomes, args.tx_modes)
    tm_metrics_plot(tm_metrics, args.ordered_genomes, biotypes, transcript_biotype_map, args.tm_coverage,
                    args.tm_identity)
    tm_para_plot(tm_data, args.ordered_genomes, biotypes, args.paralogy)
    tm_gene_family_plot(tm_data, args.ordered_genomes, biotypes, args.gene_collapse)
    consensus_metrics_plot(consensus_data, args.ordered_genomes, biotypes, args.coverage, args.identity)
    missing_rate_plot(consensus_data, args.ordered_genomes, biotypes, args.missing)
    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Annotation Support', 'Exon Annotation Support', 'Original Introns'],
                           title='Reference annotation support',
                           tgt=args.consensus_annot_support)
    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Support', 'Exon Support'],
                           title='Extrinsic support',
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


###
# Plots
###


def tm_metrics_plot(tm_metrics, ordered_genomes, biotypes, transcript_biotype_map, tm_coverage_tgt, tm_identity_tgt):
    """plots for transMap coverage, identity"""
    tm_iter = zip(*[['transMap Coverage', 'transMap Identity'],
                    [tm_coverage_tgt, tm_identity_tgt]])
    for mode, tgt in tm_iter:
        df = dict_to_df_with_biotype(tm_metrics[mode], transcript_biotype_map)
        df = pd.melt(df, id_vars='biotype', value_vars=ordered_genomes).dropna()
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


def tm_para_plot(tm_data, ordered_genomes, biotypes, para_tgt):
    """transMap paralogy plots"""
    legend_labels = ['= 1', '= 2', '= 3', u'\u2265 4']
    title_string = 'Proportion of transcripts that have multiple alignments'
    biotype_title_string = 'Proportion of {} transcripts that have multiple alignments'
    df = json_biotype_nested_counter_to_df(tm_data, 'Paralogy')
    # we want a dataframe where each row is the counts, in genome order
    # we construct the transpose first
    r = []
    df['Paralogy'] = pd.to_numeric(df['Paralogy'])
    # make sure genomes are in order
    df['genome'] = pd.Categorical(df['genome'], ordered_genomes, ordered=True)
    df = df.sort_values('genome')
    for biotype, biotype_df in df.groupby('biotype'):
        for genome, genome_df in biotype_df.groupby('genome'):
            high_para = genome_df[genome_df.Paralogy >= 4]['count'].sum()
            counts = dict(zip(genome_df['Paralogy'], genome_df['count']))
            r.append([biotype, genome, counts.get(1, 0), counts.get(2, 0), counts.get(3, 0), high_para])
    df = pd.DataFrame(r, columns=['biotype', 'genome', '1', '2', '3', u'\u2265 4'])
    sum_df = df.groupby('genome', sort=False).aggregate(sum).T

    plot_fn = generic_unstacked_barplot if len(df.columns) <= 5 else generic_stacked_barplot
    box_label = 'Number of\nalignments'
    with para_tgt.open('w') as outf, PdfPages(outf) as pdf:
        plot_fn(sum_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes, box_label)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                biotype_df = biotype_df.drop(['genome', 'biotype'], axis=1).T
                title_string = biotype_title_string.format(biotype)
                plot_fn(biotype_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes,
                        box_label)


def tm_gene_family_plot(tm_data, ordered_genomes, biotypes, gene_family_tgt):
    """transMap gene family collapse plots."""
    try:
        df = json_biotype_nested_counter_to_df(tm_data, 'Gene Family Collapse')
    except ValueError:  # no gene family collapse. probably the test set.
        with gene_family_tgt.open('w') as outf:
            pass
        return
    df['Gene Family Collapse'] = pd.to_numeric(df['Gene Family Collapse'])
    tot_df = df[['Gene Family Collapse', 'genome', 'count']].\
        groupby(['genome', 'Gene Family Collapse']).aggregate(sum).reset_index()
    tot_df = tot_df.sort_values('Gene Family Collapse')
    with gene_family_tgt.open('w') as outf, PdfPages(outf) as pdf:
        g = sns.factorplot(y='count', col='genome', x='Gene Family Collapse', data=tot_df, kind='bar',
                           col_order=ordered_genomes, col_wrap=4)
        g.fig.suptitle('Number of genes collapsed during gene family collapse')
        g.set_xlabels('Number of genes collapsed to one locus')
        g.set_ylabels('Number of genes')
        g.fig.subplots_adjust(top=0.9)
        multipage_close(pdf, tight_layout=False)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is None:
                continue
            biotype_df = biotype_df.sort_values('Gene Family Collapse')
            g = sns.factorplot(y='count', col='genome', x='Gene Family Collapse', data=biotype_df, kind='bar',
                               col_order=[x for x in ordered_genomes if x in set(biotype_df.genome)], col_wrap=4)
            g.fig.suptitle('Number of genes collapsed during gene family collapse for {}'.format(biotype))
            g.set_xlabels('Number of genes collapsed to one locus')
            g.set_ylabels('Number of genes')
            g.fig.subplots_adjust(top=0.9)
            multipage_close(pdf, tight_layout=False)


def missing_rate_plot(consensus_data, ordered_genomes, biotypes, missing_plot_tgt):
    """Missing genes/transcripts"""
    base_title = 'Number of missing orthologs in consensus set'
    gene_missing_df = json_biotype_counter_to_df(consensus_data, 'Gene Missing')
    gene_missing_df.columns = ['biotype', 'Genes', 'genome']
    transcript_missing_df = json_biotype_counter_to_df(consensus_data, 'Transcript Missing')
    transcript_missing_df.columns = ['biotype', 'Transcripts', 'genome']
    df = transcript_missing_df.merge(gene_missing_df, on=['genome', 'biotype'])
    df = pd.melt(df, id_vars=['biotype', 'genome'])
    ylabel = 'Number of genes or transcripts'
    with missing_plot_tgt.open('w') as outf, PdfPages(outf) as pdf:
        tot_df = df.groupby(['genome', 'biotype', 'variable']).aggregate(sum).reset_index()
        generic_barplot(tot_df, pdf, '', ylabel, base_title, x='genome', y='value',
                        col='variable', row_order=ordered_genomes)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is None:
                continue
            biotype_df = biotype_df.groupby(['genome', 'variable']).aggregate(sum).reset_index()
            title = base_title + ' for biotype {}'.format(biotype)
            generic_barplot(biotype_df, pdf, '', ylabel, title, x='genome', y='value',
                            col='variable', row_order=ordered_genomes)


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
            df = df.fillna(0)
            generic_stacked_barplot(df, pdf, title_string, df.index, ylabel, ordered_genomes, 'Transcript mode(s)',
                                    bbox_to_anchor=(1.25, 0.7))

        else:
            generic_barplot(pd.melt(df, id_vars='Modes'), pdf, 'Transcript mode(s)', ylabel, title_string, x='Modes',
                            y='value', order=ordered_values)


def denovo_plot(consensus_data, ordered_genomes, denovo_tgt):
    with denovo_tgt.open('w') as outf, PdfPages(outf) as pdf:
        try:
            df = json_biotype_nested_counter_to_df(consensus_data, 'denovo')
        except ValueError:
            # No de novo results. Probably the test set.
            return
        # fix column names because json_biotype_nested_counter_to_df makes assumptions
        df.columns = ['Result', 'Number of transcripts', 'Augustus mode', 'genome']
        has_pb = len(set(df['Augustus mode'])) == 2
        if len(set(df.genome)) > 1:  # if we ran in PB only, we may not have multiple genomes
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
        if len(df) == 0:
            # no support information
            return
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
            sns.kdeplot(x, y, ax=ax, cut=0, cmap='Purples_d', shade=True, shade_lowest=False, n_levels=n_levels, bw=bw,
                        rasterized=True)
        except:
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
            for ax in [ax1, ax2, ax3, ax4]:  # goodness plots are allowed to auto-set scale
                ax.set_xlim(0, 100)
                ax.set_ylim(0, 100)
            
            do_kdeplot(data['transMap original introns'], data['Original introns'], ax1, n_levels=25, bw=2)
            sns.regplot(x=data['transMap original introns'], y=data['Original introns'], ax=ax1,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)
            do_kdeplot(data['transMap intron annotation support'], data['Intron annotation support'], ax2,
                       n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron annotation support'], y=data['Intron annotation support'], ax=ax2,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)          
            do_kdeplot(data['transMap intron RNA support'], data['Intron RNA support'], ax3, n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron RNA support'], y=data['Intron RNA support'], ax=ax3,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)
            
            do_kdeplot(data['transMap alignment goodness'], data['Alignment goodness'], ax4, n_levels=20, bw=1)
            sns.regplot(x=data['transMap alignment goodness'], y=data['Alignment goodness'], ax=ax4,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)
            

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
        try:  # this is a hack to deal with weird small input datasets
            tm_df['genome'] = ordered_genomes
        except:
            return
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
        if key in d:
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
