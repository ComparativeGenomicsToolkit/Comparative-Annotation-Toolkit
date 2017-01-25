import os
import luigi
import tools.fileOps
import tools.misc
from frozendict import frozendict
from cat import RebuildableTask, Consensus, PipelineTask, GenomeFiles


class AssemblyHub(RebuildableTask, luigi.WrapperTask):
    """
    Construct an assembly hub out of all the results
    """
    def run(self):
        tools.fileOps.ensure_dir(self.out_dir)
        yield self.clone(CreateDirectoryStructure)
        yield self.clone(CreateTracks)
        yield self.clone(CreateTrackDb)


class CreateDirectoryStructure(PipelineTask):
    """
    Constructs the directory structure. Creates symlinks for all relevant files.
    """
    @staticmethod
    def get_args(pipeline_args):
        args = tools.misc.HashableNamespace()
        args.out_dir = os.path.join(pipeline_args.out_dir, 'assemblyHub')
        args.hub = os.path.join(args.out_dir, 'hub.txt')
        args.genomes = os.path.join(args.out_dir, 'genomes.txt')
        genome_files = frozendict({genome: GenomeFiles.get_args(pipeline_args, genome)
                                   for genome in pipeline_args.hal_genomes})
        sizes = {}
        twobits = {}
        for genome, genome_file in genome_files.iteritems():
            sizes[genome] = (genome_file.sizes, os.path.join(args.out_dir, genome, 'chrom.sizes'))
            twobits[genome] = (genome_file.two_bit, os.path.join(args.out_dir, genome, '{}.2bit'.format(genome)))
        args.sizes = frozendict(sizes)
        args.twobits = frozendict(twobits)
        return args

    def requires(self):
        yield self.clone(Consensus)

    def output(self):
        pipeline_args = self.get_pipeline_args()
        args = CreateDirectoryStructure.get_args(pipeline_args)
        yield luigi.LocalTarget(args.hub)
        yield luigi.LocalTarget(args.genomes)
        for local_path, hub_path in args.sizes.itervalues():
            yield luigi.LocalTarget(hub_path)
        for local_path, hub_path in args.twobits.itervalues():
            yield luigi.LocalTarget(hub_path)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        args = CreateDirectoryStructure.get_args(pipeline_args)
        hal = os.path.splitext(os.path.basename(pipeline_args.hal))[0]
        # write the hub.txt file
        with open(args.hub, 'w') as outf:
            outf.write(hub_str.format(hal))

        # write the genomes.txt file, construct a dir
        # TODO: include ancestors? Provide a flag?
        with open(args.genomes, 'w') as outf:
            for genome, (sizes_local_path, sizes_hub_path) in args.sizes.iteritems():
                outf.write(genome_str.format(genome, find_default_pos(sizes_local_path)))

        # construct a directory for each genome
        for genome, (sizes_local_path, sizes_hub_path) in args.sizes.iteritems():
            tools.fileOps.ensure_file_dir(sizes_hub_path)
            os.link(sizes_local_path, sizes_hub_path)
            twobit_local_path, twobit_hub_path = args.twobits[genome]
            os.link(twobit_local_path, twobit_hub_path)




def find_default_pos(chrom_sizes, window_size=30000):
    """
    Returns a window_size window over the beginning of the largest chromosome
    :param chrom_sizes: chrom sizes file
    :param window_size: window size to extend from
    :return: string
    """
    sizes = [x.split() for x in open(chrom_sizes)]
    sorted_sizes = sorted(sizes, key=lambda (chrom, size): -int(size))
    return '{}:{}-{}'.format(sorted_sizes[0][0], 1, window_size)


hub_str = '''hub {hal}
shortLabel {hal}
longLabel {hal}
genomesFile genomes.txt
email NoEmail

'''

genome_str = '''genome {genome}
twoBitPath {genome}/{genome}.2bit
trackDb {genome}/trackDb.txt
organism {genome}
description {genome}
scientificName {genome}
defaultPos {default_pos}

'''

