from setuptools import setup

setup(
    name='cat',
    version='2.0',
    packages=['cat', 'tools'],
    python_requires='>=3.7.0',
    install_requires=[
        'pyfasta>=0.5.2',
        'toil==5.0',
        'luigi>=2.5',
        'seaborn>=0.12.1',
        'pandas>=1.1.0',
        'frozendict',
        'configobj>=5.0',
        'sqlalchemy>=1.0',
        'ete3',
        'pysam>=0.10',
        'numpy>=1.10',
        'scipy>=0.18.1',
        'bx-python>=0.7.1',
        'parasail',
        'gffutils',
        'biopython'
    ],
    scripts=['programs/cat_to_ncbi_submit', 'programs/translate_gene_pred',
             'programs/validate_gff3', 'programs/cat_parse_ncbi_genbank',
             'programs/cat_parse_ncbi_refseq', 'programs/cat_parse_prokka_gff3'],
    author='Ian Fiddes',
    description='Comparative Annotation Toolkit',
    url='https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit',
    license='Apache 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Bioinformatics',
        'Topic :: Bioinformatics',
        'License :: Apache 2.0',
        'Programming Language :: Python :: 3.7'
    ],
    keywords='bioinformatics comparative genomics',
)
