from setuptools import setup

setup(
    name='cat',
    version='0.1',
    py_modules=['cat', 'tools'],
    install_requires=[
        'pyfasta>=0.5.2',
        'toil>=3.5',
        'luigi>=2.5',
        'seaborn>=0.7',
        'pandas>=0.18',
        'frozendict',
        'configobj>=5.0',
        'sqlalchemy>=1.0',
        'ete3',
        'pysam>=0.10',
        'numpy>=1.10',
        'scipy>=0.18.1',
        'bx-python>=0.7.1',
        'bd2k-python-lib==1.14a1.dev48'
    ],
    author='Ian Fiddes',
    description='Comparative Annotation Toolkit',
    url='https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit',
    license='Apache 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Bioinformatics',
        'Topic :: Bioinformatics',
        'License :: Apache 2.0',
        'Programming Language :: Python :: 2.7'
    ],
    keywords='bioinformatics comparative genomics',
)
