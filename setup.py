from setuptools import setup

setup(
    name='cat',
    version='0.1',
    py_modules=['cat', 'tools'],
    install_requires=[
        'pyfasta>=0.5.2',
        'toil>=3.0',
        'luigi>=2.0',
        'seaborn>=0.7',
        'pandas>=0.18',
        'frozendict',
        'configobj',
        'sqlalchemy',
        'ete3'
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
