"""
Object representation of GFF3 data and parser for GFF3 files.

See: http://www.sequenceontology.org/gff3.shtml
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pandas as pd
import urllib
import copy
import re
import collections
import fileOps
from tools import PycbioException

GFF3_HEADER = "##gff-version 3"


class GFF3Exception(PycbioException):
    """
    Exception associated with GFF3 data.
    """
    def __init__(self, message, fileName=None, lineNumber=None, cause=None):
        if fileName is not None:
            message = "{}:{}: {}".format(fileName, lineNumber, message)
        super(GFF3Exception, self).__init__(message, cause)

# characters forcing encode for columns
_encodeColRegexpStr = "\t|\n|\r|%|[\x00-\x1F]|\x7f"
_encodeColRegexp = re.compile(_encodeColRegexpStr)

# characters forcing encode for attribute names or values
_encodeAttrReStr = _encodeColRegexpStr + "|;|=|&|,"
_encodeAttrRe = re.compile(_encodeAttrReStr)


def _encodeCol(v):
    """
    Encode a column if is has special characters.
    """
    if _encodeColRegexp.search(v):
        return urllib.quote(v)
    else:
        return v


def _encodeAttr(v):
    """
    Encode a attribute name or value if is has special characters.
    """
    if _encodeAttrRe.search(v):
        return urllib.quote(v)
    else:
        return v


class Feature(object):
    """
    One feature notation from a GFF3, missing attribute, as code by `.' in
    the file are stored as None.
    """

    def __init__(self, seqname, source, type, start, end, score, strand,
                 frame, attributes, gff3Set, lineNumber=None, copy_attrs=False):
        """
        The attributes field is a dict of lists/tupes as attributes maybe
        multi-valued.
        """
        self.seqname = seqname
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        if copy_attrs is True:
            self.attributes = copy.deepcopy(attributes)
        else:
            self.attributes = attributes
        self.gff3Set = gff3Set
        self.lineNumber = lineNumber
        self.parents = set()
        self.children = set()

    @staticmethod
    def __dotIfNone(val):
        return val if val is not None else '.'

    def __attributeStr(self, name):
        """
        Return name=value for a single attribute
        """
        return "{}={}".format(
            _encodeAttr(name),
            ",".join([_encodeAttr(v) for v in self.attributes[name]]))

    def __attributeStrs(self):
        """
        Return name=value, semi-colon-separated string for attributes, including
        url-style quoting
        """
        return ";".join([self.__attributeStr(name)
                         for name in self.attributes.iterkeys()])

    def __str__(self):
        """
        Return the object as a valid GFF3 record line.
        """
        return "\t".join([self.seqname, self.source, self.type,
                          str(self.start), str(self.end), self.__dotIfNone(self.score),
                          self.__dotIfNone(self.strand), self.__dotIfNone(self.frame),
                          self.__attributeStrs()])

    @property
    def featureId(self):
        """
        ID attribute or None if records doesn't have an ID
        """
        featId = self.attributes.get("ID")
        if featId is not None:
            featId = featId[0]
        return featId

    def getRequiredAttribute(self, name):
        "get a require attributes list of values, error if not found"
        values = self.attributes.get(name)
        if values is None:
            raise GFF3Exception("required attribute not found: {}".format(name),
                                self.gff3Set.fileName, self.lineNumber)
        return values


class Gff3Set(object):
    """
    A set of GFF3 sequence annotations
    """
    def __init__(self, fileName=None):
        self.fileName = fileName
        self.roots = set()     # root nodes (those with out parents)
        # index of features by id. GFF3 allows disjoint features with
        # the same id.  None is used to store features without ids
        self.byFeatureId = collections.defaultdict(list)

    def add(self, feature):
        """
        Add a feature record
        """
        # None featureId allowed
        self.byFeatureId[feature.featureId].append(feature)

    def __linkFeature(self, feature):
        """
        Link a feature with it's parents.
        """
        parentIds = feature.attributes.get("Parent")
        if parentIds is None:
            self.roots.add(feature)
        else:
            for parentId in parentIds:
                self.__linkToParent(feature, parentId)

    def __linkToParent(self, feature, parentId):
        """
        Link a feature with it's children
        """
        parentParts = self.byFeatureId.get(parentId)
        if parentParts is None:
            raise GFF3Exception("Parent feature does not exist: {}".format(parentId),
                                self.fileName, feature.lineNumber)
        # parent maybe disjoint
        for parentPart in parentParts:
            feature.parents.add(parentPart)
            parentPart.children.add(feature)

    def finish(self):
        """
        finish loading the set, constructing the tree
        """
        # features maybe disjoint
        for featureParts in self.byFeatureId.itervalues():
            for feature in featureParts:
                self.__linkFeature(feature)

    @staticmethod
    def __recSortKey(r):
        return (r.seqname, r.start, -r.end, r.type)

    def __writeRec(self, fh, rec):
        fh.write(str(rec) + "\n")
        for child in sorted(rec.children, key=self.__recSortKey):
            self.__writeRec(fh, child)

    def write(self, fh):
        """
        Write set to a GFF3 format file.
        """
        fh.write(GFF3_HEADER + "\n")
        for root in sorted(self.roots, key=self.__recSortKey):
            self.__writeRec(fh, root)


class Gff3Parser(object):
    """
    Parser for GFF3 files.  This parse does basic validation, but does not
    fully test for conformance.
    """

    def __init__(self, fileName, copy_attrs=False):
        """
        If fileName ends with .gz or .bz2, it will decompressed.  If fh is
        specified, then parse the already opened stream, with file name still
        used for error message, otherwise the file be opened and parsed.
        """
        self.fileName = fileName
        self.lineNumber = 0
        self.copy_attrs = copy_attrs

    def __open(self):
        """
        open input file, optionally with decompression
        """
        return fileOps.opengz(self.fileName)

    # parses `attr=val'; GFF3 spec is not very specific on the allowed values.
    SPLIT_ATTR_RE = re.compile("^([a-zA-Z][^=]*)=(.*)$")

    def __parseAttrVal(self, attrStr):
        """
        Returns tuple of tuple of (attr, value), multiple are returned to
        handle multi-value attributes.
        """
        m = self.SPLIT_ATTR_RE.match(attrStr)
        if m is None:
            raise GFF3Exception("can't parse attribute/value: '" + attrStr +
                                "'", self.fileName, self.lineNumber)
        name = urllib.unquote(m.group(1))
        val = m.group(2)
        # Split by comma separate then unquote.  Commas in values must be
        # url encoded.
        return name, [urllib.unquote(v) for v in val.split(',')]

    SPLIT_ATTR_COL_RE = re.compile("; *")

    def __parseAttrs(self, attrsStr):
        """
        Parse the attributes and values
        """
        attributes = dict()
        for attrStr in self.SPLIT_ATTR_COL_RE.split(attrsStr):
            name, vals = self.__parseAttrVal(attrStr)
            if name in attributes:
                raise GFF3Exception("duplicated attribute name: {}".format(name),
                                    self.fileName, self.lineNumber)
            attributes[name] = vals
        return attributes

    GFF3_NUM_COLS = 9

    def __parseRecord(self, gff3Set, line):
        """
        Parse one record.
        """
        row = line.split("\t")
        if len(row) != self.GFF3_NUM_COLS:
            raise GFF3Exception(
                "Wrong number of columns, expected {}, got {}".format(
                    self.GFF3_NUM_COLS, len(row)),
                self.fileName, self.lineNumber)
        feature = Feature(urllib.unquote(row[0]), urllib.unquote(row[1]),
                          urllib.unquote(row[2]),
                          int(row[3]), int(row[4]), row[5], row[6], row[7],
                          self.__parseAttrs(row[8]), gff3Set, self.lineNumber, self.copy_attrs)
        gff3Set.add(feature)

    # spaces or comment line
    IGNORED_LINE_RE = re.compile("(^[ ]*$)|(^[ ]*#.*$)")

    @staticmethod
    def __isIgnoredLine(line):
        return Gff3Parser.IGNORED_LINE_RE.search(line) is not None

    def __checkHeader(self, line):
        # split to allow multiple spaces and tabs
        if line.split() != GFF3_HEADER.split():
            raise GFF3Exception(
                "First line is not GFF3 header ({}), got: {}".format(
                    GFF3_HEADER, line), self.fileName, self.lineNumber)

    def __parseLine(self, gff3Set, line):
        if self.lineNumber == 1:
            self.__checkHeader(line)
        elif not self.__isIgnoredLine(line):
            self.__parseRecord(gff3Set, line)

    def parse(self):
        """
        Parse and return a Gff3Set object in format.
        """
        fh = self.__open()
        try:
            gff3Set = Gff3Set(self.fileName)
            for line in fh:
                self.lineNumber += 1
                self.__parseLine(gff3Set, line[0:-1])
        finally:
            fh.close()
        gff3Set.finish()
        return gff3Set


def extract_attrs(gff3):
    """
    Extracts attributes table from gff3.
    The attributes we care about are gene_id, gene_name, gene_type, transcript_id, transcript_type, start_codon,
    stop_codon, five_prime_UTR, three_prime_UTR
    :param gff3: input gff3
    :returns: DataFrame
    """
    results = {}
    parser = Gff3Parser(gff3)
    tree = parser.parse()
    for gene in tree.roots:
        if 'gene_id' not in gene.attributes:  # this isn't a gene
            continue
        try:
            gene_id = gene.attributes['gene_id'][0]
            try:
                gene_biotype = gene.attributes['biotype'][0]
            except KeyError:  # attempt Gencode naming scheme
                gene_biotype = gene.attributes['gene_type'][0]
            try:
                gene_name = gene.attributes['Name'][0]
            except KeyError:  # attempt Gencode naming scheme
                try:
                    gene_name = gene.attributes['gene_name'][0]
                except KeyError:
                    # this record lacks a common name field. Probably from Ensembl
                    gene_name = gene.attributes['gene_id'][0]
            for tx in gene.children:
                if tx.type != 'transcript':
                    continue
                try:
                    tx_id = tx.attributes['transcript_id'][0]
                    try:
                        name = tx.attributes['transcript_name'][0]
                    except KeyError:
                        name = tx.attributes['Name'][0]
                    try:
                        tx_biotype = tx.attributes['biotype'][0]
                    except KeyError:  # attempt Gencode naming scheme
                        tx_biotype = tx.attributes['transcript_type'][0]
                    feature_types = {c.type for c in tx.children}
                    r = {'TranscriptBiotype': tx_biotype, 'GeneId': gene_id,
                         'GeneName': gene_name, 'GeneBiotype': gene_biotype,
                         'StartCodon': 'start_codon' in feature_types,
                         'StopCodon': 'stop_codon' in feature_types,
                         'TranscriptName': name}
                    results[tx_id] = r
                except KeyError, e:
                    raise GFF3Exception('Unable to parse field {} from the input gff3 on line {}'.format(e,
                                                                                                         tx.lineNumber))
        except KeyError, e:
            raise GFF3Exception('Unable to parse field {} from the input gff3 on line {}'.format(e, gene.lineNumber))
    df = pd.DataFrame.from_dict(results, orient='index')
    df.StartCodon = pd.to_numeric(df.StartCodon)  # force the types here for better sql dtypes
    df.StopCodon = pd.to_numeric(df.StopCodon)
    df.index.rename('TranscriptId', inplace=True)
    # make into a nice order for the sqlite database
    df = df[['TranscriptName', 'TranscriptBiotype', 'GeneId', 'GeneName', 'GeneBiotype', 'StartCodon', 'StopCodon']]
    return df
