# Basic feature finder
#
# Developed on Python 2.7.3
# By: Katie Cunningham
# Last Updated: April 11, 2013


from operator import attrgetter
from collections import defaultdict

# Class structure:
#
#                   ----------
#                   |Interval|
#                   ----------
#           ((chrom, left, right, name))
#                       |
#                  -----------
#                  |          |
#              --------   ---------
#              |Region|   |Feature|
#              --------   ---------
#                             |
#                       --------------
#                       |            |
#                     ------
#                     |Gene|
#                     ------


class Interval(object):

    def __init__(self, chrom, left, right, name):

        self.chrom = chrom
        self.left = left
        self.right = right
        self.name = name

    def __str__(self):
        return "{} {} @ {} ({},{})".format(type(self).__name__, self.name,
                                            self.chrom, self.left, self.right)

    def __repr__(self):
        return self.__str__()


class Region(Interval):

    def __init__(self, Chromosome, StartPosition, EndPosition, RegionName):

        # superclass fields
        super(Region, self).__init__(Chromosome, int(StartPosition),
                                     int(EndPosition), RegionName)

        # fields directly from file
        self.Chromosome = Chromosome
        self.StartPosition = StartPosition
        self.EndPosition = EndPosition
        self.RegionName = RegionName


class Feature(Interval):

    def __init__(self, chrom, left, right, name):

        # superclass fields
        super(Feature, self).__init__(chrom, left, right, name)


class Gene(Feature):

    def __init__(self, chrom, strand, txStart, txEnd,
                 exonCount, exonStarts, exonEnds,
                 geneSymbol, refseq):

        # superclass fields
        super(Gene, self).__init__(chrom, int(txStart), int(txEnd),
                                   geneSymbol.rstrip(','))

        # other fields from file
        # (maybe stick these in a dict?)
        self.strand = strand
        self.txStart = txStart
        self.txEnd = txEnd
        self.exonCount = exonCount
        self.exonStarts = exonStarts
        self.exonEnds = exonEnds
        self.geneSymbol = geneSymbol.rstrip(',')
        self.refseq = refseq.rstrip(',')

        # calculated fields
        if strand == '+':
            self.positive_strand = True
        else:
            self.positive_strand = False
        self.exons = zip([int(s) for s in exonStarts.rstrip(',').split(',')],
                         [int(e) for e in exonEnds.rstrip(',').split(',')])


def sort_intervals(intervals_list):
    """Sort first by chrom (a str), then by left (an int), then by right
    (an int)

    """

    intervals_list.sort(key=attrgetter('chrom', 'left', 'right'))

def create_gene_list(gene_fp):
    """
    Creates and returns a sorted list of Genes from the file with the given
    filename.

    """

    newfile = False
    # Try to read from gene_file. If we can't, assume it's something we can
    # use to open a file.
    if not hasattr(gene_fp, 'read'):
        newfile = True
        gene_fp = file(gene_fp, 'r')

    header = gene_fp.readline()
    if not header.startswith('#'):
        print 'ERROR: No header on gene file. Wrong file type?'
        print 'First line beginning with "#" expected'
        exit(1)

    header = header.lstrip('#')  # Remove any leading #
    header = header.strip()  # Remove surrounding white space
    # Create a list of the header entries
    header_entries = header.split('\t')
    # Get the string after the last "." ("hg19.refGene.strand" -> "strand")
    # We don't care about version info
    header_entries = [h.split('.')[-1] for h in header_entries]
    ###print header_entries

    # This list will hold all the Gene objects
    genes = []

    # This loop starts at the next line of the file after the header
    for index, line in enumerate(gene_fp):

        ###print index

        line = line.strip() # Remove surrounding white space
        # Create a list of the line entries
        line_entries = line.split('\t')

        # Get the item in the list of line entries that matches certain fields
        # Note that this technique allows the columns in the file to be in
        # any order as long as they have the expected names
        chrom = line_entries[header_entries.index('chrom')]
        strand = line_entries[header_entries.index('strand')]
        txStart = line_entries[header_entries.index('txStart')]
        txEnd = line_entries[header_entries.index('txEnd')]
        exonCount = line_entries[header_entries.index('exonCount')]
        exonStarts = line_entries[header_entries.index('exonStarts')]
        exonEnds = line_entries[header_entries.index('exonEnds')]
        geneSymbol = line_entries[header_entries.index('geneSymbol')]
        refSeq = line_entries[header_entries.index('refseq')]

        # Make a Gene object
        gene = Gene(chrom, strand, txStart, txEnd,
                    exonCount, exonStarts, exonEnds,
                    geneSymbol, refSeq)

        # Add it to the end of the list
        genes.append(gene)

    # If we opened a new file, clean it up.
    if newfile:
        gene_fp.close()

    sort_intervals(genes)

    return genes


def create_region_list(region_fp):
    """
    Creates and returns a sorted list of Regions from the file with the given
    filename

    """
    newfile = False
    # Try to read from gene_file. If we can't, assume it's something we can
    # use to open a file.
    if not hasattr(region_fp, 'read'):
        newfile = True
        region_fp = open(region_fp, 'r')

    # Advance until line starting with "#" is read
    header = region_fp.readline()
    while not header.startswith('#'):
        header = region_fp.readline()
        # If we're at the end of the file, something is wrong
        if header == '':
            print 'ERROR: No header on region file. Wrong file type?'
            print 'Line beginning with "#" expected before data'
            exit(1)

    header = header.lstrip('#')  # Remove any leading #
    header = header.strip()  # Remove surrounding white space
    # Create a list of the header entries
    header_entries = header.split('\t')

    #This list will hold all the Region objects
    regions = []

    for index, line in enumerate(region_fp):

        line = line.strip()
        line_entries = line.split('\t')

        Chromosome = line_entries[header_entries.index('Chromosome')]
        StartPosition = line_entries[header_entries.index('StartPosition')]
        EndPosition = line_entries[header_entries.index('EndPosition')]
        RegionName = line_entries[header_entries.index('RegionName')]

        region = Region(Chromosome, StartPosition, EndPosition, RegionName)

        regions.append(region)

    if newfile:
        region_fp.close()

    sort_intervals(regions)

    return regions

def print_comp(feature, region):

    pass
    #print 'Feat: ', feature.chrom, '(', feature.left, ',', feature.right, ')'
    #print 'Reg:  ', region.chrom, '(', region.left, ',', region.right, ')'


def overlaps(feature, region):
    """Returns True if the feature overlaps the region"""

    ###print 'In overlaps...'
    print_comp(feature, region)

    return (region.chrom == feature.chrom and
            min(region.right, feature.right)
            - max(region.left, feature.left) > 0)


def before(feature, region):
    """Returns True if the feature is completely before the region

    If the end of the feature is before the beginning of the region,
    returns True"""

    ###print 'In before...'
    ###print_comp(feature, region)

    if feature.chrom < region.chrom:
        return True
    elif feature.chrom > region.chrom:
        return False
    else: # They are on the same chromosome
        return feature.right < region.left


def after(feature, region):
    """Returns True if the feature is completely after the region

    If the beginning of the feature is after the end of the region,
    returns True"""

    ###print 'In after...'
    ###print_comp(feature, region)

    if feature.chrom > region.chrom:
        return True
    elif feature.chrom < region.chrom:
        return False
    else: # They are on the same chromosome
        return feature.left > region.right


def get_overlap(region, feature):
    """ Returns the number of bases that the Region and the Feature overlap """

    if region.chrom != feature.chrom:
        return 0

    return max(0, min(region.right, feature.right)
                     - max(region.left, feature.left) + 1)


def get_amount_before(feature, region):
    """ Returns the number of bases the the Feature is before the region """

    return -1


def get_amount_after(feature, region):
    """ Returns the number of bases the the Feature is after the region """

    return -1


def print_overlap_info(feature, region):
    """
    Prints an intelligable message about how the Feature and the Region overlap

    """

    print "{} overlaps {} by {} base pairs".format(region.name,
                                                   feature.name,
                                                   get_overlap(region,
                                                               feature))

def find_features(region_list, gene_list=None):
    """Returns a list of the features found"""

    found = {region : [] for region in region_list}

    if gene_list is None:

        #print "No gene list given"
        return

    gene_index = 0

    # For each region
    for region_index in xrange(len(region_list)):

        cur_region = region_list[region_index]
        print 'cur_region:', cur_region.name
        cur_gene = gene_list[gene_index]
        print 'cur_gene:', cur_gene.name

        # Advance features until feature overlaps or is after current region
        while ((gene_index < len(gene_list)) and
               before(cur_gene, cur_region)):
            print 'skipping:', cur_gene.name
            gene_index += 1
            if gene_index >= len(gene_list):
                break
            else:
                cur_gene = gene_list[gene_index]

        # Advance features and add them to that regions' entry in the dict
        while overlaps(cur_gene, cur_region):

            print 'adding:', cur_gene.name
            print_overlap_info(cur_region, cur_gene)
            found[cur_region].append(cur_gene)

            if cur_gene.right > cur_region.right:  # Don't move ahead when gene
                break                              # may overlap next region
            else:
                gene_index += 1
                if gene_index >= len(gene_list):
                    gene_index -= 1
                    break
                else:
                    cur_gene = gene_list[gene_index]


    return found

"""
GENE_FILENAME = 'gene_file'
REGION_FILENAME = 'region_file'


genes = create_gene_list(GENE_FILENAME)
regions = create_region_list(REGION_FILENAME)

genes = genes[:100]
regions = regions[:100]
"""



