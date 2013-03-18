
import random
import unittest

from feature_finder import *

# NOTE these tests run in the order of the class names
# So, TestOverlap(...) runs before TestRegionList(...)

class TestGeneList(unittest.TestCase):

    def setUp(self):
        GENE_FILENAME = 'gene_file'
        self.genes = create_gene_list(GENE_FILENAME)


    def test_sorted(self):
        """ Make sure that the lists of all types are aorted correctly"""

        for gene in self.genes:

            # Check that the start value is always before the end value
            self.assertTrue(int(gene.txStart) < int(gene.txEnd))
            # Check that the start value is always the same as the first of
            # the exons
            self.assertEquals(int(gene.txStart), gene.exons[0][0])
            # Check that the end value is always the same as the last of
            # the exons
            self.assertEquals(int(gene.txEnd), gene.exons[-1][-1])
            # Check that the exonCount is always the same as the number of
            # exons listed
            self.assertEquals(int(gene.exonCount), len(gene.exons))


class TestRegionList(unittest.TestCase):

    def setUp(self):
        REGION_FILENAME = 'region_file'
        self.regions = create_region_list(REGION_FILENAME)


    def test_sorted(self):
        """ Make sure that the lists of all types are aorted correctly"""

        for region in self.regions:

            # Check that the start value is always before the end value
            self.assertTrue(int(region.StartPosition)
                            < int(region.EndPosition))
            # Check that the left value is always before the right value
            self.assertTrue(int(region.left) < int(region.right))


class TestBeforeAndAfter(unittest.TestCase):

    def setUp(self):
        pass

    def test_before(self):
        """Make sure the before function works"""

        # Same chr...

        f_0_500 = Feature('chr3', 0, 500, '')
        f_50_75 = Feature('chr3', 50, 75, '')
        f_90_150 = Feature('chr3', 90, 150, '')
        f_99_201 = Feature('chr3', 99, 201, '')
        f_100_170 = Feature('chr3', 100, 170, '')
        f_100_200 = Feature('chr3', 100, 200, '')
        f_150_170 = Feature('chr3', 150, 170, '')
        f_150_200 = Feature('chr3', 150, 200, '')
        f_170_230 = Feature('chr3', 170, 230, '')
        f_300_500 = Feature('chr3', 300, 500, '')

        r_100_200 = Region('chr3', 100, 200, '')

        # feature is before region
        self.assertTrue(before(f_50_75, r_100_200))

        # feature is after region
        self.assertFalse(before(f_300_500, r_100_200))

        # feature is exactly the same as region
        self.assertFalse(before(f_100_200, r_100_200))

        # feature overlaps region to the left
        self.assertFalse(before(f_90_150, r_100_200))

        # feature overlaps region to the right
        self.assertFalse(before(f_170_230, r_100_200))

        # feature completely overlaps region
        self.assertFalse(before(f_0_500, r_100_200))
        self.assertFalse(before(f_99_201, r_100_200))

        # feature completely within region
        self.assertFalse(before(f_150_170, r_100_200))
        self.assertFalse(before(f_100_170, r_100_200))
        self.assertFalse(before(f_150_200, r_100_200))

        # different chrs
        f_ch1_500_700 = Feature('chr1', 500, 700, '')
        f_ch7_50_70 = Feature('chr7', 50, 70, '')

        # feature is before region
        self.assertTrue(before(f_ch1_500_700, r_100_200))

        # feature is after region
        self.assertFalse(before(f_ch7_50_70, r_100_200))


    def test_after(self):
        """Make sure the after function works"""

        # Same chr...

        f_0_500 = Feature('chr3', 0, 500, '')
        f_50_75 = Feature('chr3', 50, 75, '')
        f_90_150 = Feature('chr3', 90, 150, '')
        f_99_201 = Feature('chr3', 99, 201, '')
        f_100_170 = Feature('chr3', 100, 170, '')
        f_100_200 = Feature('chr3', 100, 200, '')
        f_150_170 = Feature('chr3', 150, 170, '')
        f_150_200 = Feature('chr3', 150, 200, '')
        f_170_230 = Feature('chr3', 170, 230, '')
        f_300_500 = Feature('chr3', 300, 500, '')

        r_100_200 = Region('chr3', 100, 200, '')

        # feature is before region
        self.assertFalse(after(f_50_75, r_100_200))

        # feature is after region
        self.assertTrue(after(f_300_500, r_100_200))

        # feature is exactly the same as region
        self.assertFalse(after(f_100_200, r_100_200))

        # feature overlaps region to the left
        self.assertFalse(after(f_90_150, r_100_200))

        # feature overlaps region to the right
        self.assertFalse(after(f_170_230, r_100_200))

        # feature completely overlaps region
        self.assertFalse(after(f_0_500, r_100_200))
        self.assertFalse(after(f_99_201, r_100_200))

        # feature completely within region
        self.assertFalse(after(f_150_170, r_100_200))
        self.assertFalse(after(f_100_170, r_100_200))
        self.assertFalse(after(f_150_200, r_100_200))

        # different chrs
        f_ch1_500_700 = Feature('chr1', 500, 700, '')
        f_ch7_50_70 = Feature('chr7', 50, 70, '')

        # feature is before region
        self.assertFalse(after(f_ch1_500_700, r_100_200))

        # feature is after region
        self.assertTrue(after(f_ch7_50_70, r_100_200))


class TestOverlap(unittest.TestCase):
    """Make sure the overlaps and get_overlap functions work"""

    def setUp(self):
        pass

    def test_overlaps(self):

        # Same chr...

        f_0_500 = Feature('chr3', 0, 500, '')
        f_50_75 = Feature('chr3', 50, 75, '')
        f_90_150 = Feature('chr3', 90, 150, '')
        f_99_201 = Feature('chr3', 99, 201, '')
        f_100_170 = Feature('chr3', 100, 170, '')
        f_100_200 = Feature('chr3', 100, 200, '')
        f_150_170 = Feature('chr3', 150, 170, '')
        f_150_200 = Feature('chr3', 150, 200, '')
        f_170_230 = Feature('chr3', 170, 230, '')
        f_300_500 = Feature('chr3', 300, 500, '')

        r_100_200 = Region('chr3', 100, 200, '')

        # feature is before region
        self.assertFalse(overlaps(f_50_75, r_100_200))

        # feature is after region
        self.assertFalse(overlaps(f_300_500, r_100_200))

        # feature is exactly the same as region
        self.assertTrue(overlaps(f_100_200, r_100_200))

        # feature overlaps region to the left
        self.assertTrue(overlaps(f_90_150, r_100_200))

        # feature overlaps region to the right
        self.assertTrue(overlaps(f_170_230, r_100_200))

        # feature completely overlaps region
        self.assertTrue(overlaps(f_0_500, r_100_200))
        self.assertTrue(overlaps(f_99_201, r_100_200))

        # feature completely within region
        self.assertTrue(overlaps(f_150_170, r_100_200))
        self.assertTrue(overlaps(f_100_170, r_100_200))
        self.assertTrue(overlaps(f_150_200, r_100_200))

        # different chrs
        f_ch1_500_700 = Feature('chr1', 500, 700, '')
        f_ch7_50_70 = Feature('chr7', 50, 70, '')

        f_ch3_120_180 = Feature('chr3', 120, 180, '')
        f_ch3_120_180 = Feature('chr3', 120, 180, '')
        f_ch3_150_270 = Feature('chr7', 50, 70, '')

        # feature is before region
        self.assertFalse(overlaps(f_ch1_500_700, r_100_200))

        # feature is after region
        self.assertFalse(overlaps(f_ch7_50_70, r_100_200))

    def test_get_overlap(self):

        # Same chr...

        f_0_500 = Feature('chr3', 0, 500, '')
        f_50_75 = Feature('chr3', 50, 75, '')
        f_90_150 = Feature('chr3', 90, 150, '')
        f_99_201 = Feature('chr3', 99, 201, '')
        f_100_170 = Feature('chr3', 100, 170, '')
        f_100_200 = Feature('chr3', 100, 200, '')
        f_150_170 = Feature('chr3', 150, 170, '')
        f_150_200 = Feature('chr3', 150, 200, '')
        f_170_230 = Feature('chr3', 170, 230, '')
        f_300_500 = Feature('chr3', 300, 500, '')

        r_100_200 = Region('chr3', 100, 200, '')

        # feature is before region
        self.assertEquals(0, get_overlap(f_50_75, r_100_200))

        # feature is after region
        self.assertEquals(0, get_overlap(f_300_500, r_100_200))

        # feature is exactly the same as region
        self.assertEquals(101, get_overlap(f_100_200, r_100_200))

        # feature get_overlap region to the left
        self.assertEquals(51, get_overlap(f_90_150, r_100_200))

        # feature get_overlap region to the right
        self.assertEquals(31, get_overlap(f_170_230, r_100_200))

        # feature completely get_overlap region
        self.assertEquals(101, get_overlap(f_0_500, r_100_200))
        self.assertEquals(101, get_overlap(f_99_201, r_100_200))

        # feature completely within region
        self.assertEquals(21, get_overlap(f_150_170, r_100_200))
        self.assertEquals(71, get_overlap(f_100_170, r_100_200))
        self.assertEquals(51, get_overlap(f_150_200, r_100_200))

        # different chrs
        f_ch1_500_700 = Feature('chr1', 500, 700, '')
        f_ch7_50_70 = Feature('chr7', 50, 70, '')

        # feature is before region
        self.assertEquals(0, get_overlap(f_ch1_500_700, r_100_200))

        # feature is after region
        self.assertEquals(0, get_overlap(f_ch7_50_70, r_100_200))


class TestFindFeatures(unittest.TestCase):
    """Make sure the overlaps and get_overlap functions work"""

    def setUp(self):
        self.genes = []

        realistic_file = """#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	+	66999824	67210768	2	66999824,67208755,	67000051,67210768,	SGIP1,	NM_032291,
chr1	+	33546713	33585995	1	33546713,	33585995,	ADC,	NM_052998,
chr1	+	90098643	90185094	3	90098643,90152028,90178267,	90098881,90152170,90185094,	LRRC8C	NM_032270
chr1	-	173010359	173020103	3	173010359,173013078,173019880,	173010853,173013109,173020103,	TNFSF18	NM_005092
chr1	-	231727037	231747836	2	231727037,231747748,	231727546,231747836,	n/a	n/a
chr1	+	110881944	110889303	2	110881944,110888929,	110884890,110889303,	RBM15,	NM_001201545,
chr1	+	110881944	110889303	3	110881944,110888160,110888929,	110884890,110888271,110889303,	RBM15	NM_022768
chr1	-	144300511	144340773	3	144300511,144301189,144340525,	144300755,144301536,144340773,	n/a	n/a
chr1	-	144300516	144341755	4	144300516,144301189,144340525,144341669,	144300755,144301536,144340623,144341755,	n/a	n/a
chr1	-	144300511	144340773	2	144300511,144340525,	144301536,144340773,	n/a	n/a
chr1	-	150599521	150602098	2	150599521,150601889,	150600068,150602098,	ENSA,	NM_207168,
chr1	+	149279475	149291742	1	149279475,	149291742,	n/a	n/a
chr1	+	149553002	149553662	1	149553002,	149553662,	PPIAL4A,	NM_178230,
chr1	+	149553002	149553787	1	149553002,	149553787,	n/a	n/a
chr1	+	149553002	149553787	1	149553002,	149553787,	PPIAL4B	NM_001143883"""

        fake_file = """#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	+	0010	0020	2	10,16,	18,20,	gs_10_20,	rs1,
chr1	+	0033	0099	1	33,	99,	gs33_99,	rs2,
chr1	+	0222	0444	3	222,333,400,	300,330,444,	gs222_444	rs3
chr1	-	0500	0600	3	500,530,560,	520,550,600,	gs500_600	rs4
chr1	-	1000	2000	2	1000,1042,	1040,2000,	gs1000_2000	n/a
chr1	+	3000	4000	2	3000,3500,	3501,4000,	gs3000_4000_1,	rs5_1,
chr1	+	3000	4000	3	3000,3600,3800,	3400,3700,4000,	gs3000_4000_2	rs_1
chr1	-	4511	4773	3	4511,4689,4705,	4555,4696,4773,	n/a_1	n/a
chr1	-	4516	4755	4	4516,4689,4705,4740,	4555,4696,4723,4755,	n/a	n/a
chr1	-	4511	4773	2	4511,4625,	    4536,4773,	n/a_2	n/a
chr1	-	5521	5998	2	5521,5889,	57068,5998,	ENSA,	NM_207168,"""


        fake_file = fake_file.split('\n')

        print fake_file

        header = fake_file[0]
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

        for line in fake_file[1:]:

            line = line.strip() # Remove surrounding white space
            # Create a list of the line entries
            line_entries = line.split('\t')
            print line_entries
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
            self.genes.append(gene)


        sort_intervals(self.genes)

        for gene in self.genes:
            print_interval(gene)

        #print self.genes

        self.regions = []
        self.regions.append(Region('chr1', '100', '200', 'Region_100_200'))
        self.regions.append(Region('chr1', '40', '142', 'Region_40_142'))
        self.regions.append(Region('chr1', '1030', '2005', 'Region_1030_2005'))
        self.regions.append(Region('chr1', '4200', '4242', 'Region_4200_4242'))
        self.regions.append(Region('chr1', '5600', '5900', 'Region_5600_5900'))

        sort_intervals(self.regions)

        for region in self.regions:
            print_interval(region)


    def test_find_features(self):

        print "\nin test_find_features"
        #print self.genes
        find_features(self.regions, gene_list=self.genes)


if __name__ == '__main__':
    unittest.main()

