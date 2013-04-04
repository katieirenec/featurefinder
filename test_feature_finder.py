
import random
import unittest
#import StringIO


from feature_finder import *

# NOTE these tests run in the order of the class names
# So, TestOverlap(...) runs before Testregion_list(...)

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


class Testregion_list(unittest.TestCase):

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
       pass


    def test_find_features_big_gene(self):

        print "\n*In test_find_features_big_gene"

        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #-------=====================================================----------
        #----------------------------------------------------------------------
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Regions:
        #-====------======----=====----=====----======----=====----=====---===-
        # Overlap:
        #           ======    =====    =====    ======    =====    ==

        gene_file = \
"""#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	+	8	60	1	8,	60,	G1,	G1,"""

        region_file = \
"""browser position chr1:26339-26399
track name="boygirl_12"
#Chromosome	StartPosition	EndPosition	RegionName	Score	nProbes	maxMinZvalue	maxPropDiff	avgTreatP	avgConP	
chr1	2	5	R1	4.1	1	4.1	0.2	0.7	0.5
chr1	12	17	R2	-4.0	1	-4.0	-0.2	0.0	0.3	
chr1	22	26	R3	-11.7	2	-4.2	-0.2	0.4	0.7	
chr1	31	35	R4	-4.0	1	-4.0	-0.2	0.3	0.6	
chr1	40	45	R5	-4.3	1	-4.3	-0.3	0.5	0.8	
chr1	50	54	R6	4.6	1	4.6	0.3	0.5	0.2	
chr1	59	63	R7	-4.1	1	-4.1	-0.2	0.2	0.5	
chr1	67	69	R8	4.0	1	4.0	0.2	0.6	0.4	"""

        import StringIO
        gene_fp = StringIO.StringIO(gene_file)
        gene_list = create_gene_list(gene_fp)

        region_fp = StringIO.StringIO(region_file)
        region_list = create_region_list(region_fp)

        find_features(region_list, gene_list=gene_list)

    def test_find_features_big_gene_neg(self):

        print "\n*In test_find_features_big_gene_neg"


        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #----------------------------------------------------------------------
        #-------=====================================================----------
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Regions:
        #-====------======----=====----=====----======----=====----=====---===-
        # Overlap:
        #           ======    =====    =====    ======    =====    ==

        gene_file = \
"""#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	-	8	60	1	8,	60,	G1,	G1,"""

        region_file = \
"""browser position chr1:26339-26399
track name="boygirl_12"
#Chromosome	StartPosition	EndPosition	RegionName	Score	nProbes	maxMinZvalue	maxPropDiff	avgTreatP	avgConP	
chr1	2	5	R1	4.1	1	4.1	0.2	0.7	0.5
chr1	12	17	R2	-4.0	1	-4.0	-0.2	0.0	0.3	
chr1	22	26	R3	-11.7	2	-4.2	-0.2	0.4	0.7	
chr1	31	35	R4	-4.0	1	-4.0	-0.2	0.3	0.6	
chr1	40	45	R5	-4.3	1	-4.3	-0.3	0.5	0.8	
chr1	50	54	R6	4.6	1	4.6	0.3	0.5	0.2	
chr1	59	63	R7	-4.1	1	-4.1	-0.2	0.2	0.5	
chr1	67	69	R8	4.0	1	4.0	0.2	0.6	0.4	"""

        import StringIO
        gene_fp = StringIO.StringIO(gene_file)
        gene_list = create_gene_list(gene_fp)

        region_fp = StringIO.StringIO(region_file)
        region_list = create_region_list(region_fp)

        find_features(region_list, gene_list=gene_list)


    def test_find_features_big_reg(self):

        print "\n*In test_find_features_big_reg"

        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #-====------======----=====----=====----======----=====----=====---===-
        #----------------------------------------------------------------------
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Regions:
        #-------=====================================================----------
        # Overlap:
        #           ======    =====    =====    ======    =====    ==

        gene_file = \
"""#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	+	2	5	1	2,	5,	G1,	G1,
chr1	+	12	17	1	12,	17,	G2,	G2,
chr1	+	22	26	1	22,	26,	G3,	G3,
chr1	+	31	35	1	31,	35,	G4,	G4,
chr1	+	40	45	1	40,	45,	G5,	G5,
chr1	+	50	54	1	50,	54,	G6,	G6,
chr1	+	59	63	1	59,	63,	G7,	G7,
chr1	+	67	69	1	67,	69,	G8,	G8,"""

        region_file = \
"""browser position chr1:26339-26399
track name="boygirl_12"
#Chromosome	StartPosition	EndPosition	RegionName	Score	nProbes	maxMinZvalue	maxPropDiff	avgTreatP	avgConP	
chr1	8	60	R1	-4.0	1	-4.0	-0.2	0.0	0.3	"""

        import StringIO
        gene_fp = StringIO.StringIO(gene_file)
        gene_list = create_gene_list(gene_fp)

        region_fp = StringIO.StringIO(region_file)
        region_list = create_region_list(region_fp)

        find_features(region_list, gene_list=gene_list)


    def test_find_features_big_reg_neg(self):

        print "\n*In test_find_features_big_reg_neg"

        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #----------------------------------------------------------------------
        #-====------======----=====----=====----======----=====----=====---===-
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Regions:
        #-------=====================================================----------
        # Overlap:
        #           ======    =====    =====    ======    =====    ==

        gene_file = \
"""#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	geneSymbol	refseq
chr1	+	2	5	1	2,	5,	G1,	G1,
chr1	+	12	17	1	12,	17,	G2,	G2,
chr1	+	22	26	1	22,	26,	G3,	G3,
chr1	+	31	35	1	31,	35,	G4,	G4,
chr1	+	40	45	1	40,	45,	G5,	G5,
chr1	+	50	54	1	50,	54,	G6,	G6,
chr1	+	59	63	1	59,	63,	G7,	G7,
chr1	+	67	69	1	67,	69,	G8,	G8,"""

        region_file = \
"""browser position chr1:26339-26399
track name="boygirl_12"
#Chromosome	StartPosition	EndPosition	RegionName	Score	nProbes	maxMinZvalue	maxPropDiff	avgTreatP	avgConP	
chr1	8	60	R1	-4.0	1	-4.0	-0.2	0.0	0.3	"""

        import StringIO
        gene_fp = StringIO.StringIO(gene_file)
        gene_list = create_gene_list(gene_fp)

        region_fp = StringIO.StringIO(region_file)
        region_list = create_region_list(region_fp)

        find_features(region_list, gene_list=gene_list)


        """def test_find_features_no_overlap(self):


        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #----=====-----=====-----=====-----=====-----=====-----=====-----=====-
        #----------------------------------------------------------------------
        # Regions:
        #==---------===-------==-------====-------===-----=====-------==------
        # Overlap:
        #

        gene_list = [ Gene(1, '+',  5,  9, 1,  5,  9, 'G1', 'G1'),
                      Gene(1, '+', 15, 19, 1, 15, 19, 'G2', 'G2'),
                      Gene(1, '+', 25, 29, 1, 25, 29, 'G3', 'G3'),
                      Gene(1, '+', 35, 39, 1, 35, 39, 'G4', 'G4'),
                      Gene(1, '+', 45, 49, 1, 45, 49, 'G5', 'G5'),
                      Gene(1, '+', 55, 59, 1, 55, 59, 'G6', 'G6'),
                      Gene(1, '+', 65, 69, 1, 65, 69, 'G7', 'G7') ]

        region_list = [ Region(1,  1,  2, 'R1', 'R1'),
                        Region(1, 12, 14, 'R2', 'R2'),
                        Region(1, 22, 23, 'R3', 'R3'),
                        Region(1, 31, 34, 'R4', 'R4'),
                        Region(1, 42, 44, 'R5', 'R5'),
                        Region(1, 50, 54, 'R6', 'R6'),
                        Region(1, 62, 63, 'R7', 'R7') ]

        find_features(region_list, gene_list=gene_list)

    def test_find_features_some_overlap(self):


        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #----======------================----=========--===----------=======---
        #----------------------------------------------------------------------
        # Regions:
        #======-----===-----=====----------===================-----============
        # Overlap:
        #    ==             =====            =========  ===          =======

        gene_list = [ Gene(1, '+',  5, 10, 1,  5, 10, 'G1', 'G1'),
                      Gene(1, '+', 17, 32, 1, 17, 32, 'G2', 'G2'),
                      Gene(1, '+', 37, 45, 1, 37, 45, 'G3', 'G3'),
                      Gene(1, '+', 48, 50, 1, 48, 50, 'G4', 'G4'),
                      Gene(1, '+', 61, 67, 1, 61, 67, 'G5', 'G5') ]

        region_list = [ Region(1,  1,  2, 'R1', 'R1'),
                        Region(1, 12, 14, 'R2', 'R2'),
                        Region(1, 22, 23, 'R3', 'R3'),
                        Region(1, 31, 34, 'R4', 'R4'),
                        Region(1, 42, 44, 'R5', 'R5'),
                        Region(1, 50, 54, 'R6', 'R6'),
                        Region(1, 62, 63, 'R7', 'R7') ]

        find_features(self.region_list, gene_list=self.gene_list)

    def test_find_features_realistic_genes(self):

        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #======-----===-----=====----------===================-----============
        #----------------------------------------------------------------------
        # Regions:
        #----======------================----=========--===----------=======---
        # Overlap:
        #    ==             =====            =========  ===          =======


        #    5    10   15   20   25   30   35   40   45   50   55   60   65  70
        #----|----|----|----|----|----|----|----|----|----|----|----|----|----|
        # Genes:
        #----======-----=====-----=====-----=====-----=====-----=====-----=====
        #----------------------------------------------------------------------
        # Regions:
        #=====-------=====------=====------=====------=====------=====---------
        # Overlap:
        #    =          ==        ===       ====      =====      ====



        self.region_list = []
        self.region_list.append(Region('chr1', '100', '200', 'Region_100_200'))
        self.region_list.append(Region('chr1', '40', '142', 'Region_40_142'))
        self.region_list.append(Region('chr1', '1030', '2005', 'Region_1030_2005'))
        self.region_list.append(Region('chr1', '4200', '4242', 'Region_4200_4242'))
        self.region_list.append(Region('chr1', '5600', '5900', 'Region_5600_5900'))

        sort_intervals(self.region_list)"""


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

        import StringIO

        gene_fp = StringIO.StringIO(fake_file)

        self.gene_list = create_gene_list(gene_fp)


if __name__ == '__main__':
    unittest.main()

