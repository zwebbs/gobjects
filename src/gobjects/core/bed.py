# File Name: bed.py
# Created By: ZW
# Created On: 2022-12-19
# Purpose: defines class for BED6 and BED12 style intervals inherited from
#  Intervals class for use in gobjects


# module imports
# ----------------------------------------------------------------------------
from .intervals import Interval
from dataclasses import dataclass
from typing import Union

# class definitions
# ----------------------------------------------------------------------------

# class Bed6() - child class of Interval() which adds score and strand.
# * complies with the BED6 standard found on the UCSC file format standards
#   webpage: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# 
# * the class adds two attributes to the Interval base class:
#    1. interval score (score); typically 1-1000 -required, missing denoted by '.'-
#    2. feature strandness (strand) -required, missingness denoted by '.'-
@dataclass
class Bed6(Interval):
    score: Union[int,str]
    strand: str

    # define a custom printout representation for the Bed6
    def __repr__(self):
        spec_bed6 = (
            f" {self.chrom} {self.chromStart} {self.chromEnd} {self.name}"
            f" {self.score} {self.strand} "
        )
        return repr(f"Bed6({spec_bed6})")
    

# class Bed12() - child class of Bed6() which adds several more fields.
# * complies with the BED12 standard found on the UCSC file format standards
#   webpage: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# 
# * the class adds six additional attributes to the Bed6 base class:
#    1. thickStart- The starting position at which the feature is drawn thickly (e.g. start codon)
#    2. thickEnd- The ending position at which the feature is drawn thickly (e.g. stop codon)
#    3. itemRgb- An RGB value of the form R,G,B (e.g. 255,0,0) to determine the display color
#    4. blockCount- The number of blocks (e.g. exons) in the BED line.
#    5. blockSizes- A comma-separated string of the block sizes.
#    6. blockStarts- A comma-separated string of block starts. positions should be relative to chromStart
@dataclass
class Bed12(Bed6):
    thickStart: int
    thickEnd: int
    itemRgb: str
    blockCount: int
    blockSizes: str
    blockStarts: str

    # define a custom printout representation for the Bed6
    def __repr__(self):
        spec_bed12 = (
            f" {self.chrom} {self.chromStart} {self.chromEnd} {self.name}"
            f" {self.score} {self.strand} "
            f" {self.thickStart} {self.thickEnd} {self.itemRgb}"
            f" {self.blockCount} {self.blockSizes} {self.blockStarts} "
        )
        return repr(f"Bed12({spec_bed12})")