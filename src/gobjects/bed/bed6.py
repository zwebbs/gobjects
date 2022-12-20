# File Name: intervals.py
# Created By: ZW
# Created On: 2022-12-19
# Purpose: defines class for BED6 style intervals inherited from
#  Intervals class for use in gobjects


# module imports
# ----------------------------------------------------------------------------
from .intervals import Interval
from dataclasses import dataclass
from typing import Union

# class Bed6() - child class of Interval() which adds score and strand.
# * complies with the BED6 standard found on the UCSC file format standards
#   webpage: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# 
# * the class adds two attributes to the Interval base class:
#    1. interval score (score); typically 1-1000 -required, missing denoted by '.'-
#    2. feature strandness (strand) -required, missingness denoted by '.'-
# * the class also overwrites the intersect class to include an option for strand
#    aware intersection.
@dataclass
class Bed6(Interval):
    score: Union[int,str]
    strand: str

    # define a custom printout representation for the Bed6
    def __repr__(self):
        spec_int = f" {self.chrom} {self.chromStart} {self.chromEnd} {self.name}"
        spec_bed = spec_int + f" {self.score} {self.strand} "
        return repr(f"Bed6({spec_bed})")
    
    # define a modification of Interval.__intersect__ which provides 
    # an option for strand-aware intersections
    def __intersect__(self,other, strand_aware=False):
        if strand_aware and (self.strand != other.strand): return False
        else:  # regular intersect routine
            if self.chrom != other.chrom: return True
            else:
                return not (  # define intersection conditions
                    (other.chromEnd < self.chromStart) or 
                    (other.chromStart >= self.chromEnd))


