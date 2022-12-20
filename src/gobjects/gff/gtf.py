# File Name: gtf.py
# Created By: ZW
# Created On: 2022-12-20
# Purpose: defines class for gtf style objects (identical to the GFF2 spec),
#  as defined by standards at https://genome.ucsc.edu/FAQ/FAQformat.html#format3

# module imports
# ----------------------------------------------------------------------------
from dataclasses import dataclass
from operator import lt, gt
from re import split as re_split
from re import match
from typing import Union

# function definitions
# ----------------------------------------------------------------------------

# function prep_chrom_compare() splits contig names that are a combo of
# strings and digits to use as sorting keys in sort functions and similar 
# processes. for example: 'chr1' -> ('chr', 1). this function can be 
# used in the sorted() function with the key=prep_chrom_comp argument
def prep_chrom_comp(chrom: str):
    def recode(substr: str):
        return int(substr) if substr.isdigit() else substr.lower()
    return [recode(sub) for sub in re_split('([0-9]+)',chrom)]

# function check_and_convert_nums() takes a string and converts it to 
# a floating point python number if it looks like an integer or float,
# otherwise, it returns the string back unchanged
def check_and_convert_nums(string):
    pattern = "^[+-]?((\d+(\.\d+)?)|(\.\d+))$"
    mth = match(pattern, string)
    return (float(mth.string) if mth is not None else string)

# class definitions
# ----------------------------------------------------------------------------

# TODO: Description for the GTF class and attributes!!

@dataclass
class GTF:
    chrom: str
    source: str
    feature: str
    chromStart: int
    chromEnd: int
    score: int
    strand: str
    frame: Union[str,int]
    attributes: str
    
    # define post_init routines to create attribute dictionary
    def __post_init__(self): # TODO
        self.attr_dict = {}
        pre_proc_attr = self.attributes.split("; ")
        for pair in pre_proc_attr:
            k,v = [p.strip('"') for p in pair]
            if k not in self.attr_dict:
                self.attr_dict[k] = [check_and_convert_nums(v)]
            else:
                self.attr_dict[k].append(check_and_convert_nums(v))

    # define custom representation when printing the GTF object
    def __repr__(self):
        spec = (
            f" {self.chrom} {self.source} {self.feature}"
            f" {self.chromStart} {self.chromEnd} {self.score}"
            f" {self.strand} {self.frame} {self.attributes} "
        )
        return repr(f"GTF({spec})")

        # define a custom function for the equal to (==) comparator
    # based on identical interval information on matching coordinates
    def __eq__(self, other):
        comp = ((self.chrom == other.chrom) and
                (self.chromStart == other.chromStart) and
                (self.chromEnd == other.chromEnd))
        return comp
    
    # define a custom function for the less than (<) comparator
    # based on interval algebra on matching chromosomes
    def __lt__(self, other):
        if lt(*[prep_chrom_comp(c) for c in [self.chrom, other.chrom]]): return True
        elif gt(*[prep_chrom_comp(c) for c in [self.chrom, other.chrom]]): return False
        else: # if the chromosomes names are equal by natural sort
            comp = ((self.chromEnd < other.chromStart) or
                    ((self.chromStart == other.chromStart) and
                    (self.chromEnd < other.chromEnd)))
            return comp
    
    # define a custom function for the greater than (>) comparator
    # based on interval algebra on matching chromosomes
    def __gt__(self, other):
        if lt(*[prep_chrom_comp(c) for c in [self.chrom, other.chrom]]): return False
        elif gt(*[prep_chrom_comp(c) for c in [self.chrom, other.chrom]]): return True
        else: # if the chromosomes names are equal by natural sort
            comp = (self.chromStart > other.chromStart)
            return comp
    
    # define a custom function for the less than or equal to (<=) comparator
    # based on niterval algebra on matching chromosomes
    def __le__(self,other):
        return (self.__lt__(other) or self.__eq__(other))
    
    # define a custom function for the greater than or equal to (>=) comparator
    # based on interval algebra on matching chromosomes
    def __ge__(self,other):
        return (self.__gt__(other) or self.__eq__(other))
