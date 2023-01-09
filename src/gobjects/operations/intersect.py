# File Name: intersect.py
# Created By: ZW
# Created On: 2022-12-22
# Purpose: define functions for interval algebra ( much like those found
#  in bedtools ) specifically bedtools intersect-like functionality.

# module imports
# ----------------------------------------------------------------------------
from ..core.intervals import Interval
from ..core.bed import Bed6, Bed12
from ..core.gtf import GTF

# function definitions
# ----------------------------------------------------------------------------

# define check_stranded() which takes an object and ensures it contains the
# self.strand attribute. If it does, the function returns None, if it doesn't
# the function raises an Attribute error.
def check_stranded(gobj):
    if not hasattr(gobj, "strand"):
        msg = (
            "Error! attemped strand-aware comparison"
            " for object with no strand information"
            )
        raise AttributeError(msg)
    else:
        return None

# define intersect_BEDtoBED() which takes two BED-interval-like objects
# and returns the boolean result of the interval intersect. the function
# has a strand-aware option but is not checked internally for strand attrs
def intersect_BEDtoBED(gobjA, gobjB, strand_aware=False):
    if strand_aware and (gobjA.strand != gobjB.strand): return False
    else:  # regular intersect routine
        if gobjA.chrom != gobjB.chrom: return True
        else:
            return not (  # define NOT intersection conditions
                (gobjB.chromEnd <= gobjA.chromStart) or 
                (gobjB.chromStart >= gobjA.chromEnd))

# define intersect_GTFtoGTF() which takes two GTF-like objects
# and returns the boolean result of the interval intersect. the funciton has a 
# strand-aware option but is not checked internally for strand attrs
def intersect_GTFtoGTF(gobjA, gobjB, strand_aware=False):
    if strand_aware and (gobjA.strand != gobjB.strand): return False
    else :
        return not (  # define NOT intersection conditions
                (gobjB.chromEnd < gobjA.chromStart) or 
                (gobjB.chromStart > gobjA.chromEnd)) # gtfs are right inclusive

# define intersect_BEDtoGTF() which takes an BED6/Interval-like object as its first
# argument and a GTF-like object as its second argument and returns the boolean 
# result of the interval intersect. the function has a strand aware option but is
# not checked internally for strand attrs:
def intersect_BEDtoGTF(gobjA, gobjB, strand_aware=False):
    if strand_aware and (gobjA.strand != gobjB.strand): return False
    else :
        return not (  # define NOT intersection conditions
                (gobjB.chromEnd < gobjA.chromStart) or 
                (gobjB.chromStart >= gobjA.chromEnd)) # take care of uneven boundary inclusion

# define intersect() which takes two objects and returns the boolean result of 
# interval intersect. this function wraps several subroutines which handle the various
# comparisons between and among object types, like BED vs BED or BED vs GTF, etc.
def intersect(a_obj, b_obj, strand_aware=False):
    if strand_aware:  # check for strand information on strand aware items
        check_stranded(a_obj) 
        check_stranded(b_obj)

    # handle comparisons between different object types
    if (type(a_obj) == Interval) and (type(b_obj) == Interval):
        return intersect_BEDtoBED(a_obj, b_obj, strand_aware)
    elif (type(a_obj) == GTF) and (type(b_obj) == GTF):
        return intersect_GTFtoGTF(a_obj, b_obj, strand_aware)
    elif (type(a_obj) == Interval) and (type(b_obj) == GTF):
        return intersect_BEDtoGTF(a_obj, b_obj, strand_aware)
    elif (type(a_obj) == Interval) and (type(b_obj) == GTF):
        return intersect_BEDtoGTF(b_obj, a_obj, strand_aware) # switch a and b for gtf to bed 
    else:
        raise TypeError(
            f"Could not Infer the types for intersect objects "
            "{a_obj} and/or {b_obj}"
        )