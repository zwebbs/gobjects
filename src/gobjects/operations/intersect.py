# File Name: intersect.py
# Created By: ZW
# Created On: 2022-12-22
# Purpose: define functions for interval algebra ( much like those found
#  in bedtools ) specifically bedtools intersect-like functionality.

# module imports
# ----------------------------------------------------------------------------
from ..core.bed import Bed6, Bed12
from ..core.bedpe import Bedpe
from ..core.gtf import Gtf
from ..core.intervals import Interval

# constants definitions
# ----------------------------------------------------------------------------
CONTIGUOUS = [Interval, Bed6, Bed12, Gtf]
NON_CONTIGUOUS = [Bedpe]
ALL_GOBJECTS = CONTIGUOUS + NON_CONTIGUOUS

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

# define intersect_contiguous() which takes two contiguous features, like Bed or GTF
# and returns the boolean result of the interval intersect. the function
# has a strand-aware option but is not checked internally for strand attrs
def intersect_contiguous(gobjA, gobjB, strand_aware=False):
    if strand_aware and (gobjA.strand != gobjB.strand): return False
    else:  # regular intersect routine
        if gobjA.chrom != gobjB.chrom: return False
        else:
            return not (  # define NOT intersection conditions
                ((gobjA.zero_idx_end) < gobjB.zero_idx_start) or
                 (gobjA.zero_idx_start > (gobjB.zero_idx_end)))

# define intersect_noncontiguous() which takes two non-contiguous features, like Bedpe
# and returns the boolean result of the interval intersect. the function
# has a strand-aware option but is not checked internally for strand attrs
def intersect_noncontiguous(gobjA, gobjB, strand_aware=False):
    pass

# define intersect() which takes two objects and returns the boolean result of 
# interval intersect. this function wraps several subroutines which handle the various
# comparisons between and among object types, like BED vs BED or BED vs GTF, etc.
def intersect(a_obj, b_obj, strand_aware=False):
    if strand_aware:  # check for strand information on strand aware items
        check_stranded(a_obj) 
        check_stranded(b_obj)

    # handle comparisons between contiguous objects:
    if all((type(a_obj) in CONTIGUOUS, type(b_obj) in CONTIGUOUS)):
        return intersect_contiguous(a_obj, b_obj, strand_aware) 
    # handle noncontiguous to noncontiguous
    elif all((type(a_obj) in NON_CONTIGUOUS, type(b_obj) in NON_CONTIGUOUS)):
        return intersect_noncontiguous(a_obj, b_obj, strand_aware)
    else:
        raise TypeError(
            "Could not Infer the types for intersect objects "
            f"{a_obj} and/or {b_obj}. please ensure their type is in {ALL_GOBJECTS}"
        )

