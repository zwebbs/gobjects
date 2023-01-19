# File Name: writers.py
# Created By: ZW
# Created On: 2023-01-12
# Purpose: defines functions that write lists of gobjects to output text files
#  in the proper format. these writers in general are specific to a single 
#  record type and their corresponding file format.

# module imports
# ----------------------------------------------------------------------------
from gobjects.core import Bed6, Bedpe
from pathlib import Path


# function definitions
# ----------------------------------------------------------------------------

# define write_Bed6() which takes a list of bed6 records, and a file path and 
# writes the records out to a bed6 formatted file. optionally, you can append
# the records to a pre-existing file.
def write_Bed6(bed6_objects, filepath, fappend=False):
    # check that all objects in the passed list are Bed6
    check_types = [(type(b) == Bed6) for b in bed6_objects]
    if not all(check_types):
        raise TypeError("Error! not all record in passed list are Bed6\n")
    
    # open file pointer and ensure files to append exist
    fpath = Path(filepath).resolve()
    if not fappend: fobj = open(fpath, 'w')
    elif (fappend and fpath.exists()): fobj = open(fpath, 'a')
    else:
        msg = (
            f"Cannot find existing file {fpath} for append mode..\n"
            "plese ensure the file exists or switch to write mode"
        )
        raise FileNotFoundError(msg)
    
    # write records out sequentially
    print(f"Writing Bed6 records to {fpath}..")
    for record in bed6_objects:
        line = (
            f"{record.chrom}\t{record.chromStart}\t{record.chromEnd}\t"
            f"{record.name}\t{record.score}\t{record.strand}\n"
        )
        fobj.write(line)
    print(f"Done writing Bed6 records to {fpath}..")
    fobj.close()

# define write_Bedpe() which takes a list of bedpe records, and a file path and 
# writes the records out to a bedpe formatted file. optionally, you can append
# the records to a pre-existing file.
def write_Bedpe(bedpe_objects, filepath, fappend=False):
    # check that all objects in the passed list are Bedpe
    check_types = [(type(b) == Bedpe) for b in bedpe_objects]
    if not all(check_types):
        raise TypeError("Error! not all record in passed list are Bedpe\n")
    
    # open file pointer and ensure files to append exist
    fpath = Path(filepath).resolve()
    if not fappend: fobj = open(fpath, 'w')
    elif (fappend and fpath.exists()): fobj = open(fpath, 'a')
    else:
        msg = (
            f"Cannot find existing file {fpath} for append mode..\n"
            "plese ensure the file exists or switch to write mode"
        )
        raise FileNotFoundError(msg)
    
    # write records out sequentially
    print(f"Writing Bedpe records to {fpath}..")
    for record in bedpe_objects:
        line = (
            f"{record.chrom1}\t{record.chromStart1}\t{record.chromEnd1}\t"
            f"{record.chrom2}\t{record.chromStart2}\t{record.chromEnd2}\t"
            f"{record.name}\t{record.score}\t{record.strand1}\t{record.strand2}\t"
            f"{record.attributes}\n"
        )
        fobj.write(line)
    print(f"Done writing Bedpe records to {fpath}..")
    fobj.close()