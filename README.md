# md_pocket_occlusion
SYNOPSIS

md_pocket_occlusion.py is a Python script written to measure the occlusion of a particular binding pocket on a protein interface in an MD simulation of the protein. The pocket location and volume are defined by a set of atoms in PDB format that the user uploads separately, which, in practice, would be the atoms of the bound ligand that bury themselves in the pocket. 


USAGE

md_pocket_occlusion.py [-h] [-s RESIDUES] [-l FRAME_LIST] [-ot OUTTEXT]
                       [-o OUTFILE] [-start START] [-end END] [--redo]


optional arguments:
  -h, --help     show this help message and exit
  -s RESIDUES    residues to check overlap of
  -l FRAME_LIST  list of frames 
  -ot OUTTEXT    output text file
  -o OUTFILE     output pdb file
  -start START   optional, choose a starting from to start at
  -end END       optional, choose a starting from to end at
  --redo         If specified, ignore any saved results and recompute
                 (Default: False)


INPUT

RESIDUES - The PDB file containing the buried ligand atoms that define the binding pocket. Each atom in this file will have it's overlap computed as a proxy measure of pocket occlusion.

RRAME_LIST - A text file listing the MD frame files that will be used in overlap calculations. Each frame must be saved as a separate PDB file and its path must be listed on it's own line in the FRAME_LIST file. EX:

    /path/to/frame1.pdb
    /path/to/frame2.pdb
    /path/to/frame3/pdb
    ...

OUTPUT

OUTTEXT - A text log file listing which atoms of reference RESIDUES are overlapped in each MD frame. The file is setup as a large table with rows representing the MD frames and columns representing the atoms in reference RESIDUES. The entries in the table correspond to the residue number of the MD atoms that overlapped the given reference atom in the given frame. Entries of 0 indicate that the reference atom was not overlapped.

OUTFILE - A PDB file containing the same atoms as the reference RESIDUES but where the last column in the PDB file indicates the percentage of MD frames that each atom was overlapped.
