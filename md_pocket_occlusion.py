#Python script for computing the fractional overlap of the PD-L2 Trp110 side chain.
#Part of the Supplemental information for:
#Nicolas A. Pabon & Carlos J. Camacho, CJ. "Probing protein flexibility reveals a mechanism for selective promiscuity" 

#original code written by Matthew Baumgartner
#1-11-16

import sys
import os
import subprocess as sp
import argparse
import math

#control how much debugging information is printed
DEBUG = 2

def main():
    
    parser = argparse.ArgumentParser(description = 'Compute the overlap of md frames to an input structure')
    parser.add_argument('-s', dest = 'residues', type = str, help = 'residues to check overlap of')
    parser.add_argument('-l', dest = 'frame_list',  help = 'list of frames')
    parser.add_argument('-ot', dest = 'outtext', help = 'output text file')
    parser.add_argument('-o', dest = 'outfile',  help = 'output pdb file')
    parser.add_argument('-start', dest = 'start', default = 0, type = float, help = 'optional, choose a starting from to start at')
    parser.add_argument('-end', dest = 'end', default = -1, type = float, help = 'optional, choose a starting from to end at')
    parser.add_argument('--redo', dest = 'redo', action = 'store_const', const = True, default = False, help = 'If specified, ignore any saved results and recompute (Default: False)')
   
    args = parser.parse_args()

    resfile = args.residues
    frame_list = args.frame_list
    outfile = args.outtext
    phe_outfile = args.outfile 
    start = int(args.start)
    end = int(args.end)
    REDO = args.redo

    if (resfile == None or frame_list == None or outfile == None or phe_outfile == None):
        print '\nERROR: Required arguments missing. Run with -h to show proper usage.\n'
        sys.exit()
    
    #https://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html#allatom
    vdw = {'C' : 1.7,
       'N' : 1.625,
       'O' : 1.48,
       'S' : 1.782,
       'H' : 1.00,
       'P' : 1.871,
       'F' : 1.56,
       'Cl' : 1.735,
       'CL' : 1.735,
       'Br' : 1.978,
       'BR' : 1.978,
       'I' : 2.094}
    
    if not REDO:
        if os.path.exists(outfile):
            print 'Output file:', outfile, 'already exists, not recomputing!\nRun with --redo to force recompute'
            sys.exit()
        

    res_atoms = [ map(float,[ line[30:38], line[38:46], line[46:54]]) for line in open(resfile).readlines() if (line.startswith('ATOM') or line.startswith('HETATM')) ]
    
    #proceed each atom name with the resid to avoid collisions
    res_atom_names = [ line[12:16].strip() + '_' + line[22:26].strip() for line in open(resfile).readlines() if (line.startswith('ATOM') or line.startswith('HETATM')) ] 
    
    #Get the atom tyes from the last column
    res_atom_type = [ line[76:80].strip().strip('-+0123456789') for line in open(resfile).readlines() if (line.startswith('ATOM') or line.startswith('HETATM')) ]
    
    #read in text file that lists the absolute paths of the pdb files of the frames. I used mdpocket from fpocket to split the simulation and create the input file
    frames = [ line.strip() for line in open(frame_list).readlines() ]
    
    #the frames in which an atom is over lapped
    frame_overlaps = { nm : [] for nm in res_atom_names }
    
    for frame_id, frame in enumerate(frames):
        
        #check to see if you are inside the range you want to be
        if frame_id < start:
            continue
        if end != -1:
            if frame_id >= end:
                continue
        
        
        if frame_id > 0:
            if DEBUG > 2:
                print 'frame:', frame_id, 'frame:', frame
            
            for line in open(frame).xreadlines():
                if line.startswith('ATOM'):
                    atom_name = line[12:16].strip()
                    if not atom_name.startswith('H'): #skip hydrogens
                        traj_coor = map(float, [line[30:38], line[38:46], line[46:54]])
                        
                        for i,res_coor in enumerate(res_atoms):
                            res_name = res_atom_names[i]
                            res_vdw = vdw[res_atom_type[i]] #vdw of the first letter (element)
                            
                            #distance
                            dist = math.sqrt((res_coor[0] - traj_coor[0])**2 + (res_coor[1] - traj_coor[1])**2 + (res_coor[2] - traj_coor[2])**2)
                            
                            #if it is closer than the vdw
                            if dist < res_vdw:
                                #get the residue number of the overlapping atom
                                traj_resid = int(line[22:26].strip())
                                
                                #if you haven't added one already
                                if len(frame_overlaps[res_name]) < frame_id:

                                    frame_overlaps[res_name].append(traj_resid)
                                
#                                break
                    
            #if you didn't find any overlapping atoms, add a 0 to the list
            for res_name in res_atom_names:
                if DEBUG > 2:
                    print res_name, 'overlap len:', len(frame_overlaps[res_name])
                if len(frame_overlaps[res_name]) < frame_id:
                    frame_overlaps[res_name].append(0)
                    
        
    atom_names = []
    lists = []
    
    for res_name, overlapped in frame_overlaps.iteritems():
        if DEBUG > 2:
            print res_name, ':',  overlapped   
            print len(overlapped)
        atom_names.append(res_name)
        lists.append(overlapped)
        
        
    fout = open(outfile, 'w')
    
    
    #right now there is a bug where if you aren't looking at all of the frames, the following will crash. Because I'm not using it, I'm just going to ignore it for now
    if start == 0 and end == -1:
    
        #print the header
        header = 'Frame\t' + '\t'.join(atom_names) + '\n'
        fout.write(header)
        
        for i in range(len(lists[0])): #
            str_line = str(i)
            
            for lst in lists:
                str_line += '\t' + str(lst[i])
                
            fout.write(str_line + '\n')
        fout.close()
        
        if DEBUG > 2:
            print outfile, 'written'
            
            
        
    #write out the residue again but put the perc overlapped in the b factor field
    name_perc = {}
    if DEBUG > 2:
        print 'Calculating the percent overlapped'
    for i,lst in enumerate(lists):
        non_zero = len([ x for x in lst if x != 0 ])
        perc_overlapped = float(non_zero) / len(lst)
        if DEBUG > 2:
            print atom_names[i], perc_overlapped
        name_perc[atom_names[i]] = perc_overlapped
        
    
    output = []
    
    for line in open(resfile).readlines():
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            output.append(line)
        elif line[76:80].strip().strip('-+0123456789').startswith('H'):
            output.append(line)
        else:
            atom_name = line[12:16].strip() + '_' + line[22:26].strip()
            b_fac = '  ' + '%03.2f' % round(name_perc[atom_name],2)
        
            new_line = line[:60] + b_fac + line[66:]
            output.append(new_line)
            
    fout = open(phe_outfile, 'w')
    fout.write(''.join(output))
    fout.close()
    
    print phe_outfile, 'written'
        
        
        
        
        
if __name__ == '__main__':
    main()
        
        
