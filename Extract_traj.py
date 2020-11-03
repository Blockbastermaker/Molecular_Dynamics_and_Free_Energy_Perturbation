#!/usr/bin/env python
#Convert a Q restart file to a pdb snapshot
#
#Nicolas Panel, Carlsson Group, 2019
#
#usage: get_snapshot.py [-h] -f QFILE [-l LIB] -l1 L1 -l2 L2 [-t TOP]
#                       [-q QPREP]
#Arguments:
#  -h, --help   show this help message and exit
#  -f  QFILE    restart file in Q format
#  -l  LIB      .lib used during system preparation
#  -l1 L1       .lib file of the first ligand
#  -l2 L2       .lib file of the second ligand
#  -t  TOP      topology file of the system
#  -q  QPREP    path to qprep program 


import os
import argparse

def write_input(lib,l1,l2,top,qfile):
    """
    Create Qprep input file for snapshot extraction
    """
    if os.path.splitext(qfile)[1] != ".re":
        print "ERROR: %s should be a Q restart file with .re extension"
        return 0

    #Set pdb output name
    outname = os.path.splitext(qfile)[0]+".pdb"

    #Write qprep input file
    filout = open("get_snapshot.inp","w")                                              
    filout.write("rl %s\n" % lib)                 
    filout.write("rl %s\n" % l1)
    filout.write("rl %s\n" % l2)
    filout.write("rt %s\n" % top)                                                 
    filout.write("mask none\n")                                                      
    filout.write("mask not excluded\n")                                              
    filout.write("rx %s\n" % qfile)                                             
    filout.write("wp %s n\n" % outname)                                          
    filout.write("quit\n")                                                           
    filout.close() 
    return outname

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', action='store', dest='qfile',
                    help='Restart file in Q format', required=True)
    parser.add_argument('-l', action='store', dest='lib',
                    default='/home/nourka/qligfep_Nicolas/qligfep/FF/OPLS2015.lib', 
                    help='General lib file', type=str)
    parser.add_argument('-l1', action='store', dest='l1', required=True,    
                    help='Ligand 1 lib file', type=str)
    parser.add_argument('-l2', action='store', dest='l2', required=True,    
                    help='Ligand 2 lib file', type=str)
    parser.add_argument('-t', action='store', dest='top',                   
                    default='dualtop.top',       
                    help='Topology file', type=str)
    parser.add_argument('-q', action='store', dest='qprep',                   
                    default='/home/nourka/Q5_BAR/bin/qprep5',                                           

                    help='Topology file', type=str)
    args = parser.parse_args() 


    #Create qprep input
    out = write_input(args.lib,args.l1,args.l2,args.top,args.qfile)

    #Run Qprep
    os.system("%s < get_snapshot.inp > get_snapshot.out" % args.qprep)
    print "%s has been created" % out
