Gff-Scripts
===========

daveclarke.bioinfo@gmail.com

Scripts for modifying and parsing GFF files.

These scripts are created in order to assist in the manipulation of GFF files.

--------------------------------
GfftoIsoforms.py
--------------------------------
Version 1.1 12.10.2014

Description

A script that identifies overlapping CDS entries in gff files
It creates the extra fields 'Cluster', 'ClusterName' and 'ClusterNum' in the info field (9th column).

'Cluster' identifies only the cluster number of which the CDS entry is a part

'ClusterNum' Identifies the number of the CDS entry within the cluster

'ClusterName' is a composite of the sequence identifier, Cluster and ClusterNum

Input

usage: GfftoIsoforms.py [-h] [-o OUTFILE] [-n IDFIELD] [-i INPUT] [--rename]
                        [--renameonly] [--clusterall]

Required arguments:
  -i INPUT, --input INPUT
                        Gff3 input file to be clustered.
  -o OUTFILE, --outfile OUTFILE
                        Desired name of output files (two files created .gff3
                        and .txt)

Optional arguments:
  -h, --help            show this help message and exit
  -n IDFIELD, --idfield IDFIELD
                        Field in gff3 notes column to use as unique identifier
                        for CDS annotations, defaults to 'Name'.
  --rename              If specified then the 'ClusterName' data will instead
                        be written into 'Name' field.
  --renameonly          Same as 'rename option' but will not also include
                        'Cluster' and 'ClusterNum' as additional note
                        entries.
  --clusterall          If specified then all CDS lines will be assigned to
                        clusters (clusters of one).

Output

'outputfile_without_extention.gff3'
Edited gff3 file, only the CDS lines will be edited, with the 'Cluster' and 'ClusterNum' fields added.

'outputfile_without_extention.txt'
A text file with the numbers of the clusters with multiple entries and the members of those clusters.


