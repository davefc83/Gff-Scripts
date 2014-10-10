Gff-Scripts
===========

daveclarke.bioinfo@gmail.com

Scripts for modifying and parsing gffs

These scripts are created in order to assist in the manipulation of gff files.

--------------------------------
GfftoIsoforms.1.1.py
--------------------------------
Description

A script that identifies overlapping CDS enteries in gff files and creates the extra fields 'Cluster' and 'ClusterNum' in the info field (9th cloumn).
'ClusterNum' Identifies only the cluster number of which the CDS entry is a part
'Cluster' is a composite of the ClusterNum, sequence identifier and number designated as part of the cluster.

input

python GfftoIsoforms.1.1.py -i inputfile.gff3 -o outputfile_without_extention ( -n field_in_gff_if_not_name = Name )

-i Gff3 file to be clustered

-o desired name of output files (two files created .gff3 and .txt)

-n name in info field used for unique identification of CDS regions (optional ; 'Name' is defult)

output

'outputfile_without_extention.gff3'
edited gff3 file, only the CDS lines will be eddited, with the 'Cluster' and 'ClusterNum' fields added.

'outputfile_without_extention.txt'
a text file with the numbers of the clusters with multiple enteries and the members of those clusters.

