Gff-Scripts
===========

daveclarke.bioinfo@gmail.com

Scripts for modifying and parsing gffs

These scripts are created in order to assist in the manipulation of gff files.

--------------------------------
GfftoIsoforms.py
--------------------------------
version 1.1 12.10.2014

Description

A script that identifies overlapping CDS enteries in gff files
It creates the extra fields 'Cluster', 'ClusterName' and 'ClusterNum' in the info field (9th cloumn).
'Cluster' Identifies only the cluster number of which the CDS entry is a part
'ClusterNum' Identifies the number of the CDS entry within the cluster
'ClusterName' is a composite of the sequence identifier, Cluster and ClusterNum

input

python GfftoIsoforms.1.1.py -i inputfile.gff3 -o outputfile_without_extention ( -n field_in_gff_if_not_name = Name )

-i Gff3 file to be clustered

-o desired name of output files (two files created .gff3 and .txt)

-n name in info field used for unique identification of CDS regions (optional ; 'Name' is defult)

Additional output options are --rename, renameonly and --clusterall

--rename if specified then the 'ClusterName' data will instead be written into 'Name'
--renameonly will do the same as above but also not inclusde the 'Cluster' and 'ClusterNum' enteries
--clusterall if specified then all CDS lines will be assigned to clusters (clusters of one)

output

'outputfile_without_extention.gff3'
edited gff3 file, only the CDS lines will be eddited, with the 'Cluster' and 'ClusterNum' fields added.

'outputfile_without_extention.txt'
a text file with the numbers of the clusters with multiple enteries and the members of those clusters.

