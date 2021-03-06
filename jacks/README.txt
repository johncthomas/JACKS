To install:
------------

1.  At command prompt, cd to this directory
2.  pip install .
3.  Add jacks folder to system path e.g. EXPORT PATH=$PATH:path_to_jacks/JACKS/jacks

To check installation:
----------------------

Run python, type 'import jacks' and check that no error is thrown.


To run JACKS on full screen data:
----------------------------

python run_JACKS.py countfile replicatemapfile:replicate_hdr:sample_hdr:ctrl_sample sgrnamappingfile:sgrna_hdr:gene_hdr outprefix

where

countfile  - A tab or comma (if comma, must end in '.csv') delimited file containing the raw counts of each guide in each cell line. 
The first column should contain the guide ids.
The second column can optionally contain the gene mappings, or not.
The remaining columns should contain the count data, with the first row containing the column headings, which should be the replicate identifiers, which should match those in the replicate map file.
Only those columns with an entry in the replicatemap file below will be processed.

replicatemapfile:replicate_hdr:sample_hdr:ctrl_sample_or_hdr (4 items, separated by colons)

replicatemapfile is a tab or comma (if comma, must end in '.csv') delimited file containing the mappings from replicates to samples.
The file can contain other columns, replicate_hdr specifies the column header of the column containing the replicate identifiers (matching the column headers in the count file).
sample_hdr specifies the column header of the column containing the sample mappings for each replicate (i.e. an identifier for the cell line or condition)
ctrl_sample_or_hdr specifies the sample identifier of the sample which is to be used as a control by JACKS (and which can contain multiple replicates),
or alternatively ctrl_sample_or_hdr can specify the column header of the column in the replicatemapfile which contains the sample identifiers (as used in the sample_hdr column)
of the control for each sample (e.g. see example/example_repmap_matched_ctrls.tab). For control samples, the control and sample identifiers should be identical.
    
sgrnamappingfile:sgrna_hdr:gene_hdr (3 items, separated by colons)

sgrnamappingfile is a tab or comma (if comma, must end in '.csv') delimited file containing the mappings from guides to genes.
The file can contain other columns, sgrna_hdr specifies the column header of the column containing the guide identifiers.
gene_hdr specifies the column header of the column containing the gene identifiers.
If the count file has a gene column, the count file can be reused here.

outprefix: the output prefix of the JACKS output files. Three output files will be produced.
    outprefix_gene_JACKS_results.txt contains the gene essentiality scores E(w) for each cell line
    outprefix_genestd_JACKS_results.txt contains the standard deviations of the gene essentiality scores std(w) for each cell line
    outprefix_grna_JACKS_results.txt contains the gRNA efficacy scores E(X) and E(X^2) for each guide
    outprefix_JACKS_full_data.pickle is a pickle file containing the full screen results

examples:

python run_JACKS.py example/example_count_data.tab example/example_repmap.tab:Replicate:Sample:CTRL example/example_count_data.tab:sgRNA:gene example_jacks/example_jacks

python run_JACKS.py example/example_count_data.tab example/example_repmap_matched_ctrls.tab:Replicate:Sample:Control example/example_count_data.tab:sgRNA:gene example_jacks/example_jacks

    
OR to run JACKS on new screen with previously used library:
----------------------------------------------------------

python run_JACKS_reference.py countfile replicatemapfile:replicate_hdr:sample_hdr:ctrl_sample sgrnamappingfile:sgrna_hdr:gene_hdr grnaeffiacyfile outprefix

all arguments as above except:
   
grnaeffiacyfile contains previously trained gRNA efficacy value. This is the file outprefix_grna_JACKS_results.txt returned above, or one 
    
example:   

python run_JACKS_reference.py  example/example_count_data.tab example/example_repmap.tab:Replicate:Sample:CTRL example/example_count_data.tab:sgRNA:gene example/example_grna_JACKS_results.txt example_jacks_ref/example_jacks_ref
    
    
see 2018_paper_materials/README.txt for further examples.

Then, to plot heatmap outputs for a gene of interest
-----------------------------------------------

python plot_heatmap.py picklefile gene(or "random" to randomly select one) outfile

where:

picklefile:  the pickle file output by run_JACKS.py

gene: the name of a gene to plot the output for (or "random" to select a random one)

outfile: the name of a .png output file to write the figure to

example:

python plot_heatmap.py example_jacks/example_jacks_JACKS_results_full.pickle KRAS example_jacks/KRAS.png
 