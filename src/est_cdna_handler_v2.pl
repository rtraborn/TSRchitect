#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

##-------------------------------------------------------------------------------
### Copyright (c) 2012, Krishnakumar Sridharan, Iowa State University
### Copyright (c) 2012, Brendel Group, Iowa State University and Indiana University
###
### Permission to use, copy, modify, and/or distribute this software for any
### purpose with or without fee is hereby granted, provided that the above
### copyright notice and this permission notice appear in all copies.
###
### THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
### WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
### MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
### ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
### WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
### ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
### OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
###-------------------------------------------------------------------------------

### TITLE --- est_cdna_handler.pl 
#### AUTHOR --- Krishnakumar Sridharan 09/18/2012
#### PURPOSE --- To run GeneSeqer on EST datasets and nucleotide sequence so as to associate 5'-Exon Start positions from GeneSeqer output with genes from GFF file. 
#### USAGE --- perl est_cdna_handler.pl -trans <cDNA_EST_dataset_to_be used_for_spliced_alignment> -seq <Nucleotide_sequence_to_be_used_for_spliced_alignment> -gff <GFF_annotation_file> -gsq [GeneSeqer_output_file] -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Cluster-format-file] -locn [Location_of_bin_folder_of_GeneSeqer] -species [Species_name]
#### PRE-REQUISITES --- 	1) "Getopt::Long" package available at CPAN
#				2) GeneSeqer package installed (available at http://www.plantgdb.org/)

#NOTE: Look at GeneSeqer usage to understand different Splice Site models for different species

my ($infile_trans,$infile_seq,$infile_gff,$infile_gsq,$infile_gth,$outfile,$gsq_locn,$gth_locn,$species,$dist);

GetOptions (    "trans=s" => \$infile_trans,		# cDNA/EST transcript dataset to be used for spliced alignment. Mandatory
		"seq=s" => \$infile_seq,		# Nucleotide sequence to be used for spliced alignment. Mandatory
                "gff=s" => \$infile_gff,                # GFF format annotation file. Mandatory
		"gsq=s" => \$infile_gsq,                # GeneSeqer output file. Optional
 		"gth=s" => \$infile_gth,		# GenomeThreader output file. Optional
		"dist=i" => \$dist,			# Distance from annotated gene start to be used for getting relevant transcripts. Optional
		"locn_gsq=s" => \$gsq_locn,		# Location of bin folder of GeneSeqer on system. Optional
		"locn_gth=s" => \$gth_locn,		# Location of bin folder of GenomeThreader on system. Optional
		"species=s" => \$species,		# Species from which sequence is derived for GeneSeqer. Optional
                "out=s" => \$outfile);                  # Output file name. Optional

my $usage = "\n\tUSAGE: perl est_cdna_handler.pl -trans <cDNA_EST_dataset_to_be used_for_spliced_alignment> -seq <Nucleotide_sequence_to_be_used_for_spliced_alignment> -gff <GFF_annotation_file> -gsq [GeneSeqer_output_file] -gth [GenomeThreader_output_file] -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Cluster-format-file] -locn_gsq [Location_of_bin_folder_of_GeneSeqer] -locn_gth [Location_of_bin_folder_of_GenomeThreader] -species [Species_name]\n
\t<>-Mandatory parameter\t\t[]-Optional parameter\n
NOTE: Please define all mandatory paramters and install (1) GeneSeqer package (available at http://www.plantgdb.org/) OR GenomeThreader package (available at http://www.genomethreader.org/) (2) Getopt::Long perl module from CPAN before proceeding. Also, see 0README\n\n";

if(defined $infile_gsq && $infile_gff) {}		##### If there is a readymade GSQ and GFF file, no need for all other commands
elsif(defined $infile_gth) {}
else { if(defined $infile_trans && defined $infile_seq && defined $infile_gff) {} else { print $usage; exit; } }   ##### Show usage if cmd-line parameters are undefined

if(defined $gsq_locn || defined $gth_locn) {} else { $gsq_locn = "/usr/local/src/GENESEQER/bin/"; }
if(defined $species) {} else { $species = "generic"; }
if(defined $dist) {} else { $dist=500; }

my @outname_temp;
my $outlabel;

if (defined $infile_seq) 
{ 
 @outname_temp = split (/\//,$infile_seq); 
 $outlabel = $outname_temp[$#outname_temp];
}
elsif (defined $infile_gsq)
{
 @outname_temp = split (/\//,$infile_gsq);
 $outlabel = $outname_temp[$#outname_temp];
}
elsif (defined $infile_gth)
{
 @outname_temp = split (/\//,$infile_gth);
 $outlabel = $outname_temp[$#outname_temp];
}

my $gsq_name = "gsq." . $outlabel;
my $gth_name = "gth." . $outlabel;
if(defined $outfile) {} else { $outfile = $outlabel. ".bed"}

if(defined $infile_gsq) { `perl gsq_gff_extractor.pl -gsq $infile_gsq -gff $infile_gff -dist $dist -out $outfile`; }
elsif(defined $infile_gth) { `perl gsq_gff_extractor.pl -gth $infile_gth -gff $infile_gff -dist $dist -out $outfile`; }
else 
{
 if(defined $gsq_locn)
 {
  `"$gsq_locn/"MakeArray $infile_trans`;
  `"$gsq_locn/"GeneSeqer -species $species -d $infile_trans -m 1000 -l $infile_seq -O $gsq_name`;
  `perl gsq_gff_extractor.pl -gsq $gsq_name -gff $infile_gff -dist $dist -out $outfile`;
#  `rm $gsq_name`;
 }
 elsif(defined $gth_locn)
 {
  `"$gth_locn/"gth -genomic $infile_seq -cdna $infile_trans -o $gth_name`;
  `perl gsq_gff_extractor.pl -gth $gth_name -gff $infile_gff -dist $dist -out $outfile`; 
#  `rm $gth_name`;
 }
}
