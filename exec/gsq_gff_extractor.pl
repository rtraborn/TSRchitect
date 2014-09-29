#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

#-------------------------------------------------------------------------------
## Copyright (c) 2012, Krishnakumar Sridharan, Iowa State University
## Copyright (c) 2012, Brendel Group, Iowa State University and Indiana University
##
## Permission to use, copy, modify, and/or distribute this software for any
## purpose with or without fee is hereby granted, provided that the above
## copyright notice and this permission notice appear in all copies.
##
## THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
## WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
## ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
## WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
## ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
## OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
##-------------------------------------------------------------------------------

### TITLE --- gsq_gff_extractor.pl 
### AUTHOR --- Krishnakumar Sridharan 09/16/2012
### PURPOSE --- To associate 5'-Exon Start positions from GeneSeqer output with genes from GFF file. 
### USAGE --- perl gsq_gff_extractor.pl -gsq <GeneSeqer_output_file> -gff <GFF_annotation_file> -chrnum <Chromosome_number_of_sequence> -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Cluster-format-file]
### PRE-REQUISITES --- "Getopt::Long" package available at CPAN

my ($infile_gsq,$infile_gth,$infile_gff,$chrnum,$dist,$outfile);

GetOptions (    "gsq=s" => \$infile_gsq,                # GeneSeqer output file. 
		"gth=s" => \$infile_gth,                # GenomeThreader output file.
                "gff=s" => \$infile_gff,               	# GFF format annotation file. Mandatory
		"dist=i" => \$dist,			# Distance from annotated gene start to be used for getting relevant transcripts. Optional
		"chrnum=i" => \$chrnum,			# Chromosome number of sequence input. Mandatory
                "out=s" => \$outfile);             	# Output file name. Optional

my $usage = "\nUSAGE: perl gsq_gff_extractor.pl -gsq [GeneSeqer_output_file] -gth [GenomeThreader_output_file] -gff <GFF_annotation_file> -chrnum <Chromosome_number_of_sequence> -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Cluster-format-file]\n
\t<>-Mandatory parameter\t\t[]-Optional parameter\n
NOTE: Please define all mandatory paramters and install Getopt::Long perl module from CPAN before proceeding\n";

if(defined $infile_gsq || defined $infile_gth) {} else { print $usage; exit; }
if(defined $infile_gff && defined $chrnum) {} else { print $usage; exit; }   ##### Show usage if cmd-line parameters are undefined
if(defined $dist) {} else { $dist=500; }   ####### Default distance to be used for annotation = [-500,+500]

if(defined $infile_gsq) { open (INFILE_TRANS, $infile_gsq) or die ("Could not open $infile_gsq file!!"); }
elsif(defined $infile_gth) { open (INFILE_TRANS, $infile_gth) or die ("Could not open $infile_gth file!!"); }
open (INFILE_GFF, $infile_gff) or die ("Could not open $infile_gff file!!");

####### Setting defaults for cmd-line parameters ################
my (@outname_temp,$outlabel);
if (defined $outfile) 
{ 
 @outname_temp = split (/\//,$outfile);
 $outlabel = $outname_temp[$#outname_temp];
 $outfile = $outfile . ".txt";                  ## Adding ".csv" suffix to filename to use "read.csv" command in R
}
else 
{ 
 @outname_temp = split (/\//,$infile_gsq);
 $outlabel = $outname_temp[$#outname_temp];
 $outfile = $outlabel . ".txt";  
}
open (OUTFILE, "> $outfile") or die ("Could not open $outfile!!");

my ($exon_start,$exon_end,$exon_orientation,$species_chr_id,$gene_orientation,$gene_start,$gene_end,$gene_id,$gene_chr);
my (@exon_positions,@gff_coords,%start_to_id,%id_to_allstarts,%gsq_start_orientation,%gff_start_orientation,%exon_start_end,%gff_start_chr) = ();

while(<INFILE_TRANS>)
{
 if(/^\s*Exon\s+1\s+(\d+)\s+(\d+)\s+\(\s+(\d+)\s+n\);\s+cDNA\s+(\d+)\s+(\d+)\s+\(\s+(\d+)\s+n\);\s+score:\s+(\d*\.\d+)\s*\n$/)       ##### Read in all the 5'-First Exons from the GeneSeqer output
 {
 # print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
  if($1>$2)
  {
   $exon_start = $1;
   $exon_end = $2;
   $exon_orientation = "-";
  }
  elsif ($2>$1)
  {
   $exon_start = $1;
   $exon_end = $2;
   $exon_orientation = "+";
  }
  $gsq_start_orientation{$exon_start}=$exon_orientation; 
  $exon_start_end{$exon_start} = $exon_end;
  push(@exon_positions, $exon_start);
 }
}

while(<INFILE_GFF>)
{
 if(/^[a-zA-z_]*(\d+)\s+([A-Z_a-z0-9]+)\s+(gene)\s+(\d+)\s+(\d+)\s+\.\s+([+-]*).*ID=([A-Za-z0-9_]+);.*$/)  #### GFF file format, check for compatibility 
 {												#### Read in all the gene start, end and IDs from gff file
#  print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
  $species_chr_id = $1;
  if($species_chr_id==$chrnum)
  {
   $gene_chr = $1;
   $gene_orientation = $6;
   $gene_start = $4;
   $gene_end = $5;
   $gene_id = $7;
   push(@gff_coords,$gene_start);
   push(@gff_coords,$gene_end);
   $gff_start_orientation{$gene_start}=$gene_orientation;
   $gff_start_chr{$gene_start} = $gene_chr;
   $start_to_id{$gene_start} = $gene_id;
  }
 }
}

for(my $i=0;$i<=$#exon_positions;$i++)      ####### Associate 5'-transcript ends with gene IDs based on distance from 5'-annotated gene start and orientation
{
 my $exon_min = &min($exon_positions[$i],$exon_start_end{$exon_positions[$i]});
 my $exon_max = &max($exon_positions[$i],$exon_start_end{$exon_positions[$i]});
 for (my $j=0;$j<=$#gff_coords;$j+=2)
 {
  if($gff_start_orientation{$gff_coords[$j]} eq '+')
  {
   if($exon_positions[$i]>($gff_coords[$j]-$dist-1) && $exon_positions[$i]<($gff_coords[$j]+$dist+1))
   {
    print OUTFILE "chr$gff_start_chr{$gff_coords[$j]}\t$start_to_id{$gff_coords[$j]}\t$gff_start_orientation{$gff_coords[$j]}\t$exon_min\t$exon_max\n";
    last;
   }
  }
  elsif($gff_start_orientation{$gff_coords[$j]} eq '-')
  {
   if($exon_positions[$i]>($gff_coords[$j+1]-$dist-1) && $exon_positions[$i]<($gff_coords[$j+1]+$dist+1))
   {
    print OUTFILE "chr$gff_start_chr{$gff_coords[$j]}\t$start_to_id{$gff_coords[$j]}\t$gff_start_orientation{$gff_coords[$j]}\t$exon_min\t$exon_max\n";
    last;
   }
  }
 }
} 

close INFILE_TRANS;
close INFILE_GFF;
close OUTFILE;

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }
