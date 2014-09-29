/*
--------------------------------------------------------------------------------
Copyright (c) 2011-2013, Daniel S. Standage

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
--------------------------------------------------------------------------------
*/

// System libraries
#include <unistd.h>

// Third-party libraries
#include "genometools.h"
#include "sam.h"

// Project libraries
#include "CagepNodeVisitor.h"


// Data structure for program options and arguments
typedef struct
{
  GtStrArray * cageattrs;
  GtStrArray * geneattrs;
  const char * outfile;
  FILE *       outstream;
  unsigned int region;
  int          verbose;
  const char * cagefile;
  const char * genefile;
} cagep_opts;


// Function prototypes
void             cagep_opts_init(cagep_opts *options);
int              cagep_parse_alignments(GtFeatureIndex *genes,
                                        cagep_opts *options, GtError *error);
int              cagep_parse_args(int argc, const char **argv,
                                  cagep_opts *options);
GtFeatureIndex * cagep_parse_genes(cagep_opts *options, GtError *error);
void             print_usage(FILE *outstream);


// Method implementations
void cagep_opts_init(cagep_opts *options)
{
  options->cageattrs = gt_str_array_new();
  options->geneattrs = gt_str_array_new();
  options->outfile   = NULL;
  options->outstream = stdout;
  options->region    = 50000;
  options->verbose   = 0;
  options->cagefile  = NULL;
  options->genefile  = NULL;
}

int cagep_parse_alignments(GtFeatureIndex *genes, cagep_opts *options,
                           GtError *error)
{
  samfile_t *sam = samopen(options->cagefile, "r", NULL);
  if(sam == NULL)
  {
    fprintf(stderr, "[cagep] error: unable to open CAGE file '%s'\n",
            options->cagefile);
    return 1;
  }
  gt_assert(sam->type == 2);

  bam1_t *alignment = bam_init1();
  int bytesread;
  GtRange cagerange, searchrange;
  while((bytesread = samread(sam, alignment)) > 0)
  {
    if(alignment->core.flag & 4)
    {
      // unmapped; ignore
      continue;
    }

    cagerange.start = alignment->core.pos + 1;
    cagerange.end   = bam_calend(&alignment->core, bam1_cigar(alignment)) + 1;
    char *seqid = sam->header->target_name[alignment->core.tid];
    GtStrand strand = GT_STRAND_FORWARD;
    if(alignment->core.flag & 16)
    {
      strand = GT_STRAND_REVERSE;
      searchrange.start = cagerange.start - options->region + 1;
      if(cagerange.start < options->region) searchrange.start = 1;
      searchrange.end = cagerange.end;
    }
    else
    {
      searchrange.start = cagerange.start;
      searchrange.end = cagerange.end + options->region - 1;
    }

    GtArray *downstream_genes = gt_array_new( sizeof(GtFeatureNode *) );
    gt_feature_index_get_features_for_range(genes, downstream_genes, seqid,
                                            &searchrange, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "[cagep] error: %s (seqid='%s')\n", gt_error_get(error),
              seqid);
      return 1;
    }
    if(strand == GT_STRAND_FORWARD)
    {
      gt_array_reverse(downstream_genes);
    }
    while(gt_array_size(downstream_genes) > 0)
    {
      GtFeatureNode *gene = *(GtFeatureNode **)gt_array_pop(downstream_genes);
      GtStrand genestrand = gt_feature_node_get_strand(gene);
      if(genestrand == strand)
      {
        char strc = '+';
        if(strand == GT_STRAND_REVERSE) strc = '-';
        fprintf(options->outstream, "%s\t%s\t%c\t%lu\t%lu", seqid,
                gt_feature_node_get_attribute(gene, "ID"), strc,
                cagerange.start, cagerange.end);
        
        unsigned long i;
        for(i = 0; i < gt_str_array_size(options->geneattrs); i++)
        {
          const char *attrkey   = gt_str_array_get(options->geneattrs, i);
          const char *attrvalue = gt_feature_node_get_attribute(gene, attrkey);
          if(attrvalue == NULL)
            fprintf(options->outstream, "\tgene:%s:undefined", attrkey);
          else
            fprintf(options->outstream, "\tgene:%s:%s", attrkey, attrvalue);
        }
        fputs("\n", options->outstream);
        break;
      }
    }
    gt_array_delete(downstream_genes);
  }
  bam_destroy1(alignment);
  samclose(sam);
  return 0;
}

int cagep_parse_args(int argc, const char **argv, cagep_opts *options)
{
  char c;
  cagep_opts_init(options);
  while((c = getopt(argc, (char **)argv, "a:g:ho:r:v")) >= 0)
  {
    if(c == 'a')
    {
      gt_str_array_add_cstr(options->cageattrs, optarg);
    }
    if(c == 'g')
    {
      gt_str_array_add_cstr(options->geneattrs, optarg);
    }
    else if(c == 'h')
    {
      print_usage(stdout);
      return -1;
    }
    else if(c == 'o')
    {
      options->outfile = optarg;
      options->outstream = fopen(optarg, "w");
      if(options->outstream == NULL)
      {
        fprintf(stderr, "[cagep] error: could not open output file '%s'\n",
                optarg);
        return 1;
      }
    }
    else if(c == 'r')
    {
      options->region = atoi(optarg);
    }
    else if(c == 'v')
    {
      options->verbose += 1;
    }
    else
    {
      fprintf(stderr, "[cagep] error: unknown option '%c'\n", c);
      return 1;
    }
  }

  int numfiles = argc - optind;
  if(numfiles != 2)
  {
    fprintf(stderr, "[cagep] error: must provide two input files\n");
    print_usage(stderr);
    return 1;
  }
  options->cagefile = argv[optind + 0];
  options->genefile = argv[optind + 1];

  return 0;
}

GtFeatureIndex *cagep_parse_genes(cagep_opts *options, GtError *error)
{
  GtNodeStream *gff3 = gt_gff3_in_stream_new_unsorted(1, &options->genefile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3);

  GtFeatureIndex *genes = gt_feature_index_memory_new();
  unsigned int ngenes = 0, nseqs = 0;
  GtNodeVisitor *nv = cagep_node_visitor_new(genes, &ngenes, &nseqs);
  GtGenomeNode *gn;
  bool loaderror;
  while(!(loaderror = gt_node_stream_next(gff3, &gn, error)) && gn)
  {
    gt_genome_node_accept(gn, nv, error);
  }
  gt_node_stream_delete(gff3);
  if(options->verbose)
  {
    fprintf(stderr, "[cagep] loaded %u genes (from %u sequences) into memory\n",
            ngenes, nseqs);
  }
  gt_node_visitor_delete(nv);

  return genes;
}

void print_usage(FILE *outstream)
{
  fputs("\ncagep: associate CAGE alignments with genes\n"
"Usage: cagep [options] cage-aligns.sam genes.gff3\n"
"  Options:\n"
"    -a STRING    attribute to extract from each CAGE alignment in the SAM\n"
"                 file; this option can be used multiple times to extract\n"
"                 multiple attribute values\n"
"    -g STRING    attribute to extract from each gene in the GFF3 file; this\n"
"                 option can be used multiple times to extract multiple\n"
"                 attribute values\n"
"    -h           print this help message and exit\n"
"    -o FILE      print output to file; default is the terminal (stdout)\n"
"    -r INT       how many nucleotides downstream of each CAGE alignment to\n"
"                 look for associated genes; default is 50000\n"
"    -v           print warnings in addition to errors\n\n",
        outstream);
}

int main(int argc, const char **argv)
{
  cagep_opts options;
  gt_lib_init();
  int had_parse_error = cagep_parse_args(argc, argv, &options);
  if(had_parse_error < 0) return 0;
  if(had_parse_error > 0) return 1;

  GtError *error = gt_error_new();

  GtFeatureIndex *genes = cagep_parse_genes(&options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[cagep] error parsing GFF3 file '%s': %s\n",
            options.genefile, gt_error_get(error));
    return 1;
  }

  had_parse_error = cagep_parse_alignments(genes, &options, error);
  if(had_parse_error) return 1;

  gt_feature_index_delete(genes);

  fclose(options.outstream);
  gt_str_array_delete(options.cageattrs);
  gt_str_array_delete(options.geneattrs);
  gt_error_delete(error);
  gt_lib_clean();
  return 0;
}
