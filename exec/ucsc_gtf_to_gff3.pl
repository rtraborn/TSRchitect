#!/usr/bin/env perl
use strict;

my $ignored_lines = {};
my $transcripts = {};
my $transcript_features = {};

while(my $line = <STDIN>)
{
  my @fields = split(/\t/, $line);
  
  # Only process exons and CDS segments
  unless($fields[2] eq 'exon' or $fields[2] eq 'CDS')
  {
    $ignored_lines->{$fields[2]} = 0 unless(defined($ignored_lines->{$fields[2]}));
    $ignored_lines->{$fields[2]} += 1;
    next;
  }
  
  # Grab transcript ID
  my($id) = $fields[8] =~ m/transcript_id "(.+?)";/;
  
  # If transcript has been seen before, make sure the sequence and strand match up
  if(defined($transcripts->{$id}))
  {
    if($fields[0] ne $transcripts->{$id}->[0] or $fields[6] ne $transcripts->{$id}->[6])
    {
      printf(STDERR "Error: transcript '%s' has features on different sequences and/or strands\n", $id);
      die();
    }
  }
  else # If not, store the transcript data
  {
    $transcripts->{$id} = [$fields[0], $fields[1], "mRNA", $fields[3], $fields[4], ".", $fields[6], "."];
    $transcript_features->{$id} = [];
  }
  $transcripts->{$id}->[3] = $fields[3] if($fields[3] < $transcripts->{$id}->[3]);
  $transcripts->{$id}->[4] = $fields[4] if($fields[4] > $transcripts->{$id}->[4]);
  
  # Double check for features with negative length
  if($fields[3] > $fields[4])
  {
    printf(STDERR "Feature start > end, skipping: %s", $line);
  }
  else
  {
    push(@{$transcript_features->{$id}}, [@fields]);
  }
}

# Print out features
foreach my $id(keys(%$transcripts))
{
  if($transcripts->{$id}->[3] > $transcripts->{$id}->[4])
  {
    printf(STDERR "Transcript '%s' start > end, skipping: %s\n", $id);
  }
  else
  {
    printf("%s\t%s\tgene\t%lu\t%lu\t%s\t%s\t%s\t%s\n", @{$transcripts->{$id}}[0, 1], @{$transcripts->{$id}}[3..7], "ID=$id.gene");
    printf("%s\t%s\t%s\t%lu\t%lu\t%s\t%s\t%s\t%s\n", @{$transcripts->{$id}}, "ID=$id.transcript;Parent=$id.gene");
    foreach my $feature(@{$transcript_features->{$id}})
    {
      $feature->[8] = "Parent=$id.transcript";
      printf("%s\t%s\t%s\t%lu\t%lu\t%s\t%s\t%s\t%s\n", @$feature);
    }
  }
}