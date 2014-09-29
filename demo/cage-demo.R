###Note: in progress 2/20/2014
### cagep and TSRchitect use cases are combined here, because that's the way I want the package's functionilty to work
## waiting for Linux segmentation error to be fixed before proceeding

_____________ TSRchitect v1.0 __________________

CAGE Analysis Use Case

Processing mapped CAGE data using cagep

>$ cd TSRchitect/data

Let's use the example file 'Cage_Haoec_Sample.sam' that was included in /data

How many lines does it have? Let's enter the following command to find out:

>$ wc -l Cage_Haoec_Sample.sam   
1000000 Cage_Haoec_Sample.sam

Now that we have the CAGE file, we need to assign each TSS to an annotated coding gene. 
Do do this, we'll need 1) the current gene annotation files for the appropriate genome (in this case the human genome hg19) and 2) a script that does the assignment process quickly.

For 1) We've included the current annotation for the human genome (hg19), 'knownGene.gff3'
and 2) We have included the script 'cagep', which we'll need to have installed before we can proceed.

After cagep has been successfully compiled and installed, enter the following commands to run cagep:

###To do: ######
Update the cagep commands necessary to create the .cagep file



