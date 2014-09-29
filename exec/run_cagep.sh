for i in ../data/sam_files/*.sam
do
echo $i
./cagep $i ../data/knownGene.gff3 > $i.cagep
done
