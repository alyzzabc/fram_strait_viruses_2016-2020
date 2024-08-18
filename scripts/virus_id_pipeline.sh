#Assembly of PacBio HiFi long reads, compressed into .gz files
#sample_names is a file containing sample names (or base names of filenames)
for i in `cat sample_names`
do
	hifiasm_meta -t 32 -o $i.asm $i.gz 2> $i.asm-meta.log
done

#After assembly, extract sequences from primary contigs file
for i in `cat sample_names`
do
	awk '/^S/{print ">"$2;print $3}' $i.asm.p_ctg.gfa > $i.p_ctg.fa
done

#seqtk, extract only sequences greater than 1500 bp
for i in *.asm.p_ctg.fa
do
	seqtk seq -L 1500 $i > 1.5k.$i
done

#VirSorter2, identify viruses
for i in `cat sample_names`
do
	virsorter run -w $i.vs2.out -i $i --include-groups 'dsDNAphage, ssDNA, NCLDV, RNA, lavidaviridae' --min-length 1500 -j 4 --min-score 0.5 all
done

#CheckV, assess quality of viral contigs
for i in `cat sample_names`
do
	mkdir $i.checkv.out
 	checkv end_to_end $i.vs2.out/final-viral-combined.fa $i.checkv.out -t 28 -d /path/to/checkv/database
done

#Concatenate sequences of viruses and proviruses
for i in `cat sample_names`
do
	cat $i.checkv.out/proviruses.fna $i.checkv.out/viruses.fna > $i.provirus.virus.fa
done

#Merge "vs2-pass1/final-viral-score.tsv" and "checkv/contamination.tsv" to screen for high confidence viral contigs
for i in `cat sample_names`
do
	paste $i.vs2.out/final-viral-score.tsv $i.checkv.out/contamination.tsv | cut --complement -f13 > $i.vs2-checkv.tsv
done

#Run screening.py to assign confidence of viral contigs
python screening.py

#Extract sequences of Keep1 contigs, i.e. high confidence
for i in `cat sample_names`
do
grep "keep1" $i"_screened.tsv" | cut -f1 > $i.keep1.seqIDs
done

for i in `cat sample_names`
do
grep -A1 -h -f $i.keep1.seqIDs $i.provirus.virus.fa > $i.keep1.fa
done

#DeepMicroClass, identify non-viral Keep1 contigs
for i in `cat sample_names`
do
	mkdir -p $i.dmc.out
	DeepMicroClass predict -i $i.keep1.fa -o $i.dmc.out -d cuda
done

#Remove keep1 contigs that were identified by DeepMicroClass as prokaryotic or eukaryotic
#Combine all high confidence viral contigs in one file, e.g. all_keep1_high_conf.fa

#CD-HIT, cluster viral contigs into vOTUs
cd-hit-est -M 100000 -c 0.95 -d 100 -g 1 -aS 0.85 -t 16 -i all_keep1_high_conf.fa  -o cd95_all_keep1_high_conf.fa

#Εxtract only vOTUs greater than 10kb
seqtk seq -L 10000 cd95_all_keep1_high_conf.fa > 10k.cd95_all_keep1_high_conf.fa
#VPF-Class, assign taxonomy
for i in `cat sample_names`
do
	singularity exec --bind /path/to/vpf-class-data:/opt/vpf-tools/vpf-data \
                 --bind /path/to/10k.cd95_all_keep1_high_conf.fa:/opt/vpf-tools/input-sequences/10k.cd95_all_keep1_high_conf.fa \
                 --bind /path/to/output:/opt/vpf-tools/outputs \
                 /path/to/vpf-tools_latest.sif \
                 vpf-class -i path/to/10k.cd95_all_keep1_high_conf.fa -o /opt/vpf-tools/outputs/test-classified
done

#VIBRANT, determine lifestyle
for i in `cat sample_names`
do
	VIBRANT_run.py -i 10k.cd95_all_keep1_high_conf.fa -folder /path/to/output/directory
done

#iPHoP, predict hosts
for i in `cat sample_names`
do
	iphop predict --fa_file 10k.cd95_all_keep1_high_conf.fa --db_dir /path/to/iPHoP/database --out_dir /path/to/output/directory
done

#Prodigal, predict proteins
prodigal -i 10k.cd95_all_keep1_high_conf.fa -o 10k.cd95_all_keep1_high_conf.genes -a 10k.cd95_all_keep1_high_conf.faa -p meta

#To calculate abundance, start by mapping reads (.fastq.gz files) to contigs using minimap2
for i in `cat sample_names`
do
	minimap2 -ax asm5 --sam-hit-only 10k.cd95_all_keep1_high_conf.fa $i.fastq.gz > $i.sam
done

#Convert .sam files to .bam files
for i in `cat sample_names`
do
	samtools view -bS $i.sam > $i.bam
done

#Sort .bam file
for i in `cat sample_names`
do
	samtools sort $i.bam -o $i.sorted.bam
done

#Get only primary alignments
for i in `cat sample_names`
do
	samtools view -b -F 0x800 -F 0x100 $i.sorted.bam > $i.sorted.primary.bam
done

#Run calc_cvpg.sh to actually calculate abundance. calc_cvpg.sh and calc_coverage_statistics.py should be in the same directory as the .sorted.primary.bam files
bash calc_cvpg.sh -i /path/to/.sorted.primary.bam/files -t 10



