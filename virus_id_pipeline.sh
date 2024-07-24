ssh smomw535@nesh-fe.rz.uni-kiel.de

##Convert multple fastq files to fasta using seqtk:
# seqtk seq -a file.fastq > file.fa

for i in *.fastq
do
sbatch --job-name=seqtk --cpus-per-task=1 --mem=4000 --time=4:00:00 --error=err --partition=cluster --wrap="seqtk seq -a $i > /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/$i.fa"
done

##Split multiple fasta file into several fasta files using genometools splitfasta. VirSorter2 takes a long time so splitting the fasta files can make it run faster.
# gt splitfasta -numfiles 5 fasta.fa

for i in *.fa
do
	mkdir ./$i.split
	cp $i ./$i.split
	cd ./$i.split
	gt splitfasta -numfiles 5 $i
	cp *.fa.* /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta
	cd ..
done

##VirSorter
# virsorter run -w test.out -i test.fa --min-length 1500 -j 4 all

conda activate vs2

#dsDNAphages and ssDNA only, which is the default

for i in /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/*
do
sbatch --job-name=vs2 --cpus-per-task=30 --mem=4000 --time=48:00:00 --error=err --partition=cluster --wrap="virsorter run -w $i.out -i $i --min-length 1500 -j 4 all"
done

#to include all groups
for i in /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/*.fa.*
do
sbatch --job-name=vs2 --cpus-per-task=30 --mem=4000 --time=48:00:00 --error=err --partition=cluster --wrap="virsorter run -w $i.all.out -i $i --include-groups 'dsDNAphage, ssDNA, NCLDV, RNA, lavidaviridae' --min-length 1500 -j 4 --min-score 0.5 all"
done


##Combine VirSorter2 result final-viral-combined.fa of each split fasta file using combine-fasta-files.py

##CheckV

conda activate checkv

for i in *.combined.fa
do
	mkdir ./$i.checkv
 	sbatch --job-name=checkv --cpus-per-task=30 --mem=4000 --time=48:00:00 --error=err --partition=cluster --wrap="checkv end_to_end /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/$i ./$i.checkv -t 28 -d /gxfs_home/geomar/smomw535/checkv-db-v1.4"
done

##Run VirSorter again
conda activate vs2

#Concatenate provirus and virus fasta files
for i in *.checkv
do
	cd $i
	cat proviruses.fna viruses.fna > /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus/$i.combined.fna
	cd ..
done

#Run VirSorter2 on concatenated fasta files in preparation for DRAM-v
for i in *.combined.fna
do
	sbatch --job-name=checkv-pass2 --cpus-per-task=30 --mem=4000 --time=48:00:00 --error=err --partition=cluster --wrap="virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i $i -w $i.vs2-pass2-all --include-groups 'dsDNAphage, ssDNA, NCLDV, RNA, lavidaviridae' --min-length 5000 --min-score 0.5 -j 28 all"
done

###Run prodigal then hmmscan, while I can't get DRAM-v to work.

#prodigal
conda activate prodigal

for i in *.vs2-pass2-all
do
	cd $i
	sbatch --job-name=prodigal --cpus-per-task=8 --mem=4000 --time=4:00:00 --error=err --partition=cluster --wrap="prodigal -i final-viral-combined.fa -o $i.genes -a $i.proteins.faa"
	cd ../
done

#hmmscan
conda activate hmmer

for i in *.vs2-pass2-all
do
	cd $i
	sbatch --job-name=hmmscan --cpus-per-task=10 --mem=8000 --time=48:00:00 --error=err --partition=cluster --wrap="hmmscan --tblout $i.a.tblout --cpu 1 $WORK/databases/pfam/Pfam-A.hmm *proteins.faa"
	cd ../
done

##DRAM-v

# step 1 annotate

#DRAM-v.py annotate -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab -o dramv-annotate --skip_trnascan --threads 28 --min_contig_size 1000

conda activate DRAM-py3.10

#keep a low number of CPUs/threads and a large memory
for i in ./*.vs2-pass2-all
do
	cd $i/for-dramv
	sbatch --job-name=dramv-annotate --cpus-per-task=10 --mem=192000 --time=48:00:00 --error=err --partition=cluster --wrap="DRAM-v.py annotate 
	\ -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab -o dramv-annotate --skip_trnascan --min_contig_size 1000"
	cd ../../
done

#step 2 summarize annotations
for i in ./*.vs2-pass2-all
do
	cd $i/for-dramv/dramv-annotate
	sbatch --job-name=dramv-distill --cpus-per-task=10 --mem=64000 --time=12:00:00 --error=err --partition=cluster --wrap="DRAM-v.py distill -i annotations.tsv 
	\ -o dramv-distill"
	cd ../../../
done

## Screening based on viral and host gene counts, score, hallmark gene counts, and contig length
# Merge "vs2-pass1/final-viral-score.tsv" and "checkv/contamination.tsv" using vs2-check-merge-tables.py
# Assign Keep1, Keep2, Manual, Discard using screening.py

## Manual curation
# See annotations in "dramv-annotate/annotations.tsv"


#Make a list of all Keep 1 contigs

#Get sequences of all Keep 1 contigs
grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/all.keep1s *combined.fna > all.keep1s.vs2.checkv.fa


cat /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/all.keep1s | wc -l
613809

#Make a fasta file of all contigs from VS2-pass1.
#cat all.vs2.pass1.contigs.fa | wc -l
#2691412
#cat all.vs2.pass1.contigs.fa | grep ">" | wc -l
#1345706
#Make fasta file of all keep1 contigs from all.vs2.pass1.contigs.fa
#grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/all.keep1s all.vs2.pass1.contigs.fa > all.keep1_vs2.pass1.fa

grep -c ">" all.keep1s.fa 
613809

#Make a fasta file of all contigs from checkv.combined.fna
#(base) cat all_reads.fa | wc -l
#2691412
#(base) cat all_reads.fa | grep ">" | wc -l
#1345706
#Make fasta file of all keep1 contigs from all.checkv.combined.reads.fa
#grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/all.keep1s all.checkv.combined.reads.fa > all.keep1_checkv.combined.fa

#(base) cat all.keep1_checkv.combined.fa | grep ">" | wc -l
#616995

#Make a fasta file of all proteins
#cat *.vs2-pass2-all/*.proteins.faa > all.aa.fa
#cat all.aa.fa | grep ">" | wc -l
#8045059
#Make fasta file of all keep1 proteins from all.aa.fa
#grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/all.keep1s all.aa.fa > all.keep1.aa.fa
#(base) cat all.keep1.aa.fa | grep ">" | wc -l
#5272397

#CD-HIT
#95% identity
sbatch --job-name=cdhit-95 --cpus-per-task=10 --mem=1000 --time=12:00:00 --error=err --partition=cluster --wrap="cd-hit -i all.keep1.aa.fa -o all.keep1.aa.95 -c 0.95 -n 5"

#99% identity
sbatch --job-name=cdhit-99 --cpus-per-task=10 --mem=1000 --time=12:00:00 --error=err-99 --partition=cluster --wrap="cd-hit -i all.keep1.aa.fa -o all.keep1.aa.99 -c 0.99 -n 5"

(cdhit) cat all.keep1.aa.95 | grep ">" | wc -l
378288
(cdhit) cat all.keep1.aa.99 | grep ">" | wc -l
431254

#get all NCLDV
cat final_df.tsv | grep NCLDV | cut -d , -f 2 > all.ncldv
cat final_df.tsv | grep NCLDV | cut -d , -f 2 | wc -l
475405

cat final_df.tsv | grep NCLDV | grep keep1 | cut -d , -f 2 > all.ncldv.keep1s
cat final_df.tsv | grep NCLDV | grep keep1 | cut -d , -f 2 | wc -l
51194

grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/vs2_checkv_tables/all.ncldv.keep1s all.aa.fa | grep ">" | wc -l

#vContact2

#Run prodigal -meta
conda activate prodigal

for i in *.vs2-pass2-all
do
	cd $i
	sbatch --job-name=prodigal-meta --cpus-per-task=8 --mem=4000 --time=4:00:00 --error=err --partition=cluster --wrap="prodigal -i final-viral-combined.fa -o $i.genes -a $i.meta.faa"
	cd ../
done

#Make a fasta file of all proteins from meta
cat *.vs2-pass2-all/*.meta.faa > all.aa.meta.fa
cat all.aa.meta.fa | grep ">" | wc -l
8045059
#Make fasta file of all keep1 proteins from all.aa.fa
grep -A 1 -h -f /$WORK/arctic_metag/all_data/all_data_fasta/splitfasta/all.keep1s all.aa.fa > all.keep1.aa.meta.fa
cat all.keep1.aa.meta.fa | grep ">" | wc -l
5272397


#Generate gene to genome mapping file from Prodigal output .faa
vcontact2_gene2genome -p all.keep1.aa.meta.fa -o all.keep1.meta.g2g.csv -s 'Prodigal-FAA'

#Run vContact2
sbatch --job-name=vcontact2_resume --cpus-per-task=10 --mem=16000 --time=48:00:00 --error=err-vc2 --partition=cluster --wrap="vcontact2 --raw-proteins all.keep1.aa.meta.fa --rel-mode 'Diamond' --proteins-fp all.keep1.meta.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/cluster_one-1.0.jar --output-dir vcontact2.out"
#INFO:vcontact2: Building the cluster and profiles (this may take some time...)
#If it fails, try re-running using --blast-fp flag and specifying merged.self-diamond.tab (or merged.self-blastp.tab)
#INFO:vcontact2: Saving intermediate files...
#INFO:vcontact2: Read 5276659 entries (dropped 155142 singletons) from vcontact2.out/vConTACT_profiles.csv
#/var/spool/slurmd/job6308971/slurm_script: line 4: 1153495 Killed                  vcontact2 --raw-proteins all.keep1.aa.meta.fa --rel-mode 'Diamond' --proteins-fp all.keep1.meta.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/cluster_one-1.0.jar --output-dir vcontact2.out
#slurmstepd: error: Detected 1 oom-kill event(s) in step 6308971.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.

#Rerun the whole thing with larger memory?
sbatch --job-name=vcontact2 --cpus-per-task=10 --mem=192000 --time=48:00:00 --error=err-vc2-again --partition=cluster --wrap="vcontact2 --raw-proteins all.keep1.aa.meta.fa --rel-mode 'Diamond' --proteins-fp all.keep1.meta.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/cluster_one-1.0.jar --output-dir vcontact2.rerun.out"

#Use output c1.ntw to visualize network on Cytoscape with annotations from output genome_by_genome_overview.csv

#While vContact2 is running, download iphop database.
sbatch -p data -t 24:00:00 --mem=100000 -J iphob-db -o iphop-db.out -e iphop-db.err -n 1 --wrap="echo yes | iphop download --db_dir iphop_db/ -dbv iPHoP_db_for-test"

##CD-HIT: Cluster the gene sequences (not the aa sequences)
sbatch -p cluster -t 12:00:00 --mem=192000 -J cdhit-est -o cdhit-est.out -e cdhit-est.err -n 30 --wrap="cd-hit-est -M 100000 -c 0.99 -i /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/all.keep1_vs2.pass1.fa -o cdhit.99.all.keep1_vs2.pass1.fa"
sbatch -p cluster -t 12:00:00 --mem=192000 -J cdhit95-est -o cdhit95-est.out -e cdhit95-est.err -n 30 --wrap="cd-hit-est -M 100000 -c 0.95 -i /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/all.keep1_vs2.pass1.fa -o cdhit.95.all.keep1_vs2.pass1.fa"

sbatch -p cluster -t 12:00:00 --mem=192000 -J cdhit99-est -o cdhit-est99.out -e cdhit-est.err -n 30 --wrap="cd-hit-est -M 100000 -T 30 -c 0.99 -i /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/all.keep1_vs2.pass1.fa -o cdhit.99.rep.all.keep1_vs2.pass1.fa"
sbatch -p cluster -t 12:00:00 --mem=192000 -J cdhit95-rep -o cdhit95-rep.out -e cdhit95-rep.err -n 30 --wrap="cd-hit-est -M 100000 -T 30 -c 0.95 -i /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/all.keep1_vs2.pass1.fa -o cdhit.95.rep.all.keep1_vs2.pass1.fa"


##VContact2, run on clustered proteins
#Generate gene to genome mapping file from Prodigal output .faa
vcontact2_gene2genome -p all.keep1.aa.95 -o all.keep1.aa.95.csv -s 'Prodigal-FAA'
vcontact2_gene2genome -p all.keep1.aa.99 -o all.keep1.aa.99.csv -s 'Prodigal-FAA'

#Actually run vContact2, on both 0.95 and 0.99 identity thresholds
sbatch --job-name=vc2-95 --cpus-per-task=20 --mem=512000 --time=48:00:00 --error=err-vc2-95 --partition=cluster --wrap="vcontact2 --raw-proteins all.keep1.aa.95 --rel-mode 'Diamond' --proteins-fp all.keep1.aa.95.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/cluster_one-1.0.jar --output-dir vc2.aa.95.out"
sbatch --job-name=vc2-99 --cpus-per-task=20 --mem=512000 --time=48:00:00 --error=err-vc2-99 --partition=cluster --wrap="vcontact2 --raw-proteins all.keep1.aa.99 --rel-mode 'Diamond' --proteins-fp all.keep1.aa.99.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/cluster_one-1.0.jar --output-dir vc2.aa.99.out"



sbatch --job-name=iphop-test --cpus-per-task=10 --mem=16000 --time=24:00:00 --error=err --partition=cluster --wrap="iphop predict --fa_file test_input_phages.fna --db_dir Test_db/ --out_dir iphop_test_results/test_input_phages_iphop"


##Check vContact2 out.out4 after adding java path to PATH
export PATH=/usr/bin/java:$PATH



sbatch --job-name=vc2 --cpus-per-task=20 --mem=16000 --time=48:00:00 --partition=cluster --wrap="vcontact2 --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/miniconda3/bin/cluster_one-1.0.jar --output-dir vc2.out -vvv"

sbatch --job-name=vc2 --cpus-per-task=20 --mem=16000 --time=48:00:00 --partition=cluster --wrap="vcontact2 --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode MCL --output-dir vc2.out1 --vvv"


##vContact2 in smomw535
conda activate vContact2
#/bin/sh: java: command not found
#ERROR:vcontact2: Error in contig clustering
#ERROR:vcontact2: No columns to parse from file
##even --vcs-mode MCL gives an error in contig clustering, so ClusterONE is probably not the problem


#where python = 3.7
conda activate vContact2.3
sbatch --job-name=vc2-arch --cpus-per-task=20 --mem=16000 --time=48:00:00 --partition=cluster --wrap="vcontact2 --raw-proteins test_data/VIRSorter_genome.faa --rel-mode 'Diamond' --proteins-fp test_data/VIRSorter_genome_g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/miniconda3/bin/cluster_one-1.0.jar --output-dir vc2.arch.out -v"
#it finished omg, maybe because the env has python 3.7? wait why did it fail???

#this time, try the database we want
sbatch --job-name=vc2-prok --cpus-per-task=20 --mem=16000 --time=48:00:00 --partition=cluster --wrap="vcontact2 --raw-proteins test_data/VIRSorter_genome.faa --rel-mode 'Diamond' --proteins-fp test_data/VIRSorter_genome_g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gxfs_home/geomar/smomw535/miniconda3/bin/cluster_one-1.0.jar --output-dir vc2.prok.out -v"

##vContact2 in sulky229, where python = 3.7
conda activate vContact2.4

#--rel-mode 'Diamond' --db 'ArchaeaViralRefSeq211-Merged'
sbatch --job-name=vc2 --cpus-per-task=20 --mem=16000 --time=12:00:00 --partition=cluster --wrap="vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir vc2.out24 -v"
#ERROR:vcontact2: Error in contig clustering
#ERROR:vcontact2: No columns to parse from file

#--rel-mode 'BLASTP'
vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'BLASTP' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir blastp.cli.out -vvv
sbatch --job-name=vc2 --cpus-per-task=20 --mem=16000 --time=12:00:00 --partition=cluster --wrap="vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'BLASTP' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir blastp.slurm.out -vvv"
##this works in command line, not in slurm

#--rel-mode 'Diamond'
vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.cli.out -vvv
sbatch --job-name=vc2 --cpus-per-task=20 --mem=16000 --time=12:00:00 --partition=cluster --wrap="vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.slurm.out -vvv"
#this works in command line, not in slurm

salloc --time=6:00:00 -n 32 -p cluster --mem=100000
vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'BLASTP' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir blastp.salloc.out -vvv
#this works

#how about try cli with the database we want?
vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.cli.out -vvv
vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'BLASTP' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir blastp.prok.cli.out -vvv
#both gives this error:
#ValueError: invalid mode: 'rU'

#did this again but it's working so far I AM SO CONFUSEDT
nohup vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.cli.out -vvv &


--db {None,ProkaryoticViralRefSeq85-ICTV,ProkaryoticViralRefSeq85-Merged,ProkaryoticViralRefSeq88-Merged,ProkaryoticViralRefSeq94-Merged,ProkaryoticViralRefSeq97-Merged,
ProkaryoticViralRefSeq99-Merged,ProkaryoticViralRefSeq201-Merged,ProkaryoticViralRefSeq207-Merged,ProkaryoticViralRefSeq211-Merged,ArchaeaViralRefSeq85-Merged,
ArchaeaViralRefSeq94-Merged,ArchaeaViralRefSeq97-Merged,ArchaeaViralRefSeq99-Merged,ArchaeaViralRefSeq201-Merged,ArchaeaViralRefSeq207-Merged,ArchaeaViralRefSeq211-Merged}

##try if nohup will work
nohup vcontact2 --c1-bin $HOME/conda/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.arch.nohup.out -v > nohup.out &

##try if nohup will work on smomw535
nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ArchaeaViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.arch.nohup.out -v > nohup.out &
#works!

##try if ProkaryoticViralRefSeq211-Merged will work in smomw535
salloc --time=6:00:00 -n 32 -p cluster --mem=100000
nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.nohup.out -v > nohup1.out &
##ok OMG THIS WORKED
#maybe follow with disown next time?
salloc --time=6:00:00 -n 32 -p cluster --mem=100000
conda activate vcontact2.3
nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.nohup.disown.out -v > test.out &
disown
#[1] 2390207
#allocation was relinquished huhu
#again, no disown
salloc --time=6:00:00 -n 32 -p cluster --mem=100000
conda activate vcontact2.3
nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.nohup.disown.out -v > test.out &
#[1] 2394841
#relinquished



#again, no salloc
nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.nohup.disown.out -v > test.out &
#[1] 2401535
#ok this is still working, it worked. ok so my nohup works, whew.


nohup vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir diamond.prok.test1.out -v > test1.out &
#[1] 2782822

#Now write a script.
#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --time=01:00:00
#SBATCH --job-name=vc2-test
#SBATCH --mem=16000 
#SBATCH --partition=cluster

# Activate vCotact environment
source activate vcontact2.3

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load openjdk/11.0.2

# And finally run the job​
srun vcontact2 --c1-bin $HOME/miniconda3/bin/cluster_one-1.0.jar --raw-proteins VIRSorter_genome_1-5.faa --rel-mode 'Diamond' --proteins-fp VIRSorter_genome_1-5.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --output-dir test3.out -v


##Run prodigal on cd-hit clustered sequences
conda activate prodigal
sbatch --job-name=prodigal --cpus-per-task=8 --mem=4000 --time=4:00:00 --partition=cluster --wrap="prodigal -i cdhit.95.all.keep1_vs2.pass1.fa -o cdhit.95.all.keep1_vs2.pass1.fa.genes -a cdhit.95.all.keep1_vs2.pass1.fa.proteins.faa -p meta"
sbatch --job-name=prodigal --cpus-per-task=8 --mem=4000 --time=4:00:00 --partition=cluster --wrap="prodigal -i cdhit.99.all.keep1_vs2.pass1.fa -o cdhit.99.all.keep1_vs2.pass1.fa.genes -a cdhit.99.all.keep1_vs2.pass1.fa.proteins.faa -p meta"

#Then generate gene to genome mapping files
conda activate vcontact2.3
vcontact2_gene2genome -p cdhit.95.all.keep1_vs2.pass1.fa.proteins.faa -o cdhit.95.all.keep1_vs2.pass1.fa.proteins.faa.meta.g2g.csv -s 'Prodigal-FAA'
vcontact2_gene2genome -p cdhit.99.all.keep1_vs2.pass1.fa.proteins.faa -o cdhit.99.all.keep1_vs2.pass1.fa.proteins.faa.meta.g2g.csv -s 'Prodigal-FAA'

#Create the bash scripts for running vcontact
touch vc2.cdhit95.sh
touch vc2.cdhit99.sh

#sbatch the bash scripts
sbatch vc2.cdhit95.sh
sbatch vc2.cdhit99.sh

#Run mummer
sbatch --job-name=nucmer --cpus-per-task=10 --mem=4000 --time=48:00:00 --partition=cluster --wrap="/gxfs_work1/geomar/smomw535/stampede-clustergenomes/bin/Cluster_genomes.pl -f all.keep1_vs2.pass1.fa -c 80 -i 95"

#Create blast database

#use megablast
sbatch --job-name=blastn --cpus-per-task=32 --mem=4000 --time=24:00:00 --partition=cluster --wrap="blastn -query /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus/all.keep1_vs2.pass1.fa -db /gxfs_work1/geomar/smomw535/databases/blastdb/arctic_viral_db -outfmt '6 std qlen slen' -max_target_seqs 10000 -o /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus/all.keep1_vs2.pass1.blast.tsv -num_threads 32"


##Run vContact again on unclustered keep1s
conda activate prodigal
sbatch --job-name=prodigal-meta --cpus-per-task=8 --mem=4000 --time=4:00:00 --error=err --partition=cluster --wrap="prodigal -i all.keep1_vs2.pass1.fa -o all.keep1_vs2.pass1.genes -a all.keep1_vs2.pass1.faa"

conda activate vcontact2.3
vcontact2_gene2genome -p all.keep1_vs2.pass1.faa -o all.keep1_vs2.pass1.faa.g2g.csv -s 'Prodigal-FAA'
touch vc2.unclustered.sh
sbatch vc2.unclustered.sh

##Using CHERRY for host prediction

export PATH=$WORK/databases/release207_v2/
sbatch --job-name=cherry --cpus-per-task=10 --mem=16000 --time=48:00:00 --partition=cluster --wrap="python run_Speed_up.py --contigs test_contigs.fa --len 8000 --model pretrain --topk 1"

python run_Speed_up.py --contigs test_keep1.fa --len 8000 --model pretrain --topk 1

#CHERRY works on home node, not on slurm

## Retrying iphop again — idk what I did with the LD_LIBRARY_PATH but for sure iphop_env works. Redownloaded each module separately (bitbucket to local then scp to server) because git clone or wget do not download the model completely.
## David said I can use data partition next time

#Download entire iphop db
sbatch --job-name=iphop-db --cpus-per-task=10 --mem=16000 --time=24:00:00 --partition=data --wrap="echo yes | iphop download --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/"

#All done ! The database has been put in /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub === This is the path you need to give to iPHoP (using the -d argument)
#If you need to save space, you should be able to delete /gxfs_work1/geomar/smomw535/databases/iphop_db/iPHoP_db_Sept21.tar.gz (although we would recommend doing so only after you checked that the database download indeed worked as expected)

#Now run iphop on real data
sbatch --job-name=iphop --cpus-per-task=10 --mem=32000 --time=48:00:00 --partition=cluster --wrap="iphop predict --fa_file new.all.keep1_vs2.pass1.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir iphop_results/"

#Remove "--" in this fasta file
awk '{sub(/\--$/,"")}1' all.keep1s.vs2.checkv.fa > new.all.keep1s.vs2.checkv.fa all.keep1_vs2.pass1.fa > new.all.keep1_vs2.pass1.fa

#Looks like iphop will take more than 48h to complete, so need to run it again with qos=long
sbatch --job-name=iphop-long --cpus-per-task=10 --mem=32000 --time=120:00:00 --qos=long --partition=cluster --wrap="iphop predict --fa_file new.all.keep1_vs2.pass1.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir iphop_results/"
#still running but i accidentally deleted slurm standard output ugh

#vContact issues
vc2.unclustered.latest.out/ failed due to out-of-memory. Reran this under job name vc2-latest-continued
vc2.unclustered.again.out/ still running as job name vc2-latest but will probably fail due to out-of-memory

#Try an entire vcontact2 run, just in case! Script in vc2.unclustered.latest.sh
##Job ID: 7170755

#Try a new iphop run, just in case! And make this verbose --debug!!
sbatch --job-name=iphop-final --cpus-per-task=10 --mem=128000 --time=120:00:00 --qos=long --partition=cluster --wrap="iphop predict --fa_file new.all.keep1_vs2.pass1.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir iphop_final_results/ --debug --num_threads 10"
##Job ID: 7170770

#Split multifasta file (~617k sequences) into 10 (~53k-71k sequences per file)
gt splitfasta -numfiles 10 new.all.keep1_vs2.pass1.fa

#Run iphop again
for i in *fa
do
	sbatch --job-name=iphop-split --cpus-per-task=10 --mem=64000 --qos=long --time=120:00:00 --partition=cluster --wrap="iphop predict --fa_file $i --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir $i.iphop/"	
done

#Run vcontact2 on sequences greater than either (or both) 5000bp or 10000bp

#Extract sequences >10000bp
seqkit seq -m 10000 new.all.keep1_vs2.pass1.fa -o greater.than.10000bp.all.keep1_vs2.pass1.fa

#Extract sequences >10000bp
seqkit seq -m 5000 new.all.keep1_vs2.pass1.fa -o greater.than.5000bp.all.keep1_vs2.pass1.fa

#Run prodigal first
sbatch --job-name=prodigal-meta --cpus-per-task=8 --mem=4000 --time=4:00:00 --partition=cluster --wrap="prodigal -i greater.than.10000bp.all.keep1_vs2.pass1.fa -o greater.than.10000bp.all.keep1_vs2.pass1.fa.genes -a greater.than.10000bp.all.keep1_vs2.pass1.faa"

#run cdhit -c 90, 95, 99, 1
#sbatch --job-name=cdhit90 --cpus-per-task=32 --mem=100000 --time=48:00:00  --partition=cluster --wrap="cd-hit-est -M 100000 -c 0.90 -d 100 -g 1 -aS 0.85 -t 32 -i new.all.keep1_vs2.pass1.fa -o new.cdhit.90.with.aS.all.keep1s.fa"
#tried again but still so slow??
#sbatch --job-name=cdhit90-again --cpus-per-task=32 --mem=10000 --time=120:00:00  --qos=long --partition=cluster --wrap="cd-hit-est -M 100000 -c 0.90 -d 100 -g 1 -aS 0.85 -t 32 -i all.keep1s.vs2.checkv.fa -o cdhit.90.with.aS.all.keep1s.fa"
sbatch --job-name=cdhit95 --cpus-per-task=32 --mem=10000 --time=72:00:00  --qos=long --partition=cluster --wrap="cd-hit-est -M 100000 -c 0.95 -d 100 -g 1 -aS 0.85 -t 32 -i all.keep1s.vs2.checkv.fa -o cdhit.95.with.aS.all.keep1s.fa"
sbatch --job-name=cdhit99 --cpus-per-task=32 --mem=10000 --time=72:00:00  --qos=long --partition=cluster --wrap="cd-hit-est -M 100000 -c 0.99 -d 100 -g 1 -aS 0.85 -t 32 -i all.keep1s.vs2.checkv.fa -o cdhit.99.with.aS.all.keep1s.fa"
sbatch --job-name=cdhit100 --cpus-per-task=32 --mem=10000 --time=72:00:00  --qos=long  --partition=cluster --wrap="cd-hit-est -M 100000 -c 1 -d 100 -g 1 -aS 0.85 -t 32 -i all.keep1s.vs2.checkv.fa -o cdhit.100.with.aS.all.keep1s.fa"
#try again in case 48h is not enough
sbatch --job-name=cdhit100-again --cpus-per-task=32 --mem=10000 --time=72:00:00  --qos=long --partition=cluster --wrap="cd-hit-est -M 100000 -c 1 -d 100 -g 1 -aS 0.85 -t 32 -i new.all.keep1_vs2.pass1.fa -o again.new.cdhit.1.with.aS.all.keep1s.fa"

##Run prodigal on cd-hit clustered sequences
conda activate prodigal
#sbatch --job-name=prodigal90 --cpus-per-task=32 --mem=10000 --time=12:00:00 --partition=cluster --wrap="prodigal -i new.cdhit.90.with.aS.all.keep1s.fa -o new.cdhit.90.with.aS.all.keep1s.fa.genes -a new.cdhit.90.with.aS.all.keep1s.proteins.faa -p meta"
#sbatch --job-name=prodigal95 --cpus-per-task=32 --mem=10000 --time=12:00:00 --partition=cluster --wrap="prodigal -i cdhit.95.with.aS.all.keep1s.fa -o cdhit.95.with.aS.all.keep1s.fa.genes -a cdhit.95.with.aS.all.keep1s.proteins.faa -p meta"
#sbatch --job-name=prodigal99 --cpus-per-task=32 --mem=10000 --time=12:00:00 --partition=cluster --wrap="prodigal -i cdhit.99.with.aS.all.keep1s.fa -o cdhit.99.with.aS.all.keep1s.fa.genes -a cdhit.99.with.aS.all.keep1s.proteins.faa -p meta"
#sbatch --job-name=prodigal100 --cpus-per-task=32 --mem=10000 --time=12:00:00 --partition=cluster --wrap="prodigal -i cdhit.100.with.aS.all.keep1s.fa -o cdhit.100.with.aS.all.keep1s.fa.genes -a cdhit.100.with.aS.all.keep1s.proteins.faa -p meta"

#Then generate gene to genome mapping files
conda activate vcontact2.3
#vcontact2_gene2genome -p new.cdhit.90.with.aS.all.keep1s.proteins.faa -o new.cdhit.90.with.aS.all.keep1s.meta.g2g.csv -s 'Prodigal-FAA'
#vcontact2_gene2genome -p cdhit.95.with.aS.all.keep1s.proteins.faa -o cdhit.95.with.aS.all.keep1s.meta.g2g.csv -s 'Prodigal-FAA'
#vcontact2_gene2genome -p cdhit.99.with.aS.all.keep1s.proteins.faa -o cdhit.99.with.aS.all.keep1s.meta.g2g.csv -s 'Prodigal-FAA'
#vcontact2_gene2genome -p cdhit.100.with.aS.all.keep1s.proteins.faa -o cdhit.100.with.aS.all.keep1s.meta.g2g.csv -s 'Prodigal-FAA'

#Create the bash scripts then sbatch vcontact2
#vc2-100 failed at mem=1028000 so I increased memory to 1500000. Memory couldn't be set to 2000000.

#Check for rRNAs using barrnap
cd /gxfs_work1/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus

for i in *combined.fna
do
	sbatch -p cluster -t 12:00:00 --mem=100000 -J cat -o drep.out -e drep.err -n 12 --wrap="barrnap --threads 12 $i --outseq $i.barrnap"
done

grep ">" *barrnap | cut -f4 -d ":" | sort | uniq > rRNAs.from.checkv.trimmed.reads

#Check out cleanup_cdhit_clusters.sh to see how I preprocessed cdhit clusters for subsequent time-series analysis, vcontact runs, and iphop runs (5k and 10k bp cutoffs)

#Try iphop with cluster reps >10k and >5k
conda activate iphop_new

sbatch --job-name=iphop-test --cpus-per-task=32 --mem=64000  --time=48:00:00 --partition=cluster --wrap="iphop predict --fa_file 10k.cdhit.95.with.aS.all.keep1s.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir 10k.cdhit.95.iphop/"
 
#fasta files too big, split them first
for i in 10k*.fa
do
	mkdir ./$i.split
	cp $i ./$i.split
	cd ./$i.split
	gt splitfasta -numfiles 100 10k*fa
	cd ..
done

##haven't split 5k yet because i want to run a an iphop test
for i in 5k*.fa
do
	mkdir ./$i.split
	cp $i ./$i.split
	cd ./$i.split
	gt splitfasta -numfiles 100 5k*fa
	cd ..
done

rm ./10k*/*fa
rm ./5k*/*fa

#test iphop with a split multifasta file
conda activate iphop_env ##this must be the iphop env that works

for i in ./10k.cdhit.99.with.aS.all.keep1s.fa.split/*fa*
do
sbatch --job-name=iphop --cpus-per-task=8 --mem=64000  --time=48:00:00 --partition=cluster --wrap="iphop predict --fa_file $i --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir 10k.cdhit.99.iphop/"
done

#making iphop work again ugh
#with test.fa first ugh
iphop predict --fa_file test/test_input_phages.fna --db_dir iphop_db/Test_db/ --out_dir test2_env/ --debug --num_threads 10
#omg it fucking worked in iphop_new??? i'm so confusedt

#ay wait let's run it again with the REAL database
sbatch --job-name=iphop --cpus-per-task=8 --mem=64000  --time=1:00:00 --partition=cluster --wrap="iphop predict --fa_file test/test_input_phages.fna --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir test4_env/ --debug --num_threads 10"

#ok it worked!!!
#apparently restarting your command line everytime you change variables is a good practice —— exiting the server is not enough

#now run a subset of real data in iphop_env
sbatch --job-name=iphop_test --cpus-per-task=10 --mem=64000  --time=48:00:00 --partition=cluster --wrap="iphop predict --fa_file 10k.cdhit.99.with.aS.all.keep1s.fa.26 --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir iphop.fa.26.out --num_threads 10"

#it's not splitting the fiiile. file size has nothing to do with it. spaces in fasta headers have nothing to do with it. why is BioseqIO not parsing and writing (aka splitting) my file??
#ugh characters that can't go in filenames like | and / must be deleted duh??

#so now run it without those special characters
iphop predict --fa_file test.26.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir test.out --num_threads 10
#now it successfully splits the fasta file but it fails at wish again wtf???? this doesn't happen with the test data???
#ok it finished when input file is ran in installation directory

#ok so i just found out that rafah fails when iphop is ran on slurm.
#the next is run this via slurm... maybe atdd to PATH first then to LD_LIBRARY_PATH — didn't work!

#asked help from HPC guys. meanwhile, try making soft links and see if it works.
ln -s $WORK/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus/cdhit_latest/iphop.out/10k.cdhit.99.with.aS.all.keep1s.fa.split/test.26.fa $WORK/iphop_renamed/soft.link.fa

#now try making this soft link work
iphop predict --fa_file soft.link.fa --db_dir /gxfs_work1/geomar/smomw535/databases/iphop_db/Sept_2021_pub --out_dir soft.link.out --num_threads 10
#ok soft link works!

#Karsten from HPC replied — see his email for instructions. There was a typo though, so do this instead:
source /gxfs_home/sw/spack/spack0.16.0/usr/opt/spack/linux-rhel8-x86_64/gcc-10.2.0/miniconda3-4.12.0-fnuv53yaf2s6fs2eukuipnlugl564wtw/etc/profile.d/conda.sh

#download latest version of iphop database
sbatch --job-name=iphop_db --cpus-per-task=10 --mem=32000  --time=48:00:00 --partition=data --wrap="echo yes | iphop download --db_dir iphop_db/"

#Now try it with real data
#but first fix fasta headers with:
for i in *; do sed -i 's/[\/||]/_/g; s/__/_/g' $i; done

#then run iphop_sbatch.sh

#rename "fixed" fasta headers with the original
sed 's/_ccs/\/ccs/g ; s/_full/\||full/g ; s/_lt2gene/\||lt2gene/g ; s/_partial/\||partial/g' 10k.cdhit.*.with.aS.all.keep1s.fa.split/iphop/task_*/*enome* | head

##run iphop (codes are in the iphop results directory somewhere)

#run gtdbtk on Taylor's MAGs

sbatch --mem=150000 --time=24:00:00 -n 8 -J gtdb -o gtdb.egc.out -e gtdb.egc.err --wrap="gtdbtk classify_wf --genome_dir RAS_EGC_MAGs/ --out_dir RAS_EGC_MAGs.gtdb -x fa --cpus 8"
sbatch --mem=150000 --time=24:00:00 -n 8 -J gtdb -o gtdb.fram18.out -e gtdb.fram18.err --wrap="gtdbtk classify_wf --genome_dir FRAM18_MAGs/ --out_dir FRAM18_MAGs.gtdb -x fa --cpus 8"
sbatch --mem=150000 --time=24:00:00 -n 8 -J gtdb -o gtdb.wsc20.out -e gtdb.wsc20.err --wrap="gtdbtk classify_wf --genome_dir FRAM_WSC20_MAGs/ --out_dir FRAM_WSC20_MAGs.gtdb -x fa --cpus 8"

# run checkm2 on Taylor's MAGs

sbatch --mem=10000 --time=24:00:00 -n 30 -J checkm2 -o checkm2.egc.out -e checkm2.egc.err --wrap="checkm2 predict --threads 30 --input RAS_EGC_MAGs --output-directory RAS_EGC_MAGs.checkm2 --extension .fa"
sbatch --mem=10000 --time=24:00:00 -n 30 -J checkm2 -o checkm2.fram18.out -e checkm2.fram18.err --wrap="checkm2 predict --threads 30 --input FRAM18_MAGs --output-directory FRAM18_MAGs.checkm2 --extension .fa"
sbatch --mem=10000 --time=24:00:00 -n 30 -J checkm2 -o checkm2.wsc20.out -e checkm2.wsc20.err --wrap="checkm2 predict --threads 30 --input FRAM_WSC20_MAGs --output-directory FRAM_WSC20_MAGs.checkm2 --extension .fa"


# run checkm on Taylor's MAGs
sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.egc.arc.out -n 8 --wrap="checkm taxonomy_wf domain Archaea -x fa -t 8 -f RAS_EGC_MAGs.checkm.archaea.out.tsv --tab_table RAS_EGC_MAGs/ RAS_EGC_MAGs.checkm/"
sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.egc.bac.out -n 8 --wrap="checkm taxonomy_wf domain Bacteria -x fa -t 8 -f RAS_EGC_MAGs.checkm.bacteria.out.tsv --tab_table RAS_EGC_MAGs/ RAS_EGC_MAGs.checkm/"

sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.fram18.arc.out -n 8 --wrap="checkm taxonomy_wf domain Archaea -x fa -t 8 -f FRAM18_MAGs.checkm.archaea.out.tsv --tab_table FRAM18_MAGs/ FRAM18_MAGs.checkm/"
sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.fram18.bac.out -n 8 --wrap="checkm taxonomy_wf domain Bacteria -x fa -t 8 -f FRAM18_MAGs.checkm.bacteria.out.tsv --tab_table FRAM18_MAGs/ FRAM18_MAGs.checkm/"

sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.wsc20.arc.out -n 8 --wrap="checkm taxonomy_wf domain Archaea -x fa -t 8 -f FRAM_WSC20_MAGs.checkm.archaea.out.tsv --tab_table FRAM_WSC20_MAGs/ FRAM_WSC20_MAGs.checkm/"
sbatch -p cluster -t 24:00:00 --mem=100000 -J checkm -o checkm.wsc20.bac.out -n 8 --wrap="checkm taxonomy_wf domain Bacteria -x fa -t 8 -f FRAM_WSC20_MAGs.checkm.bacteria.out.tsv --tab_table FRAM_WSC20_MAGs/ FRAM_WSC20_MAGs.checkm/"

# add taylor's MAGs to database
sbatch -p cluster -t 4:00:00 --mem=100000 -J gtdb  -n 10 --wrap="gtdbtk de_novo_wf --genome_dir $WORK/fram_mags/all_FRAM_MAGs/ --bacteria --outgroup_taxon p__Patescibacteria --out_dir $WORK/databases/fram_mags_gtdb --cpus 10 --force --extension fa"
sbatch -p cluster -t 4:00:00 --mem=100000 -J gtdb  -n 10 --wrap="gtdbtk de_novo_wf --genome_dir $WORK/fram_mags/all_FRAM_MAGs/ --archaea --outgroup_taxon p__Altarchaeota --out_dir $WORK/databases/fram_mags_gtdb --cpus 10 --force --extension fa"
#didn't work probably because time was too short and/or outgroup sequence was not in input fasta file

#try again with random outgroup from identified phyla
sbatch -p cluster -t 4:00:00 --mem=100000 -J gtdb  -n 10 --wrap="gtdbtk de_novo_wf --genome_dir $WORK/fram_mags/all_FRAM_MAGs/ --bacteria --outgroup_taxon p__Actinomycetota --out_dir $WORK/databases/fram_mags_gtdb --cpus 10 --force --extension fa"

sbatch -p cluster -t 24:00:00 --mem=100000 -J gtdb  -n 20 --wrap="gtdbtk de_novo_wf --genome_dir $WORK/fram_mags/all_FRAM_MAGs/ --bacteria --outgroup_taxon p__Patescibacteria --out_dir $WORK/databases/fram_mags_gtdb --cpus 20 --force --extension fa"
#next time use 20 CPUs

sbatch -p cluster -t 24:00:00 --mem=100000 -J gtdb  -n 20 --wrap="gtdbtk de_novo_wf --genome_dir $WORK/fram_mags/all_FRAM_MAGs/ --archaea --outgroup_taxon p__Altiarchaeota --out_dir $WORK/databases/fram_mags_gtdb --cpus 20 --force --extension fa"

# actually add the bacterial MAGs and archaeal MAGs
sbatch -p cluster -t 24:00:00 --mem=10000 -J iphop_add  -n 10 --wrap="iphop add_to_db --fna_dir $WORK/fram_mags/all_FRAM_MAGs/ --gtdb_dir $WORK/databases/fram_mags_gtdb/ --out_dir $WORK/databases/iphop_db/Sept_2021_pub_rw_FRAM_hosts --db_dir $WORK/databases/iphop_db/Sept_2021_pub_rw/"

#removed the bad MAGs

#run iphop again (for old database YOY and new database)
##split fasta file into 100 fasta files
##rerun these, i.e. unfinished iphop jobs

##remove bad symbols
##rename??

##MArVD2
MArVD2.py -i m5514_AA_1#12.provirus.virus.fa --load-model /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/marvd2_working_model/rf_model.pkl -o marvd2_test.out \
--db-pvog AllvogHMMprofiles.hmm --db-nr nr.faa --marine-jackhmmer-db pVOG_prots_ref_marine_pVOG.faa \
--viral-refseq-txt viruses.txt --pvog-dir /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs \
--db-accession2tax /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs/prot.accession2taxid.trimmed \
--dbs-dir /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs

MArVD2.py -i MArVD2_files_zip/model_build_datasets/toy_pred_dataset.fasta \
--load-model /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/marvd2_working_model/rf_model.pkl -o marvd2_test.out \
--db-pvog AllvogHMMprofiles.hmm --db-nr nr.faa --marine-jackhmmer-db pVOG_prots_ref_marine_pVOG.faa \
--viral-refseq-txt viruses.txt --pvog-dir /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs \
--db-accession2tax /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs/prot.accession2taxid.trimmed \
--dbs-dir /gxfs_work1/geomar/smomw535/marvd2/MArVD2_files_zip/other_dbs



sbatch -p base -t 24:00:00 --mem=10000 -J iphop_add  -n 10 --wrap="iphop add_to_db --fna_dir /gxfs_work/geomar/smomw535/fram_mags/all_FRAM_MAGs/ --gtdb_dir /gxfs_work/geomar/smomw535/databases/fram_mags_gtdb/ --out_dir /gxfs_work/geomar/smomw535/databases/iphop_db/Aug_2023_pub_rw_FRAM_hosts --db_dir /gxfs_work/geomar/smomw535/databases/iphop_db/Aug_2023_pub_rw/"



sbatch -p base -t 2:00:00 --mem=60000 -J iphop_test  -n 10 --wrap="iphop predict --fa_file ./iphop/test/test_input_phages.fna --db_dir $WORK/databases/iphop_db/Aug_2023_pub_rw/ --out_dir ./iphop_test_results/iphop_test_23Jan24_new_db"




total_lines=$(for i in $(cat sample.names); do wc -l "${i}_screened.tsv"; done | awk '{sum += $1} END {print sum}')
echo "Total lines: $total_lines"
echo "concat_wsc_screened.tsv lines:" $(wc -l < concat_wsc_screened.tsv)


#Run vpf tools on Keep 1 reads

for i in `cat sample.names`
sbatch -p base -t 48:00:00 --mem=10000 -J vpf -n 10 --wrap="
singularity exec --bind /gxfs_work/geomar/smomw535/vpf-class-data:/opt/vpf-tools/vpf-data \
                 --bind /gxfs_work/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combined/provirus_virus/$i_new_clean.fna:/opt/vpf-tools/input-sequences/$i_new_clean.fna \
                 --bind /gxfs_work/geomar/smomw535/arctic_metag/all_data/all_data_fasta/splitfasta/virsorter_viral_combinedprovirus_virus/$i_vpf-class.out:/opt/vpf-tools/outputs \
                 $WORK/vpf-tools_latest.sif \
                 vpf-class -i /opt/vpf-tools/input-sequences/$i_new_clean.fna -o /opt/vpf-tools/outputs/test-classified"

done


#Filter our family.tsv files
for i in `cat sample.names`
do
	echo $i
	awk 'NR==1 || ($3 > 0.5 && $5 > 0.5)' $i"_vpf-class.out"/*/family.tsv | wc -l
done
