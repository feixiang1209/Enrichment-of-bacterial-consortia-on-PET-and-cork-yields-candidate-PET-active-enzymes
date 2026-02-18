##Triiming with fastq

for file in *1.fastq.gz
do
name2=${file%1.fastq.gz}2.fastq.gz
name=$(echo $file | cut -d'_' -f1)
fastp \
-i $file \
-I $name2 \
-o trimmed/${name}_trimmed_R1.fastq.gz \
-O trimmed/${name}_trimmed_R2.fastq.gz \
--thread 30 \
-h trimmed/${name}_report.html \
-j trimmed/${name}_report.json
done



##Host removal with bowtie2

for file in *R1.fastq.gz
do
name2=${file%R1.fastq.gz}R2.fastq.gz
name=${file%_trimmed_R1.fastq.gz}
bowtie2 -x /path/to/hg38_mm10_index  -1 $file -2 $name2 -S host_removed/${name}.sam --un-conc host_removed/${name}_host_removed -p 100
done

## Assemble with megahit

for file in *1.fastq
do
megahit -1 $file -2 ${file%1.fastq}2.fastq -o megahit_output/${file%_host_removed_1.fastq}_megahit_output -t 100
done


## QC of contigs with quast

path=/path/to/megahit_output

for dir in *megahit_output
do
cd $path/$dir
quast -o ${dir%_megahit_output}_quast_output final.contigs.fa
done

##Binning and bin refinement with metawrap

for dir in *_megahit_output
do
metawrap binning -o metawrap_binning/${dir%_megahit_output} -t 100 -a $dir/final.contigs.fa \
--metabat2 --maxbin2 --concoct \
 /path/to/${dir%_megahit_output}_host_removed_1.fastq  /path/to/${dir%_megahit_output}_host_removed_2.fastq
done

path=/path/to/megahit_output
cd $path/metawrap_binning
for dir in *
do cd $path/metawrap_binning/$dir
metawrap bin_refinement -o $path/metawrap_binning/$dir/refinement_output -t 100 \
-A $path/metawrap_binning/$dir/metabat2_bins/ \
-B $path/metawrap_binning/$dir/maxbin2_bins/ \
-C $path/metawrap_binning/$dir/concoct_bins/ \
-c 70 -x 10
done




## De-rep for Cork, PET and PET_Cork
#1. rename and de-rep
# (1) PET
for dir in 10MicrobesPET 11MicrobesPET 12MicrobesPET 
do cd /path/to/$dir/refinement_output/metawrap_70_10_bins/  || exit
   for file in *.fa
   do cp $file /path/to/de-rep/PET/${dir}_${file}
   done
done 

dRep dereplicate derep_out \
  -g /path/to/de-rep/PET/*.fa \
  -p 50 \
  -comp 70 \
  -con 10 \
  --S_algorithm fastANI \
  --P_ani 0.95

cd /path/to/de-rep/PET/derep_out/dereplicated_genomes
#rename the de-reped bins 
echo -e "Original_Name\tNew_Name" > rename_mapping.txt  
counter=1  
for file in $(ls *MicrobesPET_bin.*.fa | sort -V); do
    new_name="Cork_bin.${counter}.fa"
    mv "$file" "$new_name"
    echo -e "$file\t$new_name" >> rename_mapping.txt   
    ((counter++))
done
  
# (2) Cork
for dir in 15MicrobesCork 17MicrobesCork 19MicrobesCork 
do cd /path/to/$dir/refinement_output/metawrap_70_10_bins/  || exit
   for file in *.fa
   do cp $file /path/to/de-rep/Cork/${dir}_${file}
   done
done

Rep dereplicate derep_out \
  -g /path/to/de-rep/Cork/*.fa \
  -p 50 \
  -comp 70 \
  -con 10 \
  --S_algorithm fastANI \
  --P_ani 0.95

cd /path/to/de-rep/Cork/derep_out/dereplicated_genomes
#rename the de-reped bins 
echo -e "Original_Name\tNew_Name" > rename_mapping.txt  
counter=1  
for file in $(ls *MicrobesCork_bin.*.fa | sort -V); do
    new_name="Cork_bin.${counter}.fa"
    mv "$file" "$new_name"
    echo -e "$file\t$new_name" >> rename_mapping.txt   
    ((counter++))
done 


#(3) PET_Cork
for dir in 6MicrobesPETCork 7MicrobesPETCork 8MicrobesPETCork
do cd /path/to/$dir/refinement_output/metawrap_70_10_bins/  || exit
   for file in *.fa
   do cp $file /path/to/de-rep/PET_Cork/${dir}_${file}
   done
done 

dRep dereplicate derep_out \
  -g /path/to/de-rep/PET_Cork/*.fa \
  -p 50 \
  -comp 70 \
  -con 10 \
  --S_algorithm fastANI \
  --P_ani 0.95

cd /path/to/de-rep/PET_Cork/derep_out/dereplicated_genomes
#rename the de-reped bins 
echo -e "Original_Name\tNew_Name" > rename_mapping.txt  
counter=1  
for file in $(ls *MicrobesPET_Cork_bin.*.fa | sort -V); do
    new_name="PET_Cork_bin.${counter}.fa"
    mv "$file" "$new_name"
    echo -e "$file\t$new_name" >> rename_mapping.txt   
    ((counter++))
done


#2. pool all de-reped bins and run de-rep again
for dir in Cork  PET  PET_Cork
do
cp $dir/derep_out/dereplicated_genomes/*.fa pooled_bins/
done

cd /path/to/pooled_bins/
dRep dereplicate derep_out \
  -g ./*.fa \
  -p 50 \
  -comp 70 \
  -con 10 \
  --S_algorithm fastANI \
  --P_ani 0.95


#3. rename MAGs based on de-rep result
: <<'COMMENT'
Cork_bin.1.fa   Cork_1.fa
Cork_bin.10.fa  Cork_2.fa
Cork_bin.13.fa  Cork_3.fa
Cork_bin.14.fa  Cork_4.fa
Cork_bin.15.fa  Cork_5.fa
Cork_bin.16.fa  Cork_6.fa
Cork_bin.18.fa  Cork_7.fa
Cork_bin.19.fa  Cork_8.fa
Cork_bin.20.fa  Cork_9.fa
Cork_bin.23.fa  Cork_10.fa
Cork_bin.24.fa  Cork_11.fa
Cork_bin.25.fa  Cork_12.fa
Cork_bin.27.fa  Cork_13.fa
Cork_bin.28.fa  Cork_14.fa
Cork_bin.29.fa  Cork_15.fa
Cork_bin.3.fa   Cork_16.fa
Cork_bin.30.fa  Cork_17.fa
Cork_bin.31.fa  Cork_18.fa
Cork_bin.32.fa  Cork_19.fa
Cork_bin.34.fa  Cork_20.fa
Cork_bin.35.fa  Cork_21.fa
Cork_bin.5.fa   Cork_22.fa
Cork_bin.6.fa   Cork_23.fa
Cork_bin.7.fa   Cork_24.fa
Cork_bin.9.fa   Cork_25.fa
PET_bin.10.fa   PET_1.fa
PET_bin.11.fa   PET_2.fa
PET_bin.12.fa   PET_3.fa
PET_bin.13.fa   PET_4.fa
PET_bin.16.fa   PET_5.fa
PET_bin.18.fa   PET_6.fa
PET_bin.19.fa   PET_7.fa
PET_bin.2.fa    PET_8.fa
PET_bin.4.fa    PET_9.fa
PET_bin.9.fa    PET_10.fa
PETCork_bin.10.fa       PETCork_1.fa
PETCork_bin.11.fa       PETCork_2.fa
PETCork_bin.15.fa       PETCork_3.fa
PETCork_bin.17.fa       PETCork_4.fa
PETCork_bin.19.fa       PETCork_5.fa
PETCork_bin.2.fa        PETCork_6.fa
PETCork_bin.22.fa       PETCork_7.fa
PETCork_bin.25.fa       PETCork_8.fa
PETCork_bin.26.fa       PETCork_9.fa
PETCork_bin.28.fa       PETCork_10.fa
PETCork_bin.29.fa       PETCork_11.fa
PETCork_bin.31.fa       PETCork_12.fa
PETCork_bin.34.fa       PETCork_13.fa
PETCork_bin.35.fa       PETCork_14.fa
PETCork_bin.6.fa        PETCork_15.fa
PETCork_bin.7.fa        PETCork_16.fa
PET_bin.8.fa    PET_PETCork_1.fa
PETCork_bin.8.fa        PET_PETCork_2.fa
PETCork_bin.16.fa       PET_PETCork_3.fa
PET_bin.3.fa    PET_PETCork_4.fa
PET_bin.17.fa   PET_PETCork_5.fa
PETCork_bin.14.fa       PET_PETCork_6.fa
PETCork_bin.24.fa       PET_PETCork_7.fa
PET_bin.6.fa    PET_PETCork_8.fa
Cork_bin.8.fa   Cork_PETCork_1.fa
Cork_bin.26.fa  Cork_PETCork_2.fa
PETCork_bin.21.fa       Cork_PETCork_3.fa
PETCork_bin.23.fa       Cork_PETCork_4.fa
Cork_bin.2.fa   Cork_PETCork_5.fa
PETCork_bin.3.fa        Cork_PETCork_6.fa
PETCork_bin.1.fa        Cork_PETCork_7.fa
PETCork_bin.12.fa       Cork_PETCork_8.fa
Cork_bin.36.fa  Cork_PETCork_9.fa
PETCork_bin.30.fa       Cork_PETCork_10.fa
PETCork_bin.32.fa       PET_Cork_PETCork_1.fa
COMMENT



#4. rename contigs with MAG name
mkdir -p MAGs
for file in *.fa; do
  sed "s/^>/>${file%.fa}|/" "$file" > MAGs/${file%.fa}.fna
done



## Gene/protein prediction
cd /path/to/NAGs
for file in *.fna
do
mkdir -p  prodigal_output/${file%.fna}
prodigal -i $file -o prodigal_output/${file%.fna}/${file%.fna}.gff -f gff \
-a prodigal_output/${file%.fna}/${file%.fna}.faa -d prodigal_output/${file%.fna}/${file%.fna}.fna -p single
done

## Diamond with PAZy/PETase

# Ultra
for dir in Cork* PET*
do
diamond blastp -d /path/to/PAZy_db -q  $dir/${dir}.faa -o $dir/ultra_output.txt --id 70 --query-cover 93.8
done

for dir in Cork* PET*
do
mkdir -p $dir/PETases_diamond_output
diamond blastp -d /path/to/PETases_db -q  $dir/${dir}.faa -o $dir/PETases_diamond_output/ultra_output.txt --id 70 --query-cover 93.8
done

#High
for dir in Cork* PET*
do
diamond blastp -d /path/to/PAZy_db -q  $dir/${dir}.faa -o $dir/high_output.txt --id 50.4 --query-cover 81.5
done

for dir in Cork* PET*
do
mkdir -p $dir/PETases_diamond_output
diamond blastp -d /path/to/PETases_db -q  $dir/${dir}.faa -o $dir/PETases_diamond_output/high_output.txt --id 50.4 --query-cover 81.5
done

## Functional annotation using Eggnog-mapper

for dir in Cork* PET*
do cd /path/to/$dir/
emapper.py -i  ${dir}.faa   -o emapper_output/$dir --cpu 50   --data_dir /path/to/eggnog-mapper_db/
done


##Quantify MAGs using strobe align and coverM

cd /path/to/MAGs
cat *.fna > combined_bins/combined.fna

for dir in 10MicrobesPET  12MicrobesPET  15MicrobesCork  19MicrobesCork    7MicrobesPETCork  11MicrobesPET  14MicrobesPET  17MicrobesCork  6MicrobesPETCork  8MicrobesPETCork
mkdir -p $dir
strobealign -t 50 -o $dir/${dir}_alignment.sam \
/path/to/combined.fna \
/path/to/10MicrobesPET_host_removed_1.fastq \
/path/to/10MicrobesPET_host_removed_2.fastq 

samtools view -bS $dir/${dir}_alignment.sam | samtools sort -o $dir/${dir}_alignment.sorted.bam
samtools index $dir/${dir}_alignment.sorted.bam
done

coverm genome \
    -b 10MicrobesPET/10MicrobesPET_alignment.sorted.bam \
11MicrobesPET/11MicrobesPET_alignment.sorted.bam \
12MicrobesPET/12MicrobesPET_alignment.sorted.bam \
14MicrobesPET/14MicrobesPET_alignment.sorted.bam \
15MicrobesCork/15MicrobesCork_alignment.sorted.bam \
17MicrobesCork/17MicrobesCork_alignment.sorted.bam \
19MicrobesCork/19MicrobesCork_alignment.sorted.bam \
6MicrobesPETCork/6MicrobesPETCork_alignment.sorted.bam \
7MicrobesPETCork/7MicrobesPETCork_alignment.sorted.bam \
8MicrobesPETCork/8MicrobesPETCork_alignment.sorted.bam \
    -t 50 \
    --genome-fasta-directory /path/to/MAGs \
    --methods rpkm \
    --output-file coverm_results_all.tsv


## Taxanomy annotation for MAGs
cd /path/to/MAGs
gtdbtk classify_wf --genome_dir ./  \
--out_dir ./gtdbtk_output   --extension fna   --skip_ani_screen  --cpus 50

