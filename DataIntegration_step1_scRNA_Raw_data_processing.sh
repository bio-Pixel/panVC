######################################################################################
############################1. Processing Genomics 10X data###########################

######################################################################################
#########################Download files from sra (public data)########################
prefetch -O ./srr/ --option-file ./SRR_Acc_List.txt --max-size 999999999

########################Convert sra to fastq files####################################
for i in `cat ./SRR_Acc_List.txt`;do
mkdir ./fastq/$i
/hengya/apps/sratoolkit.3.0.0-centos_linux64/bin/fasterq-dump -e 20 ./$i/$i.sra --split-files  -O ./fastq/$i/
done

#################Convert FASTQ files to the format used by 10X Genomics###############
###########Processing two FASTQ files
cat SRR_Acc_List.txt | while read i 
do 
mv ${i}_1*.gz ${i}_S1_L001_R1_001.fastq.gz
mv ${i}_2*.gz ${i}_S1_L001_R2_001.fastq.gz
done
###########Processing three FASTQ files
cat SRR_Acc_List.txt | while read i 
do 
mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz
mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz
mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz
done

###########################Run 10X Genomics data by CellRanger##########################
cellranger count --id h2_TCR --transcriptome $index_file_folder --fastqs $fastq_folder --sample $sample_name --nosecondary


########################################################################################
############################2. Processing SeekOne data##################################
seekonetools rna run --fq1 $fastq1_path --fq2 $fastq2_path --samplename $sample_name --outdir $out_dir --chemistry DDV2 --genomeDir $index_file_folder --gtf $gtf_file --star_path $star_path --region exon


########################################################################################
############################3. Processing Singleron data################################
celescope rna sample --outdir .//${sample_name}/00.sample --sample $sample_name --thread 40 --chemistry auto  --fq1 $fastq1_path
celescope rna barcode --outdir .//${sample_name}/01.barcode --sample $sample_name --thread 40 --chemistry auto --lowNum 2  --fq1 $fastq1_path --fq2 $fastq2_path
celescope rna cutadapt --outdir .//${sample_name}/02.cutadapt --sample $sample_name --thread 40 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq .//${sample_name}/01.barcode/${sample_name}_2.fq 
celescope rna star --outdir .//${sample_name}/03.star --sample $sample_name --thread 40 --genomeDir $genomedir --outFilterMultimapNmax 1 --starMem 30  --fq .//${sample_name}/02.cutadapt/${sample_name}_clean_2.fq 
celescope rna featureCounts --outdir .//${sample_name}/04.featureCounts --sample $sample_name --thread 40 --gtf_type gene --genomeDir $genomedir  --input .//${sample_name}/03.star/${sample_name}_Aligned.sortedByCoord.out.bam 
celescope rna count --outdir .//${sample_name}/05.count --sample $sample_name --thread 40 --genomeDir $genomedir --expected_cell_num 3000 --cell_calling_method EmptyDrops_CR  --bam .//${sample_name}/04.featureCounts/${sample_name}_name_sorted.bam --force_cell_num None 
celescope rna analysis --outdir .//${sample_name}/06.analysis --sample $sample_name --thread 40 --genomeDir $genomedir  --matrix_file .//${sample_name}/05.count/${sample_name}_filtered_feature_bc_matrix 


#########################################################################################
#############################4. Processing BD data#######################################
cd $path
outDir=$path
temPrefix=${path}/temp/tempfile/
temDir=${path}/tmp/tmpfile
cwL=${path}/rhapsody_wta_1.9.1.cwl
ymL=${path}/template_wta_1.9.1.yml

echo "Begin:" `date` && \

cwl-runner \
        --parallel \
        --outdir $outDir \
        --tmpdir-prefix $temPrefix \
        --tmp-outdir-prefix $temDir \
        --rm-tmpdir \
        $cwL $ymL

echo "End: " `date`
