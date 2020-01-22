#!/bin/bash
#BSUB -J STARwrap_symC1
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o STARwrap_symC1%J.out
#BSUB -e STARwrap_symC1%J.err
#BSUB -n 8
#BSUB -u dxc947@miami.edu

#specify variable containing sequence file prefixes and directory paths
dcs="/scratch/projects/transcriptomics/demi"
coldir="/scratch/projects/transcriptomics/demi/sequences/EAPSI/"
exp="heat"
EAPSIsamples="Wt1-4a Wt1-4b Wt1-4c Wt1-6a Wt1-6b Wt1-6c Wt2-4a Wt2-4b Wt2-4c Wt2-6a Wt2-6b Wt2-6c Hw1-4a Hw1-4b Hw1-4c Hw1-6a Hw1-6b Hw1-6c Hw2-4a Hw2-4b Hw2-4c Hw2-6b Hw2-6c"

#lets me know which files are being processed
echo "These are the reads to be aligned to the Cladocopium goreaui reference genome: $EAPSIsamples"

#loop to automate generation of scripts to direct sequence file trimming
for EAPSIsample in $EAPSIsamples
do \
echo "Aligning ${EAPSIsample}"

#   input BSUB commands
echo '#!/bin/bash' > "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -q general' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -J '"${EAPSIsample}"_staralign_symC1'' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -o '"${coldir}"/"${exp}"/logfiles/"$EAPSIsample"staralign_SymC1%J.out'' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -e '"${coldir}"/"${exp}"/errorfiles/"$EAPSIsample"staralign_SymC1%J.err'' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -n 8' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo '#BSUB -W 4:00' >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job

#   input command to run STAR aligner
echo ${dcs}/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--quantMode TranscriptomeSAM \
--runThreadN 16 \
--readFilesIn ${coldir}/${exp}/STARalign_Pdam/${EAPSIsample}_Pdam-Unmapped.out.mate1 \
--genomeDir ${dcs}/sequences/genomes/symbiodinium/STARindex \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${dcs}/sequences/genomes/symbiodinium/symC1_genome.gff \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${coldir}/${exp}/STARalign_SymC1/${EAPSIsample}_SymC1 >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job

#lets me know file is done
echo 'echo' "STAR alignment of $EAPSIsample complete" >> "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
echo "STAR alignment script of $EAPSIsample submitted"
#   submit generated trimming script to job queue
bsub < "${coldir}"/"${exp}"/scripts/"${EAPSIsample}"_staralign_symC1.job
done
