#!/bin/bash
#~/scripts/EAPSI-symITS2/STARalign_ITS2.sh
#purpose: align trimmed non-coral RNAseq reads against the GeoSymbio ITS2 dataset using STAR

#BSUB -J staralign_its2
#BSUB -q bigmem
#BSUB -P transcriptomics
#BSUB -o star_its2%J.out
#BSUB -e star_its2%J.err
#BSUB -n 8
#BSUB -W 6:00
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/EAPSI_symITS2"
EAPSIsamples="Wt1-6a Wt1-6b Wt1-6c Wt2-6a Wt2-6b Wt2-6c Hw1-6a Hw1-6b Hw1-6c Hw2-6b Hw2-6c"

echo "These are the reads to be aligned to the Symbiodatabaceae reference: $EAPSIsamples"

#loop to automate generation of scripts to direct sequence file trimming
for EAPSIsample in $EAPSIsamples
do \
echo "#Aligning ${EAPSIsample}"

#   input BSUB commands
echo '#!/bin/bash' > "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -q bigmem' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -J '"${EAPSIsample}"_staralign_its2'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -o '"${prodir}"/outputs/logfiles/"$EAPSIsample"staralign_its2%J.out'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -e '"${prodir}"/outputs/errorfiles/"$EAPSIsample"staralign_its2%J.err'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -n 8' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo '#BSUB -W 4:00' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job

#   input command to run STAR aligner
echo "Aligning ${EAPSIsample}" >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo ${mcs}/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn ${prodir}/data/reads/${EAPSIsample}_PdamUnmapped.out.mate1.fastq \
--genomeDir ${prodir}/data/indices/STARindex \
--sjdbGTFtagExonParentTranscript Parent \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${prodir}/outputs/${EAPSIsample}_ITS2 >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
#lets me know file is done
echo "STAR alignment of $EAPSIsample complete"
#lets me know file is done
echo 'echo' "STAR alignment of $EAPSIsample complete"'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
echo "STAR alignment script of $EAPSIsample submitted"
#   submit generated trimming script to job queue
bsub < "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_its2.job
done
