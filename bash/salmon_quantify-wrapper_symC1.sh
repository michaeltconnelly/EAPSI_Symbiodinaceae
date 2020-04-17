#!/bin/bash
#./bash/salmon_quantify-wrapper_symC1.sh
#purpose: quantify non-coral trimmed RNAseq reads against the Cladocopium C1 reference transcriptome (Levin et al. 2016)
#To start this job from the EAPSI_Symbiodinaceae directory, use:
#bsub -P transcriptomics < /scratch/projects/transcriptomics/mikeconnelly/projects/EAPSI_Pocillopora_AxH/bash/salmon_quantify-wrapper_symC1.sh

#BSUB -J salmonwrap_symC1
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o salmonwrap_symC1%J.out
#BSUB -e salmonwrap_symC1%J.err
#BSUB -n 8
#BSUB -u m.connelly1@umiami.edu

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/EAPSI_Symbiodinaceae"
EAPSIsamples="Wt1-1a Wt1-1b Wt1-1c Wt1-4a Wt1-4b Wt1-4c Wt1-5a Wt1-5b Wt1-5c Wt1-6a Wt1-6b Wt1-6c Wt2-1a Wt2-1b Wt2-1c Wt2-4a Wt2-4b Wt2-4c Wt2-5a Wt2-5b Wt2-5c Wt2-6a Wt2-6b Wt2-6c Hw1-1a Hw1-1b Hw1-1c Hw1-4a Hw1-4b Hw1-4c Hw1-5a Hw1-5b Hw1-5c Hw1-6a Hw1-6b Hw1-6c Hw2-1a Hw2-1b Hw2-1c Hw2-4a Hw2-4b Hw2-4c Hw2-5a Hw2-5b Hw2-5c Hw2-6b Hw2-6c"

#lets me know which files are being processed
echo "These are the reads to be quantified against the Cladocopium C1 MI reference transcriptome: $EAPSIsamples"

#loop to automate generation of scripts to direct sequence file trimming
for EAPSIsample in $EAPSIsamples
do \
echo "Quantifying ${EAPSIsample}"

#   input BSUB commands
echo '#!/bin/bash' > "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -q general' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -J '"${EAPSIsample}"_salmonquant_symC1'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -o '"${prodir}"/logfiles/"$EAPSIsample"_salmonquant_symC1%J.out'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -e '"${prodir}"/errorfiles/"$EAPSIsample"_salmonquant_symC1%J.err'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -n 8' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo '#BSUB -W 4:00' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job

#   input command to run salmon
echo "${mcs}"/programs/salmon-latest_linux_x86_64/bin/salmon quant \
-i ${prodir}/data/indices/symC1_MI_index \
-l SR \
-r ${prodir}/outputs/STARalign_Pdam/${EAPSIsample}_PdamUnmapped.out.mate1 \
--validateMappings -o ${prodir}/outputs/${EAPSIsample}transcripts_quant >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job

#lets me know file is done
echo 'echo' "Salmon quantification of $EAPSIsample complete"'' >> "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
echo "Salmon quantification script of $EAPSIsample submitted"
#   submit generated trimming script to job queue
bsub < "${prodir}"/bash/jobs/"${EAPSIsample}"_salmonquant_symC1.job
done
