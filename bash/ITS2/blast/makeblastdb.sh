#!/bin/bash
#./bash/ITS2/blast/makeblastdb.sh
#purpose: Create BLAST databases

#Make nucleotide blast databases from Cladocopium genome fasta files
makeblastdb \
-dbtype nucl \
-parse_seqids \
-in data/refs/SymbC1.Genome.Scaffolds.fasta \
-input_type 'fasta' \
-out data/refs/SymC1 \
-max_file_sz "4GB"

#Make nucleotide blast databases from Symbiodatabaceae fasta files
makeblastdb \
-dbtype nucl \
-parse_seqids \
-in data/GeoSymbio/Symbiodatabaceae.fasta \
-input_type 'fasta' \
-out data/refs/Symbiodatabaceae \
-max_file_sz "4GB"

#Make nucleotide blast databases from Symbiodatabaceae fasta files
makeblastdb \
-dbtype nucl \
-parse_seqids \
-in data/GeoSymbio/GeoSymbio_ITS2_LocalDatabase.fasta \
-input_type 'fasta' \
-out data/refs/GeoSymbio \
-max_file_sz "4GB"
