### Running Step2
DATA_DIR=/diazlab/data3/abhinav/resource/alphafold_database/                                                               # Please replace YourPath with the exact path on your machine {Your downloaded data directory}. # Path to the directory containing the AlphaFold 2 downloaded data.
FASTA_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/fasta_2/ # Please replace YourPath with the exact path on your machine. # Path to the directory containing the input FASTA files.
OUTPUT_DIR=/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/DeepAIR/preprocessing_structure_feature/sampledata/        # Please replace YourPath with the exact path on your machine # Path to the directory where the results will be saved.

for i in $(ls $FASTA_DIR); do
  echo $FASTA_DIR$i
  bash /diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/alphafold-2.3.1/run_alphafold.sh -f $FASTA_DIR$i \
    -t 2024-12-06 \
    -d $DATA_DIR \
    -o $OUTPUT_DIR
done
