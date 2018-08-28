BridGE (Bridging Genes with Epistasis)

1. Set up environment variables
 
To make BridGE run properly, BridGE and its sub directories need to be added to PATH and MATLABPATH. Depending on the shell, this can be done with the following commands:
source setup_env_csh.sh # c shell
source setup_env_bash.sh # bash shell

2. BridGE components and examples

BridGE has four main components: 1) data processing (DataProcess); 2) compute SNP-SNP interaction (ComputeInteraction); 3) sample permutation (SamplePermutation); 4) post analysis (Analysis). "bridge.sh" is the main script that is used to call each component (with the --job flag). All input parameters and their default values can be found by typing the following command:
bridge.sh --help

2.1. Processing data (DataProcess)
(1) In addition to the standard gwas data quality control procedures, BridGE's data processing procedure also removes sex chromosomes, filters out related samples, matches case and control samples (optional), excludes SNPs that can't be mapped to the given gene sets, converts data from plink format to matlab, and binarizes the SNP data if it is running with recessive or dominant assumptions.

Example 2.1.1, run data processing using example yeast genetic data: 
bridge.sh --projectDir=$BRIDGEPATH/example_project --job=DataProcess --plinkFile=gwas_example --geno=0.05 --mind=0.05 --hwe=0 --matchCC=0 --genesets=$BRIDGEPATH/refdata/Costanzo2016_S12_complexes.mat --geneAnnotation=$BRIDGEPATH/refdata/glist_yeast_orfs.txt --minPath=10 --maxPath=300 --mappingDistance=500 --ldR2=0.2 --pihat=1 

Example 2.1.2, run data processing with all default parameters (don't test this with given example file):
bridge.sh --job=DataProcess --projectDir=<your project directory> --plinkFile=<your data in plink format without the extension>

Example 2.1.3, directly call the data processing script processdata.sh (don't test this with given example file):
processdata.sh --plinkFile=<your data in plink format without the extension>

Example 2.1.4, use pre-processed data that is already in MATLAB format and map SNPs to pathways:
bridge.sh --projectDir=$BRIDGEPATH/example_project --job=DataProcess --matlabFile=gwas_example.mat --genesets=$BRIDGEPATH/refdata/Costanzo2016_S12_complexes.mat --geneAnnotation=$BRIDGEPATH/refdata/glist_yeast_orfs.txt --minPath=10 --maxPath=300 --mappingDistance=500

(2) "--projectDir=<your project directory>" can be skipped if you are running these commands under your project directory. 

(3) Many input parameters used in this data processing step, such as "--mind", "--geno", "--maf", and "--hwe" are not used in later steps. Also, "--genesets" and "--geneAnnotation" are used to map SNPs to pathways. Once the mapping file is generated, these parameters are also no longer needed in the pilot run and sample permutation steps, but they will be used in the post analysis step to map SNPs back to genes.

2.2 Compute SNP-SNP interaction (ComputeInteraction)
(1) SNP-SNP interaction can be computed independently or in SamplePermutation. One important note for this step is that, to improve the computational efficiency for hygeSSI interaction, a dictionary table will be computed during first time SNP-SNP hygeSSI interaction computation. Don't launch several hygeSSI computation in the same folder before the dictionary table is computed.

(2) Supported interaction type and corresponding interaction measurement:
AA: additive-additive interaction, measured by logistic regression using software cassi
RR: recessive-recessive interaction, measured by hygeSSI
DD: dominant-dominant interaction, measured by hygeSSI
RD: recessive-dominant interaction, measured by hygeSSI
combined: RR/DD/RD combined interaction 

Example 2.2.1, compute hygeSSI based RR/DD/RD combined interaction using example yeast genetic data
bridge.sh --projectDir=$BRIDGEPATH/example_project --job=ComputeInteraction --nWorker=10 --model=combined --samplePerms=10

Example 2.2.2, compute lrSSI based AA interaction using example yeast genetic data
bridge.sh --projectDir=$BRIDGEPATH/example_project --job=ComputeInteraction --model=AA --samplePerms=10

2.3 Run sample permutation (SamplePermutation)
(1) BridGE can be used to test pathway-level enrichment based on binarized network or weighted network. 

Example 2.3.1, run BridGE with weighted combined interaction network using example yeast genetic data
bridge.sh --job=SamplePermutation --projectDir=$BRIDGEPATH/example_project --interaction=hygeSSI --model=combined --binaryNetwork=0 --samplePerms=10 --snpPerms=1000 --minPath=10

Example 2.3.2, run BridGE with binary AA interaction network using example yeast genetic data
bridge.sh --job=SamplePermutation --projectDir=$BRIDGEPATH/example_project --model=AA --binaryNetwork=1 --samplePerms=10 --snpPerms=1000 --minPath=10

(2) Indiviudal sample permutation runs can be called by "runsampleperm.sh" instead of using the master BridGE script (bridge.sh).

Example 2.3.2, run sample permutation with RR interaction using example yeast genetic data, and with real pheynotype data (randRun=0):
runsampleperm.sh --projectDir=$BRIDGEPATH/example_project --minPath=10 --model=RR  --randRun=0

2.4 Run analysis based on permutation results
This step will compute false discovery rates, and write significant discoveries (with detailed statistical information and gene drivers) to an excel file.

Example 2.4.1, run analysis on permutation results based on combined interaction: 
bridge.sh --job=Analysis --projectDir=$BRIDGEPATH/example_project --model=combined  --samplePerms=10 --minPATH=10 --pvalueCutoff=0.05 --mappingDistance=500

Inputs for --model --samplePerms --minPATH --mappingDistance  should also be specified (unless they are the same as the default values).
