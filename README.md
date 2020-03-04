# cgMLST-Pa

cgMLST to *Pseudomonas aeruginosa*.

## chewBBACA

* https://github.com/B-UMMI/chewBBACA_tutorial

## Workflow

* Step 1: Schema Creation
* Step 2: Allele calling
* Step 3: Schema Validation (Allele call)
* Step 4: Extracting the Matrix loci
* Step 5: Evaluation of the schema cgMLST
* Step 6: Analyze the proteins in the genes of the wgMLST
 
## Softwares and Downloads (Main dependencies)

* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* Prodigal 2.6.0 https://github.com/hyattpd/prodigal/releases/ or above

## Other dependencies (for schema evaluation only):

* ClustalW2 http://www.clustal.org/download/current/
* MAFFT https://mafft.cbrc.jp/alignment/software/linux.html
* DATAMASH https://www.gnu.org/software/datamash/
* MLST https://github.com/sanger-pathogens/mlst_check

## Step 1: Schema Creation

## Selection of complete genomes 

For *Pseudomonas aeruginosa* select the option RefSeq from GenBank at https://www.ncbi.nlm.nih.gov/assembly. RefSeq corresponds to a comprehensive, non-redundant, well-annotated set of reference sequences. A set of 142 complete genomes, 25 sequences at chromosome level, 1709 at contig level, and 1025 scaffold sequences of P. aeruginosa were publicly available in GenBank (https://www.ncbi.nlm.nih.gov/assembly) in September 2018. 

Genomes that could not be assigned MLST type sequence (ST) (https://github.com/sanger-pathogens/mlst_check) were not included for the construction of the cgMLST scheme. Of the 2901 available genomes, it was possible to assign STs to 2828 genomes. A second filter was added to remove unfinished genomes that had ≥200 contigs and 502 genomes were removed. 

In the end, 73 genomes were removed due to the absence of MLST *loci* and 502 were removed because the available sequences consisted of ≥200 contigs. Thus, of the 2901 genomes obtained from NCBI, 2326 genomes were used for the construction and validation of the cgMLST scheme, (142 complete sequences and 2184 unfinished genomes). 

Among the 142 genomes, *Pseudomonas aeruginosa* PAO1 reference genome (GCF_000006765.1) was included so that the Prodigal algorithm could use it as reference to recognize coding sequences (CDs). Prodigal generated the PAO1.trn file at this step. 

**The PAO1 genome was then removed from further analysis**.

## Step 1: Definition of CD sequences 

## Command: 

```bash
# create schema
chewBBACA.py CreateSchema -i Genomes20181011Sem_Plasmideo141/ --cpu 15 -o schema_seed --ptf PAO1.trn
```

**Note:** Genomes20181011Sem_Plasmideo141: Folder containing the 141 complete genomes that created the scheme.

## Step 2: Allele calling

In this step the allele calling is performed using the resulting set of *loci* found in step 1.

## Command: 

```bash
# run allelecall
chewBBACA.py AlleleCall -i Genomes20181011Sem_Plasmideo141 -g schema_seed/ -o results_cg --cpu 15 --ptf PAO1.trn
```

## Step 2.1: Paralog detection

In this step *loci* considered paralogous from result of the allelecall (see above) are removed

## Command: 

```bash
# run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_20190921T183955/results_alleles.tsv -g results_cg/results_20190921T183955/RepeatedLoci.txt -o alleleCallMatrix_cg
```

## Step 2.2: Genomes Quality Control

In this step we define a *Threshold* of the scheme that limits the loss of *loci* targets defined in the previous steps while excluding genomes considered to be of poor quality due to significant *loci* absence. 

We then define the percentage of *loci* that will constitute the scheme based on how many targets we want to keep in this phase:  **100%, 99.5%, 99% and 95%** of the *loci* present in the set of high quality genomes. This is one of the main steps in defining the cgMLST schema targets.

## Commands:

```bash
# run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```

The list of low qualiy genomes will then be removed from the original list using

```bash
# run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_120 -p 1.0 -g GenomeRemoved120thr.txt
```

This script selects all * loci * present in the selected * Threshold *. The value * p * is associated with the percentage of * loci * that has been set, for example: * p1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_120 folder is created which receives the result of the allelic profile for each of the cgMLST candidate * loci * (allelic profile matrix). The file in this folder (cgMLST.tsv) contains the allelic profile of each selected * loci * and will be used to create the core gene list.


## Step 2.3: Creating the core gene list

This command selects all target genes from the "cgMLST.tsv" spreadsheet.

```bash
# 10 list
head -n 1 cgMLST.tsv > Genes_100%_Core_120.txt
```
This list needs to be transposed so that each core gene name is reported in a single line:

```bash
# transpose table
datamash -W transpose < Genes_100%_Core_120.txt > Genes_Core_Al.txt 
```

This step generated the file> Genes_Core_Al.txt 

This list was then modified so that each name was preceeded by *schema_seed* :

```bash
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("schema_seed/"$1)}' > listgenes_core_100_120ca%.txt
```

## Step 3: Schema Validation (Allele calling)

In this step we repeat the allele calling using only the selected candidate *loci* for each of the genomes selected for validation (2184 genome drafts)

```bash
chewBBACA.py AlleleCall -i GenomasValidacao210919 -g listgenes_core_100_120ca%.txt -o results --cpu 15 --ptf PAO1.trn
```

This folder generated from this step **GenomasValidacao210919** has all 2184 validation draft genomes acquired from the NCBI that has ST and and they had less than 200 contigs.

## Step 3.1: Concatenate the allelic

Concatenate the allelic profile matrix obtained from the creation of the scheme with the matrix obtained for the validation genomes

To concatenate the matrix of the *loci* that defined the scheme and matrix of the loci of the validation genomes was used the following scripts:

```bash
# create header
head -n 1 cgMLST_120/cgMLST.tsv > cgMLST_all.tsv
```

```bash
# concatenate
grep -v ^FILE cgMLST_120/cgMLST.tsv results/ results_20190922T222448/results_alleles.tsv >> cgMLST_all.tsv
```

## Step 3.2: Evaluation of genome quality

After concatenation, the *TestGenomeQuality* to assess the impact of each validation genome in relation to the *loci* candidates in order to exclude low quality validation genomes. In this step you need to define a new *Threshold* for this dataset, as well as a new value of the parameter *p*, because loci that remain after the filters are the ones that constituted the final scheme.

```bash
 chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```

In order to exclude validation genomes that have left the scheme it is necessary to follow the steps described in **Step 2.2**

## Step 4: Extracting the Matrix loci

In this stage we chose to choose the *loci* present in 99% (*p0.99*) of the validation genomes and the *Threshold* 200 that limited the loss of the *loci* in the genomes. In this *Threshold* (200) 5 unfinished genomes were removed due to loss of *loci* targets.

To transpose (put the names of the validation genomes one in each line) I used the datamash and created the file > removedGenomes200thr.txt.

```bash
# transpose
datamash -W transpose < removedGenomes200.txt > removedGenomes200thr.txt
``` 
The genomes that were excluded in the Threshold 200 have been placed in the **removedGenomes200thr.txt**

```bash
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_200 -p0.99 -g removedGenomes200thr.txt 
```

This script selects *loci* and genomes that remained in the *Threshold* 200 and excludes the validation and *loci* genomes that were excluded in this *Threshold*.

## Step 5: Evaluation of the schema cgMLST

To assess the variability of the *loci* targets of cgMLST as well as the quality of the *loci* we run this script and graphically visualize the data.

```bash
chewBBACA.py SchemaEvaluator -i schema_seed/ -l rms/RmS.html -ta 11 --title "cgMLST custom r sales" --cpu 6
```

## Step 6: Analyze the proteins in the genes of the wgMLST

To check which protein encodes each loci found in the wg/cgMLST.

```bash
chewBBACA.py UniprotFinder -i schema_seed/ -t proteinID_Genome.tsv --cpu 10
```
