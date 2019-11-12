# cgMLST-Pa

Script utilizado para realizar o esquema de cgMLST com o chewBBACA



## chewBBACA

* https://github.com/B-UMMI/chewBBACA_tutorial

 
## Softwares e Downloads

* NCBI complete genomes?


## Step 1: Schema creation

To train Prodigal in recognizing the coding sequences (CDs) of the *Pseudomonas aeruginosa* genome it was trained with the PAO1 reference genome (GCF_000006765.1) that generated the PAO1.trn file. 

The genome used for training was removed from the analysis.


```bash
# create schema
chewBBACA.py CreateSchema -i Genomes20181011Sem_Plasmideo141/ --cpu 15 -o schema_seed --ptf PAO1.trn
```

**Note:** Genomes20181011Nem_Plasmideo141: Folder containing the complete NCBI genomes that created the schema.



## Step 2: Allele calling

In this step the alleles are assigned to each of the CDs found in the breeding genomes.

```bash
# run allelecall
chewBBACA.py AlleleCall -i Genomes20181011Sem_Plasmideo141 -g schema_seed/ -o results_cg --cpu 15 --ptf PAO1.trn
```


## Step 2.1: Paralog detection

In this step the script removes the *loci* considered parallogues from the result of the previous script.

```bash
# run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_20190921T183955/results_alleles.tsv -g results_cg/results_20190921T183955/RepeatedLoci.txt -o alleleCallMatrix_cg
```


## Step 2.2: Genomes Quality Control

In this step it is necessary to define the **Threshold** of the Scheme that limits the loss of *loci* targets of the Schema Creation genomes and excludes genomes considered to be of poor quality due to the *loci*  loss. 

The next step is to define the percentage of *loci* that will constitute the Scheme, this can be **100%, 99.5%, 99% and 95%** of the *loci* present in the breeding genomes. This is one of the main steps in defining the cgMLST schema targets.



```bash
# run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```



## Step 2.3: Extracting Matrix Loci

In this step we use the list of **removedGenomes** which is output from the **TestGenomeQuality**. In this step we select the Threshold chosen from the removedGenomes.txt file, for example, Threshold 120 and select all the genomes that came out of this Threshold and create a new file named, for example, GenomeRemoved120thr.txt. The name of the genomes that were obtained from the removedGenomes.txt file must be one on each line for use in the next step.



```bash
# run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_120 -p 1.0 -g GenomeRemoved120thr.txt
```


This script selects all * loci * present in the selected * Threshold *. The value * p * is associated with the percentage of * loci * that has been set, for example: * p1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_120 folder is created which receives the result of the allelic profile for each of the cgMLST scheme * loci * that will be in the allelic profile matrix. The file in this folder (cgMLST.tsv) contains the allelic profile of each selected * loci * and will be used to create the core gene list.


## Step 2.4 :Creating the core gene list

This command selects all target genes from the "cgMLST.tsv" spreadsheet that corresponds to the first line of this file. Then we perform the script below:

```bash
# 10 list
head -n 1 cgMLST.tsv > Genes_100%_Core_120.txt
```



To transpose (put the names of the core genes on each line) I used the datamash and created the file> Genes_Core_Al.txt which is the previous script.

```bash
# transpose table
datamash -W transpose <Genes_100%_Core_120.txt > Genes_Core_Al.txt 
```

This print schema_seed / command on each cgMLST target gene for use in the validation step.

```bash
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("schema_seed/"$1)}' > listgenes_core_100_120ca%.txt
```

# Schema Validation


## Step 3: Schema Validation (Allele call)



```bash
chewBBACA.py AlleleCall -i GenomasValidacao210919 -g listgenes_core_100_120ca%.txt -o results --cpu 15 --ptf PAO1.trn
```

 

\##**Comentário:** Este comando identifica os *loci* candidatos em cada um dos genomas de validação e atribui um perfil alélico para cada um deles pelo AlleleCall.



\##**Comentário:** Esta pasta “GenomasValidacao210919” possui todos os genomas de validação adquiridos no NCBI que possui ST.



\##**Parte 3.1 -** **Concatenar a matrix do perfil alélico obtido da criação do esquema com a matrix obtida para os genomas de validação**



\##**Comentário:** Para concatenar a matrix dos *loci* que definiu o esquema e a matrix dos loci dos genomas de validação foi utilizado os scripts abaixo:

```bash
head -n 1 cgMLST_120/cgMLST.tsv > cgMLST_all.tsv
```



```bash
grep -v ^FILE cgMLST_120/cgMLST.tsv results/ results_20190922T222448/results_alleles.tsv >> cgMLST_all.tsv
```



## Step 3.2:  Avaliação da qualidade dos genomas





```bash
 chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```



\##**Comentário:** Após a concatenação é realizado o *TestGenomeQuality* para avaliar o impacto de cada genoma de validação em relação aos *loci* candidatos visando excluir genomas de validação de baixa qualidade. Nesta etapa é necessário definir um novo *Threshold* para este conjunto de dados, bem como um novo valor do parâmetro *p*, pois os loci que ficarem após os filtros são aqueles que constituíram o esquema final. 



\##**Comentário: Para excluir os genomas de validação que saíram do esquema é necessário seguir os passos descritos na Parte 2.2**



## Step 4:  Extraindo os loci da Matrix



\##**Comentário:** Nesta etapa optamos em escolher os *loci* presentes em 99% (*p0.99*) dos genomas de validação e o *Threshold* 200 que limitou a perda dos *loci* nos genomas.



\##**Comentário:** Os genomas que foram excluídos no *Threshold* 200 foram colocados no arquivo removedGenomes200thr.txt 



```bash
datamash -W transpose < removedGenomes200.txt > removedGenomes200thr.txt
```



\##**Comentário:** Para transpor (colocar os nomes dos genomas de validação um em cada linha) utilizei o datamash e criei o arquivo > removedGenomes200thr.txt.

 

```bash
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_200 -p0.99 -g removedGenomes200thr.txt 
```



\##**Comentário:** Este script seleciona os *loci* e genomas que permaneceram no *Threshold* 200 e exclui os genomas de validação e loci que foram excluídos neste *Threshold*.



\##**Avaliação do esquema de cgMLST** (**Schema Evaluator was run on the cgMLST schema):**



```bash
chewBBACA.py SchemaEvaluator -i schema_seed/ -l rms/RmS.html -ta 11 --title "cgMLST schema GBS tutorial schema evaluator" --cpu 6
```



\##**Comentário:** Para avaliar a variabilidade dos *loci* alvos do cgMLST bem como a qualidade dos *loci* selecionados rodamos este script e visualizamos graficamente os dados.



\##**Analisar as proteínas dos genes do wgMLST**

```bash
chewBBACA.py UniprotFinder -i schema_seed/ -t proteinID_Genome.tsv --cpu 10
```



\##**Comentário:** Para verificar qual proteína codifica cada loci encontrado do wg/cgMLST.
