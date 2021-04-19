# Files for POSTAR3 database

### CLIPdb (RBP)

| Table name in MySQL          | content                                                      | position                                                     | species                                                      | note                 |
| ---------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------------------- |
| {species}_RBPanno            | RBP annotation, with motif from Pfam                         | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/RBP/RBP_anno/{species}_RBP_anno.txt` | all 7 (arabidopsis, fly, human, mouse, worm, yeast, zebrafish) |                      |
| {species}_GO                 | GO annotation for RBP                                        | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/RBP/RBP_GO/{species}_RBP_GO.txt` | all 7                                                        |                      |
| {species}_gene_map           | translate different gene name into what is used in database  | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/RBP/RBP_gene_map/{species}_gene_map.txt` | all 7                                                        |                      |
| GO                           | GO annotation for RBP                                        | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/RBP/GO/GO.txt` | all 7 in one file                                            | not used in database |
| pathway                      | pathway information      for RBP                             | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/RBP/pathway/pathway.txt` | all 7 in one file                                            | not used in database |
| {species}_RBP_clipdb_circRNA | RBP binding sites on circRNA junction from ENCODE and CircInteractome | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/circRNA/circRNA_CI_summary.txt`<br>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/circRNA/circRNA_ENCODE_summary.txt` | human                                                        |                      |



### RBS (RBP Binding Sites)

| Table name in MySQL                                          | content                                                      | position                                                     | species      | note                       |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------ | -------------------------- |
| {species}_clipdb_exp                                         | binding site records from CLIP-seq experiments (raw CLIP-seq data, ENCODE eCLIP, and PIP-seq from PMID 24393486) | `/data/zhaoweihao/project/POSTAR3/human_clipdb_exp.txt` (pcluster, human)<br>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/clipdb_eclip/mouse_clipdb_exp.txt` (BioII, mouse)<br>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/clipdb_eclip/{species}_clipdb.txt` (BioII, other 5 species) | all 7        | stored in different places |
| {species}_clipdb_pred                                        | predicted binding site records in POSTAR (FIMO, TESS, DeepBind) | `/data/zhaoweihao/project/POSTAR3/human_clipdb_pred.txt` (pcluster, human)<br>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/clipdb_eclip/mouse_clipdb_pred.txt` (BioII, mouse) | human, mouse |                            |
| {species}_gene                                               | gene information for RNAs in POSTAR                          | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/RNA_gene_info/{species}_gene.txt` | all 7        |                            |
| CancerGene<br>CoreTF<br>Disease<br>DiseaseGene<br>Drug<br>SpecificGene | gene annotation from various resource                        | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/RNA_anno/human_*.txt` | human        |                            |
| {species}_geneExp                                            | gene expression value (fpkm) of different species from ExpressionAtlas | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/geneExp/{species}_geneExp.txt` | all 7        |                            |
| {species}_RBP_hotspot                                        | RBP binding records on 20nt bin across genome                | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/binding_sites/RBP_hotspot/{species}_RBP_hotspot.txt` | all 7        |                            |
| {species}_circRNAanno                                        | circRNA annotation information from circBase                 | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/circRNA/anno/human_circRNAanno.txt` | human        |                            |
| {species}_circRNAexp                                         | circRNA expression information from circBase                 | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/circRNA/anno/human_circRNAexp.txt` | Human        |                            |



### RNA crosstalk

| Table name in MySQL                    | content                                                      | position                                                     | species                | note |
| -------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ---------------------- | ---- |
| {species}_RBP_clipdb_exp_miRNApredict  | RBS records (experimental) overlapped with miRNA binding predicted by different algorithms (miRanda, RNAhybrid, psRobot, psRNAtarget) | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/miRNA/{species}_RBP_clipdb_exp_miRNApredict.txt` (human, mouse)<br>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/miRNA/{species}_RBP_clipdb_miRNApredict.txt` (other species) | all 7 except zebrafish |      |
| {species}_RBP_clipdb_exp_miRNAvalidate | RBS records (experimental) overlapped with miRNA binding validated by AGO2 CLIP | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/miRNA/{species}_RBP_clipdb_exp_miRNAvalidate.txt` | human, mouse           |      |
| {species}_RBP_clipdb_exp_RNAediting    | RBS records overlapped with RNA editing sites from RADAR, DARNED, PMID 25373143, and GTEx | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/editing/{species}_RBP_clipdb_exp_RNAediting.txt` (human, mouse)<br/>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/editing/{species}_RBP_clipdb_RNAediting.txt` (other species) | all 7                  |      |
| {species}_RBP_clipdb_exp_RNAmod        | RBS records overlapped with RNA modification sites from DMBase2 | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/modification/{species}_RBP_clipdb_exp_RNAmod.txt` (human, mouse)<br/>`/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/crosstalk/modification/{species}_RBP_clipdb_RNAmod.txt` (other species) | all 7                  |      |



### Genomic Variants

| Table name in MySQL          | content                                                      | position                                                     | species | note |
| ---------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------- | ---- |
| {species}_RBP_clipdb_exp_SNP | RBS records overlapped with SNP annotation from dbSNP        | `/data/zhaoweihao/project/POSTAR3/database/variation/SNP/{species}_RBP_clipdb_exp_SNP.txt` (pcluster) | all 7   |      |
| human_RBP_clipdb_exp_1KG     | RBS records overlapped with SNP annotation from 1000 Genomes | `/data/zhaoweihao/project/POSTAR3/database/variation/1KG/human_RBP_clipdb_exp_1KG.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_gnomAD  | RBS records overlapped with SNP annotation from Genome Aggregation database (gnomAD) | `/data/zhaoweihao/project/POSTAR3/database/variation/gnomAD/human_RBP_clipdb_exp_gnomAD.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_eQTL    | RBS records overlapped with eQTL from GTEx                   | `/data/zhaoweihao/project/POSTAR3/database/variation/eQTL/human_RBP_clipdb_exp_eQTL.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_sQTL    | RBS records overlapped with sQTL from GTEx                   | `/data/zhaoweihao/project/POSTAR3/database/variation/sQTL/human_RBP_clipdb_exp_sQTL.txt` (pcluster) | human   |      |



### Disease Mutations

| Table name in MySQL            | content                                             | position                                                     | species | note |
| ------------------------------ | --------------------------------------------------- | ------------------------------------------------------------ | ------- | ---- |
| human_RBP_clipdb_exp_CancerWES | RBS records overlapped with PanTCGA exome mutations | `/data/zhaoweihao/project/POSTAR3/database/disease/CancerWES/human_RBP_clipdb_exp_CancerWES.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_CancerWGS | RBS records overlapped with PCAWG genome mutations  | `/data/zhaoweihao/project/POSTAR3/database/disease/CancerWGS/human_RBP_clipdb_exp_CancerWGS.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_CCLE      | RBS records overlapped with mutations from CCLE     | `/data/zhaoweihao/project/POSTAR3/database/disease/CCLE/human_RBP_clipdb_exp_CCLE.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_ClinVar   | RBS records overlapped with mutations from ClinVar  | `/data/zhaoweihao/project/POSTAR3/database/disease/ClinVar/human_RBP_clipdb_exp_ClinVar.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_COSMIC    | RBS records overlapped with mutations from COSMIC   | `/data/zhaoweihao/project/POSTAR3/database/disease/COSMIC/human_RBP_clipdb_exp_COSMIC.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_denovodb  | RBS records overlapped with mutations from denovodb | `/data/zhaoweihao/project/POSTAR3/database/disease/denovodb/human_RBP_clipdb_exp_denovodb.txt` (pcluster) | human   |      |
| human_RBP_clipdb_exp_GWASdb    | RBS records overlapped with mutations from GWASdb2  | `/data/zhaoweihao/project/POSTAR3/database/disease/GWASdb2/human_RBP_clipdb_exp_GWASdb.txt` (pcluster) | human   |      |



### Structurome

| Table name in MySQL         | content                                                      | position                                                     | species                      | note                                                         |
| --------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ---------------------------- | ------------------------------------------------------------ |
| structurome_{species}       | structurome information of each RBS, including its binding site information, structurome reactivity, sequence, and secondary structure model predicted by RNAfold/Fold | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/Structure/Fold/{species}/{species}_structurome_summary.txt` | all 7 (human in calculation) | file only contains structurome information, data table in MySQL database was generated by joining this table with `{species}_clipdb_exp` table using the following MySQL command |
| strcuturome\_{species}\_all | like above, but contain every information (structurome_human only has one record for one RBP binding site) |                                                              | human                        |                                                              |

```mysql
CREATE TABLE temp
(
datasetName char(40),
myIDRBP char(60),
protocolStructurome char(30),
softwareStructurome char(10),
sequence text,
reactivity text,
structure text,
mfe double,
tissueTypeStructurome char(20),
conditionStructurome char(30),
INDEX myIDRBP (myIDRBP)
) ENGINE=MyISAM;

LOAD DATA LOCAL INFILE 'temp.txt' INTO TABLE temp;

CREATE TABLE structurome_human
(
datasetName char(40),
chr char(5),
start int(11),
end int(11),
myIDRBP char(60),
sthRBP int(5),
strandRBP char(3),
RBP char(10),
protocolCLIP char(30),
tissueTypeCLIP char(20),
sourceCLIP char(200),
score double,
chrGene char(5),
startGene int(11),
endGene int(11),
geneID char(30),
geneType char(20),
geneName char(30),
PhastCons double,
Phylop double,
l1 char(40),
protocolStructurome char(30),
softwareStructurome char(10),
sequence text,
reactivity text,
structure text,
mfe double,
tissueTypeStructurome char(20),
conditionStructurome char(30),
INDEX datasetName (datasetName),
INDEX protocolCLIP (protocolCLIP),
INDEX geneName (geneName),
INDEX protocolStructurome (protocolStructurome),
INDEX softwareStructurome (softwareStructurome)
) ENGINE=MyISAM;

insert into structurome_human select b.datasetName as datasetName, a.chr as chr, a.start as start, a.end as end, b.myIDRBP as myIDRBP, a.sthRBP as sthRBP, a.strandRBP as strandRBP, a.RBP as RBP, a.protocol as protocolCLIP, a.tissueType as tissueTypeCLIP, a.source as sourceCLIP, a.score as score, a.chrGene as chrGene, a.startGene as startGene, a.endGene as endGene, a.geneID as geneID, a.geneType as geneType, a.geneName as geneName, a.PhastCons as PhastCons, a.Phylop as Phylop, a.l1 as l1, b.protocolStructurome as protocolStructurome, b.softwareStructurome as softwareStructurome, b.sequence as sequence, b.reactivity as reactivity, b.structure as structure, b.mfe as mfe, b.tissueTypeStructurome as tissueTypeStructurome, b.conditionStructurome as conditionStructurome from human_clipdb_exp as a, temp as b where a.myIDRBP = b.myIDRBP;

```





### Translatome

| Table name in MySQL                       | content                                                      | position                                                     | species                | note                     |
| ----------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ---------------------- | ------------------------ |
| translatome\_{species}\_DE                | Ribo-seq density on each ORF under every condition           | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/translatome/translatome_density/translatome_{species}_density.txt` | all 7 except zebrafish |                          |
| translatome\_{species}\_TE                | translation efficiency of each ORF under every condition     | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/translatome/translatome_TE/translatome_{species}_TE.txt` | all 7 except zebrafish |                          |
| translatome\_{species}\_TS                | translation score (Ribotaper, ORFscore, RibORF) of each ORF under every condition | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/translatome/translatome_TS/translatome_{species}_TS.txt` | all 7 except zebrafish |                          |
| translatome\_{species}\_heatmap_den       | Ribo-seq density on each ORF under every condition (raw, used to draw heatmap on web page) | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/translatome/translatome_heatmap_den/translatome_{species}_heatmap_den.txt` | all 7 except zebrafish |                          |
| translatome\_{species}\_name              | ORF information, including transcriptID, geneID, ORFID, gene name, and ORF category | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/translatome/translatome_name/translatome_{species}_name.txt` | all 7 except zebrafish |                          |
| translatome\_{species}\_bedgraph          | Ribo-seq signal of each ORF (raw)                            |                                                              |                        | not calculated by Ziyuan |
| translatome\_{species}\_denoised_bedgraph | Ribo-seq signal of each ORF (denoised)                       |                                                              |                        | not calculated by Ziyuan |



### Degradome

| Table name in MySQL | content                                                      | position                                                     | species                        | note |
| ------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------ | ---- |
| degradome_{species} | degradome information, including binding information of miRNA and target RNA | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/database/degradome/degradome_{species}.txt` | arabidopsis, fly, human, mouse |      |
| {species}_miRNAanno | miRNA annotation from miRBase                                | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/miRNA/{species}_miRNAanno.txt` | arabidopsis, fly, human, mouse |      |
| {species}_miRNAGO   | miRNA GO information from Ensembl                            | `/BioII/lulab_b/zhaoweihao/project/POSTAR3/miRNA/{species}_miRNAGO.txt` | human                          |      |

