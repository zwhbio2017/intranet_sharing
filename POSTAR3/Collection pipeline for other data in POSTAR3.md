# Collection pipeline for other data in POSTAR3

### CLIPdb (RBP)

#### {species}_RBPanno

**Note**: 

1. RBP geneName needs to be unified (in some paper they use specific names that do not fit ensembl, please refer to ensembl for standard geneName)
2. geneID was also from ensembl
3. protein domain (motif): obtain sequence from ncbi protein, then input protein sequence into http://pfam.xfam.org/ -> SEQUENCE SEARCH, only curate significant match



#### {species}_GO

**Note**:

1. Input RBP geneName/geneID into ensembl BioMart, search Gene name, Gene stable ID, GO term accession, GO term name, and GO domain



#### {species}_gene_map

**Note**:

1. gene synonyms from various resource (does not need to change during collection)



#### {species}_RBP_clipdb_circRNA

**Note**:

1. Collected from PMID33287884 (ENCODE) and PMID26669964 (CircInteractome), only contains human circRNA binding data, could be updated if corresponding resource is updated, or new resource, especially with new species, have emerged.
2. CircInteractome used hg19 annotation, and some of the genome position was relative to the circRNA, remember to use liftOver to change the position into hg38 coordinate and translate relative position to absolute position



#### Motif

**Software**: Homer, MEME (+ FIMO), RNAcontext, RNApromo

**Steps**:

1. Prepare data and folders
2. Split Piranha binding sites for each RBP
3. Random select binding site records as training set and test set
4. Prepare fasta file for motif discovery
5. Run corresponding software
6. Summarize plot for each motif

**Code**:

```bash
# 1. Prepare binding site called by Piranha ({species}_Piranha_peak.txt) and the total list of RBP (gene_list.txt) that needs to calculate motif

# Folders: 

# 1.peak_all
# 2.peak_training
# 3.peak_test
# 4.MEME
# 5.RNApromo
# 6.RNAcontext
# 7.Homer_training
# 8.Homer_test
# 9.FIMO
# 10.PWM_MEME
# 11.PWM_Homer
# 12.PWM_RNAcontext
# 13.figure_RNApromo

# 2. split Piranha binding sites for each RBP

perl split.pl
```

```perl
# split.pl (split each binding site into corresponding file for each RBP in 1.peak_all/)

use strict;

open INPUT1,"gene_list.txt" or die;

while(<INPUT1>){
        chomp;
        my $RBP = $_;
        open INPUT2,"Fly_Piranha_peak.txt" or die;
        open OUTPUT,">./1.peak_all/$RBP.all_peak.txt" or die;
        while(<INPUT2>){
                chomp;
                my @line = split /\t/,$_;
                if($RBP eq $line[6]){
                        print OUTPUT "$_\n";
                }
        }
}
```

```bash
# 3. Select training set and test set for motif discovery

for i in `ls | grep txt`
do
		perl run.pl $i
done		
```

```perl
# run.pl (split into training and test set)

use strict;

#Usage: for i in `ls | grep txt`;do perl run.pl $i;done

my $total_line = `wc -l $ARGV[0]`;

chomp($total_line);

$ARGV[0] =~ /(.*)\.all_peak\.txt/;

my $file_name = $1;

if($total_line>=1000){
        `sort -k 11rn $ARGV[0] | head -n 500 > ../2.peak_training/$file_name.training_peak.txt`;
        `sort -k 11rn $ARGV[0] | head -n 1000 | tail -n 500 > ../3.peak_test/$file_name.test_peak.txt`;
}
else{
        my $cutoff = int($total_line/2);
        `sort -k 11rn $ARGV[0] | head -n $cutoff > ../2.peak_training/$file_name.training_peak.txt`;
        `sort -k 11rn $ARGV[0] | tail -n $cutoff > ../3.peak_test/$file_name.test_peak.txt`;
}
```

```bash
# 4. prepare fasta file for motif discovery

# MEME (Homer)

for i in `ls | grep txt`
do
		perl run_MEME.pl $i
done

# RNAcontext

for i in `ls | grep txt`
do
		perl run_RNAcontext.pl $i
done

# RNApromo

for i in `ls | grep txt`
do
		perl run_RNApromo.pl $i
done
```

```perl
# run_MEME.pl

use strict;

#usage: for i in `ls | grep txt`;do perl run.pl $i;done

open INPUT1,"$ARGV[0]" or die;

$ARGV[0] =~ /(.*)\.training_peak.txt/;

my $file_name = $1;

open TEMP,">rename.txt" or die;

my $i;

while(<INPUT1>){
        chomp;
        $i++;
        my @line = split /\t/,$_;
        print TEMP "$line[0]\t$line[1]\t$line[2]\t$i\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
}

open INPUT2,"rename.txt" or die;

my %strand;

while(<INPUT2>){
        chomp;
        my @line = split /\t/,$_;
        $strand{$line[3]} = $line[5];
}

# bedtools getfasta (need to change -fi according to the species)
`fastaFromBed -s -name -fi /BioII/lulab_b/zhaoweihao/database/Fly/dmel-all-chromosome-r6.18.fasta -bed rename.txt -fo temp.fa`;

open INPUT3,"temp.fa" or die;

open OUTPUT,">MEME/$file_name.training_peak.fa" or die;

my $j;

my $seq_name;

while(<INPUT3>){
        chomp;
        $j++;
        if(($j-1)%2==0){
                $_ =~ />(.*)/;
                $seq_name = $1;
                print OUTPUT ">$seq_name\n";
        }
        elsif(($j-2)%2==0){
                $_ =~ tr/ACGTacgt/ACGTACGT/;
                print OUTPUT "$_\n";
        }
}

unlink("rename.txt");
unlink("temp.fa");
```

```perl
# run_RNAcontext.pl

use strict;

#usage: for i in `ls | grep txt`;do perl run.pl $i;done

open INPUT1,"$ARGV[0]" or die;

$ARGV[0] =~ /(.*)\.training_peak.txt/;

my $file_name = $1;

open TEMP,">rename.txt" or die;

my $i;

# extend to 60nt
while(<INPUT1>){
        chomp;
        $i++;
        my @line = split /\t/,$_;
        if($line[2]-$line[1]<60){
                $line[1] = $line[1] - 20;
                if($line[1] < 0){
                        $line[1] = 0;
                }
                $line[2] = $line[2] + 20;
        }
        print TEMP "$line[0]\t$line[1]\t$line[2]\t$i\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
}

open INPUT2,"rename.txt" or die;

my %strand;

while(<INPUT2>){
        chomp;
        my @line = split /\t/,$_;
        $strand{$line[3]} = $line[5];
}

# bedtools getfasta (need to change -fi according to the species)
`fastaFromBed -s -name -fi /BioII/lulab_b/zhaoweihao/database/Fly/dmel-all-chromosome-r6.18.fasta -bed rename.txt -fo temp.fa`;

open INPUT3,"temp.fa" or die;

open OUTPUT,">RNAcontext/$file_name.training_peak.fa" or die;

my $j;

my $seq_name;

while(<INPUT3>){
        chomp;
        $j++;
        if(($j-1)%2==0){
                $_ =~ />(.*)/;
                $seq_name = $1;
                print OUTPUT ">$seq_name\n";
        }
        elsif(($j-2)%2==0){
                $_ =~ tr/ACGTacgt/ACGUACGU/; # replace T with U
                print OUTPUT "$_\n";
        }
}

unlink("rename.txt");
unlink("temp.fa");
```

```perl
# run_RNApromo.pl

use strict;

#usage: for i in `ls | grep txt`;do perl run.pl $i;done

open INPUT1,"$ARGV[0]" or die;

$ARGV[0] =~ /(.*)\.training_peak.txt/;

my $file_name = $1;

open TEMP,">rename.txt" or die;

my $i;

# extend to 60nt
while(<INPUT1>){
        chomp;
        $i++;
        my @line = split /\t/,$_;
        if($line[2]-$line[1]<=60){
                $line[1] = $line[1] - 20;
                if($line[1]<0){
                        $line[1] = 0;
                }
                $line[2] = $line[2] + 20;
        }
        print TEMP "$line[0]\t$line[1]\t$line[2]\t$i\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
}

open INPUT2,"rename.txt" or die;

my %strand;

while(<INPUT2>){
        chomp;
        my @line = split /\t/,$_;
        $strand{$line[3]} = $line[5];
}

# bedtools getfasta (need to change -fi according to the species)
`fastaFromBed -s -name -fi /BioII/lulab_b/zhaoweihao/database/Fly/dmel-all-chromosome-r6.18.fasta -bed rename.txt -fo temp.fa`;

open INPUT3,"temp.fa" or die;

open OUTPUT,">RNApromo/$file_name.training_peak.fa" or die;

my $j;

my $seq_name;

while(<INPUT3>){
        chomp;
        $j++;
        if(($j-1)%2==0){
                $_ =~ />(.*)/;
                $seq_name = $1;
                print OUTPUT ">$seq_name\n";
        }
        elsif(($j-2)%2==0){
                $_ =~ tr/ACGTacgt/ACGTACGT/;
                print OUTPUT "$_\n";
        }
}

unlink("rename.txt");
unlink("temp.fa");
```

```bash
# 5. Run algorithm

# MEME + FIMO

# MEME

cd /Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/4.MEME/AUB
source /BioII/lulab_b/containers/singularity/wrappers/bashrc

/Share2/home/zhaoweihao/app/meme_4.11.4/bin/meme \
/Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/2.peak_training/MEME/AUB.training_peak.fa \ # input training fasta file
-o output \ # output directory
-dna \ # input is DNA sequence
-minw 4 \ # shortest reported motif length
-maxw 10 \ # longest reported motif length
-nmotifs 25 # max number of reported motif

# FIMO

# replace with proper numbering system
sed 's/10.0e+000/1.0e+001/g' ../4.MEME/AUB/output/meme.txt | sed 's/10.0e+001/1.0e+002/g' | sed 's/10.0e+002/1.0e+003/g' | sed 's/10.0e+003/1.0e+004/g' | sed 's/10.0e+004/1.0e+005/g' | sed 's/10.0e+005/1.0e+006/g' | sed 's/10.0e+006/1.0e+007/g' | sed 's/10.0e+007/1.0e+008/g' | sed 's/10.0e+008/1.0e+009/g' | sed 's/10.0e+009/1.0e+010/g' | sed 's/10.0e+010/1.0e+011/g' > ./AUB/meme.txt

cd /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/9.FIMO/AUB/

/Share2/home/lulab/zhuyumin/share/zhuyumin/apps/meme_4.11.4/bin/fimo \
--thresh 0.01 \ # p value threshold
-o output \ # output directory
meme.txt \ # input MEME result file (modified)
/BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/3.peak_test/MEME/AUB.test_peak.fa # input test fasta file

# Homer

# Training

cd /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/7.Homer_training/AUB/

/BioII/lulab_b/zhaoweihao/app/homer/bin/findMotifs.pl \ /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/2.peak_training/MEME/AUB.training_peak.fa \ # input training fasta file
fasta \ # input format
output \ # output directory
-len 4,5,6,7,8,9,10 \ # possible motif length
-rna # look for RNA motif

# Test

cd /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/8.Homer_test/AUB/

/BioII/lulab_b/zhaoweihao/app/homer/bin/findMotifs.pl \
/BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/3.peak_test/MEME/AUB.test_peak.fa \ # input test fasta file
fasta \ # input format
output \ # output directory
-rna \ # look for RNA motif 
-find /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/7.Homer_training/AUB/output/homerMotifs.all.motifs \ # input training motif result file
> count.txt

# RNAcontext

# random.sh (prepare random background sequence for motif discovery)

# select random sequence to generate bed file
bedtools random \
-n 500 \ # total number of sequence
-l 60 \ # length for each sequence
-g /BioII/lulab_b/zhaoweihao/database/Fly/dmel-6.18.chrom.sizes \ # input genome size file (chromosome name \t chromosome size)
> random1.bed

bedtools random \
-n 500 \ # total number of sequence
-l 60 \ # length for each sequence
-g /BioII/lulab_b/zhaoweihao/database/Fly/dmel-6.18.chrom.sizes \ # input genome size file (chromosome name \t chromosome size)
> random2.bed

# get the random sequence fasta file
fastaFromBed \
-s \ # consider strand
-name \ # use bed file name as fasta header name
-fi /BioII/lulab_b/zhaoweihao/database/Fly/dmel-all-chromosome-r6.18.fasta \ # input genome fasta file
-bed random1.bed \ # input random bed file
-fo random1.fa # output sequence fasta file

fastaFromBed \
-s \ # consider strand
-name \ # use bed file name as fasta header name
-fi /BioII/lulab_b/zhaoweihao/database/Fly/dmel-all-chromosome-r6.18.fasta \ # input genome fasta file
-bed random2.bed \ # input random bed file
-fo random2.fa # output sequence fasta file

# pre.sh (prepare sequence file, 0 \t sequence)

for i in `cat ../gene_list.txt`
do
        mkdir ${i}/pre
        cat ../2.peak_training/RNAcontext/${i}.training_peak.fa | grep -v '>' | awk -v FS='\t' -v OFS='\t' '{print "0\t"$1}' > ${i}/pre/${i}.training_peak.sequences.txt
        cat ../3.peak_test/RNAcontext/${i}.test_peak.fa | grep -v '>' | awk -v FS='\t' -v OFS='\t' '{print "0\t"$1}' > ${i}/pre/${i}.test_peak.sequences.txt
done

# run.sh (run analysis for each RBP)

for i in `cat ../gene_list.txt`
do
        echo $i start at `date`
        cd /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/6.RNAcontext
        cp run.bsub ${i}/
        bash ${i}/run.bsub $i
        echo $i end at `date`
done    

# run.bsub (run analysis for a single RBP)

RBP=$1

# use sfold to prepare RNAcontext input
prepare_rnacontext(){
    n_threads=${n_threads:=16}
    if [ "$#" -ne 3 ];then
        echo "Usage: prepare_rnacontext sequence_file sprofile_file prefix"
        return 1
    fi
    sequence_file=${1:?}
    sprofile_file=${2:?}
    prefix=${3:?}
    rm -rf $prefix/sprofile
    [ -d $prefix/seq ] || mkdir -p $prefix/seq
    [ -d $prefix/sfold ] || mkdir  -p $prefix/sfold
    [ -d $prefix/sprofile ] || mkdir -p $prefix/sprofile
    # generate sfold input sequence
    awk -v output_dir=$prefix/seq \
        '{fname=output_dir "/" NR ".fa"; print ">" NR > fname; print $2 > fname}' $sequence_file
    n_seqs=$(ls $prefix/seq | wc -l)
    # run sfold
    {
    for seq_id in $(seq $n_seqs);do
        echo /BioII/lulab_b/zhaoweihao/app/sfold-2.2/bin/sfold -a 0 -o $prefix/sfold/$seq_id $prefix/seq/${seq_id}.fa
    done
    } | parallel -j $n_threads
    # run sprofile
    {
    for seq_id in $(seq $n_seqs);do
        sequence=$(sed '2 !d' $prefix/seq/${seq_id}.fa)
        echo /Share2/home/lulab/zhuyumin/share/zhuyumin/apps/sfold-2.2/sprofile \
            $prefix/sfold/$seq_id/sstrand.out \
            $prefix/sfold/$seq_id/loopr.out \
            $prefix/sprofile/${seq_id}.txt \
            $sequence
    done
    } | parallel -j $n_threads
    # cat sprofile output files
    {
    for seq_id in $(seq $n_seqs);do
        echo $prefix/sprofile/${seq_id}.txt
    done
    } | xargs cat > $sprofile_file
}

cd /BioII/lulab_b/zhaoweihao/project/POSTAR3/motif/Fly/update/6.RNAcontext/${RBP}

# merge training with random1, test with random2 to generate input sequence file
cat <(awk 'BEGIN{FS="\t";OFS="\t"}{print "1",$2}' pre/${RBP}.training_peak.sequences.txt) <(grep -v '>' ../random1.fa|awk 'BEGIN{FS="\t";OFS="\t"}{print "0",$0}'|sed 's/T/U/g'|tr a-z A-Z)|grep -v 'N' > ${RBP}.training_peak.sequences.txt
cat <(awk 'BEGIN{FS="\t";OFS="\t"}{print "1",$2}' pre/${RBP}.test_peak.sequences.txt) <(grep -v '>' ../random2.fa|awk 'BEGIN{FS="\t";OFS="\t"}{print "0",$0}'|sed 's/T/U/g'|tr a-z A-Z)|grep -v 'N' > ${RBP}.test_peak.sequences.txt

# prepare sprofile for training and test set
prepare_rnacontext "${RBP}.training_peak.sequences.txt" "${RBP}.training_peak.sprofile.txt" RBP.training_peak
prepare_rnacontext "${RBP}.test_peak.sequences.txt" "${RBP}.test_peak.sprofile.txt" RBP.test_peak

mkdir outputs
/BioII/lulab_b/zhaoweihao/app/RNAcontext/bin/rnacontext \
-w 4-10 \ # output motif length
-a ACGU \ # sequence base (RNA)
-e PLMU \ # annotation alphabet for structural motif
-s 5 \ # number of initializations or restarts
-c ${RBP}.training_peak.sequences.txt \ # input training sequence file
-h ${RBP}.training_peak.sprofile.txt  \ # input training sprofile file
-d ${RBP}.test_peak.sequences.txt \ # input test sequence file
-n ${RBP}.test_peak.sprofile.txt \ # input test sprofile file
-o ${RBP} # output directory

# RNApromo

cd /Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/5.RNApromo/AUB

# rename input fasta header
sed 's/>/>train_/g' /Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/2.peak_training/RNApromo/AUB.training_peak.fa > AUB.training_peak.fa

sed 's/>/>test_/g' /Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/3.peak_test/RNApromo/AUB.test_peak.fa > AUB.test_peak.fa

# merge input fasta
cat AUB.training_peak.fa AUB.test_peak.fa > AUB.fa
rm AUB.training_peak.fa
rm AUB.test_peak.fa

# need to specify global parameters
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Share2/home/zhaoweihao/lib64/ # need libstdc++.so.5 in the directory
export SEGAL_LAB_DEBUG=1

/Share2/home/zhaoweihao/app/RNApromo/rnamotifs08_motif_finder.pl \
-positive_seq AUB.fa \ # input fasta file
-output_dir /Share2/home/zhaoweihao/project/POSTAR3/motif/Fly/update/5.RNApromo/AUB/output \ # output directory
-bg 0.1 \ # p value cutoff
-n 3 # learn 3 motif for given sequence

# 6. summarize plot for each motif

# MEME

# run grep_PWM.pl to get the Position-Weight Matrix

perl grep_PWM.pl AUB \ # input RBP name
CTGAGGGTTT GGTCGACCRW CCTGGCRG # select three most significant motif
```

```perl
# grep_PWM.pl (for MEME)

use strict;

my @id;

push(@id, $ARGV[1]);
push(@id, $ARGV[2]);
push(@id, $ARGV[3]);

my $i;

foreach my $id (@id){
        $i++;
        open OUTPUT,">fly_CLIPdb_MEME_$ARGV[0]_$i.txt" or die;
        print OUTPUT "ID fly_CLIPdb_MEME_$ARGV[0]_$i\n";
        print OUTPUT "BF Drosophila melanogaster\n";
        print OUTPUT "P0\tA\tC\tG\tU\n";
        my $motif_anno = `grep -A 2 "position-specific probability matrix" ../4.MEME/$ARGV[0]/output/meme.txt |grep -A 2 "Motif $id"| tail -n 1`;
        chomp($motif_anno);
        $motif_anno =~ /letter-probability matrix: alength= (.*) w= (.*) nsites= (.*) E= (.*)/;
        my $w = $2;
        my $temp = 2 + $w;
        for(my $j=1;$j<=$w;$j++){
                my $line = `grep -A $temp "position-specific probability matrix" ../4.MEME/$ARGV[0]/output/meme.txt |grep -A $temp "Motif $id"| tail -n $w | head -n $j | tail -n 1`;
                chomp($line);
                my @line = split /\s+/,$line;
                if($j<10){
                        $j = "0".$j;
                }
                if($j==$w){
                        print OUTPUT "$j\t$line[1]\t$line[2]\t$line[3]\t$line[4]\nXX\n//\n";
                }
                else{
                        print OUTPUT "$j\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
                }
        }
}
```

```bash
# run WebLogo to plot MEME PWM

for i in `ls | grep txt | cut -d '.' -f 1`
do
        echo $i
        /Share2/home/lulab/zhuyumin/share/zhuyumin/apps/weblogo-master/weblogo \
        -f ${i}.txt \ # input file
        -D transfac \ # input file format
        -o ${i}.png \ # output file name
        -F png \ # output format
        -A rna \ # sequence type
        --errorbars False \ # no error bar
        --yaxis 1.0 \ # height of y axis
        --size medium \ # the logo size
        --fineprint 'WebLogo 3.4' \ # print on the lower right area
        --resolution 600 \ # resolution
        -c 'classic' # color scheme
done

# Homer

perl grep_PWM.pl AUB \ # input RBP name
3-GACCATCA 5-ATGGTCTG 1-GTCTGAC # select three most significant motif
```

```perl
# grep_PWM.pl (for Homer)

use strict;

my @id;

push(@id, $ARGV[1]);
push(@id, $ARGV[2]);
push(@id, $ARGV[3]);

my $i;

foreach my $id (@id){
        $i++;
        open OUTPUT,">fly_CLIPdb_HOMER_$ARGV[0]_$i.txt" or die;
        print OUTPUT "ID fly_CLIPdb_HOMER_$ARGV[0]_$i\n";
        print OUTPUT "BF Drosophila melanogaster\n";
        print OUTPUT "P0\tA\tC\tG\tU\n";
        $id =~ /(.*)-(.*)/;
        my $w = length($2);
        for(my $j=1;$j<=$w;$j++){
                my $line = `grep -A $w "$id" ../7.Homer_training/$ARGV[0]/output/homerMotifs.all.motifs | tail -n $w | head -n $j | tail -n 1`;
                chomp($line);
                my @line = split /\t/,$line;
                if($j<10){
                        $j = "0".$j;
                }
                if($j==$w){
                        print OUTPUT "$j\t$line[0]\t$line[1]\t$line[2]\t$line[3]\nXX\n//\n";
                }
                else{
                        print OUTPUT "$j\t$line[0]\t$line[1]\t$line[2]\t$line[3]\n";
                }
        }
}
```

```bash
# run WebLogo to plot Homer PWM (same as EME)

for i in `ls | grep txt | cut -d '.' -f 1`
do
        echo $i
        /Share2/home/lulab/zhuyumin/share/zhuyumin/apps/weblogo-master/weblogo -f ${i}.txt -D transfac -o ${i}.png -F png -A rna --errorbars False --yaxis 1.0 --size medium --fineprint 'WebLogo 3.4' --resolution 600 -c 'classic'
done

# RNAcontext

perl grep_PWM.pl AUB \ # input RBP name
8 9 10 # input three most significant motif length
```

```perl
# grep_PWM.pl (for RNAcontext)

use strict;

my @id;

push(@id, $ARGV[1]);
push(@id, $ARGV[2]);
push(@id, $ARGV[3]);

my $i;

foreach my $id (@id){
        $i++;
        open OUTPUT1,">fly_CLIPdb_RNAcontext_$ARGV[0]_$i.txt" or die;
        print OUTPUT1 "ID fly_CLIPdb_RNAcontext_$ARGV[0]_$i\n";
        print OUTPUT1 "BF Drosophila melanogaster\n";
        print OUTPUT1 "P0\tA\tC\tG\tU\n";
        my $line1 = `grep -A 2 "Motif width $id" ../6.RNAcontext/$ARGV[0]/outputs/params_$ARGV[0].txt | tail -n 1`;
        chomp($line1);
        my @line1 = split /\t/,$line1;
        my $line2 = `grep -A 3 "Motif width $id" ../6.RNAcontext/$ARGV[0]/outputs/params_$ARGV[0].txt | tail -n 1`;
        chomp($line2);
        my @line2 = split /\t/,$line2;
        my $line3 = `grep -A 4 "Motif width $id" ../6.RNAcontext/$ARGV[0]/outputs/params_$ARGV[0].txt | tail -n 1`;
        chomp($line3);
        my @line3 = split /\t/,$line3;
        my $line4 = `grep -A 5 "Motif width $id" ../6.RNAcontext/$ARGV[0]/outputs/params_$ARGV[0].txt | tail -n 1`;
        chomp($line4);
        my @line4 = split /\t/,$line4;
        my $w = $#line1;
        for(my $j=1;$j<=$w;$j++){
                if($j<10){
                        $j = "0".$j;
                }
                if($j==$w){
                        print OUTPUT1 "$j\t$line1[$j]\t$line2[$j]\t$line3[$j]\t$line4[$j]\nXX\n//\n";
                }
                else{
                        print OUTPUT1 "$j\t$line1[$j]\t$line2[$j]\t$line3[$j]\t$line4[$j]\n";
                }
        }
        open OUTPUT2,">fly_CLIPdb_RNAcontext_$ARGV[0]_$i\_structure.txt" or die;
        my $PLMU = `grep -A 11 "Motif width $id" ../6.RNAcontext/$ARGV[0]/outputs/params_$ARGV[0].txt | tail -n 4`;
        print OUTPUT2 "$PLMU";
}
```

```bash
# plot RNAcontext structure motif using R code

for i in `ls | grep structure.txt | cut -d '.' -f 1`
do
        Rscript plot.R ${i}.txt ${i}.png
done
```

```R
# plot.R
library(ggplot2)

Args <- commandArgs(T)
input <- Args[1]
output <- Args[2]

mx <- read.table(input, sep = '\t', header = F)
hah <- ggplot(mx, aes(x=V1,y=V2)) + 
  geom_bar(stat="identity",position=position_dodge(),fill="#999999",colour="black") + 
  ylab("Relative probability of binding") + 
  theme_bw() + 
  theme(axis.title.x=element_text(size=0),axis.title.y=element_text(size=20),axis.text.x=element_text(angle=45,hjust=1,size=22,family='sans'),axis.text.y=element_text(size=17,family='sans'))

png(output)
print(hah)
dev.off()
```

```bash
# RNApromo

# copy plot to summarize

cp ../5.RNApromo/AUB/output/model_cons_1.png ./mouse_CLIPdb_RNApromo_AUB_1_structure.png
```





### RBS

#### {species}_clipdb_exp

**Note**:

1. Contained peaks from collected CLIP-seq data, ENCODE eCLIP, and PIP-seq from PMID 24393486
2. When updating, split file into two parts: ENCODE eCLIP (replace with newest version) and others (merge with new peaks)
3. This table was further split into {species}\_clipdb\_{protocol_name} to ensure search efficiency even the table is large. Protocol naming convention is listed as follows:

| number | hits (HITS-CLIP)    | par (PAR-CLIP) | i (iCLIP) | a (iCLAP) | f (Fr-iCLIP) | su (4SU-iCLIP) | br (BrdU-CLIP) | ur (urea-iCLIP) | ce (eCLIP from experiment) | eclip (eCLIP from ENCODE) | pip (PIP-seq from 24393486) |
| ------ | ------------------- | -------------- | --------- | --------- | ------------ | -------------- | -------------- | --------------- | -------------------------- | ------------------------- | --------------------------- |
| 1      | CIMS                | PARalyzer      | CITS      | PureCLIP  | PureCLIP     | PureCLIP       | CTK            | PureCLIP        | PureCLIP                   | -                         | -                           |
| 2      | Piranha             | Piranha        | Piranha   | Piranha   | Piranha      | Piranha        | Piranha        | Piranha         | Piranha                    | -                         | -                           |
| 3      | CLIPper (human)     | MiClip         | PureCLIP  | -         | -            | -              | -              | -               | -                          | -                         | -                           |
| 4      | CTK (other species) | -              | -         | -         | -            | -              | -              | -               | -                          | -                         | -                           |



#### {species}\_clipdb\_pred

**Note**:

1. The pred table was built during POSTAR, with predicted binding sites from collected CLIP-seq data in POSTAR using FIMO, TESS, and DeepBind, only contains human and mouse
2. Could be updated in the next version (or discard?)



#### {species}_gene

**Note**:

1. Gene information for each RNA in POSTAR3
2. Search method: ensembl BioMart: Chromosome/scaffold name, Gene start, Gene end, Gene stable ID, Strand, Gene type, Gene name, transform into bed-like file (add 0 before Strand, change Strand into +/-)



#### CancerGene/CoreTF/Disease/DiseaseGene/Drug/SpecificGene

**Note**: 

1. Curated in previous versions, no need to update unless new data resource is available



#### {species}_geneExp

**Note**: 

1. Download from ExpressionAltas (no need to update unless new species is added)



#### {species}_RBP_hotspot

**Note**: 

1. Calculate using {species}_clipdb_exp on every 20nt bin for every transcript across genome



#### {species}_circRNAanno

**Note**:

1. circRNA annotation was obtained from circBase: list search -> input circRNAID -> export results (txt), then further formatted
2. Genome coordinate of circBase was hg19 for human, so liftOver to hg38



#### {species}_circRNAexp

**Note**: 

1. circRNA expression value was obtained from circBase (in the previous downloaded circRNA annotation file)
2. The expression value was total reads at circular junction in each study



### RNA Crosstalk

#### {species}_RBP_clipdb_exp_miRNApredict

**Note**: 

1. miRNA binding sites was predicted using miRanda, RNAhybrid, psRobot, psRNAtarget (script not found)
2. No need to update predicted miRNA binding sites. When updating this table, could extract miRNA information from existing table and overlap with new binding site records to generate new version of this table



#### {species}_RBP_clipdb_exp_miRNAvalidate

**Note**:

1. miRNA binding sites was calculated using AGO2 CLIP-seq data (script not found)
2. Like in previous part, no need to update experiment miRNA binding sites



#### {species}_RBP_clipdb_exp_RNAediting

**Note**: 

1. The editing site was curated from RADAR, DARNED (other species), PMID25373143 (worm), and GTEx (human, mouse)
2. Could be further updated if new RNA editing site database is available, otherwise update like in previous parts



#### {species}_RBP_clipdb_exp_RNAmod

**Note**: 

1. The modification site was curated from DMBase2
2. Could be further updated if new RNA modification site database is available, otherwise update like in previous parts



### Genomic Variants & Disease Mutations

**Note**:

1. Download data from corresponding database and overlap with binding site records to update

```bash
bedtools intersect \
-a ${other_data_bed} \ # bed file for other types of data (crosstalk, variation, and disease)
-b ${binding_site_bed} \ # bed file for RBP binding sites
-sorted \ # input bed file is sorted by "sort -k1,1 -k2,2n"
-wa \ # write -a file in output
-wb \ # write -b file in output
> ${output_file} # output file name
```



### Structurome

**Note**: 

1. Processed reactivity could be obtained from both analysis (following the original pipeline in the paper) and direct download from GEO processed files and process them into the uniform format



### Translatome

**Note**:

1. heatmap_den table was the raw density for each ORF
2. Bedgraph did not update (maybe discard in future versions?)



### Degradome

#### {species}_miRNAanno

**Note**:

1. The data was curated from miRBase (manually collected, did not find batch search and download)



#### {species}_miRNAGO

**Note**:

1. The data was curated from ensembl BioMart (only human has data)

