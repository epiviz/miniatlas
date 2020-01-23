

# /usr/local/packages/r-3.6.0/bin/R


.libPaths('/local/projects-t3/idea/bherb/software/R_lib/R_3_6')


library(RcppParallel)
library(Seurat)
library(RJSONIO)
library(ShortRead)
library(Rsamtools)
library(readxl)
library(R.utils)


flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}


###########################
### Tracks for BICCN - MiniAtlas

### working off of Eran's integrated analysis 

intclust = read.csv('/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/Mini-Atlas_Publication_2019/Mukamel_IntegratedMOpAnalysis/miniatlas_joint_clusterings.tsv',sep='\t')



###########################
### ATAC data 


brain/biccn/grant/huang/macosko/transcriptome/sncell/processed/align/pBICCNsMMrMOpRAiF003d190318.bam.tar



BAMind <- dataSum$bamInd[tmpind]
        tmpfiles <- untar(totBAM[BAMind],list=TRUE) ### good - but even if I find a bai file, is it useful? 
        
        if(length(tmpfiles)==1){
            BAMfile <- totBAM[BAMind]
### skip if file size too large?
            if(file.info(BAMfile)$size>50000000000) next
    


BAMfile= '/local/projects-t3/NEMO/dmz/brain/biccn/grant/huang/macosko/transcriptome/sncell/processed/align/pBICCNsMMrMOpRMiF007d190314.bam.tar'



tmpdir = '/local/scratch/bherb/bamtest'

untar(BAMfile,exdir=tmpdir)

### for methylation             
#bam <- scanBam(paste(tmpdir,tmpfiles[2],sep="/"))

#tmpbam <- scanBam('/local/scratch/bherb/bamtest/pBICCNsMMrMOpRMiF007d190314/pBICCNsMMrMOpRMiF007d190314.bam') - too big to load? 

### smaller file - pBICCNsMMrPAGi2d180603


CEMBA171206_3C

### in command line - copy, then open in /local/scratch/bherb/bamtest


/local/projects-t3/NEMO/dmz/brain/biccn/grant/cemba/ecker/chromatin/scell/processed/align/3C/CEMBA171206_3C/3C1.nsrt.bw

track type=bigWig name="Check Big Wig" description="Samples?" bigDataUrl=http://data.nemoarchive.org/nemoHub/mm10/3C1.nsrt.bw

tmpbam2 <- scanBam('/local/projects-t3/NEMO/dmz/brain/biccn/grant/cemba/ecker/chromatin/scell/processed/counts/3C/CEMBA171206_3C/CEMBA171206_3C.sorted.bed.gz')



#tmpbed = read.table('/local/scratch/bherb/bamtest/CEMBA171206_3C.sorted.bed')


122338974 


### test with CEMBA171206_3C

int3C = intclust[grep('CEMBA171206_3C',intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot)) ## try 2 and 3

for(i in cells){
tmpnames = data.frame(names=flexsplit(int3C$sample[which(int3C$joint_cluster_round4_annot==i)],'_')[,5])
subname=gsub(' ','_',i)
subname=gsub('/','-',subname) ## rerun 
write.table(tmpnames,file=paste('/local/scratch/bherb/bamtest/CEMBA171206_3C_',subname,'.txt',sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

}

### ephatic 
#for i in ATTCAGAATACGCTGCCTCTCTATCCTATCCT CGCTCATTATGCGCAGTTATGCGACAGGACGT GAGATTCCCTCTGGTATTCTAGCTTATAGCCT; do 
#grep -w "$i" CEMBA171206_3C.sorted.bed
#done > test3.bed

sink(file='/local/scratch/bherb/bamtest/parse_CEMBA171206_3C.sh')
cat('#/bin/bash/')
cat('\n')
for(i in cells){
subname=gsub(' ','_',i)
cat(paste('for i in $(cat /local/scratch/bherb/bamtest/CEMBA171206_3C_',subname,'.txt); do grep -w "$i" /local/scratch/bherb/bamtest/CEMBA171206_3C.sorted.bed; done > /local/scratch/bherb/bamtest/CEMBA171206_3C_',subname,'.bed','\n',sep=""))
}
sink()

### and then just submit - qsub -P "sament-lab" -l mem_free=100G parse_CEMBA171206_3C.sh 

### looping 



samples = unique(apply(flexsplit(intclust$sample[grep('snatac_gene_CEMBA',intclust$sample)],'_')[,3:4],1,paste,collapse='_'))

for(k in samples){
dir.create(paste('/local/scratch/bherb/bamtest/',k,sep=""))
int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot))

for(i in cells){
tmpnames = data.frame(names=flexsplit(int3C$sample[which(int3C$joint_cluster_round4_annot==i)],'_')[,5])
subname=gsub(' ','_',i,fixed=TRUE)
subname=gsub('/','-',subname,fixed=TRUE)
write.table(tmpnames,file=paste('/local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.txt',sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}
}

### ephatic 
#for i in ATTCAGAATACGCTGCCTCTCTATCCTATCCT CGCTCATTATGCGCAGTTATGCGACAGGACGT GAGATTCCCTCTGGTATTCTAGCTTATAGCCT; do 
#grep -w "$i" CEMBA171206_3C.sorted.bed
#done > test3.bed

for(k in samples){
int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot))

sink(file=paste('/local/scratch/bherb/bamtest/',k,'/parse_',k,'.sh',sep=""))
cat('#/bin/bash/')
cat('\n')
for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
cat(paste('for i in $(cat /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.txt); do grep -w "$i" /local/scratch/bherb/bamtest/CEMBA171206_3C.sorted.bed; done > /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bed','\n',sep=""))
}
sink()
}



sink(file='/local/scratch/bherb/bamtest/qsub_calls.sh')
for(k in samples){
cat('\n\n')
cat(paste('qsub -P "sament-lab" -l mem_free=100G -o /local/scratch/bherb/bamtest/',k,'/',k,'_output.txt -e /local/scratch/bherb/bamtest/',k,'/',k,'_error.txt /local/scratch/bherb/bamtest/',k,'/parse_',k,'.sh','\n\n',sep=""))
}
sink()


### sort beds 
### testing with
/local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bed

### sort - works 
sort -k 1,1 /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bed > /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bed.sorted

### form bedgraph
genomeCoverageBed -bg -i /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bed.sorted -g /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes > /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bedgraph


### looping - CEMBA171207_3C will be test



for(k in samples){
int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot))

sink(file=paste('/local/scratch/bherb/bamtest/',k,'/sortBed_',k,'.sh',sep=""))
cat('#/bin/bash/')
cat('\n')
for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
cat(paste('sort -k 1,1 /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bed > /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bed.sorted','\n',sep=""))
}
sink()

sink(file=paste('/local/scratch/bherb/bamtest/',k,'/makeBedGraph_',k,'.sh',sep=""))
cat('#/bin/bash/')
cat('\n')
for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
cat(paste('genomeCoverageBed -bg -i /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bed.sorted -g /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes > /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bedgraph','\n',sep=""))
}
sink()

}


#qsub -P "sament-lab" -l mem_free=100G /local/scratch/bherb/bamtest/CEMBA171207_3C/sortBed_CEMBA171207_3C.sh

#qsub -P "sament-lab" -l mem_free=100G /local/scratch/bherb/bamtest/CEMBA171207_3C/makeBedGraph_CEMBA171207_3C.sh


### command line

### use UCSC chr sizes 

fetchChromSizes mm10 > hg18.chrom.sizes - should be ok, probably downloaded from UCSC before 




bedGraphToBigWig file.bg hg19_chrom_info.txt file.bw

bedGraphToBigWig /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bedgraph /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bw

### check

cp /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_ET_1.bw /local/projects-t3/NEMO/dmz/nemoHub/mm10

 http://data.nemoarchive.org/nemoHub/mm10/CEMBA171207_3C_L5_ET_1.bw

### good.

### do these need headers? 


### looping





for(k in samples){
int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot))

sink(file=paste('/local/scratch/bherb/bamtest/',k,'/makeBigWig_',k,'.sh',sep=""))
cat('#/bin/bash/')
cat('\n')
for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
cat(paste('bedGraphToBigWig /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bedgraph /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/scratch/bherb/bamtest/',k,'/',k,'_',subname,'.bw','\n',sep=""))
}
sink()

}


qsub -P "sament-lab" -l mem_free=100G /local/scratch/bherb/bamtest/CEMBA171207_3C/makeBigWig_CEMBA171207_3C.sh


#### copy files 

cp /local/scratch/bherb/bamtest/CEMBA171207_3C/*.bw /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC


### make trackDb.txt file 

for(k in samples){


bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/',k,'_ATAC/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

### maybe add a catch that elimiates files that do not end in .bw ... for now it works 

### also need to eliminate empty bw... 

int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round4_annot))

sink(file=paste('/local/scratch/bherb/bamtest/',k,'/trackDb.txt',sep=""))
cat(paste('track ',k,'\n','container multiWig','\n','shortLabel ',k,'\n','longLabel ',k,'\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(subname,bwFiles)]
if(length(tmpfile)==0) next
intInd = intersect(grep(k,intclust$sample),grep(i,intclust$joint_cluster_round4_annot))
hexcolor = unique(intclust$joint_cluster_round4_color[intInd])
rgbcolor = col2rgb(hexcolor[1])
rgbcol = paste(rgbcolor[,1],collapse=', ')
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'\n','longLabel ',k,'-',subname,'\n','parent ',k,'\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()
} #samples 



### command line

cp /local/scratch/bherb/bamtest/CEMBA171207_3C/trackDb.txt /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/mm10

#http://data.nemoarchive.org/nemoHub/CEMBA171207_3C_ATAC/hub.txt

#http://data.nemoarchive.org/nemoHub/CEMBA171207_3C_ATAC/mm10/

#http://data.nemoarchive.org/nemoHub/testHub/hub.txt 



### try to combine bigwig files to use joint_cluster_round2_annot clustering  - need to filter out empty bigwigs...




#for(k in samples){
k='CEMBA171207_3C'

bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/',k,'_ATAC/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

int3C = intclust[grep(k,intclust$sample),]
cells4 = unique(as.character(int3C$joint_cluster_round4_annot))
cells2 = unique(as.character(int3C$joint_cluster_round2_annot))

for(i in cells2){
round4cells = unique(int3C$joint_cluster_round4_annot[which(int3C$joint_cluster_round2_annot==i)])
### eliminate ones without bw
badind=NA
for(n in round4cells){
subname=gsub(' ','_',n)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(subname,bwFiles)]
if(length(tmpfile)==0){
badind = c(badind,which(round4cells==n))
}
}
badind=na.omit(badind)
if(length(badind)>0) round4cells = round4cells[-badind]

if(length(round4cells)==0){
next
} else if(length(round4cells)==1){
subname=gsub(' ','_',round4cells[1])
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(subname,bwFiles)]
system(command=paste('cp ',tmpfile,' /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C',sep=""))
} else {

for(j in round4cells){
subname=gsub(' ','_',as.character(j))
subname=gsub('/','-',subname,fixed=TRUE)
if(j==round4cells[1]){
tmpfile=bwFiles[grep(subname,bwFiles)]
} else {
tmpfile=c(tmpfile,bwFiles[grep(subname,bwFiles)])
}
}
tmpfile2 = paste(tmpfile,collapse=" ")
subname2=gsub(' ','_',i)
subname2=gsub('/','-',subname2,fixed=TRUE)

system(command=paste('bigWigMerge ',tmpfile2,' /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_',subname2,'.bedgraph',sep=""))

system(command=paste('bedGraphToBigWig /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_',subname2,'.bedgraph /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_',subname2,'.bw',sep=""))
}
}




#### copy files  - command line 

mkdir /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2

mkdir /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2/mm10

cp /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/*.bw /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2/mm10

cp /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2

cp /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/hub.txt /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2 #edit hub.txt 


### make trackDb.txt file 

#for(k in samples){

bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/',k,'_ATAC_round2/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

### maybe add a catch that elimiates files that do not end in .bw ... for now it works 

### also need to eliminate empty bw... 

int3C = intclust[grep(k,intclust$sample),]
cells = unique(as.character(int3C$joint_cluster_round2_annot))

sink(file=paste('/local/scratch/bherb/bamtest/round2cluster/',k,'/trackDb.txt',sep=""))
cat(paste('track ',k,'_round2','\n','container multiWig','\n','shortLabel ',k,'_round2','\n','longLabel ',k,'_round2','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(subname,bwFiles)]
if(length(tmpfile)==0) next
intInd = intersect(grep(k,intclust$sample),grep(i,intclust$joint_cluster_round2_annot))
hexcolor = unique(intclust$joint_cluster_round4_color[intInd])
rgbcolor = col2rgb(hexcolor[1])
rgbcol = paste(rgbcolor[,1],collapse=', ')
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_round2','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'_round2','\n','longLabel ',k,'-',subname,'_round2','\n','parent ',k,'_round2','\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()


} #samples 


### command line 
cp /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/trackDb.txt /local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC_round2/mm10

http://data.nemoarchive.org/nemoHub/CEMBA171207_3C_ATAC_round2/hub.txt



### testing 

system(command='bigWigMerge /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_IT_Rspo1_1.bw /local/scratch/bherb/bamtest/CEMBA171207_3C/CEMBA171207_3C_L5_IT_Rspo2_1.bw /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_L5_IT_Rspo1_1.bedgraph')

system(command='bedGraphToBigWig /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_L5_IT_Rspo1_1.bedgraph /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/scratch/bherb/bamtest/round2cluster/CEMBA171207_3C/CEMBA171207_3C_L5_IT_Rspo1_1.bw')

table(intclust$joint_cluster_round2_annot[grep('L5 ET_',intclust$joint_cluster_round4_annot)])



### add colors above? have colors? 
## -trackline -trackopts name='name="CEMBA171207_3C_L5_ET_1" visibility=2 color=#3CBC78'


Convert bed files to bedgraph, then bigwig - might need to sort bed first

example:

bedtools genomecov -ibam -bg -split -strand + -i accepted_hits.bam -g ChromInfo.txt > accepted_hits.plus.bedGraph

bedGraphToBigWig accepted_hits.plus.bedGraph ChromInfo.txt accepted_hits.plus.bw

need ChromInfo....



library(XML)
masterColor = xmlParse('/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/EpigenomeClustering/Total.hue_by_name.xml') 

masterColor = xmlToList(masterColor)

masterColor = unlist(masterColor[[2]])

masterColor = data.frame(cellType = masterColor[which(names(masterColor)=='Track.name')], color=masterColor[which(names(masterColor)=='Track.color')])

MajCellTypes = flexsplit(masterColor$cellType,' - ')[,2]

cellind = match(unique(MajCellTypes),MajCellTypes)

masterColor2 = data.frame(cellType = MajCellTypes[cellind], color = masterColor$color[cellind])

rgb2hex <- function(a){
hexcol=NA
for(i in c(1:length(a))){
tmp = as.numeric(flexsplit(a[i],",")[1,])
hexcol[i] = rgb(tmp[1], tmp[2], tmp[3], maxColorValue = 255)
}
return(hexcol)
}

masterColor2$hex = rgb2hex(masterColor2$color)
masterColor2$subname = gsub(" ","_",masterColor2$cellType,fixed=TRUE)

### can be further reduced.... common color per cell type  


###############
#### DNA  methylation files  - update - new clustering from Ecker group
 ## old file  = /local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster_depreciated

## CG methylation 

#/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/Mini-Atlas_Publication_2019/EckerRen_Mouse_MOp_methylation+ATAC/bigwig/MajorCluster


dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster')
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster')
dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/mm10')
file.copy('/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/EpigenomeClustering/mCG/*.rate.bw','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/mm10') ## didn't work - just used cp 



### make trackDb.txt file 

#for(k in samples){

bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

### maybe add a catch that elimiates files that do not end in .bw ... for now it works 

### also need to eliminate empty bw... 

#int3C = intclust[grep('snmcseq',intclust$sample),]
#cells = unique(as.character(int3C$joint_cluster_round2_annot))

cells=list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/mm10',sep=""),pattern='.bw',full.names=FALSE)
cells = flexsplit(cells,'.')[,1]

colorind = match(cells,masterColor2$subname)
hexcolors = masterColor2$hex[colorind]
rgbcolors = as.character(masterColor2$color[colorind])

## get color from master color 

## here 
#colors=sample(unique(intclust$joint_cluster_round3_color),length(cells)) #random, since nothing matches

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/mm10/trackDb.txt',sep=""))
cat(paste('track Ecker_CG_Methylation_MajorCluster','\n','container multiWig','\n','shortLabel Ecker_CG_Methylation_MajorCluster','\n','longLabel Ecker_CG_Methylation_MajorCluster','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(i,bwFiles)]
if(length(tmpfile)==0) next
#intInd = intersect(grep(k,intclust$sample),grep(i,intclust$joint_cluster_round2_annot))
#hexcolor = unique(intclust$joint_cluster_round4_color[intInd])
#hexcolor = colors[which(cells==i)]
#rgbcolor = col2rgb(hexcolor[1])
#rgbcol = paste(rgbcolor[,1],collapse=', ')
rgbcol = rgbcolors[which(cells==i)]
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'_MajorCluster','\n','longLabel Ecker_CG_Methylation_MajorCluster-',subname,'_MajorCluster','\n','parent Ecker_CG_Methylation_MajorCluster','\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()

## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CG_Methylation_MajorCluster/hub.txt',sep=""))
cat(paste('hub Ecker_CG_Methylation_MajorCluster','\n','shortLabel Ecker_CG_Methylation_MajorCluster','\n','longLabel Ecker_CG_Methylation_MajorCluster','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()

http://data.nemoarchive.org/nemoHub/Ecker_CG_Methylation_MajorCluster/hub.txt



################
### CG Hyper peak beds  (moved files to mm10)  ## here   - need to convert beds to big beds - wait for Jayaram
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/hyperDMR')

bbFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/hyperDMR/mm10',sep=""),pattern='.bb',full.names=TRUE)

sizes = file.info(bbFiles)$size

bbFiles = bbFiles[which(sizes>400)]

cells=list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/hyperDMR/mm10',sep=""),pattern='.bb',full.names=FALSE)
cells = flexsplit(cells,'.')[,1]

colorind = match(cells,masterColor2$subname)
hexcolors = masterColor2$hex[colorind]
rgbcolors = as.character(masterColor2$color[colorind])

#colors=sample(unique(intclust$joint_cluster_round3_color),length(cells)) #random, since nothing matches - and using same colors from CG - cells in same order

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/hyperDMR/mm10/trackDb.txt',sep=""))
#cat(paste('track Peak','\n','container multiWig','\n','shortLabel Peak','\n','longLabel Peak','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bbFiles[grep(i,bbFiles)]
if(length(tmpfile)==0) next
rgbcol = rgbcolors[which(cells==i)]
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster_Peaks','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'\n','longLabel ',subname, '_ATAC','\n','visibility dense','\n','type bigBed','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()


## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Peak/hub.txt',sep=""))
cat(paste('hub Ecker_ATAC_MajorCluster_Peaks','\n','shortLabel Ecker_ATAC_MajorCluster_Peaks','\n','longLabel Ecker_ATAC_MajorCluster_Peaks','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()

http://data.nemoarchive.org/nemoHub/Peak/hub.txt







#############
### CH methylation 

dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster')
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster')
dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/mm10')
file.copy('/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/EpigenomeClustering/mCH/*.bw','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/mm10') ## didn't work - just used cp 



### make trackDb.txt file 

#for(k in samples){

bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

### maybe add a catch that elimiates files that do not end in .bw ... for now it works 

### also need to eliminate empty bw... 

#int3C = intclust[grep('snmcseq',intclust$sample),]
#cells = unique(as.character(int3C$joint_cluster_round2_annot))

cells=list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/mm10',sep=""),pattern='.bw',full.names=FALSE)
cells = flexsplit(cells,'.')[,1]

colorind = match(cells,masterColor2$subname)
hexcolors = masterColor2$hex[colorind]
rgbcolors = as.character(masterColor2$color[colorind])

#colors=sample(unique(intclust$joint_cluster_round3_color),length(cells)) #random, since nothing matches - and using same colors from CG - cells in same order

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/mm10/trackDb.txt',sep=""))
cat(paste('track Ecker_CH_Methylation_MajorCluster','\n','container multiWig','\n','shortLabel Ecker_CH_Methylation_MajorCluster','\n','longLabel Ecker_CH_Methylation_MajorCluster','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(i,bwFiles)]
if(length(tmpfile)==0) next
#intInd = intersect(grep(k,intclust$sample),grep(i,intclust$joint_cluster_round2_annot))
#hexcolor = unique(intclust$joint_cluster_round4_color[intInd])
#hexcolor = colors[which(cells==i)]
#rgbcolor = col2rgb(hexcolor[1])
#rgbcol = paste(rgbcolor[,1],collapse=', ')
rgbcol = rgbcolors[which(cells==i)]
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'_MajorCluster','\n','longLabel Ecker_CH_Methylation_MajorCluster-',subname,'_MajorCluster','\n','parent Ecker_CH_Methylation_MajorCluster','\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()

## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_CH_Methylation_MajorCluster/hub.txt',sep=""))
cat(paste('hub Ecker_CH_Methylation_MajorCluster','\n','shortLabel Ecker_CH_Methylation_MajorCluster','\n','longLabel Ecker_CH_Methylation_MajorCluster','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()

http://data.nemoarchive.org/nemoHub/Ecker_CH_Methylation_MajorCluster/hub.txt



#############
### Ecker ATAC

dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster')
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster')
dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/mm10')
file.copy('/local/projects-t3/NEMO/incoming/biccn/miniatlas/BICCN_MiniBrain_JointDataAnalysis/EpigenomeClustering/ATAC/*.bw','/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/mm10') ## didn't work - just used cp 

bwFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/mm10',sep=""),pattern='.bw',full.names=TRUE)

sizes = file.info(bwFiles)$size

bwFiles = bwFiles[which(sizes>400)]

cells=list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/mm10',sep=""),pattern='.bw',full.names=FALSE)
cells = flexsplit(cells,'.')[,1]

colorind = match(cells,masterColor2$subname)
hexcolors = masterColor2$hex[colorind]
rgbcolors = as.character(masterColor2$color[colorind])

#colors=sample(unique(intclust$joint_cluster_round3_color),length(cells)) #random, since nothing matches - and using same colors from CG - cells in same order

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/mm10/trackDb.txt',sep=""))
cat(paste('track Ecker_ATAC_MajorCluster','\n','container multiWig','\n','shortLabel Ecker_ATAC_MajorCluster','\n','longLabel Ecker_ATAC_MajorCluster','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(i,bwFiles)]
if(length(tmpfile)==0) next
rgbcol = rgbcolors[which(cells==i)]
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'_MajorCluster','\n','longLabel Ecker_ATAC_MajorCluster-',subname,'_MajorCluster','\n','parent Ecker_ATAC_MajorCluster','\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()

## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Ecker_ATAC_MajorCluster/hub.txt',sep=""))
cat(paste('hub Ecker_ATAC_MajorCluster','\n','shortLabel Ecker_ATAC_MajorCluster','\n','longLabel Ecker_ATAC_MajorCluster','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()

http://data.nemoarchive.org/nemoHub/Ecker_ATAC_MajorCluster/hub.txt




################
### ATAC peak beds  (moved files to mm10)
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Peak')

bbFiles = list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Peak/mm10',sep=""),pattern='.bb',full.names=TRUE)

sizes = file.info(bbFiles)$size

bbFiles = bbFiles[which(sizes>400)]

cells=list.files(path=paste('/local/projects-t3/NEMO/dmz/nemoHub/Peak/mm10',sep=""),pattern='.bb',full.names=FALSE)
cells = flexsplit(cells,'.')[,1]

colorind = match(cells,masterColor2$subname)
hexcolors = masterColor2$hex[colorind]
rgbcolors = as.character(masterColor2$color[colorind])

#colors=sample(unique(intclust$joint_cluster_round3_color),length(cells)) #random, since nothing matches - and using same colors from CG - cells in same order

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Peak/mm10/trackDb.txt',sep=""))
#cat(paste('track Peak','\n','container multiWig','\n','shortLabel Peak','\n','longLabel Peak','\n','type bigWig 0 30000','\n','viewLimits 0:100','\n','visibility full','\n','maxHeightPixels 150:30:11','\n','aggregate transparentOverlay','\n','showSubtrackColorOnUi on','\n','windowingFunction mean','\n','priority 1.4','\n','configurable on','\n','autoScale on','\n\n',sep=""))

for(i in cells){
subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bbFiles[grep(i,bbFiles)]
if(length(tmpfile)==0) next
rgbcol = rgbcolors[which(cells==i)]
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster_Peaks','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'\n','longLabel ',subname, '_ATAC','\n','visibility dense','\n','type bigBed','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()


## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Peak/hub.txt',sep=""))
cat(paste('hub Ecker_ATAC_MajorCluster_Peaks','\n','shortLabel Ecker_ATAC_MajorCluster_Peaks','\n','longLabel Ecker_ATAC_MajorCluster_Peaks','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()

http://data.nemoarchive.org/nemoHub/Peak/hub.txt





##############
### smartseq cells 

mouseGenes = read.csv('/local/projects-t3/idea/bherb/NeMO_work/MouseGene_Pos.csv',row.names=1)

tot = read.csv('https://obj.umiacs.umd.edu/nemo-miniatlas/by_cell_type/smater_cells_collapsed_exon_counts.csv')

data = tot[,-c(1:2)]

geneLoc = data.frame(gene=tot$gene,chr="NONE",start=NA,end=NA,stringsAsFactors=FALSE)

for(i in c(1:nrow(geneLoc))){
tmpind = which(mouseGenes$mgi_symbol==as.character(geneLoc$gene[i]))
if(length(tmpind)>0){
geneLoc$chr[i]=as.character(mouseGenes$chromosome_name[tmpind[1]])
geneLoc$start[i]=mouseGenes$start_position[tmpind[1]]
geneLoc$end[i]=mouseGenes$end_position[tmpind[1]]
}
if(i%%100==0) cat(paste(i,', ',sep=""))
}

### 12603 not present

naind = which(geneLoc$chr=='NONE')

tot = tot[-naind,]
data = data[-naind,]
geneLoc = geneLoc[-naind,]

geneLoc$ensembl = mouseGenes$ensembl_gene_id[match(geneLoc$gene,mouseGenes$mgi_symbol)]

geneLoc$chr=paste('chr',as.character(geneLoc$chr),sep="")

badChr=which(geneLoc$chr=='chrCHR_MG4264_PATCH'|geneLoc$chr=='chrJH584293.1') #

tot = tot[-badChr,]
data = data[-badChr,]
geneLoc = geneLoc[-badChr,]

#### genes out of bounds 

chrSizes = read.delim('/local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes',header=FALSE)

badSize=NA
for(i in unique(geneLoc$chr)){
tmpind=which(geneLoc$chr==i)
tmpbad=which(geneLoc$end[tmpind]>chrSizes[which(chrSizes[,1]==i),2])

} ## ok no


### round data

data=round(data,2)

### order cell names 

data = data[,order(colnames(data))]

cells = gsub('.',' ',colnames(data),fixed=TRUE)
subname=gsub(' ','_',cells)
subname=gsub('/','-',subname,fixed=TRUE)

### looking to create a BED6+ bedfile - example:

colors = intclust$joint_cluster_round3_color[match(cells,intclust$joint_cluster_round3_annot)]

colors[1] = '#6F836B'  ## had to do manually, not perfect match   

rgbColors = col2rgb(colors)


dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster')
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster')
dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10')

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.bed',sep=""))
#cat(paste('chrom','\t','chromStart','\t','chromEnd','\t','name','\t','score','\t','strand','\t','name2','\t','expCount','\t','expScores','\t','_dataOffset','\t','_dataLen','\n',sep=""))

for(i in c(1:nrow(geneLoc))){
cat(paste(geneLoc$chr[i],'\t',geneLoc$start[i],'\t ',geneLoc$end[i],'\t ',geneLoc$gene[i],'\t','999','\t','-','\t',geneLoc$ensembl[i],'\t',length(subname),'\t',paste(data[i,],collapse=','),'\t','0','\t','0','\n',sep=""))
}
sink()

system(command='bedSort /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.bed /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster_sorted.bed')

file.copy('/local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/Test_10XnucleiMO.as','/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.as')

system(command='bedToBigBed -as=/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.as -type=bed6+5 /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster_sorted.bed /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.bb')

## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/hub.txt',sep=""))
cat(paste('hub Zeng_SMARTer_cells_MajorCluster','\n','shortLabel Zeng_SMARTer_cells_MajorCluster','\n','longLabel Zeng_SMARTer_cells_MajorCluster','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()


bbFile = list.files('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10',pattern='.bb',full.names=TRUE)

bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',bbFile)

## trackDb.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/trackDb.txt',sep=""))
cat(paste('track Zeng_SMARTer_cells_MajorCluster','\n','type bigBarChart','\n','visibility full','\n','shortLabel Gene Expression','\n','longLabel Zeng_SMARTer_cells_MajorCluster','\n','barChartBars ',paste(subname,collapse=' '),'\n','barChartColors ',paste(colors,collapse=' '),'\n','bigDataUrl ',bigurl,'\n\n',sep=""))
sink()


http://data.nemoarchive.org/nemoHub/Zeng_SMARTer_cells_MajorCluster/hub.txt








##############
### smartseq nuclei

mouseGenes = read.csv('/local/projects-t3/idea/bherb/NeMO_work/MouseGene_Pos.csv',row.names=1)

tot = read.csv('https://obj.umiacs.umd.edu/nemo-miniatlas/by_cell_type/smater_nuclei_collapsed_exon_counts.csv')

data = tot[,-c(1:2)]

geneLoc = data.frame(gene=tot$gene,chr="NONE",start=NA,end=NA,stringsAsFactors=FALSE)

for(i in c(1:nrow(geneLoc))){
tmpind = which(mouseGenes$mgi_symbol==as.character(geneLoc$gene[i]))
if(length(tmpind)>0){
geneLoc$chr[i]=as.character(mouseGenes$chromosome_name[tmpind[1]])
geneLoc$start[i]=mouseGenes$start_position[tmpind[1]]
geneLoc$end[i]=mouseGenes$end_position[tmpind[1]]
}
if(i%%100==0) cat(paste(i,', ',sep=""))
}

### 12603 not present

naind = which(geneLoc$chr=='NONE')

tot = tot[-naind,]
data = data[-naind,]
geneLoc = geneLoc[-naind,]

geneLoc$ensembl = mouseGenes$ensembl_gene_id[match(geneLoc$gene,mouseGenes$mgi_symbol)]

geneLoc$chr=paste('chr',as.character(geneLoc$chr),sep="")

badChr=which(geneLoc$chr=='chrCHR_MG4264_PATCH'|geneLoc$chr=='chrJH584293.1') ##

tot = tot[-badChr,]
data = data[-badChr,]
geneLoc = geneLoc[-badChr,]

#### genes out of bounds 

chrSizes = read.delim('/local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes',header=FALSE)

badSize=NA
for(i in unique(geneLoc$chr)){
tmpind=which(geneLoc$chr==i)
tmpbad=which(geneLoc$end[tmpind]>chrSizes[which(chrSizes[,1]==i),2])
} ## ok none...


### round data

data=round(data,2)

### order cell names 

data = data[,order(colnames(data))]

cells = gsub('.',' ',colnames(data),fixed=TRUE)
subname=gsub(' ','_',cells)
subname=gsub('/','-',subname,fixed=TRUE)

### looking to create a BED6+ bedfile - example:

colors = intclust$joint_cluster_round3_color[match(cells,intclust$joint_cluster_round3_annot)]

colors[1] = '#6F836B'  ## had to do manually, not perfect match   

rgbColors = col2rgb(colors)

dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster')
file.copy('/local/projects-t3/NEMO/dmz/nemoHub/CEMBA171207_3C_ATAC/genomes.txt','/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster')
dir.create('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10')

sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster.bed',sep=""))
#cat(paste('chrom','\t','chromStart','\t','chromEnd','\t','name','\t','score','\t','strand','\t','name2','\t','expCount','\t','expScores','\t','_dataOffset','\t','_dataLen','\n',sep=""))

for(i in c(1:nrow(geneLoc))){
cat(paste(geneLoc$chr[i],'\t',geneLoc$start[i],'\t ',geneLoc$end[i],'\t ',geneLoc$gene[i],'\t','999','\t','-','\t',geneLoc$ensembl[i],'\t',length(subname),'\t',paste(data[i,],collapse=','),'\t','0','\t','0','\n',sep=""))
}
sink()

system(command='bedSort /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster.bed /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster_sorted.bed')

file.copy('/local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/Test_10XnucleiMO.as','/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster.as')

system(command='bedToBigBed -as=/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster.as -type=bed6+5 /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster_sorted.bed /local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub/mm10.chrom.sizes /local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/Zeng_SMARTer_nuclei_MajorCluster.bb')

## hub.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/hub.txt',sep=""))
cat(paste('hub Zeng_SMARTer_nuclei_MajorCluster','\n','shortLabel Zeng_SMARTer_nuclei_MajorCluster','\n','longLabel Zeng_SMARTer_nuclei_MajorCluster','\n','genomesFile genomes.txt','\n','email  bherb@som.umaryland.edu','\n\n',sep=""))
sink()


bbFile = list.files('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10',pattern='.bb',full.names=TRUE)

bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',bbFile)

## trackDb.txt file 
sink(file=paste('/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/mm10/trackDb.txt',sep=""))
cat(paste('track Zeng_SMARTer_nuclei_MajorCluster','\n','type bigBarChart','\n','visibility full','\n','shortLabel Gene Expression','\n','longLabel Zeng_SMARTer_nuclei_MajorCluster','\n','barChartBars ',paste(subname,collapse=' '),'\n','barChartColors ',paste(colors,collapse=' '),'\n','bigDataUrl ',bigurl,'\n\n',sep=""))
sink()


http://data.nemoarchive.org/nemoHub/Zeng_SMARTer_nuclei_MajorCluster/hub.txt







### track info - no? just copied 

sink(file='/local/projects-t3/NEMO/dmz/nemoHub/Zeng_SMARTer_cells_MajorCluster/mm10/Zeng_SMARTer_cells_MajorCluster.as')
cat(paste('table bigBarChart','\n','"bigBarChart bar graph display"','\n','(
    string chrom;               "Reference sequence chromosome or scaffold"','\n','uint chromStart;            "Start position in chromosome"','\n','uint chromEnd;              "End position in chromosome"','\n','string name;                "Name or ID of item"','\n','uint score;                 "Score (0-1000)"','\n','char[1] strand;             "'+','-' or '.'. Indicates whether the query aligns to the + or - strand on the reference"','\n','string name2;               "Alternate name of item"','\n','uint expCount;              "Number of bar graphs in display, must be <= 100"','\n','float[expCount] expScores;  "Comma separated list of category values."','\n','bigint _dataOffset;         "Offset of sample data in data matrix file, for boxplot on details page, optional only for barChart format"','\n','int _dataLen;               "Length of sample data row in data matrix file, optional only for barChart format"','\n',')',sep=""))








bedSort Test_10XnucleiMO.bed Test_10XnucleiMO_sorted.bed ## had to delete header first 
## add back header to Test_10XnucleiMO_sorted.bed? - nope

bedToBigBed -as=Test_10XnucleiMO.as -type=bed6+5 Test_10XnucleiMO_sorted.bed mm10.chrom.sizes Test_10XnucleiMO.bb

bigBedInfo Test_10XnucleiMO.bb


BED6+5 with additional fields for category count and median values, and sample matrix fields

http://data.nemoarchive.org/nemoHub/hub.txt

/local/projects-t3/NEMO/dmz/nemoHub/mm10
/local/projects-t3/idea/bherb/NeMO_work/NemoTrackHub



















 95086227        95158010        DICER1  999     -       ENSG00000100697.10      5       2.94,11.60,38.00,6.69,4.89)

# chrom chromStart chromEnd string name score strand name2 expCount expScores _dataOffset _dataLen
chr14   95086227        95158010        DICER1  999     -       ENSG00000100697.10      5       2.94,11.60,38.00,6.69,4.89





subname=gsub(' ','_',i)
subname=gsub('/','-',subname,fixed=TRUE)
tmpfile=bwFiles[grep(i,bwFiles)]
if(length(tmpfile)==0) next
#intInd = intersect(grep(k,intclust$sample),grep(i,intclust$joint_cluster_round2_annot))
#hexcolor = unique(intclust$joint_cluster_round4_color[intInd])
hexcolor = colors[which(cells==i)]
rgbcolor = col2rgb(hexcolor[1])
rgbcol = paste(rgbcolor[,1],collapse=', ')
bigurl = gsub('/local/projects-t3/NEMO/dmz/nemoHub','http://data.nemoarchive.org/nemoHub',tmpfile)
cat(paste('track ',subname,'_MajorCluster','\n','bigDataUrl ',bigurl,'\n','shortLabel ',subname,'_MajorCluster','\n','longLabel Ecker_CH_Methylation_MajorCluster-',subname,'_MajorCluster','\n','parent Ecker_CH_Methylation_MajorCluster','\n','type bigWig','\n','color ',rgbcol,'\n\n',sep=""))
}
sink()






track type=barChart name="barChart Example One" description="A barChart file" barChartBars="adiposeSubcut breastMamTissue colonTransverse muscleSkeletal wholeBlood" visibility=pack
browser position chr14:95,081,796-95,436,280
# chrom chromStart chromEnd string name score strand name2 expCount expScores _dataOffset _dataLen
chr14   95086227        95158010        DICER1  999     -       ENSG00000100697.10      5       2.94,11.60,38.00,6.69,4.89
chr14   95181939        95319906        CLMN    999     -       ENSG00000165959.7       5       7.08,69.53,9.32,1.38,1.68
chr14   95417493        95475836        SYNE3   999     -       ENSG00000176438.8       5       7.29,3.73,0.74,20.35,1.39

### Junbque 

cat(paste('track type=barChart name="Zeng_SMARTer_cells_MajorCluster" description="A barChart file" barChartBars="', paste(subname,collapse=' '),'", visibility=pack
browser position chr14:95,081,796-95,436,280','\n\n',sep=""))




