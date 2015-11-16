
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # file name
n=as.numeric(options[2]) # number of PCs to output
LOD=as.numeric(options[3]) # lower limit of detection
Norm=options[4] # whether perform normalization; if "T" is specified, median-by-ratio normalization will be performed.
Plot=options[5] # whether plot

if(length(options)<2)n=5
if(length(options)<3)LOD=0
if(length(options)<4)Norm="F"
if(length(options)<5)Plot="T"

if(Plot=="T") X11()

# csv or txt
tmp=strsplit(File, split="\\.")[[1]]
FileType=tmp[length(tmp)]

if(FileType=="csv"){
	cat("\n Read in csv file \n")
	prefix=strsplit(File,split="\\.csv")[[1]][1]
	In=read.csv(File,stringsAsFactors=F,row.names=1)
}
if(FileType!="csv"){
	cat("\n Read in tab delimited file \n")
	prefix=strsplit(File,split=paste0("\\.",FileType))[[1]][1]
	In=read.table(File,stringsAsFactors=F,row.names=1,sep="\t",header=T,quote="\"")
}



Matraw=data.matrix(In)

Max=apply(Matraw,1,max)
WhichRM=which(Max<LOD)
print(paste(length(WhichRM),"genes with max expression < ", LOD, "are removed"))

Mat=Matraw
if(length(WhichRM)>0)Mat=Matraw[-WhichRM,]
print(str(Mat))

if(Norm=="T"){
cat("\n ==== Performing normalization ==== \n")
library(EBSeq)
Sizes=MedianNorm(Mat)
if(is.na(Sizes[1]))cat("\n Warning: all genes have 0(s), normalization is not performed \n")
else Mat=GetNormalizedMat(Mat, MedianNorm(Mat))
}

#Rescale
MatSC=t(apply(Mat,1,scale))
rownames(MatSC)=rownames(Mat)
colnames(MatSC)=colnames(Mat)
MatSC[which(is.na(MatSC))]=0

print(n)

PCAres=prcomp(t(MatSC))

if(Plot=="T"){
	pairs(PCAres$x[,1:n])
}
pdf(paste0(prefix,"_PC_pairs.pdf"),width=15,height=15)
pairs(PCAres$x[,1:n])
dev.off()

if(Plot=="T"){
library(rgl)
plot3d(PCAres$x[,1:3])
rgl.texts(PCAres$x[,1],PCAres$x[,2],PCAres$x[,3],colnames(Mat),col="darkgrey")
}
PCA_sort=sapply(1:n,function(j){
			 tmp=abs(PCAres$rotation[,j])
			 t2=names(sort(tmp,decreasing=T))
		})
colnames(PCA_sort)=paste0("PC",1:n)

Perc=PCAres$sdev/sum(PCAres$sdev)

write.csv(PCAres$rotation[,1:n],file=paste0(prefix,"_loading.csv"))
write.csv(Perc,file=paste0(prefix,"_perc_sdev.csv"))
write.csv(PCA_sort[,1:n],file=paste0(prefix,"_sort_by_absloading.csv"))

if(Plot=="T")Sys.sleep(1e30)


