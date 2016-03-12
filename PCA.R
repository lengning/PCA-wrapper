
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # file name (output pca plot will show transformed data using PC calculated from this file)
n=as.numeric(options[2]) # number of PCs to output
LOD=as.numeric(options[3]) # lower limit of detection
Norm=options[4] # whether perform normalization; if "T" is specified, median-by-ratio normalization will be performed.
Plot=options[5] # whether plot
Projected = options[6]
FileSC = options[7] # file name (this file will be used for projection )

if(length(options)<2)n=5
if(length(options)<3)LOD=0
if(length(options)<4)Norm="F"
if(length(options)<5)Plot="T"
if(length(options)<6)Projected="F"

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

if(Projected=="T"){
  # csv or txt
  tmp=strsplit(FileSC, split="\\.")[[1]]
  FileType=tmp[length(tmp)]
  
  if(FileType=="csv"){
    cat("\n Read in csv file \n")
    prefix=strsplit(FileSC,split="\\.csv")[[1]][1]
    InSC=read.csv(FileSC,stringsAsFactors=F,row.names=1)
  }
  if(FileType!="csv"){
    cat("\n Read in tab delimited file \n")
    prefix=strsplit(FileSC,split=paste0("\\.",FileType))[[1]][1]
    InSC=read.table(FileSC,stringsAsFactors=F,row.names=1,sep="\t",header=T,quote="\"")
  }  
  MatSCraw=data.matrix(InSC)
  MaxSC=apply(MatSCraw,1,max)
  WhichSCRM=which(MaxSC<LOD)
  print(paste(length(WhichSCRM),"File2: genes with max expression < ", LOD, "are removed"))
  
  MatSC=MatSCraw
  if(length(WhichSCRM)>0)MatSC=MatSCraw[-WhichSCRM,]
  print(str(MatSC))
  
  if(Norm=="T"){
    cat("\n ==== Performing normalization ==== \n")
    library(EBSeq)
    Sizes=MedianNorm(MatSC)
    if(is.na(Sizes[1]))cat("\n Warning: all genes have 0(s), normalization is not performed \n")
    else MatSC=GetNormalizedMat(MatSC, MedianNorm(MatSC))
  }
  
  #Rescale
  MatSCscale=t(apply(MatSC,1,scale))
  rownames(MatSCscale)=rownames(MatSC)
  colnames(MatSCscale)=colnames(MatSC)
  MatSCscale[which(is.na(MatSCscale))]=0
}


Matraw=data.matrix(In)
Max=apply(Matraw,1,max)
WhichRM=which(Max<LOD)
print(paste(length(WhichRM),"File1: genes with max expression < ", LOD, "are removed"))


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
Matscale=t(apply(Mat,1,scale))
rownames(Matscale)=rownames(Mat)
colnames(Matscale)=colnames(Mat)
Matscale[which(is.na(Matscale))]=0

print(n)

PCAres=prcomp(t(Matscale))

if(Projected=="T"){
  if(length(all.equal(rownames(Matscale),rownames(MatSCscale)))>1) print("For PCA, use genes that are in both File1 (bulk) and File2 (sc)")
  Commg = intersect(rownames(Matscale),rownames(MatSCscale))
  Matscalecommg = Matscale[Commg,]
  PCAres=bkPCAres=prcomp(t(Matscalecommg))
  MatSCscalecommg = MatSCscale[Commg,]
  scPCAres=prcomp(t(MatSCscalecommg))  
  
  SConbulkPCA=t(PCAres$rotation) %*% MatSCscalecommg # bulk projected PCA
  PCAres = list(x=t(SConbulkPCA))
}


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

if(Projected!="T"){
PCA_sort=sapply(1:n,function(j){
  tmp=abs(PCAres$rotation[,j])
  t2=names(sort(tmp,decreasing=T))
})
colnames(PCA_sort)=paste0("PC",1:n)
}

if(Projected!="T")Perc=PCAres$sdev/sum(PCAres$sdev)
if(Projected=="T"){
  SCvar=apply(SConbulkPCA,1,var)
  Perc=SCvar/nrow(bkPCAres$rotation)
}
write.csv(Perc,file=paste0(prefix,"_perc_sdev.csv"))
if(Projected!="T") {
	write.csv(PCA_sort[,1:n],file=paste0(prefix,"_sort_by_absloading.csv"))
	write.csv(PCAres$rotation[,1:n],file=paste0(prefix,"_loading.csv"))
}
if(Plot=="T")Sys.sleep(1e30)


