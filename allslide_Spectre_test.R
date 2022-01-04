#Set1 slide 5개 합쳐서 clustering

#Asinh 0.001 (숫자 작아질수록 density 분리됨)
##MHCI이 0인 cell은 없앰
##area에서 50보다 작거나 920보다 큰 cell은 제거
setwd("F:/MICSSS_Analysis/20211229 spectre/")
library(dplyr)
library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages



### Import data

list.files('F:/MICSSS_Analysis/flowjo/export_Set1_flowjo/', '.csv')

data.list <- Spectre::read.files(file.loc = 'F:/MICSSS_Analysis/flowjo/export_Set1_flowjo/',
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)

###########
check <- do.list.summary(data.list)

check$name.table # Review column names and their subsequent values
check$ncol.check # Review number of columns (features, markers) in each sample
check$nrow.check # Review number of rows (cells) in each sample

data.list[[1]]

### Merge data

cell.dat <- Spectre::do.merge.files(dat = data.list)
cell.dat

cell.dat <- cell.dat[!(cell.dat$`MHC class I` <= 0), ]
cell.dat2 <- cell.dat[!(cell.dat$`Area (px^2)` < 50), ]
cell.dat <- cell.dat2[!(cell.dat2$`Area (px^2)` > 920), ]
summary(cell.dat)
##########################################################################################################
#### Import meta data  ###########################################################
################################################################################

#metadata는 flowjo에서 export한 파일 list 액셀로 붙여넣기해서 만듦       
meta.dat <- fread('F:/MICSSS_Analysis/20211229 spectre/metadata/sample.details.csv')
meta.dat
getwd()
fwrite(cell.dat, "set1_merge_data.csv")
##########################################################################################################
#### Asinh transformation
##########################################################################################################

### Co-factor targets

cf.low <- 0.001


### Transformation settings

as.matrix(names(cell.dat))
cell.dat <- select(cell.dat, c(1, 2:6, 8:23, 7)) #CellID가 제일앞에 오게 정렬
asinh.low <- names(cell.dat)[c(3:12)] #AEC빼고 단백질들 선택


### Test transformation settings on subsampled data

sub <- do.subsample(cell.dat, 100000)
sub <- do.asinh(sub, use.cols = asinh.low, cofactor = cf.low)


summary(sub)
as.matrix(names(sub))

transf.cols <- names(sub)[grepl('_asinh', names(sub))]
as.matrix(transf.cols)

### Make plots of transformed columns from the subsampled data
##(전체데이터사용하면 엄청 오래걸림!)

plot.against <- 'Pax5_asinh'
which(names(sub) == plot.against)

for(i in transf.cols){
  make.colour.plot(sub, i, plot.against)
} 

### Apply transformation to full dataset

cell.dat <- do.asinh(cell.dat, use.cols = asinh.low, cofactor = cf.low)


as.matrix(names(cell.dat))

transf.cols <- names(cell.dat)[grepl('_asinh', names(cell.dat))]
as.matrix(transf.cols)

##########################################################################################################
#### Add metadata
##########################################################################################################

meta.dat

sample.info <- meta.dat
sample.info

### Add sample metadata to primary data.table

cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)


##########################################################################################################
#### Write data to disk
##########################################################################################################

##false positive data 정리하기(CD3+PAX5+, CD8+PAX5+)
JN <- cell.dat[!(cell.dat$Pax5_asinh > 0 & cell.dat$CD3_asinh > 0), ]
JN_ <- JN[!(JN$Pax5_asinh > 0 & JN$CD8_asinh > 0), ]
Jsub <- do.subsample(JN_, 50000)
make.colour.plot(Jsub, 'CD3_asinh', 'Pax5_asinh')
make.colour.plot(Jsub, 'CD8_asinh', 'Pax5_asinh')

### Write cellular data and analysis  preferences to disk

fwrite(JN_, "cell.data.csv") # data

### Save session info to disk

setwd(OutputDirectory)
dir.create("Output - info", showWarnings = FALSE)
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()
######여기까지 데이터 전처리 이제 clustering


##########################################################################################################
#### Setup preferences
##########################################################################################################
### Sample preferences

as.matrix(names(JN_))

sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Downsample preferences
sub.by <- "Group"
as.matrix(unique(JN_[[sub.by]]))

data.frame(table(JN_[[group.col]]))

sub.targets <- c(10000, 10000, 10000, 10000, 10000) #slide group이 5개라서
sub.targets

### Clustering preferences

## Cellular cols
as.matrix(names(JN_))

cellular.cols <- names(JN_)[c(24:33)] #모든 asinh값 
cellular.cols

## Columns for clustering
clustering.cols <- names(JN_)[c(25,26,32,29,27)]#CD3, CD8, Pax5, MHCI, Gal-9
clustering.cols

## Cluster numbers etc
metak <- 'auto'

##########################################################################################################
#### Run clustering and dimensionality reduction
##########################################################################################################

### Run clustering
#column정리(안해도됨)
names(JN_)
cell.JN <- select(JN_, c(1, 14, 15, 21:36))
cell.JN <- run.flowsom(JN_, clustering.cols, meta.k = metak)

fwrite(cell.JN, "Clustered.csv")

### Run DR
head(cell.JN)
cell.sub <- do.subsample(cell.JN, sub.targets, sub.by)
#run UMAP
cell.sub <- Spectre::run.umap(dat = cell.sub, use.cols = clustering.cols)

fwrite(cell.sub, "RD_umap.csv")

### Make expression heatmap

exp <- do.aggregate(cell.JN, cellular.cols, by = "FlowSOM_metacluster")
exp

make.pheatmap(exp, "FlowSOM_metacluster", plot.cols = cellular.cols)



### Make expression plots
make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor')
make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor')

make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)

for(i in cellular.cols){
  make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", i, group.col, figure.title = paste0('Multiplot - ', i, '.png'))
}




##########################################################################################################
#### Annotate clusters and write summary data
##########################################################################################################
getwd()
setwd(OutputDirectory)
dir.create("reannot")
setwd("reannot")

### Identify cellular populations
#우선 임의로 annotation하지 않기. 20211228 여기까지만 했다.
annots <- list('T cells' = c(4),
               'B cells' = c(2),
               'Tumor' = c(1),
               'APC' = 5)

annots <- do.list.switch(annots)
setorderv(annots, 'Population')
annots

### Add population names to datasets

names(annots) <- c('Values', "Population")
annots

cell.dat <- do.add.cols(cell.dat.f, "FlowSOM_metacluster", annots, "Values")
cell.dat

cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
cell.sub

### Save annotated data and make population plots

fwrite(cell.dat, "Clustered.annotated.csv") 
fwrite(cell.sub, "DR.sub.annotated.csv")

### Population plots

annot.exp <- do.aggregate(cell.dat, cellular.cols, by = "Population")
make.pheatmap(annot.exp, "Population", plot.cols = cellular.cols)

make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Values", col.type = 'factor', add.label = TRUE)
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')


##########################################################################################################
#### Percent positive plots and cutoffs
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 3.4 - percent positive")
setwd("Output 3.4 - percent positive")

dir.create("Percent positive plots")
setwd("Percent positive plots")

### Plots for cytoff

as.matrix(names(cell.dat))

perc.pos.markers <- names(cell.dat)[c(13, 20)] #CD3랑 Pax5
perc.pos.markers

plot.against <- 'CD3_asinh'

for(i in perc.pos.markers){
  make.multi.plot(cell.sub, i, plot.against, plot.by = clustering.cols)
}

### Percent positive cutoffs

perc.pos.markers
perc.pos.cutoffs <- c(1, 0) #Pax5는 1에서 cut, CD3는 0에서 cut

as.matrix(perc.pos.markers)
as.matrix(perc.pos.cutoffs)

for(i in c(1:length(perc.pos.markers))){
  a <- perc.pos.markers[[i]]
  b <- perc.pos.cutoffs[[i]]
  
  x <- cell.sub[cell.sub[[a]] >= b]
  x$Pos <- TRUE
  y <- cell.sub[cell.sub[[a]] < b]
  y$Pos <- FALSE
  z <- rbind(x,y)
  
  make.multi.plot(z, a, plot.against, 'Pos', group.col)
}

### Generate data.table

perc.pos <- data.table(perc.pos.markers, perc.pos.cutoffs)
names(perc.pos) <- c('Columns', 'Cutoff')
perc.pos

##########################################################################################################
#### Write summary data
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 3.5 - summary data")
setwd("Output 3.5 - summary data")

### Select columns to measure MFI

as.matrix(cellular.cols)
dyn.cols <- cellular.cols[c(1,5)] #Pax5, CD3 
dyn.cols

### Setup cell count data

as.matrix(unique(cell.dat[[sample.col]]))
meta.dat

# counts <- meta.dat[,c(sample.col, 'Cells per sample'), with = FALSE]
# counts

### Create summary tables

sum.dat <- create.sumtable(dat = cell.dat, 
                           sample.col = sample.col,
                           pop.col = "Population",
                           use.cols = dyn.cols, 
                           annot.cols = c(group.col, batch.col), 
                           #counts = counts, 
                           perc.pos = perc.pos)

as.matrix(names(sum.dat))
sum.dat    

### Write summary data

fwrite(sum.dat, 'sum.dat.csv')

##########################################################################################################
#### Output session info
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output - info")
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()


