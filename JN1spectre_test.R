#이렇게 저렇게 해보는 test용

##MHCI이 0인 cell은 없앰
library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages
library(dplyr)

#"F:/MICSSS_Analysis/giotto/JN1_spectregiotto"
##area에서 50보다 작거나 920보다 큰 cell은 제거
head(raw)
raw[(raw$area < 50), ] = 0
raw[(raw$area > 920), ] = 0
raw2 <- raw[!(raw$area == 0), ]
summary(raw2)

raw3 <- raw2[!(raw2$MHCI < 0), ]
summary(raw3)

names(raw3)
rawE <- select(raw3, c(2, 4:13)) #expression file
names(rawE)

rawL <- select(raw3, c(2, 14, 15)) #location file
names(rawL)

rawM <- select(raw3, c(2, 23, 24)) #metadata file
names(rawM)

fwrite(rawE, 'expression.csv')
fwrite(rawL, 'location.csv')
fwrite(rawM, 'metadata.csv')

##########################################################################################################
#### Asinh transformation
##########################################################################################################

setwd(OutputDirectory)
dir.create("/Output 1.1 - transformation plots4")
setwd("/Output 1.1 - transformation plots4")
getwd()
### Co-factor targets

cf.low <- 0.1
cf.mid <- 500
cf.high <- 1000

cf.low
cf.mid
cf.high

### Transformation settings

as.matrix(names(rawE))

asinh.low <- names(rawE)[c(2:11)]
asinh.mid <- names(cell.dat)[c(1:11)]
asinh.high <- names(cell.dat)[c()]

asinh.low
asinh.mid
asinh.high

### Test transformation settings on subsampled data

sub <- do.subsample(rawE, 20000)

sub <- do.asinh(sub, use.cols = asinh.low, cofactor = cf.low)
sub <- do.asinh(sub, use.cols = asinh.mid, cofactor = cf.mid)
sub <- do.asinh(sub, use.cols = asinh.high, cofactor = cf.high)

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

cell.dat <- do.asinh(rawE, use.cols = asinh.low, cofactor = cf.low)
cell.dat <- do.asinh(cell.dat, use.cols = asinh.mid, cofactor = cf.mid)
cell.dat <- do.asinh(cell.dat, use.cols = asinh.high, cofactor = cf.high)

as.matrix(names(cell.dat))

transf.cols <- names(cell.dat)[grepl('_asinh', names(cell.dat))]
as.matrix(transf.cols)

##########################################################################################################
#### Add metadata
##########################################################################################################

meta.dat

sample.info <- rawM
sample.info

### Add sample metadata to primary data.table

cell.dat <- do.add.cols(cell.dat, "CellID", sample.info, "CellID", rmv.ext = TRUE)
head(cell.dat)

##########################################################################################################
#### Write data to disk
##########################################################################################################
getwd()
setwd("F:/hjikim/flowjo/211223 JN1/")
dir.create("Output 1.2 - transformed data")
setwd("Output 1.2 - transformed data")
getwd()
### Write cellular data and analysis  preferences to disk

fwrite(cell.dat, "cell.dat.csv") # data

### Save session info to disk

setwd(OutputDirectory)
dir.create("Output - info", showWarnings = FALSE)
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()

######여기까지 데이터 전처리 이제 clustering
### Set output directory
getwd()
setwd(PrimaryDirectory)
dir.create("Output 3 - clustering and DR", showWarnings = FALSE)
setwd("Output 3 - clustering and DR")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

##########################################################################################################
#### Setup preferences
##########################################################################################################
cell.datori <- cell.dat
cell.dat <- cell.dat3
### Sample preferences

as.matrix(names(cell.dat))

sample.col <- "core"
group.col <- "Group"
batch.col <- "Batch"

### Downsample preferences
x <- rep("JN1", 503671)
x <- rep(1, 503671)
cell.dat2 <- cbind(cell.dat2, x)

sub.by <- "Group"
names(cell.dat2)
names(cell.dat2)[25] <- "Batch"
as.matrix(unique(cell.dat2[[sub.by]]))

data.frame(table(cell.dat2[[group.col]]))

sub.targets <- c(3632, 3632)
sub.targets

### Clustering preferences

## Cellular cols
as.matrix(names(cell.dat))

cellular.cols <- names(cell.dat)[c(13,15,16,19,20)]
cellular.cols

## Columns for clustering
as.matrix(names(cell.dat))

clustering.cols <- names(cell.dat)[c(13,15,16,19,20)]
clustering.cols

## Cluster numbers etc
metak <- 'auto'

##########################################################################################################
#### Run clustering and dimensionality reduction
##########################################################################################################
#밑에서 구했던 cutoff값으로 자르기 -> cell 수 3600개 남았고.. 제대로 clustering안됨
cell.off <- cell.dat[!(cell.dat$CD3_asinh < 0), ]
cell.off <- cell.off[!(cell.off$Pax5_asinh <1), ]

setwd(OutputDirectory)
dir.create("Output 3.1 - clustered2")
setwd("Output 3.1 - clustered2")
getwd()
### Run clustering

cell.dat.f2 <- run.flowsom(cell.off, clustering.cols, meta.k = metak)
cell.dat

fwrite(cell.dat.f2, "Clustered_cutoff.csv")

### Run DR


cell.sub2 <- do.subsample(cell.dat.f2, sub.targets, sub.by)
cell.sub2 <- run.umap(cell.sub2, clustering.cols)
cell.sub2

fwrite(cell.sub, "RD.cutoff.csv")

### Make expression heatmap

exp <- do.aggregate(cell.dat.f2, cellular.cols, by = "FlowSOM_metacluster")
exp

make.pheatmap(exp, "FlowSOM_metacluster", plot.cols = cellular.cols)

### Make expression plots

make.colour.plot(cell.sub2, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.colour.plot(cell.sub2, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor')
make.colour.plot(cell.sub2, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor')

make.multi.plot(cell.sub2, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
make.multi.plot(cell.sub2, "UMAP_X", "UMAP_Y", cellular.cols)

for(i in cellular.cols){
  make.multi.plot(cell.sub2, "UMAP_X", "UMAP_Y", i, group.col, figure.title = paste0('Multiplot - ', i, '.png'))
}

##########################################################################################################
#### Annotate clusters and write summary data
##########################################################################################################
getwd()
setwd(OutputDirectory)
dir.create("reannot")
setwd("reannot")

### Identify cellular populations

annots <- list('T cells' = c(3,4),
               'B cells' = c(6,7),
               'Tumor' = c(1,2),
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

fwrite(cell.dat, "Clustered.annotated.csv") #이 데이터로 giotto에 적용해보기
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


