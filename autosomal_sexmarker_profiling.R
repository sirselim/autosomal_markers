#### 
## date created: 06-06-2014
## last modified: 10-06-2014

# load required packages
library('marmalaid')
library('minfi')
library('FactoMineR')
library('glmnet')

load('cd14_sexmarker_profiling.RData')
load('cd14_betas.RData')
load('cd3_betas.RData')

# get sample annotation from public database
# data(annotation_v1.1)
# read in chip annotation
# meth.anno <- read.csv('humanmethylation450_15017482_v1-2.csv', head = T, skip = 7, as.is = T)
# remove controls and SNPs from annotation
meth.anno2 <- meth.anno[c(1:485462),]
# look at probe distribution
# probe.distribution <- table(meth.anno2$CHR)
# barplot(probe.distribution, ylab = 'Number of probes', xlab = 'Chromosome', ylim = c(0, 50000))

# get probes that map multiple times and filter them out
multi.map <- read.csv('HumanMethylation450_15017482_v.1.1_hg19.mult.txt', head = F, as.is = T)
multi.map <- as.character(as.matrix(multi.map))
# added in cross-reactive probes from recent study
probes.crossReact <- read.csv('48639-non-specific-probes-Illumina450k.csv', head = T, as.is = T)
probes.crossReact <- as.character(probes.crossReact$TargetID)
# check numbers
table(meth.anno2$IlmnID %in% multi.map)
# remove probes that map multiple times
probe.filter <- meth.anno2$IlmnID %in% multi.map
meth.anno2 <- meth.anno2[!probe.filter,]
# ADDITIONAL: remove probes that have been shown to cross-react
probe.filter <- meth.anno2$IlmnID %in% probes.crossReact
meth.anno2 <- meth.anno2[!probe.filter,]


cd14.samples <- annotation[grep('CD14', annotation$TISSUE_SUBTYPE),]
all.probes <- as.character(meth.anno2$IlmnID)
write.csv(all.probes, 'filteredProbes_450k.csv', row.names = F)

## the below getbeta function is really slow if you want to get all probes,
## better to just use it for getting smaller sets of probes (still not ideal).
## to get around this I'll get all the information from the local database
## on the ESR server methead
# cd14.beta <- getbeta(cd14.samples$Id, probes = all.probes)

# now loaded up front
# load('cd14_betas.RData')

# clean up sample details
cd14.samples$SEX <- gsub(' female', 'Female', cd14.samples$SEX)
cd14.samples$SEX <- gsub(' male', 'Male', cd14.samples$SEX)
cd14.samples$SEX <- gsub(' Female', 'Female', cd14.samples$SEX)
cd14.samples$SEX <- gsub(' NA', 'NA', cd14.samples$SEX)

cd14.samples$SEX <- gsub('NA', NA, cd14.samples$SEX)
cd14.samples[is.na(cd14.samples$SEX),]$SEX <- 'unknown'

cd14.samples.sex <- cd14.samples[c('Id', 'SEX')]

mdsPlot(cd14.beta, numPositions = 5000, sampGroups = cd14.samples$SEX)

# the above looks a bit strange, check with just X chr markers
table(meth.anno2$CHR == 'X')

x.probes <- meth.anno2[meth.anno2$CHR == 'X',]$IlmnID
cd14.beta.x <- cd14.beta[rownames(cd14.beta) %in% x.probes,]

mdsPlot(cd14.beta[rownames(cd14.beta) %in% x.probes,], numPositions = 10000, sampGroups = cd14.samples$SEX)
mdsPlot(cd14.beta.x, sampGroups = cd14.samples$SEX, legendPos = "topleft", numPositions = 10000)

cd14.x.pca <- PCA(t(cd14.beta.x))

cd14.pca.coords <- data.frame(Id = rownames(cd14.x.pca$ind$coord), cd14.x.pca$ind$coord)

merge(cd14.samples.sex, cd14.pca.coords, by = 'Id', sort = F)

# plot.PCA(cd14.x.pca, axes = c(1,2))
# plot.PCA(cd14.x.pca, axes = c(1,3))
# plot.PCA(cd14.x.pca, axes = c(2,3))
# plot.PCA(cd14.x.pca, axes = c(3,4))

## NOTE: it looks like there are issues of mixed gender within the CD14 data set.
## We can use the samples that have a correct call based on X chr markers, and then
## once we have a 'trained' model could try and determine the 'unknowns'/mixed samples.


## trying CD3 data, there should be 6 male and 6 female
cd3.samples <- annotation[grep('CD3', annotation$TISSUE_SUBTYPE),]
cd3.samples <- tail(cd3.samples, n = 12)
cd3.samples$SEX <- gsub(' ', '', cd3.samples$SEX)

# below is now loaded at the front end
# load('cd3_betas.RData')

mdsPlot(cd3.beta, numPositions = 5000, sampGroups = cd3.samples$SEX, legendPos = "topright")
mdsPlot(cd3.beta[rownames(cd3.beta) %in% x.probes,], numPositions = 10000, sampGroups = cd3.samples$SEX, legendPos = "topright")


table(meth.anno2$CHR)
sexChr.probes <- meth.anno2[meth.anno2$CHR == 'X' | meth.anno2$CHR == 'Y',]$IlmnID

autosomal.filter <- meth.anno2$IlmnID %in% sexChr.probes
autosomal.probes <- meth.anno2[!autosomal.filter,]$IlmnID

mdsPlot(cd3.beta[rownames(cd3.beta) %in% autosomal.probes,], numPositions = 100000, 
        sampGroups = cd3.samples$SEX, legendPos = "topright")


colnames(cd3.beta) <- paste(colnames(cd3.beta), cd3.samples$SEX, sep = '_')

## testing glmnet
cd3.samples$SEX2 <- ifelse(cd3.samples$SEX == 'female', 1, 2)

x <- t(cd3.beta)
y <- cd3.samples$SEX2

fit1 <- glmnet(x, y)
fit1
plot(fit1)
plot(fit1, xvar="lambda")

cv.fit1 <- cv.glmnet(x, y)  # Perform cross validation on the fited model
plot(cv.fit1)

tpred <- predict(fit1, x)
coef(fit1, s=0.01)

mte <- apply((tpred - y)^2, 2, mean)  # Compute mse for the predictions
points(log(fit1$lambda), mte, col="blue", pch="*")	# overlay the mse predictions on the plot
legend("topleft", legend = c("10 fold CV", "Test"), pch="*", col = c("red", "blue"))


# fit1m <- glmnet(x, y, family="mgaussian")
# plot(fit1m, type.coef="2norm")

head(coef(fit1, s=0.01), n = 100)
cd3.beta[grep('cg03278611', rownames(cd3.beta)),]
meth.anno[grep('cg03278611', meth.anno$IlmnID),]

cd3.beta[grep('cg05964935', rownames(cd3.beta)),]
meth.anno[grep('cg05964935', meth.anno$IlmnID),]

test <- coef(fit1, s=0.01)
test@i
cd3.beta[c(58, 1109),]

meth.anno[grep('cg01804836', meth.anno$IlmnID),]
cd3.beta[test@i,]
rownames(cd3.beta[test@i,])

meth.anno[meth.anno$IlmnID %in% rownames(cd3.beta[test@i,]),]$CHR
PCA(t(cd3.beta[test@i,]))
## awesome! The above model appears to have worked and has pulled out sex chromosome markers! good start.

cd3.autosomal <- cd3.beta[rownames(cd3.beta) %in% autosomal.probes,]
cd3.samples$SEX2 <- ifelse(cd3.samples$SEX == 'female', 1, 0)

x <- t(cd3.autosomal)
y <- cd3.samples$SEX2

fit2 <- glmnet(x, y)
fit2
plot(fit2)
plot(fit2, xvar="lambda")

cv.fit2 <- cv.glmnet(x, y)  # Perform cross validation on the fited model
plot(cv.fit2)

tpred <- predict(fit2, x)
test3 <- coef(fit2, s=0.01)

mte <- apply((tpred - y)^2, 2, mean)  # Compute mse for the predictions
points(log(fit2$lambda), mte, col="blue", pch="*")  # overlay the mse predictions on the plot
legend("topleft", legend = c("10 fold CV", "Test"), pch="*", col = c("red", "blue"))

cd3.autosomal[test3@i,]
PCA(t(cd3.autosomal[test3@i,]))

cd3.beta[grep('cg00341297', rownames(cd3.beta)),]
meth.anno[grep('cg00341297', meth.anno$IlmnID),]

# meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal[test3@i,]),]
# genes <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal[test3@i,]),]$UCSC_RefGene_Name
# 
# genes <- unique(unlist(strsplit(genes, split = ';')))


cd3.autosome.markers <- data.frame(IlmnID = rownames(cd3.autosomal[test3@i,]),
                                   female.mean = rowMeans(cd3.autosomal[test3@i, grep('female', colnames(cd3.autosomal))]),
                                   male.mean = rowMeans(cd3.autosomal[test3@i, grep('_male', colnames(cd3.autosomal))]))

cd3.autosome.markers$abs_diff <- abs(cd3.autosome.markers$female.mean - cd3.autosome.markers$male.mean)
cd3.autosome.markers$percent_diff <- round(cd3.autosome.markers$abs_diff * 100, 2)
cd3.autosome.markers$CHR <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal[test3@i,]),]$CHR
cd3.autosome.markers$gene <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal[test3@i,]),]$UCSC_RefGene_Name

gene.filter <- matrix(NA,18)

for(i in c(1:nrow(cd3.autosome.markers))) {
  
  gene.filter[i,] <- paste(unique(strsplit(cd3.autosome.markers$gene, split = ';')[[i]]), collapse = "; ")
    
}

cd3.autosome.markers$gene <- gene.filter

write.csv(cd3.autosome.markers, 'cd3_autosomalMarkers_20140606.csv', row.names = F)



# ## can I use the model to predict sex in the cd14+ samples?
# cd14.autosomal <- cd14.beta[rownames(cd14.beta) %in% autosomal.probes,]
# colnames(cd14.autosomal) <- paste(colnames(cd14.autosomal), cd14.samples$SEX, sep = '_')
# 
# cd14.predict <- predict(fit2, t(cd14.autosomal))
# 
# cd14.autosomal[test3@i,]
# 
# cd14.autosomal[grep('cg16218221', rownames(cd14.autosomal)),]
# 
# table(rownames(cd14.autosomal[test3@i,]) %in% probes.crossReact$TargetID)
# 
# ### rerun the cd3 model with cross-reactive probes removed
# cd3.autosomal <- cd3.beta[rownames(cd3.beta) %in% autosomal.probes,]
# cd3.samples$SEX2 <- ifelse(cd3.samples$SEX == 'female', 1, 0)
# filter.crossReact <- rownames(cd3.autosomal) %in% probes.crossReact$TargetID
# cd3.autosomal.2 <- cd3.autosomal[!filter.crossReact,]
# 
# 
# x <- t(cd3.autosomal.2)
# y <- cd3.samples$SEX2
# 
# fit2 <- glmnet(x, y)
# fit2
# plot(fit2)
# plot(fit2, xvar="lambda")
# 
# cv.fit2 <- cv.glmnet(x, y)  # Perform cross validation on the fited model
# plot(cv.fit2)
# 
# tpred <- predict(fit2, x)
# test3 <- coef(fit2, s=0.01)
# 
# mte <- apply((tpred - y)^2, 2, mean)  # Compute mse for the predictions
# points(log(fit2$lambda), mte, col="blue", pch="*")  # overlay the mse predictions on the plot
# legend("topleft", legend = c("10 fold CV", "Test"), pch="*", col = c("red", "blue"))
# 
# cd3.autosomal.2[test3@i,]
# PCA(t(cd3.autosomal.2[test3@i,]))
# 
# cd3.beta[grep('cg00341297', rownames(cd3.beta)),]
# meth.anno[grep('cg00341297', meth.anno$IlmnID),]
# 
# meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal[test3@i,]),]
# genes <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal.2[test3@i,]),]$UCSC_RefGene_Name
# 
# # genes <- unique(unlist(strsplit(genes, split = ';')))
# 
# 
# cd3.autosome.markers <- data.frame(IlmnID = rownames(cd3.autosomal.2[test3@i,]),
#                                    female.mean = rowMeans(cd3.autosomal.2[test3@i, grep('female', colnames(cd3.autosomal.2))]),
#                                    male.mean = rowMeans(cd3.autosomal.2[test3@i, grep('_male', colnames(cd3.autosomal.2))]))
# 
# cd3.autosome.markers$abs_diff <- abs(cd3.autosome.markers$female.mean - cd3.autosome.markers$male.mean)
# cd3.autosome.markers$percent_diff <- round(cd3.autosome.markers$abs_diff * 100, 2)
# cd3.autosome.markers$CHR <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal.2[test3@i,]),]$CHR
# cd3.autosome.markers$gene <- meth.anno[meth.anno$IlmnID %in% rownames(cd3.autosomal.2[test3@i,]),]$UCSC_RefGene_Name
# 
# gene.filter <- matrix(NA,length(genes))
# 
# for(i in c(1:nrow(cd3.autosome.markers))) {
#   
#   gene.filter[i,] <- paste(unique(strsplit(cd3.autosome.markers$gene, split = ';')[[i]]), collapse = "; ")
#   
# }
# 
# cd3.autosome.markers$gene <- gene.filter
# 
# # write.csv(cd3.autosome.markers, 'cd3_autosomalMarkers_20140606.csv', row.names = F)
# 
# cd14.autosomal[grep('cg11643285', rownames(cd14.autosomal)),]
# 
# 
# cd14.test <- data.frame(beta = cd14.autosomal[grep('cg11643285', rownames(cd14.autosomal)),], sex = cd14.samples$SEX)
# boxplot(cd14.test$beta ~ cd14.test$sex)
# 
# 
# 
# cd3.test <- data.frame(beta = cd3.autosomal[grep('cg11643285', rownames(cd3.autosomal)),], sex = cd3.samples$SEX)
# boxplot(cd3.test$beta ~ cd3.test$sex)

# # liver.samples <- annotation[grep('Liver', annotation$TISSUE),]
# # liver.samples <- liver.samples[liver.samples$DISEASE == 'Healthy',]
# # 
# # as.character(na.omit(liver.samples$Id))
# # get more data
# load('cd34_beta.RData')
# load('cd19_beta.RData')
# load('liver_beta.RData')
# 
# 
# par(mfrow = c(1,3))
# probe <- 'cg13056653'
# cd14.test <- data.frame(beta = cd14.autosomal[grep(probe, rownames(cd14.autosomal)),], sex = cd14.samples$SEX)
# boxplot(cd14.test$beta ~ cd14.test$sex, main = paste('CD14+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd3.test <- data.frame(beta = cd3.autosomal[grep(probe, rownames(cd3.autosomal)),], sex = cd3.samples$SEX)
# boxplot(cd3.test$beta ~ cd3.test$sex, main = paste('CD3+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd19.test <- data.frame(beta = cd19.beta[grep(probe, rownames(cd19.beta)),], 
#                         sex = annotation[grep('CD19', annotation$TISSUE_SUBTYPE),]$SEX)
# cd19.test$sex <- as.character(cd19.test$sex)
# cd19.test[is.na(cd19.test)] <- 'unknown'
# boxplot(cd19.test$beta ~ cd19.test$sex, main = paste('CD19+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# 
# 
# par(mfrow = c(1,3))
# probe <- 'cg11232448'
# cd14.test <- data.frame(beta = cd14.autosomal[grep(probe, rownames(cd14.autosomal)),], sex = cd14.samples$SEX)
# boxplot(cd14.test$beta ~ cd14.test$sex, main = paste('CD14+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd3.test <- data.frame(beta = cd3.autosomal[grep(probe, rownames(cd3.autosomal)),], sex = cd3.samples$SEX)
# boxplot(cd3.test$beta ~ cd3.test$sex, main = paste('CD3+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd19.test <- data.frame(beta = cd19.beta[grep(probe, rownames(cd19.beta)),], 
#                         sex = annotation[grep('CD19', annotation$TISSUE_SUBTYPE),]$SEX)
# cd19.test$sex <- as.character(cd19.test$sex)
# cd19.test[is.na(cd19.test)] <- 'unknown'
# boxplot(cd19.test$beta ~ cd19.test$sex, main = paste('CD19+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# 
# 
# par(mfrow = c(1,4))
# probe <- 'cg11643285'
# cd14.test <- data.frame(beta = cd14.autosomal[grep(probe, rownames(cd14.autosomal)),], sex = cd14.samples$SEX)
# boxplot(cd14.test$beta ~ cd14.test$sex, main = paste('CD14+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd3.test <- data.frame(beta = cd3.autosomal[grep(probe, rownames(cd3.autosomal)),], sex = cd3.samples$SEX)
# boxplot(cd3.test$beta ~ cd3.test$sex, main = paste('CD3+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd19.test <- data.frame(beta = cd19.beta[grep(probe, rownames(cd19.beta)),], 
#                         sex = annotation[grep('CD19', annotation$TISSUE_SUBTYPE),]$SEX)
# cd19.test$sex <- as.character(cd19.test$sex)
# cd19.test[is.na(cd19.test)] <- 'unknown'
# boxplot(cd19.test$beta ~ cd19.test$sex, main = paste('CD19+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))
# 
# cd34.test <- data.frame(beta = cd34.beta[grep(probe, rownames(cd34.beta)),], 
#                         sex = annotation[grep('CD34', annotation$TISSUE_SUBTYPE),]$SEX)
# cd34.test$sex <- as.character(cd34.test$sex)
# cd34.test[is.na(cd34.test)] <- 'unknown'
# cd34.test[cd34.test$sex == ' NA',]$sex <- 'unknown'
# cd34.test$predict.sex <- ifelse(cd34.test$beta >= 0.88, 'female', 'male')
# boxplot(cd34.test$beta ~ cd34.test$predict.sex, main = paste('CD19+:', probe, sep = ' '), ylab = 'beta', ylim = c(0,1))


