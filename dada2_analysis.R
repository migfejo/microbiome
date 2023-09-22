
## Carga de paquetes necesarios
library(dada2)
library(phyloseq)

## Lectura de archivos y primera evaluación
fnFs <- sort(list.files(path = "fastq", pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path = "fastq", pattern="_R2_001.fastq.gz", full.names = TRUE))
# Se guardan los identificadores de las muestras
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(sample(fnFs, 1))
plotQualityProfile(sample(fnRs, 1))


## Filtrado y truncado
filtFs <- file.path("myfilt", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("myfilt", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,240), maxN=0,
                     maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE,
                     multithread=FALSE)
head(out)

sort(sapply(fnFs, function(x){sum(getUniques(x))}))  # Lecturas por muestra originales

plotQualityProfile(sample(filtFs, 1))
plotQualityProfile(sample(filtRs, 1))


## Aprender tasas de error
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## Inferencia de ASVs
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaFs[[1]]
dadaRs[[1]]



## Unión de secuencias
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspeccion del objeto
head(mergers[[1]])



## Construccion de la tabla de abundancias (Desprelicacion)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Distribución de la longitud de secuencias
table(nchar(getSequences(seqtab)))



## Eliminar quimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

# Ratio de secuencias supervivientes
ncol(seqtab.nochim)/ncol(seqtab)
# Ratio de lecturas supervivientes
sum(seqtab.nochim)/sum(seqtab)


## Asignar taxonomía
taxa <- assignTaxonomy(seqtab.nochim, 
                       paste0("/tax/silva_nr_v132_train_set.fa.gz"),
                       multithread=FALSE)




