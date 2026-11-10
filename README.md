# Bioc-ERBS-ChIP-Seq-Analysis
# Harvard Introduction to Bioconductor eğitimi ile birlikte yapmış olduğum analiz. Bu çalışmanın amacı genomik veri analizi ve anotasyon yapmak üzerine, dir. Hücre hattı olarak HepG2 ve GM12878 kullanılmaktadır. HepG2; hepatosellüler karsinomundan izole edilmiş, epitel benzeri morfolojiye sahiptir. GM12878;  Epstein-Barr virüsü (EBV) ile dönüştürülmüş bir insan lenfoblastoid hücre hattıdır. Meme kanseri tanısında önemli rol alan Östradiol Reseptörünün DNA üzerindeki bağlanma bölgelerini (ChIP-seq peak) Bioconductor paketleri kullanarak karşılaştırmak, ortak/hücreye özgü bölgeleri tespit etmek ve gen anotasyon veritabanları yardımıyla hedef genleri ortaya çıkarmaktır.
# Eğitimle birlikte yaparken hatalar çıkmaktadır. Bu çalışma hatalardan temizlenmiş kısım ile oluşturulmaktadır.


if( !requireNamespace("BiocManager", quietly = TRUE ))
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "Biobase", "Biostrings"))
if (!requireNamespace("remotes", quietly = TRUE))
install.packages("remotes")
remotes::install_github("genomicsclass/tissuesGeneExpression")
  library(tissuesGeneExpression) #Bununla geri çağırıyoruz
BiocManager::install(c("GenomicRanges", "IRanges", "Biobase", "SummarizedExperiment",
                       "genefilter", "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "org.Hs.eg.db", "airway"), update = FALSE)
remotes::install_github("genomicsclass/ERBS")
remotes::install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression); data(tissuesGeneExpression)
library(ERBS)
library(GenomicRanges)
library(Biobase)
data(HepG2); data(GM12878)
 #Tepe noktasını bulmak için
 length(HepG2); length(GM12878)
 #İkisinin ortak tepelerini bulmak için
 ov <- findOverlaps(HepG2,GM12878)
 length(unique(queryHits(ov)))
 mat <- matrix(rnorm(20), nrow = 5,
              dimnames = list(paste0("gene", 1:5), paste0("s", 1:4)))
pheno <- data.frame(group = c("ctrl", "ctrl", "case", "case"),
                    row.names = colnames(mat))
eset <- ExpressionSet(assayData = mat,
                      phenoData = AnnotatedDataFrame(pheno))
eset
all(colnames(exprs(eset)) == rownames(pData(eset)))
BiocManager::install("airway")
library(airway)
data(airway)
airway
assay(airway)[1:3, 1:4] #okuma sayıları
colData(airway)        # hangi örnek tedavi gördü, hangisi görmedi
#ExperimenSet'te konumlar verilemz sadece isimleri vardır ancak SummarizedExperiment her satıra bir genomik aralık eklenebilir.
#Aynı genin birden fazla adı vardır bu yüzden bunları tek isime çevirmek için aşağıdakiler kullanılır.
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)         # elindeki anahtar türleri
columns(org.Hs.eg.db)          # isteyebileceğin bilgiler
# birden fazla bilgi tek tabloda
select(org.Hs.eg.db, keys = c("TP53", "BRCA1"),
       columns = c("ENTREZID", "ENSEMBL", "GENENAME"),
       keytype = "SYMBOL")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

g  <- genes(txdb)                       # gen aralıkları (GRanges)
tx <- transcripts(txdb)                 # transkriptler
ex <- exonsBy(txdb, by = "gene")        # her genin eksonları (GRangesList)
pr <- promoters(g, upstream = 2000, downstream = 200)   # promotör bölgeleri
library(ERBS); data(HepG2)

# tepelerin promotöre denk gelenleri
subsetByOverlaps(HepG2, pr)

# her tepenin en yakın geni
d <- distanceToNearest(HepG2, g)
d
mcols(d)$distance           # uzaklıklar (baz çifti)

# en yakın genin Entrez numarasından sembolüne
nearest <- g[subjectHits(d)]
mapIds(org.Hs.eg.db, keys = nearest$gene_id,
       column = "SYMBOL", keytype = "ENTREZID")
library(genefilter)
library(tissuesGeneExpression); data(tissuesGeneExpression)

ind <- tissue %in% c("kidney", "liver")
tt  <- rowttests(e[, ind], factor(tissue[ind]))   # her probe için t-testi

sum(tt$p.value < 0.05)                            # düzeltmesiz
sum(p.adjust(tt$p.value, "bonferroni") < 0.05)   # Bonferroni
sum(p.adjust(tt$p.value, "BH") < 0.05)           # FDR
# Beklenen: 100 * 500 / 20000 = 2.5
tab <- matrix(c(15, 485, 85, 19415), nrow = 2,
              dimnames = list(c("kumede", "kumede_degil"),
                              c("anlamli", "anlamsiz")))
tab
fisher.test(tab)
