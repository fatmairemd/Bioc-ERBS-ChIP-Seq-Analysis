# Bioc-ERBS-ChIP-Seq-Analysis
# Harvard Introduction to Bioconductor eğitimi ile birlikte yapmış olduğum analiz. Bu çalışmanın amacı genomik veri analizi ve anotasyon yapmak üzerine, dir. Hücre hattı olarak HepG2 ve GM12878 kullanılmaktadır. HepG2; hepatosellüler karsinomundan izole edilmiş, epitel benzeri morfolojiye sahiptir. GM12878;  Epstein-Barr virüsü (EBV) ile dönüştürülmüş bir insan lenfoblastoid hücre hattıdır. Meme kanseri tanısında önemli rol alan Östradiol Reseptörünün DNA üzerindeki bağlanma bölgelerini (ChIP-seq peak) Bioconductor paketleri kullanarak karşılaştırmak, ortak/hücreye özgü bölgeleri tespit etmek ve gen anotasyon veritabanları yardımıyla hedef genleri ortaya çıkarmaktır.
# Eğitimle birlikte yaparken hatalar çıkmaktadır. Bu çalışma hatalardan temizlenmiş kısım ile oluşturulmaktadır.

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

gerekli_paketler <- c("GenomicRanges", "IRanges", "ERBS", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
for (pkg in gerekli_paketler) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, force = TRUE)
  }
}

library(GenomicRanges)
library(IRanges)
library(ERBS)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
cikti_klasoru <- "proje_ciktilari"
if (!dir.exists(cikti_klasoru)) {
  dir.create(cikti_klasoru)
}
data(HepG2)
data(GM12878)
hepg2_gr <- GRanges(seqnames = seqnames(HepG2), ranges = ranges(HepG2), strand = strand(HepG2))
gm_gr     <- GRanges(seqnames = seqnames(GM12878), ranges = ranges(GM12878), strand = strand(GM12878))
cat("--- Veri Özeti ---\n")
cat("HepG2 Tepe Sayısı:", length(hepg2_gr), "\n")
cat("GM12878 Tepe Sayısı:", length(gm_gr), "\n")
ozet_df <- data.frame(
  Hucere_Hatti = c("HepG2", "GM12878"),
  Tepe_Sayisi = c(length(hepg2_gr), length(gm_gr)),
  Ortalama_Genislik = c(mean(width(hepg2_gr)), mean(width(gm_gr))),
  Genome_Surumu = c("hg19", "hg19")
)
write.csv(ozet_df, file.path(cikti_klasoru, "01_genel_ozet.csv"), row.names = FALSE)
chr_hepg2 <- seqnames(hepg2_gr)
chr_tablo <- table(chr_hepg2)

png(file.path(cikti_klasoru, "02_hepg2_kromozom_dagilimi.png"), width = 800, height = 600)
barplot(chr_tablo[1:24], las = 2, main = "HepG2 Kromozom Başına Tepe Dağılımı", col = "steelblue", ylab = "Tepe Sayısı")
dev.off()
png(file.path(cikti_klasoru, "03_tepe_genisligi_karsilastirma.png"), width = 600, height = 600)
boxplot(width(hepg2_gr), width(gm_gr), 
        names = c("HepG2", "GM12878"), 
        main = "Tepe Genişlikleri Dağılımı", 
        ylab = "Genişlik (bp)", 
        col = c("orange", "purple"))
dev.off()
ov <- findOverlaps(hepg2_gr, gm_gr)
ortak_tepe_sayisi <- length(unique(queryHits(ov)))
ozgu_hepg2_sayisi <- length(hepg2_gr) - ortak_tepe_sayisi
cat("Ortak tepe sayısı:", ortak_tepe_sayisi, "\n")
cat("Sadece HepG2'ye özgü tepe sayısı:", ozgu_hepg2_sayisi, "\n")
ortak_ozgu_df <- data.frame(
  Kategori = c("Ortak Tepeler", "Sadece HepG2 Ozgu"),
  Sayi = c(ortak_tepe_sayisi, ozgu_hepg2_sayisi)
)
write.csv(ortak_ozgu_df, file = "proje_ciktilari/04_ortak_ve_ozgu_tepeler.csv", row.names = FALSE)
cat("Anotasyon ve en yakın gen hesaplanıyor...\n")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
g <- suppressMessages(genes(txdb))
n_sinir <- min(500, length(hepg2_gr))
ornek_hepg2 <- hepg2_gr[1:n_sinir] 
d <- distanceToNearest(ornek_hepg2, g)
en_yakin_genler <- g[subjectHits(d)]
gen_sembolleri <- mapIds(org.Hs.eg.db, 
                         keys = en_yakin_genler$gene_id, 
                         column = "SYMBOL", 
                         keytype = "ENTREZID", 
                         multiVals = "first")
anotasyon_sonuc <- data.frame(
  Tepe_ID = paste0("Peak_", 1:n_sinir),
  Mesafe_bp = mcols(d)$distance,
  En_Yakin_Gen = gen_sembolleri
)

write.csv(anotasyon_sonuc, file.path(cikti_klasoru, "05_tepe_en_yakin_genler.csv"), row.names = FALSE)

cat("\nTebrikler! Tüm analiz adımları hatasız tamamlandı. Dosyalar 'proje_ciktilari' klasörüne kaydedildi.\n")

# Vol2
if( !requireNamespace("BiocManager", quietly = TRUE ))
+ install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "Biobase", "Biostrings"))
if (!requireNamespace("remotes", quietly = TRUE))
+     install.packages("remotes")
 remotes::install_github("genomicsclass/tissuesGeneExpression")
  library(tissuesGeneExpression) #Bununla geri çağırıyoruz
 BiocManager::install(c("GenomicRanges", "IRanges", "Biobase", "SummarizedExperiment", "ERBS", "tissuesGeneExpression", "genefilter","TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db", "airway"), update = FALSE, forse = TRUE)
 library(tissueGeneExpression)
 data(tissuesGeneExpression
 library(GenomicRanges)
 library(IRanges)
 data(HepG2)
 data(GM12878)
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
              dimnames = list(c("kümede", "kümede değil"),
                              c("anlamlı", "anlamsız")))
tab
fisher.test(tab)
