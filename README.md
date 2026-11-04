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
