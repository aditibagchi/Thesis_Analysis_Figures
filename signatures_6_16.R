
MBREC6_R_Sigs_input <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MBREC6_R_Sigs_input.csv")
MBREC16_P_Sigs_input <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MBREC16_P_Sigs_input.csv")

MBREC6_P_DF_Sigs_input <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MBREC6_P_DF_Sigs_input.csv", header=TRUE)
sigs.input.MB16_P = mut.to.sigs.input(mut.ref = MBREC16_P_Sigs_input, sample.id = "sample",chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = BSgenome.Hsapiens.UCSC.hg19)
MB6_R <- whichSignatures(tumor.ref = sigs.input.MB6_R, sample.id = "MB-REC-06_recurrent_g1.crossmap.hg38_to_hg19",
                signatures.ref = signatures.nature2013, associated = c(),
                signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                tri.counts.method = "default")
plotSignatures(MB6_R)
mut.to.sigs.input(MBREC16_P_Sigs_input, sample.id = "sample", chr = "chr", pos = "pos",
                  ref = "ref", alt = "alt", bsg = NULL)