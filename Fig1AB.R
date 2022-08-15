
library(maftools)

laml.maf <- "output/mutation.data.type.txt" 
laml.clin <- "output/clinical.data.spe.txt"
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

length(unique(laml@data[laml@data$Type == 'SNV', 'Hugo_Symbol']$Hugo_Symbol))
length(unique(laml@data$Hugo_Symbol))

unique(laml@data[laml@data$Type == 'SNV' & laml@data$Tumor_Sample_Barcode == 'T983', 'Hugo_Symbol']$Hugo_Symbol)

library(data.table)
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
clinical <- fread("output/clinical.data.spe.txt")

Gender <- cols[3:4]
names(Gender) <- as.character(unique(na.omit(clinical$Gender)))

Liver_cirrhosis <- c("gray70","black")
names(Liver_cirrhosis) <- sort(as.character(unique(na.omit(clinical$Liver_cirrhosis))))

Lymph_node_metastasis <- c("gray70","black")
names(Lymph_node_metastasis) <- sort(as.character(unique(na.omit(clinical$Lymph_node_metastasis))))

Tumor_thrombus <- c("gray70","black")
names(Tumor_thrombus) <- sort(as.character(unique(na.omit(clinical$Tumor_thrombus))))

TNM_stage <- c("gray90","gray46","gray36","black")
names(TNM_stage) <- sort(as.character(unique(na.omit(clinical$TNM_stage))))


ann_colors = list(Gender = Gender, TNM_stage=TNM_stage, Liver_cirrhosis=Liver_cirrhosis,
                  Lymph_node_metastasis=Lymph_node_metastasis,Tumor_thrombus=Tumor_thrombus)


pdf("figs/Oncoplot.pdf", width = 10, height = 8)
oncoplot(maf = subsetMaf(maf = laml, query = "Type == 'SNV'"), top=22,
         clinicalFeatures = c('Gender', 'TNM_stage','Liver_cirrhosis','Tumor_thrombus'),
         sortByAnnotation = FALSE, keepGeneOrder=FALSE,
         anno_height=1.5, gene_mar = 6, bgCol="gray90",
         annotationColor = ann_colors
)
dev.off()


pdf("figs/mut_interaction.pdf", width = 8, height = 6)
somaticInteractions(maf = subsetMaf(maf = laml, query = "Type == 'SNV'"), top = 22, pvalue = c(0.05, 0.1), fontSize = 0.75)
dev.off()
