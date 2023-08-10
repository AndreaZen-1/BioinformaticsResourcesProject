# Zen Andrea Bioinformatics Resources Code

# 1. Load the data
load("Breast_Cancer.RData")



# 2. Extract protein-coding genes only
library("biomaRt")

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
att <- c("ensembl_gene_id", "gene_biotype")
flt <- c("ensembl_gene_id")
vals <- r_anno_df$gene_id

# get a list of all the genes starting from the ids
query <- getBM(attributes = att, filters = flt, values = vals, mart = ensembl)

# find the ids of the protein coding genes
query_pc <- query[which(query$gene_biotype == "protein_coding"),]

# create new dataframes with only the protein coding genes
r_anno_pc <- r_anno_df[which(r_anno_df$gene_id %in% query_pc$ensembl_gene_id),]
raw_counts_pc <- raw_counts_df[query_pc$ensembl_gene_id,]  # This gives us NAs, we need to remove them!
raw_counts_pc <- na.omit(raw_counts_pc)

# making subplots of raw counts in different conditions with shared Y axis
fig1 <- plot_ly(x = colSums(raw_counts_df)/1000000,
                y = colnames(raw_counts_df),
                type = 'bar',
                name = "All genes",
                marker = list(color = '#9ECAE1', 
                              line = list(color = '#3A8EBC',
                                          width = 1.5)))

fig2 <- plot_ly(x = colSums(raw_counts_pc)/1000000,
                y = colnames(raw_counts_pc),
                type = 'bar',
                name = "Protein-coding genes",
                marker = list(color = '#E19F9D',
                              line = list(color = '#C74E4A',
                                          width = 1.5)))

fig <- subplot(fig1, fig2, shareY = TRUE) %>% layout(title = "Millions of reads per sample")
fig



# 3. Differential Expression Analysis
library("limma")
library("edgeR")
library("ggplot2")
library("tidyr")

cnt_thr <- 20
repl_thr <- 2
filter_vec <- apply(raw_counts_pc,1,function(y) max(by(y, c_anno_df$condition, function(x) sum(x>=cnt_thr))))
filter_counts_df <- raw_counts_pc[filter_vec>=repl_thr,]

# I would expect it to work, but this only returns NAs
#filter_anno_df <- r_anno_pc[rownames(filter_counts_df),]
# but we can always use %in%
filter_anno_df <- r_anno_pc[which(r_anno_pc$gene_id %in% rownames(filter_counts_df)),]

# Same plot as before but added also the filtered genes
fig1 <- plot_ly(x = colSums(raw_counts_df)/1000000,
                y = colnames(raw_counts_df),
                type = 'bar',
                name = "All genes",
                marker = list(color = '#9ECAE1', 
                              line = list(color = '#3A8EBC',
                                          width = 1.5)),
                width = 800, height = 500)

fig2 <- plot_ly(x = colSums(raw_counts_pc)/1000000,
                y = colnames(raw_counts_pc),
                type = 'bar',
                name = "Protein-coding genes",
                marker = list(color = '#E19F9D',
                              line = list(color = '#C74E4A',
                                          width = 1.5)),
                width = 800, height = 500)

fig3 <- plot_ly(x = colSums(filter_counts_df)/1000000,
                y = colnames(filter_counts_df),
                type = 'bar',
                name = "Filtered genes",
                marker = list(color = '#9DE1CC',
                              line = list(color = '#33A783',
                                          width = 1.5)),
                width = 800, height = 500)

fig <- subplot(fig1, fig2, fig3, shareY = TRUE) %>% layout(title = "Millions of reads per sample")
fig

# print the amounts of genes after the different filtering steps
cat(" All, raw genes: ", dim(raw_counts_df)[1], "\n",
    "Protein-coding genes: ", dim(raw_counts_pc)[1], "\n",
    "Filtered, protein-coding genes: ", dim(filter_counts_df)[1])


### PCA analysis of the data to check differences between Control and Cases

data.matrix <- filter_counts_df
# rank. sets the maximum rank to compute, useful if we have way more
prin_comp <- prcomp(t(data.matrix), rank. = 2, scale. = TRUE)
#prin_comp <- prcomp(t(data.matrix), rank. = 2)
components <- prin_comp[["x"]]              # get the components values
components <- data.frame(components)        # make a data.frame obj of them
components <- cbind(components, c_anno_df)  # add to the df sample and condition

# make a palette to have red -> Case,  blue -> Control
pal <- c("#E19F9D", "#9ECAE1")
pal <- setNames(pal, c("Case", "Control"))

d_pca <- plot_ly(components, x = ~PC1, y = ~PC2, color = ~condition,
                 colors = pal, type = 'scatter', mode = 'markers',
                 # Hover text:
                 text = ~paste("Sample: ", sample, "<br>Condition:",
                               condition),
                 width = 800)

# title with variance explained
tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio[1:2])
tit = paste("Total Explained Variance = ", tot_explained_variance_ratio)
d_pca <- d_pca %>% layout(title = tit, legend = list(title = list(text = "Condition")))


## Test of PCA plot with 3 Principal components, in 3D
data.matrix <- filter_counts_df
# rank. sets the maximum rank to compute, useful if we have way more
prin_comp <- prcomp(t(data.matrix), rank. = 3, scale. = TRUE)
#prin_comp <- prcomp(t(data.matrix), rank. = 3)
components <- prin_comp[["x"]]              # get the components values
components <- data.frame(components)        # make a data.frame obj of them
components <- cbind(components, c_anno_df)  # add to the df sample and condition

# make a palette to have red -> Case,  blue -> Control
pal <- c("#E19F9D", "#9ECAE1")
pal <- setNames(pal, c("Case", "Control"))

t_pca <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition,
                 colors = pal,
                 # note that type and mode are different for the 3D scatter
                 type = 'scatter3d', mode = 'markers',
                 # Hover text:
                 text = ~paste("Sample: ", sample, "<br>Condition:",
                               condition),
                 width = 800)

# title with variance explained
tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio[1:3])
tit = paste("Total Explained Variance = ", tot_explained_variance_ratio)

# changing camera angle so that PC1 and PC2 are as in the 2d plot
scene = list(camera = list(eye = list(x = -1, y = -1, z = 1)))

t_pca <- t_pca %>% layout(title = tit, legend = list(title = list(text = "Condition")), 
                          scene = scene)

d_pca
t_pca

# Explained variance plot.
# Notice the difference when using scale = TRUE in the amount of variance explained.
data.matrix <- filter_counts_df

prin_comp <- prcomp(t(data.matrix), scale.=TRUE)
#prin_comp <- prcomp(t(data.matrix))

exp_var_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
cumsum <- cumsum(exp_var_ratio[1:30]) # I'm limiting it to the first 30 components..
data <- data.frame(cumsum, seq(1, length(cumsum), 1))
colnames(data) <- c("Explained_Variance", "Components")

fig <- plot_ly(data = data, x = ~Components, y = ~Explained_Variance, type='scatter', mode='lines', fill='tozeroy') %>%
  layout(xaxis = list(title = "# Components", tickvals = seq(1, length(cumsum), 1)),
         yaxis = list(title = "Explained Variance"))

fig

## NORMALIZATION
# create a DGEList object
# this object has information about counts, samples and genes.
edge_c <- DGEList(counts=filter_counts_df, group=c_anno_df$condition, samples=c_anno_df, genes=filter_anno_df) 

# this does intra and between sample normalization keeping into account the dimensions and number of genes, amount of samples etc.
edge_n <- calcNormFactors(edge_c, method="TMM")

# checking the normalization 

# we can keep the normalized values into a new table with cpm values
# create a cpm-rpkm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2))

# checking the normalization factors we obtained in the last block
norm_factors <- mean(edge_n$samples$lib.size * edge_n$samples$norm.factors)/(edge_n$samples$lib.size*edge_n$samples$norm.factors)
names(norm_factors) <- edge_n$samples$sample

# check the same PCA plots as before but after normalization
data.matrix <- cpm_table
prin_comp <- prcomp(t(data.matrix), rank. = 2, scale. = TRUE)
components <- prin_comp[["x"]]              # get the components values
components <- data.frame(components)        # make a data.frame obj of them
components <- cbind(components, c_anno_df)  # add to the df the sample and condition information

# make a palette to have red -> Case,  blue -> Control
pal <- c("#E19F9D", "#9ECAE1")
pal <- setNames(pal, c("Case", "Control"))

d_pca <- plot_ly(components, x = ~PC1, y = ~PC2, color = ~condition,
                 colors = pal, type = 'scatter', mode = 'markers',
                 # Hover text:
                 text = ~paste("Sample: ", sample, "<br>Condition:",
                               condition),
                 width = 800)

# title with variance explained
tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio[1:2])
tit = paste("Total Explained Variance = ", tot_explained_variance_ratio)

d_pca <- d_pca %>% layout(title = tit, legend = list(title = list(text = "Condition")))

d_pca

# and the 3D version
data.matrix <- cpm_table
# rank. sets the maximum rank to compute, useful if we have way more
prin_comp <- prcomp(t(data.matrix), rank. = 3, scale. = TRUE)
components <- prin_comp[["x"]]              # get the components values
components <- data.frame(components)        # make a data.frame obj of them
components <- cbind(components, c_anno_df)  # add to the df sample and condition

# make a palette to have red -> Case,  blue -> Control
pal <- c("#E19F9D", "#9ECAE1")
pal <- setNames(pal, c("Case", "Control"))

t_pca <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition,
                 colors = pal,
                 # note that type and mode are different for the 3D scatter
                 type = 'scatter3d', mode = 'markers',
                 # Hover text:
                 text = ~paste("Sample: ", sample, "<br>Condition:",
                               condition),
                 width = 800)

# title with variance explained
tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio[1:3])
tit = paste("Total Explained Variance = ", tot_explained_variance_ratio)

# changing camera angle so that PC1 and PC2 are as in the 2d plot
scene = list(camera = list(eye = list(x = -1, y = -1, z = 1)))

t_pca <- t_pca %>% layout(title = tit, legend = list(title = list(text = "Condition")),
                          scene = scene)

t_pca

# I also want to print the old bargraph with the now normalized data.
fig1 <- plot_ly(x = colSums(raw_counts_df)/1000000,
                y = colnames(raw_counts_df),
                type = 'bar',
                name = "All genes",
                marker = list(color = '#9ECAE1', 
                              line = list(color = '#3A8EBC',
                                          width = 1.5)),
                width = 800, height = 500)

fig2 <- plot_ly(x = colSums(raw_counts_pc)/1000000,
                y = colnames(raw_counts_pc),
                type = 'bar',
                name = "Protein-coding genes",
                marker = list(color = '#E19F9D',
                              line = list(color = '#C74E4A',
                                          width = 1.5)),
                width = 800, height = 500)

fig3 <- plot_ly(x = colSums(filter_counts_df)/1000000,
                y = colnames(filter_counts_df),
                type = 'bar',
                name = "Filtered genes",
                marker = list(color = '#9DE1CC',
                              line = list(color = '#33A783',
                                          width = 1.5)),
                width = 800, height = 500)

fig4 <- plot_ly(x = colSums(cpm_table)/1000000,
                y = colnames(cpm_table),
                type = 'bar',
                name = "Normalized genes",
                marker = list(color = '#9DA6E0',
                              line = list(color = '#3B4CBB',
                                          width = 1.5)),
                width = 800, height = 500)

four_bgraph <- subplot(fig1, fig2, fig3, fig4, shareY = TRUE) %>% layout(title = "Millions of reads per sample")
four_bgraph

## Define the experimental design matrix
design <- model.matrix(~0+group, data=edge_c$samples)
colnames(design) <- levels(edge_c$samples$group)
rownames(design) <- edge_c$samples$sample
design

## Calculate dispersion and fit with edgeR (necessary for differential expression analysis)
# estimate the dispersion for each gene in the data
edge_d <- estimateDisp(edge_n, design)
edge_f <- glmQLFit(edge_d, design)

# Function to perform the differential analysis
edgeRglmQLF <- function(mat=edge_f,contro,cpm_mat=edge_n,label="",sig_thr=0.5,sig_col="CPM",fc_thr=0.5,pval_col="p_val",pval_thr=0.05,names=FALSE)
{
  degs <- glmQLFTest(edge_f,contrast=contro)$table[,-3]
  colnames(degs) <- c("log2_FC","log2_CPM","p_val")
  a_levels <- rownames(contro)[which(contro!=0)]
  a_samples <- which(cpm_mat$samples$group%in%a_levels)
  cpm_sele <- cpm(cpm_mat,log=T)[,a_samples]
  degs$log2_CPM <- apply(cpm_sele,1,function(x) mean(x))
  #degs<-exactTest(edge_c, pair=cond, dispersion=bcv^2)$table
  degs$p_adj <- p.adjust(degs$p_val, method ="BH")
  degs$class <- "="
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC>=fc_thr & degs[,pval_col]<=pval_thr),"class"] <- "+"
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC<=(-fc_thr) & degs[,pval_col]<=pval_thr),"class"] <- "-"
  degs$class <- as.factor(degs$class)
  degs$comp <- label
  degs$id <- rownames(degs)
  degs <- degs[,c("id","comp","log2_FC","log2_CPM","p_val","p_adj","class")]
  if(names=="TRUE"){
    newnames <- paste(label,colnames(degs),sep="_")
    colnames(degs) <- newnames
  }
  return(degs)
}

contro <- makeContrasts("Case-Control", levels=design)
DEGs <- edgeRglmQLF(mat=edge_f, cpm_mat=edge_n, contro=contro, label="CasevsControl", sig_thr=1, sig_col="log2_CPM", fc_thr=1.5, pval_thr=0.25, pval_col="p_adj",names=F)

# -- VOLCANO PLOT --
# to visualize divergency between expressions
input_df <- DEGs
xlabel <- "log2 FC Control vs Cases"
ylabel <- "-log10 adj_pvalue (FDR)"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$log2_FC, -log(input_df$p_adj, base=10), xlab=xlabel, ylab=ylabel, 
     col=ifelse(input_df$class=="=", "grey80", "grey50"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot")

# add lines for the used p-values
abline(v=0, lty=2, col="grey20")
abline(v=1.5, lty=2, col="grey20")
abline(v=-1.5, lty=2, col="grey20")
abline(h=-log10(0.25), lty=2, col="grey20")

# color up- and down- regulated
# red: D47774, blue: 74B2D4
ups = input_df[which(input_df$class=="+"),]
with(subset(ups, log2_FC > 1.5 & -log(p_adj, base=10) > 0.25), points(ups$log2_FC, -log(ups$p_adj, base=10), pch=20, cex=0.8, col="#D47774"))
dns = input_df[which(input_df$class=="-"),]
with(subset(dns, log2_FC < 1.5 & -log(p_adj, base=10) > 0.25), points(dns$log2_FC, -log(dns$p_adj, base=10), pch=20, cex=0.8, col="#74B2D4"))
rm(ups, dns)
# -- VOLCANO PLOT --

# -- HEATMAP --
# red: D47774, blue: 74B2D4
pal <- c("#74B2D4", "white", "#D47774") 
pal <- colorRampPalette(pal)(100)
clean_cpm_table <- cpm_table[which(rownames(cpm_table) %in% DEGs$id[which(DEGs$class!="=")]), ]
heatmap(as.matrix(clean_cpm_table), cexCol = 0.5, margins = c(4,4), col=pal, cexRow = 0.2)
# -- HEATMAP --


# Some cleanup commands
all_obj <- as.list(x = ls())
# keep c_anno, the filtered df and the DEGs
part_obj <- all_obj[which(all_obj %in% c("c_anno_df", "filter_anno_df", "filter_counts_df", "DEGs") == FALSE )]
rm(list = as.character(part_obj), part_obj)



# 4. Gene set enrichment analysis with clutserProfiler

# make a map of ids for the different packages
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
convert <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), filters=c("ensembl_gene_id"), values=DEGs$id, mart = ensembl)

# now we expand our table to contain this information
DEGs <- merge(DEGs, convert,by.x="id", by.y="ensembl_gene_id")

# Excluding all genes that don't have an entrezgene_id (from 17841 to 17421 genes!)
DEGs <- DEGs[which(!is.na(DEGs$entrezgene_id)), ]

# Removing entries that have the same entrezgene_id (duplicates) (17421 to 17374)
DEGs <- DEGs[-which(duplicated(DEGs$entrezgene_id)), ]

# then we can split again into the two up and downregulated lists
up_DEGs <- DEGs[which(DEGs$class=="+"), ]
down_DEGs <- DEGs[which(DEGs$class=="-"), ]

library("clusterProfiler")
library("org.Hs.eg.db")

### -- GO analysis --
## Biological Process
# up-regulated
ego_BP_up <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
# down-regulated
ego_BP_down <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

## Molecular Function
# up-regulated
ego_MF_up <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
# down-regulated
ego_MF_down <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

## plotting the results
## Biological Process
# upregulated
#barplot(ego_BP_up, showCategory=15)
dotplot(ego_BP_up, showCategory=10, title = "Biological Process up-regulated")
# downregulated
#barplot(ego_BP_down, showCategory=15)
dotplot(ego_BP_down, showCategory=10, title = "Biological Process down-regulated")
## Molecular Function
# upregulated
#barplot(ego_MF_up, showCategory=15)
dotplot(ego_MF_up, showCategory=10, title = "Molecular Function up-regulated")
# downregulated
#barplot(ego_MF_down, showCategory=15)
dotplot(ego_MF_down, showCategory=10, title = "Molecular Function down-regulated")


## -- KEGG enrichment analysis --
# this function automatically gets the needed KEGG pathways from the database
# upregulated
ekegg_up <- enrichKEGG(gene = up_DEGs$entrezgene_id,
                       organism = 'human',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
# downregulated
ekegg_down <- enrichKEGG(gene = down_DEGs$entrezgene_id,
                         organism = 'human',
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

# but we can look at them just by plotting the top 10 pathways, like for the GO analysis
#barplot(ekegg_up, showCategories = 10)
dotplot(ekegg_up, showCategory=10, title = "KEGG pathways up-regulated")
dotplot(ekegg_down, showCategory=10, title = "KEGG pathways down-regulated")



# 5. Pathway visualization
library("pathview")

## Visualizing the Alcoholism pathways

# First I need to find the relative id of the pathway
ekegg_up@result[["Description"]][1:10]
ekegg_up@result[["ID"]][6]  # --> "hsa05034"

# Then I use the up-regulated genes
logFC <- up_DEGs$log2_FC
names(logFC) <- up_DEGs$entrezgene_id

# and the pathway.id I found, into this function that visualizes the graphical pathway
pathview(gene.data = logFC, 
         pathway.id = "hsa05034", 
         species = "human")



# 6. Enriched Transcription Factors in the promoters of DE genes
library("MotifDb")
library("seqLogo")
library("PWMEnrich")
library("PWMEnrich.Hsapiens.background")

# connect to BioMart
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# get promoter sequences 500nt upstream
up_promoters <- getSequence(id = up_DEGs$id, 
                            type="ensembl_gene_id",
                            seqType="gene_flank",
                            upstream=500,
                            mart=ensembl)

data(PWMLogn.hg19.MotifDb.Hsap)
sequences <- sapply(up_promoters$gene_flank, DNAString)

#enriched_TFs <- motifEnrichment(sequences, PWMLogn.hg19.MotifDb.Hsap, score = "affinity")
# save the TFs because it took a lot to make them
#save(enriched_TFs, file="eTF.RData")

# so we can directly load them now
load("eTF.RData")

report <- groupReport(enriched_TFs)

plot(report[1:5], fontsize=15, id.fontsize=15)



# 7. Computation of Empirical distributions of scores for a TF
# NRL gene
report$id[2]

# to get the seqLogo of this TF
# seqLogo(apply(report@pwms[[2]]@pfm,2,function(x) x/sum(x)))

nrl_tfs_motifs <- subset(MotifDb, organism=='Hsapiens' & geneSymbol==report$target[2])
nrl_PWMs <- toPWM(as.list(nrl_tfs_motifs))
ecdf <- motifEcdf(nrl_PWMs, organism = "hg19", quick=TRUE)
thresholds_995 <- lapply(ecdf, function(x) log2(quantile(x, 0.995)))



# 8. Genes with over-threshold binding scores
scores_995 <- motifScores(sequences, nrl_PWMs, raw.score=FALSE, cutoff = unlist(thresholds_995))
#ecdf <- motifEcdf(nrl_PWMs, organism = "hg19", quick=TRUE)

# rows where the sum over the row is bigger than 0
# if it was 0 --> not over threshold
over_threshold_995 <- scores_995[which(apply(scores_995, 1, sum) > 0), ]

# fraction of the upregulated
fraction_995 <- c(length(which(apply(scores_995, 1, sum)>0))/808)
fraction_995

# same analysis at 0.999 threshold

# we already have the old ecfd, only compute new thresholds and calculate the scores
thresholds_999 <- lapply(ecdf, function(x) log2(quantile(x, 0.999)))
scores_999 <- motifScores(sequences, nrl_PWMs, raw.score=FALSE, cutoff = unlist(thresholds_999))

over_threshold_999 <- scores_999[which(apply(scores_999, 1, sum) > 0), ]
fraction_999 <- c(length(which(apply(scores_999, 1, sum)>0))/808)
fraction_999 # -> 0.9715 (785 genes)

# finding names of the 0.999 genes
names_999 <- up_promoters$ensembl_gene_id[which(up_promoters$gene_flank %in% rownames(over_threshold_999))]
# here we get the ensembl ids, but we can get the external names:
ext_names_999 <- up_DEGs$external_gene_name[which(up_DEGs$id %in% names_999)]



# 9. List of up-regulated gene ids to do an analysis on STRING
library("STRINGdb")

# I will use R to download the interaction .tsv file
# make a list of all the ids
#list_id_upreg <- up_DEGs$id

# make an obj to connect to the STRING database 
#string_db <- STRINGdb$new(species = 9606, score_threshold = 400, version = '11')

# map the upDEGs to the STRING_id
#up_DEGs_mapped <- string_db$map(up_DEGs, 'id', removeUnmappedRows = TRUE)
# luckily there were no unmapped rows

# we can then obtain all the interactions
#up_DEGs_interactions <- string_db$get_interactions(up_DEGs_mapped$STRING_id)

# the code above worked, but I couldn't get the same results as the list downloaded
# manually from the STRINGdb webpage.. I guess the tool is just old..
# so I wrote the upDEGS on a file and used that to compute the graph with the webapp
# and download the interactions to analyze them with igraph
#write.table(up_DEGs_mapped$STRING_id, file="list_id_DEGs.txt", row.names=F, col.names=T, sep="\t", quote=F)



# 10. Create network from STRING PPI
library("igraph")

links <- read.delim("string_interactions.tsv")

# Create nodes annotations using biomaRt
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
nodes <- getBM(attributes=c("external_gene_name","ensembl_gene_id","description","gene_biotype","start_position","end_position","chromosome_name","strand"),
               filters=c("ensembl_gene_id"), 
               values=up_DEGs$id,
               mart = ensembl)
# dropping ensembl_gene_id, chromosome_name, and strand
nodes <- unique(nodes[,c(1,3:6)])

# Some of the gene products in the links are not present in the nodes list, and must therefore be removed
symbols <- union(links$X.node1,links$node2)
links <- links[which(links$X.node1%in%intersect(nodes$external_gene_name,symbols)),]
links <- links[which(links$node2%in%intersect(nodes$external_gene_name,symbols)),]

## Create the network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=FALSE) 

## Plot the PPI network
# for the report I actually decided to use the one from STRING directly
plot(net, 
     edge.width=3,
     vertex.color="orange",
     vertex.size=5,
     vertex.frame.color="darkgray",
     vertex.label.color="black", 
     vertex.label.cex=0.2,
     edge.curved=0.1) 

# Analysis of the network
# Average diameter
mean_distance(net, directed=F)

# Clustering coefficient
transitivity(net, type="global")

# Find largest connected component
comp <- components(net)

# here we can see that we have
comp$no # 131 clusters

# the largest connected cluster is of 672 genes, and it is the cluster '1'
max(comp$csize)

# we can then check which genes are in the cluster '1'
names(comp$membership[which(comp$membership == 1)])

# the most connected gene can be found by looking at the degree of the nodes
res <- degree(net)
print(res[which(res == max(res))]) # -> CDK1 with degree 370


