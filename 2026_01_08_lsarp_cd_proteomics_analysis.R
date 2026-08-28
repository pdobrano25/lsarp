### 2024_10_21  Host Proteomics + MLI ASV Processing + Analysis

# Goals (not in this order)
# 1: process proteomics data feature table for LSARP-CD analysis
# 2: process MLI ASV data

# 3: Combined proteomics analysis (PCA, complex LM + enrichment)
#.     RS and Placebo; Responders and Non-Responders (RS and Placebo)
# 4: Combined MLI analysis (PCA, complex LM + enrichment)
#.     RS and Placebo; Responders and Non-Responders (RS and Placebo)

# Sensitivity
# : Combined proteomics with GO instead of MitoCarta
# : Enrichment, p value thresholds
# : Proteomics per location
# : MLI per location
# : Reduced proteomics combined
# : Reduced proteomics per location
# : Reduced MLI combined
# : Reduced MLI per location



# :: save -----------------------------------------------------------------

v.date = "2026_03_04"
# save.image("./2026_03_04_lsarp_big_save.Renv")

# load("./2026_03_04_lsarp_big_save.Renv")

# :: load packages --------------------------------------------------------

# TOGGLE TO RE-RUN LENGTHY ANALYSES
rerun=T

library("dplyr"); library("ggplot2")

setwd("~/Documents/PhD/git_oars_archfolder")

# load processed R objects
# load(file = "./2025_12_02_lsarp_16s_data_meta.Renv")

# :: functions ------------------------------------------------------------


# function to make volcano plots have an intelligible y axis
neg_log10_trans <- scales::new_transform(
  name = "neglog10",
  transform = function(x) -log10(x),
  inverse = function(x) 10^(-x),
  format = function(x) format(x, scientific = FALSE)
)

# :: load metadata ------------------------------------------------------------

proteomics.metadata = read.csv("./host_proteomics/plate_design.csv", sep=",")

proteomics.metadata$Time %>% unique()
proteomics.metadata$Location %>% unique()
proteomics.metadata$Status %>% unique()
# fix typo
colnames(proteomics.metadata)[colnames(proteomics.metadata)=="Stasus"] = "Status"

## create interpretable code for samples based on HM-Time-Location-Status
proteomics.metadata$HM = ifelse(nchar(proteomics.metadata$SampleID) == 3, paste("HM0",proteomics.metadata$SampleID, sep=""),
                                paste("HM", proteomics.metadata$SampleID, sep=""))
proteomics.metadata$code = paste(proteomics.metadata$HM,
                                 proteomics.metadata$Time,
                                 proteomics.metadata$Location,
                                 proteomics.metadata$Status,
                                 sep="_")
## create code matching feature table columns
proteomics.metadata$sample = paste(
  paste("Plate", proteomics.metadata$Plate, sep=""),
  proteomics.metadata$Well, sep="_")

## add clinical data

# unblinding: treatment group, clinical data, flares, compliance
unblinding = read.csv("./2024_07_31_unblinding_DM.csv")
unblinding$group = ifelse(unblinding$group == "(RS)", "RS",
                          ifelse(unblinding$group == "(Plac)", "Plac",
                                 unblinding$group))
unblinding$Group = ifelse(unblinding$group == "RS", "RS",
                          ifelse(unblinding$group == "Plac", "Placebo", NA))

# merge
proteomics.metadata$study_id = paste(proteomics.metadata$HM, ".00", sep="")
proteomics.metadata = merge(proteomics.metadata,
                            unblinding[,c("study_id", "Group", "flare.call_verified", "days.taken", "non.compliant")],
                            by="study_id")
unique(proteomics.metadata$flare.call_verified)
# n = 40 total

unique(subset(proteomics.metadata, flare.call_verified == "completed (7)")$HM) %>% length()
# 22 completed study
unique(subset(proteomics.metadata, flare.call_verified == "flare (8)")$HM) %>% length()

# subset to only completed / flare (i.e. no non-compliant / excluded)

proteomics.metadata = subset(proteomics.metadata, 
                             flare.call_verified %in% c("completed (7)", "flare (8)"))

table(distinct(proteomics.metadata[,c("HM", "Group")])$Group)
# 12 Placebo + 18 RS = n 30 total

table(distinct(proteomics.metadata[,c("HM", "Group")]))

# Exclude additional participants based on compliance (keep ~70% placebo)
proteomics.metadata = subset(proteomics.metadata, !HM %in% c("HM0899", "HM0966", "HM0970", "HM0978", "HM1030", "HM1039"))
proteomics.metadata = subset(proteomics.metadata, !HM %in% c("HM1045")) # note: 1045 (appendicitis) wasn't even there

unique(proteomics.metadata$HM)
unique(proteomics.metadata$HM) %>% length()

# n = 24 (good; missing n=2 (HM0978 and HM0998)
distinct(proteomics.metadata[,c("Group", "HM")])[,c("Group")] %>% table()
# n = 11 placebo + 13 RS

# add response from other scripts
proteomics.metadata = merge(proteomics.metadata,
      lsarp.metadata.responders[,c("HM", "flare.group")]%>% distinct(),
      by="HM")

# :: visualize samples ----------------------------------------------------

ggplot(proteomics.metadata %>% subset(Time != 2) %>%
         mutate(Status = factor(Status, levels=c("N", "A")))%>%
         arrange((Status)),
       aes(x=as.factor(Time), y=HM))+
  ggbeeswarm::geom_beeswarm(aes(shape=Status, color=flare.group), 
                            alpha=0.5, size=2.5, cex=5)+
  scale_color_manual(values=c("red", "blue"))+
  #scale_shape_manual(values=c(24,25))+
  scale_shape_manual(values=c(17,16))+
  
  facet_grid(Group~Location, scales="free")+
  theme_classic()+
  labs(x="Scope", y="")



# :: OPTIONAL: subset samples -------------------------------------------------------

# for sensitivity analyses:
proteomics.subset = proteomics.metadata %>% 
  # remove time 2 (12M scope)
  subset(Time != 2) %>%
  # select most inflamed per patient/timepoint/location
  mutate(Status = factor(Status, levels=c("N", "A")))%>%
  group_by(HM, Time, Location) %>%
  #mutate(status.num = as.numeric(Status)) %>% as.data.frame() %>%
  #dplyr::select(HM, Time, Location, Status, status.num)
  slice_max(order_by=as.numeric(Status), n=1)

ggplot(proteomics.subset %>% subset(Time != 2) %>%
         arrange((Status)),
       aes(x=as.factor(Time), y=HM))+
  ggbeeswarm::geom_beeswarm(aes(shape=Status, color=flare.group), 
                            alpha=0.5, size=2.5, cex=5)+
  scale_color_manual(values=c("red", "blue"))+
  #scale_shape_manual(values=c(24,25))+
  scale_shape_manual(values=c(17,16))+
  
  facet_grid(Group~Location, scales="free")+
  theme_classic()+
  labs(x="Scope", y="")

# Use these for sensitivity analyses, if necessary

proteomics.metadata # full data
proteomics.subset # subset data

# :: load data ------------------------------------------------------------

proteomics.data.raw = read.csv("./host_proteomics/biopsy_report_proteinGroups_allsamples.csv", sep=",")
colnames(proteomics.data.raw) # no QCs

# replace names on feature table
colnames(proteomics.data.raw)[-c(1:5)]  = proteomics.metadata[match(colnames(proteomics.data.raw)[-c(1:5)] , proteomics.metadata$sample),]$code
# remove samples that are not from patients
proteomics.data.raw = proteomics.data.raw[,!is.na(colnames(proteomics.data.raw))]

# good

# check annotations
table(proteomics.data.raw$Genes) %>% data.frame() %>% arrange(-Freq)

length(unique(proteomics.data.raw$Genes))
# 7279 genes

# :: process data ---------------------------------------------------------

# replace unannotated Genes with "unannotated"
proteomics.data.raw$Genes = ifelse(proteomics.data.raw$Genes == "", "unannotated", proteomics.data.raw$Genes)

# sum up protein intensities of the same Protein.Name name
proteomics.data.df = reshape2::melt(proteomics.data.raw[,-c(1:2,4:5)])
proteomics.data.df = proteomics.data.df %>%
  # use Protein.Names, not Genes (and then remove "_HUMAN" tag. Otherwise, Excel transforms to date)
  # Note, Protein names do not always = Gene name.
  mutate(Protein.Names = gsub("_HUMAN", "", Protein.Names)) %>%
  group_by(Protein.Names, variable) %>%
  mutate(value = sum(value)) %>% 
  dplyr::select(Protein.Names, variable, value) %>% distinct() %>% data.frame()
proteomics.data = reshape2::acast(proteomics.data.df, variable ~ Protein.Names, value.var="value")

# replace NA with 0
proteomics.data[is.na(proteomics.data)] <- 0

# calculate universal pseudocount
proteomics.pseudocount = min(proteomics.data[proteomics.data!=0])/2

# log and add pseudocount
proteomics.data.processed = log2(proteomics.data + proteomics.pseudocount)
dim(proteomics.data.processed)

# save Gene=Protein map
proteomics.gene.map = proteomics.data.raw[,c("Protein.Names", "Genes")] %>% distinct()

colnames(proteomics.gene.map) = c("protein", "gene")
# remove "_HUMAN"
proteomics.gene.map$protein = gsub("_HUMAN", "", proteomics.gene.map$protein)
# if the gene is a date, convert to protein name

proteomics.gene.map$gene = ifelse(grepl(paste(c("-Jan", "-Feb", "-Mar", "-Apr", "-May", "-Jun", "-Jul", "-Aug", "-Sep", "-Oct", "-Nov", "-Dec"), collapse="|"),proteomics.gene.map$gene),
                                        proteomics.gene.map$protein, proteomics.gene.map$gene)
# add make.names version of protein for merging after data.frame()
proteomics.gene.map$protein.x = make.names(proteomics.gene.map$protein)

proteomics.gene.map %>% arrange(protein)

# :: analyze --------------------------------------------------------------

# reduce to samples to keep (based on full data or subset)
# note: samples have already been pruned to those non-compliant / early flared / appendicitis
proteomics.data.raw = proteomics.data[rownames(proteomics.data) %in% proteomics.metadata$code,]

# filter proteins by 80% prevalence + abundance
dim(proteomics.data.raw) # 7,484 proteins
proteins.to.keep = proteomics.data.raw
proteins.to.keep[proteins.to.keep != 0] <- 1
proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
length(proteins.to.keep) # 4,813 proteins present > 80% of samples

proteomics.data = proteomics.data.raw[,colnames(proteomics.data.raw) %in% proteins.to.keep]
dim(proteomics.data) # 4,813 proteins, 177 samples

# ensure only T0 and T1 are kept
proteomics.data = proteomics.data[rownames(proteomics.data) %in% subset(proteomics.metadata, Time %in% c(0,1))$code,]
nrow(proteomics.data)

# :: mitocarta ------------------------------------------------------------

mitocarta.proteins = read.csv("./host_proteomics/mitocarta3_proteins.csv")
mitocarta.pathways = read.csv("./host_proteomics/mitocarta3_pathways.csv")

# first, unpack all protein names (and synonyms) identified as mito-proteins
mito.proteins = c(mitocarta.proteins$Symbol, unlist(strsplit(mitocarta.proteins$Synonyms, split = "\\|")))
# these are actually gene names
nrow(mitocarta.proteins) # 1157 mitochondrial genes
length(mito.proteins) # with 4726 total potential gene names

# assign mito-pathway per gene
mitocarta.pathways.expanded <- mitocarta.pathways %>%
  tidyr::separate_rows(Genes, sep = ",") %>% dplyr::select(-X2)
# fix colname to match other data
colnames(mitocarta.pathways.expanded)[colnames(mitocarta.pathways.expanded) == "Genes"] = "gene"

# delete extraneous spaces
mitocarta.pathways.expanded$gene = gsub(" ", "", mitocarta.pathways.expanded$gene)
# remember to use this to link mito to proteins:
proteomics.gene.map

# :: GO analysis ----------------------------------------------------------

proteomics.gene.map

# BiocManager::install("biomaRt")

ensembl <- biomaRt::useMart("ensembl", 
                            dataset = "hsapiens_gene_ensembl",
                            verbose = T)
# version 0.7

go_annotations <- biomaRt::getBM(
  attributes = c("hgnc_symbol",          # Gene name (input)
                 "go_id",                # GO term ID
                 "name_1006",            # GO term name (description)
                 "namespace_1003",       # GO namespace (BP, MF, CC)
                 "go_linkage_type"),     # Evidence code (e.g., IEA for inferred)
  filters = "hgnc_symbol",              # Filter by HGNC gene symbols
  values = proteomics.gene.map$gene,
  mart = ensembl
)

# clean up
go_annotations

saveRDS(go_annotations, "./2026_03_03_go_annotations.Rds")

length(unique(go_annotations$hgnc_symbol)) # 6849 annotated, expressed genes

# clean up (only keep GO term and gene; unique)
go_annotations = go_annotations %>% 
  dplyr::select(hgnc_symbol, name_1006, namespace_1003) %>% 
  subset(name_1006 != "") %>% distinct()
length(unique(go_annotations$hgnc_symbol)) # 6821 kept

# keep only biological process (default in GO enrichment)
go_annotations_bio = subset(go_annotations, namespace_1003 == "biological_process")
length(unique(go_annotations_bio$hgnc_symbol)) # 6462 kept

go_mitochondria_annotations = subset(go_annotations, grepl("itoch", name_1006))
unique(go_mitochondria_annotations$hgnc_symbol) %>% length() # 1153 flagged as mitochondrial (cellular component OR biological process)
unique(go_mitochondria_annotations$name_1006)

sum(unique(go_mitochondria_annotations$hgnc_symbol) %in% unique(mitocarta.pathways.expanded$gene)) %>% sum() 
# 666 genes overlap

# create final mito map

# add mitochondrial labels
proteomics.mito.map = merge(proteomics.gene.map,
                            mitocarta.pathways.expanded,
                            by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# note that feature names (in feature table) are proteins, not genes

# save these, because biomart is unreliable (transitioning servers?)

# >>> 1. MPX ------------------------------------------------------------------

# Note: here, MPX refers to "Mucosal Proteomics", not metaproteomics

# :: Sample selection -----------------------------------------------------

# CHECK SAMPLES
proteomics.metadata$HM %>% unique()
proteomics.subset$HM %>% unique()

# write.csv(proteomics.metadata[,c("SampleID", "Time")] %>% distinct(), "./host_proteomics/2025_09_10_host_biopsy_timing.csv")
# NOTE: I added dates from redcap
biopsy.timings = read.csv("./host_proteomics/2025_09_10_host_biopsy_timing.csv")
colnames(biopsy.timings) = c("x", "HM", "Time", "scope", "rs", "x1")
biopsy.timings = biopsy.timings %>% dplyr::select(-x, -x1)
biopsy.timings = biopsy.timings %>%
  mutate(diff = as.Date(scope,"%d_%m_%Y")-as.Date(rs, "%d-%m-%Y"))

# HM0878 does not have 4-5M scope
# others look fine
# HM0960 had a scope, but no biopsies

# but, it seems like some had scopes during the RS treatment period?
# let's visualize

biopsy.timings = biopsy.timings %>%
  mutate(HM = ifelse(nchar(HM) == 3, paste("HM0", HM, sep=""), paste("HM", HM, sep="")))
biopsy.timings = reshape2::acast(biopsy.timings, HM ~ Time, value.var="scope") %>% as.data.frame()
colnames(biopsy.timings) = c("scope_0", "scope_1")
biopsy.timings$HM = rownames(biopsy.timings)

biopsy.timings = merge(biopsy.timings,
                       metadata.lsarp.stool, by="HM", all.y=T)
biopsy.timings = biopsy.timings %>%
  mutate(scope_0_days = as.Date(scope_0,"%d_%m_%Y") - as.Date(rs_start_date,"%Y-%m-%d"),
         scope_1_days = as.Date(scope_1, "%d_%m_%Y") - as.Date(rs_start_date,"%Y-%m-%d"))

# keep only samples to be analyzed
biopsy.timings = subset(biopsy.timings, HM %in% unique(proteomics.metadata$HM))

# visualize stools + timings
ggplot(biopsy.timings %>% arrange(lsarp.days),
                                   aes(x=lsarp.days, y=reorder(HM, (exit.day))))+
  annotate("rect", xmin=0, xmax=max(subset(biopsy.timings, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  # add background point to look neater
  geom_point(fill = "white", shape=21, size=5)+
  # then plot real points (they overlap in the window, but look good when saved at full size)
  geom_point(aes(fill=lsarp.on.rs, 
                 alpha = ifelse(lsarp.days <= exit.day, "1", "0")),
             color="white", 
             shape=21, size=5)+
  # add exit reason
  geom_point(data=biopsy.timings[,c("exit.day","exit.reason", "lsarp.days","HM", "Group")] %>% distinct()%>% arrange(lsarp.days),
             aes(x=exit.day, y=reorder(HM, (lsarp.days)), shape=exit.reason),  size=4)+
  # add scope dates
  geom_point(data=biopsy.timings[,c("HM","Group", "scope_0_days", "lsarp.days")] %>% group_by(HM) %>% slice_max(order_by=lsarp.days, n =1)%>% arrange(lsarp.days),
             aes(x=scope_0_days, y=reorder(HM, (lsarp.days))),  shape=1, size=3)+
  geom_point(data=biopsy.timings[,c("HM","Group", "scope_1_days", "lsarp.days")] %>% group_by(HM) %>% slice_max(order_by=lsarp.days, n =1) %>%arrange(lsarp.days),
             aes(x=scope_1_days, y=reorder(HM, (lsarp.days))), shape=1, size=3)+
  scale_fill_manual(values=c(2,"grey","grey"))+
  scale_alpha_manual(values=c(0.25,1))+
  scale_shape_manual(values=c(4))+
  # geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), color="white", nudge_y=-0.1, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-1, nchar(standard.name))), size=2)+
  theme_minimal()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~Group,nrow=1, scales="free")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        panel.grid.major.y=element_line(),
                        strip.background = element_rect(
                          color="black"))
# perfect: diamonds = scope + proteomic data
# Missing Baseline: HM0978, 998
# Missing Re-scope: HM0960, 1051, 978, 906, 998
# + = "excluded"
# x = "flare"

table(proteomics.metadata[,c("HM", "Time")])

# :: Outlier Identification -----------------------------------------------

# Note: RPCA for outlier detection has been used here:  DOI: 10.1186/s12859-020-03608-0 
# Unfortunately, their "outlier map" method is opaque
# So we'll use 95% confidence interval

dim(proteomics.data)

# remove 0 var proteins
proteomics.data = proteomics.data[,apply(proteomics.data, 2, sd)>0]
# remove non-patient samples // already done
proteomics.data = proteomics.data[grepl("_", rownames(proteomics.data)),]
# run PCA
proteomics.pca = rrcov::PcaHubert(log2(proteomics.data+(min(proteomics.data[proteomics.data!=0])/2)), 
                                  scale=T, center=T)
# extract data
proteomics.pca.df = proteomics.pca$scores[,c(1:4)] %>% data.frame()
# extract var explained
#proteomics.pca.var <- ((proteomics.pca$sdev^2) / sum(proteomics.pca$sdev^2))[1:2]
# merge with meta
proteomics.pca.df$code = rownames(proteomics.data)
proteomics.pca.df = merge(proteomics.pca.df,
                          proteomics.metadata, by="code")

# identify samples within confidence ellipse (loop through locations)
proteomics.pca.outliers = do.call(rbind, lapply(c("DC", "PC", "TI"), function(x){
  data.subset = subset(proteomics.pca.df, Location == x)
  confidence_level = 0.95
  # on PC1 and PC2
  ellipse = car::dataEllipse(data.subset[,c(2)], data.subset[,c(3)], 
                              levels = confidence_level, draw = FALSE)
  
  # Extract ellipse boundaries
  mahal_dist = mahalanobis(data.subset[,c(2:3)],
                            center = colMeans(data.subset[,c(2:3)]), 
                            cov = cov(data.subset[,c(2:3)]))

  # Determine critical value for % confidence (chi-squared distribution)
  critical_value = qchisq(confidence_level, df = 2)  # df = 2 for PC1 and PC2
  
  # Identify outliers (samples with Mahalanobis distance > critical value)
  outliers = which(mahal_dist > critical_value)
  outlier_samples = data.subset$code[outliers]
  # if no outliers, input NA
  if(length(outlier_samples) == 0){
    outlier_samples = NA
  }
  # output dataframe
  data.frame(outliers = outlier_samples,
             location = x)
}))
proteomics.pca.outliers

proteomics.pca.df$outlier = ifelse(proteomics.pca.df$code %in% proteomics.pca.outliers$outliers, "outlier", NA)

# use this one to show outliers
ggplot(proteomics.pca.df,
       aes(x=PC1, y=PC2))+
  #geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(aes(shape = Location, fill=Status))+
  #stat_ellipse(level = 0.95, linetype=2, type="t")+
  ggrepel::geom_text_repel(aes(label = ifelse(outlier == "outlier", code, NA)),size=2)+
  scale_shape_manual(values=c(21,22,23))+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~Location)

table(proteomics.pca.df$outlier)
# 9 outliers

# :: PCA -------------------------------------------------------

# 
proteomics.data.clean = proteomics.data[rownames(proteomics.data) %in% subset(proteomics.pca.df, is.na(outlier))$code,]
dim(proteomics.data.clean) # 168 samples remain

# check if calprotectin proteins correlate
plot(log10(as.data.frame(proteomics.data.clean)$S10A8),
     log10(as.data.frame(proteomics.data.clean)$S10A9))
# good

# any Excel funny business?
colnames(proteomics.data.clean)[grepl("Sep", colnames(proteomics.data.clean))] %>% dim()
colnames(proteomics.data.clean)[grepl("SEPT", colnames(proteomics.data.clean))]

proteomics.clean.pca = prcomp(log2(proteomics.data.clean+(min(proteomics.data.clean[proteomics.data.clean!=0])/2)))
# extract data
proteomics.clean.pca.df = proteomics.clean.pca$x[,c(1:4)] %>% data.frame()
# extract var explained
proteomics.clean.pca.var <- ((proteomics.clean.pca$sdev^2) / sum(proteomics.clean.pca$sdev^2))[1:2]
# merge with meta
proteomics.clean.pca.df$code = rownames(proteomics.clean.pca.df)
proteomics.clean.pca.df = merge(proteomics.clean.pca.df,
                                proteomics.metadata, by="code")

# check status!
proteomics.clean.pca.plot = ggplot(proteomics.clean.pca.df %>%
                                     mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
       aes(x=PC1, y=PC2))+
  #geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(aes(shape = Location, fill=Status))+
  stat_ellipse(aes(color=Status), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = code),size=2)+
  scale_shape_manual(values=c(21,22,23), guide=F)+
  theme_classic()+theme(#legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~Location, scales="free")+
  labs(x=paste("PC1:", round(proteomics.clean.pca.var, digits=4)[1]*100, "%"),
       y=paste("PC2:", round(proteomics.clean.pca.var, digits=4)[2]*100, "%"))
proteomics.clean.pca.plot

# numbers
proteomics.clean.pca.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()
proteomics.clean.pca.df[,c("HM", "Group")] %>% table()
proteomics.clean.pca.df[,c("HM", "Group", "Location", "Time")] %>% table()

# PERMANOVA
set.seed(25)
proteomics.clean.pca.permanova = vegan::adonis2(dist(proteomics.clean.pca.df[,c("PC1", "PC2")]) ~
                  Location + Status + Plate + Time*Group, 
               by="margin",
               # control for repeat measures of HM
               strata=proteomics.clean.pca.df$HM,
               # data
                proteomics.clean.pca.df %>% mutate(Plate = as.factor(Plate)))
proteomics.clean.pca.permanova = as.data.frame(proteomics.clean.pca.permanova)[c(1,2,4), c(3,5)] %>%
  mutate(Variable = rownames(.)) %>%
  mutate(R2 = round(R2, digits=3))
proteomics.clean.pca.permanova = proteomics.clean.pca.permanova[,c(3,1,2)] %>%
  reshape2::melt()
proteomics.clean.pca.permanova.plot = ggplot(proteomics.clean.pca.permanova %>% 
         mutate(Variable = gsub("Time:Group", "Time:Group", Variable))%>%
         mutate(variable = factor(variable, levels=c("R2", "Pr(>F)")),
                Variable = factor(Variable, levels=c("Location", "Status", "Time:Group"))),
       aes(x=variable, y=Variable))+
  geom_tile(fill="white", color="black")+
  geom_text(aes(label=value, fontface = ifelse(variable == "Pr(>F)" & `value` < 0.05, "bold", "plain")))+
  scale_x_discrete(position = "top") +
  theme_minimal()+theme(axis.text = element_text(size=12),
                        panel.grid.major=element_blank())+
  labs(x="", y="")
proteomics.clean.pca.permanova.plot

# make a plot, not a table
proteomics.clean.pca.permanova.r2.plot = proteomics.clean.pca.permanova %>%
  reshape2::acast(Variable ~ variable, value.var="value") %>% data.frame() %>%
  mutate(`𝘱 value` = `Pr..F.`) %>%
  mutate(FDR = p.adjust(`𝘱 value`, method="bonferroni")) %>%
  mutate(sig = ifelse(`FDR` < 0.05, "***",
                      ifelse(`FDR` < 0.20, "*",
                      ifelse(`𝘱 value` < 0.05, "+", NA)))) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = factor(Variable, levels=c("Time:Group", "Status", "Location"))) %>%
  ggplot(aes(x=R2, y=Variable))+
  geom_bar(stat="identity", color="white", fill="black", width=0.5)+
 # scale_color_manual(values=c("*" = "grey", "***" = "black"))+
  geom_text(aes(x=R2+0.015, label=sig, vjust=ifelse(sig == "+", 1, 0.8)), size=6)+
  scale_x_continuous(limits=c(0,0.2))+
  theme_classic()+
  labs(x=expression(R^2), y=NULL)+
  theme(panel.grid.major.x = element_line(color="grey", linewidth=0.2, linetype=2))

# "When controlling for inflammation (and region),
# RS has no impact on host proteome"
# (but it's a small effect size)

# plot groups

proteomics.clean.pca.plot = ggplot(proteomics.clean.pca.df %>%
                                          #subset(Group == "RS") %>%
                                          mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
                                        aes(x=PC1, y=PC2))+
  # geom_path(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(shape = Location, fill=flare.group), color="white", size=2.5)+
  stat_ellipse(aes(color=Group), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = code),size=2)+
  scale_shape_manual(values=c(21,22,23))+
  guides(shape="none", fill="none")+
  theme_classic()+theme(#legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))+
  facet_wrap(~Location, scales="free")+
  labs(x=paste("PC1:", round(proteomics.clean.pca.var, digits=4)[1]*100, "%"),
       y=paste("PC2:", round(proteomics.clean.pca.var, digits=4)[2]*100, "%"),
       fill="Group", color="Group")
proteomics.clean.pca.plot



# :: PCA Responders  --------------------------------------------


# numbers
proteomics.clean.pca.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% distinct() %>% dplyr::select(flare.group) %>% table()
proteomics.clean.pca.df[,c("HM", "Group")] %>%  subset(Group == "RS") %>% table()
proteomics.clean.pca.df[,c("HM", "Group", "Location", "Time", "flare.group")] %>%  subset(Group == "RS") %>% table()

# 
set.seed(25)
proteomics.clean.pca.resp.permanova = vegan::adonis2(dist(subset(proteomics.clean.pca.df, Group == "RS")[,c("PC1", "PC2")]) ~
                 Location + Status + Plate + Time*flare.group, 
               by="margin",
               # control for repeat measures
               strata=as.factor(subset(proteomics.clean.pca.df, Group == "RS")$HM),
               # data
               subset(proteomics.clean.pca.df, Group == "RS") %>%
                 mutate(Plate = as.factor(Plate)))
proteomics.clean.pca.resp.permanova = as.data.frame(proteomics.clean.pca.resp.permanova)[c(1,2,4), c(3,5)] %>%
  mutate(Variable = rownames(.)) %>%
  mutate(R2 = round(R2, digits=3))
proteomics.clean.pca.resp.permanova = proteomics.clean.pca.resp.permanova[,c(3,1,2)] %>%
  reshape2::melt()
proteomics.clean.pca.permanova.resp.plot = ggplot(proteomics.clean.pca.resp.permanova %>% 
                                               mutate(Variable = gsub("Time:flare.group", "Time:Response", Variable))%>%
                                               mutate(variable = factor(variable, levels=c("R2", "Pr(>F)")),
                                                      Variable = factor(Variable, levels=c("Location", "Status", "Time:Response"))),
                                             aes(x=variable, y=Variable))+
  geom_tile(fill="white", color="black")+
  geom_text(aes(label=value, fontface = ifelse(variable == "Pr(>F)" & `value` < 0.05, "bold", "plain")))+
  scale_x_discrete(position = "top") +
  theme_minimal()+theme(axis.text = element_text(size=12),
                        panel.grid.major=element_blank())+
  labs(x="", y="")
proteomics.clean.pca.permanova.resp.plot

# make plot, not table
proteomics.clean.pca.permanova.resp.r2.plot = proteomics.clean.pca.resp.permanova %>%
  mutate(Variable = gsub("Time:flare.group", "Time:Response", Variable))%>%
  reshape2::acast(Variable ~ variable, value.var="value") %>% data.frame() %>%
  mutate(`𝘱 value` = `Pr..F.`) %>%
  mutate(FDR = p.adjust(`𝘱 value`, method="bonferroni")) %>%
  mutate(sig = ifelse(`FDR` < 0.05, "***",
                      ifelse(`FDR` < 0.20, "*",
                             ifelse(`𝘱 value` < 0.05, "+", NA)))) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = factor(Variable, levels=c("Time:Response", "Status", "Location"))) %>%
  ggplot(aes(x=R2, y=Variable))+
  geom_bar(stat="identity", color="white", fill="black", width=0.5)+
  # scale_color_manual(values=c("*" = "grey", "***" = "black"))+
  geom_text(aes(x=R2+0.015, label=sig, vjust=ifelse(sig == "+", 1, 0.8)), size=6)+
  scale_x_continuous(limits=c(0,0.2))+
  theme_classic()+
  labs(x=expression(R^2), y=NULL)+
  theme(panel.grid.major.x = element_line(color="grey", linewidth=0.2, linetype=2))

# Among RS treated patients, there is still no significant difference in host proteome
# when controlling for inflammation status

proteomics.clean.pca.resp.plot = ggplot(proteomics.clean.pca.df %>%
         subset(Group == "RS") %>%
         mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
       aes(x=PC1, y=PC2))+
 # geom_path(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(shape = Location, fill=flare.group), color="white", size=2.5)+
  stat_ellipse(aes(color=flare.group), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = code),size=2)+
  scale_shape_manual(values=c(21,22,23))+
  guides(shape="none", fill="none")+
  theme_classic()+theme(#legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))+
  facet_wrap(~Location, scales="free")+
  labs(x=paste("PC1:", round(proteomics.clean.pca.var, digits=4)[1]*100, "%"),
       y=paste("PC2:", round(proteomics.clean.pca.var, digits=4)[2]*100, "%"),
       fill="Response", color="Response")
proteomics.clean.pca.resp.plot

(proteomics.clean.pca.plot+
    proteomics.clean.pca.permanova.plot+
    patchwork::plot_layout(widths=c(2,1)))/
(proteomics.clean.pca.resp.plot+
   proteomics.clean.pca.permanova.resp.plot+
   patchwork::plot_layout(widths=c(2,1)))

# Conclusion: clinical responders have no altered 
# host proteome relative to non-responders after the RS treatment

# Next question: Which proteins are different?


# :: MPX Interactions -------------------------------------------------------

# Question: are there significant differences before and after RS relative to Placebo

# loop through locations and perform LMER
# work with "proteomics.data.clean", which has been QC'ed

if(rerun==T){
  t1 = Sys.time()
  proteomics.data.group.interaction.lm = do.call(rbind, lapply(c("combined"), function(location){
    # subset to location & time
    metadata.subset = subset(proteomics.metadata, Time %in% c(0,1))
    data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
    # keep intersecting samples (shared by metadata and data)
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location (or, for combined, all locations)
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      # this should not be necessary // defunct, since we filtered to 80% prevalence above
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 0){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        metadata.subset$Plate = as.factor(metadata.subset$Plate) # likely overfits
        
        lmer.results = lmerTest::lmer(protein ~ Group*Time + Location + Status + (1|study_id) + (1|Plate), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 3 min
  # calculate padj
  proteomics.data.group.interaction.lm = proteomics.data.group.interaction.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.group.interaction.lm, "./host_proteomics/proteomics.data.group.interaction.lm.Rds")
}
proteomics.data.group.interaction.lm = readRDS("./host_proteomics/proteomics.data.group.interaction.lm.Rds")

# add mito linker
colnames(proteomics.data.group.interaction.lm)[colnames(proteomics.data.group.interaction.lm) == "protein"] = "protein.x"
proteomics.data.group.interaction.lm = merge(proteomics.data.group.interaction.lm,
                                       proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.group.interaction.lm.mito = merge(proteomics.data.group.interaction.lm,
                                            mitocarta.pathways.expanded,
                                            by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# arrange by significance
proteomics.data.group.interaction.lm.mito %>% arrange(padj)

subset(proteomics.data.group.interaction.lm.mito, !is.na(proteomics.data.group.interaction.lm.mito))$protein %>% unique() %>% length()

subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)$protein %>% unique() %>% length()
subset(proteomics.data.group.interaction.lm.mito, pval < 0.05 & coef > 0)$protein %>% unique() %>% length() # 56 up
subset(proteomics.data.group.interaction.lm.mito, pval < 0.05 & coef < 0)$protein %>% unique() %>% length() # 118 down

# :: MPX Volcano ----------------------------------------------------------

# volcano plots
proteomics.data.group.interaction.lm.mito.volcano = ggplot(proteomics.data.group.interaction.lm.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% 
                                                             mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                     aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                               color=mito),
                           size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Treatment Group")+
  labs(x="Interaction Coefficient",
       y="𝘱 value")
proteomics.data.group.interaction.lm.mito.volcano

# no proteins with FDR-adjusted p value < 0.20 (or even 0.80)

# treatment - decreased: 
# PFKAL = ATP-dependent 6-phosphofructokinase
# TBB3 = tubulin; cytoskeleton
# F6XY72 = nucleoside-disphosphate kinase; phosphate transfer
# ACOX1 = acyl-CoA oxidase; very LCFA oxidation
# S2512 = SLC25A12; malate-aspartate NADH shuttle; mitochondrial
# PSDE = proteosome
# MYOM1 = muscle structure
# FABP5 = fatty acid binding
# ATPO = mitochondrial
# ELAF = serine protease inhibitor; targets neutrophil elastase
# IF5A1 = translation initiation factor
# C1QBP = mitochondrial regulator
# PGM1 = phosphoglucomutase 1; glucose metabolism
# DX39B = spliceosome
# ETFB = electron transfer
# SUCB2 = succinate-CoA ligase; TCA cycle

subset(proteomics.data.group.interaction.lm.mito, protein == "SUCB2")
# S2512 = not present
# ATPO = decreased (consistent)
# C1QBP = increased (different)
# ETFB = increased (different)
# SUCB2 = increased (different)


# :: MPX: Mito Enrich -----------------------------------------------------

subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table()%>%
  fisher.test()
# p-value = 0.5885
# OR = 0.68
# not sig with GO mito either
subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito2) %>% table()%>%
  fisher.test()

proteomics.data.group.interaction.lm.mito.data = subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

proteomics.data.group.interaction.lm.mito.plot = proteomics.data.group.interaction.lm.mito.data%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
 # geom_point(shape=21, aes(size=perc), fill="white")+
 # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(proteomics.data.group.interaction.lm.mito.data, mito == "Mitochondrial"),
            aes(
              x=direction, y = .97*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="black", fill="white", size=4)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
proteomics.data.group.interaction.lm.mito.plot



# :: __ mito sensitivity ----------------------------------------------------------

# goal: see whether:
# A) p val threshold matters
# B) MitoCarta vs GO matters

# A: p val sensitivity
mito.enrich.sens.group = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.group.interaction.lm.mito, pval < x)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito.enrich.sens.group.plot = mito.enrich.sens.group %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Group ~ MitoCarta3.0")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

# B) GO sensitivity

# A: p val sensitivity
mito2.enrich.sens.group = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.group.interaction.lm.mito, pval < x)[,c("gene", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito2) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito2.enrich.sens.group.plot = mito2.enrich.sens.group %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Group ~ GO")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

mito.enrich.sens.group.plot+
mito2.enrich.sens.group.plot

# robust to p values


# :: __ oxphos sensitivity ----------------------------------------------------------

# goal: see if oxphos is also enriched
unique(proteomics.data.group.interaction.lm.mito$MitoPathways.Hierarchy)

subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)$protein %>% unique() $%>% length()
# n = 174 total sig proteins
data.frame(subset(proteomics.data.group.interaction.lm.mito, pval < 0.05) %>%
  mutate(oxphos = ifelse(grepl("OXPHOS", MitoPathways.Hierarchy), "oxphos", "")))[,c("protein", "coef", "pval", "oxphos")] %>% distinct() %>%
  # need to delete duplicated genes
  group_by(protein) %>%
  slice_max(order_by=oxphos, n = 1) %>% data.frame() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, oxphos) %>% table() %>%
  fisher.test()
# not sig; OR 0, p = 0.3066


# what's happening?
data.frame(subset(proteomics.data.group.interaction.lm.mito, pval < 0.05) %>%
             mutate(oxphos = ifelse(grepl("OXPHOS", MitoPathways.Hierarchy), "oxphos", "")))[,c("protein", "coef", "pval", "oxphos")] %>% distinct() %>%
  subset(oxphos=="oxphos")
# oxphos proteins are ONLY decreased in RS-treated (n=4 proteins)

# :: MPX: Heatmap ---------------------------------------------------------

# plot a one dimensional heatmap (with yellow-green-purple color)

proteomics.data.group.interaction.lm.mito.heatmap.data = subset(proteomics.data.group.interaction.lm.mito, pval < 0.05) %>%
  mutate(mito = ifelse(mito == "mito", " ", "  ")) %>%
  dplyr::select(mito, protein, coef)  %>% distinct() %>%
  mutate(pro.rank = rank(coef))

# add mito gene annotations
proteomics.data.group.interaction.lm.mito.heatmap.data.proteins = proteomics.data.group.interaction.lm.mito.heatmap.data %>%
  # subset to sig mito proteins
  subset(mito == " ") %>% dplyr::select(protein, pro.rank, coef) %>% distinct() %>%
  # give them a position along the x axis
  arrange(coef) %>%
  mutate(rank.pos = rank(coef))%>%
  mutate(rel.pos = min(rank.pos) + (((rank.pos - min(rank.pos)) * ((length(unique(proteomics.data.group.interaction.lm.mito.heatmap.data$protein)))-1))/
                                      (max(rank.pos)-1)))


proteomics.data.group.interaction.lm.mito.heatmap = ggplot(proteomics.data.group.interaction.lm.mito.heatmap.data,
       aes(x=reorder(protein, coef), y=1))+
  geom_tile(aes(fill=coef), color="black")+
  geom_point(aes(x=reorder(protein, coef), y=2, color=mito), shape=73, size=5)+
  geom_tile(aes(y=5), fill="white")+
  geom_text(data=proteomics.data.group.interaction.lm.mito.heatmap.data.proteins,
            aes(x=rel.pos, y=4, label=protein), angle=90, size=3, hjust=0, nudge_y=0)+
  scale_color_manual(values=c(" " = "red"), na.value=NA)+
  geom_segment(data = proteomics.data.group.interaction.lm.mito.heatmap.data.proteins,
               aes(x=pro.rank, xend=rel.pos, y=2.5, yend= 3.7), linetype=1, linewidth=0.2)+
  #guides(color="none")+
  viridis::scale_fill_viridis()+
  scale_x_discrete(expand = c(0.01, 1))+
  theme_minimal()+
  theme(legend.position="top",
        panel.grid=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())+
  labs(fill="Coefficient", title=NULL, x=NULL, y=NULL, color="Mitochondrial")+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
proteomics.data.group.interaction.lm.mito.heatmap



# :: ----------------------------------------------------------------------


# :: MPX Interactions: Responders ------------------------------------------
rerun=T
if(rerun==T){
  t1 = Sys.time()
  proteomics.data.group.interaction.resp.lm = do.call(rbind, lapply("combined", function(location){
    # subset to location (or combined) & RS group
    metadata.subset = subset(proteomics.metadata, Group == "RS" & Time %in% c(0,1))
    data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location (or combined)
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      # // defunct
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 0){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        metadata.subset$Plate = as.factor(metadata.subset$Plate)
        
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Location + Status + (1|study_id) + (1|Plate), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 7 min
  # calculate padj
  proteomics.data.group.interaction.resp.lm = proteomics.data.group.interaction.resp.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.group.interaction.resp.lm, "./host_proteomics/proteomics.data.group.interaction.resp.lm.Rds")
}
proteomics.data.group.interaction.resp.lm = readRDS("./host_proteomics/proteomics.data.group.interaction.resp.lm.Rds")

# add mito linker
colnames(proteomics.data.group.interaction.resp.lm)[colnames(proteomics.data.group.interaction.resp.lm) == "protein"] = "protein.x"
proteomics.data.group.interaction.resp.lm = merge(proteomics.data.group.interaction.resp.lm,
                                            proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.group.interaction.resp.lm.mito = merge(proteomics.data.group.interaction.resp.lm,
                                                 mitocarta.pathways.expanded,
                                                 by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# arrange by significance
proteomics.data.group.interaction.resp.lm.mito %>% arrange(pval)

# if we had only assessed mitochondrial proteins:
proteomics.data.group.interaction.resp.lm.mito %>%
  subset(mito == "mito") %>%
  dplyr::select(protein, coef, pval) %>%
  distinct() %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)
# still not passing FDR < 0.20

# :: MPX Volcano: Responders ----------------------------------------------------------

# volcano plots
proteomics.data.group.interaction.resp.lm.mito.volcano = ggplot(proteomics.data.group.interaction.resp.lm.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                          aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Therapy Response")+
  labs(x="Interaction Coefficient",
       y="𝘱 value")
proteomics.data.group.interaction.resp.lm.mito.volcano

# no longer < 0.20
# S10A9 protein is S100A9 gene (calprotectin)
# so, clinical responders have a reduction in calprotectin
# makes sense

proteomics.data.group.interaction.resp.lm.mito %>%
  arrange(pval)

# MCT1 ?
subset(proteomics.data.group.interaction.resp.lm.mito, gene == "SLC16A1") # p = 0.02

# H2S detox?
subset(proteomics.data.group.interaction.resp.lm.mito, grepl("SQR", gene)) # p = 0.04 (increased)
subset(proteomics.data.group.interaction.resp.lm.mito, grepl("ETHE", gene)) # p = 0.98
subset(proteomics.data.group.interaction.resp.lm.mito, grepl("TST", gene)) # TSTD1 p = 0.25 (reduced)


# "underpowered, but shows a trend towards increased oxidative phosphorylation and reduced inflammation"

# confirm this:

# :: MPX: Responders Mito Enrich -----------------------------------------------------

# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito = subset(proteomics.data.group.interaction.resp.lm.mito, mito == "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.other = subset(proteomics.data.group.interaction.resp.lm.mito, mito != "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.mito$mito = "mito"
proteomics.resp.other$mito = "other"
proteomics.resp.mito.other = rbind(proteomics.resp.other,
                                   proteomics.resp.mito)
# any edge cases:
table(proteomics.resp.mito$protein) %>% range()
table(proteomics.resp.other$protein) %>% range()
# no; perfect, continue

subset(proteomics.resp.mito.other, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito = factor(mito, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito) %>% table() %>%
  fisher.test()
# significantly enriched for mitochondrial genes
# p-value = 3.661e-07
# OR = 7.05
# same! good

# CHECK GO MITO (also sig)
# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito2 = subset(proteomics.data.group.interaction.resp.lm.mito, mito2 == "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.other2 = subset(proteomics.data.group.interaction.resp.lm.mito, mito2 != "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.mito2$mito2 = "mito"
proteomics.resp.other2$mito2 = "other"
proteomics.resp.mito.other2 = rbind(proteomics.resp.other2,
                                   proteomics.resp.mito2)
subset(proteomics.resp.mito.other2, pval < 0.05)[,c("protein", "coef", "pval", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito2 = factor(mito2, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito2) %>% table() %>%
  fisher.test()
# significantly enriched for mitochondrial genes
# p-value = 1.43e-07
# OR = 6.33
# similar! good



# :: __ mito2 sensitivity -------------------------------------------------------

# goal: see whether:
# A) p val threshold matters
# B) MitoCarta vs GO matters

# A: p val sensitivity
mito.enrich.sens.resp = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.group.interaction.resp.lm.mito, pval < x)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito.enrich.sens.resp.plot = mito.enrich.sens.resp %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ MitoCarta3.0")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

# B) GO sensitivity
mito2.enrich.sens.resp = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.group.interaction.resp.lm.mito, pval < x)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito2) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito2.enrich.sens.resp.plot = mito2.enrich.sens.resp %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ GO")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

mito.enrich.sens.resp.plot+
  mito2.enrich.sens.resp.plot
# Robust across p value thresholds

proteomics.mito.gene.enrichment.data = subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

proteomics.mito.gene.enrichment.plot =proteomics.mito.gene.enrichment.data%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(proteomics.mito.gene.enrichment.data, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="black",fill="white",size=4)+
  geom_label(data=subset(proteomics.mito.gene.enrichment.data, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .9*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="black",fill="white", size=5)+
  
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
proteomics.mito.gene.enrichment.plot


# :: __ oxphos sensitivity ----------------------------------------------------------

# goal: see if oxphos is also enriched
unique(proteomics.data.group.interaction.resp.lm.mito$MitoPathways.Hierarchy)

data.frame(subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05) %>%
             mutate(oxphos = ifelse(grepl("OXPHOS", MitoPathways.Hierarchy), "oxphos", "")))[,c("protein", "coef", "pval", "oxphos")] %>% distinct() %>%
  # need to delete duplicated genes
  group_by(protein) %>%
  slice_max(order_by=oxphos, n = 1) %>% data.frame() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, oxphos) %>% table()%>%
  fisher.test()
# sig enriched by RS

# what's happening?
data.frame(subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05) %>%
             mutate(oxphos = ifelse(grepl("OXPHOS", MitoPathways.Hierarchy), "oxphos", "")))[,c("protein", "coef", "pval", "oxphos")] %>% distinct() %>%
  subset(oxphos=="oxphos")
# oxphos proteins are ONLY increased in responders (n=8 proteins)

# :: MPX: Heatmap ---------------------------------------------------------

# plot a one dimensional heatmap (with yellow-green-purple color)

proteomics.resp.mito.other.heatmap = ggplot(subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05)%>%
                                              mutate(mito = ifelse(mito == "mito", " ", "  ")),
                                                           aes(x=reorder(protein, coef), y=1))+
  geom_tile(aes(fill=coef), color="black")+
  geom_point(aes(x=reorder(protein, coef), y=2, color=mito), shape=73, size=5)+
  scale_color_manual(values=c(" " = "red"), na.value=NA)+
  #guides(color="none")+
  viridis::scale_fill_viridis()+
  theme_void()+
  theme(legend.position="top")+
  labs(fill="Coefficient", title=NULL, x=NULL, y=NULL, color="Mitochondrial")
proteomics.resp.mito.other.heatmap


proteomics.data.group.interaction.lm.mito.resp.heatmap.data = subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05) %>%
  mutate(mito = ifelse(mito == "mito", " ", "  ")) %>%
  dplyr::select(mito, protein, coef)  %>% distinct() %>%
  mutate(pro.rank = rank(coef))

# add mito gene annotations
proteomics.data.group.interaction.lm.mito.resp.heatmap.data.proteins = proteomics.data.group.interaction.lm.mito.resp.heatmap.data %>%
  # subset to sig mito proteins
  subset(mito == " ") %>% dplyr::select(protein, pro.rank, coef) %>% distinct() %>%
  # give them a position along the x axis
  arrange(coef) %>%
  mutate(rank.pos = rank(coef))%>%
  mutate(rel.pos = min(rank.pos) + (((rank.pos - min(rank.pos)) * ((length(unique(proteomics.data.group.interaction.lm.mito.resp.heatmap.data$protein)))-1))/
                                      (max(rank.pos)-1)))


proteomics.data.group.interaction.lm.mito.resp.heatmap = ggplot(proteomics.data.group.interaction.lm.mito.resp.heatmap.data,
                                                           aes(x=reorder(protein, coef), y=1))+
  geom_tile(aes(fill=coef), color="black")+
  geom_point(aes(x=as.factor(reorder(protein, coef)), y=2, color=mito), shape=73, size=5)+
  geom_tile(aes(y=5), fill="white")+
  geom_text(data=proteomics.data.group.interaction.lm.mito.resp.heatmap.data.proteins,
            aes(x=rel.pos, y=4, label=protein), angle=90, size=3, hjust=0, nudge_y=0)+
  scale_color_manual(values=c(" " = "red"), na.value=NA)+
  geom_segment(data = proteomics.data.group.interaction.lm.mito.resp.heatmap.data.proteins,
               aes(x=pro.rank, xend=rel.pos, y=2.5, yend= 3.7), linetype=1, linewidth=0.2)+
  #guides(color="none")+
  viridis::scale_fill_viridis()+
  scale_x_discrete(expand = c(0.01 , 1))+
  theme_minimal()+
  theme(legend.position="top",
        panel.grid=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())+
  labs(fill="Coefficient", title=NULL, x=NULL, y=NULL, color="Mitochondrial")+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
proteomics.data.group.interaction.lm.mito.resp.heatmap



# :: MPX: Interpret -------------------------------------------------------

subset(proteomics.data.group.interaction.resp.lm.mito, 
       grepl("OXPHOS", MitoPathway)) %>%
  arrange(-abs(coef))
# SDHA, SURF1, SDHB are enriched and == OXPHOS
# NDUF == OXPHOS, but several NDUF are decreased, too

subset(proteomics.data.group.interaction.resp.lm.mito, 
       grepl("TCA", MitoPathway))
# SDHA, SDHB,  IDH3G are enriched and == TCA cycle

subset(proteomics.data.group.interaction.resp.lm.mito, 
       grepl("Fatty acid", MitoPathway))
# ACADM is enriched and == Fatty acid oxidation

# SDH = Succinate dehydrogenase
# ACADM = Medium-chain Acyl-CoA dehydrogenase

proteomics.data.group.interaction.resp.lm.mito %>% subset(mito == "mito") %>% arrange(coef)
# PISD is depleted and == Phospholipid metabolism
# phosphatidylserine decarboxylase (related to ethanolamine)



# :: :: PLACEBO ----------------------------------------------------------------------

# repeat with Placebo group to see if the same trends hold, or if the associations are specific to RS
# note: the responder/non-responder analysis is underpowered! only n=1 Placebo responder with biopsy data

# :: :: MPX Interactions: Responders ------------------------------------------


  t1 = Sys.time()
  proteomics.data.placebo.interaction.resp.lm = do.call(rbind, lapply("combined", function(location){
    # subset to location (or combined) & RS group
    metadata.subset = subset(proteomics.metadata, Group == "Placebo" & Time %in% c(0,1))
    data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location (or combined)
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      # // defunct, because of 80% prevalence filter
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 0){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 7 min
  # calculate padj
  proteomics.data.placebo.interaction.resp.lm = proteomics.data.placebo.interaction.resp.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
# add mito linker
colnames(proteomics.data.placebo.interaction.resp.lm)[colnames(proteomics.data.placebo.interaction.resp.lm) == "protein"] = "protein.x"
proteomics.data.placebo.interaction.resp.lm = merge(proteomics.data.placebo.interaction.resp.lm,
                                                  proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.placebo.interaction.resp.lm.mito = merge(proteomics.data.placebo.interaction.resp.lm,
                                                       mitocarta.pathways.expanded,
                                                       by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# arrange by significance
proteomics.data.placebo.interaction.resp.lm.mito %>% arrange(pval)

# :: :: MPX Volcano: Responders ----------------------------------------------------------

# volcano plots
proteomics.data.placebo.interaction.resp.lm.mito.volcano = ggplot(proteomics.data.placebo.interaction.resp.lm.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                                aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Therapy Response")+
  labs(x="Interaction Coefficient",
       y="𝘱 value")
proteomics.data.placebo.interaction.resp.lm.mito.volcano


proteomics.data.placebo.interaction.resp.lm.mito %>%
  arrange(pval)

# :: :: MPX: Responders Mito Enrich -----------------------------------------------------

# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito.placebo = subset(proteomics.data.placebo.interaction.resp.lm.mito, mito == "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.other.placebo = subset(proteomics.data.placebo.interaction.resp.lm.mito, mito != "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.mito.placebo$mito = "mito"
proteomics.resp.other.placebo$mito = "other"
proteomics.resp.mito.other.placebo = rbind(proteomics.resp.other.placebo,
                                   proteomics.resp.mito.placebo)
# any edge cases:
table(proteomics.resp.mito.placebo$protein) %>% range()
table(proteomics.resp.other.placebo$protein) %>% range()
# no; perfect, continue

subset(proteomics.resp.mito.other.placebo, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito = factor(mito, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito) %>% table() %>%
  fisher.test()
# not sig

# CHECK GO MITO (also sig)
# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito2.placebo = subset(proteomics.data.placebo.interaction.resp.lm.mito, mito2 == "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.other2.placebo = subset(proteomics.data.placebo.interaction.resp.lm.mito, mito2 != "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.mito2.placebo$mito2 = "mito"
proteomics.resp.other2.placebo$mito2 = "other"
proteomics.resp.mito.other2.placebo = rbind(proteomics.resp.other2.placebo,
                                    proteomics.resp.mito2.placebo)
subset(proteomics.resp.mito.other2.placebo, pval < 0.05)[,c("protein", "coef", "pval", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito2 = factor(mito2, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito2) %>% table() %>%
  fisher.test()
# 

proteomics.data.resp.interaction.lm.mito.placebo = subset(proteomics.resp.mito.other.placebo, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

proteomics.data.resp.interaction.placebo.lm.mito.plot = proteomics.data.resp.interaction.lm.mito.placebo%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(proteomics.data.resp.interaction.lm.mito.placebo, mito == "Mitochondrial"),
            aes(
              x=direction, y = .97*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
proteomics.data.resp.interaction.placebo.lm.mito.plot
# same results as RS (7% vs 7%)


# :: :: __ mito2 sensitivity -------------------------------------------------------

# goal: see whether:
# A) p val threshold matters
# B) MitoCarta vs GO matters

# A: p val sensitivity
mito.enrich.sens.resp.placebo = do.call(rbind, lapply(seq(from=0.05, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.placebo.interaction.resp.lm.mito, pval < x)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito.enrich.sens.resp.placebo.plot = mito.enrich.sens.resp.placebo %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ MitoCarta3.0 (Placebo)")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw p value)",
       y="Enrichment 𝘱 value")

# B) GO sensitivity
mito2.enrich.sens.resp.placebo = do.call(rbind, lapply(seq(from=0.05, to=0.25, by=0.01), function(x){
  fisher.result = subset(proteomics.data.placebo.interaction.resp.lm.mito, pval < x)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito2) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
mito2.enrich.sens.resp.placebo.plot = mito2.enrich.sens.resp.placebo %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.05, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ GO (Placebo)")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw p value)",
       y="Enrichment 𝘱 value")

mito.enrich.sens.resp.placebo.plot+
  mito2.enrich.sens.resp.placebo.plot
# 

# :: all plots ------------------------------------------------------------
(proteomics.data.group.interaction.lm.mito.volcano+
    proteomics.data.group.interaction.lm.mito.plot+
    proteomics.data.group.interaction.resp.lm.mito.volcano+
    proteomics.mito.gene.enrichment.plot)/
  (proteomics.data.placebo.interaction.resp.lm.mito.volcano+
   mito2.enrich.sens.resp.placebo.plot)



# >>> 2. ASV MLI --------------------------------------------------------------


# :: process ASV ----------------------------------------------------------

# load disease status
trials.disease.activity = read.csv("./master_metadata_trials_241030.csv") %>% subset(grepl("ASP", standard.name))
trials.disease.activity = trials.disease.activity[,c("standard.name", "sample.inflammation.status")]

# load 16S data
load("./all_trials_16S_min50k_mappingFile_250530_gg2species_allData_250910.Rdata")
lsarp.mli.ps <- physeq.pooled.wTree
# subset to patients with biopsy data and ASP
lsarp.mli.ps = phyloseq::subset_samples(lsarp.mli.ps, grepl("ASP", standard.name) & 
                                          study_id %in% unique(proteomics.metadata$study_id) & 
                                          trial.samples == "Resistant Starch Trial")

# merge duplicates (add reads together, which will be rarefied later)
lsarp.mli.ps.otu = speedyseq::psmelt(lsarp.mli.ps)
lsarp.mli.ps.otu = lsarp.mli.ps.otu %>% reshape2::acast(standard.name ~ OTU, value.var="Abundance", fun.aggregate=sum)
lsarp.mli.ps.meta = phyloseq::sample_data(lsarp.mli.ps) %>% data.frame()
lsarp.mli.ps.meta = lsarp.mli.ps.meta[,c("study_id", "standard.name", "trial.patient.status", "baseline.diagnosis",
                                         "fecalcal_res")] %>% distinct() %>% data.frame()
rownames(lsarp.mli.ps.meta) = lsarp.mli.ps.meta$standard.name
lsarp.mli.ps.otu.ps = phyloseq::phyloseq(phyloseq::otu_table(lsarp.mli.ps.otu, taxa_are_rows=F),
                                         phyloseq::sample_data(lsarp.mli.ps.meta),
                                         phyloseq::tax_table(lsarp.mli.ps))
# rarefy
set.seed(25)
lsarp.mli.ps.otu.ps = phyloseq::rarefy_even_depth(lsarp.mli.ps.otu.ps, sample.size=50000)
# edit taxonomy
lsarp.mli.ps.df = speedyseq::psmelt(lsarp.mli.ps.otu.ps)
lsarp.mli.ps.df
# analyze
lsarp.mli.ps.df = lsarp.mli.ps.df %>%
  mutate(LCA = 
  ifelse(!is.na(Species), paste(Genus, Species, sep="_"),
         ifelse(!is.na(Genus) & is.na(Species), Genus,
                ifelse(!is.na(Family), Family,
                       ifelse(!is.na(Order), Order,
                              ifelse(!is.na(Class), Class,
                                     ifelse(!is.na(Phylum), Phylum, Kingdom)))))))
lsarp.mli.ps.df$Abundance = lsarp.mli.ps.df$Abundance / 50000                    
lsarp.pseudo = min(lsarp.mli.ps.df$Abundance[lsarp.mli.ps.df$Abundance != 0])/2
# make ASV table (%)
lsarp.mli.ps.mat = reshape2::acast(lsarp.mli.ps.df,
                                   standard.name ~ OTU, value.var = "Abundance",
                                   fun.aggregate=sum)
# make LCA table (%)
lsarp.mli.ps.glom.mat = reshape2::acast(lsarp.mli.ps.df,
                                   standard.name ~ LCA, value.var = "Abundance",
                                   fun.aggregate=sum)

dim(lsarp.mli.ps.mat)
# make mapping file
lsarp.mli.ps.meta = lsarp.mli.ps.df[,c("Sample","study_id","standard.name")] %>% distinct() %>% data.frame()    
# unpack meta from standard.name
lsarp.mli.ps.meta = lsarp.mli.ps.meta %>%
  tidyr::separate(Sample, into=c("HM", "type", "Location"), sep="-", remove=F) %>%
  tidyr::separate(HM, into=c("HM", "Time"), sep="\\.0", remove=F)
# merge with other meta
lsarp.mli.ps.meta = merge(lsarp.mli.ps.meta,
                          trials.disease.activity, by="standard.name")
colnames(lsarp.mli.ps.meta)[colnames(lsarp.mli.ps.meta) == "sample.inflammation.status"] = "Status"         
lsarp.mli.ps.meta$Status = ifelse(lsarp.mli.ps.meta$Status == "Non-inflamed", "N", "A")
# add group
lsarp.mli.ps.meta = merge(lsarp.mli.ps.meta,
                          subset(unblinding, !grepl("excluded", flare.call_verified))[,c("study_id", "Group")], by="study_id")
lsarp.mli.ps.meta$sample_name = paste(lsarp.mli.ps.meta$HM, lsarp.mli.ps.meta$Time, sep=".0")
nrow(distinct(lsarp.mli.ps.meta)) # n = 137

# match to proteomics samples
lsarp.mli.ps.meta = subset(lsarp.mli.ps.meta, sample_name %in% paste(proteomics.clean.pca.df$HM, proteomics.clean.pca.df$Time, sep=".0")) %>% distinct()
nrow(distinct(lsarp.mli.ps.meta)) # n = 124

table(lsarp.mli.ps.meta[,c("HM", "Time")])

# :: Butyrogens I -----------------------------------------------------------

# Use the original method to annotate butyrogens
# (so we can use a Fisher test on enriched/depleted)

# requires reassigning taxonomy with gg138, but can use original ASV-level data
library("dada2")
set.seed(25)
lsarp.mli.gg138.tax = dada2::assignTaxonomy((phyloseq::refseq(lsarp.mli.ps)), "~/Documents/PhD/16S_databases/gg_13_8_train_set_97.fa.gz")
lsarp.mli.gg138.tax.df = data.frame(lsarp.mli.gg138.tax)
lsarp.mli.gg138.tax.df$OTU = names(phyloseq::refseq(lsarp.mli.ps)) 

# add LCA
lsarp.mli.gg138.tax.df <- data.frame(lsarp.mli.gg138.tax.df) %>%
  mutate(lca.gg138 = 
           ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                  ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                         ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                       ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                              ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                     ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))

# flag whether the ASV is a Butyrogen
lsarp.mli.gg138.tax.df = lsarp.mli.gg138.tax.df %>%
  mutate(Species = ifelse(is.na(Species), "", Species),
         Genus =ifelse(is.na(Genus), "", Genus),
         Family = ifelse(is.na(Family), "", Family)) %>%
  mutate(butyrogen.i = ifelse(Family=="f__Lachnospiraceae" | Genus=="g__Blautia" | Genus=="g__Roseburia" | Genus=="g__Eubacterium" | Genus=="g__Ruminococcus" | Genus=="g__Clostridium" | Genus=="g__Faecalibacterium", "butyrogen", "other"))

lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")]

# save this
saveRDS(lsarp.mli.gg138.tax.df, "./2025_09_14_lsarp_mli_butyrogens_i.Rds")


lsarp.mli.gg138.tax.df = readRDS("./2025_09_14_lsarp_mli_butyrogens_i.Rds")



# :: Butyrogens II -----------------------------------------------------------

# Use the Kircher method to annotate butyrogens
# (so we can use a Fisher test on enriched/depleted)

# requires Kircher method

# identify ASVs present in 1+ sample
subset(lsarp.mli.ps.df, Abundance >0)$OTU %>% unique() %>% length()
# n = 3594

lsarp.mli.refseq = lsarp.mli.ps@refseq[names(lsarp.mli.ps@refseq) %in% rownames(phyloseq::tax_table(lsarp.mli.ps.otu.ps)),]
# reduced to 3594 ASVs present

run.picrust = FALSE
if(run.picrust == T) {
  # note: make new picrust2 instance
  # use this data:
  
  # prepare seq.table for picrust2
  fasta_seqs <- Biostrings::DNAStringSet(data.frame(lsarp.mli.refseq)[,1])
  names(fasta_seqs) <- rownames(data.frame(lsarp.mli.refseq))  # assign ASV IDs as sequence names
  # prepare asv.table for picrust2
  lsarp.mli.picrust2.abuntable = data.frame(phyloseq::otu_table(lsarp.mli.ps))
  lsarp.mli.picrust2.abuntable = lsarp.mli.picrust2.abuntable[,colnames(lsarp.mli.picrust2.abuntable) %in%
                                                                names(fasta_seqs)]
  # export
  Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/lsarp.mli.picrust2.seqtable.fasta")
  Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/predict_SCFA_producers/lsarp.mli.picrust2.seqtable.fasta")
  write.table(t(lsarp.mli.picrust2.abuntable), "oars_picrust2/lsarp.mli.picrust2.abuntable.tsv", 
              sep="\t", quote = F, col.names = NA)
  
  # STEP 1: run through default picrust2 (to ensure picrust2 works)
  cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2
  conda activate oars_picrust2
  # In R, may need to install Rcpp, jsonlite, lattice, Matrix, RSpectra, castor
  rm -r lsarp_mli_picrust2_out
  picrust2_pipeline.py \
  -s lsarp.mli.picrust2.seqtable.fasta \
  -i lsarp.mli.picrust2.abuntable.tsv \
  -o lsarp_mli_picrust2_out \
  -p 1
  
  # optionally perform stratification to source functions to taxa
  # use EC's (not preferable; 32 million rows)
  metagenome_pipeline.py -i lsarp.mli.picrust2.abuntable.tsv -m lsarp_mli_picrust2_out/marker_predicted_and_nsti.tsv.gz -f lsarp_mli_picrust2_out/EC_predicted.tsv.gz \
  -o lsarp_mli_picrust2_out/EC_metagenome_out --strat_out
  
  # from: 
  
  # STEP 2: run through Vitals' predict_SCFA_producers to identify butyrogens
  # source: https://github.com/ag-vital/predict_SCFA_producers/tree/master
  # first: git clone https://github.com/ag-vital/predict_SCFA_producers.git
  # then, rename "picrust" folder as "SCFA"
  cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2/predict_SCFA_producers
  # when rerunning, delete placement_working first
  rm -r placement_working
  place_seqs.py -s lsarp.mli.picrust2.seqtable.fasta -o placed_seqs.tre -p 1 --intermediate placement_working --ref_dir SCFA
  # 97 sequences failed to align
  hsp.py -t placed_seqs.tre --observed_trait_table SCFA/SCFA_pathwaydata.txt -o SCFA_predicted.tsv -p 1 -m emp_prob -n
  
  # save files to new folder
  cp ./SCFA_predicted.tsv ../picrust2_saved/lsarp_mli_SCFA_predicted.tsv
  cp ../lsarp_mli_picrust2_out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz ../picrust2_saved/lsarp_mli_pred_metagenome_contrib.tsv
  
}
# import predicted SCFA data
lsarp.mli.vital.butyrogens = read.csv("./oars_picrust2/picrust2_saved/lsarp_mli_SCFA_predicted.tsv", sep="\t")
# Butyrogens = "Possesses a) AcetylCoA pathway, and at least one of: b) but or c) buk"
lsarp.mli.vital.butyrogens.acetylcoa = subset(lsarp.mli.vital.butyrogens, acetylcoa == 1)
lsarp.mli.vital.butyrogens.acetylcoa = subset(lsarp.mli.vital.butyrogens.acetylcoa, but == 1 | buk == 1)
# these will be "butyrogens"
lsarp.mli.vital.butyrogens.acetylcoa

# identify whether glom'ed taxa is butyrogen
lsarp.mli.butyrogens = lsarp.mli.ps.df %>%
  mutate(butyrogen.ii = ifelse(OTU %in% lsarp.mli.vital.butyrogens.acetylcoa$sequence,
                            "butyrogen", "other"))
lsarp.mli.butyrogens = lsarp.mli.butyrogens[,c("OTU", "LCA", "butyrogen.ii")]
colnames(lsarp.mli.butyrogens) = c("OTU", "lca.gg2", "butyrogen.ii")

saveRDS(lsarp.mli.butyrogens, "./2025_09_14_lsarp_mli_butyrogens_ii.Rds")

lsarp.mli.butyrogens = readRDS("./2025_09_14_lsarp_mli_butyrogens_ii.Rds")


# :: save -----------------------------------------------------------------

saveRDS(lsarp.mli.ps.mat, "./lsarp.mli.ps.mat.Rds")
saveRDS(lsarp.mli.ps.glom.mat, "./lsarp.mli.ps.glom.mat.Rds")
saveRDS(lsarp.mli.ps.meta, "./lsarp.mli.ps.meta.Rds")



# :: load -----------------------------------------------------------------


# load
# Glom'ed
lsarp.mli.ps.glom.mat = readRDS("./lsarp.mli.ps.glom.mat.Rds")
# ASVs
lsarp.mli.ps.mat = readRDS("./lsarp.mli.ps.mat.Rds")
lsarp.mli.ps.meta = readRDS("./lsarp.mli.ps.meta.Rds")

lsarp.mli.gg138.tax.df = readRDS("./2025_09_14_lsarp_mli_butyrogens_i.Rds")
lsarp.mli.butyrogens = readRDS("./2025_09_14_lsarp_mli_butyrogens_ii.Rds")

colnames(lsarp.mli.ps.meta) = gsub("Stasus", "Status", colnames(lsarp.mli.ps.meta))

distinct(lsarp.metadata.responders[,c("HM","Group", "flare.group")]) %>%
  arrange(Group, flare.group)
# missing HM0978 (the other placebo remit sample, therefore we're down to n=1 placebo remit)

ggplot(lsarp.mli.ps.meta %>% subset(Time != 2) %>%
         merge(distinct(lsarp.metadata.responders[,c("HM", "flare.group")]), by="HM") %>%
         mutate(Status = factor(Status, levels=c("N", "A")))%>%
         arrange((Status)),
       aes(x=as.factor(Time), y=HM))+
  ggbeeswarm::geom_beeswarm(aes(shape=Status, color=flare.group), 
                            alpha=0.5, size=2.5, cex=5)+
  scale_color_manual(values=c("red", "blue"))+
  #scale_shape_manual(values=c(24,25))+
  scale_shape_manual(values=c(17,16))+
  
  facet_grid(Group~Location, scales="free")+
  theme_classic()+
  labs(x="Scope", y="")

# :: ASV PCoA -----------------------------------------------------------------

lsarp.mli.ps.mat = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat)%in% lsarp.mli.ps.meta$standard.name,]

# 20% prev
ncol(lsarp.mli.ps.mat) # 3594
lsarp.mli.ps.mat.pa = lsarp.mli.ps.mat
lsarp.mli.ps.mat.pa[lsarp.mli.ps.mat.pa > 0] = 1
lsarp.mli.ps.mat = lsarp.mli.ps.mat[,colSums(lsarp.mli.ps.mat.pa) > nrow(lsarp.mli.ps.mat.pa)*0.2]
ncol(lsarp.mli.ps.mat) # 509

# calculate Bray-Curtis dissimilarities
lsarp.mli.ps.bray = vegan::vegdist(lsarp.mli.ps.mat)
# perform PCoA
lsarp.mli.ps.pcoa = ape::pcoa(lsarp.mli.ps.bray)
# extract data from pcoa
lsarp.mli.ps.pcoa.df = data.frame(lsarp.mli.ps.pcoa$vectors[,c(1:2)])
lsarp.mli.ps.pcoa.df$standard.name = rownames(lsarp.mli.ps.pcoa.df)
# add metadata
lsarp.mli.ps.pcoa.df = merge(lsarp.mli.ps.pcoa.df,
                          distinct(lsarp.mli.ps.meta), by="standard.name")
# extract variance explained
lsarp.mli.ps.pcoa.var = lsarp.mli.ps.pcoa$values[c(1:2),2]
lsarp.mli.ps.pcoa.df$var1 = round(lsarp.mli.ps.pcoa.var[1]*100, digits=2)
lsarp.mli.ps.pcoa.df$var2 = round(lsarp.mli.ps.pcoa.var[2]*100, digits=2)
rownames(lsarp.mli.ps.pcoa.df) = lsarp.mli.ps.pcoa.df$standard.name
# clean up "lsarp.on.rs" variable

# numbers
lsarp.mli.ps.pcoa.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()
lsarp.mli.ps.pcoa.df[,c("HM", "Group")] %>%  table()
lsarp.mli.ps.pcoa.df[,c("HM", "Group", "Location", "Time")] %>%
  dplyr::select(HM, Location, Time)  %>% table()

set.seed(25)
t1 = Sys.time()
lsarp.mli.ps.permanova = vegan::adonis2(lsarp.mli.ps.bray ~ Location + Status + Time*Group,
                                     lsarp.mli.ps.pcoa.df,
                                     strata = lsarp.mli.ps.pcoa.df$HM,
                                     by="margin")
t2 = Sys.time()
t2 - t1
lsarp.mli.ps.permanova 
# Time:Group interaction is very sig
# trends hold with glom'ed data

lsarp.mli.ps.permanova = as.data.frame(lsarp.mli.ps.permanova)[1:3, c(3,5)] %>%
  mutate(Variable = rownames(.)) %>%
  mutate(R2 = round(R2, digits=3))
lsarp.mli.ps.permanova = lsarp.mli.ps.permanova[,c(3,1,2)] %>%
  reshape2::melt()
lsarp.mli.ps.permanova.plot = ggplot(lsarp.mli.ps.permanova %>% 
                                               mutate(Variable = gsub("Time:Group", "Time:Group", Variable))%>%
                                               mutate(variable = factor(variable, levels=c("R2", "Pr(>F)")),
                                                      Variable = factor(Variable, levels=c("Location", "Status", "Time:Group"))),
                                             aes(x=variable, y=Variable))+
  geom_tile(fill="white", color="black")+
  geom_text(aes(label=value, fontface = ifelse(variable == "Pr(>F)" & `value` < 0.05, "bold", "plain")))+
  scale_x_discrete(position = "top") +
  theme_minimal()+theme(axis.text = element_text(size=12),
                        panel.grid.major=element_blank())+
  labs(x="", y="")
lsarp.mli.ps.permanova.plot


# make a plot, not a table
lsarp.mli.ps.permanova.r2.plot = lsarp.mli.ps.permanova %>%
  reshape2::acast(Variable ~ variable, value.var="value") %>% data.frame() %>%
  mutate(`𝘱 value` = `Pr..F.`) %>%
  mutate(FDR = p.adjust(`𝘱 value`, method="bonferroni")) %>%
  mutate(sig = ifelse(`FDR` < 0.05, "***",
                      ifelse(`FDR` < 0.20, "*",
                             ifelse(`𝘱 value` < 0.05, "+", NA)))) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = factor(Variable, levels=c("Time:Group", "Status", "Location"))) %>%
  ggplot(aes(x=R2, y=Variable))+
  geom_bar(stat="identity", color="white", fill="black", width=0.5)+
  # scale_color_manual(values=c("*" = "grey", "***" = "black"))+
  geom_text(aes(x=R2+0.015, label=sig, vjust=ifelse(sig == "+", 1, 0.8)), size=6)+
  scale_x_continuous(limits=c(0,0.2))+
  theme_classic()+
  labs(x=expression(R^2), y=NULL)+
  theme(panel.grid.major.x = element_line(color="grey", linewidth=0.2, linetype=2))

# check status
lsarp.mli.ps.pcoa.plot = ggplot(lsarp.mli.ps.pcoa.df %>% mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
       aes(x=Axis.1, y=Axis.2))+
  #geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(aes(shape = Location, fill=Group), color="white", size=2.5)+
  stat_ellipse(aes(color=Group), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = HM),size=2)+
  scale_shape_manual(values=c(21,22,23))+
  guides(shape = "none", fill="none")+
  theme_classic()+theme(#legend.position="none",
    plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~Location, scales="free")+
  labs(x=paste("Axis 1:", round(lsarp.mli.ps.pcoa.df$var1, digits=4)[1], "%"),
       y=paste("Axis 2:", round(lsarp.mli.ps.pcoa.df$var2, digits=4)[2], "%"))
lsarp.mli.ps.pcoa.plot


# :: PCoA Responders  --------------------------------------------


lsarp.metadata.responders # = readRDS("./2025_07_14_lsarp_cd_response_list.Rds") # already loaded

lsarp.mli.ps.pcoa.df = merge(lsarp.mli.ps.pcoa.df,
                                lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                                by="HM")
# redraw Bray-Curtis

# calculate Bray-Curtis dissimilarities
dim(lsarp.mli.ps.mat) # 509 taxa; 124 samples
lsarp.mli.ps.bray.resp = vegan::vegdist(lsarp.mli.ps.mat[grepl(paste(c(subset(lsarp.metadata.responders, Group == "RS")$HM), collapse="|"), rownames(lsarp.mli.ps.mat)),])
# perform PCoA
lsarp.mli.ps.pcoa.resp = ape::pcoa(lsarp.mli.ps.bray.resp)
# extract data from pcoa
lsarp.mli.ps.pcoa.resp.df = data.frame(lsarp.mli.ps.pcoa.resp$vectors[,c(1:2)])
lsarp.mli.ps.pcoa.resp.df$standard.name = rownames(lsarp.mli.ps.pcoa.resp.df)
# add metadata
lsarp.mli.ps.pcoa.resp.df = merge(lsarp.mli.ps.pcoa.resp.df,
                             distinct(lsarp.mli.ps.meta), by="standard.name")
lsarp.mli.ps.pcoa.resp.df = merge(lsarp.mli.ps.pcoa.resp.df,
                                  lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                                  by="HM")
# extract variance explained
lsarp.mli.ps.pcoa.resp.var = lsarp.mli.ps.pcoa.resp$values[c(1:2),2]
lsarp.mli.ps.pcoa.resp.df$var1 = round(lsarp.mli.ps.pcoa.resp.var[1]*100, digits=2)
lsarp.mli.ps.pcoa.resp.df$var2 = round(lsarp.mli.ps.pcoa.resp.var[2]*100, digits=2)
rownames(lsarp.mli.ps.pcoa.resp.df) = lsarp.mli.ps.pcoa.resp.df$standard.name
# clean up "lsarp.on.rs" variable

# numbers
lsarp.mli.ps.pcoa.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% distinct() %>% dplyr::select(flare.group) %>% table()
lsarp.mli.ps.pcoa.df[,c("HM", "Group")] %>%  subset(Group == "RS") %>% table()
lsarp.mli.ps.pcoa.df[,c("HM", "Group", "Location", "Time", "flare.group")] %>% subset(Group == "RS") %>% 
  dplyr::select(HM, Location, Time)  %>% table()


set.seed(25)
t1 = Sys.time()
lsarp.mli.ps.resp.permanova = vegan::adonis2(lsarp.mli.ps.bray.resp ~ Location + Status + Time*flare.group,
                                        lsarp.mli.ps.pcoa.resp.df,
                                        strata = lsarp.mli.ps.pcoa.resp.df$HM,
                                        by="margin")
t2 = Sys.time()
t2 - t1
lsarp.mli.ps.resp.permanova # Only Time:flare.group is significant!

lsarp.mli.ps.resp.permanova = as.data.frame(lsarp.mli.ps.resp.permanova)[1:3, c(3,5)] %>%
  mutate(Variable = rownames(.)) %>%
  mutate(R2 = round(R2, digits=3))
lsarp.mli.ps.resp.permanova = lsarp.mli.ps.resp.permanova[,c(3,1,2)] %>%
  reshape2::melt()
lsarp.mli.ps.resp.permanova.plot = ggplot(lsarp.mli.ps.resp.permanova %>% 
                                       mutate(Variable = gsub("Time:flare.group", "Time:Response", Variable))%>%
                                       mutate(variable = factor(variable, levels=c("R2", "Pr(>F)")),
                                              Variable = factor(Variable, levels=c("Location", "Status", "Time:Response"))),
                                     aes(x=variable, y=Variable))+
  geom_tile(fill="white", color="black")+
  geom_text(aes(label=value, fontface = ifelse(variable == "Pr(>F)" & `value` < 0.05, "bold", "plain")))+
  scale_x_discrete(position = "top") +
  theme_minimal()+theme(axis.text = element_text(size=12),
                        panel.grid.major=element_blank())+
  labs(x="", y="")
lsarp.mli.ps.resp.permanova.plot


# make a plot, not a table
lsarp.mli.ps.resp.permanova.r2.plot = lsarp.mli.ps.resp.permanova %>%
  mutate(Variable = gsub("Time:flare.group", "Time:Response", Variable))%>%
  reshape2::acast(Variable ~ variable, value.var="value") %>% data.frame() %>%
  mutate(`𝘱 value` = `Pr..F.`) %>%
  mutate(FDR = p.adjust(`𝘱 value`, method="bonferroni")) %>%
  mutate(sig = ifelse(`FDR` < 0.05, "***",
                      ifelse(`FDR` < 0.20, "*",
                             ifelse(`𝘱 value` < 0.05, "+", NA)))) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = factor(Variable, levels=c("Time:Response", "Status", "Location"))) %>%
  ggplot(aes(x=R2, y=Variable))+
  geom_bar(stat="identity", color="white", fill="black", width=0.5)+
  # scale_color_manual(values=c("*" = "grey", "***" = "black"))+
  geom_text(aes(x=R2+0.015, label=sig, vjust=ifelse(sig == "+", 1, 0.8)), size=6)+
  scale_x_continuous(limits=c(0,0.2))+
  theme_classic()+
  labs(x=expression(R^2), y=NULL)+
  theme(panel.grid.major.x = element_line(color="grey", linewidth=0.2, linetype=2))


# check status
lsarp.mli.ps.pcoa.resp.plot = ggplot(lsarp.mli.ps.pcoa.resp.df %>% mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
                                aes(x=Axis.1, y=Axis.2))+
  #geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(aes(shape = Location, fill=flare.group), color="white", size=2.5)+
  stat_ellipse(aes(color=flare.group), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = code),size=2)+
  scale_shape_manual(values=c(21,22,23))+
  guides(shape = "none", fill="none")+
  theme_classic()+theme(#legend.position="none",
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))+
  facet_wrap(~Location, scales="free")+
  labs(x=paste("Axis 1:", round(lsarp.mli.ps.pcoa.df$var1, digits=4)[1], "%"),
       y=paste("Axis 2:", round(lsarp.mli.ps.pcoa.df$var2, digits=4)[2], "%"),
       color="Response")
lsarp.mli.ps.pcoa.resp.plot

lsarp.mli.ps.pcoa.plot/
lsarp.mli.ps.pcoa.resp.plot

# :: ASV Interaction -------------------------------------------------------

# Question: are there significant differences before and after RS relative to Placebo

dim(lsarp.mli.ps.mat)

# loop through locations and perform LMER
asv.data.rs.interaction = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat) %in% lsarp.mli.ps.meta$Sample,]
# n = 124

if(rerun==T){
  lsarp.pseudo = min(asv.data.rs.interaction[asv.data.rs.interaction != 0])/2
  
  t1 = Sys.time()
  asv.data.group.rs.interaction.lm = do.call(rbind, lapply(c("combined"), function(location){
    # subset to location (or combined)
    metadata.subset = subset(lsarp.mli.ps.meta, Time %in% c(0,1)) %>% distinct()
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, (Sample) %in% rownames(data.subset))
    # apply filtering per location (or combined)
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = (Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      # // should be defunct
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(asv ~ Group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1
  # calculate padj
  asv.data.group.rs.interaction.lm = asv.data.group.rs.interaction.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.group.rs.interaction.lm, "./host_proteomics/asv.data.group.interaction.lm.Rds")
}
asv.data.group.rs.interaction.lm = readRDS("./host_proteomics/asv.data.group.interaction.lm.Rds")

# number of ASVs present in 20% of samples:
subset(asv.data.group.rs.interaction.lm, !is.na(coef))$feature %>% unique() %>% length()
# 509

# number of sig ASVs (FDR < 0.20):
subset(asv.data.group.rs.interaction.lm, padj < 0.20)$feature %>% unique() %>% length()
# 109

# :: ASV Volcano ----------------------------------------------------------

# append butyrogen (vital) annotation and LCA
asv.data.group.rs.interaction.lm = merge(asv.data.group.rs.interaction.lm %>% mutate(OTU = feature),
                                         lsarp.mli.butyrogens[,c("OTU", "lca.gg2", "butyrogen.ii")],
                                         by="OTU") %>% distinct()
# append butyrogen (original) annotation and LCA
asv.data.group.rs.interaction.lm = merge(asv.data.group.rs.interaction.lm %>% mutate(OTU = feature),
                                         lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")],
                                         by="OTU") %>% distinct()

asv.data.group.rs.interaction.lm = asv.data.group.rs.interaction.lm %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                             paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))

# volcano plots
asv.data.group.rs.interaction.lm.volcano = ggplot(asv.data.group.rs.interaction.lm[,c("coef","butyrogen.i","lca.gg2", "pval", "padj", "feature", "location", "clean.name")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                            aes(x=as.numeric(coef), y=padj))+
  geom_point(shape = 21, aes(fill=butyrogen.i, alpha=ifelse(padj < 0.20, 1, 0.5)))+
  geom_hline(yintercept=0.20, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", clean.name), NA),
                                     color=butyrogen.i),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("butyrogen" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("butyrogen" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Treatment Group")+
  labs(x="Interaction Coefficient",
       y="FDR")
asv.data.group.rs.interaction.lm.volcano
# Not Kircher Butyrogens!

# :: ASV: But I Enrich -------------------------------------------

# need to prevent double-counting proteins
butyrogens.i.but = subset(asv.data.group.rs.interaction.lm, butyrogen.i == "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.other = subset(asv.data.group.rs.interaction.lm, butyrogen.i != "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.but$butyrogen.i = "butyrogen"
butyrogens.i.other$butyrogen.i = "other"
butyrogens.i.but.other = rbind(butyrogens.i.other,
                             butyrogens.i.but)
# any edge cases:
table((butyrogens.i.but$feature)) %>% range()
table(butyrogens.i.other$feature) %>% range()
# no; perfect, continue

subset(butyrogens.i.but.other, padj < 0.2)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
  dplyr::select(direction, butyrogen.i) %>% table() %>%
  fisher.test()
# not enriched for butyrogen ASVs
# p-value = 0.5483
# OR = 0.7386477

mli.butyrogen.i.enrichment.data = subset(butyrogens.i.but.other, padj < 0.20)[,c("feature", "coef", "pval", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.i) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.i = ifelse(butyrogen.i == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq))

mli.butyrogen.i.enrichment.plot = mli.butyrogen.i.enrichment.data %>%
  ggplot()+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.i),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(mli.butyrogen.i.enrichment.data, butyrogen.i == "Butyrogen"),
            aes(
    x=direction, y = .90*100,
    label = paste(round(perc, digits=2)*100, "%", sep="")), color="black",fill="white", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen ASV Enrichment")+
  labs(x="", y="")
mli.butyrogen.i.enrichment.plot


# :: __ but sensitivity -------------------------------------------------------


# goal: see whether:
# A) p val threshold matters
but.enrich.sens.group = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(butyrogens.i.but.other, padj < x)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
    dplyr::select(direction, butyrogen.i) %>% table() %>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
but.enrich.sens.group.plot = but.enrich.sens.group %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.2, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Group ~ Butyrogens")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (adj 𝘱 value)",
       y="Enrichment 𝘱 value")


# :: ASV: Heatmap ---------------------------------------------------------


asv.data.group.rs.interaction.lm.heatmap.data = subset(asv.data.group.rs.interaction.lm, padj < 0.20) %>%
  mutate(butyrogen.i = ifelse(butyrogen.i == "butyrogen", " ", "  ")) %>%
  mutate(lca.gg2 = make.unique(lca.gg2))%>%
  dplyr::select(butyrogen.i, lca.gg2, coef)  %>% distinct() %>%
  mutate(but.rank = rank(coef))

# add butyrogen annotations
asv.data.group.rs.interaction.lm.heatmap.data.taxa = asv.data.group.rs.interaction.lm.heatmap.data %>%
  # subset to sig mito proteins
  subset(butyrogen.i == " ") %>% dplyr::select(lca.gg2, but.rank, coef, butyrogen.i) %>% distinct() %>%
  # give them a position along the x axis
  arrange(coef) %>%
  mutate(rank.pos = rank(coef))%>%
  mutate(rel.pos = min(rank.pos) + (((rank.pos - min(rank.pos)) * ((length(unique(asv.data.group.rs.interaction.lm.heatmap.data$lca.gg2)))-1))/
                                      (max(rank.pos)-1)))
# clean names slightly
asv.data.group.rs.interaction.lm.heatmap.data = asv.data.group.rs.interaction.lm.heatmap.data %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                                   paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))
asv.data.group.rs.interaction.lm.heatmap.data.taxa = asv.data.group.rs.interaction.lm.heatmap.data.taxa %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                             paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))
         
asv.data.group.rs.interaction.lm.heatmap.plot = ggplot(asv.data.group.rs.interaction.lm.heatmap.data,
                                                           aes(x=reorder(lca.gg2, coef), y=1))+
  geom_tile(aes(fill=coef), color="black")+
  geom_point(aes(x=reorder(lca.gg2, coef), y=2, color=butyrogen.i), shape=73, size=5)+
  geom_tile(aes(y=5), fill="white")+
  geom_text(data=asv.data.group.rs.interaction.lm.heatmap.data.taxa,
            aes(x=rel.pos, y=4, label=gsub("\\..*", "", clean.name)), angle=90, size=3, hjust=0, nudge_y=0)+
  geom_text(data=asv.data.group.rs.interaction.lm.heatmap.data.taxa,
            aes(x=rel.pos, y=20, label=" "), angle=90, size=3, hjust=0, nudge_y=0)+
  scale_color_manual(values=c(" " = "red"), na.value=NA)+
  geom_segment(data = asv.data.group.rs.interaction.lm.heatmap.data.taxa,
               aes(x=but.rank, xend=rel.pos, y=2.5, yend= 3.7), linetype=1, linewidth=0.2)+
  #guides(color="none")+
  viridis::scale_fill_viridis()+
  scale_x_discrete(expand = c(0.01, 1))+
  theme_minimal()+
  theme(legend.position="top",
        panel.grid=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())+
  labs(fill="Coefficient", title=NULL, x=NULL, y=NULL, color=" Butyrogen")+
  theme(plot.margin = unit(c(0, 1, 1, 1), "cm"))
asv.data.group.rs.interaction.lm.heatmap.plot



# :: ASV: But II Enrich -------------------------------------------

# need to prevent double-counting proteins
butyrogens.ii.but = subset(asv.data.group.rs.interaction.lm, butyrogen.ii == "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.ii.other = subset(asv.data.group.rs.interaction.lm, butyrogen.ii != "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.ii.but$butyrogen.ii = "butyrogen"
butyrogens.ii.other$butyrogen.ii = "other"
butyrogens.ii.but.other = rbind(butyrogens.ii.other,
                                butyrogens.ii.but)
# any edge cases:
table(butyrogens.ii.but$feature) %>% range()
table(butyrogens.ii.other$feature) %>% range()
# no; perfect, continue

subset(butyrogens.ii.but.other, padj < 0.2)[,c("feature", "coef", "padj", "butyrogen.ii")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(butyrogen.ii = factor(butyrogen.ii, levels=c("other", "butyrogen")))%>%
  dplyr::select(direction, butyrogen.ii) %>% table() %>%
  fisher.test()
# significantly enriched for butyrogen ASVs
# p-value = 0.02189
# OR = 2.878781

mli.butyrogen.ii.enrichment.data = subset(butyrogens.ii.but.other, padj < 0.20)[,c("feature", "coef", "pval", "butyrogen.ii")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.ii) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.ii = ifelse(butyrogen.ii == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq))

mli.butyrogen.ii.enrichment.plot = mli.butyrogen.ii.enrichment.data %>%
  ggplot()+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.ii),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(mli.butyrogen.ii.enrichment.data, butyrogen.ii == "Butyrogen"),
            aes(
              x=direction, y = .90*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="black", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen ASV Enrichment")+
  labs(x="", y="")
mli.butyrogen.ii.enrichment.plot



# ::  -----------------------------------------------------------


# :: ASV Interactions: Responders -----------------------------------------------------------

# goal: conduct analysis using Response/Non-Response paradigm


if(rerun==T){
  t1 = Sys.time()
  asv.data.group.rs.interaction.lm.resp = do.call(rbind, lapply(c("combined"), function(location){
    # subset to location (or combined)
    metadata.subset = subset(lsarp.mli.ps.meta, Group == "RS" & Time %in% c(0, 1)) %>% distinct()
    # add Flare.group
    metadata.subset = merge(metadata.subset,
                            lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(), by="HM")
    
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, (Sample) %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[192]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = (Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      # 20% prevalence filter
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))

        lmer.results = lmerTest::lmer(asv ~ flare.group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 1 min
  # calculate padj
  asv.data.group.rs.interaction.lm.resp = asv.data.group.rs.interaction.lm.resp %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.group.rs.interaction.lm.resp, "./host_proteomics/asv.data.group.interaction.resp.lm.Rds")
}
asv.data.group.rs.interaction.lm.resp = readRDS("./host_proteomics/asv.data.group.interaction.resp.lm.Rds")

# controlling for inflammation and location

# :: ASV Volcano: Responders ----------------------------------------------------------

# append butyrogen annotation and LCA
asv.data.group.rs.interaction.lm.resp = merge(asv.data.group.rs.interaction.lm.resp %>% mutate(OTU = feature),
                                              lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")],
                                         by="OTU")
asv.data.group.rs.interaction.lm.resp = merge(asv.data.group.rs.interaction.lm.resp %>% mutate(OTU = feature),
                                              lsarp.mli.butyrogens[,c("OTU", "lca.gg2", "butyrogen.ii")],
                                              by="OTU")


asv.data.group.rs.interaction.lm.resp = asv.data.group.rs.interaction.lm.resp %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                             paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))


# volcano plots
asv.data.group.rs.interaction.lm.resp.volcano = ggplot(asv.data.group.rs.interaction.lm.resp[,c("coef","lca.gg2", "butyrogen.i", "pval", "padj", "feature", "location", "clean.name")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                            aes(x=as.numeric(coef), y=padj))+
  geom_point(shape = 21, aes(fill=butyrogen.i, alpha=ifelse(padj < 0.20, 1, 0.5)))+
  geom_hline(yintercept=0.20, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", clean.name), NA),
                                     color=butyrogen.i),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("butyrogen" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("butyrogen" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Therapy Response")+
  labs(x="Interaction Coefficient",
       y="FDR")
asv.data.group.rs.interaction.lm.resp.volcano

    
# HM0865 and HM0938 have A in TI

# HM0865, Time 01 (rescope)
# large ulcers in ileum; affecting 50-75% of surface
# aphthous ulcers in transverse and distal colon; affecting < 50% of surface
# total SES-CD = 12 (6 in ileum), 3/3 in transverse and distal colon
# versus 31 (10, 8/5) at baseline
# versus 10 at emergency scope after re-scope

# HM0938, Time 01 (rescope)
# SES = 1 in ileum (and total); but no presence of ulcer; only affected surface
# versus 4 at baseline (in ileum)



# :: ASV: Responders But I Enrich -------------------------------------------

# need to prevent double-counting proteins
butyrogens.i.resp.but = subset(asv.data.group.rs.interaction.lm.resp, butyrogen.i == "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.other = subset(asv.data.group.rs.interaction.lm.resp, butyrogen.i != "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.but$butyrogen.i = "butyrogen"
butyrogens.i.resp.other$butyrogen.i = "other"
butyrogens.i.resp.but.other = rbind(butyrogens.i.resp.other,
                                   butyrogens.i.resp.but)
# any edge cases:
table((butyrogens.i.resp.but$feature)) %>% range()
table((butyrogens.i.resp.other$feature)) %>% range()
# no; perfect, continue

subset(butyrogens.i.resp.but.other, padj < 0.2)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
  dplyr::select(direction, butyrogen.i) %>% table() %>%
  fisher.test()
# sig enriched for butyrogen ASVs
# p-value = 9.65e-05
# OR = 8.135093

mli.resp.butyrogen.i.enrichment.data = subset(butyrogens.i.resp.but.other, padj < 0.2)[,c("feature", "coef", "pval", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.i) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.i = ifelse(butyrogen.i == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq))

mli.resp.butyrogen.i.enrichment.plot = mli.resp.butyrogen.i.enrichment.data %>%
  ggplot(aes(x=direction, y=butyrogen.i, fill=direction))+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.i),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(mli.resp.butyrogen.i.enrichment.data, butyrogen.i == "Butyrogen"),
            aes(
              x=direction, y = .90*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="black",fill="white", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen ASV Enrichment")+
  labs(x="", y="")
mli.resp.butyrogen.i.enrichment.plot
# Original Butyrogens


# :: ASV: Heatmap ---------------------------------------------------------


asv.data.group.rs.interaction.lm.resp.heatmap.data = subset(asv.data.group.rs.interaction.lm.resp %>% distinct(), padj < 0.20) %>%
  mutate(butyrogen.i = ifelse(butyrogen.i == "butyrogen", " ", "  ")) %>%
  mutate(lca.gg2 = make.unique(lca.gg2))%>%
  dplyr::select(butyrogen.i, lca.gg2, coef)  %>% distinct() %>%
  mutate(but.rank = rank(coef))

# add butyrogen annotations
asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa = asv.data.group.rs.interaction.lm.resp.heatmap.data %>%
  # subset to sig mito proteins
  subset(butyrogen.i == " ") %>% dplyr::select(lca.gg2, but.rank, coef, butyrogen.i) %>% distinct() %>%
  # give them a position along the x axis
  arrange(coef) %>%
  mutate(rank.pos = rank(coef))%>%
  mutate(rel.pos = min(rank.pos) + (((rank.pos - min(rank.pos)) * ((length(unique(asv.data.group.rs.interaction.lm.resp.heatmap.data$lca.gg2)))-1))/
                                      (max(rank.pos)-1)))


asv.data.group.rs.interaction.lm.resp.heatmap.data = asv.data.group.rs.interaction.lm.resp.heatmap.data %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                             paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))
asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa = asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa %>%
  mutate(clean.name = ifelse(grepl("s__", lca.gg2), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep=""),
                             paste("(", substr(lca.gg2, 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", lca.gg2))), sep="")))


asv.data.group.rs.interaction.lm.resp.heatmap.plot = ggplot(asv.data.group.rs.interaction.lm.resp.heatmap.data,
                                                       aes(x=reorder(lca.gg2, coef), y=1))+
  geom_tile(aes(fill=coef), color="black")+
  geom_point(aes(x=reorder(lca.gg2, coef), y=2, color=butyrogen.i), shape=73, size=5)+
  geom_tile(aes(y=5), fill="white")+
  geom_text(data=asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa,
            aes(x=rel.pos, y=4, label=gsub("\\..*", "", clean.name)), angle=90, size=3, hjust=0, nudge_y=0)+
  geom_text(data=asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa,
            aes(x=rel.pos, y=15, label=" "), angle=90, size=3, hjust=0, nudge_y=0)+
  scale_color_manual(values=c(" " = "red"), na.value=NA)+
  geom_segment(data = asv.data.group.rs.interaction.lm.resp.heatmap.data.taxa,
               aes(x=but.rank, xend=rel.pos, y=2.5, yend= 3.7), linetype=1, linewidth=0.2)+
  #guides(color="none")+
  viridis::scale_fill_viridis()+
  scale_x_discrete(expand = c(0.01, 1))+
  theme_minimal()+
  theme(legend.position="top",
        panel.grid=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())+
  labs(fill="Coefficient", title=NULL, x=NULL, y=NULL, color=" Butyrogen")+
  theme(plot.margin = unit(c(0, 1, 1, 1), "cm"))
asv.data.group.rs.interaction.lm.resp.heatmap.plot

# :: __ but sensitivity -------------------------------------------------------


# goal: see whether:
# A) p val threshold matters

# A: p val sensitivity
but.enrich.sens.resp = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(asv.data.group.rs.interaction.lm.resp, padj < x)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
    dplyr::select(direction, butyrogen.i) %>% table() %>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
but.enrich.sens.resp.plot = but.enrich.sens.resp %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.20, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ Butyrogens")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

but.enrich.sens.group.plot+
but.enrich.sens.resp.plot

# :: ASV: Responders But II Enrich -------------------------------------------

# need to prevent double-counting proteins
butyrogens.ii.resp.but = subset(asv.data.group.rs.interaction.lm.resp, butyrogen.ii == "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.ii.resp.other = subset(asv.data.group.rs.interaction.lm.resp, butyrogen.ii != "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.ii.resp.but$butyrogen.ii = "butyrogen"
butyrogens.ii.resp.other$butyrogen.ii = "other"
butyrogens.ii.resp.but.other = rbind(butyrogens.ii.resp.other,
                                  butyrogens.ii.resp.but)
# any edge cases:
table(butyrogens.ii.resp.but$feature) %>% range()
table(butyrogens.ii.resp.other$feature) %>% range()
# no; perfect, continue

subset(butyrogens.ii.resp.but.other, padj < 0.2)[,c("feature", "coef", "padj", "butyrogen.ii")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(butyrogen.ii = factor(butyrogen.ii, levels=c("other", "butyrogen")))%>%
  dplyr::select(direction, butyrogen.ii) %>% table() %>%
  fisher.test()
# not enriched for butyrogen ASVs
# p-value = 1
# OR = 1.064232


mli.resp.butyrogen.ii.enrichment.data = subset(butyrogens.ii.resp.but.other, padj < 0.2)[,c("feature", "coef", "pval", "butyrogen.ii")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.ii) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.ii = ifelse(butyrogen.ii == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq))

mli.resp.butyrogen.ii.enrichment.plot = mli.resp.butyrogen.ii.enrichment.data %>%
  ggplot(aes(x=direction, y=butyrogen.ii, fill=direction))+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.ii),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_label(data=subset(mli.resp.butyrogen.ii.enrichment.data, butyrogen.ii == "Butyrogen"),
            aes(
              x=direction, y = .90*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="black",fill="white", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen ASV Enrichment")+
  labs(x="", y="")
mli.resp.butyrogen.ii.enrichment.plot


# :: :: PLACEBO -----------------------------------------------------------

# n = 1 placebo Remitter

if(rerun==T){
  t1 = Sys.time()
  asv.data.placebo.rs.interaction.lm.resp = do.call(rbind, lapply(c("combined"), function(location){
    # subset to location
    metadata.subset = subset(lsarp.mli.ps.meta, Group == "Placebo" & Time %in% c(0, 1)) %>% distinct()
    # add Flare.group
    metadata.subset = merge(metadata.subset,
                            lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(), by="HM")
    
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, (Sample) %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[192]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = (Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      # 20% prevalence filter
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        
        lmer.results = lmerTest::lmer(asv ~ flare.group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 1 min
  # calculate padj
  asv.data.placebo.rs.interaction.lm.resp = asv.data.placebo.rs.interaction.lm.resp %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.placebo.rs.interaction.lm.resp, "./host_proteomics/asv.data.placebo.interaction.resp.lm.Rds")
}
asv.data.placebo.rs.interaction.lm.resp = readRDS("./host_proteomics/asv.data.placebo.interaction.resp.lm.Rds")

# controlling for inflammation and location

# :: :: ASV Volcano: Responders ----------------------------------------------------------

# append butyrogen annotation and LCA
asv.data.placebo.rs.interaction.lm.resp = merge(asv.data.placebo.rs.interaction.lm.resp %>% mutate(OTU = feature),
                                              lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")],
                                              by="OTU")
asv.data.placebo.rs.interaction.lm.resp = merge(asv.data.placebo.rs.interaction.lm.resp %>% mutate(OTU = feature),
                                              lsarp.mli.butyrogens[,c("OTU", "lca.gg2", "butyrogen.ii")],
                                              by="OTU")
# volcano plots
asv.data.placebo.rs.interaction.lm.resp.volcano = ggplot(asv.data.placebo.rs.interaction.lm.resp[,c("coef","lca.gg2", "butyrogen.i", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                       aes(x=as.numeric(coef), y=padj))+
  geom_point(shape = 21, aes(fill=butyrogen.i, alpha=ifelse(padj < 0.20, 1, 0.5)))+
  geom_hline(yintercept=0.20, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", lca.gg2), NA),
                                     color=butyrogen.i),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("butyrogen" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("butyrogen" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Therapy Response")+
  labs(x="Interaction Coefficient",
       y="FDR")
asv.data.placebo.rs.interaction.lm.resp.volcano

# enrich
# need to prevent double-counting proteins
butyrogens.i.resp.but.plac = subset(asv.data.placebo.rs.interaction.lm.resp, butyrogen.i == "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.other.plac = subset(asv.data.placebo.rs.interaction.lm.resp, butyrogen.i != "butyrogen")[,c("feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.but.plac$butyrogen.i = "butyrogen"
butyrogens.i.resp.other.plac$butyrogen.i = "other"
butyrogens.i.resp.but.other.plac = rbind(butyrogens.i.resp.other.plac,
                                    butyrogens.i.resp.but.plac)
# any edge cases:
table(butyrogens.i.resp.but.plac$feature) %>% range()
table(butyrogens.i.resp.other.plac$feature) %>% range()
# no; perfect, continue

subset(butyrogens.i.resp.but.other.plac, padj < 0.2)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
  dplyr::select(direction, butyrogen.i) %>% table() %>%
  fisher.test()
# sig enriched for butyrogen ASVs
# p-value = 0.004548
# OR = 3.579375

mli.resp.butyrogen.i.enrichment.data.plac = subset(butyrogens.i.resp.but.other.plac, padj < 0.2)[,c("feature", "coef", "pval", "butyrogen.i")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.i) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.i = ifelse(butyrogen.i == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq))

mli.resp.butyrogen.i.enrichment.plac.plot = mli.resp.butyrogen.i.enrichment.data.plac %>%
  ggplot(aes(x=direction, y=butyrogen.i, fill=direction))+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.i),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(mli.resp.butyrogen.i.enrichment.data.plac, butyrogen.i == "Butyrogen"),
            aes(
              x=direction, y = .90*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="white", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen ASV Enrichment")+
  labs(x="", y="")
mli.resp.butyrogen.i.enrichment.plac.plot
# Original Butyrogens


# :: :: __ but sensitivity -------------------------------------------------------


# goal: see whether:
# A) p val threshold matters

# A: p val sensitivity
but.enrich.sens.resp.placebo = do.call(rbind, lapply(seq(from=0.01, to=0.25, by=0.01), function(x){
  fisher.result = subset(asv.data.placebo.rs.interaction.lm.resp, padj < x)[,c("feature", "coef", "padj", "butyrogen.i")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    mutate(butyrogen.i = factor(butyrogen.i, levels=c("other", "butyrogen")))%>%
    dplyr::select(direction, butyrogen.i) %>% table() %>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             threshold = x)
}))
but.enrich.sens.resp.placebo.plot = but.enrich.sens.resp.placebo %>%
  ggplot(aes(x=threshold, y=pval))+
  geom_line()+
  geom_hline(yintercept = 0.05, color="red", linetype=2, linewidth=0.2)+
  geom_vline(xintercept = 0.20, color="black", linetype=2, linewidth=0.2)+
  theme_classic()+
  facet_wrap(~"Responders ~ Butyrogens (Placebo)")+
  theme(strip.text = element_text(size=10))+
  labs(x="Threshold for Significance (raw 𝘱 value)",
       y="Enrichment 𝘱 value")

but.enrich.sens.resp.placebo.plot


# >> Threshold Sensitivity ----------------------------------------------------------

# show all sensitivity analyses

(mito.enrich.sens.group.plot+
  mito.enrich.sens.resp.plot)/
(mito2.enrich.sens.group.plot+
  mito2.enrich.sens.resp.plot)/
(but.enrich.sens.group.plot+
  but.enrich.sens.resp.plot)
# For RS: mito are only enriched in responders, not RS treated (across low thresholds)
# For RS: but are only enriched in responders, not RS treated (across low thresholds)

(mito.enrich.sens.resp.placebo.plot)/
     (mito2.enrich.sens.resp.placebo.plot)
# For Placebo: favors significance with GO, but not sig with MitoCarta
but.enrich.sens.resp.placebo.plot
# For Placebo: butyrogens are also enriched in responders

# >>> 3. Correlations ---------------------------------------------------------

## goal: correlate proteins with ASVs

# data:
proteomics.data.clean
proteomics.pseudo = min(proteomics.data.clean[proteomics.data.clean!=0])/2

lsarp.mli.ps.mat
asv.pseudo = min(lsarp.mli.ps.mat[lsarp.mli.ps.mat!=0])/2


# :: subset to sig ---------------------------------------------------------------

# subset to sig proteins (resp/nonresp comparisons)
protein.sig.protein = subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05)$protein
proteomics.data.clean.sig = proteomics.data.clean[,colnames(proteomics.data.clean) %in%
                                                protein.sig.protein]
dim(proteomics.data.clean.sig)
# 265

# subset to sig asv (resp/nonresp comparisons)
asv.sig.asv = subset(asv.data.group.rs.interaction.lm, padj < 0.20)$OTU
asv.data.clean = lsarp.mli.ps.mat[,colnames(lsarp.mli.ps.mat) %in%
                                    asv.sig.asv]
dim(asv.data.clean)
# 109 ASVs

# how many total correlations calculated?
265 * 109
# 28,885


# :: regress out covar ----------------------------------------------------

# first, regress out individual, location, and inflammation status

# proteomics:
proteomics.data.regressed = do.call(cbind, lapply(1:ncol(proteomics.data.clean.sig), function(x){
  print(x)
  data.subset = data.frame(protein = proteomics.data.clean.sig[,x]) %>%
    mutate(sample = rownames(.)) %>%
    tidyr::separate(sample, into=c("HM", "Time", "Location", "Status"), remove=F)%>%
    mutate(code = sample)
  # remove plate effects
  data.subset = merge(data.subset, proteomics.metadata[,c("code", "Plate")], by="code")%>%
    mutate(plate = as.factor(Plate))
  # remove 0
  data.subset.no0 = subset(data.subset, protein != 0)
  #
  lm.model = lmerTest::lmer(log2(protein) ~ 1 + Time + Location + Status + (1|HM) + (1|plate),
                            data.subset.no0)
  # use model to predict residuals of full data (keeping order intact and preserving 0 as NA)
  new.data = data.frame(resid = log2(data.subset$protein) - predict(lm.model, data.subset %>% mutate(protein = log2(protein))))
  new.data[new.data == "-Inf"] = NA
  colnames(new.data) = colnames(proteomics.data.clean.sig)[x]
  rownames(new.data) = rownames(proteomics.data.clean.sig)
  new.data
}))

# confirm good quality (randomly select 9 features)

do.call(rbind, lapply(1:9, function(x){
  set.seed(x)
  x = sample(1:ncol(proteomics.data.clean.sig), size = 1)
  data.subset = data.frame(protein = proteomics.data.clean.sig[,x]) %>%
    mutate(sample = rownames(.)) %>%
    tidyr::separate(sample, into=c("HM", "Time", "Location", "Status"), remove=F)%>%
    mutate(code = sample)
  # remove plate effects
  data.subset = merge(data.subset, proteomics.metadata[,c("code", "Plate")], by="code")%>%
    mutate(plate = as.factor(Plate))
  # remove 0 counts
  data.subset = subset(data.subset, protein != 0)
  # build model
  lm.resid = lmerTest::lmer(log2(protein) ~ 1 + Time + Location + Status + (1|HM) + (1|plate),
                            data.subset)
  # check fit
  data.frame(resid = resid(lm.resid), 
             fitted = fitted(lm.resid),
             # add predict to check that predict will fit data properly
             pred = predict(lm.resid, data.subset %>% mutate(protein = log2(protein))),
             pred.resid = log2(data.subset$protein) - predict(lm.resid, data.subset %>% mutate(protein = log2(protein)),
                                                              allow.new.levels=T),
             # it's the same
             feature = colnames(proteomics.data.clean.sig)[x])
})) %>%
  ggplot(aes(x=pred.resid, y=pred))+
  geom_point(shape=21, fill="black", color="white")+
  theme_classic()+
  facet_wrap(~feature, scales="free")

# asv
asv.data.regressed = do.call(cbind, lapply(1:ncol(asv.data.clean), function(x){
  print(x)
  data.subset = data.frame(asv = asv.data.clean[,x]) %>%
    mutate(sample = rownames(.)) %>%
    tidyr::separate(sample, into=c("HM","ASP", "Location"), sep="-") %>%
    tidyr::separate(HM, into=c("HM", "Time"), sep="\\.") %>%
    mutate(Time = substr(Time, 2,2)) %>%
    mutate(sample = paste(HM, Time, Location, sep="_")) %>%
    mutate(HM = as.factor(HM))
  # remove 0
  data.subset.no0 = subset(data.subset, asv != 0)
  #
  lm.model = lmerTest::lmer(log2(asv) ~ 1 + Time + Location + (1|HM),
                            data.subset.no0)
  # use model to predict residuals of full data (keeping order intact and preserving 0 as NA)
  new.data = data.frame(resid = log2(data.subset$asv) - predict(lm.model, data.subset %>% mutate(asv = log2(asv)),
                                                                allow.new.levels=T))
  new.data[new.data == "-Inf"] = NA
  colnames(new.data) = colnames(asv.data.clean)[x]
  rownames(new.data) = data.subset$sample
  new.data
}))

# confirm good quality (randomly select 9 features)

do.call(rbind, lapply(1:9, function(x){
  set.seed(x)
  x = sample(1:ncol(asv.data.clean), size = 1)
  data.subset = data.frame(asv = asv.data.clean[,x]) %>%
    mutate(sample = rownames(.)) %>%
    tidyr::separate(sample, into=c("HM","ASP", "Location"), sep="-") %>%
    tidyr::separate(HM, into=c("HM", "Time"), sep="\\.") %>%
    mutate(Time = substr(Time, 2,2)) %>%
    mutate(sample = paste(HM, Time, Location, sep="_"))
  # remove 0 counts
  data.subset = subset(data.subset, asv != 0)
  # build model
  lm.resid = (lmerTest::lmer(log2(asv) ~ 1 + Time + Location + (1|HM),
                                  data.subset))
  # check fit
  resid.df = data.frame(resid = resid(lm.resid), 
             fitted = fitted(lm.resid),
             feature = colnames(proteomics.data.clean.sig)[x])
  # remove 0 values; replace with NA
 # resid.df$resid[asv.data.clean[,x] == (asv.pseudo*2)] = NA
  # collect data
  resid.df
})) %>%
  ggplot(aes(x=resid, y=fitted))+
  geom_point(shape=21, fill="black", color="white")+
  theme_classic()+
  facet_wrap(~feature, scales="free")


# now take mean values per HM, Time, Location
proteomics.data.regressed.mean = proteomics.data.regressed %>% as.matrix()%>%
  reshape2::melt() %>%
  tidyr::separate(col=Var1, into=c("HM", "Time", "Location", "Status")) %>%
  group_by(HM, Time, Location, Var2) %>%
  mutate(mean.value = mean(na.omit(value))) %>%
  dplyr::select(HM, Time, Location, Var2, mean.value) %>% distinct() %>% as.data.frame() %>%
  mutate(sample = paste(HM, Time, Location, sep="_")) %>%
  dplyr::select(sample, Var2, mean.value)

asv.data.regressed.mean = asv.data.regressed %>% as.matrix()%>%
  reshape2::melt() %>%
  tidyr::separate(col=Var1, into=c("HM", "Time", "Location")) %>%
  group_by(HM, Time, Location, Var2) %>%
  mutate(mean.value = mean(na.omit(value))) %>%
  dplyr::select(HM, Time, Location, Var2, mean.value) %>% distinct() %>% as.data.frame() %>%
  mutate(sample = paste(HM, Time, Location, sep="_")) %>%
  dplyr::select(sample, Var2, mean.value)

# merge datasets
proteomics.data.regressed.mean.mat = proteomics.data.regressed.mean %>%
  reshape2::acast(sample ~ Var2, value.var="mean.value") %>%
  as.data.frame() %>%
  mutate(sample = rownames(.))
asv.data.regressed.mean.mat = asv.data.regressed.mean %>%
  reshape2::acast(sample ~ Var2, value.var="mean.value") %>%
  as.data.frame() %>%
  mutate(sample = rownames(.))

proteomics.asv.data.regressed.merge = 
  merge(proteomics.data.regressed.mean.mat,
        asv.data.regressed.mean.mat, by="sample") 

rownames(proteomics.asv.data.regressed.merge) = proteomics.asv.data.regressed.merge$sample
proteomics.asv.data.regressed.merge$sample = NULL


grepl("ASV", colnames(proteomics.asv.data.regressed.merge)) %>% sum()
(!grepl("ASV", colnames(proteomics.asv.data.regressed.merge))) %>% sum()


# :: subset to RS responders ----------------------------------------------

rs.responders = subset(proteomics.metadata, Group == "RS" & flare.group == "Remit")[,c("HM")] %>% unique()

proteomics.asv.data.regressed.merge.rs = proteomics.asv.data.regressed.merge[grepl(paste(rs.responders, collapse="|"), rownames(proteomics.asv.data.regressed.merge)),]
dim(proteomics.asv.data.regressed.merge.rs)
# 27 samples x 444 features

# subset to present in >20%
proteomics.asv.data.regressed.merge.rs.pa = proteomics.asv.data.regressed.merge.rs
proteomics.asv.data.regressed.merge.rs.pa[is.na(proteomics.asv.data.regressed.merge.rs.pa)] = 0
proteomics.asv.data.regressed.merge.rs.pa[(proteomics.asv.data.regressed.merge.rs.pa)!=0] = 1
proteomics.asv.data.regressed.merge.rs = proteomics.asv.data.regressed.merge.rs[,colSums(proteomics.asv.data.regressed.merge.rs.pa) >
                                                                                  (nrow(proteomics.asv.data.regressed.merge.rs.pa)*0.2)]
dim(proteomics.asv.data.regressed.merge.rs)
# 27 samples x 427 samples

# however, some end up having R = 1
plot(proteomics.asv.data.regressed.merge.rs$ASV00105,
     proteomics.asv.data.regressed.merge.rs$CIR1A)
# which is fine for pearson

# :: correlate ------------------------------------------------------------

# calculate spearman between ASVs and proteins
proteomics.asv.data.regressed.cor = Hmisc::rcorr(as.matrix(proteomics.asv.data.regressed.merge.rs), type="spearman")
# subset to important comparisons
proteomics.asv.data.regressed.cor = cbind(reshape2::melt(proteomics.asv.data.regressed.cor$r),
                                          data.frame(p= reshape2::melt(proteomics.asv.data.regressed.cor$P)[,3])) %>%
  subset(Var1 != Var2) %>%
  arrange(-abs(value)) %>%
  mutate(type1 = ifelse(Var1 %in% protein.sig.protein, "protein", "ASV")) %>%
  mutate(type2 = ifelse(Var2 %in% protein.sig.protein, "protein", "ASV")) #%>%
  #subset(type1 != type2) %>%
  #subset(type1 == "ASV")

# add taxa
proteomics.asv.data.regressed.cor$OTU = proteomics.asv.data.regressed.cor$Var1
proteomics.asv.data.regressed.cor = merge(proteomics.asv.data.regressed.cor,
                                          asv.data.group.rs.interaction.lm[,c("OTU", "lca.gg2","lca.gg138", "butyrogen.i")],
                                          by="OTU")
proteomics.asv.data.regressed.cor$Var1 = NULL
# add mito (MitoCarta)
proteomics.asv.data.regressed.cor$protein = proteomics.asv.data.regressed.cor$Var2
proteomics.asv.data.regressed.cor = merge(proteomics.asv.data.regressed.cor,
                                          proteomics.data.group.interaction.resp.lm.mito[,c("protein", "mito")],
                                          by="protein")
proteomics.asv.data.regressed.cor$Var2 = NULL

# subset to butyrogens & mitochondrial proteins
proteomics.asv.data.regressed.cor.sig = proteomics.asv.data.regressed.cor %>%
  subset(butyrogen.i == "butyrogen" & mito == "mito") %>%
  distinct() %>%
  mutate(padj = p.adjust(p, method="BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "***", ifelse(padj < 0.20, "*", ""))) %>%
  arrange(p)

proteomics.asv.data.regressed.cor.sig %>%
  arrange(padj)

# among strongest correlation = Agathobacter rectalis == ACADM

# :: plot -----------------------------------------------------------------

library("network")
cor_data <- proteomics.asv.data.regressed.cor.sig %>%
  filter(p < 0.05 & value > 0) %>% # or p < 0.01 etc.
  select(protein, OTU, value, lca.gg2)

# Create the bipartite network object -----------------------------------
net <- network::network(as.matrix(cor_data[,1:2]), matrix.type = "edgelist", directed = FALSE)

# Add vertex attributes -------------------------------------------------
set.vertex.attribute(net, "type",      ifelse(network.vertex.names(net) %in% cor_data$protein, "Protein", "Taxon"))
set.vertex.attribute(net, "label",     ifelse(network.vertex.names(net) %in% cor_data$protein,
                                              network.vertex.names(net),
                                              cor_data$lca.gg2[match(network.vertex.names(net), cor_data$OTU)]))
set.edge.attribute(net, "weight", abs(cor_data$value))
set.edge.attribute(net, "sign",   ifelse(sign(cor_data$value) == 1, "Positive", "Negative"))  # +1 or –1

# Plot with ggnetwork ---------------------------------------------------
library("ggnetwork")

set.seed(2)
network.plot = ggnetwork(net)
# fix ames
network.plot$label = ifelse(network.plot$type == "Taxon", ifelse(grepl("s__", (network.plot$label)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (network.plot$label)))), sep=""),
                                   paste("(", substr((network.plot$label), 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (network.plot$label)))), sep="")),
                            network.plot$label)

asv.protein.network.plot = ggplot(network.plot, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(size = weight, color = sign),
             curvature = 0, alpha = 0.7,
             show.legend = FALSE) +
  scale_size_continuous(range = c(0.3, 2.5)) +
  geom_nodes(aes(color = type, shape = type),
             size = 5) +
  scale_color_manual(values = c("Protein" = "lightblue", "Taxon" = "plum",
                                "Negative" = "#1f78b4", "Positive" = "grey")) +  # green vs orange #e31a1c
  scale_shape_manual(values = c("Protein" = 19, "Taxon" = 15)) +  # circle vs square
  geom_nodetext_repel(data = subset(network.plot, type == "Protein"),
                       aes(label = label, x=x, y=y), size = 4, fontface = "bold", color = "black") +
  geom_nodetext_repel(data = subset(network.plot, type == "Taxon"),
                       aes(label = label, x=x,y=y), size = 3, color = "black", max.overlaps=5) +
  theme_blank() +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  labs(title = "Mucosal Proteome ↔ Microbiome Correlation Network",
       subtitle = "Spearman | 𝘱 < 0.05, edges sized by |ρ|") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

asv.protein.network.plot

proteomics.asv.data.regressed.cor.sig %>%
  subset(p < 0.05) %>%
  subset(grepl("bromii", lca.gg2))


proteomics.asv.data.regressed.cor.sig %>% 
  # add oxphos label
  mutate(oxphos = ifelse(protein %in% unique(subset(proteomics.data.group.interaction.lm.mito, grepl("OXPHOS", MitoPathways.Hierarchy))$protein), "oxphos", ""))%>%
  arrange(-value) %>% 
  head(n=20) %>%
  subset(oxphos == "oxphos")
 
# plot select correlations
mpx.16s.cor.plot.1 = proteomics.asv.data.regressed.merge.rs %>%
  dplyr::select(unique(subset(network.plot, label == "(s) Anaerostipes_hadrus")$vertex.names),
                SURF1) %>%
  ggplot(aes(x=ASV00296,
             y=SURF1)) +
  geom_point(fill="black", color="white", shape=21, size=3)+
  geom_smooth(method="lm", fill="black", color="white")+
  ggpubr::stat_cor(method="spearman", size=3)+
  theme_classic()+
  labs(x="Anaerostipes hadrus",
       y="SURF1")
  
mpx.16s.cor.plot.2 = proteomics.asv.data.regressed.merge.rs %>%
  dplyr::select(unique(subset(network.plot, label == "(s) Ruminococcus_E_bromii_B")$vertex.names),
                DHX30) %>%
  ggplot(aes(x=ASV00105,
             y=DHX30)) +
  geom_point(fill="black", color="white", shape=21, size=3)+
  geom_smooth(method="lm", fill="black", color="white")+
  ggpubr::stat_cor(method="spearman", size=3)+
  theme_classic()+
  labs(x="Ruminococcus bromii",
       y="DHX30")

mpx.16s.cor.plot.3 = proteomics.asv.data.regressed.merge.rs %>%
  dplyr::select(unique(subset(network.plot, label == "(g) Faecalibacterium")$vertex.names),
                SDHB) %>%
  ggplot(aes(x=ASV00081,
             y=SDHB)) +
  geom_point(fill="black", color="white", shape=21, size=3)+
  geom_smooth(method="lm", fill="black", color="white", size=1)+
  ggpubr::stat_cor(method="spearman", size=3)+
  theme_classic()+
  labs(x="Faecalibacterium",
       y="SDHB")



# >>> FIGURES  -----------------------------------------------------------------

  # Host-Proteome PCA
library("patchwork")
((proteomics.clean.pca.plot+
   proteomics.clean.pca.permanova.r2.plot+
   patchwork::plot_layout(widths=c(2,1)))/
  (proteomics.clean.pca.resp.plot+
     proteomics.clean.pca.permanova.resp.r2.plot+
     patchwork::plot_layout(widths=c(2,1)))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_3_proteome_pca.pdf",
         width=14, height=6,device = cairo_pdf)

# ASV Bacteriome PCoA
((lsarp.mli.ps.pcoa.plot+
    lsarp.mli.ps.permanova.r2.plot+
    patchwork::plot_layout(widths=c(2,1)))/
    (lsarp.mli.ps.pcoa.resp.plot+
       lsarp.mli.ps.resp.permanova.r2.plot+
       patchwork::plot_layout(widths=c(2,1)))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_3_bacteriome_pca.pdf",
         width=14, height=6,device = cairo_pdf)

# new version of PCoAs
pca1 = (lsarp.mli.ps.pcoa.plot+
          lsarp.mli.ps.permanova.r2.plot+
          patchwork::plot_layout(widths=c(3,1)))
pca2 = (lsarp.mli.ps.pcoa.resp.plot+
          lsarp.mli.ps.resp.permanova.r2.plot+
          patchwork::plot_layout(widths=c(3,1)))
pca3 = (proteomics.clean.pca.plot+
          proteomics.clean.pca.permanova.r2.plot+
          patchwork::plot_layout(widths=c(3,1)))
pca4 = (proteomics.clean.pca.resp.plot+
          proteomics.clean.pca.permanova.resp.r2.plot+
          patchwork::plot_layout(widths=c(3,1)))

(pca1/pca2/pca3/pca4) %>%
  ggsave(filename="./lsarp_plots/2026_03_03_lsarp_3_pca.pdf",
         width=14, height=12,device = cairo_pdf)

# Volcanos + Enrichment (Original Butyrogens)
proteomics.heatmaps = (proteomics.data.group.interaction.lm.mito.heatmap + proteomics.data.group.interaction.lm.mito.resp.heatmap + patchwork::plot_layout(guides = "collect") & 
    viridis::scale_fill_viridis(limits = range(c(proteomics.data.group.interaction.lm.mito.heatmap$data$coef, proteomics.resp.mito.other.heatmap$data$coef)))) &
  theme(legend.position = "top")
butyrogen.heatmaps = (asv.data.group.rs.interaction.lm.heatmap.plot + asv.data.group.rs.interaction.lm.resp.heatmap.plot + patchwork::plot_layout(guides = "collect") & 
                         viridis::scale_fill_viridis(limits = range(c(asv.data.group.rs.interaction.lm.heatmap.plot$data$coef, asv.data.group.rs.interaction.lm.resp.heatmap.plot$data$coef)))) &
  theme(legend.position = "top")

proteomics.volcano.enrichments = (
  proteomics.data.group.interaction.lm.mito.volcano+
    proteomics.data.group.interaction.lm.mito.plot+
  proteomics.data.group.interaction.resp.lm.mito.volcano+
    proteomics.mito.gene.enrichment.plot)+
  patchwork::plot_layout(nrow=1, widths=c(1.5,1,1.5,1))

butyrogen.volcano.enrichments = (
  asv.data.group.rs.interaction.lm.volcano+
    mli.butyrogen.i.enrichment.plot+
  asv.data.group.rs.interaction.lm.resp.volcano+
    mli.resp.butyrogen.i.enrichment.plot)+
    patchwork::plot_layout(nrow=1, widths=c(1.5,1,1.5,1))

((
  butyrogen.heatmaps / butyrogen.volcano.enrichments /
    proteomics.heatmaps / proteomics.volcano.enrichments)+
  patchwork::plot_layout(heights=c(2, 1.7, 1, 1.5))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_3_volcanos_flagged_i.pdf",
         width=20, height=20,device = cairo_pdf)

# Volcanos + Enrichment (Kircher Butyrogens)
(proteomics.data.group.interaction.lm.mito.plot+proteomics.data.group.interaction.lm.mito.volcano+
    proteomics.data.group.interaction.resp.lm.mito.volcano+proteomics.mito.gene.enrichment.plot+
    mli.butyrogen.ii.enrichment.plot+asv.data.group.rs.interaction.lm.volcano+
    asv.data.group.rs.interaction.lm.resp.volcano+mli.resp.butyrogen.ii.enrichment.plot+
    patchwork::plot_layout(nrow=2))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_3_volcanos_flagged_ii.pdf",
         width=20, height=8,device = cairo_pdf)
  
(asv.protein.network.plot/
  (mpx.16s.cor.plot.1+
     mpx.16s.cor.plot.2+
     mpx.16s.cor.plot.3+patchwork::plot_layout(nrow=1))+
    patchwork::plot_layout(nrow=2, heights=c(3,1))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_3_network.pdf",
         width=8, height=9,device = cairo_pdf)


# correlations

((but.enrich.sens.group.plot+
    but.enrich.sens.resp.plot)/
    (mito.enrich.sens.group.plot+
    mito.enrich.sens.resp.plot)/
  (mito2.enrich.sens.group.plot+
     mito2.enrich.sens.resp.plot)) %>%
  ggsave(filename="./lsarp_plots/2026_03_04_lsarp_supp_enrich_sens.pdf",
         width=7, height=8, device = cairo_pdf)

# :: ----------------------------------------------------------------------


# >>> Tables --------------------------------------------------------------

# group mpx
table1 = proteomics.data.group.interaction.lm.mito %>%
  dplyr::select(protein, gene, coef, pval, padj, MitoPathway, MitoPathways.Hierarchy, mito, mito2) %>%
  mutate(mito.mitocarta = mito,
         mito.go = mito2) %>%
  dplyr::select(-mito, -mito2)

# group asv
table2 = asv.data.group.rs.interaction.lm %>%
  dplyr::select(OTU, lca.gg138, lca.gg2, coef, pval, padj, butyrogen.i) %>%
  mutate(butyrogen = butyrogen.i) %>%
  dplyr::select(-butyrogen.i) %>% distinct()

# resp mpx (treatment)
table3 = proteomics.data.group.interaction.resp.lm.mito %>%
  dplyr::select(protein, gene, coef, pval, padj, MitoPathway, MitoPathways.Hierarchy, mito, mito2) %>%
  mutate(mito.mitocarta = mito,
         mito.go = mito2) %>%
  dplyr::select(-mito, -mito2)

# resp asv (treatment)
table4 = asv.data.group.rs.interaction.lm.resp %>%
  dplyr::select(OTU, lca.gg138, lca.gg2, coef, pval, padj, butyrogen.i) %>%
  mutate(butyrogen = butyrogen.i) %>%
  dplyr::select(-butyrogen.i) %>% distinct()

# correlations (butyrogens sig correlated with mitochondrial proteins)
table5 = proteomics.asv.data.regressed.cor.sig %>%
  dplyr::select(-butyrogen.i, -mito, -type1, -type2, -sig)


openxlsx::write.xlsx(
  list(
    `st 30 mli_asv` = table2,
    `st 31 mli_asv_resp` = table4,
    `st 32 mli_mpx` = table1,
    `st 33 mli_mpx_resp` = table3,
    `st 34 mli_cor` = table5
  ),
  
  "lsarp_paper/lsarp_heatmap_source_mli.xlsx")


# >>>  SENSITIVITY ----------------------------------------------------------

# Plan
# ✓ Combined proteomics with GO instead of MitoCarta
# ✓ Enrichment, p value thresholds
# ✓ Proteomics per location
# ✓ Proteomics reduced, combined
# ✓ Proteomics reduced, per location
# ✓ MLI per location

# :: Sens: GO -------------------------------------------------------------
# volcano + enrichment
proteomics.data.group.interaction.lm.mito2.volcano = ggplot(proteomics.data.group.interaction.lm.mito[,c("coef", "pval", "padj", "protein", "location", "mito2")] %>% distinct() %>% 
                                                             mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                           aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito2, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito2),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Treatment Group")+
  labs(x="Interaction Coefficient",
       y="p value")
proteomics.data.group.interaction.lm.mito2.volcano


proteomics.data.group.interaction.lm.mito2.data = subset(proteomics.data.group.interaction.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito2) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito2 == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

proteomics.data.group.interaction.lm.mito2.plot = proteomics.data.group.interaction.lm.mito2.data%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(proteomics.data.group.interaction.lm.mito2.data, mito == "Mitochondrial"),
            aes(
              x=direction, y = .97*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
proteomics.data.group.interaction.lm.mito2.plot


# volcano plots
proteomics.data.group.interaction.resp.lm.mito2.volcano = ggplot(proteomics.data.group.interaction.resp.lm.mito[,c("coef", "pval", "padj", "protein", "location", "mito2")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                                aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito2, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito2),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Clinical Response")+
  labs(x="Interaction Coefficient",
       y="p value")
proteomics.data.group.interaction.resp.lm.mito2.volcano


proteomics.mito2.gene.enrichment.data = subset(proteomics.data.group.interaction.resp.lm.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito2) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito2 == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

proteomics.mito2.gene.enrichment.plot =proteomics.mito2.gene.enrichment.data%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(proteomics.mito2.gene.enrichment.data, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  geom_text(data=subset(proteomics.mito2.gene.enrichment.data, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .9*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=5)+
  
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
proteomics.mito2.gene.enrichment.plot


# :: Sens: Thresholds  ----------------------------------------------------

(mito.enrich.sens.group.plot/
   mito.enrich.sens.resp.plot)|
  (mito2.enrich.sens.group.plot/
     mito2.enrich.sens.resp.plot)|
  (but.enrich.sens.group.plot/
     but.enrich.sens.resp.plot)|
# For RS: mito are only enriched in responders, not RS treated (across low thresholds)
# For RS: but are only enriched in responders, not RS treated (across low thresholds)

(mito.enrich.sens.resp.placebo.plot/
  mito2.enrich.sens.resp.placebo.plot)|
# For Placebo: favors significance with GO, but not sig with MitoCarta
(but.enrich.sens.resp.placebo.plot/
   patchwork::plot_spacer())
# For Placebo: butyrogens are also enriched in responders



(mito.enrich.sens.group.plot|
  mito2.enrich.sens.group.plot|
    mito.enrich.sens.resp.plot|
     mito2.enrich.sens.resp.plot)

# :: Sens: MPX - Location ------------------------------------------

t1 = Sys.time()
lsarp.proteomics.location.mito = do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
  print(paste(l))
  
  metadata.subset.2 = subset(proteomics.metadata, Location == l)
  
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset.2$code,]
  
  # loop through proteins
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset.2$code,] %>% data.frame()
  metadata.subset.2 = subset(metadata.subset.2, code %in% rownames(data.subset))
  
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset.3 = merge(metadata.subset.2,
                              data.subset[,colnames(data.subset) %in% c("code", protein)],
                              by="code")
    colnames(metadata.subset.3)[colnames(metadata.subset.3) == protein] = "protein"
    
    # 10% of samples
    if(sum(metadata.subset.3$protein != 0) <= nrow(metadata.subset.3)*0.10){
      data.frame(
        location = l,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset.3$protein = scale(log2(metadata.subset.3$protein+proteomics.pseudocount)) 
      metadata.subset.3$Status = factor(metadata.subset.3$Status, levels=c("N", "A"))
      metadata.subset.3$Location = factor(metadata.subset.3$Location, levels=c("TI", "PC", "DC"))
      metadata.subset.3$Time = factor(metadata.subset.3$Time, levels=c("0", "1"))
      # metadata.subset.3$Plate = as.factor(metadata.subset.3$Plate)
      
      lmer.results = lmerTest::lmer(protein ~ Group*Time + Status + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      
      # extract data
      data.frame(
        location = l,
        feature = protein,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))
}))
t2 = Sys.time()
t2 - t1 # 8 min

lsarp.proteomics.location.mito = lsarp.proteomics.location.mito %>%
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# enrichment
lsarp.proteomics.location.mito.enrich = do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
  data.subset = subset(lsarp.proteomics.location.mito, location == x)
  fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             location = x)
}))
# TI is significantly depleted
# but not corroborated by GO

# volcano
lsarp.proteomics.location.mito.volcano = ggplot(lsarp.proteomics.location.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                                 aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~location)+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.location.mito.volcano

# enrichment plots
lsarp.proteomics.location.mito.enrich = subset(lsarp.proteomics.location.mito, pval < 0.05)[,c("location", "protein", "coef", "pval", "mito")] %>% distinct() %>%
  group_by(location) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(location, direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(location, direction) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.location.mito.enrich.plot =lsarp.proteomics.location.mito.enrich%>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC"))) %>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.location.mito.enrich, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  geom_text(data=subset(lsarp.proteomics.location.mito.enrich, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~location)+
  labs(x="", y="")
lsarp.proteomics.location.mito.enrich.plot
# no differences


# :: Sens: MPX - Location Resp --------------------------------------------

t1 = Sys.time()
lsarp.proteomics.location.mito.resp = 
  do.call(rbind, lapply(c("RS", "Placebo"), function(y){
    do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
      print(paste(y, " ", l))
      
      metadata.subset = subset(proteomics.metadata, Location == l & Group == y)
      
      data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
      
      # loop through proteins
      # keep intersecting samples
      data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
      metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
      
      # apply filtering per location
      proteins.to.keep = data.subset
      proteins.to.keep[proteins.to.keep != 0] <- 1
      proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
      
      data.subset = data.subset[,proteins.to.keep]
      
      # apply variance filter per location
      proteins.to.keep.variance = data.frame(sd = apply((log2(data.subset+proteomics.pseudocount)), 2, sd)) 
      # select top 10%
      proteins.to.keep.variance = proteins.to.keep.variance %>%
        arrange(sd) %>%
        slice_max(order_by=sd, n=round(nrow(proteins.to.keep.variance) *0.10, digits=0))
      
      proteins.to.keep = rownames(proteins.to.keep.variance)
      
      # loop through lmer
      do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
        print(paste0(y, " ", l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
        # protein = proteins.to.keep[1]
        # add protein abundance data by merging
        data.subset$code = rownames(data.subset)
        
        metadata.subset.2 = merge(metadata.subset,
                                  data.subset[,colnames(data.subset) %in% c("code", protein)],
                                  by="code")
        colnames(metadata.subset.2)[colnames(metadata.subset.2) == protein] = "protein"
        
        
        metadata.subset.2$protein = scale(log2(metadata.subset.2$protein+proteomics.pseudocount)) 
        metadata.subset.2$Status = factor(metadata.subset.2$Status, levels=c("N", "A"))
        metadata.subset.2$Location = factor(metadata.subset.2$Location, levels=c("TI", "PC", "DC"))
        metadata.subset.2$Time = factor(metadata.subset.2$Time, levels=c("0", "1"))
        metadata.subset.2$Plate = as.factor(metadata.subset.2$Plate)
        metadata.subset.2$flare.group = factor(metadata.subset.2$flare.group, levels=c("Relapse", "Remit"))
        
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Status + (1|study_id), 
                                      metadata.subset.2) %>% summary() %>% coef() %>% data.frame()
        
        # extract data
        data.frame(
          location = l,
          feature = protein,
          Group = y,
          nfeatures = length(proteins.to.keep),
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
        
      }))
    }))
  }))
t2 = Sys.time()
t2 - t1 # 2 min

lsarp.proteomics.location.mito.resp = lsarp.proteomics.location.mito.resp %>%
  group_by(location, Group) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# enrichment
lsarp.proteomics.location.mito.resp.enrich.fisher = 
  do.call(rbind, lapply(c("Placebo", "RS"), function(y){
    do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
      data.subset = subset(lsarp.proteomics.location.mito.resp, location == x & Group == y)
      if(nrow(subset(data.subset, pval < 0.05 & mito == "mito"))==0){
        data.frame(pval = NA,
                   or = NA,
                   Group = y,
                   location = x)
      }else{
      fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
        mutate(direction = sign(coef)+2) %>%
        dplyr::select(direction, mito) %>% table()%>%
        fisher.test()
      data.frame(pval = fisher.result$p.value,
                 or = fisher.result$estimate,
                 Group = y,
                 location = x)
      }
    }))
  }))
# not enough values to assess RS-PC and RS-DC

# volcano
lsarp.proteomics.location.mito.resp.volcano = ggplot(lsarp.proteomics.location.mito.resp[,c("coef", "pval", "padj", "protein", "location", "mito", "Group", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                                  aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_grid(Group~location, scales="free")+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.location.mito.resp.volcano

# enrich

# enrichment plots
lsarp.proteomics.location.mito.resp.enrich = subset(lsarp.proteomics.location.mito.resp, pval < 0.05)[,c("location", "protein", "coef", "pval", "mito", "location", "Group")] %>% distinct() %>%
  group_by(location, Group) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(location, direction, mito, Group) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(location, direction, Group) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.location.mito.resp.enrich.plot =lsarp.proteomics.location.mito.resp.enrich%>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC")))%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.location.mito.resp.enrich, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  geom_text(data=subset(lsarp.proteomics.location.mito.resp.enrich, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_grid(Group~location)+
  labs(x="", y="")
lsarp.proteomics.location.mito.resp.enrich.plot
# more depletions, but due to low numbers?


# :: Sens: MPX - Reduced --------------------------------------------------

# use only MOST INFLAMED out of pairs of samples

t1 = Sys.time()
lsarp.proteomics.combined.red.mito = do.call(rbind, lapply(c("combined"), function(location){
  # subset to location
  metadata.subset = subset(proteomics.subset, Time %in% c(0,1))
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
  metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset = merge(metadata.subset,
                            data.subset[,colnames(data.subset) %in% c("code", protein)])
    colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
    
    if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
      data.frame(
        location = location,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
      metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
      metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
      metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
      #metadata.subset$Plate = as.factor(metadata.subset$Plate)
      
      lmer.results = lmerTest::lmer(protein ~ Group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
      # extract data
      data.frame(
        location = location,
        protein = protein,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))}))
t2 = Sys.time()
t2 - t1 # 3 min
# calculate padj
lsarp.proteomics.combined.red.mito = lsarp.proteomics.combined.red.mito %>% 
  mutate(padj = p.adjust(pval, method="BH")) %>% 
  arrange(pval)
# Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  

# add mito linker
colnames(lsarp.proteomics.combined.red.mito)[colnames(lsarp.proteomics.combined.red.mito) == "protein"] = "protein.x"
lsarp.proteomics.combined.red.mito = merge(lsarp.proteomics.combined.red.mito,
                                             proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
lsarp.proteomics.combined.red.mito = merge(lsarp.proteomics.combined.red.mito,
                                                  mitocarta.pathways.expanded,
                                                  by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# arrange by significance
lsarp.proteomics.combined.red.mito %>% arrange(padj)

# volcano plots
lsarp.proteomics.combined.red.mito.volcano = ggplot(lsarp.proteomics.combined.red.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% 
                                                             mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                           aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Treatment Group")+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.combined.red.mito.volcano

subset(lsarp.proteomics.combined.red.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table()%>%
  fisher.test()

subset(lsarp.proteomics.combined.red.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito2) %>% table()%>%
  fisher.test()
# neither are sig

lsarp.proteomics.combined.red.mito.data = subset(lsarp.proteomics.combined.red.mito, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.combined.red.mito.enrich.plot = lsarp.proteomics.combined.red.mito.data %>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.combined.red.mito.data, mito == "Mitochondrial"),
            aes(
              x=direction, y = .97*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
lsarp.proteomics.combined.red.mito.enrich.plot

# :: Sens: MPX - Reduced Resp --------------------------------------------------

lsarp.proteomics.combined.red.mito.resp = do.call(rbind, lapply("combined", function(location){
  # subset to location & RS group
  metadata.subset = subset(proteomics.subset, Group == "RS" & Time %in% c(0,1))
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
  metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset = merge(metadata.subset,
                            data.subset[,colnames(data.subset) %in% c("code", protein)])
    colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
    
    # must be detected in more than 3 samples to enter linear model (note: this doesn't actually activate)
    if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
      data.frame(
        location = location,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
      metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
      metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
      metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
      metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
      #metadata.subset$Plate = as.factor(metadata.subset$Plate)
      
      lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Location + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
      # extract data
      data.frame(
        location = location,
        protein = protein,
        coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
    }
  }))}))
t2 = Sys.time()
t2 - t1 # 7 min
# calculate padj
lsarp.proteomics.combined.red.mito.resp = lsarp.proteomics.combined.red.mito.resp %>% 
  mutate(padj = p.adjust(pval, method="BH")) %>% 
  arrange(pval)
# Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  

# add mito linker
colnames(lsarp.proteomics.combined.red.mito.resp)[colnames(lsarp.proteomics.combined.red.mito.resp) == "protein"] = "protein.x"
lsarp.proteomics.combined.red.mito.resp = merge(lsarp.proteomics.combined.red.mito.resp,
                                                  proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
lsarp.proteomics.combined.red.mito.resp = merge(lsarp.proteomics.combined.red.mito.resp,
                                                       mitocarta.pathways.expanded,
                                                       by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""),
         mito2 = ifelse(gene %in% go_mitochondria_annotations$hgnc_symbol, "mito", ""))
# arrange by significance
lsarp.proteomics.combined.red.mito.resp %>% arrange(pval)

# if we had only assessed mitochondrial proteins:
lsarp.proteomics.combined.red.mito.resp %>%
  subset(mito == "mito") %>%
  dplyr::select(protein, coef, pval) %>%
  distinct() %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)
# still not passing FDR < 0.20

# volcano plots
lsarp.proteomics.combined.red.mito.resp.volcano = ggplot(lsarp.proteomics.combined.red.mito.resp[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                                aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Associations with Clinical Response")+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.combined.red.mito.resp.volcano

# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito.2 = subset(lsarp.proteomics.combined.red.mito.resp, mito == "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.other.2 = subset(lsarp.proteomics.combined.red.mito.resp, mito != "mito")[,c("protein", "coef", "pval")] %>% distinct()
proteomics.resp.mito.2$mito = "mito"
proteomics.resp.other.2$mito = "other"
proteomics.resp.mito.other.2 = rbind(proteomics.resp.other.2,
                                   proteomics.resp.mito.2)


subset(proteomics.resp.mito.other.2, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito = factor(mito, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito) %>% table() %>%
  fisher.test()
# not sig

# CHECK GO MITO (also sig)
# need to ensure proteins aren't duplicated (because of some proteins have multiple mito annotations)
proteomics.resp.mito2.2 = subset(lsarp.proteomics.combined.red.mito.resp, mito2 == "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.other2.2 = subset(lsarp.proteomics.combined.red.mito.resp, mito2 != "mito")[,c("protein", "coef", "pval")] %>% distinct() %>% arrange(coef)
proteomics.resp.mito2.2$mito2 = "mito"
proteomics.resp.other2.2$mito2 = "other"
proteomics.resp.mito.other2.2 = rbind(proteomics.resp.other2.2,
                                    proteomics.resp.mito2.2)
subset(proteomics.resp.mito.other2.2, pval < 0.05)[,c("protein", "coef", "pval", "mito2")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  mutate(mito2 = factor(mito2, levels=c("other", "mito")))%>%
  dplyr::select(direction, mito2) %>% table() %>%
  fisher.test()
# not sig

# enrich plots

lsarp.proteomics.combined.red.mito.resp.enrich = subset(lsarp.proteomics.combined.red.mito.resp, pval < 0.05)[,c("protein", "coef", "pval", "mito")] %>% distinct() %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(direction) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.combined.red.mito.resp.enrich.plot =lsarp.proteomics.combined.red.mito.resp.enrich%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.combined.red.mito.resp.enrich, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .95*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=4)+
  geom_text(data=subset(lsarp.proteomics.combined.red.mito.resp.enrich, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .95*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=4)+
  
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~"Mitochondrial Protein Enrichment")+
  labs(x="", y="")
lsarp.proteomics.combined.red.mito.resp.enrich.plot


# :: Sens: MPX - Reduced Location--------------------------------------------------


t1 = Sys.time()
lsarp.proteomics.location.red.mito = do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
  print(paste(l))
  
  metadata.subset.2 = subset(proteomics.subset, Location == l)
  
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset.2$code,]
  
  # loop through proteins
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset.2$code,] %>% data.frame()
  metadata.subset.2 = subset(metadata.subset.2, code %in% rownames(data.subset))
  
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset.3 = merge(metadata.subset.2,
                              data.subset[,colnames(data.subset) %in% c("code", protein)],
                              by="code")
    colnames(metadata.subset.3)[colnames(metadata.subset.3) == protein] = "protein"
    
    # 10% of samples
    if(sum(metadata.subset.3$protein != 0) <= nrow(metadata.subset.3)*0.10){
      data.frame(
        location = l,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset.3$protein = scale(log2(metadata.subset.3$protein+proteomics.pseudocount)) 
      metadata.subset.3$Status = factor(metadata.subset.3$Status, levels=c("N", "A"))
      metadata.subset.3$Location = factor(metadata.subset.3$Location, levels=c("TI", "PC", "DC"))
      metadata.subset.3$Time = factor(metadata.subset.3$Time, levels=c("0", "1"))
      # metadata.subset.3$Plate = as.factor(metadata.subset.3$Plate)
      
      lmer.results = lmerTest::lmer(protein ~ Group*Time + Status + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      
      # extract data
      data.frame(
        location = l,
        feature = protein,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))
}))
t2 = Sys.time()
t2 - t1 # 8 min

lsarp.proteomics.location.red.mito = lsarp.proteomics.location.red.mito %>%
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# enrichment
lsarp.proteomics.location.red.mito.enrich = do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
  data.subset = subset(lsarp.proteomics.location.red.mito, location == x)
  fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             location = x)
}))
lsarp.proteomics.location.red.mito.enrich
# TI is significantly depleted
# but not corroborated by GO

# volcano
lsarp.proteomics.location.red.mito.volcano = ggplot(lsarp.proteomics.location.red.mito[,c("coef", "pval", "padj", "protein", "location", "mito")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~location)+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.location.red.mito.volcano

# enrichment plots
lsarp.proteomics.location.red.mito.enrich = subset(lsarp.proteomics.location.red.mito, pval < 0.05)[,c("location", "protein", "coef", "pval", "mito")] %>% distinct() %>%
  group_by(location) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(location, direction, mito) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(location, direction) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.location.red.mito.enrich.plot =lsarp.proteomics.location.red.mito.enrich%>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC"))) %>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.location.red.mito.enrich, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  geom_text(data=subset(lsarp.proteomics.location.red.mito.enrich, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~location)+
  labs(x="", y="")
lsarp.proteomics.location.red.mito.enrich.plot
# no differences


# :: Sens: MPX - Reduced Location Resp --------------------------------------------

t1 = Sys.time()
lsarp.proteomics.location.red.mito.resp = 
  do.call(rbind, lapply(c("RS", "Placebo"), function(y){
    do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
      print(paste(y, " ", l))
      
      metadata.subset = subset(proteomics.subset, Location == l & Group == y)
      
      data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
      
      # loop through proteins
      # keep intersecting samples
      data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
      metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
      
      # apply filtering per location
      proteins.to.keep = data.subset
      proteins.to.keep[proteins.to.keep != 0] <- 1
      proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
      
      data.subset = data.subset[,proteins.to.keep]
      
      # apply variance filter per location
      proteins.to.keep.variance = data.frame(sd = apply((log2(data.subset+proteomics.pseudocount)), 2, sd)) 
      # select top 10%
      proteins.to.keep.variance = proteins.to.keep.variance %>%
        arrange(sd) %>%
        slice_max(order_by=sd, n=round(nrow(proteins.to.keep.variance) *0.10, digits=0))
      
      proteins.to.keep = rownames(proteins.to.keep.variance)
      
      # loop through lmer
      do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
        print(paste0(y, " ", l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
        # protein = proteins.to.keep[1]
        # add protein abundance data by merging
        data.subset$code = rownames(data.subset)
        
        metadata.subset.2 = merge(metadata.subset,
                                  data.subset[,colnames(data.subset) %in% c("code", protein)],
                                  by="code")
        colnames(metadata.subset.2)[colnames(metadata.subset.2) == protein] = "protein"
        
        
        metadata.subset.2$protein = scale(log2(metadata.subset.2$protein+proteomics.pseudocount)) 
        metadata.subset.2$Status = factor(metadata.subset.2$Status, levels=c("N", "A"))
        metadata.subset.2$Location = factor(metadata.subset.2$Location, levels=c("TI", "PC", "DC"))
        metadata.subset.2$Time = factor(metadata.subset.2$Time, levels=c("0", "1"))
        metadata.subset.2$Plate = as.factor(metadata.subset.2$Plate)
        metadata.subset.2$flare.group = factor(metadata.subset.2$flare.group, levels=c("Relapse", "Remit"))
        
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Status + (1|study_id), 
                                      metadata.subset.2) %>% summary() %>% coef() %>% data.frame()
        
        # extract data
        data.frame(
          location = l,
          feature = protein,
          Group = y,
          nfeatures = length(proteins.to.keep),
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
        
      }))
    }))
  }))
t2 = Sys.time()
t2 - t1 # 2 min

lsarp.proteomics.location.red.mito.resp = lsarp.proteomics.location.red.mito.resp %>%
  group_by(location, Group) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# enrichment
lsarp.proteomics.location.red.mito.resp.enrich.fisher = 
  do.call(rbind, lapply(c("Placebo", "RS"), function(y){
    do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
      data.subset = subset(lsarp.proteomics.location.red.mito.resp, location == x & Group == y)
      if(nrow(subset(data.subset, pval < 0.05 & mito == "mito"))==0){
        data.frame(pval = NA,
                   or = NA,
                   Group = y,
                   location = x)
      }else{
        fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
          mutate(direction = sign(coef)+2) %>%
          dplyr::select(direction, mito) %>% table()%>%
          fisher.test()
        data.frame(pval = fisher.result$p.value,
                   or = fisher.result$estimate,
                   Group = y,
                   location = x)
      }
    }))
  }))
# not enough values to assess RS-PC and RS-DC

# volcano
lsarp.proteomics.location.red.mito.resp.volcano = ggplot(lsarp.proteomics.location.red.mito.resp[,c("coef", "pval", "padj", "protein", "location", "mito", "Group", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                     aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=mito, alpha=ifelse(pval < 0.05, 1, 0.5)))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA),
                                     color=mito),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("mito" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("mito" = "red", "other" = "black"), na.value="black")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_grid(Group~location, scales="free")+
  labs(x="Interaction Coefficient",
       y="p value")
lsarp.proteomics.location.red.mito.resp.volcano

# enrich

# enrichment plots
lsarp.proteomics.location.red.mito.resp.enrich = subset(lsarp.proteomics.location.red.mito.resp, pval < 0.05)[,c("location", "protein", "coef", "pval", "mito", "location", "Group")] %>% distinct() %>%
  group_by(location, Group) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(location, direction, mito, Group) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Downregulated", "Upregulated"),
         mito = ifelse(mito == "mito", "Mitochondrial", "Other")) %>%
  group_by(location, direction, Group) %>%
  mutate(perc = Freq / sum(Freq)) 

lsarp.proteomics.location.red.mito.resp.enrich.plot =lsarp.proteomics.location.red.mito.resp.enrich%>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC")))%>%
  ggplot(aes(x=direction, y=mito, fill=direction))+
  # geom_point(shape=21, aes(size=perc), fill="white")+
  # geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=mito),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(lsarp.proteomics.location.red.mito.resp.enrich, mito == "Mitochondrial" & direction == "Downregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  geom_text(data=subset(lsarp.proteomics.location.red.mito.resp.enrich, mito == "Mitochondrial" & direction == "Upregulated"),
            aes(x=direction, y = .97*100,
                label = paste(round(perc, digits=2)*100, "%", sep="")), 
            color="white", size=3)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Mitochondrial" = "red", "Other" = "white"))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_grid(Group~location)+
  labs(x="", y="")
lsarp.proteomics.location.red.mito.resp.enrich.plot
# more depletions, but _due to low numbers_


# :: Sens: MLI - Location -------------------------------------------------


t1 = Sys.time()
asv.lm.group.location = do.call(rbind, lapply(c("DC", "PC", "TI"), function(location){
  # subset to location
  metadata.subset = subset(lsarp.mli.ps.meta, Time %in% c(0,1) & Location == location) %>% distinct()
  data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
  # keep intersecting samples
  data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
  metadata.subset = subset(metadata.subset, (Sample) %in% rownames(data.subset))
  # apply filtering per location
  asvs.to.keep = data.subset
  asvs.to.keep[asvs.to.keep != 0] <- 1
  asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
  # loop through lmer
  do.call(rbind, lapply(asvs.to.keep, function(asv){
    print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
    # asv = asvs.to.keep[1]
    # add asv abundance data by merging
    data.subset$Sample = rownames(data.subset)
    
    metadata.subset = merge(metadata.subset %>% mutate(Sample = (Sample)),
                            data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
    colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
    
    if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
      data.frame(
        location = location,
        feature = asv,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
      metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
      metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
      metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
      lmer.results = lmerTest::lmer(asv ~ Group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
      # extract data
      data.frame(
        location = location,
        feature = asv,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))}))
t2 = Sys.time()
t2 - t1
# calculate padj
asv.lm.group.location = asv.lm.group.location %>% 
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>% 
  mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
  arrange(pval)
# Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  


# append butyrogen annotation and LCA
asv.lm.group.location = merge(asv.lm.group.location %>% mutate(OTU = feature),
                                         lsarp.mli.butyrogens[,c("OTU", "lca.gg2", "butyrogen.ii")],
                                         by="OTU") %>% distinct()
# append butyrogen annotation and LCA
asv.lm.group.location = merge(asv.lm.group.location %>% mutate(OTU = feature),
                                         lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")],
                                         by="OTU") %>% distinct()
# volcano plots
asv.lm.group.location.volcano = ggplot(asv.lm.group.location[,c("location", "coef","butyrogen.i","lca.gg2", "pval", "padj", "feature")] %>% distinct() %>% 
                                         mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                  aes(x=as.numeric(coef), y=padj))+
  geom_point(shape = 21, aes(fill=butyrogen.i, alpha=ifelse(padj < 0.20, 1, 0.5)))+
  geom_hline(yintercept=0.20, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", lca.gg2), NA),
                                     color=butyrogen.i),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("butyrogen" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("butyrogen" = "red", "other" = "black"), na.value="black")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~location)+
  labs(x="Interaction Coefficient",
       y="FDR")
asv.lm.group.location.volcano

# need to prevent double-counting proteins
butyrogens.i.but.2 = subset(asv.lm.group.location, butyrogen.i == "butyrogen")[,c("location", "feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.other.2 = subset(asv.lm.group.location, butyrogen.i != "butyrogen")[,c("location", "feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.but.2$butyrogen.i = "butyrogen"
butyrogens.i.other.2$butyrogen.i = "other"
butyrogens.i.but.other.2 = rbind(butyrogens.i.other.2,
                               butyrogens.i.but.2)
# any edge cases:
table(butyrogens.i.but.2$feature) %>% range()
table(butyrogens.i.other.2$feature) %>% range()
# no; perfect, continue


# enrichment
asv.lm.group.location.enrich = do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
  data.subset = subset(asv.lm.group.location, location == x)
  fisher.result = subset(data.subset, pval < 0.05)[,c("OTU", "butyrogen.i", "coef", "pval")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, butyrogen.i) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             location = x)
}))
# not sig

asv.lm.group.location.data = subset(asv.lm.group.location, pval < 0.05)[,c("location", "feature", "coef", "pval", "butyrogen.i")] %>% distinct() %>%
  group_by(location) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.i, location) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.i = ifelse(butyrogen.i == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction, location) %>%
  mutate(perc = Freq / sum(Freq))

asv.lm.group.location.plot = asv.lm.group.location.data %>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC"))) %>%
  ggplot()+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.i),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(asv.lm.group.location.data, butyrogen.i == "Butyrogen"),
            aes(
              x=direction, y = .90*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="white", size=5)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_wrap(~location)+
  labs(x="", y="")
asv.lm.group.location.plot

# none passed FDR < 0.20
# no enrichments were sig

# :: Sens: MLI - Location Resp -------------------------------------------------

# only 1 Placebo Remit, so cannot run stats
# only use RS group

  t1 = Sys.time()
asv.lm.resp.location = 
  do.call(rbind, lapply(c("RS"), function(gr){
    
  do.call(rbind, lapply(c("DC", "PC", "TI"), function(location){
    # subset to location
    metadata.subset = subset(lsarp.mli.ps.meta, Group == gr & Time %in% c(0, 1) & Location == location) %>% distinct()
    # add Flare.group
    metadata.subset = merge(metadata.subset,
                            lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(), by="HM")
    
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, (Sample) %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(gr, " ", location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[192]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = (Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      # 20% prevalence filter
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          Group = gr,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Location = factor(metadata.subset$Location, levels=c("TI", "PC", "DC"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        
        lmer.results = lmerTest::lmer(asv ~ flare.group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          Group = gr,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))
    }))
  }))
  t2 = Sys.time()
  t2 - t1 # 1 min
  # calculate padj
  asv.lm.resp.location = asv.lm.resp.location %>% 
    group_by(Group, location) %>%
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
# append butyrogen annotation and LCA
  asv.lm.resp.location = merge(asv.lm.resp.location %>% mutate(OTU = feature),
                                              lsarp.mli.gg138.tax.df[,c("OTU", "lca.gg138", "butyrogen.i")],
                                              by="OTU")
  asv.lm.resp.location = merge(asv.lm.resp.location %>% mutate(OTU = feature),
                                              lsarp.mli.butyrogens[,c("OTU", "lca.gg2", "butyrogen.ii")],
                                              by="OTU")
# volcano plots
  asv.lm.resp.location.volcano = ggplot(asv.lm.resp.location[,c("Group", "coef","lca.gg2", "butyrogen.i", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                       aes(x=as.numeric(coef), y=padj))+
  geom_point(shape = 21, aes(fill=butyrogen.i, alpha=ifelse(padj < 0.20, 1, 0.5)))+
  geom_hline(yintercept=0.20, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", lca.gg2), NA),
                                     color=butyrogen.i),
                                 size=2.5)+
  # scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=c("butyrogen" = "red", "other" = "white"), na.value="white")+
  scale_color_manual(values=c("butyrogen" = "red", "other" = "black"), na.value="black")+
  
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_grid(Group ~ location)+
  labs(x="Interaction Coefficient",
       y="FDR")
  asv.lm.resp.location.volcano

# need to prevent double-counting proteins
butyrogens.i.resp.but.2 = subset(asv.lm.resp.location, butyrogen.i == "butyrogen")[,c("Group", "location", "feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.other.2 = subset(asv.lm.resp.location, butyrogen.i != "butyrogen")[,c("Group", "location", "feature", "coef", "pval","padj")] %>% distinct()
butyrogens.i.resp.but.2$butyrogen.i = "butyrogen"
butyrogens.i.resp.other.2$butyrogen.i = "other"
butyrogens.i.resp.but.other.2 = rbind(butyrogens.i.resp.other.2,
                                    butyrogens.i.resp.but.2)
# any edge cases:
table(butyrogens.i.resp.but.2$feature) %>% range()
table(butyrogens.i.resp.other.2$feature) %>% range()
# no; perfect, continue


# enrichment
butyrogens.i.resp.but.other.2.fisher = 
  do.call(rbind, lapply(c("Placebo", "RS"), function(y){
    do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
      data.subset = subset(butyrogens.i.resp.but.other.2, location == x & Group == y)
      if(nrow(subset(data.subset, pval < 0.05 & butyrogen.i == "butyrogen"))==0){
        data.frame(pval = NA,
                   or = NA,
                   Group = y,
                   location = x)
      }else{
        fisher.result = subset(data.subset, pval < 0.05)[,c("butyrogen.i", "coef", "pval", "feature")] %>% distinct() %>%
          mutate(direction = sign(coef)+2) %>%
          dplyr::select(direction, butyrogen.i) %>% table()%>%
          fisher.test()
        data.frame(pval = fisher.result$p.value,
                   or = fisher.result$estimate,
                   Group = y,
                   location = x)
      }
    }))
  }))
# not enough values to assess RS-PC and RS-DC


asv.lm.resp.location.enrich = subset(asv.lm.resp.location, pval < 0.05)[,c("Group", "location", "feature", "coef", "pval", "butyrogen.i")] %>% distinct() %>%
  group_by(Group, location) %>%
  mutate(direction = sign(coef)+2) %>%
  dplyr::select(direction, butyrogen.i, Group, location) %>% table() %>% as.data.frame() %>%
  mutate(direction = ifelse(direction == 1, "Depleted", "Enriched"),
         butyrogen.i = ifelse(butyrogen.i == "butyrogen", "Butyrogen", "Other")) %>%
  group_by(direction, Group, location) %>%
  mutate(perc = Freq / sum(Freq))

asv.lm.resp.location.enrich.plot = asv.lm.resp.location.enrich %>%
  mutate(location = factor(location, levels=c("TI", "PC", "DC"))) %>%
  ggplot(aes(x=direction, y=butyrogen.i, fill=direction))+
  geom_bar(position="stack", stat="identity",
           aes(x=direction, y=perc*100, fill=butyrogen.i),
           color="black")+
  #geom_point(shape=21, aes(size=perc), fill="white")+
  #geom_point(shape=21, aes(size=perc), alpha=0.6)+
  geom_text(data=subset(asv.lm.resp.location.enrich, butyrogen.i == "Butyrogen"),
            aes(
              x=direction, y = .95*100,
              label = paste(round(perc, digits=2)*100, "%", sep="")), color="white", size=4)+
  #scale_fill_manual(values=c("Depleted" = "blue", "Enriched" = "red"))+
  scale_fill_manual(values=c("Butyrogen" = "red", "Other" = "white"))+
  # scale_size_continuous(range=c(20,60))+
  theme_minimal()+
  theme(legend.position="none",
        strip.background = element_rect(color="black"),
        strip.text=element_text(size=10))+
  facet_grid(Group ~ location) +
  labs(x="", y="")
asv.lm.resp.location.enrich.plot

# very noisy


# :: Sens PLOTS -----------------------------------------------------------

# MPX, placebo - RESPONDERS (note; n = 1 'placebo responder')
(proteomics.data.placebo.interaction.resp.lm.mito.volcano+
   proteomics.data.resp.interaction.placebo.lm.mito.plot)|
# MLI, placebo - RESPONDERS (note; n = 1 'placebo responder')
(asv.data.placebo.rs.interaction.lm.resp.volcano+
   mli.resp.butyrogen.i.enrichment.plac.plot )
# Summary: looks similar to RS, but likely underpowered

# MPX, default combined - GROUP - GO
(proteomics.data.group.interaction.lm.mito2.volcano+
proteomics.data.group.interaction.lm.mito2.plot)|
# MPX, default combined - RESPONDERS - GO
(proteomics.data.group.interaction.resp.lm.mito2.volcano+
proteomics.mito2.gene.enrichment.plot)
# Summary: GO terms recapitulate MitoCarta3.0


# MPX, default location - GROUP
(lsarp.proteomics.location.mito.volcano+
lsarp.proteomics.location.mito.enrich.plot)/
# MPX, default location - RESPONDERS
(lsarp.proteomics.location.mito.resp.volcano+
lsarp.proteomics.location.mito.resp.enrich.plot)+
  patchwork::plot_layout(heights=c(1,3))
# Summary: A depletion of Mito in TI (Group level)
# and a depletion in RS responders, but far fewer p < 0.05


# MPX, reduced combined - GROUP
(lsarp.proteomics.combined.red.mito.volcano+
lsarp.proteomics.combined.red.mito.enrich.plot)|
# MPX, reduced combined - RESPONDERS
(lsarp.proteomics.combined.red.mito.resp.volcano+
lsarp.proteomics.combined.red.mito.resp.enrich.plot)
# Summary: Failed to recapitulate original analysis
# no differences between groups; group or responders


# MPX, reduced location - GROUP
(lsarp.proteomics.location.red.mito.volcano+
lsarp.proteomics.location.red.mito.enrich.plot)/
# MPX, reduced location - RESPONDERS
(lsarp.proteomics.location.red.mito.resp.volcano+
lsarp.proteomics.location.red.mito.resp.enrich.plot)+
  patchwork::plot_layout(heights=c(1,3))
# Summary: noisy results; TI is depleted, DC is enriched
# far fewer p < 0.05 in responders/non-responders


# MLI, location - GROUP
(asv.lm.group.location.volcano+
asv.lm.group.location.plot)/
# MLI, location - RESPONDERS
(asv.lm.resp.location.volcano+
asv.lm.resp.location.enrich.plot)
# Summary: 1 taxa achieves FDR < 0.20
# enriched butyrogens in DC and PC responders
# depleted in TI
# underpowered to compare to Placebo

# :: ----------------------------------------------------------------------


# >>> TEMPORARY -----------------------------------------------------------

metadata.subset = subset(proteomics.metadata, Time %in% c(0,1))
data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]

lsarp.proteomics.analysis.2.0 = do.call(rbind, lapply(c("DC", "PC", "TI", "combined"), function(l){
  print(paste(l))
   if(l == "combined"){
    metadata.subset.2 = metadata.subset
    # reduce to MOST INFLAMED out of pairs of samples
    # i.e. if A and N were taken, use only A
    metadata.subset.2 = metadata.subset.2 %>%
      group_by(HM, Time, Location) %>%
      mutate(status.factor = factor(Status, levels=c("A", "N"))) %>%
      subset(Status == min(Status))
  }else{
  metadata.subset.2 = subset(metadata.subset, Location == l)
  # reduce to MOST INFLAMED out of pairs of samples
  # i.e. if A and N were taken, use only A
  metadata.subset.2 = metadata.subset.2 %>%
    group_by(HM, Time, Location) %>%
    mutate(status.factor = factor(Status, levels=c("A", "N"))) %>%
    subset(Status == min(Status))
  }
  
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset.2$code,]
  
  # loop through proteins
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset.2$code,] %>% data.frame()
  metadata.subset.2 = subset(metadata.subset.2, code %in% rownames(data.subset))
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset.3 = merge(metadata.subset.2,
                            data.subset[,colnames(data.subset) %in% c("code", protein)],
                            by="code")
    colnames(metadata.subset.3)[colnames(metadata.subset.3) == protein] = "protein"
    
    if(sum(metadata.subset.3$protein != 0) <= 3){
      data.frame(
        location = l,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset.3$protein = scale(log2(metadata.subset.3$protein+proteomics.pseudocount)) 
      metadata.subset.3$Status = factor(metadata.subset.3$Status, levels=c("N", "A"))
      metadata.subset.3$Location = factor(metadata.subset.3$Location, levels=c("TI", "PC", "DC"))
      metadata.subset.3$Time = factor(metadata.subset.3$Time, levels=c("0", "1"))
      if(l == "combined"){
      lmer.results = lmerTest::lmer(protein ~ Group*Time + Location + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      }else{
      lmer.results = lmerTest::lmer(protein ~ Group*Time + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      }
      # extract data
      data.frame(
        location = l,
        feature = protein,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))
  }))
  
  
  
lsarp.proteomics.analysis.2.0 %>%
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.data.group.interaction.lm.mito[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# AT GROUP LEVEL
# only PC has sig proteins (FDR < 0.20), and only 2 are mito
  

lsarp.proteomics.analysis.resp.2.0 = do.call(rbind, lapply(c("DC", "PC", "TI", "combined"), function(l){
  print(paste(l))
  if(l == "combined"){
    metadata.subset.2 = metadata.subset %>% subset(Group == "RS")
    # reduce to MOST INFLAMED out of pairs of samples
    # i.e. if A and N were taken, use only A
    metadata.subset.2 = metadata.subset.2 %>%
      group_by(HM, Time, Location) %>%
      mutate(status.factor = factor(Status, levels=c("A", "N"))) %>%
      subset(Status == min(Status))
  }else{
    metadata.subset.2 = subset(metadata.subset, Location == l) %>% subset(Group == "RS")
    # reduce to MOST INFLAMED out of pairs of samples
    # i.e. if A and N were taken, use only A
    metadata.subset.2 = metadata.subset.2 %>%
      group_by(HM, Time, Location) %>%
      mutate(status.factor = factor(Status, levels=c("A", "N"))) %>%
      subset(Status == min(Status))
  }
  
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset.2$code,]
  
  # loop through proteins
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset.2$code,] %>% data.frame()
  metadata.subset.2 = subset(metadata.subset.2, code %in% rownames(data.subset))
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset.3 = merge(metadata.subset.2,
                              data.subset[,colnames(data.subset) %in% c("code", protein)],
                              by="code")
    colnames(metadata.subset.3)[colnames(metadata.subset.3) == protein] = "protein"
    
    if(sum(metadata.subset.3$protein != 0) <= 3){
      data.frame(
        location = l,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset.3$protein = scale(log2(metadata.subset.3$protein+proteomics.pseudocount)) 
      metadata.subset.3$Status = factor(metadata.subset.3$Status, levels=c("N", "A"))
      metadata.subset.3$Location = factor(metadata.subset.3$Location, levels=c("TI", "PC", "DC"))
      metadata.subset.3$flare.group = factor(metadata.subset.3$flare.group, levels=c("Relapse", "Remit"))
       metadata.subset.3$Time = factor(metadata.subset.3$Time, levels=c("0", "1"))
      if(l == "combined"){
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Location + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      }else{
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      }
      # extract data
      data.frame(
        location = l,
        feature = protein,
        coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
    }
  }))
}))



lsarp.proteomics.analysis.resp.2.0 %>%
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.data.group.interaction.lm.mito[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# AT GROUP LEVEL
# only PC has sig proteins (FDR < 0.20), and only 2 are mito


# :: MPX LM ---------------------------------------------------------------

# create function to perform interaction analysis per location


lsarp.proteomics.analysis.3.0 = do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
  print(paste(l))
  
  metadata.subset.2 = subset(proteomics.metadata, Location == l)
  
  data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset.2$code,]
  
  # loop through proteins
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset.2$code,] %>% data.frame()
  metadata.subset.2 = subset(metadata.subset.2, code %in% rownames(data.subset))
  
  # apply filtering per location
  proteins.to.keep = data.subset
  proteins.to.keep[proteins.to.keep != 0] <- 1
  proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  
  # loop through lmer
  do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
    print(paste0(l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset.3 = merge(metadata.subset.2,
                              data.subset[,colnames(data.subset) %in% c("code", protein)],
                              by="code")
    colnames(metadata.subset.3)[colnames(metadata.subset.3) == protein] = "protein"
    
    # 10% of samples
    if(sum(metadata.subset.3$protein != 0) <= nrow(metadata.subset.3)*0.10){
      data.frame(
        location = l,
        feature = protein,
        coef = NA,
        pval = NA)
    }else{
      metadata.subset.3$protein = scale(log2(metadata.subset.3$protein+proteomics.pseudocount)) 
      metadata.subset.3$Status = factor(metadata.subset.3$Status, levels=c("N", "A"))
      metadata.subset.3$Location = factor(metadata.subset.3$Location, levels=c("TI", "PC", "DC"))
      metadata.subset.3$Time = factor(metadata.subset.3$Time, levels=c("0", "1"))
      metadata.subset.3$Plate = as.factor(metadata.subset.3$Plate)
      
      lmer.results = lmerTest::lmer(protein ~ Group*Time + Status + Plate + (1|study_id), metadata.subset.3) %>% summary() %>% coef() %>% data.frame()
      
      # extract data
      data.frame(
        location = l,
        feature = protein,
        coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
        pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
    }
  }))
}))

lsarp.proteomics.analysis.3.0 = lsarp.proteomics.analysis.3.0 %>%
  group_by(location) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

# enrichment
mito.enrich.location.group = do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
  data.subset = subset(lsarp.proteomics.analysis.3.0, location == x)
  fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
    mutate(direction = sign(coef)+2) %>%
    dplyr::select(direction, mito2) %>% table()%>%
    fisher.test()
  data.frame(pval = fisher.result$p.value,
             or = fisher.result$estimate,
             location = x)
}))
# none are significant


# :: :: Plots ---------------------------------------------------------

# volcano plot


# heatmap


# enrichment bars



# :: MPX LM Responders ----------------------------------------------------


lsarp.proteomics.analysis.responders.3.0 = 
  do.call(rbind, lapply(c("RS", "Placebo"), function(y){
    do.call(rbind, lapply(c("DC", "PC", "TI"), function(l){
      print(paste(y, " ", l))
      
      metadata.subset = subset(proteomics.metadata, Location == l & Group == y)
      
      data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
      
      # loop through proteins
      # keep intersecting samples
      data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
      metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
      
      # apply filtering per location
      proteins.to.keep = data.subset
      proteins.to.keep[proteins.to.keep != 0] <- 1
      proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
      
      data.subset = data.subset[,proteins.to.keep]
      
      # apply variance filter per location
      proteins.to.keep.variance = data.frame(sd = apply((log2(data.subset+proteomics.pseudocount)), 2, sd)) 
      # select top 10%
      proteins.to.keep.variance = proteins.to.keep.variance %>%
        arrange(sd) %>%
        slice_max(order_by=sd, n=round(nrow(proteins.to.keep.variance) *0.10, digits=0))
      
      proteins.to.keep = rownames(proteins.to.keep.variance)
      
      # loop through lmer
      do.call(rbind, lapply(proteins.to.keep, function(protein){
        print(paste0(y, " ", l, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
        # protein = proteins.to.keep[1]
        # add protein abundance data by merging
        data.subset$code = rownames(data.subset)
        
        metadata.subset.2 = merge(metadata.subset,
                                  data.subset[,colnames(data.subset) %in% c("code", protein)],
                                  by="code")
        colnames(metadata.subset.2)[colnames(metadata.subset.2) == protein] = "protein"
        
        
        metadata.subset.2$protein = scale(log2(metadata.subset.2$protein+proteomics.pseudocount)) 
        metadata.subset.2$Status = factor(metadata.subset.2$Status, levels=c("N", "A"))
        metadata.subset.2$Location = factor(metadata.subset.2$Location, levels=c("TI", "PC", "DC"))
        metadata.subset.2$Time = factor(metadata.subset.2$Time, levels=c("0", "1"))
        metadata.subset.2$Plate = as.factor(metadata.subset.2$Plate)
        metadata.subset.2$flare.group = factor(metadata.subset.2$flare.group, levels=c("Relapse", "Remit"))
        
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Status + Plate + (1|study_id), 
                                      metadata.subset.2) %>% summary() %>% coef() %>% data.frame()
        
        # extract data
        data.frame(
          location = l,
          feature = protein,
          Group = y,
          nfeatures = length(proteins.to.keep),
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
        
      }))
    }))
  }))

lsarp.proteomics.analysis.responders.3.0 = lsarp.proteomics.analysis.responders.3.0 %>%
  group_by(location, Group) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  #subset(padj < 0.20) %>%
  mutate(protein = feature) %>%
  merge(distinct(proteomics.mito.map[,c("protein", "mito", "mito2")]),
        by="protein") %>%
  arrange(pval) %>% as.data.frame()

lsarp.proteomics.analysis.responders.3.0 %>%
  subset(Group == "Placebo") %>%
  arrange(padj)

# enrichment
mito.enrich.location.resp = 
  do.call(rbind, lapply(c("Placebo", "RS"), function(y){
    do.call(rbind, lapply(c("TI", "PC","DC"), function(x){
      data.subset = subset(lsarp.proteomics.analysis.responders.3.0, location == x & Group == y)
      fisher.result = subset(data.subset, pval < 0.05)[,c("protein", "coef", "pval", "mito", "mito2")] %>% distinct() %>%
        mutate(direction = sign(coef)+2) %>%
        dplyr::select(direction, mito) %>% table()%>%
        fisher.test()
      data.frame(pval = fisher.result$p.value,
                 or = fisher.result$estimate,
                 Group = y,
                 location = x)
    }))
  }))
# none are significant


# :: MPX limma -----------------------------------------------------------

library(limma)
library(edgeR)

proteomics.data.clean
proteomics.metadata

limma.metadata = proteomics.metadata %>%
  subset(Time %in% c(0,1))
limma.metadata = limma.metadata[limma.metadata$code %in% rownames(proteomics.data.clean),]

# Create DESeqDataSet (counts is your count matrix, genes x samples)
dge <- DGEList(counts = t(proteomics.data.clean))   # yes, even if already normalized

# 2. Set up the design matrix with participant pairing and interaction
# Relevel factors
limma.metadata$Group <- relevel(factor(limma.metadata$Group), ref = "Placebo")
limma.metadata$Time  <- relevel(factor(limma.metadata$Time),  ref = "0")
limma.metadata$Location  <- relevel(factor(limma.metadata$Location),  ref = "TI")
limma.metadata$Status  <- relevel(factor(limma.metadata$Status),  ref = "N")

limma.metadata$HM <- factor(limma.metadata$HM)

design <- model.matrix(
  ~ Group * Time + Location + Status,
  data = limma.metadata
)# This creates: groupControl, groupTreatment, timePre, timePost, groupTreatment:timePost, and participant terms
colnames(design) <- make.names(colnames(design))

# 3. voom: normalize (if needed) and transform to logCPM with precision weights
v <- voom(dge, design = design, plot = TRUE)  

# Better duplicateCorrelation for true repeated measures (recommended!)
# This explicitly models the within-participant correlation

cont.matrix <- makeContrasts(
  Interaction        = GroupRS.Time1,
  GroupRS_main       = GroupRS,
  TimeEffect_Placebo = Time1,
  TimeEffect_RS      = Time1 + GroupRS.Time1,
  levels = design
)

corfit <- duplicateCorrelation(v, design, block = limma.metadata$HM)
v2 <- voom(dge, design = design, block = limma.metadata$HM, correlation = corfit$consensus)

fit_dup <- lmFit(v2, design, block = limma.metadata$HM, correlation = corfit$consensus)
fit_dup <- contrasts.fit(fit_dup, cont.matrix)
fit_dup <- eBayes(fit_dup, robust = TRUE)

# Now use topTable(fit_dup, ...) — this is the gold-standard for paired designs
res_interaction_dup <- topTable(fit_dup, coef = "Interaction", number = Inf, sort.by = "P")
res_interaction_dup %>% head()

with(res_interaction_dup,
     plot(logFC, -log10(adj.P.Val), pch = 16))

# check
df <- data.frame(
  expr = v2$E["DOB", ],
  Group = limma.metadata$Group,
  Time  = limma.metadata$Time,
  HM    = limma.metadata$HM
)

ggplot(df, aes(Time, expr, group = HM, color = Group)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group),
               fun = mean, geom = "line", linewidth = 1.2)
# looks good


# :: MPX limma resp -------------------------------------------------------

limma.metadata.resp = subset(limma.metadata, Group == "RS")

# Create DESeqDataSet (counts is your count matrix, genes x samples)
dge.resp <- DGEList(counts = t(proteomics.data.clean[rownames(proteomics.data.clean) %in% limma.metadata.resp$code,]))   # yes, even if already normalized

# 2. Set up the design matrix with participant pairing and interaction
# Relevel factors
limma.metadata.resp$flare.group <- relevel(factor(limma.metadata.resp$flare.group), ref = "Relapse")
limma.metadata.resp$Time  <- relevel(factor(limma.metadata.resp$Time),  ref = "0")
limma.metadata.resp$Location  <- relevel(factor(limma.metadata.resp$Location),  ref = "TI")
limma.metadata.resp$Status  <- relevel(factor(limma.metadata.resp$Status),  ref = "N")

limma.metadata.resp$HM <- factor(limma.metadata.resp$HM)

design.resp <- model.matrix(
  ~ flare.group * Time + Location + Status,
  data = limma.metadata.resp
)# This creates: groupControl, groupTreatment, timePre, timePost, groupTreatment:timePost, and participant terms
colnames(design.resp) <- make.names(colnames(design.resp))

# 3. voom: normalize (if needed) and transform to logCPM with precision weights
v.resp <- voom(dge.resp, design = design.resp, plot = TRUE)  

# Better duplicateCorrelation for true repeated measures (recommended!)
# This explicitly models the within-participant correlation

cont.matrix.resp <- makeContrasts(
  Interaction        = flare.groupRemit.Time1,
  FlareGroup_main       = flare.groupRemit,
  TimeEffect_Relapse = Time1,
  TimeEffect_Remit      = Time1 + flare.groupRemit.Time1,
  levels = design.resp
)

corfit.resp <- duplicateCorrelation(v.resp, design.resp, block = limma.metadata.resp$HM)
v2.resp <- voom(dge.resp, design = design.resp, block = limma.metadata.resp$HM, correlation = corfit.resp$consensus)

fit_dup.resp <- lmFit(v2.resp, design.resp, block = limma.metadata.resp$HM, correlation = corfit.resp$consensus)
fit_dup.resp <- contrasts.fit(fit_dup.resp, cont.matrix.resp)
fit_dup.resp <- eBayes(fit_dup.resp, robust = TRUE)

# Now use topTable(fit_dup, ...) — this is the gold-standard for paired designs
res_interaction_dup.resp <- topTable(fit_dup.resp, coef = "Interaction", number = Inf, sort.by = "P")
res_interaction_dup.resp %>% head()
# even less significant

with(res_interaction_dup.resp,
     plot(logFC, -log10(adj.P.Val), pch = 16))

# check
df <- data.frame(
  expr = v2$E["SMDC1", ],
  Group = limma.metadata$Group,
  Time  = limma.metadata$Time,
  HM    = limma.metadata$HM
)

ggplot(df, aes(Time, expr, group = HM, color = Group)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = Group),
               fun = mean, geom = "line", linewidth = 1.2)
# looks bad

# >>> Old -----------------------------------------------------------------

# :: MPX: GO Enrichment //defunct ---------------------------------------------------

proteomics.upregulated = subset(proteomics.data.group.interaction.lm.mito, pval < 0.05 & coef > 0)$gene %>% unique()
proteomics.downregulated = subset(proteomics.data.group.interaction.lm.mito, pval < 0.05 & coef < 0)$gene %>% unique()

# note: clusterProfiler does not work
# do it manually
proteomics.background = table(go_annotations$name_1006) %>% data.frame() #%>% arrange(-Freq) %>% head(n=20)
colnames(proteomics.background) = c("name_1006", "bg")

proteomics.upregulated.go = subset(go_annotations, hgnc_symbol %in% proteomics.upregulated) %>% dplyr::select(name_1006) %>% table() %>% data.frame()
proteomics.downregulated.go = subset(go_annotations, hgnc_symbol %in% proteomics.downregulated)  %>% dplyr::select(name_1006) %>% table() %>% data.frame()
proteomics.total.go = subset(go_annotations, hgnc_symbol %in% c(proteomics.upregulated,
                                                                proteomics.downregulated)) %>% dplyr::select(name_1006) %>% table() %>% data.frame()

# loop through processes
proteomics.upregulated.enrich = do.call(rbind, lapply(unique(proteomics.total.go$name_1006), function(go){
  observed.go = subset(proteomics.total.go, name_1006 == go)
  observed.total = length(proteomics.upregulated) # total genes "sig" affected
  bg.go = subset(proteomics.background, as.character(name_1006) == go)
  bg.total = length(unique(go_annotations_bio$hgnc_symbol)) # total genes observed (annotated and eligible for GO enrichment)
  prob = phyper(q = observed.go$Freq-1, # sig genes involved in GO X
                m = observed.total, # total sig genes
                n = bg.total, # total genes
                k = bg.go$bg, # total genes involved in GO X
                lower.tail = FALSE) # testing for enrichment, not depletion, too
  data.frame(pval = 1 - prob,
             GO = go)
}))
proteomics.upregulated.enrich$padj = p.adjust(proteomics.upregulated.enrich$pval, method="BH")
proteomics.upregulated.enrich$sig = ifelse(proteomics.upregulated.enrich$padj< 0.05, "*", "")
proteomics.upregulated.enrich %>% arrange(pval) %>% head(n=20)

# nothing interpretable


# :: MPX ------------------------------------------------------------------


# :: Question 1: Baseline -------------------------------------------------

# Question: are there significant differences between N and A in groups

asv.metadata.baseline = subset(lsarp.mli.ps.meta, Time == "0") %>% distinct()


# loop through locations and perform LMER
asv.data.baseline = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat) %in% asv.metadata.baseline$Sample,]

rerun=T
if(rerun==T){
  
  t1 = Sys.time()
  asv.data.baseline.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(asv.metadata.baseline, Location == location) %>% data.frame()
    data.subset = asv.data.baseline[rownames(asv.data.baseline) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,] %>% as.data.frame()
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add protein abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)], by="Sample")
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        lmer.results = lm(asv ~ Status, metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "StatusA",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "StatusA",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 20 sec
  # calculate padj
  asv.data.baseline.lm = asv.data.baseline.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(padj)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.baseline.lm, "./host_proteomics/asv.data.baseline.lm.Rds")
}
asv.data.baseline.lm = readRDS("./host_proteomics/asv.data.baseline.lm.Rds")

# arrange by significance
asv.data.baseline.lm %>% arrange(padj)

# volcano plots
ggplot(asv.data.baseline.lm[,c("coef", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
       aes(x=as.numeric(coef), y=-log10(pval)))+
  geom_point(shape = 21, aes(fill=coef))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient")+
  facet_wrap(~location)


# :: Question 2: RS -------------------------------------------------------

# Question: are there significant differences before and after RS

asv.metadata.rs = subset(lsarp.mli.ps.meta, Group == "RS") %>% distinct()

# loop through locations and perform LMER
asv.data.rs = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat) %in% asv.metadata.rs$Sample,]

if(rerun==T){
  t1 = Sys.time()
  asv.data.rs.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(asv.metadata.rs, Location == location)
    data.subset = asv.data.rs[rownames(asv.data.rs) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,] %>% data.frame()
    metadata.subset = subset(metadata.subset, Sample %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add protein abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)])
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(asv ~ Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time() # ~20 min
  t2 - t1
  # calculate padj
  asv.data.rs.lm = asv.data.rs.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.rs.lm, "./host_proteomics/asv.data.rs.lm.Rds")
}
asv.data.rs.lm = readRDS("./host_proteomics/asv.data.rs.lm.Rds")

# volcano plots
ggplot(asv.data.rs.lm[,c("coef", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
       aes(x=as.numeric(coef), y=-log10(pval)))+
  geom_point(shape = 21, aes(fill=coef))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient")+
  facet_wrap(~location)

# :: Question 3: Placebo -------------------------------------------------------

# Question: are there significant differences before and after Placebo

asv.metadata.placebo = subset(lsarp.mli.ps.meta, Group == "Placebo") %>% distinct()

# loop through locations and perform LMER
asv.data.placebo = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat) %in% asv.metadata.placebo$Sample,]

if(rerun==T){
  t1 = Sys.time()
  asv.data.placebo.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(asv.metadata.placebo, Location == location)
    data.subset = asv.data.placebo[rownames(asv.data.placebo) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, Sample %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)])
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(asv ~ Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1
  # calculate padj
  asv.data.placebo.lm = asv.data.placebo.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.placebo.lm, "./host_proteomics/asv.data.placebo.lm.Rds")
}
asv.data.placebo.lm = readRDS("./host_proteomics/asv.data.placebo.lm.Rds")


# volcano plots
ggplot(asv.data.placebo.lm[,c("coef", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
       aes(x=as.numeric(coef), y=-log10(pval)))+
  geom_point(shape = 21, aes(fill=coef))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient")+
  facet_wrap(~location)


# :: MPX Interactions (per location) -------------------------------------------------------

# Question: are there significant differences before and after RS relative to Placebo

# loop through locations and perform LMER

if(rerun==T){
  t1 = Sys.time()
  proteomics.data.interaction.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(proteomics.metadata, Location == location)
    data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(protein ~ Group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1
  # calculate padj
  proteomics.data.interaction.lm = proteomics.data.interaction.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.interaction.lm, "./host_proteomics/proteomics.data.interaction.lm.Rds")
}
proteomics.data.interaction.lm = readRDS("./host_proteomics/proteomics.data.interaction.lm.Rds")

# add mito linker
colnames(proteomics.data.interaction.lm)[colnames(proteomics.data.interaction.lm) == "protein"] = "protein.x"
proteomics.data.interaction.lm = merge(proteomics.data.interaction.lm,
                                       proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.interaction.lm.mito = merge(proteomics.data.interaction.lm,
                                            mitocarta.pathways.expanded,
                                            by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""))
# arrange by significance
proteomics.data.interaction.lm.mito %>% arrange(pval)

# :: MPX Volcano (per location) ----------------------------------------------------------

# volcano plots
proteomics.data.interaction.lm.mito.volcano = ggplot(proteomics.data.interaction.lm.mito[,c("coef", "pval", "padj", "protein", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                     aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=coef))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Interaction Coefficient",
       y="Unadjusted p value")+
  facet_wrap(~location)
proteomics.data.interaction.lm.mito.volcano

# no proteins with FDR-adjusted p value < 0.20 (or even 0.80)


# odd shape; I suspect on account of the normalization strategy; can this be accounted for in the lm?
proteomics.data.interaction.lm.mito$MitoPathway %>% unique()

subset(proteomics.data.interaction.lm.mito, mito=="mito" & padj < 0.20) %>% arrange(padj)

subset(proteomics.data.interaction.lm.mito, padj < 0.20) %>% arrange(padj)

proteomics.data.interaction.lm.mito %>% arrange(pval)


# :: MPX Interactions: Responders (per location) ------------------------------------------

if(rerun==T){
  t1 = Sys.time()
  proteomics.data.interaction.resp.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location & RS group
    metadata.subset = subset(proteomics.metadata, Location == location & Group == "RS")
    # add Flare.group
    metadata.subset = merge(metadata.subset,
                            lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(), by="HM")
    data.subset = proteomics.data.clean[rownames(proteomics.data.clean) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(protein ~ flare.group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 7 min
  # calculate padj
  proteomics.data.interaction.resp.lm = proteomics.data.interaction.resp.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.interaction.resp.lm, "./host_proteomics/proteomics.data.interaction.resp.lm.Rds")
}
proteomics.data.interaction.resp.lm = readRDS("./host_proteomics/proteomics.data.interaction.resp.lm.Rds")

# add mito linker
colnames(proteomics.data.interaction.resp.lm)[colnames(proteomics.data.interaction.resp.lm) == "protein"] = "protein.x"
proteomics.data.interaction.resp.lm = merge(proteomics.data.interaction.resp.lm,
                                            proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.interaction.resp.lm.mito = merge(proteomics.data.interaction.resp.lm,
                                                 mitocarta.pathways.expanded,
                                                 by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""))
# arrange by significance
proteomics.data.interaction.resp.lm.mito %>% arrange(pval)

# :: MPX Volcano: Responders (per location)  ----------------------------------------------------------

# volcano plots
proteomics.data.interaction.resp.lm.mito.volcano = ggplot(proteomics.data.interaction.resp.lm.mito[,c("coef", "pval", "padj", "protein", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                          aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=coef))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", protein), NA)),
                           size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Interaction Coefficient",
       y="Unadjusted p value")+
  facet_wrap(~location)
proteomics.data.interaction.resp.lm.mito.volcano

# :: MLI ------------------------------------------------------------------


# :: Question 1: Baseline -------------------------------------------------

# Question: are there significant differences between N and A in groups
# Solid question to ask

proteomics.metadata.baseline = subset(proteomics.metadata, Time == "0")

# loop through locations and perform LMER
proteomics.data.baseline = proteomics.data.clean[rownames(proteomics.data.clean) %in% proteomics.metadata.baseline$code,]


if(rerun==T){
  
  t1 = Sys.time()
  proteomics.data.baseline.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(proteomics.metadata.baseline, Location == location)
    data.subset = proteomics.data.baseline[rownames(proteomics.data.baseline) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, )
    # apply filtering per location (present in 80% of samples)
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, lapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      # note: this isn't necessary with 80% prevalence filter
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = (log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        lmer.results = lmerTest::lmer(protein ~ Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # build separate model for beta-coef
        lmer.results.scaled = lmerTest::lmer(scale(protein) ~ Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          log2fc = lmer.results[rownames(lmer.results) == "StatusA",]$Estimate,
          coef = lmer.results.scaled[rownames(lmer.results.scaled) == "StatusA",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "StatusA",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # ~20 min
  # calculate padj
  proteomics.data.baseline.lm = proteomics.data.baseline.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.baseline.lm, "./host_proteomics/proteomics.data.baseline.lm.Rds")
}
proteomics.data.baseline.lm = readRDS("./host_proteomics/proteomics.data.baseline.lm.Rds")

# add mito linker
colnames(proteomics.data.baseline.lm)[colnames(proteomics.data.baseline.lm) == "protein"] = "protein.x"
proteomics.data.baseline.lm = merge(proteomics.data.baseline.lm,
                                    proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.baseline.lm.mito = merge(proteomics.data.baseline.lm,
                                         mitocarta.pathways.expanded,
                                         by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""))
# arrange by significance
proteomics.data.baseline.lm.mito %>% arrange(padj)

# volcano plots
proteomics.data.baseline.lm.mito.volcano = ggplot(proteomics.data.baseline.lm.mito[,c("coef", "pval", "padj","mito", "protein", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                  aes(x=as.numeric(coef), y=-log10(padj)))+
  geom_point(aes(fill=coef, shape=mito))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, protein, NA)), size=3)+
  geom_hline(yintercept=-log10(0.05), linetype=2, color="black", alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("black", "red"))+
  scale_shape_manual(values=c(21,24))+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient")+ # note: Log2FC forms clusters of values; scaled is smooth
  facet_wrap(~location)
proteomics.data.baseline.lm.mito.volcano

# many proteins have padj < 0.05

# odd shape; I suspect on account of the normalization strategy; can this be accounted for in the lm?
proteomics.data.baseline.lm.mito$MitoPathway %>% unique()

subset(proteomics.data.baseline.lm.mito, mito=="mito" & padj < 0.20) %>% arrange(padj)

# many proteins make sense and can be commented on (e.g. SLC's, ZG16, HSP's)


# :: Question 2: RS -------------------------------------------------------


ggplot(subset(proteomics.metadata, Time %in% c(0,1)) %>%
         subset(code %in% proteomics.clean.pca.df$code) %>%
         mutate(Location = factor(Location, levels=c("TI", "PC", "DC"))),
       aes(x=as.factor(Time), y=HM))+
  ggbeeswarm::geom_beeswarm(shape=21, size=3,aes(fill=Status), cex=2)+
  facet_wrap(~Location)+
  theme_minimal()

# Question: are there significant differences before and after RS
# control for inflammation status

proteomics.metadata.rs = subset(proteomics.metadata, Group == "RS")

# loop through locations and perform LMER
proteomics.data.rs = proteomics.data.clean[rownames(proteomics.data.clean) %in% proteomics.metadata.rs$code,]

if(rerun==T){
  t1 = Sys.time()
  proteomics.data.rs.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(proteomics.metadata.rs, Location == location)
    data.subset = proteomics.data.rs[rownames(proteomics.data.rs) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = (log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(protein ~ Status + Time + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time() # 10 min with parallel
  t2 - t1
  # calculate padj
  proteomics.data.rs.lm = proteomics.data.rs.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Baseline, controlling for inflammation status (so, upregulated or downregulated at scope)  
  
  saveRDS(proteomics.data.rs.lm, "./host_proteomics/proteomics.data.rs.lm.Rds")
}
proteomics.data.rs.lm = readRDS("./host_proteomics/proteomics.data.rs.lm.Rds")


# add mito linker
colnames(proteomics.data.rs.lm)[colnames(proteomics.data.rs.lm) == "protein"] = "protein.x"
proteomics.data.rs.lm = merge(proteomics.data.rs.lm,
                              proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.rs.lm.mito = merge(proteomics.data.rs.lm,
                                   mitocarta.pathways.expanded,
                                   by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""))
# arrange by significance
proteomics.data.rs.lm.mito %>% arrange(padj)

# volcano plots
proteomics.data.rs.lm.mito.volcano = ggplot(proteomics.data.rs.lm.mito[,c("coef", "pval", "padj", "protein","mito", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                            aes(x=as.numeric(coef), y=-log10(padj)))+
  geom_point(aes(fill=coef, shape=mito))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, protein, NA)), size=3)+
  geom_hline(yintercept=-log10(0.05), linetype=2, color="black", alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("black", "red"))+
  scale_shape_manual(values=c(21,24))+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Log2FC")+
  facet_wrap(~location)
proteomics.data.rs.lm.mito.volcano

# odd shape; I suspect on account of the normalization strategy; can this be accounted for in the lm?
proteomics.data.rs.lm.mito$MitoPathway %>% unique()

subset(proteomics.data.rs.lm.mito, mito=="mito" & padj < 0.20) %>% arrange(padj)

# After controlling for inflammation status, the RS group had these sig impacted features compared to baseline

# :: Question 3: Placebo -------------------------------------------------------

# Question: are there significant differences before and after Placebo

proteomics.metadata.placebo = subset(proteomics.metadata, Group == "Placebo")

# loop through locations and perform LMER
proteomics.data.placebo = proteomics.data.clean[rownames(proteomics.data.clean) %in% proteomics.metadata.placebo$code,]

if(rerun==T){
  t1 = Sys.time()
  proteomics.data.placebo.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(proteomics.metadata.placebo, Location == location)
    data.subset = proteomics.data.placebo[rownames(proteomics.data.placebo) %in% metadata.subset$code,]
    # keep intersecting samples
    data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
    metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
    # apply filtering per location
    proteins.to.keep = data.subset
    proteins.to.keep[proteins.to.keep != 0] <- 1
    proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
    # loop through lmer
    do.call(rbind, parallel::mclapply(proteins.to.keep, function(protein){
      print(paste0(location, " ", protein, " ", round(which(proteins.to.keep == protein) / length(proteins.to.keep) * 100, digits=4)))
      # protein = proteins.to.keep[1]
      # add protein abundance data by merging
      data.subset$code = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset,
                              data.subset[,colnames(data.subset) %in% c("code", protein)])
      colnames(metadata.subset)[colnames(metadata.subset) == protein] = "protein"
      
      if(sum(metadata.subset$protein != (proteomics.pseudocount)) <= 3){
        data.frame(
          location = location,
          feature = protein,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$protein = scale(log2(metadata.subset$protein+proteomics.pseudocount)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(protein ~ Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          protein = protein,
          coef = lmer.results[rownames(lmer.results) == "Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 15 min with parallel
  # calculate padj
  proteomics.data.placebo.lm = proteomics.data.placebo.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(proteomics.data.placebo.lm, "./host_proteomics/proteomics.data.placebo.lm.Rds")
}
proteomics.data.placebo.lm = readRDS("./host_proteomics/proteomics.data.placebo.lm.Rds")

# add mito linker
colnames(proteomics.data.placebo.lm)[colnames(proteomics.data.placebo.lm) == "protein"] = "protein.x"
proteomics.data.placebo.lm = merge(proteomics.data.placebo.lm,
                                   proteomics.gene.map, by="protein.x", all.x=T)

# add mitochondrial labels
proteomics.data.placebo.lm.mito = merge(proteomics.data.placebo.lm,
                                        mitocarta.pathways.expanded,
                                        by="gene", all.x=T) %>%
  mutate(mito = ifelse(gene %in% mito.proteins, "mito", ""))
# arrange by significance
proteomics.data.placebo.lm.mito %>% arrange(padj)

# volcano plots
proteomics.data.placebo.lm.mito.volcano = ggplot(proteomics.data.placebo.lm.mito[,c("coef", "pval", "padj","mito", "protein", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                 aes(x=as.numeric(coef), y=-log10(pval)))+
  geom_point(shape = 21, aes(fill=coef))+
  ggrepel::geom_text_repel(aes(label = ifelse(padj < 0.05, protein, NA),
                               color = ifelse(mito == "mito", "red", "black")), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("black", "red"))+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient")+
  facet_wrap(~location)
proteomics.data.placebo.lm.mito.volcano

# odd shape; I suspect on account of the normalization strategy; can this be accounted for in the lm?
proteomics.data.placebo.lm.mito$MitoPathway %>% unique()

subset(proteomics.data.placebo.lm.mito, mito=="mito" & padj < 0.20) %>% arrange(padj)
# only Succinate Dehydrogenase in DC

# note: fewer sig features with Placebo likely because of lower power than RS group



# :: ASV Interaction (per location) -------------------------------------------------------

# Question: are there significant differences before and after RS relative to Placebo

# loop through locations and perform LMER
asv.data.rs.interaction = lsarp.mli.ps.mat[rownames(lsarp.mli.ps.mat) %in% lsarp.mli.ps.meta$Sample,]

if(rerun==T){
  lsarp.pseudo = min(asv.data.rs.interaction[asv.data.rs.interaction != 0])/2
  
  t1 = Sys.time()
  asv.data.rs.interaction.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(lsarp.mli.ps.meta, Location == location)
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, make.names(Sample) %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = make.names(Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(asv ~ Group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "GroupRS:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1
  # calculate padj
  asv.data.rs.interaction.lm = asv.data.rs.interaction.lm %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.rs.interaction.lm, "./host_proteomics/asv.data.interaction.lm.Rds")
}
asv.data.rs.interaction.lm = readRDS("./host_proteomics/asv.data.interaction.lm.Rds")


# :: ASV Volcano (per location) ----------------------------------------------------------

# volcano plots
asv.data.rs.interaction.lm.volcano = ggplot(asv.data.rs.interaction.lm[,c("coef", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                            aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=coef))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Interaction Coefficient",
       y="Unadjusted p value")+
  facet_wrap(~location)
asv.data.rs.interaction.lm.volcano

# :: ASV Interactions: Responders (per location) -----------------------------------------------------------

# goal: conduct analysis using Response/Non-Response paradigm


if(rerun==T){
  t1 = Sys.time()
  asv.data.rs.interaction.lm.resp = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
    # subset to location
    metadata.subset = subset(lsarp.mli.ps.meta, Location == location & Group == "RS")
    # add Flare.group
    metadata.subset = merge(metadata.subset,
                            lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(), by="HM")
    
    data.subset = asv.data.rs.interaction[rownames(asv.data.rs.interaction) %in% metadata.subset$Sample,]
    # keep intersecting samples
    data.subset = data.subset[metadata.subset$Sample,]  %>% data.frame()
    metadata.subset = subset(metadata.subset, make.names(Sample) %in% rownames(data.subset))
    # apply filtering per location
    asvs.to.keep = data.subset
    asvs.to.keep[asvs.to.keep != 0] <- 1
    asvs.to.keep = colnames(asvs.to.keep[,colSums(asvs.to.keep) > nrow(asvs.to.keep)*0.2])
    # loop through lmer
    do.call(rbind, lapply(asvs.to.keep, function(asv){
      print(paste0(location, " ", asv, " ", round(which(asvs.to.keep == asv) / length(asvs.to.keep) * 100, digits=4)))
      # asv = asvs.to.keep[1]
      # add asv abundance data by merging
      data.subset$Sample = rownames(data.subset)
      
      metadata.subset = merge(metadata.subset %>% mutate(Sample = make.names(Sample)),
                              data.subset[,colnames(data.subset) %in% c("Sample", asv)]) %>% distinct()
      colnames(metadata.subset)[colnames(metadata.subset) == asv] = "asv"
      
      if(sum(metadata.subset$asv != (0)) <= nrow(metadata.subset)*0.2){
        data.frame(
          location = location,
          feature = asv,
          coef = NA,
          pval = NA)
      }else{
        metadata.subset$asv = (log2(metadata.subset$asv+lsarp.pseudo)) 
        metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
        metadata.subset$flare.group = factor(metadata.subset$flare.group, levels=c("Relapse", "Remit"))
        metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
        lmer.results = lmerTest::lmer(asv ~ flare.group*Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
        # extract data
        data.frame(
          location = location,
          feature = asv,
          coef = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$Estimate,
          pval = lmer.results[rownames(lmer.results) == "flare.groupRemit:Time1",]$`Pr...t..`)
      }
    }))}))
  t2 = Sys.time()
  t2 - t1 # 1 min
  # calculate padj
  asv.data.rs.interaction.lm.resp = asv.data.rs.interaction.lm.resp %>% 
    mutate(padj = p.adjust(pval, method="BH")) %>% 
    mutate(sig = ifelse(padj < 0.05, "***",ifelse(padj < 0.20, "*",  "")))%>%
    arrange(pval)
  # Note: coefficient is relative to Not-inflamed (so, upregulated or downregulated in inflamed)  
  
  saveRDS(asv.data.rs.interaction.lm.resp, "./host_proteomics/asv.data.interaction.resp.lm.Rds")
}
asv.data.rs.interaction.lm.resp = readRDS("./host_proteomics/asv.data.interaction.resp.lm.Rds")


# :: ASV Volcano: Responders (per location) ----------------------------------------------------------

# volcano plots
asv.data.rs.interaction.lm.resp.volcano = ggplot(asv.data.rs.interaction.lm.resp[,c("coef", "pval", "padj", "feature", "location")] %>% distinct() %>% mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                                 aes(x=as.numeric(coef), y=pval))+
  geom_point(shape = 21, aes(fill=coef))+
  geom_hline(yintercept=0.05, color="black", alpha=0.5, linetype=2)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Interaction Coefficient",
       y="Unadjusted p value")+
  facet_wrap(~location)
asv.data.rs.interaction.lm.resp.volcano

# :: 16S Plots ---------------------------------------------------------

# Plot heatmaps of padj < 0.20

# first, lists of mitochondrial features to include; per location
padj20.features.ti.asv = c(unique(subset(asv.data.baseline.lm, pval < 0.05 & location == "TI")$feature),
                           unique(subset(asv.data.rs.lm, pval < 0.05 & location == "TI")$feature),
                           unique(subset(asv.data.placebo.lm, pval < 0.05 & location == "TI")$feature))
padj20.features.pc.asv = c(unique(subset(asv.data.baseline.lm, pval < 0.05 & location == "PC")$feature),
                           unique(subset(asv.data.rs.lm, pval < 0.05 & location == "PC")$feature),
                           unique(subset(asv.data.placebo.lm, pval < 0.05 & location == "PC")$feature))
padj20.features.dc.asv = c(unique(subset(asv.data.baseline.lm, pval < 0.05 & location == "DC")$feature),
                           unique(subset(asv.data.rs.lm, pval < 0.05 & location == "DC")$feature),
                           unique(subset(asv.data.placebo.lm, pval < 0.05 & location == "DC")$feature))

# next, merge dataframes with LFC and pval per location (label comparisons, too)
asv.data.lm.ti = rbind(
  asv.data.baseline.lm %>% mutate(condition = "Baseline"),
  asv.data.rs.lm %>% mutate(condition = "RS"),
  asv.data.placebo.lm %>% mutate(condition = "Placebo")
) %>% subset(feature %in% padj20.features.ti.asv & location == "TI") 
asv.data.lm.pc = rbind(
  asv.data.baseline.lm %>% mutate(condition = "Baseline"),
  asv.data.rs.lm %>% mutate(condition = "RS"),
  asv.data.placebo.lm %>% mutate(condition = "Placebo")
) %>% subset(feature %in% padj20.features.pc.asv & location == "PC")
asv.data.lm.dc = rbind(
  asv.data.baseline.lm %>% mutate(condition = "Baseline"),
  asv.data.rs.lm %>% mutate(condition = "RS"),
  asv.data.placebo.lm %>% mutate(condition = "Placebo")
) %>% subset(feature %in% padj20.features.dc.asv & location == "DC") 

# merge all
asv.data.lm.all = rbind(asv.data.lm.ti,
                        asv.data.lm.pc,
                        asv.data.lm.dc) %>% data.frame()

# now plot
asv.data.lm.ti.plot = ggplot(asv.data.lm.ti %>%
                               # reorder proteins based on baseline coef
                               mutate(feature.order = as.numeric(factor(feature, levels=arrange(subset(asv.data.lm.ti, condition == "Baseline"), coef)$feature))),
                             aes(x=condition, y=reorder(feature, feature.order)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(pval < 0.05, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")
asv.data.lm.pc.plot = ggplot(asv.data.lm.pc %>%
                               # reorder proteins based on baseline coef
                               mutate(feature.order = as.numeric(factor(feature, levels=arrange(subset(asv.data.lm.pc, condition == "Baseline"), coef)$feature))),
                             aes(x=condition, y=reorder(feature, feature.order)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(pval < 0.05, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")
asv.data.lm.dc.plot = ggplot(asv.data.lm.dc %>%
                               # reorder proteins based on baseline coef
                               mutate(feature.order = as.numeric(factor(feature, levels=arrange(subset(asv.data.lm.dc, condition == "Baseline"), coef)$feature))),
                             aes(x=condition, y=reorder(feature, feature.order)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(pval < 0.05, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")

asv.data.lm.ti.plot+asv.data.lm.pc.plot+asv.data.lm.dc.plot


# nothing is significant after FDR



# >> MitoInflammation Signature -------------------------------------------

# Goal: Calculate an inflammation score (use lasso) and see if RS decreases relative to Placebo

# train using baseline samples (inflamed vs uninflamed)
# test on ALL samples (therefore, baseline will be overfit)

t1 = Sys.time()
proteomics.data.lasso.results = lapply(c("TI", "PC", "DC"), function(location){
  # subset to location at time 0 (to train inflammation score)
  metadata.subset = subset(proteomics.metadata.baseline, Location == location)
  data.subset = proteomics.data.baseline[rownames(proteomics.data.baseline) %in% metadata.subset$code,] %>% data.frame()
  
  ### subset to MitoProteins
  data.subset = data.subset[,colnames(data.subset) %in% 
                              # subset the mapping file to just mito proteins, then filter to make.names proteins
                              subset(proteomics.gene.map, gene %in% mito.proteins)$protein.x]
  ###
  ## NOTE: Results are not significant if you use ALL proteins; only with Mito Proteins
  
  # keep intersecting samples
  data.subset = data.subset[rownames(data.subset) %in% metadata.subset$code,] %>% data.frame()
  metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
  # create data structure for lasso
  lasso.data = data.subset
  lasso.data$code = rownames(lasso.data)
  lasso.data = merge(lasso.data, 
                     subset(metadata.subset, Time == "0")[,c("code", "Status")] %>% distinct(), by="code")
  rownames(lasso.data) = lasso.data$code
  lasso.data$code = NULL
  
  #
  # prepare downstream (testing) datasets, and ensure all have the same variables
  proteomics.data.rs.lasso = proteomics.data.rs[rownames(proteomics.data.rs) %in% subset(proteomics.metadata.rs, Time %in% c("0", "1") & Location == location)$code,]
  proteomics.data.placebo.lasso = proteomics.data.placebo[rownames(proteomics.data.placebo) %in% subset(proteomics.metadata.placebo, Time %in% c("0", "1") & Location == location)$code,]
  
  proteomics.data.rs.lasso = proteomics.data.rs.lasso[,colnames(proteomics.data.rs.lasso) %in% colnames(lasso.data)]
  proteomics.data.placebo.lasso = proteomics.data.placebo.lasso[,colnames(proteomics.data.placebo.lasso) %in% colnames(lasso.data)]
  
  # for training data, keep features present in the testing data
  lasso.data = lasso.data[,colnames(lasso.data) %in% c(colnames(proteomics.data.rs.lasso), colnames(proteomics.data.placebo.lasso), "Status")]
  
  # build caret lasso
  table(lasso.data$Status) # class imbalance --> use ROC
  set.seed(25)
  library("caret")
  lasso_model <- caret::train(
    x = as.matrix(lasso.data[,colnames(lasso.data) != "Status"]),
    y = as.factor(lasso.data$Status),
    method = "glmnet",
    metric = "ROC",
    trControl = caret::trainControl(method = "cv", number = 5,
                                    summaryFunction = twoClassSummary,   # enables AUC
                                    classProbs = TRUE,                   # needed for ROC
                                    savePredictions = "final"),
    tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-4, 1, length = 100)),
    preProcess = c("center", "scale")
  )
  # apply model to RS and Placebo Time 2
  proteomics.data.rs.lasso.preds = data.frame(
    Group = "RS",
    code = rownames(proteomics.data.rs.lasso),
    prob = predict(lasso_model, proteomics.data.rs.lasso, type="prob")$A)
  proteomics.data.placebo.lasso.preds = data.frame(
    Group = "Placebo",
    code = rownames(proteomics.data.placebo.lasso),
    prob = predict(lasso_model, proteomics.data.placebo.lasso, type="prob")$A)
  # merge
  proteomics.data.lasso.preds = rbind(proteomics.data.rs.lasso.preds,
                                      proteomics.data.placebo.lasso.preds) %>% data.frame() %>%
    tidyr::separate(code, into=c("HM", "Time", "Location", "Status"))
  # record AUC of predictor
  proteomics.data.lasso.preds$AUC = subset(lasso_model$results, lambda == lasso_model$bestTune$lambda)$ROC
  
  ### extract feature coefficients
  proteomics.data.lasso.coefs = coef(lasso_model$finalModel, s =  lasso_model$bestTune$lambda) %>% as.matrix() %>% as.data.frame() %>% subset(s0 != 0)
  proteomics.data.lasso.coefs = data.frame(protein = rownames(proteomics.data.lasso.coefs)[-1],
                                           coef = proteomics.data.lasso.coefs[-1,]) %>%
    mutate(location = location)
  ###
  
  # return
  return(list(proteomics.data.lasso.preds,
              proteomics.data.lasso.coefs))
})
t2 = Sys.time()
t2 - t1

# extract predictions/scores
proteomics.data.lasso.scores = do.call(rbind, lapply(c(1:3), function(x){
  purrr::pluck(proteomics.data.lasso.results[[x]][1])[[1]]
}))

# Are the scores (predicted probability of being inflamed) lower in RS?
mito.lasso.ti = lmerTest::lmer(prob ~ Group  + (1|HM), proteomics.data.lasso.scores %>%
                                 subset(Location == "TI" & Time == 1)) %>% summary() %>% coef()
mito.lasso.pc = lmerTest::lmer(prob ~ Group   + (1|HM), proteomics.data.lasso.scores %>%
                                 subset(Location == "PC" & Time == 1)) %>% summary()%>% coef()
mito.lasso.dc = lmerTest::lmer(prob ~ Group   + (1|HM), proteomics.data.lasso.scores %>%
                                 subset(Location == "DC" & Time == 1)) %>% summary()%>% coef()

# or wilcoxon
mito.lasso.ti.tt = wilcox.test(subset(proteomics.data.lasso.scores, Location == "TI" & Time == 1 & Group == "RS")$prob,
                               subset(proteomics.data.lasso.scores, Location == "TI" & Time == 1 & Group == "Placebo")$prob)
mito.lasso.pc.tt = wilcox.test(subset(proteomics.data.lasso.scores, Location == "PC" & Time == 1 & Group == "RS")$prob,
                               subset(proteomics.data.lasso.scores, Location == "PC" & Time == 1 & Group == "Placebo")$prob)
mito.lasso.dc.tt = wilcox.test(subset(proteomics.data.lasso.scores, Location == "DC" & Time == 1 & Group == "RS")$prob,
                               subset(proteomics.data.lasso.scores, Location == "DC" & Time == 1 & Group == "Placebo")$prob)

mito.lasso = data.frame(Location = c("TI", "PC", "DC"),
                        pval = c(mito.lasso.ti[2,5], mito.lasso.pc[2,5], mito.lasso.dc[2,5])) %>%
  mutate(padj = p.adjust(pval, method="bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  mutate(Location = factor(Location, levels=c("TI", "PC","DC")))


# :: MitoInflammation Plot ------------------------------------------------

proteomics.data.lasso.scores.plot = ggplot(proteomics.data.lasso.scores %>%
                                             subset(Time == 1) %>%
                                             mutate(Location = factor(Location, levels=c("TI", "PC","DC"))),
                                           aes(x=Group, y=prob))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill = Status), size=2.5)+
  geom_text(data=mito.lasso, aes(x=1.5, y=Inf, label=paste("p = ", round(pval, digits=3))), vjust=2, size=4)+
  theme_classic()+theme(#legend.position="none",
    strip.text = element_text(size=12),
    strip.background = element_rect(
      color="black"))+
  labs(x="", y="Mitochondrial Inflammation Score", fill="Status")+
  facet_wrap(~Location)
proteomics.data.lasso.scores.plot
# RS-treated have a lower mitochondrial inflammation score 
# compared to placebo, in DC, but not TI or PC

# Interpret the features driving this score

# extract coefs
proteomics.data.lasso.coefs = do.call(rbind, lapply(c(1:3), function(x){
  purrr::pluck(proteomics.data.lasso.results[[x]][2])[[1]]
}))
proteomics.data.lasso.coefs.plot = ggplot(proteomics.data.lasso.coefs %>%
                                            mutate(location = factor(location, levels=c("TI", "PC", "DC"))),
                                          aes(x=location, y=reorder(protein, coef)))+
  geom_tile(aes(fill=coef))+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(axis.text.x=element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="", fill="Lasso\nCoefficient")+
  facet_wrap(~location, scales="free_x")
proteomics.data.lasso.coefs.plot

# plot both
proteomics.data.lasso.scores.plot + 
  proteomics.data.lasso.coefs.plot+
  patchwork::plot_layout(widths=c(3,1))

# HMGCL == inflammation (ketone body production); 3-hydroxy-3-methylglutaryl-CoA (HMG-CoA) into acetyl-CoA and acetoacetate
# PDP1 != inflammation (?); pyruvate into acetyl-CoA


# :: WGCNA ----------------------------------------------------------------

# identify protein clusters in RS and placebo groups
# see where the MitoProteins are situated

# follow tutorial: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#Preliminaries
library("WGCNA")

# clean proteomic.data
proteomics.data.clean = proteomics.data[proteomics.metadata$code,]

# identify potential outliers
gsg <-goodSamplesGenes(proteomics.data.clean)
gsg$allOK
# all ok

# pick soft threshold for weighted adjacency matrix
spt <- pickSoftThreshold(proteomics.data.clean) 
# plot soft thresholds
ggplot(spt$fitIndices %>% data.frame(),
       aes(x=Power, y=SFT.R.sq))+
  geom_text(aes(label = Power))+
  theme_minimal()
# plot mean connectivity
ggplot(spt$fitIndices %>% data.frame(),
       aes(x=Power, y=mean.k.))+
  geom_text(aes(label = Power))+
  theme_minimal()
# correlate the two
ggplot(spt$fitIndices %>% data.frame(),
       aes(x=mean.k., y=SFT.R.sq))+
  geom_text(aes(label = Power))+
  theme_minimal()
# pick 3; Highest R2 and highest Mean Connectivity
softPower <- 3
adjacency <- adjacency(proteomics.data.clean, power = softPower)
# create modules:
# first, similarity matrix
TOM <- TOMsimilarity(adjacency)
# into dissimilarity matrix
TOM.dissimilarity <- 1-TOM
# heirarchical clustering
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
# cut tree (per WGCNA's suggestion)
Modules <- cutreeDynamic(dendro = geneTree, 
                         distM = TOM.dissimilarity, 
                         deepSplit = 2, pamRespectsDendro = FALSE, 
                         minClusterSize = 30)
# check table
table(Modules)
# add colors
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors)
# plot dendrogram
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# eigengene
MElist <- moduleEigengenes(proteomics.data.clean, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)
# merge similar modules
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.5, col = "red")
merge <- mergeCloseModules(proteomics.data.clean, ModuleColors, cutHeight = .5)
# fix colors
mergedColors = merge$colors
mergedMEs = merge$newMEs

# assign trait data
allTraits <- proteomics.metadata[,c("study_id", "Time", "Location", "Status", "HM", "code", "Group","flare.call_verified", "days.taken")] %>% distinct()
rownames(allTraits) = allTraits$code
# match rows
Samples <- rownames(proteomics.data.clean)
traitRows <- match(Samples, allTraits$code)
datTraits <- allTraits[traitRows, c("Time", "Location", "Status", "Group", "flare.call_verified")]
# correlate with clinical data
nGenes = ncol(proteomics.data.clean)
nSamples = nrow(proteomics.data.clean)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# check module membership of genes of interest
Modules[which(colnames(proteomics.data.clean) == "HMGCL")] # 2
Modules[which(colnames(proteomics.data.clean) == "PDP1")] # 4

# stats on modules
mergedMEs


t1 = Sys.time()
proteomics.data.me.rs.lm = do.call(rbind, lapply(c("TI", "PC", "DC"), function(location){
  # subset to location
  metadata.subset = subset(proteomics.metadata.rs, Location == location & Group == "RS")
  data.subset = mergedMEs[rownames(mergedMEs) %in% metadata.subset$code,]
  # keep intersecting samples
  data.subset = data.subset[metadata.subset$code,] %>% data.frame()
  metadata.subset = subset(metadata.subset, code %in% rownames(data.subset))
  # apply filtering per location
  #proteins.to.keep = data.subset
  #proteins.to.keep[proteins.to.keep != 0] <- 1
  #proteins.to.keep = colnames(proteins.to.keep[,colSums(proteins.to.keep) > nrow(proteins.to.keep)*0.8])
  # loop through lmer
  do.call(rbind, lapply(colnames(mergedMEs), function(me){
    print(paste0(location, " ", me, " ", round(which(colnames(mergedMEs) == me) / length(colnames(mergedMEs)) * 100, digits=4)))
    # protein = proteins.to.keep[1]
    # add protein abundance data by merging
    data.subset$code = rownames(data.subset)
    
    metadata.subset = merge(metadata.subset,
                            data.subset[,colnames(data.subset) %in% c("code", me)])
    colnames(metadata.subset)[colnames(metadata.subset) == me] = "me"
    
    metadata.subset$me = scale(metadata.subset$me)
    metadata.subset$Status = factor(metadata.subset$Status, levels=c("N", "A"))
    metadata.subset$Time = factor(metadata.subset$Time, levels=c("0", "1"))
    lmer.results = lmerTest::lmer(me ~ Time + Status + (1|study_id), metadata.subset) %>% summary() %>% coef() %>% data.frame()
    # extract data
    data.frame(
      location = location,
      me = me,
      coef = lmer.results[rownames(lmer.results) == "Time1",]$Estimate,
      pval = lmer.results[rownames(lmer.results) == "Time1",]$`Pr...t..`)
    
  }))}))
t2 = Sys.time() # ~20 min
t2 - t1
proteomics.data.me.rs.lm %>% mutate(padj = p.adjust(pval, method="BH")) %>% arrange(padj)


# :: STRING ---------------------------------------------------------------

library("STRINGdb")

string_db <- STRINGdb$new(
  version = "11.5",       # Or latest
  species = 9606,         # Human
  score_threshold = 900,  # Very high confidence
  input_directory = ""
)

# merge all data
proteomics.significant = rbind(proteomics.data.baseline.lm.mito %>% mutate(comparison = "baseline") %>% subset(padj < 0.20),
                               proteomics.data.rs.lm.mito %>% mutate(comparison = "RS") %>% subset(padj < 0.20),
                               proteomics.data.placebo.lm.mito %>% mutate(comparison = "Placebo") %>% subset(padj < 0.20)) %>% data.frame()

mapped <- string_db$map(proteomics.significant, "protein", removeUnmappedRows = TRUE)

# plot only connected proteins
ppi_df <- string_db$get_interactions(mapped$STRING_id)

# Step 2: Identify STRING IDs with at least one interaction
connected_ids <- unique(c(ppi_df$from, ppi_df$to))

# Step 3: Keep only mapped proteins that are connected
mapped_connected <- mapped[mapped$STRING_id %in% connected_ids, ]

# Step 4: Plot only connected proteins
string_db$plot_network(mapped_connected$STRING_id)

# work with igraph
ppi_df$protein1 = mapped_connected$protein[match(ppi_df$from, mapped_connected$STRING_id)]
ppi_df$protein2 = mapped_connected$protein[match(ppi_df$to, mapped_connected$STRING_id)]

subset(ppi_df, protein1 == "HMGCL")$protein2 %>% unique()
# AUH = leucine degradation
# ECI2 = beta-oxidation of unsaturated fatty acids
# BDH1 = interconversion of acetoacetate and (R)-3-hydroxybutyrate
# NUDT19 = Fatty acyl-coenzyme A (CoA) diphosphatase 
# HSD17B4 = Peroxisomal fatty acid beta-oxidation
subset(proteomics.data.rs.lm.mito, protein == "AUH") # not sig
subset(proteomics.data.rs.lm.mito, protein == "ECI2") # not sig
subset(proteomics.data.rs.lm.mito, protein == "BDH1") # not sig
subset(proteomics.data.rs.lm.mito, protein == "NUDT19") # not sig
subset(proteomics.data.rs.lm.mito, protein == "HSD17B4") # not sig
subset(proteomics.data.rs.lm.mito, protein == "HMGCL") # not sig


subset(ppi_df, protein1 == "PDP1")$protein2 %>% unique()


# :: MPX Plots ---------------------------------------------------------

# Plot heatmaps of padj < 0.20

# first, lists of mitochondrial features to include; per location
padj20.features.ti = c(unique(subset(proteomics.data.baseline.lm.mito, padj < 0.20 & location == "TI" & mito == "mito")$protein),
                       unique(subset(proteomics.data.rs.lm.mito, padj < 0.20 & location == "TI" & mito == "mito")$protein),
                       unique(subset(proteomics.data.placebo.lm.mito, padj < 0.20 & location == "TI" & mito == "mito")$protein))
padj20.features.pc = c(unique(subset(proteomics.data.baseline.lm.mito, padj < 0.20 & location == "PC" & mito == "mito")$protein),
                       unique(subset(proteomics.data.rs.lm.mito, padj < 0.20 & location == "PC" & mito == "mito")$protein),
                       unique(subset(proteomics.data.placebo.lm.mito, padj < 0.20 & location == "PC" & mito == "mito")$protein))
padj20.features.dc = c(unique(subset(proteomics.data.baseline.lm.mito, padj < 0.20 & location == "DC" & mito == "mito")$protein),
                       unique(subset(proteomics.data.rs.lm.mito, padj < 0.20 & location == "DC" & mito == "mito")$protein),
                       unique(subset(proteomics.data.placebo.lm.mito, padj < 0.20 & location == "DC" & mito == "mito")$protein))

# next, merge dataframes with LFC and pval per location (label comparisons, too)
proteomics.data.lm.mito.ti = rbind(
  proteomics.data.baseline.lm.mito %>% mutate(condition = "Baseline"),
  proteomics.data.rs.lm.mito %>% mutate(condition = "RS"),
  proteomics.data.placebo.lm.mito %>% mutate(condition = "Placebo")
) %>% subset(protein %in% padj20.features.ti & location == "TI") %>%
  dplyr::select(-MitoPathway, -MitoPathways.Hierarchy) %>% distinct()

proteomics.data.lm.mito.pc = rbind(
  proteomics.data.baseline.lm.mito %>% mutate(condition = "Baseline"),
  proteomics.data.rs.lm.mito %>% mutate(condition = "RS"),
  proteomics.data.placebo.lm.mito %>% mutate(condition = "Placebo")
) %>% subset(protein %in% padj20.features.pc & location == "PC")  %>%
  dplyr::select(-MitoPathway, -MitoPathways.Hierarchy) %>% distinct()

proteomics.data.lm.mito.dc = rbind(
  proteomics.data.baseline.lm.mito %>% mutate(condition = "Baseline"),
  proteomics.data.rs.lm.mito %>% mutate(condition = "RS"),
  proteomics.data.placebo.lm.mito %>% mutate(condition = "Placebo")
) %>% subset(protein %in% padj20.features.dc & location == "DC")  %>%
  dplyr::select(-MitoPathway, -MitoPathways.Hierarchy) %>% distinct()

# merge all
proteomics.data.lm.mito.all = rbind(proteomics.data.lm.mito.ti,
                                    proteomics.data.lm.mito.pc,
                                    proteomics.data.lm.mito.dc) %>% data.frame()

# now plot
proteomics.data.lm.mito.ti.plot = ggplot(proteomics.data.lm.mito.ti %>%
                                           # reorder proteins based on baseline coef
                                           mutate(protein = factor(protein, levels=arrange(subset(proteomics.data.lm.mito.ti, condition == "Baseline"), coef)$protein)),
                                         aes(x=condition, y=(protein)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(padj < 0.20, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")
proteomics.data.lm.mito.pc.plot = ggplot(proteomics.data.lm.mito.pc %>%
                                           # reorder proteins based on baseline coef
                                           mutate(protein = factor(protein, levels=arrange(subset(proteomics.data.lm.mito.pc, condition == "Baseline"), coef)$protein)),
                                         aes(x=condition, y=(protein)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(padj < 0.20, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")
proteomics.data.lm.mito.dc.plot = ggplot(proteomics.data.lm.mito.dc %>%
                                           # reorder proteins based on baseline coef
                                           mutate(protein = factor(protein, levels=arrange(subset(proteomics.data.lm.mito.dc, condition == "Baseline"), coef)$protein)),
                                         aes(x=condition, y=(protein)))+
  geom_tile(aes(fill=coef), color="black")+
  geom_text(aes(label = ifelse(padj < 0.05, "***", NA)), size=6, vjust=0.75)+
  geom_text(aes(label = ifelse(padj < 0.20, "*", NA)), size=6, vjust=0.75)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="")+
  facet_wrap(~location, scales="free")

proteomics.data.lm.mito.ti.plot+proteomics.data.lm.mito.pc.plot+proteomics.data.lm.mito.dc.plot


