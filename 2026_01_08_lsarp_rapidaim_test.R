### 2025_09_20  LSARP RapidAIM Theory

# Background: RapidAIM-RS Selection Algorithm worked for some, but not others
# i.e. some participants on RS clinically responded, others did not
# Hypothesis: Wrong bacteria contributed to the butyrogen score in nonresponders

# Experiment: 
# 1. Re-run 16S processing, compare differences in Log2FC of selected RS vs PBS
# 2. See if differences are linked to selection algorithm (e.g. butyrogen list),
# and whether it's appropriate to remove them from algorithm
# 3. Re-run optimized selection algorithm, see how RS rankings change
# Hypothesis: Optimizing RS selection algorithm will change nonresponders' RS, but not change responders

# limitation: some HM's used metaproteomic data to select RS

# save.image("./2025_09_21_lsarp_rapidaim_reanalysis.Renv")
load("./2026_01_07_lsarp_rapidaim_reanalysis.Renv")

# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse"); library("phyloseq")

# load processed R objects
load(file = "./2025_12_02_lsarp_16s_data_meta.Renv")


# :: color vectors ---------------------------------------------------------------

rs.names.pbs = c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")
rs.names = c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors = data.frame(cols = c(gg_color_hue(5)[c(1,1,1,2,3,3,4,4,5)], "#000000"))

# used for nearly all figures
rs.colors = c("Authentic" = "#F8766D", "BobsRedMill" = "#F8766D", "MSPrebiotic" = "#F8766D",
              "LetsDoOrganic" = "#A3A500", "HiMaize260" = "#00BF7D", "Novelose330" = "#00BF7D",
              "ActistarRT" = "#00B0F6", "FibersymRW" = "#00B0F6", "Versafibe1490" = "#E76BF3")
# 


# :: load data ------------------------------------------------------------

# loop through original OTU files used for algorithm

hms = lsarp.metadata.responders$HM %>% unique()

lsarp.rapidaim.dir = "lsarp_rapidaim_reprocess"
CONTROLS <- readRDS("~/Documents/PhD/Resistant_Starch/controls_nonibd-dc-cloud_n83.rds")
TREE <- read_tree( "~/Documents/PhD/Resistant_Starch/97_otus_gg13v5.tree")

# select OTU or ASV
data.type = "otu"#  | asv | otu

# toggle .otu and .asv

lsarp.rapidaim.asv = do.call(merge_phyloseq, lapply(hms, function(hm){
  print(hm)
  #hm = hms[1]
  # select files
  hm.map =  list.files(lsarp.rapidaim.dir)[grepl(hm, x = list.files(lsarp.rapidaim.dir))]
  # load map
  hm.map = read.csv(paste(lsarp.rapidaim.dir, hm.map, sep="/"))
  # fix PCR
  hm.map$pcr <- as.numeric(lapply(hm.map$Post.PCR.DNA.Yield..ng.uL., function(x) replace(x, x == 0, min(hm.map$Post.PCR.DNA.Yield..ng.uL.[hm.map$Post.PCR.DNA.Yield..ng.uL.>0], na.rm = TRUE)/2))) # replace 0's with half the smallest value
  
  # extract seq file
  asv.file = colnames(hm.map)[1]
  asv.file = gsub("\\.", "_", asv.file)
  # select asv/otu file
  if(data.type == "otu"){
    asv.file = list.files(lsarp.rapidaim.dir)[grepl(asv.file, x = list.files(lsarp.rapidaim.dir))]
    asv.file = asv.file[grepl("no_doubletons", x=asv.file)]
    # now load OTU table
    BIOM.RapidAIM <- import_biom(paste(lsarp.rapidaim.dir, asv.file, sep="/"), parseFunction = parse_taxonomy_greengenes) #import biom file
    row.names(hm.map) <- hm.map[,1]
    hm.map <- sample_data(hm.map, errorIfNULL=T) #create sample data phyloseq object from metadata file
    hm.map$Sample <- rownames(hm.map) #add a new column with sample names (for downstream analyses)
    RapidAIM.ps <- merge_phyloseq(BIOM.RapidAIM, hm.map) #Create phyloseq object (with tree)
    RapidAIM.ps #check phyloseq object
    
    # OTU Filtering: Decontam  ---------------------------------------------------------------
    RapidAIM.ps.df <- as.data.frame(sample_data(RapidAIM.ps)) # Put sample_data into a ggplot-friendly data.frame
    RapidAIM.ps.df$LibrarySize <- sample_sums(RapidAIM.ps)
    RapidAIM.ps.df <- RapidAIM.ps.df[order(RapidAIM.ps.df$LibrarySize),]
    RapidAIM.ps.df$Index <- seq(nrow(RapidAIM.ps.df))
    ggplot(data=RapidAIM.ps.df, aes(x=Index, y=LibrarySize, color=Sample_or_Control))+
      geom_point()+geom_hline(yintercept=120000, color="red")+
      theme_bw()+theme(legend.position=c(0.80, 0.15))+
      labs(title="Decontam Depth")
    
    summary(taxa_sums(RapidAIM.ps)) # check if there are taxa absent in all samples
    RapidAIM.ps <- filter_taxa(RapidAIM.ps, function(x) mean(x) != 0, TRUE) #remove taxa that aren't actually present
    
    # Frequency method
    contamdf.freq <- decontam::isContaminant(RapidAIM.ps, method="frequency", conc="pcr")
    RapidAIM.ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, RapidAIM.ps)
    
    # OTU Filtering: Depth --------------------------------------------------------
    RapidAIM.counts <- data.frame(Depth=sample_sums(RapidAIM.ps.noncontam),
                                  Sample=sample_data(RapidAIM.ps.noncontam)$Sample,
                                  RS_Name=sample_data(RapidAIM.ps.noncontam)$RS_Name,
                                  Replicate=sample_data(RapidAIM.ps.noncontam)$Replicate,
                                  stringsAsFactors = FALSE)
    RapidAIM.seq.depth.plot <- ggplot(RapidAIM.counts, aes(x=reorder(RS_Name, Depth), Depth)) + 
      geom_boxplot() +
      geom_point(aes(color=Replicate))+
      geom_hline(yintercept=120000, linetype="solid", color="red",size=0.5)+
      coord_flip()+
      theme_bw()+
      labs(x="", y="Depth (Reads)", title="Sequencing Depth")
    RapidAIM.seq.depth.plot
    RapidAIM.filtered <- prune_samples(sample_sums(RapidAIM.ps.noncontam)>=120000, RapidAIM.ps.noncontam) # Include only samples with "sufficient" depth (120k reads)

    # OTU Filtering: Taxa ---------------------------------------------------------
    
    source("~/Documents/PhD/Resistant_Starch/filter_stats_190509.R")
    RapidAIM.filter.stats <- filter_stats(RapidAIM.filtered,minimum.reads = 0:5)
    RapidAIM.filtered.2 <- filter_taxa(RapidAIM.filtered,function(x) sum(x>=2)>=length(x)*.04,prune = TRUE)
    RapidAIM.filtered.3 <- subset_taxa(RapidAIM.filtered.2, Family!="mitochondria") # remove reads classified as mitochondria
    RapidAIM.filtered.4 <- subset_taxa(RapidAIM.filtered.3, Phylum!="Chloroplast") # remove reads classified as mitochondria
    sample_data(RapidAIM.filtered.4) = sample_data(RapidAIM.filtered.4) %>% data.frame() %>% dplyr::select(RS_Name, HM, pH, Qubit, Sample)
    
    # remove TREE
    
    # now we can merge
    return(RapidAIM.filtered.4)
  }
  # or, do ASVs
  if(data.type == "asv"){
    asv.file = list.files(lsarp.rapidaim.dir)[grepl(asv.file, x = list.files(lsarp.rapidaim.dir))]
    asv.file = asv.file[grepl("pooled", x=asv.file)]
    asv.file = readRDS(paste(lsarp.rapidaim.dir, asv.file, sep="/"))
    
    # fix rownames
    rownames(asv.file) = gsub("_filtered.fastq.gz", "", rownames(asv.file))
    hm.map$Sample = hm.map[,1]
    rownames(hm.map) = gsub("gg13v5.", "", hm.map$Sample)
    rownames(hm.map) = gsub("\\.", "_", rownames(hm.map))
    
    # ASV Chimera Removal -----------------------------------------------------

    
    # chimera removal
    set.seed(25)
    nonchimeras.to.be.kept = dada2::removeBimeraDenovo(
      dada2::getUniques(asv.file), 
      method="consensus", 
      multithread=TRUE, 
      verbose=FALSE)
    
    # make phyloseq
    RapidAIM.ps = merge_phyloseq(otu_table(asv.file, taxa_are_rows=F),
                                 sample_data(hm.map))
    
    # keep only non-bimeras
    RapidAIM.ps <- prune_taxa(colnames(asv.file)[colnames(asv.file) %in% names(nonchimeras.to.be.kept)], RapidAIM.ps)
    
    # ASV Filtering: Decontam  ---------------------------------------------------------------
    RapidAIM.ps.df <- as.data.frame(sample_data(RapidAIM.ps)) # Put sample_data into a ggplot-friendly data.frame
    RapidAIM.ps.df$LibrarySize <- sample_sums(RapidAIM.ps)
    RapidAIM.ps.df <- RapidAIM.ps.df[order(RapidAIM.ps.df$LibrarySize),]
    RapidAIM.ps.df$Index <- seq(nrow(RapidAIM.ps.df))
    ggplot(data=RapidAIM.ps.df, aes(x=Index, y=LibrarySize, color=Sample_or_Control))+
      geom_point()+geom_hline(yintercept=120000, color="red")+
      theme_bw()+theme(legend.position=c(0.80, 0.15))+
      labs(title="Decontam Depth")
    
    summary(taxa_sums(RapidAIM.ps)) # check if there are taxa absent in all samples
    RapidAIM.ps <- filter_taxa(RapidAIM.ps, function(x) mean(x) != 0, TRUE) #remove taxa that aren't actually present
    
    # Frequency method
    contamdf.freq <- decontam::isContaminant(RapidAIM.ps, method="frequency", conc="pcr")
    RapidAIM.ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, RapidAIM.ps)
    
    # ASV Filtering: Depth --------------------------------------------------------
    RapidAIM.counts <- data.frame(Depth=sample_sums(RapidAIM.ps.noncontam),
                                  Sample=sample_data(RapidAIM.ps.noncontam)$Sample,
                                  RS_Name=sample_data(RapidAIM.ps.noncontam)$RS_Name,
                                  Replicate=sample_data(RapidAIM.ps.noncontam)$Replicate,
                                  stringsAsFactors = FALSE)
    RapidAIM.seq.depth.plot <- ggplot(RapidAIM.counts, aes(x=reorder(RS_Name, Depth), Depth)) + 
      geom_boxplot() +
      geom_point(aes(color=Replicate))+
      geom_hline(yintercept=120000, linetype="solid", color="red",size=0.5)+
      coord_flip()+
      theme_bw()+
      labs(x="", y="Depth (Reads)", title="Sequencing Depth")
    RapidAIM.seq.depth.plot
    RapidAIM.filtered <- prune_samples(sample_sums(RapidAIM.ps.noncontam)>=120000, RapidAIM.ps.noncontam) # Include only samples with "sufficient" depth (120k reads)
    
    # ASV Filtering: Taxonomy --------------------------------------------------------
    
    # gg138
    RapidAIM.ps.tax = dada2::assignTaxonomy(colnames(phyloseq::otu_table(RapidAIM.filtered)), "~/Documents/PhD/16S_databases/gg_13_8_train_set_97.fa.gz")
    RapidAIM.ps.tax.df = as.data.frame(RapidAIM.ps.tax)
    RapidAIM.ps.tax.df$ASV = (colnames(phyloseq::otu_table(RapidAIM.filtered))) 
    
    RapidAIM.filtered.tax = merge_phyloseq(RapidAIM.filtered,
                                       tax_table(RapidAIM.ps.tax))
    
    # filter
    source("~/Documents/PhD/Resistant_Starch/filter_stats_190509.R")
    RapidAIM.filter.stats <- filter_stats(RapidAIM.filtered.tax,minimum.reads = 0:5)
    RapidAIM.filtered.2 <- filter_taxa(RapidAIM.filtered.tax,function(x) sum(x>=2)>=length(x)*.04,prune = TRUE)
    RapidAIM.filtered.3 <- subset_taxa(RapidAIM.filtered.2, Family!="mitochondria") # remove reads classified as mitochondria
    RapidAIM.filtered.4 <- subset_taxa(RapidAIM.filtered.3, Phylum!="Chloroplast") # remove reads classified as mitochondria
    sample_data(RapidAIM.filtered.4) = sample_data(RapidAIM.filtered.4) %>% data.frame() %>% dplyr::select(RS_Name, HM, pH, Qubit, Sample)
    
    # now we can merge
    return(RapidAIM.filtered.4)

  }
  }))

 # good
lsarp.rapidaim.otu # ~9000 OTUs x 1127 samples
lsarp.rapidaim.asv # ~7000 ASVs x 1127 samples

lsarp.rapidaim.otu = lsarp.rapidaim.asv

# >> Tax tables -----------------------------------------------------------


lsarp.rapidaim.otu.tax <- data.frame(tax_table(lsarp.rapidaim.otu)) %>%
  mutate(LCA = 
           ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                  ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                         ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                       ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                              ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                     ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))

lsarp.rapidaim.asv.tax <- data.frame(tax_table(lsarp.rapidaim.asv)) %>%
  mutate(LCA = 
           ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                  ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                         ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                       ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                              ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                     ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))


# >> OTU Tables -----------------------------------------------------------

# rarefy
set.seed(25)
lsarp.rapidaim.otu.rare = rarefy_even_depth(lsarp.rapidaim.otu, replace=F, 120000)
(lsarp.rapidaim.otu.rare)
# 9030 taxa remaining

sample_data(lsarp.rapidaim.otu.rare)$HM %>% unique()

# Create full OTU tables df's (median of 3 replicates)
lsarp.rapidaim.otu.mat = speedyseq::psmelt(lsarp.rapidaim.otu.rare) %>%
  # remove extra HMs
  subset(!HM %in% c("808.03", "874.04", "717.03")) %>%
  # fix HMs
  mutate(HM = gsub("\\..*.", "", HM)) %>%
  # add HM0
  mutate(HM = ifelse(nchar(HM) == 3, paste("HM0", HM, sep=""), paste("HM", HM, sep=""))) %>%
  # make unit variable
  mutate(hm_sample = paste(HM, RS_Name, sep="_")) %>%
  # remove unused samples
  subset(RS_Name %in% rs.names.pbs) %>%
  # add LCA
  mutate(LCA = ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                      ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                             ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                    ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                           ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                                  ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                         ifelse(is.na(Phylum), "undefined", paste(Phylum))))))))) %>%
  # calculate sum of LCA reads per sample
  group_by(Sample, LCA) %>%
  mutate(sum.abun = sum(Abundance/120000)) %>% # 
  # make matrix
  dplyr::select(hm_sample, LCA, sum.abun) %>% distinct() %>%
  reshape2::acast(hm_sample ~ LCA, value.var="sum.abun", fun.aggregate=median)

dim(lsarp.rapidaim.otu.mat)
range(lsarp.rapidaim.otu.mat)
rownames(lsarp.rapidaim.otu.mat)
# good


saveRDS(lsarp.rapidaim.otu.mat, "./lsarp_rapidaim_otu_mat.Rds")
saveRDS(lsarp.rapidaim.otu.tax, "./lsarp_rapidaim_otu_tax_mat.Rds")


# :: ANALYSIS -------------------------------------------------------------

lsarp.rapidaim.otu.mat = readRDS("./lsarp_rapidaim_otu_mat.Rds")
lsarp.rapidaim.otu.tax = readRDS("./lsarp_rapidaim_otu_tax_mat.Rds")


# >> Butyrogens -----------------------------------------------------------

# List taxa considered butyrogens
lsarp.rapidaim.otu.butyrogens = subset(lsarp.rapidaim.otu.tax, 
                                       Family=="Lachnospiraceae" | Genus=="Blautia" | Genus=="Roseburia" | Genus=="Eubacterium" | Genus=="Ruminococcus" | Genus=="Clostridium" | Genus=="Faecalibacterium")$LCA %>% unique() # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
lsarp.rapidaim.asv.butyrogens = subset(lsarp.rapidaim.asv.tax, 
                                       Family=="f__Lachnospiraceae" | Genus=="g__Blautia" | Genus=="g__Roseburia" | Genus=="g__Eubacterium" | Genus=="g__Ruminococcus" | Genus=="g__Clostridium" | Genus=="g__Faecalibacterium")$LCA %>% unique() # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
# sum up butyrogens
lsarp.rapidaim.otu.mat = as.data.frame(lsarp.rapidaim.otu.mat)
lsarp.rapidaim.otu.mat$butyrogens = rowSums(lsarp.rapidaim.otu.mat[,colnames(lsarp.rapidaim.otu.mat)%in%lsarp.rapidaim.otu.butyrogens])

lsarp.rapidaim.asv.mat = as.data.frame(lsarp.rapidaim.asv.mat)
lsarp.rapidaim.asv.mat$butyrogens = rowSums(lsarp.rapidaim.asv.mat[,colnames(lsarp.rapidaim.asv.mat)%in%lsarp.rapidaim.asv.butyrogens])

# >> H2S -----------------------------------------------------------

# List taxa considered butyrogens
lsarp.rapidaim.otu.h2s = subset(lsarp.rapidaim.otu.tax, 
                                       Genus=="Veillonella" | Genus=="Atopobium" | Genus=="Fusobacterium" | Genus=="Leptotrichia" | Genus=="Prevotella" | Genus=="Streptococcus")
lsarp.rapidaim.asv.h2s = subset(lsarp.rapidaim.asv.tax, 
                                Genus=="g__Veillonella" | Genus=="g__Atopobium" | Genus=="g__Fusobacterium" | Genus=="g__Leptotrichia" | Genus=="g__Prevotella" | Genus=="g__Streptococcus")
# sum up h2s
lsarp.rapidaim.otu.mat$h2s = rowSums(lsarp.rapidaim.otu.mat[,colnames(lsarp.rapidaim.otu.mat)%in%lsarp.rapidaim.otu.h2s$LCA])

lsarp.rapidaim.asv.mat$h2s = rowSums(lsarp.rapidaim.asv.mat[,colnames(lsarp.rapidaim.asv.mat)%in%lsarp.rapidaim.asv.h2s$LCA])

# >> LFC ------------------------------------------------------------------

pseudo = 0.5/120000

# followed by Log2FC to PBS
# select RS that was selected among RS-treated only

## :: OTU
lsarp.rapidaim.otu.df = reshape2::melt(as.matrix(lsarp.rapidaim.otu.mat)) %>% as.data.frame() %>%
  tidyr::separate(Var1, into=c("HM", "RS_Name"), sep="_", remove=F) %>%
  group_by(HM, Var2) %>%
  mutate(lfc = log2(value+pseudo) - log2(value[RS_Name == "PBS"]+pseudo)) %>%
  subset(RS_Name != "PBS")
# subset to selected RS
lsarp.rapidaim.otu.df = subset(lsarp.rapidaim.otu.df,
                               Var1 %in% unique(paste(lsarp.metadata.responders$HM, lsarp.metadata.responders$RS_Name, sep="_")))
unique(lsarp.rapidaim.otu.df$HM)
# make matrix
lsarp.rapidaim.otu.df.mat = reshape2::acast(lsarp.rapidaim.otu.df, 
                                            HM ~ Var2, value.var="lfc") %>% as.data.frame()
dim(lsarp.rapidaim.otu.df.mat) # 515 taxa
# remove sd = 0
lsarp.rapidaim.otu.df.mat = lsarp.rapidaim.otu.df.mat[,apply(lsarp.rapidaim.otu.df.mat, 2, sd)>0]
dim(lsarp.rapidaim.otu.df.mat) # 297 taxa
# add response
lsarp.rapidaim.otu.df.mat$response = lsarp.metadata.responders[match(rownames(lsarp.rapidaim.otu.df.mat), lsarp.metadata.responders$HM),]$flare.group


# >> Predictors -----------------------------------------------------------

# subset to RS group
lsarp.rapidaim.otu.rs.mat = lsarp.rapidaim.otu.df.mat[rownames(lsarp.rapidaim.otu.df.mat) %in% subset(lsarp.metadata.responders, Group == "RS")$HM,]

# add Phasco - Dialister for OTU
#lsarp.rapidaim.otu.rs.mat$Phasco_Dialister_2 = lsarp.rapidaim.otu.rs.mat$Phascolarctobacterium - lsarp.rapidaim.otu.rs.mat$Dialister
# and But/H2S
#lsarp.rapidaim.otu.rs.mat$bh = lsarp.rapidaim.otu.rs.mat$butyrogens - lsarp.rapidaim.otu.rs.mat$h2s
# fix order
lsarp.rapidaim.otu.rs.mat = lsarp.rapidaim.otu.rs.mat[, c(setdiff(names(lsarp.rapidaim.otu.rs.mat), "response"), "response")]

# identify responses that are different between Responders and Non-Responders
lsarp.rapidaim.otu.rf = ranger::ranger(as.factor(response) ~., data.frame(lsarp.rapidaim.otu.rs.mat), importance="permutation")

# examine variable importances (top 10)
lsarp.rapidaim.otu.rf.imp = data.frame(imp = lsarp.rapidaim.otu.rf$variable.importance) %>% 
  slice_max(order_by=imp, n = 10)

# plots
lsarp.rapidaim.otu.rf.imp.plot = ggplot(lsarp.rapidaim.otu.rf.imp %>% mutate(feature = gsub("\\.", "", gsub("\\X", "", rownames(.)))),
       aes(x=imp, y=reorder(feature, imp)))+
  geom_point(aes(fill=scale(imp)),color="black", shape=21, size=2.5)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  facet_wrap(~"OTU Importance")+
  theme_classic()+theme(legend.position="none")+
  labs(x="Mean Decrease in Accuracy", y="")
lsarp.rapidaim.otu.rf.imp.plot


# OTU: t-test
lsarp.rapidaim.otu.ttest = do.call(rbind, lapply(1:(ncol(lsarp.rapidaim.otu.rs.mat)-1), function(x){
  data.subset = lsarp.rapidaim.otu.rs.mat[,colnames(lsarp.rapidaim.otu.rs.mat)%in%
                                           c(colnames(lsarp.rapidaim.otu.rs.mat)[x], "response")]
  wilcox.output = wilcox.test(subset(data.subset, response == "Relapse")[,1],
                              subset(data.subset, response == "Remit")[,1])
  data.frame(pval = wilcox.output$p.value,
             delta = mean(subset(data.subset, response == "Remit")[,1]) - mean(subset(data.subset, response == "Relapse")[,1]),
             feature = colnames(lsarp.rapidaim.otu.rs.mat)[x])
})) %>% 
  subset(!is.na(pval)) %>%
  arrange(pval)
lsarp.rapidaim.otu.ttest %>% arrange(-pval)
# OTUs: Holdemania, Morganella, Phasco, ~Blautia

# OTU: plot
lsarp.rapidaim.otu.rs.mat.features = reshape2::melt(lsarp.rapidaim.otu.rs.mat) %>%
  subset(variable %in% slice_min(lsarp.rapidaim.otu.ttest, order_by=pval, n=12)$feature)
lsarp.rapidaim.otu.rs.mat.features.plot = ggplot(lsarp.rapidaim.otu.rs.mat.features,
                                                 aes(x=response, y=value))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response), color="white", size=2.5)+
  geom_text(data=slice_min(lsarp.rapidaim.otu.ttest, order_by=pval, n=12) %>%
              mutate(variable = feature),
            aes(label=round(pval, digits=3)), x=1.5, y=Inf,vjust = 1.5, size=3)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~variable, scale="free", ncol=3)+
  labs(x="", y="Log2FC to PBS")
lsarp.rapidaim.otu.rs.mat.features.plot
# Phasco / Dialister ?
ggplot(lsarp.rapidaim.otu.rs.mat,
       aes(x=response, y=Phascolarctobacterium - Dialister))+
  geom_boxplot()+
  geom_point(shape=21, color="white", aes(fill=response), size=2.5)+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Phascolarctobacterium - Dialister\nLog2FC to PBS")
# interesting, the most sig difference is Phasco/Dialister
# Butyrogens:
# Phasco / Dialister ?
ggplot(lsarp.rapidaim.otu.rs.mat,
       aes(x=response, y=butyrogens))+
  geom_boxplot()+
  geom_point(shape=21, color="white", aes(fill=response), size=2.5)+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Butyrogen Log2FC to PBS")


(lsarp.rapidaim.otu.rf.imp.plot+
    lsarp.rapidaim.otu.rs.mat.features.plot+
    patchwork::plot_layout(widths=c(1,3)))+
  ggsave(filename="./lsarp_plots/TEMP_lsarp_phasco_dialister.pdf",
         width=11, height=4,device = cairo_pdf)

# Well that's interesting

# Let's see if the RS had a positive impact on Phasco/Dialister in vivo
# There is no sig association with response + RS group

# No association with Butyrogens, H2S, or But:H2S ratio

# :: Final Analysis ----------------------------------------------------------

lsarp.rapidaim.otu.df.mat = readRDS("./lsarp_rapidaim_otu_mat.Rds")


# :: :: OTU ---------------------------------------------------------------

# subset to RS group
lsarp.rapidaim.otu.rs.mat = lsarp.rapidaim.otu.df.mat[rownames(lsarp.rapidaim.otu.df.mat) %in% subset(lsarp.metadata.responders, Group == "RS")$HM,]
dim(lsarp.rapidaim.otu.rs.mat)

# add Phasco - Dialister for OTU
#lsarp.rapidaim.otu.rs.mat$Phasco_Dialister_2 = lsarp.rapidaim.otu.rs.mat$Phascolarctobacterium - lsarp.rapidaim.otu.rs.mat$Dialister
# fix order
lsarp.rapidaim.otu.rs.mat = lsarp.rapidaim.otu.rs.mat[, c(setdiff(names(lsarp.rapidaim.otu.rs.mat), "response"), "response")]

# identify responses that are different between Responders and Non-Responders

# remove 0 variance predictors
lsarp.rapidaim.otu.rs.mat.responses = lsarp.rapidaim.otu.rs.mat$response

lsarp.rapidaim.otu.rs.mat = 
  lsarp.rapidaim.otu.rs.mat[,c(1:(ncol(lsarp.rapidaim.otu.rs.mat)-1))][,apply(lsarp.rapidaim.otu.rs.mat[,c(1:(ncol(lsarp.rapidaim.otu.rs.mat)-1))],
        2, sd)>0]
lsarp.rapidaim.otu.rs.mat$response = lsarp.rapidaim.otu.rs.mat.responses

dim(lsarp.rapidaim.otu.rs.mat)
# 14 samples, 260 predictors

# :: :: PCA --------------------------------------------------------------


lsarp.rapidaim.otu.rs.pca = prcomp((lsarp.rapidaim.otu.rs.mat%>%dplyr::select(-response)))
# extract data
lsarp.rapidaim.otu.rs.pca.df = lsarp.rapidaim.otu.rs.pca$x[,c(1:4)] %>% data.frame()
# extract var explained
lsarp.rapidaim.otu.rs.pca.var <- ((lsarp.rapidaim.otu.rs.pca$sdev^2) / sum(lsarp.rapidaim.otu.rs.pca$sdev^2))[1:2]
# merge with meta
lsarp.rapidaim.otu.rs.pca.df$code = rownames(lsarp.rapidaim.otu.rs.pca.df)
lsarp.rapidaim.otu.rs.pca.df$response =  lsarp.rapidaim.otu.rs.mat$response

# plot PCA
lsarp.rapidaim.otu.rs.pca.plot = ggplot(lsarp.rapidaim.otu.rs.pca.df,
                                   aes(x=PC1, y=PC2))+
  #geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(aes(shape = response, fill=response), color="white", size=2.5)+
  stat_ellipse(aes(color=response), alpha = 0.5, level = 0.95)+
  #ggrepel::geom_text_repel(aes(label = code),size=2)+
  scale_shape_manual(values=c(21,22,23), guide=F)+
  theme_classic()+theme(#legend.position = "none",
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))+
  labs(x=paste("PC1:", round(lsarp.rapidaim.otu.rs.pca.var, digits=4)[1]*100, "%"),
       y=paste("PC2:", round(lsarp.rapidaim.otu.rs.pca.var, digits=4)[2]*100, "%"))
lsarp.rapidaim.otu.rs.pca.plot

# PCA doesn't discriminate


# :: :: LOOCV ---------------------------------------------------------------

# Note: OOB malfunctions; systematically predicts the wrong class

# repeat 15 times and calculate AUC; LOOCV
lsarp.rapidaim.otu.rf.caret.15.loocv = lapply(1:15, function(iter){
  print(iter)
  set.seed(iter)
  lapply(1:nrow(lsarp.rapidaim.otu.rs.mat), function(hm){
    train.data = lsarp.rapidaim.otu.rs.mat[-hm,]
    test.data = data.frame((lsarp.rapidaim.otu.rs.mat[hm,]))
    
    # balance dataset (to n=4-5 each)
    # smaller value
    sample.size = min(table(train.data$response))
    
    set.seed(iter)
    train.data = train.data %>%
      group_by(response) %>%
      slice_sample(n=sample.size, replace=F) %>%
      mutate(response = factor(response, levels=c("Relapse", "Remit")))
    
  lsarp.rapidaim.otu.rf = randomForestSRC::rfsrc((response) ~., data.frame(train.data), 
                                         #importance="permutation",
                                         #num.trees = 500,
                                         #min.node.size = 1,
                                         probability=T)
  
  rf.predictions = data.frame(prob.relapse = predict(lsarp.rapidaim.otu.rf, test.data %>% dplyr::select(-response))$predicted[1],
             true.relapse = test.data$response,
             hm = hm,
             iter = iter)
  
  # now extract feature importances
  lsarp.rapidaim.otu.rf.imp = randomForestSRC::vimp(lsarp.rapidaim.otu.rf, method = "permute", joint=F) 
  
  lsarp.rapidaim.otu.rf.imp = data.frame(imp = lsarp.rapidaim.otu.rf.imp$importance[,1]) %>%
    mutate(iter = iter) %>%
    mutate(feature = rownames(lsarp.rapidaim.otu.rf.imp$importance)) %>%
    arrange(imp)
  return(list(rf.predictions,
              lsarp.rapidaim.otu.rf.imp))
  })
})

lsarp.rapidaim.otu.rf.caret.15.loocv.preds = 
  do.call(rbind, lapply(1:14, function(hm){
  do.call(rbind, lapply(1:15, function(iter){
  purrr::pluck(purrr::pluck(lsarp.rapidaim.otu.rf.caret.15.loocv[[iter]][[hm]])[[1]])
}))}))

lsarp.rapidaim.otu.rf.caret.15.loocv.auc = lsarp.rapidaim.otu.rf.caret.15.loocv.preds %>%
  group_by(iter) %>%
  mutate(auc = pROC::auc(true.relapse, prob.relapse, levels=c("Relapse", "Remit"))[1]) %>%
  dplyr::select(iter, auc) %>% distinct() %>%
  ungroup() %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(-iter) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high) %>% distinct()
lsarp.rapidaim.otu.rf.caret.15.loocv.auc

# 0.704

# next: do scrambled
# then: do feature importances


# :: :: LOOCV Scrambled ---------------------------------------------------------------

# repeat 15 times and calculate AUC; LOOCV
lsarp.rapidaim.otu.rf.caret.15.scramble.loocv = lapply(1:15, function(iter){
  print(iter)
  set.seed(iter)
  
  # randomly shuffle class labels
  lsarp.rapidaim.otu.rs.mat = lsarp.rapidaim.otu.rs.mat %>%
    mutate(response = as.factor(sample(response, nrow(lsarp.rapidaim.otu.rs.mat), replace=F))) 
    
    
  lapply(1:nrow(lsarp.rapidaim.otu.rs.mat), function(hm){
    train.data = lsarp.rapidaim.otu.rs.mat[-hm,]
    test.data = data.frame((lsarp.rapidaim.otu.rs.mat[hm,]))
    
    # balance dataset (5:5)
    # smaller value
    sample.size = min(table(train.data$response))
    
    set.seed(iter)
    train.data = train.data %>%
      group_by(response) %>%
      slice_sample(n=sample.size, replace=F) %>%
      data.frame() %>%
      mutate(response = factor(response, levels=c("Relapse", "Remit")))
    
    lsarp.rapidaim.otu.rf = randomForestSRC::rfsrc((response) ~., data.frame(train.data), 
                                                   #importance="permutation",
                                                   #num.trees = 500,
                                                   #min.node.size = 1,
                                                   probability=T)
    
    rf.predictions = data.frame(prob.relapse = predict(lsarp.rapidaim.otu.rf, test.data %>% dplyr::select(-response))$predicted[1],
                                true.relapse = test.data$response,
                                hm = hm,
                                iter = iter)
    
    # now extract feature importances
    lsarp.rapidaim.otu.rf.imp = randomForestSRC::vimp(lsarp.rapidaim.otu.rf, method = "permute", joint=F) 
    
    lsarp.rapidaim.otu.rf.imp = data.frame(imp = lsarp.rapidaim.otu.rf.imp$importance[,1]) %>%
      mutate(iter = iter) %>%
      mutate(feature = rownames(lsarp.rapidaim.otu.rf.imp$importance)) %>%
      arrange(imp)
    return(list(rf.predictions,
                lsarp.rapidaim.otu.rf.imp))
  })
})

lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.preds = 
  do.call(rbind, lapply(1:14, function(hm){
    do.call(rbind, lapply(1:15, function(iter){
      purrr::pluck(purrr::pluck(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv[[iter]][[hm]])[[1]])
    }))}))

lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.auc = lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.preds %>%
  group_by(iter) %>%
  mutate(auc = pROC::auc(true.relapse, prob.relapse, 
                         levels=c("Relapse", "Remit"))[1]) %>%
  dplyr::select(iter, auc) %>% distinct() %>%
  ungroup() %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(-iter) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high) %>% distinct()
lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.auc
# Scramble AUC = 0.599

# check importances


# :: :: ROC ---------------------------------------------------------------

# plot
lsarp.rapidaim.otu.rf.15.loocv.roc = do.call(rbind, lapply(1:15, function(iter){
  data.subset = subset(lsarp.rapidaim.otu.rf.caret.15.loocv.preds, iter == iter)
  # calibrate pred; not necessary
  #data.subset$pred = 1 / (1 + exp(-data.subset$pred))
  data.frame(sens = pROC::roc(response = data.subset$true.relapse, predictor = data.subset$prob.relapse,  # define case = "Relapse"
                              levels=c("Relapse", "Remit"))$sensitivities,
             spec = pROC::roc(response = data.subset$true.relapse, predictor = data.subset$prob.relapse,# define case = "Relapse"
                              levels=c("Relapse", "Remit"))$specificities,
             iter = iter)
}))

# Step 2: Compute mean and standard error for sensitivities across iterations
lsarp.rapidaim.otu.rf.15.loocv.roc_summary <- 
  lsarp.rapidaim.otu.rf.15.loocv.roc %>%
  mutate(sens = round(sens, digits=1),
         spec = round(spec, digits=1))%>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


lsarp.rapidaim.otu.rf.15.loocv.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_smooth(data = lsarp.rf.selected.output.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter),
  #          se=T, linewidth=0.3, color="grey")+
  # Error ribbon (95% CI)
  geom_ribbon(data = lsarp.rapidaim.otu.rf.15.loocv.roc_summary, 
              aes(x = fpr, ymin = lower, ymax = upper),
              alpha=0.2)+
  # Mean ROC curve
  geom_path(data = lsarp.rapidaim.otu.rf.15.loocv.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black", linewidth=1.5)+
  geom_path(data = lsarp.rapidaim.otu.rf.15.loocv.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="white", linewidth=1)+
  
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("LOOCV\nAUC: ", round(lsarp.rapidaim.otu.rf.caret.15.loocv.auc$mean.auc, digits=2),
                       "\n(", round(lsarp.rapidaim.otu.rf.caret.15.loocv.auc$auc.low, digits=2), 
                       ", ", round(lsarp.rapidaim.otu.rf.caret.15.loocv.auc$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position="none")+
  facet_wrap(~"OTU RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
lsarp.rapidaim.otu.rf.15.loocv.roc.plot


# :: :: ROC Scrambled ---------------------------------------------------------------

# plot
lsarp.rapidaim.otu.rf.15.scramble.loocv.roc = do.call(rbind, lapply(1:15, function(iter){
  data.subset = subset(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.preds, iter == iter)
  # calibrate pred; not necessary
  #data.subset$pred = 1 / (1 + exp(-data.subset$pred))
  data.frame(sens = pROC::roc(response = data.subset$true.relapse, predictor = data.subset$prob.relapse,  # define case = "Relapse"
                              levels=c("Relapse", "Remit"))$sensitivities,
             spec = pROC::roc(response = data.subset$true.relapse, predictor = data.subset$prob.relapse,# define case = "Relapse"
                              levels=c("Relapse", "Remit"))$specificities,
             iter = iter)
}))

# Step 2: Compute mean and standard error for sensitivities across iterations
lsarp.rapidaim.otu.rf.15.scramble.loocv.roc_summary <- lsarp.rapidaim.otu.rf.15.scramble.loocv.roc %>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


lsarp.rapidaim.otu.rf.15.scramble.loocv.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_smooth(data = lsarp.rf.selected.output.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter),
  #          se=T, linewidth=0.3, color="grey")+
  # Error ribbon (95% CI)
  geom_ribbon(data = lsarp.rapidaim.otu.rf.15.scramble.loocv.roc_summary, 
              aes(x = fpr, ymin = lower, ymax = upper),
              alpha=0.2)+
  # Mean ROC curve
  geom_path(data = lsarp.rapidaim.otu.rf.15.scramble.loocv.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black", linewidth=1.5)+
  geom_path(data = lsarp.rapidaim.otu.rf.15.scramble.loocv.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="white", linewidth=1)+
  
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("LOOCV\nAUC: ", round(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.auc$mean.auc, digits=2),
                       "\n(", round(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.auc$auc.low, digits=2), 
                       ", ", round(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.auc$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position="none")+
  facet_wrap(~"Scrambled OTU RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
lsarp.rapidaim.otu.rf.15.scramble.loocv.roc.plot


# :: :: Importances ---------------------------------------------------------------

lsarp.rapidaim.otu.rf.caret.15.loocv.importances = 
  do.call(rbind, lapply(1:14, function(hm){
    do.call(rbind, lapply(1:15, function(iter){
      purrr::pluck(purrr::pluck(lsarp.rapidaim.otu.rf.caret.15.loocv[[iter]][[hm]])[[2]])
    }))}))%>%
  group_by(feature) %>%
  summarize(mean.imp = mean(imp),
            imp.low = mean(imp) - (sd(imp)/sqrt(n()) * 1.96), # 95% CI
            imp.high = mean(imp) + (sd(imp)/sqrt(n()) * 1.96)) %>% 
  as.data.frame() %>%
  arrange(-mean.imp)

# take top 15 mean importance
lsarp.rapidaim.otu.rf.caret.15.loocv.importances.top = lsarp.rapidaim.otu.rf.caret.15.loocv.importances %>%
  slice_max(order_by = mean.imp, n = 15)

# plots
lsarp.rapidaim.otu.rf.imp.plot = ggplot(lsarp.rapidaim.otu.rf.caret.15.loocv.importances.top %>% mutate(feature = gsub("\\.", "", gsub("\\X", "", feature))),
                                        aes(x=mean.imp, y=reorder(feature, mean.imp)))+
  geom_segment(aes(x = imp.low, xend = imp.high, y=feature, yend=feature))+
  geom_point(aes(fill=scale(mean.imp)), shape=21, size=2.5)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  facet_wrap(~"OTU Importance")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Mean Decrease in Accuracy", y="")
lsarp.rapidaim.otu.rf.imp.plot


# :: :: Importances Scrambled ---------------------------------------------------------------

lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.importances = 
  do.call(rbind, lapply(1:14, function(hm){
    do.call(rbind, lapply(1:15, function(iter){
      purrr::pluck(purrr::pluck(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv[[iter]][[hm]])[[2]])
    }))}))%>%
  group_by(feature) %>%
  summarize(mean.imp = mean(imp),
            imp.low = mean(imp) - (sd(imp)/sqrt(n()) * 1.96), # 95% CI
            imp.high = mean(imp) + (sd(imp)/sqrt(n()) * 1.96)) %>% 
  as.data.frame() %>%
  arrange(-mean.imp)

# take top 15 mean importance
lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.importances.top = lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.importances %>%
  slice_max(order_by = mean.imp, n = 15)

# plots
lsarp.rapidaim.otu.rf.scrambled.imp.plot = ggplot(lsarp.rapidaim.otu.rf.caret.15.scramble.loocv.importances.top %>% mutate(feature = gsub("\\.", "", gsub("\\X", "", feature))),
                                        aes(x=mean.imp, y=reorder(feature, mean.imp)))+
  geom_segment(aes(x = imp.low, xend = imp.high, y=feature, yend=feature))+
  geom_point(aes(fill=scale(mean.imp)), shape=21, size=2.5)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  facet_wrap(~"OTU Importance")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Mean Decrease in Accuracy", y="")
lsarp.rapidaim.otu.rf.scrambled.imp.plot


# :: :: Phasco/Dialister --------------------------------------------------

lsarp.rapidaim.otu.rs.pd = data.frame(
  response = lsarp.rapidaim.otu.rs.mat$response,
  pd = lsarp.rapidaim.otu.rs.mat$Phascolarctobacterium - lsarp.rapidaim.otu.rs.mat$Dialister)
# wilcox
lsarp.rapidaim.otu.rs.pd.wilcox = wilcox.test(subset(lsarp.rapidaim.otu.rs.pd, response == "Relapse")$pd,
                                              subset(lsarp.rapidaim.otu.rs.pd, response == "Remit")$pd)$p.value

lsarp.rapidaim.otu.rs.pd.plot = ggplot(lsarp.rapidaim.otu.rs.pd,
                                                 aes(x=response, y=pd))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response), color="white", size=3)+
  annotate(geom="text", label=round(lsarp.rapidaim.otu.rs.pd.wilcox, digits=3), x=1.5, y=Inf,vjust = 1.5, size=3.5)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Phascolarctobacterium / Dialister Log2FC")+
  labs(x="", y="Log2FC to PBS")
lsarp.rapidaim.otu.rs.pd.plot

# AUC

# repeat 0 times and calculate AUC; LOOCV (don't need to repeat because LR is deterministic)
lsarp.rapidaim.otu.lr.15.loocv = do.call(rbind, lapply(1, function(iter){
  print(iter)
  set.seed(iter)
  do.call(rbind, lapply(1:nrow(lsarp.rapidaim.otu.rs.pd), function(hm){
    pd.data = lsarp.rapidaim.otu.rs.pd
    train.data = pd.data[-hm,]
    test.data = data.frame((pd.data[hm,]))
    
    lr.model = glm(response ~ (scale(pd)),
                         family=binomial(link='logit'),
                         data=train.data%>%mutate(response = as.factor(ifelse(response == "Relapse", 0, 1))))
    data.frame(pred.value = predict(lr.model, test.data, probability=T),
               true.value = test.data$response,
               hm = hm,
               iter = iter)
  }))
}))
lsarp.rapidaim.otu.lr.15.loocv
pROC::auc(lsarp.rapidaim.otu.lr.15.loocv$true.value,
          (lsarp.rapidaim.otu.lr.15.loocv$pred.value))
# AUC 0.711 using P/D alone
ggplot(lsarp.rapidaim.otu.lr.15.loocv %>%
         # convert pred.value to probability
         mutate(pred.value = sigmoid::sigmoid(pred.value)),
       aes(x=true.value, y=pred.value))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, color="white", fill="black", size=3)+
  theme_classic()+
  labs(x="True label", y="Probability")

# :: :: Add Phasco/Dialister ----------------------------------------------


# repeat 15 times and calculate AUC; LOOCV
lsarp.rapidaim.otu.rf.caret.pd.15.loocv = lapply(1:15, function(iter){
  print(iter)
  set.seed(iter)
  # add PD
  lsarp.rapidaim.otu.rs.mat$pd = lsarp.rapidaim.otu.rs.mat$Phascolarctobacterium - lsarp.rapidaim.otu.rs.mat$Dialister
  
  lapply(1:nrow(lsarp.rapidaim.otu.rs.mat), function(hm){
    train.data = lsarp.rapidaim.otu.rs.mat[-hm,]
    test.data = data.frame((lsarp.rapidaim.otu.rs.mat[hm,]))
    
    # balance dataset (5:5)
    # smaller value
    sample.size = min(table(train.data$response))
    
    set.seed(iter)
    train.data = train.data %>%
      group_by(response) %>%
      slice_sample(n=sample.size, replace=F) %>%
      mutate(response = factor(response, levels=c("Relapse", "Remit")))
    
    lsarp.rapidaim.otu.rf = randomForestSRC::rfsrc((response) ~., data.frame(train.data), 
                                                   #importance="permutation",
                                                   #num.trees = 500,
                                                   #min.node.size = 1,
                                                   probability=T)
    
    rf.predictions = data.frame(prob.relapse = predict(lsarp.rapidaim.otu.rf, test.data %>% dplyr::select(-response))$predicted[1],
                                true.relapse = test.data$response,
                                hm = hm,
                                iter = iter)
    
    # now extract feature importances
    lsarp.rapidaim.otu.rf.imp = randomForestSRC::vimp(lsarp.rapidaim.otu.rf, method = "permute", joint=F) 
    
    lsarp.rapidaim.otu.rf.imp = data.frame(imp = lsarp.rapidaim.otu.rf.imp$importance[,1]) %>%
      mutate(iter = iter) %>%
      mutate(feature = rownames(lsarp.rapidaim.otu.rf.imp$importance)) %>%
      arrange(imp)
    return(list(rf.predictions,
                lsarp.rapidaim.otu.rf.imp))
  })
})

lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances = 
  do.call(rbind, lapply(1:14, function(hm){
    do.call(rbind, lapply(1:15, function(iter){
      purrr::pluck(purrr::pluck(lsarp.rapidaim.otu.rf.caret.pd.15.loocv[[iter]][[hm]])[[2]])
    }))}))%>%
  group_by(feature) %>%
  summarize(mean.imp = mean(imp),
            imp.low = mean(imp) - (sd(imp)/sqrt(n()) * 1.96), # 95% CI
            imp.high = mean(imp) + (sd(imp)/sqrt(n()) * 1.96)) %>% 
  as.data.frame() %>%
  arrange(-mean.imp)

# take top 15 mean importance
lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances.top = lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances %>%
  slice_max(order_by = mean.imp, n = 15)

# plots
lsarp.rapidaim.otu.rf.pd.imp.plot = ggplot(lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances.top %>% 
                                             mutate(feature = gsub("\\.", "", gsub("\\X", "", feature))) %>%
                                             mutate(feature = ifelse(feature == "pd", "Phascolarctobacterium/Dialister", feature)),
                                        aes(x=mean.imp, y=reorder(feature, mean.imp)))+
  geom_segment(aes(x = imp.low, xend = imp.high, y=feature, yend=feature))+
  geom_point(aes(fill=scale(mean.imp)), shape=21, size=2.5)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  facet_wrap(~"OTU Importance")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Mean Decrease in Accuracy", y="")
lsarp.rapidaim.otu.rf.pd.imp.plot

# PD = 0.04
# Phasco = 0.01437
subset(lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances.top, feature == "pd")$mean.imp /
  subset(lsarp.rapidaim.otu.rf.caret.pd.15.loocv.importances.top, feature == "Phascolarctobacterium")$mean.imp
# P/D is 2.748x more important than Phasco

# :: Plots ----------------------------------------------------------------

(lsarp.rapidaim.otu.rf.15.loocv.roc.plot|
lsarp.rapidaim.otu.rf.imp.plot|
lsarp.rapidaim.otu.rf.pd.imp.plot|
lsarp.rapidaim.otu.rs.pd.plot) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_phasco_plots.pdf",
         width=18, height=4,device = cairo_pdf)
  

