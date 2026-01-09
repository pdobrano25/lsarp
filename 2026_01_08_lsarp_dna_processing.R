### 2026_01_06  LSARP-CD DNA Processing (including metadata)

# Prepares FFQ data (with minimal analyses), 16S data variables (e.g. butyrogens, predicted microbial load, functional redundancy), metagenomic data, metaproteomic data (e.g. proteins collapsed to functions), and metabolomic data. Output serves as input for most analysis scripts.



# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse"); library("patchwork")


# :: necessary files ------------------------------------------------------

# LSARP-CD patient list
"2025_05_29_lsarp_cd_rapidaim_list_2.csv"

# Master mapping file (for MBX; revised June 2025)
"2025_06_15_mbx_mapping.Rds"

# James' ASV table
"all_trials_16S_min50k_mappingFile_241030_gg2species_physeq_pooled_noTree_241115.Rds"
# updated:
"2025_09_08_full_mapping_rapidaim.jb_gg2species_physeq_pooled_wTree_250909.Rds"

# James' MGX manifest
"humann3_outputs_summary_250502.csv"

# James' MGX files (in several directories)
"metaphlan_bugs_list"

# RapidAIM data
"2025_03_11_full_mapping_rapidaim.csv"


# set vectors for RS_Name, with and without PBS
rs.names.pbs = c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")
rs.names = c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")


# >> 16S ------------------------------------------------------------------


# :: process ALL 16S data --------------------------------------------------------

# NOTE: this was performed once for OARS; we can use the same processed data (since we know the tax assignments will be held constant)


# :: prepare LSARP mapping ----------------------------------------------------

# obtain list of LSARP-CD patients for analysis (Dave's list)
lsarp.patient.list <- readxl::read_excel("./2024_07_31_unblinding_DM.xlsx") %>% as.data.frame()
# remove non-compliant
lsarp.patient.list = subset(lsarp.patient.list, is.na(non.compliant)) # based on compliance %
lsarp.patient.list = subset(lsarp.patient.list, !hm %in% c(966, 970, 974)) # based on Dave's comment "Flare before starting RS"
# make other HM column (short form, without .00)
lsarp.patient.list$HM = substr(lsarp.patient.list$study_id, 1, 6)
# remove unnecessary rows
lsarp.patient.list = subset(lsarp.patient.list, !is.na(HM))
# remove unnecessary columns
lsarp.patient.list = lsarp.patient.list[,c("HM", "study_id","scope.day", "flare.call", "group", "flare.day", "days.taken")]

# check that stools are LSARP and not other study; load main mapping file
metadata.lsarp.stool <- readRDS("./2025_06_15_mbx_mapping.Rds") %>% subset(!HM %in% c("Pos", "Neg") & grepl("LSARP", stl_study))
metadata.lsarp.stool[,c("standard.name")] %>% unique() %>% length() # 355 unique stools
# subset to those intended for analysis (from Dave's list)
metadata.lsarp.stool = subset(metadata.lsarp.stool, HM %in% lsarp.patient.list$HM)
metadata.lsarp.stool$study_id = paste(metadata.lsarp.stool$HM, ".00", sep="")
length(unique(metadata.lsarp.stool$standard.name)) # 196 unique stools

# remove columns that mess things up due to MBX being different (i.e. replicates)
metadata.lsarp.stool = metadata.lsarp.stool[,c("HM", "standard.name", "study_id",
                                               "stool_date_rec_v2",
                                               "rs_start_date", "product_end",
                                               "lsarp.timing","lsarp.days","lsarp.off", "lsarp.on.rs",
                                               "lsarp.rs", "diag", "gender", "ph_response", "fecalcal_res")] %>% distinct()

# add "phase"; phase 1 = on RS + baseline sample; phase 2 = off RS (washout)
# note: for statistics, baseline sample is added to phase 2 as comparison group
metadata.lsarp.stool$phase = ifelse(metadata.lsarp.stool$lsarp.off <= 0, "treatment", "washout")
# note: this way is more accurate than going by month 0-5, 6-12
# slightly more complicated than for OARS analysis, because it's not discretely 0,3,6 = treatment

metadata.lsarp.stool$baseline = ifelse(grepl("baseline", metadata.lsarp.stool$lsarp.timing), "baseline", "not_baseline")

# good
metadata.lsarp.stool = metadata.lsarp.stool %>% distinct() %>% data.frame()
dim(metadata.lsarp.stool)
# 196 unique stools

## add other notes
metadata.lsarp.stool = merge(metadata.lsarp.stool,
                             lsarp.patient.list[,c("HM",
                                                 "flare.call", "group", "flare.day", "days.taken", "scope.day")], by="HM")
# fix names
metadata.lsarp.stool$group = ifelse(metadata.lsarp.stool$group == "Plac", "Placebo", "RS")

metadata.lsarp.stool$flare.day

# flesh out flare.day and days.taken [i.e. there is 1 value each per HM, so replace NA with that value]
metadata.lsarp.stool = metadata.lsarp.stool %>%
  group_by(HM) %>%
  mutate(flare.day = ifelse(all(is.na(flare.day)), NA, flare.day[!is.na(flare.day)]),
         days.taken = days.taken[!is.na(days.taken)]) %>% as.data.frame()

# Goal: load 16S data (rarefied) and replace metadata
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds")
# subset to LSARP-CD samples
amplicon.data.gg.rare.lsarp = phyloseq::subset_samples(amplicon.data.gg.rare, standard.name %in% metadata.lsarp.stool$standard.name)
# merge old sample_data with new metadata
amplicon.data.gg.rare.lsarp.new.meta = merge(data.frame(phyloseq::sample_data(amplicon.data.gg.rare.lsarp)),
                                            metadata.lsarp.stool, by="standard.name")
amplicon.data.gg.rare.lsarp.new.meta$SampleID = gsub("\\.", "_", gsub("gg13v5.", "", amplicon.data.gg.rare.lsarp.new.meta$SampleID))
rownames(amplicon.data.gg.rare.lsarp.new.meta) = amplicon.data.gg.rare.lsarp.new.meta$SampleID
# merge in new metadata
amplicon.data.gg.rare.lsarp = phyloseq::merge_phyloseq(phyloseq::otu_table(amplicon.data.gg.rare.lsarp),
                                                      phyloseq::tax_table(amplicon.data.gg.rare.lsarp),
                                                      phyloseq::sample_data(amplicon.data.gg.rare.lsarp.new.meta))
# good
lsarp.asv.meta = phyloseq::sample_data(amplicon.data.gg.rare.lsarp) %>% data.frame()

lsarp.asv.data = speedyseq::psmelt(amplicon.data.gg.rare.lsarp)

# calculate median value of taxa across technical replicates (e.g. 3 barcodes per standard.name (stool))
lsarp.asv.data.median = reshape2::acast(lsarp.asv.data,
                                       standard.name ~ Taxa, value.var="Abundance", fun.aggregate=median)
dim(lsarp.asv.data.median) # 192 x 10843 taxa

# delete taxa not present
lsarp.asv.data.median = lsarp.asv.data.median[,colSums(lsarp.asv.data.median)>0]
dim(lsarp.asv.data.median) # 192 x 4902 taxa left

# use this (^) for Alpha-diversity

# prepare Glom, too (sum up LCA)
lsarp.asv.data.glom = speedyseq::psmelt(amplicon.data.gg.rare.lsarp) %>%
  group_by(Sample, LCA) %>%
  mutate(Abundance = sum(Abundance)) %>%
  dplyr::select(Sample, standard.name, LCA, Abundance) %>% distinct()
# then cast while taking median per standard.name (if median was 0, but reads in others, this will collapse taxa to not present)
lsarp.asv.data.glom = reshape2::acast(lsarp.asv.data.glom,
                                     standard.name ~ LCA, value.var="Abundance", 
                                     fun.aggregate=median)
dim(lsarp.asv.data.glom)
range(lsarp.asv.data.glom)/50000*100
# good; 192 samples x 1554 glommed taxa; range is below 100%

# examine missing stools
metadata.lsarp.stool$missing.16s = ifelse(!metadata.lsarp.stool$standard.name %in% data.frame(phyloseq::sample_data(amplicon.data.gg.rare.lsarp))$standard.name, "missing", "")
sum(metadata.lsarp.stool$missing.16s == "missing")
# missing 4 stools; proceed

subset(metadata.lsarp.stool, missing.16s == "missing")
# these have been re-attempted several times now (DNA extraction + sequencing)

# fix edge cases (timing was mislabelled)
metadata.lsarp.stool$lsarp.timing = ifelse(metadata.lsarp.stool$standard.name == "HM0999-STL-11", "Month 12 post start of product", 
                                           metadata.lsarp.stool$lsarp.timing)
metadata.lsarp.stool$lsarp.timing = ifelse(metadata.lsarp.stool$standard.name == "HM0933-STL-04", "Month 1 post start of product", 
                                           metadata.lsarp.stool$lsarp.timing)
metadata.lsarp.stool$lsarp.timing = ifelse(metadata.lsarp.stool$standard.name == "HM1035-STL-05", "Month 3 post start of product", 
                                           metadata.lsarp.stool$lsarp.timing)


# >> Calculate Fibre + RS ------------------------------------------------------

# 2025_07_08  Goal: add supplemented RS to RS intake

# load REDCap for RS dose
lsarp.redcap.data.all.loaded <- read.csv("~/Downloads/redCap_sampleData_uniques_240724.tsv", sep="\t")
# unpack monash date
lsarp.redcap.data.loaded = lsarp.redcap.data.all.loaded %>%
  # subset to HM's in the study
  subset(study_id %in% metadata.lsarp.stool$study_id) %>%
  # clean up ffq.date name
  mutate(ffq.date = as.Date(substr(monash.ffq.date, 
                                   nchar(monash.ffq.date)-6,
                                   nchar(monash.ffq.date)), "%d%b%y")) %>% 
  # remove NA entries
  subset(!is.na(ffq.date)) %>%
  # since these are new onset, we can simply subset to before withdrawn from study; diagnosis and any others will be ~within study period
  subset(ffq.date <= rs.trial.withdraw.date)

# note: we can use FFQ as predictors for ML models:
redcap.data.ffq = lsarp.redcap.data.loaded %>% select(apple:other_take_away_meals)
# change to ordinal
redcap.data.ffq.map = data.frame(entry = unique(unlist(redcap.data.ffq)),
                                 value = rank(c(2/7, 1/7, .5/30, 1/30, 0, 5/7, 1/1, 2/1, 4/1, 6/1))-1) %>%
  arrange(value)
# apply replacement
redcap.data.ffq = do.call(cbind, lapply(1:ncol(redcap.data.ffq), function(x){
  data.subset = data.frame(feature = redcap.data.ffq[,x])
  # apply replacement
  data.subset$feature = redcap.data.ffq.map$value[match(data.subset$feature,redcap.data.ffq.map$entry)]
  colnames(data.subset) = colnames(redcap.data.ffq)[x]
  data.subset
}))
# add standard.name
redcap.data.ffq$standard.name = lsarp.redcap.data.loaded$standard.name
redcap.data.ffq$HM = substr(redcap.data.ffq$standard.name, 1, 6)
# keep study stools
redcap.data.ffq = subset(redcap.data.ffq, standard.name %in% subset(metadata.lsarp.stool, baseline=="baseline")$standard.name)
# only 12 have baseline FFQ

# now for calculated nutrients, fibre, RS, etc
# clean up to extract out necessary diet data (including BSA (m2))
lsarp.redcap.data.loaded = subset(lsarp.redcap.data.loaded, !is.na(dose))[,c("study_id", "fecalcal_res","ffq.date",
                                                                             "energy", "resistant_starch","dietary_fibre", "dose", "dose_2", "dose_3",
                                                                             "percent_product", "percent_product_2", "percent_product_3",
                                                                             "lot_number", "lot_number_2", "lot_number_3", "m2")] %>% distinct()
unique(lsarp.redcap.data.loaded$study_id) # all 27
metadata.lsarp.stool$study_id %>% unique() # all 27


## calculate adj.fibre and adj.rs (adjusting for caloric intake)
lsarp.redcap.data.loaded$fibre_adj = lsarp.redcap.data.loaded$dietary_fibre / lsarp.redcap.data.loaded$energy * 4.184 * 1000
lsarp.redcap.data.loaded$rs_adj = lsarp.redcap.data.loaded$resistant_starch / lsarp.redcap.data.loaded$energy * 4.184 * 1000

# note: we need to consider the supplemented RS intake; so add to dietary RS intake
# but also consider that RS lots vary by RS% across lots, and body surface area changes
# simple solution: take average of doses/lots, since dates were not recorded with high precision

# take average of RS doses (since body size changes)
lsarp.redcap.data.loaded = lsarp.redcap.data.loaded %>%
  group_by(study_id, ffq.date) %>%
  mutate(rs_dose = mean(na.omit(c(dose,dose_2,dose_3)))) %>%
  mutate(rs_perc = mean(na.omit(c(percent_product,percent_product_2,percent_product_3)))) %>%
  mutate(rs_intake = rs_dose * rs_perc/100) %>% data.frame()

# check m2 == caloric intake
cor.test(distinct(subset(lsarp.redcap.data.loaded, study_id != "HM0960.00")[,c("study_id", "m2", "energy")])$m2,
         distinct(subset(lsarp.redcap.data.loaded, study_id != "HM0960.00")[,c("study_id", "m2", "energy")])$energy, method="spearman")
# Rho = 0.31, p = 0.00562 (repeat measures, though)
# insensitive to deleting outlier

# check dose = 7.5 / m2
# we have: m2 and dose (RS product minus digestible starch)
rs.bsa = lsarp.redcap.data.loaded$m2
rs.dose = lsarp.redcap.data.loaded$rs_intake # - (lsarp.redcap.data.loaded$rs_dose*lsarp.redcap.data.loaded$percent_product/100)
plot(rs.bsa*7.5,
     rs.dose)
cor.test(rs.bsa*7.5,
     rs.dose)
# good

# therefore, amioca is dosed at 4 * BSA

# clean and merge with other data (keep only what is needed)
lsarp.redcap.data.loaded$HM = substr(lsarp.redcap.data.loaded$study_id, 1, 6)
# add dates and group
lsarp.redcap.data.loaded = merge(lsarp.redcap.data.loaded,
                                 metadata.lsarp.stool[,c("HM", "rs_start_date", "product_end", "group")] %>% distinct(),
                                    by="HM")

# add calories from RS (2.7 kcal/g) or placebo (4 kcal/g)
# note: amioca was dosed at 4 * BSA, versus 7.5 * BSA. So this can be back-calculated
lsarp.redcap.data.loaded$kcal = lsarp.redcap.data.loaded$energy / 4.184
lsarp.redcap.data.loaded$rs.kcal = 
  ifelse(lsarp.redcap.data.loaded$group == "RS",
         # add calories from RS fraction
         (lsarp.redcap.data.loaded$rs_dose * lsarp.redcap.data.loaded$rs_perc/100 * 2.7) + # cite: Miketinas, 2020 https://pubmed.ncbi.nlm.nih.gov/32840627/
           # plus calories from non-RS fraction
           (lsarp.redcap.data.loaded$rs_dose * (1-lsarp.redcap.data.loaded$rs_perc/100) * 4), # 4 kcal/g carbohydrates
         ## FOR PLACEBO, back-calculate Amioca dose
         # note: RS dose = 7.5 g RS / m2 (BSA)
         # therefore: 7.5 / RS Dose = BSA
         # 4 / BSA = Amioca dose
         # Amioca dose * 4 = Amioca calories
         4 * 4 * lsarp.redcap.data.loaded$m2) # 4 kcal/g carbohydrates for amioca
ggplot(lsarp.redcap.data.loaded,
       aes(x=group, y=rs.kcal))+
  geom_boxplot(width=0.5)+
  geom_point(shape=21, color="white", fill="black", size=3)+
  theme_classic()+
  labs(x=NULL, y="Energy from supplement (kcal)")
# this makes sense because RS was dosed higher than amioca, such that calories from RS were additional
# (BSA also contributes to variance)

# calculate days since starting RS
lsarp.redcap.data.loaded$lsarp.days = as.Date(lsarp.redcap.data.loaded$ffq.date) - as.Date(lsarp.redcap.data.loaded$rs_start_date)

# check energy values
hist(lsarp.redcap.data.loaded$kcal)
# some pretty ridiculous values from mostly 1 participant, but consistent
lsarp.redcap.data.loaded[lsarp.redcap.data.loaded$kcal == max(lsarp.redcap.data.loaded$kcal),]
subset(lsarp.redcap.data.loaded, HM == "HM0960")
# HM0960, all entries are consistent

# replace RS value with 0 for days outside of treatment period
# (add a buffer of 7 days)
lsarp.redcap.data.loaded$rs_intake = ifelse(lsarp.redcap.data.loaded$group == "RS" & lsarp.redcap.data.loaded$ffq.date <= as.Date(lsarp.redcap.data.loaded$product_end)+7 &
                                                 lsarp.redcap.data.loaded$group == "RS" & lsarp.redcap.data.loaded$ffq.date >= as.Date(lsarp.redcap.data.loaded$rs_start_date)-7,
                                               lsarp.redcap.data.loaded$rs_intake, 0)

# add dietary energy + energy from RS/placebo supplement
lsarp.redcap.data.loaded$energy = lsarp.redcap.data.loaded$kcal + lsarp.redcap.data.loaded$rs.kcal
# total RS from diet and supplement
lsarp.redcap.data.loaded$total_rs = lsarp.redcap.data.loaded$resistant_starch + lsarp.redcap.data.loaded$rs_intake
# re-normalize to energy
lsarp.redcap.data.loaded$total_rs_adj = lsarp.redcap.data.loaded$total_rs / lsarp.redcap.data.loaded$energy * 1000
# re-normalize fiber (to total energy intake)
lsarp.redcap.data.loaded$fibre_adj = lsarp.redcap.data.loaded$dietary_fibre / lsarp.redcap.data.loaded$energy * 1000


# numbers
lsarp.redcap.data.loaded[,c("HM", "group")] %>% table()
lsarp.redcap.data.loaded[,c("HM", "group")] %>% distinct() %>% dplyr::select("group") %>% table()


# :: ``` longitudinal fiber plot -----------------------------------------------------------

metadata.lsarp.fiber.plot = ggplot(lsarp.redcap.data.loaded %>% subset(!is.na(energy)),
                                      aes(x=as.numeric(lsarp.days), 
                                          y=fibre_adj))+
  annotate("rect", xmin=0, xmax=max(na.omit(subset(metadata.lsarp.stool, lsarp.on.rs == "on.rs")$lsarp.days)),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="black",  se=F, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.days > 0 & lsarp.days < max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Fiber Intake (g per day / 1000 kcal)")
metadata.lsarp.fiber.plot


# :: ``` average fiber plot -----------------------------------------------------------

metadata.lsarp.average.fiber.plot = ggplot()+
  geom_segment(data=lsarp.redcap.data.loaded %>%
                 group_by(HM) %>%
                 mutate(min.fiber = min(na.omit(fibre_adj)),
                        max.fiber = max(na.omit(fibre_adj)),
                        mean.fiber = mean(na.omit(fibre_adj))), 
               aes(x=min.fiber, y=reorder(HM, mean.fiber),
                   xend=max.fiber, yend=reorder(HM, mean.fiber)),
               linewidth=0.3)+
  geom_point(data=lsarp.redcap.data.loaded %>% group_by(HM) %>% 
               mutate(mean.fib.adj = mean(na.omit(fibre_adj))) %>%
               dplyr::select(HM, mean.fib.adj) %>% distinct(),
             aes(x=mean.fib.adj,
                 y=reorder(HM, mean.fib.adj), 
                 fill=scale(mean.fib.adj)), color="black", shape=21, size=3)+
  coord_flip()+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        axis.text.x=element_blank())+
  labs(x="Average Fiber Intake (g per day / 1000 kcal)", y="Participant", fill = "")
metadata.lsarp.average.fiber.plot


# >> Fiber vs RS ----------------------------------------------------------

lsarp.nutrient.data.fiber.rs.plot = ggplot(lsarp.redcap.data.loaded,
                                              aes(x=fibre_adj, rs_adj))+
  geom_path(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, fill="black", color="white", size=2.5)+
  geom_smooth(method="lm", color="black", se=F)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+
  labs(x="Fiber Intake g / 1000 kcal per day",
       y="RS Intake g / 1000 kcal per day")
lsarp.nutrient.data.fiber.rs.plot


# :: ``` longitudinal rs plot ------------------------------------------------------------

# dietary RS intake
metadata.lsarp.rs.plot = ggplot(lsarp.redcap.data.loaded%>%
                                  # remove outlier
                                  subset(HM != "HM0960"),
                                   aes(x=lsarp.days, 
                                       y=resistant_starch))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="black",  se=F, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.days > 0 & lsarp.days < max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Dietary Resistant Starch Intake (g per day)")
metadata.lsarp.rs.plot

# supplemented RS intake
metadata.lsarp.rs.supp.plot = ggplot(lsarp.redcap.data.loaded %>%
                                     # remove outlier
                                       subset(HM != "HM0960"),
                                        aes(x=lsarp.days, 
                                            y=total_rs))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="black",  se=F, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.days > 0 & lsarp.days < max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.days <= 0, "1", "2")), color = "white", size=3)+  
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  #geom_text(aes(label=study_id))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Total Resistant Starch Intake (g per day)")
metadata.lsarp.rs.supp.plot



# :: impute average fiber -------------------------------------------------

# calculate average (energy-adjusted) fiber
lsarp.redcap.data.loaded = lsarp.redcap.data.loaded %>%
  group_by(HM) %>%
  mutate(ave.fiber = mean(na.omit(fibre_adj)))


# use average (energy-adjusted) fiber intake to impute missing (NA) values
lsarp.redcap.data.loaded = lsarp.redcap.data.loaded %>%
  group_by(HM) %>%
  mutate(adj.fiber = ifelse(is.na(fibre_adj), ave.fiber, fibre_adj))
# double check

ggplot(lsarp.redcap.data.loaded,
       aes(x=lsarp.days, 
           y=adj.fiber))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white",  se=T)+
  geom_point(aes(fill=ifelse(lsarp.days > 0 & lsarp.days < max(subset(metadata.lsarp.stool, lsarp.on.rs=="on.rs")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.days <= 0, "1", "2")), color = "white", size=3)+  
  scale_fill_manual(values=c(2,"grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  facet_wrap(~group)+
  labs(x="Days since starting RS", 
       y="Fiber Intake (g per day / 1000 kcal)")


# add average fiber intake to main mapping
metadata.lsarp.stool = merge(metadata.lsarp.stool, 
                                lsarp.redcap.data.loaded[,c("HM", "ave.fiber")] %>% distinct(), 
                                by="HM", all.x = T)
# do not use adj.fiber, since this value is EITHER ave.fiber or adj.fiber; 
# better to keep consistent across HMs at the expense of higher precision for only some HMs

# note: ave.fiber is used as covariate in omic analyses

# :: ICC ------------------------------------------------------------------

# what is the ICC of dietary fiber intake for RS and placebo?
lmerTest::lmer(na.omit(fibre_adj) ~ lsarp.days + (1|HM), 
               subset(lsarp.redcap.data.loaded, group == "RS")) %>%
          performance::icc()
# ICC = 0.471
lmerTest::lmer(na.omit(fibre_adj) ~ lsarp.days + (1|HM), 
               subset(lsarp.redcap.data.loaded, group == "Placebo")) %>%
  performance::icc()
# ICC = 0.546

# together
lmerTest::lmer(fibre_adj ~ lsarp.days + group + (1|HM), 
               lsarp.redcap.data.loaded) %>%
  performance::icc()
# ICC = 0.477, we'll use this in this the paper

# :: Wilcox ---------------------------------------------------------------


lsarp.fiber.mean.pval = 
  wilcox.test(distinct(subset(metadata.lsarp.stool, group == "RS")[,c("HM", "ave.fiber")])$ave.fiber,
              distinct(subset(metadata.lsarp.stool, group == "Placebo")[,c("HM", "ave.fiber")])$ave.fiber)$p.value
# no difference in fiber intake between groups (on average)

lsarp.fiber.comparison.plot = ggplot(metadata.lsarp.stool[,c("HM", "ave.fiber", "group")] %>% distinct(),
                                        aes(x=group, y=ave.fiber))+
  geom_boxplot(width=0.5, alpha=0.5)+
  geom_boxplot(width=0.5, alpha=0.5, outlier.shape = NA, coef = 0)+
  ggbeeswarm::geom_beeswarm(aes(fill=group), color="white", shape=21, size=3)+
  ylim(c(5,20))+
  annotate(geom="text", y=Inf, x=1.5, vjust=1.5, size=4,
           label=paste("p: ", round(lsarp.fiber.mean.pval, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Average Fiber Intake")+
  labs(x="", y="Average Fiber Intake g / 1000 kcal per day")
lsarp.fiber.comparison.plot

# :: save -----------------------------------------------------------------

write.csv(metadata.lsarp.stool, "./2025_07_24_lsarp_fiber.csv")

# final version
(metadata.lsarp.rs.plot+metadata.lsarp.rs.supp.plot+
    metadata.lsarp.fiber.plot+lsarp.fiber.comparison.plot) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_ffq_fiber_rs.pdf",
         width=10, height=6.5)

# Brief methods write up:
# Fiber intakes were calculated from Monash FFQ and adjusted for calculated energy intake. 
# Values were calculated per stool. Missing values were imputed using the average of 3+ FFQs.

# Resistant starch doses & percentages were averaged per RS phase
# Total energy was adjusted for the RS (2.7 kcal/g) and digestible starch (4 kcal/g)

###

# :: add stool water ----------------------------------------------------------

water.lsarp.stool = read.csv("./2024_11_08_metabolomics_samples_prioritized_PD.csv")[,c("standard.name", "HM", "lsarp.days","stool_date_rec_v2", "stool_water_perc", "lsarp.timing")] %>% distinct()
# take unique
water.lsarp.stool = water.lsarp.stool[,c("standard.name", "stool_water_perc")] %>% distinct()

# append to original mapping, to ensure correct phase label
metadata.lsarp.stool = merge(metadata.lsarp.stool,
                            water.lsarp.stool, by="standard.name")

# note: this is referred to as "stool moisture" in text, not "stool water percent", since other liquids were possibly vaporized during freeze drying

# :: add fecalcal ---------------------------------------------------------

# already added, but let's change the name to "fcal"
colnames(metadata.lsarp.stool)[colnames(metadata.lsarp.stool) == "fecalcal_res"] = "fcal"

# add simpler timing variable
metadata.lsarp.stool$timing = 
  ifelse(grepl("baseline", metadata.lsarp.stool$lsarp.timing), "0M",
         ifelse(grepl("\\.5", metadata.lsarp.stool$lsarp.timing), "0.5M",
                ifelse(grepl(" 1 ", metadata.lsarp.stool$lsarp.timing), "1M",
                       ifelse(grepl(" 2 ", metadata.lsarp.stool$lsarp.timing), "2M",
                              ifelse(grepl(" 3 ", metadata.lsarp.stool$lsarp.timing), "3M",
                                     ifelse(grepl(" 4 ", metadata.lsarp.stool$lsarp.timing), "4M",
                                            ifelse(grepl(" 5 ", metadata.lsarp.stool$lsarp.timing), "5M",
                                                   ifelse(grepl(" 6 ", metadata.lsarp.stool$lsarp.timing), "6M",
                                                          ifelse(grepl(" 7 ", metadata.lsarp.stool$lsarp.timing), "7M",
                                                                 ifelse(grepl(" 8 ", metadata.lsarp.stool$lsarp.timing), "8M",
                                                                        ifelse(grepl(" 9  ", metadata.lsarp.stool$lsarp.timing), "9M",
                                                                               ifelse(grepl(" 10 ", metadata.lsarp.stool$lsarp.timing), "10M",
                                                                                      ifelse(grepl(" 11 ", metadata.lsarp.stool$lsarp.timing), "11M",
                              ifelse(grepl(" 12 ", metadata.lsarp.stool$lsarp.timing), "12M", NA))))))))))))))
# make that a factor (since leading 0's were not used)
metadata.lsarp.stool$timing = factor(metadata.lsarp.stool$timing, levels=c("0M", "0.5M", "1M", "2M", "3M", "4M", "5M", "6M", "7M", "8M", "9M" ,"10M", "11M", "12M"))

# RS_Name is already present; just need to change column title (to match most other code)
colnames(metadata.lsarp.stool)[colnames(metadata.lsarp.stool) == "lsarp.rs"] = "RS_Name"

# add missing fcal (checked on REDCap)
# HM0860-STL-02 fcal not performed
# HM0865-STL-02 fcal not performed
# HM0868-STL-02 fcal not performed
# HM0878-STL-02 fcal not performed
# HM0883-STL-02 fcal not performed

# :: fix edge cases -------------------------------------------------------

# edge case fixes (HM0999-STL-11 is actually Month 12)
metadata.lsarp.stool$phase = ifelse(metadata.lsarp.stool$standard.name == "HM0999-STL-11", "washout",  metadata.lsarp.stool$phase)

# if "lsarp.off" is positive, but lsarp timing is Month 5, consider this "on.rs"
metadata.lsarp.stool$lsarp.on.rs = ifelse(grepl("Month 5", metadata.lsarp.stool$lsarp.timing), "on.rs", metadata.lsarp.stool$lsarp.on.rs)
# metadata.lsarp.stool.asv$lsarp.on.rs = ifelse(grepl("Month 5", metadata.lsarp.stool.asv$lsarp.timing), "on.rs", metadata.lsarp.stool.asv$lsarp.on.rs)
# this only happened once; for HM1035, which was collected 3 days after stopping RS

# create second Group (factored in reverse order; for certain plots)
metadata.lsarp.stool = metadata.lsarp.stool %>%
  mutate(Group = factor(ifelse(group == "Placebo", "Placebo", "RS"), levels=c("RS", "Placebo")))

# fill out NA values
metadata.lsarp.stool = metadata.lsarp.stool %>%
  dplyr::select(-ph_response) %>% 
  group_by(HM) %>%
  mutate(flare.day = ifelse(all(is.na(flare.day)), NA, unique(flare.day[!is.na(flare.day)])),
         days.taken = unique(days.taken[!is.na(days.taken)]))
metadata.lsarp.stool$flare.day = as.numeric(metadata.lsarp.stool$flare.day)
metadata.lsarp.stool$days.taken = as.numeric(metadata.lsarp.stool$days.taken)

# add special marker for 1045 (appendicitis)
metadata.lsarp.stool$exit.day = ifelse(metadata.lsarp.stool$HM == "HM1045", 71, metadata.lsarp.stool$flare.day)
metadata.lsarp.stool$exit.reason = ifelse(metadata.lsarp.stool$HM == "HM1045", "excluded", 
                                          ifelse(!is.na(metadata.lsarp.stool$flare.day), "flare", NA))

# :: ``` stool collections plot -----------------------------------------------

# visualize stools + timings
metadata.lsarp.stool.plot = ggplot(metadata.lsarp.stool,
                                   aes(x=lsarp.days, y=reorder(HM, (exit.day))))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(aes(fill=lsarp.on.rs, 
                 alpha = ifelse(lsarp.days <= (flare.day), "1", "0")),
             shape=21, size=3)+
  geom_point(data=metadata.lsarp.stool[,c("exit.day", "lsarp.days","HM", "group", "exit.reason")] %>% distinct(),
             aes(x=exit.day, y=reorder(HM, (lsarp.days)), shape=exit.reason), size=3)+
  scale_fill_manual(values=c(2,"grey","black"))+
  scale_alpha_manual(values=c(0.25,1))+
  scale_shape_manual(values=c(3,4))+
  #geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), nudge_y=-0.1, size=6)+
  theme_minimal()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group,nrow=1, scales="free")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.plot

# It would appear that analyses conducted during the washout would be underpowered
# (only n = 3 for placebo)

# X = flare.day
# faded = stools collected after flare day (remove; analysis would be confounded by new treatment)

# + = appendicitis

# not all of these are to be included in the analysis

# in the main analysis file, I remove the faded points to avoid confusion (i.e. raising the question "why weren't these analyzed?")

# :: 16S Alpha ------------------------------------------------------------

# calculate alpha measures
lsarp.asv.data.median.richness = vegan::specnumber(lsarp.asv.data.median)
lsarp.asv.data.median.shannon = vegan::diversity(lsarp.asv.data.median, index="shannon")

lsarp.asv.data.median.alpha = data.frame(
  standard.name = names(lsarp.asv.data.median.richness),
  richness = lsarp.asv.data.median.richness,
  shannon = lsarp.asv.data.median.shannon
)
# append to mapping file; and create new asv mapping file (keep samples missing 16S)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool,
                                lsarp.asv.data.median.alpha, by="standard.name", all.x=T) %>% data.frame()
# make rownames standard.name
rownames(metadata.lsarp.stool.asv) = metadata.lsarp.stool.asv$standard.name
metadata.lsarp.stool = metadata.lsarp.stool %>% data.frame()
rownames(metadata.lsarp.stool) = metadata.lsarp.stool$standard.name

nrow(metadata.lsarp.stool.asv)
# 196 (4 missing 16S)
length(which(is.na(metadata.lsarp.stool.asv$shannon)))

# :: 16S Butyrogens I -----------------------------------------------------

# 
amplicon.data.gg138.tax.df = readRDS("./2025_09_10_rs_trial_16s_data_tax_table_gg2.Rds") # note: misnamed "gg2", should be "gg138"

# Use classic definition of butyrogens (Lachno et al); more precisely defined as "potential butyrogens"
# Identify ASV's annotated as these butyrogens, and sum them using gg2 data
# reclassify using rarefied data and gg13.8
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds") %>%
  phyloseq::subset_samples(standard.name %in% metadata.lsarp.stool$standard.name)
# subset to butyrogens
amplicon.data.gg138.tax.df = subset(data.frame(amplicon.data.gg138.tax.df), Family=="f__Lachnospiraceae" | Genus=="g__Blautia" | Genus=="g__Roseburia" | Genus=="g__Eubacterium" | Genus=="g__Ruminococcus" | Genus=="g__Clostridium" | Genus=="g__Faecalibacterium") # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
# apply this filter to phyloseq data
lsarp.phyloseq.butyrogens.i = speedyseq::psmelt(amplicon.data.gg.rare) # gg2 data
lsarp.phyloseq.butyrogens.i = lsarp.phyloseq.butyrogens.i %>%
  subset(OTU %in% amplicon.data.gg138.tax.df$ASV) %>% # reduces from ~10,000 to ~3,000 (all ASVs to butyrogen ASVs)
  group_by(Sample) %>% # group by barcode
  mutate(but = sum(Abundance/50000)) %>% # sum of all ASVs annotated as butyrogen per barcode / rarefaction depth
  group_by(standard.name) %>% # group by stool (i.e. so we can take median of extraction replicates)
  # take median
  mutate(but = median(but)) %>% # take median
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but) %>% distinct() %>% data.frame()
lsarp.phyloseq.butyrogens.i %>% arrange(standard.name)

lsarp.phyloseq.butyrogens.i = lsarp.phyloseq.butyrogens.i[,c("standard.name", "but")]
colnames(lsarp.phyloseq.butyrogens.i) = c("standard.name", "but.i")
# rename but to but.i to distinguish it from but.ii (generated below)

# :: 16S PICRUSt2 -------------------------------------------------------

# Note: can be reperformed, but some code must be run in terminal
run.picrust = F
if(run.picrust == T){
# note: make new picrust2 instance
# use this data (gg2):
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds") %>%
  phyloseq::subset_samples(standard.name %in% metadata.lsarp.stool$standard.name)
# at some point, remove taxa that are not actually present

# prepare seq.table for picrust2
fasta_seqs <- Biostrings::DNAStringSet(data.frame(phyloseq::refseq(amplicon.data.gg.rare))[,1])
names(fasta_seqs) <- rownames(data.frame(phyloseq::refseq(amplicon.data.gg.rare)))  # assign ASV IDs as sequence names
# prepare asv.table for picrust2
lsarp.picrust2.abuntable = data.frame(phyloseq::otu_table(amplicon.data.gg.rare))
# export
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/lsarp.picrust2.seqtable.fasta")
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/predict_SCFA_producers/lsarp.picrust2.seqtable.fasta")
write.table(t(lsarp.picrust2.abuntable), "oars_picrust2/lsarp.picrust2.abuntable.tsv", 
            sep="\t", quote = F, col.names = NA)

# STEP 1: run through default picrust2 (to ensure picrust2 works)
 ##cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2
conda activate oars_picrust2
# In R, may need to install Rcpp, jsonlite, lattice, Matrix, RSpectra, castor
rm -r lsarp_picrust2_out
picrust2_pipeline.py \
-s lsarp.picrust2.seqtable.fasta \
-i lsarp.picrust2.abuntable.tsv \
-o lsarp_picrust2_out \
-p 1

# optionally perform stratification to source functions to taxa
# use EC's
metagenome_pipeline.py -i lsarp.picrust2.abuntable.tsv -m lsarp_picrust2_out/marker_predicted_and_nsti.tsv.gz -f lsarp_picrust2_out/EC_predicted.tsv.gz \
-o lsarp_picrust2_out/EC_metagenome_out --strat_out


# STEP 2: run through Vitals' predict_SCFA_producers to identify butyrogens
# source: https://github.com/ag-vital/predict_SCFA_producers/tree/master
# first: git clone https://github.com/ag-vital/predict_SCFA_producers.git
# then, rename "picrust" folder as "SCFA"
cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2/predict_SCFA_producers
# IF REDO: delete placed_seqs.tre
rm placed_seqs.tre
rm -r placement_working
place_seqs.py -s ../lsarp.picrust2.seqtable.fasta -o placed_seqs.tre -p 1 --intermediate placement_working --ref_dir SCFA
# 97 sequences failed to align
hsp.py -t placed_seqs.tre --observed_trait_table SCFA/SCFA_pathwaydata.txt -o lsarp_SCFA_predicted.tsv -p 1 -m emp_prob -n

# save files to new folder
cp ./lsarp_SCFA_predicted.tsv ../picrust2_saved/lsarp_cd_SCFA_predicted.tsv
cp ../lsarp_picrust2_out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz ../picrust2_saved/lsarp_cd_pred_metagenome_contrib.tsv

}

# :: 16S Butyrogens II (Vital) -------------------------------------------------------

# note: definition of butyrogen = contains acetyl-CoA pathway + (but OR buk)

# import
lsarp.vital.butyrogens = read.csv("./oars_picrust2/picrust2_saved/lsarp_cd_SCFA_predicted.tsv", sep="\t")
# subset to MUST HAVE: "acetylcoa" & MUST HAVE ONE OR BOTH OF: "but" or "buk"
lsarp.vital.butyrogens.acetylcoa = subset(lsarp.vital.butyrogens, acetylcoa == 1)
lsarp.vital.butyrogens.acetylcoa = subset(lsarp.vital.butyrogens.acetylcoa, but == 1 | buk == 1)
# these will be considered "butyrogens"
lsarp.vital.butyrogens.acetylcoa

# return to phyloseq object; subset ASVs to those identified as vital's butyrogens
lsarp.phyloseq.butyrogens = phyloseq::otu_table(amplicon.data.gg.rare)[,lsarp.vital.butyrogens.acetylcoa$sequence]
# remake phyloseq object
lsarp.phyloseq.butyrogens = phyloseq::merge_phyloseq(lsarp.phyloseq.butyrogens, 
                                                     phyloseq::tax_table(amplicon.data.gg.rare), 
                                                     phyloseq::sample_data(amplicon.data.gg.rare))
# format to dataframe
lsarp.phyloseq.butyrogens = speedyseq::psmelt(lsarp.phyloseq.butyrogens)
lsarp.phyloseq.butyrogens$Abundance = lsarp.phyloseq.butyrogens$Abundance / 50000
# visualize which butyrogens are most abundant
lsarp.phyloseq.butyrogens %>%
  group_by(LCA) %>%
  mutate(sum.abun = sum(Abundance)) %>% dplyr::select(LCA, sum.abun) %>% distinct() %>% data.frame() %>%
  slice_max(n=20, sum.abun) %>%
  ggplot(aes(x=reorder(LCA, sum.abun), y=sum.abun)) +
  coord_flip()+scale_y_log10()+
  geom_point(shape=21, aes(fill=sum.abun), size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none")+labs(x="", y="Total Abundance (log10)")
# top are Faecalibacterium, Gemmiger, Agathobacter, Anaerostipes, Anaerobutyricum, Roseburia

# collapse replicates
lsarp.phyloseq.butyrogens.df = lsarp.phyloseq.butyrogens %>%
  group_by(Sample) %>% # group by barcode (extraction replicate)
  mutate(but = sum(Abundance)) %>% # sum up abundance of all butyrogens per barcode
  group_by(standard.name) %>% # group by stool sample
  mutate(but = median(but)) %>% # calculate median of extraction replicates
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but) %>% distinct() %>% data.frame()
# take median
lsarp.phyloseq.butyrogens.df %>% arrange(standard.name)

lsarp.phyloseq.butyrogens.ii = lsarp.phyloseq.butyrogens.df[,c("standard.name", "but")]
colnames(lsarp.phyloseq.butyrogens.ii) = c("standard.name", "but.ii")
# rename but to but.ii, to distinguish from original butyrogen definition (above)

# :: 16S Functional Redundancy ---------------------------------------------

# Note: definition of functional redundancy, using ECs:
# mean(n taxa with function i) / n functions

# i.e. "average number of taxa sharing a function / functional richness"


# SLOW
run.fd = F
if(run.fd == T){
  # using PICRUSt2
  amplicon.data.gg.rare = readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds")
  
  picrust_data = read.csv("./oars_picrust2/picrust2_saved/lsarp_cd_pred_metagenome_contrib.tsv", sep="\t")
  #
  
  dim(picrust_data)
  
  # calculate redundancy per function
  picrust_data_redun = picrust_data %>%
    group_by(sample, function.) %>% # group by barcode and EC function
    mutate(ntaxa = length(unique(taxon))) %>% # calculate the number of ASVs with that function (per barcode)
    dplyr::select(sample, function., ntaxa) %>% distinct()
  
  # average functional redundancy (across functions) per sample
  picrust_data_redun_mean <- picrust_data_redun %>%
    group_by(sample) %>% # group by barcode
    # count n unique functions per sample
    mutate(nfunctions = length(unique(function.))) %>% # calculate the number of ECs (per barcode)
    # nfunctions = functional richness
    mutate(meantaxa = mean(ntaxa)) %>% # calculate the mean number of taxa possessing a function, across all functions
    # meantaxa = average taxa richness per function
    dplyr::select(sample, meantaxa, nfunctions) %>% distinct() %>%
    # fun redundancy = "average taxa richness per function / functional richness"
    mutate(fd = nfunctions) %>%
    mutate(fr = meantaxa) %>%
    dplyr::select(sample, fd, fr)
  # fr = functional redundancy = on average, each function has ~meantaxa contributing to it
  # fd = functional diversity (richness)
  
  colnames(picrust_data_redun_mean)[1] = "dada2.sampleNames"
  
  # link to STL and take median
  lsarp.cd.asv.functionalredundancy = merge(picrust_data_redun_mean,
                                        data.frame(phyloseq::sample_data(amplicon.data.gg.rare))[,c("standard.name", "dada2.sampleNames")],
                                        by="dada2.sampleNames")
  lsarp.cd.asv.functionalredundancy = lsarp.cd.asv.functionalredundancy %>%
    group_by(standard.name) %>%
    mutate(fd = median(fd)) %>% # median fd/fr per sample (across extraction replicates)
    mutate(fr = median(fr)) %>% # median fd/fr per sample (across extraction replicates)
    dplyr::select(standard.name, fd, fr) %>% distinct() %>% data.frame()
  
  # interestingly, fr and fd are normally distributed
  
  saveRDS(lsarp.cd.asv.functionalredundancy, "./2026_01_05_lsarp_cd_functional_redundancy.Rds")
}
lsarp.cd.asv.functionalredundancy = readRDS("./2026_01_05_lsarp_cd_functional_redundancy.Rds")


# :: 16S MLP --------------------------------------------------------------

run.16smlp = F
if(run.16smlp == T){
# load RDP data
amplicon.data.rdp.tax.df = readRDS("./2025_09_10_rs_trial_16s_data_tax_table_rdp.Rds")
# use RDP tax to apply MLP (check internal validation performance)
amplicon.data.rdp.tax.df$OTU = amplicon.data.rdp.tax.df$ASV

# Apply (but also, validate) Bork's Microbial Load Predictor
# subset to desired samples
lsarp.cd.asv.data.for.mlp = speedyseq::psmelt(phyloseq::subset_samples(amplicon.data.gg.rare, 
                                                                       standard.name %in% metadata.lsarp.stool$standard.name))
lsarp.cd.asv.data.for.mlp = lsarp.cd.asv.data.for.mlp[,c("OTU", "Abundance", "standard.name")]
lsarp.cd.asv.data.for.mlp$Abundance = lsarp.cd.asv.data.for.mlp$Abundance / 50000
# take median % of ASVs per sample (across extraction replicates)
lsarp.cd.asv.data.for.mlp = lsarp.cd.asv.data.for.mlp %>%
  group_by(standard.name, OTU) %>%
  mutate(Abundance = median(Abundance)) %>% distinct()
# add taxonomy
lsarp.cd.asv.data.for.mlp = merge(lsarp.cd.asv.data.for.mlp[,c("OTU", "standard.name", "Abundance")],
                              amplicon.data.rdp.tax.df, by="OTU") 

# format lsarp.cd taxa (e.g. all genera, or uc_o/c/f) to match MLP
lsarp.cd.asv.data.for.mlp = lsarp.cd.asv.data.for.mlp %>%
  mutate(RDP_taxa = ifelse(!is.na(Genus), paste("g_", Genus,sep=""),
                           ifelse(!is.na(Family), paste("uc_f_", Family, sep=""),
                                  ifelse(!is.na(Order), paste("uc_o_", Order, sep=""),
                                         ifelse(!is.na(Class), paste("uc_c_", Class, sep=""),
                                                ifelse(!is.na(Phylum), paste("uc_p_", Phylum, sep=""),
                                                       "Bacteria"))))))

# make matrix; sum up glom'ed taxa using RDP annotation (range = 0-68%)
lsarp.cd.asv.data.for.mlp = reshape2::acast(lsarp.cd.asv.data.for.mlp,
                                        standard.name ~ RDP_taxa, value.var="Abundance",
                                        fun.aggregate=sum) %>% as.data.frame()
# Note: I'm doing this in a different order than other 16S processing
# e.g. here: median --> sum (if median is 0, lose that taxa)
# elsewhere: sum --> median (if median is 0, sum may rescue)

# add Shannon
lsarp.cd.asv.data.for.mlp.shannon = vegan::diversity(lsarp.cd.asv.data.for.mlp)
# pseudocount
asv.mlp.min = min(lsarp.cd.asv.data.for.mlp[lsarp.cd.asv.data.for.mlp!=0])/2

# log transform with pseudocount later* (because we need to use the smallest pseudo)

# re-train MLP model
mlp.rdp = readRDS("./oars_mlp/model.16S_rRNA.rds")
# extract training data
mlp.rdp.original.data = mlp.rdp$trainingData
# and, we'll use the source data to reprocess
mlp.rdp.retrain.data = read.csv("./oars_mlp/2025_03_05_vandeputte_otu_unrarefied.csv")
mlp.rdp.retrain.data.load = read.csv("./oars_mlp/2025_03_05_vandeputte_load.csv")
# 1. match samples (add rowname, so we can match later, after rarefying)
rownames(mlp.rdp.retrain.data) = paste("sample", mlp.rdp.retrain.data$X)
rownames(mlp.rdp.retrain.data.load) = paste("sample", mlp.rdp.retrain.data.load$X)
# 2. rarefy data
mlp.rdp.retrain.data = mlp.rdp.retrain.data[!is.na(rowSums(mlp.rdp.retrain.data)),]
set.seed(25)
mlp.rdp.retrain.data = phyloseq::rarefy_even_depth(phyloseq::otu_table(mlp.rdp.retrain.data,taxa_are_rows=F), 
                                                   sample.size=10000, replace=F)
mlp.rdp.retrain.data = mlp.rdp.retrain.data / 10000
mlp.rdp.retrain.data = as.data.frame(mlp.rdp.retrain.data)
# add Shannon (calculated myself)
mlp.rdp.retrain.data.shannon = vegan::diversity(mlp.rdp.retrain.data)

# 3. log2 transform + pseudo
mlp.rdp.retrain.data = as.data.frame(mlp.rdp.retrain.data)
mlp.rdp.pseudo = min(mlp.rdp.retrain.data[mlp.rdp.retrain.data!=0])/2
if(asv.mlp.min < mlp.rdp.pseudo){
  asv.mlp.min.to.use = asv.mlp.min
  print(paste("Use LSARP pseudocount:", asv.mlp.min.to.use))
}else{
  asv.mlp.min.to.use = mlp.rdp.pseudo
  print(paste("Use MLP pseudocount:", asv.mlp.min.to.use))
}
# since my pseudo is smaller, we'll use it instead
# --
# briefly return to lsarp.cd data
lsarp.cd.asv.data.for.mlp = log2(lsarp.cd.asv.data.for.mlp+asv.mlp.min.to.use)
lsarp.cd.asv.data.for.mlp$Shannon = lsarp.cd.asv.data.for.mlp.shannon
# --
# now back to MLP data
mlp.rdp.retrain.data = log2(mlp.rdp.retrain.data+asv.mlp.min.to.use)
mlp.rdp.retrain.data$Shannon = mlp.rdp.retrain.data.shannon
# 4. add bacterial load data
mlp.rdp.retrain.data = merge(mlp.rdp.retrain.data,
                             mlp.rdp.retrain.data.load[,c("Cell_count_per_gram", "X")], by="row.names")
mlp.rdp.retrain.data$X.x = NULL
mlp.rdp.retrain.data$X.y = NULL
rownames(mlp.rdp.retrain.data) = mlp.rdp.retrain.data$Row.names
mlp.rdp.retrain.data$Row.names = NULL
mlp.rdp.retrain.data.load = log10(mlp.rdp.retrain.data$Cell_count_per_gram) # Note: Log10
mlp.rdp.retrain.data$Cell_count_per_gram = NULL
# 5. take overlapping features
dim(mlp.rdp.retrain.data) # 190 taxa
mlp.rdp.retrain.data = mlp.rdp.retrain.data[,colnames(mlp.rdp.retrain.data) %in% colnames(lsarp.cd.asv.data.for.mlp)]
mlp.rdp.retrain.data$load = mlp.rdp.retrain.data.load
dim(mlp.rdp.retrain.data) # 150 taxa remain
mlp.rdp.retrain.data = subset(mlp.rdp.retrain.data, !is.na(load))

# 6. build/validate model

# I will need to re-validate their results using their validation data.

t1 = Sys.time()
set.seed(25)
library("caret")
# xgboost package version errors:
packageVersion("xgboost") # 1.7.11.1 (errors occur if you use old (~1.7.5) or up-to-date (~ 3.0)
#remove.packages("xgboost")
#install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.11.1.tar.gz", repos = NULL, type = "source")

mlp.rdp.model.xgb = caret::train(load ~.,
                                 mlp.rdp.retrain.data,
                                 method = "xgbTree",
                                 trControl = caret::trainControl(method = "cv", number = 10,
                                                                 savePredictions=TRUE),
                                 verbose = FALSE)
t2 = Sys.time()
t2-t1 # 2 min

# compare correlation with original
mlp.rdp.model.xgb.preds = subset(mlp.rdp.model.xgb$pred, 
                                 eta == mlp.rdp.model.xgb$finalModel$tuneValue$eta &
                                   gamma  ==  mlp.rdp.model.xgb$finalModel$tuneValue$gamma &
                                   nrounds  ==  mlp.rdp.model.xgb$finalModel$tuneValue$nrounds &
                                   max_depth  ==  mlp.rdp.model.xgb$finalModel$tuneValue$max_depth &
                                   colsample_bytree ==  mlp.rdp.model.xgb$finalModel$tuneValue$colsample_bytree &
                                   min_child_weight == mlp.rdp.model.xgb$finalModel$tuneValue$min_child_weight &
                                   subsample == mlp.rdp.model.xgb$finalModel$tuneValue$subsample)[,c("pred", "obs")]
cor.test(mlp.rdp.model.xgb.preds$pred, mlp.rdp.model.xgb.preds$obs, method="pearson")
mlp.rdp.retrain.xgb.plot = ggplot(mlp.rdp.model.xgb.preds,
                                  aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Retrained")+labs(x="Prediction", y="Observed")
mlp.rdp.retrain.xgb.plot

mlp.rdp.original.preds = subset(mlp.rdp$pred, 
                                eta == mlp.rdp$finalModel$tuneValue$eta &
                                  gamma  ==  mlp.rdp$finalModel$tuneValue$gamma &
                                  nrounds  ==  mlp.rdp$finalModel$tuneValue$nrounds &
                                  max_depth  ==  mlp.rdp$finalModel$tuneValue$max_depth &
                                  colsample_bytree ==  mlp.rdp$finalModel$tuneValue$colsample_bytree &
                                  min_child_weight == mlp.rdp$finalModel$tuneValue$min_child_weight &
                                  subsample == mlp.rdp$finalModel$tuneValue$subsample)[,c("pred", "obs", "rowIndex")] %>%
  # need to take average of preds across 10x
  group_by(rowIndex) %>%
  mutate(pred = mean(pred)) %>%
  distinct() %>% data.frame()
cor.test(mlp.rdp.original.preds$pred, mlp.rdp.original.preds$obs, method="pearson")
# Pearson Cor = 0.79

mlp.rdp.original.plot = ggplot(mlp.rdp.original.preds,
                               aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Original")+labs(x="Prediction", y="Observed")
mlp.rdp.original.plot
# Cor = 0.79
# Retrained is the same! Let's use it


# apply to lsarp.cd
lsarp.cd.asv.mlp.preds.xgb = predict(mlp.rdp.model.xgb, lsarp.cd.asv.data.for.mlp)
lsarp.cd.asv.mlp.preds.xgb = data.frame(standard.name = gsub("\\.", "-", rownames(lsarp.cd.asv.data.for.mlp)),
                                    load = lsarp.cd.asv.mlp.preds.xgb)

lsarp.cd.asv.mlp.preds.xgb.plot = ggplot(lsarp.cd.asv.mlp.preds.xgb,
                                     aes(x=load, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained XGBoost")+
  labs(x="Predicted Load", y="Density")
lsarp.cd.asv.mlp.preds.xgb.plot

# :: 16S MLP Validation ---------------------------------------------------

# select models
models = c("glmnet", "svmRadial", "xgbTree", "ranger", "pls")

mlp.rdp.retrain.models = do.call(rbind, lapply(models, function(model){
  print(model)
  t1 = Sys.time()
  set.seed(25)
  mlp.rdp.model = caret::train(load ~.,
                               mlp.rdp.retrain.data,
                               method = model,
                               trControl = caret::trainControl(method = "cv", number = 10,
                                                               savePredictions=TRUE),
                               verbose = FALSE,
                               preProcess = c("center", "scale")
  )
  t2 = Sys.time()
  t2-t1 # 2 min
  
  # extract preds per model
  if(model == "glmnet"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, alpha == mlp.rdp.model$bestTune$alpha &
                                   lambda == mlp.rdp.model$bestTune$lambda)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "svmRadial"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, sigma == mlp.rdp.model$bestTune$sigma &
                                   C == mlp.rdp.model$bestTune$C)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "xgbTree"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, eta == mlp.rdp.model$bestTune$eta &
                                   gamma  ==  mlp.rdp.model$bestTune$gamma &
                                   nrounds  ==  mlp.rdp.model$bestTune$nrounds &
                                   max_depth  ==  mlp.rdp.model$bestTune$max_depth &
                                   colsample_bytree ==  mlp.rdp.model$bestTune$colsample_bytree &
                                   min_child_weight == mlp.rdp.model$bestTune$min_child_weight &
                                   subsample == mlp.rdp.model$bestTune$subsample)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "ranger"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, mtry == mlp.rdp.model$bestTune$mtry &
                                   splitrule == mlp.rdp.model$bestTune$splitrule &
                                   min.node.size == mlp.rdp.model$bestTune$min.node.size)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "pls"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, ncomp == mlp.rdp.model$bestTune$ncomp)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  # need to take average of preds across 10x
  
  mlp.rdp.model.preds = mlp.rdp.model.preds %>%
    group_by(rowIndex) %>%
    mutate(pred = mean(pred)) %>%
    distinct() %>% data.frame()
  mlp.rdp.model.preds$model = model
  mlp.rdp.model.preds$time = t2-t1
  mlp.rdp.model.preds
}))

mlp.rdp.retrain.models.plots = ggplot(mlp.rdp.retrain.models %>%
                                        rbind(data.frame(mlp.rdp.original.preds[,c("rowIndex", "pred", "obs")]) %>%
                                                mutate(model = "Original", time = NA)) %>% data.frame() %>%
                                        mutate(model = factor(model, levels=c("Original", "glmnet","pls","svmRadial","ranger","xgbTree"))),
                                      aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~model)+labs(x="Prediction", y="Observed")
mlp.rdp.retrain.models.plots
# RandomForest is comparable

# Correlate predictions
Hmisc::rcorr(reshape2::acast(mlp.rdp.retrain.models,
                             model ~ rowIndex, value.var="pred")%>%t() ,
             type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))
# indeed, the most similar to xgBoost predictions

# build final model (RandomForest) and apply
set.seed(25)
mlp.rdp.model.final = caret::train(load ~.,
                                   mlp.rdp.retrain.data,
                                   method = "ranger",
                                   trControl = caret::trainControl(method = "cv", number = 10,
                                                                   savePredictions=TRUE),
                                   verbose = FALSE,
                                   preProcess = c("center", "scale")
)

# apply to lsarp.cd
lsarp.cd.asv.mlp.preds = predict(mlp.rdp.model.final, lsarp.cd.asv.data.for.mlp)
lsarp.cd.asv.mlp.preds = data.frame(standard.name = gsub("\\.", "-", rownames(lsarp.cd.asv.data.for.mlp)),
                                load.asv = lsarp.cd.asv.mlp.preds)

lsarp.cd.asv.mlp.preds.rf.plot = ggplot(lsarp.cd.asv.mlp.preds,
                                    aes(x=load.asv, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained Ranger")+
  labs(x="Predicted Load", y="Density")
lsarp.cd.asv.mlp.preds.rf.plot
# Much smoother than xgBoost
# Good

# save (so we don't need to keep redoing it)
saveRDS(lsarp.cd.asv.mlp.preds, "./2025_09_15_lsarp_cd_mlp_asv.Rds")
}

lsarp.cd.asv.mlp.preds = readRDS("./2025_09_15_lsarp_cd_mlp_asv.Rds")


# >> MGX ------------------------------------------------------------------

# :: process MGX ----------------------------------------------------------

mgx.lsarp.manifest = read.csv("~/Downloads/humann3_outputs_summary_250502.csv")

# subset manifest to LSARP-CD data samples
mgx.lsarp.manifest = subset(mgx.lsarp.manifest, standard.name %in% metadata.lsarp.stool$standard.name)

# create variable based on MGX completion
mgx.lsarp.manifest$mgx.tax.complete = ifelse(mgx.lsarp.manifest$standard.name %in% subset(mgx.lsarp.manifest, metaphlan.bugs == "completed")$standard.name, "yes", "")
mgx.lsarp.manifest$mgx.path.complete = ifelse(mgx.lsarp.manifest$standard.name %in% subset(mgx.lsarp.manifest, pathabundance == "completed")$standard.name, "yes", "")
mgx.lsarp.manifest$mgx.sequenced = ifelse(mgx.lsarp.manifest$standard.name %in% subset(mgx.lsarp.manifest, !grepl("failed", sample.id))$standard.name, "yes", "")

sum(mgx.lsarp.manifest$mgx.tax.complete == "yes")/ nrow(mgx.lsarp.manifest) # 92%
sum(mgx.lsarp.manifest$mgx.path.complete == "yes") / nrow(mgx.lsarp.manifest) # 91%
sum(mgx.lsarp.manifest$mgx.sequenced == "yes") / nrow(mgx.lsarp.manifest) # 100%

# :: read in MGX tax data -----------------------------------------------------

library("tidyverse")
# list all files
all_files = list.files(path = "~/Downloads/humann3_main_outputs/metaphlan_bugs_list", full.names = TRUE)
# reduce to desired samples
lsarp_desired_samples = subset(mgx.lsarp.manifest, mgx.tax.complete == "yes")$standard.name
lsarp_file_list = all_files[stringr::str_detect(basename(all_files), 
                                           paste(lsarp_desired_samples, collapse = "|"))]

read_metaphlan_file <- function(file_path) {
  # Read the file (assuming tab-delimited, adjust if needed)
  data = read.csv(file_path, sep="\t", header=F) # Read all columns as character
  # remove first few rows (blank + titles); and select columns (NCBI tax ID + "additional species")
  data = data[-c(1:5),c(1,3)]
  # add colnames
  colnames(data) = c("taxa", "abundance")
  data$sample = paste("HM", gsub(".*HM", "", gsub("\\_merged.*", "", file_path)), sep="")
  # optionally, but importantly, filter out to just species level
  data = subset(data, grepl("s__", taxa))
  data = subset(data, !grepl("t__", taxa))
  data$taxa = gsub(".*s__", "", data$taxa)
  
  return(data)
}

# load all metagenomics samples, subset to species
lsarp.mgx.taxa <- do.call(rbind, lapply(lsarp_file_list, read_metaphlan_file))
lsarp.mgx.taxa$abundance = as.numeric(lsarp.mgx.taxa$abundance)

# make matrix (and take average of REPEATEDLY sequenced samples)
lsarp.mgx.taxa = reshape2::acast(lsarp.mgx.taxa, sample ~ taxa, 
                                value.var="abundance", fun.aggregate=mean)
# replace NA with 0
lsarp.mgx.taxa[is.na(lsarp.mgx.taxa)] = 0

dim(lsarp.mgx.taxa) # 181 samples x 1224 taxa

hist(rowSums(lsarp.mgx.taxa))
# thanks for averaging repeatedly sequenced samples, some samples have >100%
# so, rescale to 100
lsarp.mgx.taxa = t(apply(lsarp.mgx.taxa, 1, function(x) {x / sum(x) * 100})) %>% as.matrix()
hist(rowSums(lsarp.mgx.taxa))
# good

# ``` mgx completeness ----------------------------------------------------

ggplot(metadata.lsarp.stool %>% mutate(mgx = ifelse(standard.name %in% rownames(lsarp.mgx.taxa), "mgx", NA)),
       aes(x=lsarp.days, y=reorder(HM, (lsarp.days))))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(color="black", size=3)+
  geom_point(color="white", size=2)+
  geom_point(aes(color=lsarp.on.rs))+
  geom_point(aes(x=lsarp.days, 
                 y=reorder(HM, (lsarp.days)), shape=mgx), size=3)+
  scale_color_manual(values=c(2,"grey","black"))+
  #geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), nudge_y=-0.1, size=6)+
  theme_minimal()+
  labs(x="Day since starting product", y="", color="")+
  facet_wrap(~group, scales="free")

# looks good; very few samples missing

# >> MGX MLP --------------------------------------------------------------

run.mgxmlp = T
if(run.mgxmlp == T){
# now apply MLP to MGX data
lsarp.cd.mgx.data.for.mlp = as.data.frame(lsarp.mgx.taxa)
# rescale to %
lsarp.cd.mgx.data.for.mlp <- lsarp.cd.mgx.data.for.mlp / 100

# calculate Shannon diversity
lsarp.cd.mgx.data.for.mlp.shannon = vegan::diversity((lsarp.cd.mgx.data.for.mlp), index="shannon")

# extract MGX MLP data
mlp.mp3 = readRDS("./oars_mlp/model.metaphlan3.rds")
mlp.mp3.data = mlp.mp3$trainingData
# undo log transform (they say log10 in their methods)
mlp.mp3.data = 10^mlp.mp3.data
mlp.mp3.data.load = log10(mlp.mp3.data$.outcome) # undo the unlog on bacterial load column
mlp.mp3.data.shannon = log10(mlp.mp3.data$`Shannon diversity`) # undo the log on shannon column and save
mlp.mp3.data$`.outcome` = NULL
mlp.mp3.data$`Shannon diversity` = NULL

# replace min value (which was calculated per taxa) with 0; but ignore Shannon and Load columns
mlp.mp3.data <- do.call(cbind, lapply(1:ncol(mlp.mp3.data), function(x) {
  col = data.frame(mlp.mp3.data[,x])
  col[col == min(col, na.rm = TRUE)] <- 0
  # and convert to %
  col = col/100
  col.df = data.frame(feature = col)
  colnames(col.df) = colnames(mlp.mp3.data)[x]
  col.df
})) %>% data.frame()
rowSums(mlp.mp3.data)
# these values are way off of what'd be expected
# I would therefore consider it inappropriate to use these data
# alternatively, proceed by rescaling to %
mlp.mp3.data <- sweep(mlp.mp3.data, 1, rowSums(mlp.mp3.data), FUN = "/") * 1

# clean up taxa names
colnames(mlp.mp3.data) = gsub(".*s__", "", colnames(mlp.mp3.data))

# take only overlapping taxa
dim(mlp.mp3.data) # 156 taxa
mlp.mp3.data = mlp.mp3.data[,colnames(mlp.mp3.data) %in% c(colnames(lsarp.cd.mgx.data.for.mlp), ".outcome")] 
dim(mlp.mp3.data) # 98 taxa (100 - 2)

# log2 transform + pseudocount
mgx.mlp.min = min(mlp.mp3.data[mlp.mp3.data!=0])/2
mgx.lsarp.cd.mlp.min = min(lsarp.cd.mgx.data.for.mlp[lsarp.cd.mgx.data.for.mlp!=0])/2
if(mgx.lsarp.cd.mlp.min < mgx.mlp.min){
  mgx.mlp.min.to.use = mgx.lsarp.cd.mlp.min
  print(paste("Use lsarp.cd pseudocount:", mgx.mlp.min.to.use))
}else{
  mgx.mlp.min.to.use = mgx.mlp.min
  print(paste("Use MLP pseudocount:", mgx.mlp.min.to.use))
}
mlp.mp3.data = log2(mlp.mp3.data+mgx.mlp.min.to.use)
# add Shannon and Load cols back
mlp.mp3.data$load = mlp.mp3.data.load
mlp.mp3.data$Shannon = mlp.mp3.data.shannon
# and finish off lsarp.cd data
lsarp.cd.mgx.data.for.mlp = log2(lsarp.cd.mgx.data.for.mlp+mgx.mlp.min.to.use)
lsarp.cd.mgx.data.for.mlp$Shannon = lsarp.cd.mgx.data.for.mlp.shannon

# train model
t1 = Sys.time()
mlp.mp3.model.xgb = caret::train(load ~.,
                                 mlp.mp3.data,
                                 method = "xgbTree",
                                 trControl = caret::trainControl(method = "cv", number = 10, # note: they perform 5x 10x
                                                                 savePredictions=TRUE),
                                 verbose = FALSE
)
t2 = Sys.time()
t2-t1 # 5 min

# compare correlation with original
mlp.mp3.model.xgb.preds = subset(mlp.mp3.model.xgb$pred, 
                                 eta == mlp.mp3.model.xgb$finalModel$tuneValue$eta &
                                   gamma  ==  mlp.mp3.model.xgb$finalModel$tuneValue$gamma &
                                   nrounds  ==  mlp.mp3.model.xgb$finalModel$tuneValue$nrounds &
                                   max_depth  ==  mlp.mp3.model.xgb$finalModel$tuneValue$max_depth &
                                   colsample_bytree ==  mlp.mp3.model.xgb$finalModel$tuneValue$colsample_bytree &
                                   min_child_weight == mlp.mp3.model.xgb$finalModel$tuneValue$min_child_weight &
                                   subsample == mlp.mp3.model.xgb$finalModel$tuneValue$subsample)[,c("pred", "obs")]
cor.test(mlp.mp3.model.xgb.preds$pred, mlp.mp3.model.xgb.preds$obs, method="pearson")
mlp.mgx.retrain.xgb.plot = ggplot(mlp.mp3.model.xgb.preds,
                                  aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Retrained")+labs(x="Prediction", y="Observed")
mlp.mgx.retrain.xgb.plot# from 0.61 to 0.53 after re-processing (including reducing features)

# calculate original (from paper)
mlp.mp3.original.preds = subset(mlp.mp3$pred, 
                                eta == mlp.mp3$finalModel$tuneValue$eta &
                                  gamma  ==  mlp.mp3$finalModel$tuneValue$gamma &
                                  nrounds  ==  mlp.mp3$finalModel$tuneValue$nrounds &
                                  max_depth  ==  mlp.mp3$finalModel$tuneValue$max_depth &
                                  colsample_bytree ==  mlp.mp3$finalModel$tuneValue$colsample_bytree &
                                  min_child_weight == mlp.mp3$finalModel$tuneValue$min_child_weight &
                                  subsample == mlp.mp3$finalModel$tuneValue$subsample)[,c("pred", "obs", "rowIndex")] %>%
  # need to take average of preds across 10x
  group_by(rowIndex) %>%
  mutate(pred = mean(pred)) %>%
  distinct() %>% data.frame()
cor.test(mlp.mp3.original.preds$pred, mlp.mp3.original.preds$obs, method="pearson")
mlp.mgx.original.plot = ggplot(mlp.mp3.original.preds,
                               aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Original")+labs(x="Prediction", y="Observed")
mlp.mgx.original.plot

# apply to lsarp.cd
lsarp.cd.mgx.mlp.preds = predict(mlp.mp3.model.xgb, lsarp.cd.mgx.data.for.mlp)
lsarp.cd.mgx.mlp.preds = data.frame(standard.name = gsub("\\.", "-", rownames(lsarp.cd.mgx.data.for.mlp)),
                                load.mgx = lsarp.cd.mgx.mlp.preds)

hist(lsarp.cd.mgx.mlp.preds$load.mgx)


# :: MGX MLP Validation ---------------------------------------------------

# Like with ASVs, validate other models for MGX MLP

# select models
models = c("glmnet", "svmRadial", "xgbTree", "ranger", "pls")

mlp.mgx.retrain.models = do.call(rbind, lapply(models, function(model){
  print(model)
  t1 = Sys.time()
  set.seed(25)
  mlp.mgx.model = caret::train(load ~.,
                               mlp.mp3.data %>% as.data.frame(),
                               method = model,
                               trControl = caret::trainControl(method = "cv", number = 10,
                                                               savePredictions=TRUE),
                               verbose = FALSE,
                               preProcess = c("center", "scale")
  )
  t2 = Sys.time()
  t2-t1 # 2 min
  
  # extract preds per model
  if(model == "glmnet"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, alpha == mlp.mgx.model$bestTune$alpha &
                                   lambda == mlp.mgx.model$bestTune$lambda)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "svmRadial"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, sigma == mlp.mgx.model$bestTune$sigma &
                                   C == mlp.mgx.model$bestTune$C)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "xgbTree"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, eta == mlp.mgx.model$finalModel$tuneValue$eta &
                                   gamma  ==  mlp.mgx.model$finalModel$tuneValue$gamma &
                                   nrounds  ==  mlp.mgx.model$finalModel$tuneValue$nrounds &
                                   max_depth  ==  mlp.mgx.model$finalModel$tuneValue$max_depth &
                                   colsample_bytree ==  mlp.mgx.model$finalModel$tuneValue$colsample_bytree &
                                   min_child_weight == mlp.mgx.model$finalModel$tuneValue$min_child_weight &
                                   subsample == mlp.mgx.model$finalModel$tuneValue$subsample)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "ranger"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, mtry == mlp.mgx.model$bestTune$mtry &
                                   splitrule == mlp.mgx.model$bestTune$splitrule &
                                   min.node.size == mlp.mgx.model$bestTune$min.node.size)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "pls"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, ncomp == mlp.mgx.model$bestTune$ncomp)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  # need to take average of preds across 10x
  
  mlp.mgx.model.preds = mlp.mgx.model.preds %>%
    group_by(rowIndex) %>%
    mutate(pred = mean(pred)) %>%
    distinct() %>% data.frame()
  mlp.mgx.model.preds$model = model
  mlp.mgx.model.preds$time = t2-t1
  mlp.mgx.model.preds
}))

mlp.mgx.retrain.models.plots = ggplot(mlp.mgx.retrain.models %>%
                                        rbind(data.frame(mlp.mp3.original.preds[,c("rowIndex", "pred", "obs")]) %>%
                                                mutate(model = "Original", time = NA)) %>% data.frame() %>%
                                        mutate(model = factor(model, levels=c("Original", "glmnet","pls","svmRadial","ranger","xgbTree"))),
                                      aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~model)+labs(x="Prediction", y="Observed")
mlp.mgx.retrain.models.plots

# Correlate predictions
Hmisc::rcorr(reshape2::acast(mlp.mgx.retrain.models,
                             model ~ rowIndex, value.var="pred")%>%t() ,
             type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))
# ranger again is most comparable


# :: ```plot: MLP comparison -------------------------------------------------

# comparison of predictions
mlp.rdp.mgx.comparison.plot = ggplot(merge(lsarp.cd.asv.mlp.preds,
                                           lsarp.cd.mgx.mlp.preds, by="standard.name"),
                                     aes(x=load.asv, y=load.mgx))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"ASV vs MGX")+labs(x="ASV Load", y="MGX Load")

# ASV performs better, and they correlate
library("patchwork")
(mlp.rdp.retrain.models.plots / 
    mlp.mgx.retrain.models.plots) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_mlp_asv_mgx.pdf",
         width=8, height=8,device = cairo_pdf)

# plus
(lsarp.cd.asv.mlp.preds.xgb.plot/
    lsarp.cd.asv.mlp.preds.rf.plot/
    mlp.rdp.mgx.comparison.plot)%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_mlp_retrain.pdf",
         width=4, height=8,device = cairo_pdf)

# save (so we don't need to keep redoing it)
saveRDS(lsarp.cd.mgx.mlp.preds, "./2025_09_15_lsarp_cd_mlp_mgx.Rds")
}

lsarp.cd.mgx.mlp.preds = readRDS("./2025_09_15_lsarp_cd_mlp_mgx.Rds")


# >> MPX ----------------------------------------------------------------

lsarp.mpx = read.csv("./metaproteomics/MetaLab_iterative/final_proteins.tsv", sep="\t")
dim(lsarp.mpx)
# clean and filter to LSARP samples
colnames(lsarp.mpx)[c(1:9)] # keep 1 (protein ID), 8 and onwards (samples)
colnames(lsarp.mpx)[c(499:504)] # remove last 5 (non-samples)
lsarp.mpx = lsarp.mpx[,c(1, 8:499)]

# clean colnames (to be consistent with one another)
colnames(lsarp.mpx) = gsub("Intensity\\.\\d{4}_\\d{2}_\\d{2}_Sara_StintziLab_", "", colnames(lsarp.mpx))
colnames(lsarp.mpx) = gsub("Intensity\\.\\d{4}_\\d{2}_\\d{2}_Sara_Stintzi_", "", colnames(lsarp.mpx))

# filter to lsarp samples
lsarp.mpx.lsarp.only = lsarp.mpx[,colnames(lsarp.mpx) %in% c("Protein.IDs", "QC",
                                                         gsub("-", "_", metadata.lsarp.stool$standard.name))]
colnames(lsarp.mpx.lsarp.only) = gsub("\\.1", "",colnames(lsarp.mpx.lsarp.only))
dim(lsarp.mpx.lsarp.only)

# load functional annotations
mpx.functions = read.csv("./metaproteomics/MetaLab_iterative/functional_annotation/functions.tsv", sep="\t", header=F)
colnames(mpx.functions) = mpx.functions[1,] # make first rows into colnames
mpx.functions = mpx.functions[-1,] # and then delete them (there's a more efficient way of doing this)
# keep useful annotations
mpx.functions = mpx.functions[,colnames(mpx.functions) %in% c("Group_ID", "Name", "Protein name",
                                                              "Description", "Taxonomy Id", "Taxonomy name",
                                                              "Preferred name", "Gene_Ontology_id", "Gene_Ontology_name",
                                                              "Gene_Ontology_namespace", "EC_id", "EC_de", "EC_an",
                                                              "EC_ca", "KEGG_ko", "KEGG_Pathway_Entry", "KEGG_Pathway_Name",
                                                              "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
                                                              "CAZy", "BiGG_Reaction", "PFAMs", "COG accession", "COG category",
                                                              "COG name", "NOG accession", "NOG category", "NOG name")] %>% distinct()
# standardize spaces in colnames (use "_", not " ")
colnames(mpx.functions) = gsub(" ", "_", colnames(mpx.functions))

# create function to collapse proteins by annotation (sum up proteins sharing functional ID)
protein.collapser = function(data,
                             protein.information,
                             annotation,
                             add_tax = FALSE,
                             delimiter = ";"){
  # annotation = bquote(annotation) # dplyr cannot handle quotes
  data.melt <- reshape2::melt(data)
  # remove rows with 0 abundance
  data.melt <- subset(data.melt, value != 0)
  # set sample column name
  colnames(data.melt)[2] <- "code"
  
  # select annotation of choice, and lengthen annotation file, such that X;Y;Z becomes X, Y, Z on separate rows
  # optionally append Taxonomy_name; //defunct
  if(add_tax == T & annotation == "CAZy"){
    # paste Taxonomy_name to each annotation
    protein.information <- protein.information %>%
      dplyr::select(Name, CAZy, Taxonomy_name)%>%
      rowwise() %>%
      mutate(
        CAZy = paste0(str_split(CAZy, ",")[[1]], "_", Taxonomy_name) %>%
          paste(collapse = ",")
      ) %>%
      ungroup()
  }
  if(add_tax == T & annotation == "Preferred_name"){
    # paste Taxonomy_name to each annotation
    protein.information <- protein.information %>%
      dplyr::select(Name, Preferred_name, Taxonomy_name)%>%
      rowwise() %>%
      mutate(
        Preferred_name = paste0(str_split(Preferred_name, ",")[[1]], "_", Taxonomy_name) %>%
          paste(collapse = ",")
      ) %>%
      ungroup()
  }
  # select annotation to use
  protein.information.long = protein.information[,colnames(protein.information) %in% c("Name", annotation)]
  
  # identify the maximum number of splits (i.e. max number of functions a protein is annotated with)
  max_splits <- paste("split", 1:(max(str_count(protein.information.long[,2], delimiter)+1)))
  # separate these functions (per sample) into n columns, where n is max_splits
  protein.information.long = tidyr::separate(protein.information.long, col={{annotation}}, into=max_splits, sep=delimiter, remove=T)
  # melt to df
  protein.information.long = reshape2::melt(protein.information.long, id="Name")
  # replace blank data with NA
  protein.information.long$value = ifelse(protein.information.long$value == "", NA, protein.information.long$value)
  protein.information.long$variable = NULL
  colnames(protein.information.long)[2] <- "annotation"
  # remove NA values
  protein.information.long = subset(protein.information.long, !is.na(annotation))
  # append to data
  colnames(data.melt)[1] <- "Name"
  data.melt <- merge(data.melt, 
                     protein.information.long, by="Name")
  # EDGECASE: fix Purine"-"/" "nucleoside (sum them together)
  data.melt$annotation <- ifelse(data.melt$annotation == "Purine-nucleoside phosphorylase",
                                 "Purine nucleoside phosphorylase", data.melt$annotation)
  
  colnames(data.melt)[2] = "code"
  # sum proteins by annotation per sample (code)
  data.melt.meta.summed <- data.melt %>%
    group_by(code, annotation) %>% # group by stool sample
    mutate(sum.intensity = sum(value)) %>%
    dplyr::select(code, annotation, sum.intensity) %>% distinct() %>% data.frame()
  
  # done; ready for downstream analyses
  data.melt.meta.summed
}
# checked


# :: MPX KEGG PATHWAY --------------------------------------------------------------

# apply function to collapse protein data by KEGG pathway
lsarp.cd.mpx.kegg = protein.collapser(data=lsarp.mpx.lsarp.only, 
                                            protein.information = mpx.functions,
                                            annotation="KEGG_Pathway_Name",
                                            delimiter = ";")
# make a matrix; take average of extraction replicates
lsarp.cd.mpx.kegg = reshape2::acast(data.frame(lsarp.cd.mpx.kegg) %>% 
                                      mutate(code = as.character(code)) %>% distinct(),
                                              code ~ annotation, value.var="sum.intensity",
                                              fun.aggregate = mean)
# return "_" to "-" to make rows into standard.names
rownames(lsarp.cd.mpx.kegg) = gsub("_", "-", rownames(lsarp.cd.mpx.kegg))
# save
saveRDS(lsarp.cd.mpx.kegg, "./metaproteomics/2025_06_28_lsarp.cd_mpx_kegg.Rds")

lsarp.cd.mpx.kegg.mat = readRDS("./metaproteomics/2025_06_28_lsarp.cd_mpx_kegg.Rds")
dim(lsarp.cd.mpx.kegg.mat) # 183 samples x 181 KEGG Pathways


# :: MPX COG --------------------------------------------------------------

# apply function to collapse protein data by COG annotation
lsarp.cd.mpx.cog = protein.collapser(data=lsarp.mpx.lsarp.only, 
                                     protein.information = mpx.functions,
                                     annotation="COG_name")
# make a matrix; take average of extraction replicates
lsarp.cd.mpx.cog = reshape2::acast(data.frame(lsarp.cd.mpx.cog) %>% 
                                     mutate(code = as.character(code)) %>% distinct(),
                                   code ~ annotation, value.var="sum.intensity",
                                   fun.aggregate = mean)
# return "_" to "-" to make rows into standard.names
rownames(lsarp.cd.mpx.cog) = gsub("_", "-", rownames(lsarp.cd.mpx.cog))
# save
saveRDS(lsarp.cd.mpx.cog, "./metaproteomics/2025_06_28_lsarp.cd_mpx_cog.Rds")

lsarp.cd.mpx.cog.mat = readRDS("./metaproteomics/2025_06_28_lsarp.cd_mpx_cog.Rds")
dim(lsarp.cd.mpx.cog.mat) # 183 samples x 2696 COGs


# :: MPX CAZy -------------------------------------------------------------

# apply function to collapse protein data by CAZy annotation
lsarp.cd.mpx.cazy = protein.collapser(data = lsarp.mpx.lsarp.only, 
                                      protein.information = mpx.functions,
                                      annotation = "CAZy",
                                      delimiter = ",")

# make a matrix; take average of extraction replicates
lsarp.cd.mpx.cazy = reshape2::acast(data.frame(lsarp.cd.mpx.cazy) %>% 
                                      mutate(code = as.character(code)) %>% distinct(),
                                    code ~ annotation, value.var="sum.intensity",
                                    fun.aggregate = mean)

# remove "-" in CAZy names
lsarp.cd.mpx.cazy = lsarp.cd.mpx.cazy[,colnames(lsarp.cd.mpx.cazy) != "-"]
# return "_" to "-" to make rows into standard.names
rownames(lsarp.cd.mpx.cazy) = gsub("_", "-", rownames(lsarp.cd.mpx.cazy))

saveRDS(lsarp.cd.mpx.cazy, "./metaproteomics/2025_06_28_lsarp.cd_mpx_cazy.Rds")

lsarp.cd.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_lsarp.cd_mpx_cazy.Rds")
dim(lsarp.cd.mpx.cazy.mat) # 183 samples x 66 CAZy


# :: MPX CAZy Starch:Mucin ----------------------------------------------------

lsarp.cd.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_lsarp.cd_mpx_cazy.Rds")

# https://pmc.ncbi.nlm.nih.gov/articles/PMC9120202/ for GH's involved in mucin degradation

mucin.cazy = c("GH101", "GH20", "GH29", "GH33", "GH84", "GH95")
starch.cazy = c("GH13", "CBM48", "CBM20", "GH77")

# calculate Starch:Mucin ratio per sample
lsarp.cd.mpx.starch.mucin.mat = lsarp.cd.mpx.cazy.mat %>%
  reshape2::melt() %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>% # replace NA with 0
  mutate(value = ifelse(value == 0, min(value[value!=0])/2, value)) %>% # replace 0 with pseudocount
  mutate(type = ifelse(Var2 %in% starch.cazy, "starch",
                       ifelse(Var2 %in% mucin.cazy, "mucin", "other"))) %>%
  subset(type != "other") %>%
  group_by(Var1, type) %>% # group by CAZy type (starch- or mucin- active)
  mutate(sum.cazy = sum(value)) %>% # sum up proteins per CAZy group
  dplyr::select(Var1, sum.cazy, type) %>% distinct() %>%
  reshape2::acast(Var1 ~ type, value.var="sum.cazy") %>% # make matrix (using summed CAZy groups)
  data.frame() %>%
  mutate(starch.mucin = log2(starch / mucin)) %>%  # calculate log ratio
  mutate(standard.name = rownames(.)) %>%
  dplyr::select(standard.name, starch, mucin, starch.mucin)
# fix standard name, if necessary
lsarp.cd.mpx.starch.mucin.mat$standard.name = gsub("\\_", "\\-", lsarp.cd.mpx.starch.mucin.mat$standard.name)


# >> MBX ------------------------------------------------------------------


# :: load data ------------------------------------------------------------

# feature table
lsarp.mbx.data = read.csv("./metabolomics/lsarp_mbx_data_raw.csv")
dim(lsarp.mbx.data) # 273 samples, 6863 features

# feature map
lsarp.mbx.feature.data = read.csv("./metabolomics/lsarp_mbx_feature_data.csv")

# mbx map
lsarp.mbx.metadata <- read.csv("./2024_11_08_metabolomics_samples_prioritized_PD.csv")

# mbx sample map
lsarp.mbx.samples = read.csv("~/Documents/PhD/git_oars_archfolder/metabolomics/LSARP_naming_HM.csv")

# main map
metadata.lsarp.stool 


# :: process --------------------------------------------------------------

# how many samples
ncol(lsarp.mbx.data) - 6
# n = 267 samples

# add rownames (bucket = features)
rownames(lsarp.mbx.data) = lsarp.mbx.data$Bucket
# remove feature ID data (Bucket, m.z, M, RT min, Ions, MS/MS)
lsarp.mbx.data[,c(1:6)]
lsarp.mbx.data = lsarp.mbx.data[,-c(1:6)]

# fix sample names
colnames(lsarp.mbx.data) = lsarp.mbx.samples$standard.name[match(colnames(lsarp.mbx.data), make.names(lsarp.mbx.samples[,1]))]

# subset to LSARP
lsarp.mbx.metadata = subset(lsarp.mbx.metadata, standard.name %in% metadata.lsarp.stool$standard.name)
length(unique(lsarp.mbx.metadata$standard.name))
# n=196 samples

# annotations
lsarp.mbx.feature.data %>% nrow()
subset(lsarp.mbx.feature.data, Annotation.choisie != "")$Annotation.choisie %>% unique() %>% length()
# 1987 total unique features
subset(lsarp.mbx.feature.data, Annotation.choisie != "")$Annotation.choisie %>% unique() %>% sort()
# 1987 unique annotated features

# Save 2 versions:
# 1: annotated (DO NOT collapse (sum) features annotated the same, per AS advice)
# 2: raw (keep all features, annotate where possible)

lsarp.mbx.raw = lsarp.mbx.data %>% as.matrix() %>% reshape2::melt()
colnames(lsarp.mbx.raw) = c("Bucket", "standard.name", "value")
lsarp.mbx.raw = merge(lsarp.mbx.raw,
                     lsarp.mbx.feature.data[,c("Bucket", "Annotation.choisie")], by="Bucket")

# to keep unannotated features, annotate as the RT:m/z
lsarp.mbx.raw$Annotation = ifelse(as.character(lsarp.mbx.raw$Annotation.choisie) != "", 
                                 paste(as.character(lsarp.mbx.raw$Annotation.choisie), as.character(lsarp.mbx.raw$Bucket), sep=" | "),
                                 as.character(lsarp.mbx.raw$Bucket))
# good

unique(lsarp.mbx.raw$Annotation) %>% length()
# 6863 unique 

# now, slice off unannotated (only annotated have "|", so remove those without)
lsarp.mbx.annotated = subset(lsarp.mbx.raw, grepl(" \\| ", Annotation))
unique(lsarp.mbx.annotated$Annotation) %>% length()
# 3054 features

# make matrix
lsarp.mbx.raw.mat = reshape2::acast(lsarp.mbx.raw,
                                   standard.name ~ Annotation, value="value") %>% as.data.frame()
dim(lsarp.mbx.raw.mat) # 267 x 6863

lsarp.mbx.annotated.mat = reshape2::acast(lsarp.mbx.annotated,
                                         standard.name ~ Annotation, value="value") %>% as.data.frame()
dim(lsarp.mbx.annotated.mat) # 267 x 3054

saveRDS(lsarp.mbx.raw.mat, "./metabolomics/2025_11_16_lsarp_mbx_raw.Rds")
saveRDS(lsarp.mbx.annotated.mat, "./metabolomics/2025_11_16_lsarp_mbx_annotated.Rds")

# :: filter ---------------------------------------------------------------

# note: prevalence filtration must be applied! But, it is applied in analysis script
# since we'll use 80% prevalence filtration,
# we want to add flexibility for removing non-compliant samples first (or not)


# >> RapidAIM pH ----------------------------------------------------------

# :: load pH data ---------------------------------------------------------

# read in rapidaim scores
all.rapidaim.ph = read.csv("~/Documents/PhD/Dissertation/RS_mapping_files/2025_03_11_full_mapping_rapidaim.csv")

# calculate delta.ph
lsarp.cd.rapidaim.scores <- all.rapidaim.ph %>%
  mutate(HM_no0 = substr(HM, 1, 6)) %>% # clean up HM names
  subset(HM_no0 %in% lsarp.patient.list$HM) %>% # subset to LSARP patients
  subset(!grepl(paste(c("\\.01","\\.02", "\\.03", "\\.04", "\\.05", "\\.06"), collapse="|"), HM)) %>% # remove subsequent+ collections
  # subset to RS + PBS
  subset(RS_Name %in% rs.names.pbs) %>%
  # calculate median pH of up to 3 culture replicates
  group_by(HM, RS_Name) %>%
  mutate(med.ph = median(pH)) %>%
  dplyr::select(HM, RS_Name, med.ph) %>% distinct() %>%
  # calculate delta pH relative to PBS
  group_by(HM)%>%
  mutate(delta.ph = med.ph - med.ph[RS_Name == "PBS"]) %>%
  subset(RS_Name != "PBS") %>%
  mutate(study_id = HM) 

unique(lsarp.cd.rapidaim.scores$HM)
# good

# >> RS Selections ---------------------------------------------------------

# load original selection data (rather than re-calculate again, like with OARS)

# CD
unique(metadata.lsarp.stool$HM)

lsarp.cd.scores.files = list.files("./lsarp_rs_scores/")
lsarp.cd.scores.files = lsarp.cd.scores.files[grep("scores", lsarp.cd.scores.files)]

# note: I copied the scores files into a new folder
lsarp.cd.rs.selections = do.call(rbind, lapply(lsarp.cd.scores.files, function(x){
  data.loaded = readRDS(paste("./lsarp_rs_scores/", x, sep=""))
  data.loaded$sample = substr(x, 1, 6)
  data.loaded
}))

lsarp.cd.rs.selections
# identify RS selected
lsarp.cd.rs.selected = lsarp.cd.rs.selections %>%
  group_by(sample) %>%
  subset(Z_score >=1) %>%
  group_by(sample) %>%
  filter(Distance==min(Distance))
colnames(lsarp.cd.rs.selected)[1] = "Selected"

# note: 860, 865 and 878 used MPX as well, so fix these exceptions!
lsarp.cd.rs.selections.exceptions.1 = lsarp.cd.rs.selections %>%
  subset(sample == "HM0860" & RS == "HiMaize260")
lsarp.cd.rs.selections.exceptions.2 = lsarp.cd.rs.selections %>%
  subset(sample == "HM0878" & RS == "Novelose330")
  
lsarp.cd.rs.selections.exceptions = rbind(lsarp.cd.rs.selections.exceptions.1,
                                          lsarp.cd.rs.selections.exceptions.2) %>% as.data.frame()
colnames(lsarp.cd.rs.selections.exceptions)[1] = "Selected"

lsarp.cd.rs.selections$selected_rs = lsarp.cd.rs.selections$Selected
lsarp.cd.rs.selections$Selected = NULL

# for mapping file
lsarp.cd.rs.selected = subset(lsarp.cd.rs.selected, !sample %in% c("HM0860","HM0878")) %>%
  rbind(lsarp.cd.rs.selections.exceptions) %>% as.data.frame()

# for RS selections plot
lsarp.cd.rs.selections.all = merge(lsarp.cd.rs.selections, 
                                   lsarp.cd.rs.selected[,c("sample", "Selected")], by="sample", all.x=T)

saveRDS(lsarp.cd.rs.selections.all, "./lsarp_rs_scores/2025_06_29_lsarp.cd.rs.selections.RDS")
    


# >> RapidAIM ASV -------------------------------------------------------------

# :: Fermentation Response ------------------------------------------------

# ifelse > 1.27 (per OARS)

lsarp.cd.rapidaim.scores$HM = substr(lsarp.cd.rapidaim.scores$HM, 1, 6)
lsarp.cd.rapidaim.scores$sample = lsarp.cd.rapidaim.scores$HM

lsarp.cd.rs.selections

lsarp.cd.fermentation.response = merge(lsarp.cd.rapidaim.scores,
                                       lsarp.cd.rs.selected[,c("sample", "Selected")], by="sample")

lsarp.cd.fermentation.response = subset(lsarp.cd.fermentation.response,
                                        RS_Name == Selected)

lsarp.cd.fermentation.response$Response = ifelse(lsarp.cd.fermentation.response$delta.ph > -1.27, "Weak Response", "Strong Response")

lsarp.cd.fermentation.response = lsarp.cd.fermentation.response %>% dplyr::select(HM, Selected, delta.ph, med.ph, Response) %>% distinct()


# :: final data -----------------------------------------------------------

# add to main omics mapping file
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.phyloseq.butyrogens.i, by="standard.name", all.x=T)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.phyloseq.butyrogens.ii, by="standard.name", all.x=T)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.cd.asv.functionalredundancy, by="standard.name", all.x=T)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.cd.asv.mlp.preds, by="standard.name", all.x=T)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.cd.mgx.mlp.preds, by="standard.name", all.x=T)
metadata.lsarp.stool.asv = merge(metadata.lsarp.stool.asv,
                                 lsarp.cd.mpx.starch.mucin.mat, by="standard.name", all.x=T)

# in case you need to fix RapidAIM results: 
#metadata.lsarp.stool = metadata.lsarp.stool[,!colnames(metadata.lsarp.stool) %in%
#                                              c("Selected", "delta.ph", "med.ph", "Response")]
metadata.lsarp.stool = merge(metadata.lsarp.stool,
                             lsarp.cd.fermentation.response, by="HM", all.x=T)

## check dataset completion (n=192)
metadata.lsarp.stool.asv$standard.name %>% unique() %>% length()
# fcal
subset(metadata.lsarp.stool.asv, is.na(fcal)) # 5 samples missing fcal
# water
subset(metadata.lsarp.stool.asv, is.na(stool_water_perc)) # no samples missing water
# load
subset(metadata.lsarp.stool.asv, is.na(load.asv)) # 0 samples missing mlp.load with ASVs
subset(metadata.lsarp.stool.asv, is.na(load.mgx)) # many samples missing mlp.load with MGX

# rs.col = RS that was actually given to the participant (since some HM's used 16S + MPX, will not always equate to 16S-based approach)
metadata.lsarp.stool$rs.col = ifelse(!is.na(metadata.lsarp.stool$RS_Name), metadata.lsarp.stool$RS_Name, "grey")
metadata.lsarp.stool$rs.col = factor(metadata.lsarp.stool$rs.col, levels=c("grey", rs.names[rs.names %in% unique(metadata.lsarp.stool$RS_Name)]))
metadata.lsarp.stool.asv$rs.col = ifelse(!is.na(metadata.lsarp.stool.asv$RS_Name), metadata.lsarp.stool.asv$RS_Name, "grey")
metadata.lsarp.stool.asv$rs.col = factor(metadata.lsarp.stool.asv$rs.col, levels=c("grey", rs.names[rs.names %in% unique(metadata.lsarp.stool.asv$RS_Name)]))

# remove HM0933-STL-02 sample (duplicate; use repeat, STL-03)
metadata.lsarp.stool = subset(metadata.lsarp.stool, standard.name != "HM0933-STL-02")
metadata.lsarp.stool.asv = subset(metadata.lsarp.stool.asv, standard.name != "HM0933-STL-02")
# remove from mgx
lsarp.mgx.taxa = lsarp.mgx.taxa[rownames(lsarp.mgx.taxa)!= "HM0933-STL-02",]

# fix 0.5 m samples
metadata.lsarp.stool$phase = ifelse(metadata.lsarp.stool$lsarp.timing == "Month .5", "treatment", 
                                    metadata.lsarp.stool$phase)
# Add flag for flare (so we can remove from statistical analyses)
metadata.lsarp.stool$flare = ifelse(is.na(metadata.lsarp.stool$flare.day), "no_flare",
                                    ifelse(metadata.lsarp.stool$flare.day < metadata.lsarp.stool$lsarp.days, 
                                           "flare","no_flare"))
metadata.lsarp.stool.asv$flare = ifelse(is.na(metadata.lsarp.stool.asv$flare.day), "no_flare",
                                    ifelse(metadata.lsarp.stool.asv$flare.day < metadata.lsarp.stool.asv$lsarp.days, 
                                           "flare","no_flare"))
metadata.lsarp.stool$phase = ifelse(metadata.lsarp.stool$standard.name == "HM0860-STL-03", "treatment", 
                                    metadata.lsarp.stool$phase)
metadata.lsarp.stool$phase = ifelse(metadata.lsarp.stool$standard.name == "HM0860-STL-03", "treatment", 
                                    metadata.lsarp.stool$phase)
metadata.lsarp.stool.asv$phase = ifelse(metadata.lsarp.stool.asv$standard.name == "HM0860-STL-03", "treatment", 
                                        metadata.lsarp.stool.asv$phase)
metadata.lsarp.stool.asv$phase = ifelse(metadata.lsarp.stool.asv$standard.name == "HM0860-STL-03", "treatment", 
                                        metadata.lsarp.stool.asv$phase)

# remove duplicate
metadata.lsarp.stool = distinct(metadata.lsarp.stool)
metadata.lsarp.stool.asv = distinct(metadata.lsarp.stool.asv)

# save these, load into analysis script
save(
  # 16S data
  metadata.lsarp.stool, 
  metadata.lsarp.stool.asv, 
  lsarp.asv.data.median,
  lsarp.asv.data.glom,
  lsarp.phyloseq.butyrogens.i,
  lsarp.phyloseq.butyrogens.ii,
  lsarp.cd.asv.functionalredundancy,
  lsarp.cd.asv.mlp.preds,
  # MGX data
  lsarp.mgx.taxa,
  lsarp.cd.mgx.mlp.preds,
  # MPX data
  lsarp.cd.mpx.kegg.mat,
  lsarp.cd.mpx.cog.mat,
  lsarp.cd.mpx.cazy.mat,
  # MBX data
  lsarp.mbx.raw.mat,
  lsarp.mbx.annotated.mat,
  # RapidAIM pH
  lsarp.cd.rapidaim.scores,
  # RapidAIM selections
  lsarp.cd.rs.selections,
  # destination
  file = "./2025_12_02_lsarp_16s_data_meta.Renv")









