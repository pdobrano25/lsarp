### 2025_07_10  LSARP-UC Data processing + analysis

# Goal: all data processing and analysis for LSARP-UC


# :: save -----------------------------------------------------------------

v.date = "2026_01_08"
# save.image("./2025_09_12_lsarp_uc_analysis.Renv")

# load("./2025_09_12_lsarp_uc_analysis.Renv")

# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse")


# >> Load data ------------------------------------------------------------


# :: colors and vectors ---------------------------------------------------------------


rs.names.pbs = c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")
rs.names = c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors = data.frame(cols = c(gg_color_hue(5)[c(1,1,1,2,3,3,4,4,5)], "#000000"),
                         RS_Name = rs.names.pbs[c(2:10,1)])

# used for nearly all figures
rs.colors = c("Authentic" = "#F8766D", "BobsRedMill" = "#F8766D", "MSPrebiotic" = "#F8766D",
              "LetsDoOrganic" = "#A3A500", "HiMaize260" = "#00BF7D", "Novelose330" = "#00BF7D",
              "ActistarRT" = "#00B0F6", "FibersymRW" = "#00B0F6", "Versafibe1490" = "#E76BF3")
# 
# this function calculates within/between sample beta-diversity
# into a format compatible with standard variable analyses
beta.trajectory = function(data=lsarp.asv.bray){
  # melt
  data.subset = reshape2::melt(as.matrix(data))
  # note, I'm keeping the duplicate values (eg. X-Y and Y-X)
  # because attempting to remove them could introduce new errors (hard to pin down)
  # and the average will be the same if all values are duplicated
  
  # but we will remove self comparisons
  data.subset = subset(data.subset, Var1 != Var2)
  
  # now calculate between sample beta
  
  # unpack sample names
  data.subset = tidyr::separate(data.subset, col=Var1, into=c("HM.A", "STL.A", "Number.A"), sep="-", remove=F)
  data.subset = tidyr::separate(data.subset, col=Var2, into=c("HM.B", "STL.B", "Number.B"), sep="-", remove=F)
  # calculate between sample diversity
  between.beta = data.subset %>% 
    group_by(Var1) %>%
    mutate(between.beta = mean(value)) %>%
    dplyr::select(Var1, between.beta) %>% distinct()
  # calculate within patient diversity
  within.beta = data.subset %>% 
    subset(HM.A == HM.B) %>%
    group_by(Var1) %>%
    mutate(within.beta = mean(value)) %>%
    dplyr::select(Var1, within.beta) %>% distinct()
  # merge
  beta.data = merge(between.beta,
                    within.beta, by="Var1")
  colnames(beta.data)[1] = "standard.name"
  return(beta.data)
}


# function to make volcano plots have an intelligible y axis
neg_log10_trans <- scales::new_transform(
  name = "neglog10",
  transform = function(x) -log10(x),
  inverse = function(x) 10^(-x),
  format = function(x) format(x, scientific = FALSE)
)
# :: load patient list ----------------------------------------------------


# select patients (note; this version has lsarp.days added)
lsarp.uc = read.csv("./2025_06_03_lsarp_uc_unblinding_PD.csv")
lsarp.uc[,c("HM", "rs_selected", "flare_day", "lsarp_days", "group", "compliance")]
# subset to PP (>80% compliance)
lsarp.uc = subset(lsarp.uc, compliance == "good")

# simplify dataframe
lsarp.uc = lsarp.uc[,c("HM", "rs_selected", "flare_day", "lsarp_days", "group")]
lsarp.uc$flare = ifelse(is.na(lsarp.uc$flare_day), 0, 1)
lsarp.uc = lsarp.uc %>% arrange(HM)

# :: exclusions -----------------------------------------------------------

# patients excluded:
excluded.lsarp.uc = c("HM0908",  # poor compliance
                      "HM0923", # did not start
                      "HM0929", # did not start
                      "HM0990") # poor compliance
# already excluded from this list

# :: load RS selections ---------------------------------------------------

# uc
unique(lsarp.uc$HM)

lsarp.uc.scores.files = list.files("./lsarp_uc_rs_scores/")
lsarp.uc.scores.files = lsarp.uc.scores.files[grep("scores", lsarp.uc.scores.files)]

# note: I copied the scores files into a new folder
lsarp.uc.rs.selections = do.call(rbind, lapply(lsarp.uc.scores.files, function(x){
  data.loaded = readRDS(paste("./lsarp_uc_rs_scores/", x, sep=""))
  data.loaded$sample = substr(x, 1, 6)
  data.loaded
}))

# identify RS selected
lsarp.uc.rs.selected = lsarp.uc.rs.selections %>%
  group_by(sample) %>%
  subset(Z_score >=1) %>%
  group_by(sample) %>%
  filter(Distance==min(Distance))
colnames(lsarp.uc.rs.selected)[1] = "Selected"


# : CHECK -----------------------------------------------------------------


# check this exception:
#lsarp.cd.rs.selections.exceptions.2 = lsarp.cd.rs.selections %>%
#  subset(sample == "HM0869" & RS == "Versafibe1490")

# add HM1033, which had Fibersym pull all Z-scores < 1
lsarp.uc.rs.selected.1033 = lsarp.uc.rs.selections %>%
  subset(sample == "HM1033") %>%
  subset(Z_score == max(Z_score))
colnames(lsarp.uc.rs.selected.1033)[1] = "Selected"

lsarp.uc.rs.selected = rbind(lsarp.uc.rs.selected,
                             lsarp.uc.rs.selected.1033) %>% data.frame()

# stitch back on
lsarp.uc.rs.selections = merge(lsarp.uc.rs.selections,
                               lsarp.uc.rs.selected[,c("Selected", "sample")], by="sample")

# final clean

lsarp.uc.rs.selections$RS_Name = factor(lsarp.uc.rs.selections$RS, levels=rs.names)

lsarp.uc.rs.selections$HM = lsarp.uc.rs.selections$sample

lsarp.uc.rs.selections = merge(lsarp.uc.rs.selections,
                               distinct(lsarp.uc[,c("HM", "group")]), by="HM")
lsarp.uc.rs.selections$Selected = factor(lsarp.uc.rs.selections$Selected, levels=rs.names)

# good
lsarp.uc.rs.selections

unique(lsarp.uc.rs.selections$HM) # n=15

# :: load REDCap data ----------------------------------------------------------

lsarp.uc.redcap =  read.csv("~/Downloads/redCap_sampleData_uniques_240724.tsv", sep="\t")

# subset to LSARP-UC patients
lsarp.uc.redcap = subset(lsarp.uc.redcap, study_id %in% paste(lsarp.uc$HM, ".00", sep=""))

# reduce to important columns
lsarp.uc.redcap = lsarp.uc.redcap[,c("study_id", "fecalcal_res", "standard.name","stool_date_rec_v2",
                                     "trial.stool.timing",
                                     "rs_start_date","product_end", "rs.trial.withdraw.date")] %>% 
  subset(grepl("STL", standard.name)) %>% distinct()

# format dates
lsarp.uc.redcap$stool_date_rec_v2 = as.Date(lsarp.uc.redcap$stool_date_rec_v2)
lsarp.uc.redcap$rs_start_date = as.Date(lsarp.uc.redcap$rs_start_date)
lsarp.uc.redcap$product_end = as.Date(lsarp.uc.redcap$product_end)
lsarp.uc.redcap$rs.trial.withdraw.date = as.Date(lsarp.uc.redcap$rs.trial.withdraw.date)
# calculate durations
lsarp.uc.redcap$lsarp.days = lsarp.uc.redcap$stool_date_rec_v2 - lsarp.uc.redcap$rs_start_date
lsarp.uc.redcap$lsarp.off = lsarp.uc.redcap$stool_date_rec_v2 - lsarp.uc.redcap$product_end
# add HM
lsarp.uc.redcap$HM = substr(lsarp.uc.redcap$study_id, 1, 6)

# eliminate non-LSARP samples
lsarp.uc.redcap = subset(lsarp.uc.redcap,
                         grepl(paste(c("baseline", "Month 1", "Month 2", "Month 3",
                                       "Month 4", "Month 5", "Month 6", "Month 8",
                                       "Month 10", "Month 12"), collapse="|"), trial.stool.timing))
# eliminate samples collected after withdrawing
# add HM1053 withdrawal date: 07-05-2024
lsarp.uc.redcap$rs.trial.withdraw.date = ifelse(lsarp.uc.redcap$HM == "HM1053", "2024-05-07", as.character(lsarp.uc.redcap$rs.trial.withdraw.date))
lsarp.uc.redcap$rs.trial.withdraw.date = as.Date(lsarp.uc.redcap$rs.trial.withdraw.date, "%Y-%m-%d")
lsarp.uc.redcap = subset(lsarp.uc.redcap, stool_date_rec_v2 <= rs.trial.withdraw.date)

# assign phase
lsarp.uc.redcap$phase = ifelse(lsarp.uc.redcap$lsarp.days < 0, "preRS",
                               ifelse(lsarp.uc.redcap$lsarp.off <= 3 & lsarp.uc.redcap$lsarp.days >= 0, "onRS", "postRS"))
lsarp.uc.redcap$phase = factor(lsarp.uc.redcap$phase, levels=c("preRS", "onRS", "postRS"))
# confirmed with Dave's data; same result
# add flare day
lsarp.uc.redcap = merge(lsarp.uc.redcap,
                        lsarp.uc[,c("HM","group", "flare_day")], by="HM")
# fix group name
lsarp.uc.redcap$group = ifelse(lsarp.uc.redcap$group == "Plac", "Placebo", "RS")

# plot

# ``` plot: stool timings -------------------------------------------------

# visualize stools + timings

lsarp.uc.stool.plot = ggplot(lsarp.uc.redcap %>% subset(!HM %in% excluded.lsarp.uc) %>%
                               # remove post-flare
                               subset(lsarp.days <= flare_day | is.na(flare_day)),
                             aes(x=lsarp.days, y=reorder(HM, flare_day)))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase%in% c("preRS", "onRS"))$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  # background point
  geom_point(fill="white", shape=21, size=3)+
  # actual point
  geom_point(aes(fill=phase), 
                 #alpha = ifelse(lsarp.days <= flare_day, "1", "0")),
             color="white",
             shape=21, size=5)+
  # add flare date
  geom_point(data=lsarp.uc.redcap[,c("flare_day", "lsarp.days","HM", "group")] %>% 
               subset(!HM %in% excluded.lsarp.uc) %>% distinct(),
             aes(x=flare_day, y=reorder(HM, (lsarp.days))), shape=4, size=5)+
  scale_fill_manual(values=c("black",2, "grey"))+
  scale_alpha_manual(values=c(0.25,1))+
  #geom_text(aes(label=missing.16s) ,color="white", nudge_y=-0.1, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-2, nchar(standard.name))), vjust=-1, size=2)+
  theme_classic()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group,nrow=1, scales="free")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.major.y=element_line(color="grey", linewidth=0.2),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.uc.stool.plot


# >> Calculate Fibre + RS ------------------------------------------------------

# 2025_07_08  Goal: add supplemented RS to RS intake

# load REDCap for RS dose
lsarp.uc.redcap.data.all.loaded <- read.csv("~/Downloads/redCap_sampleData_uniques_240724.tsv", sep="\t")
# unpack monash date
lsarp.uc.redcap.data.loaded = lsarp.uc.redcap.data.all.loaded %>%
  # subset to HM's in the study
  subset(study_id %in% lsarp.uc.redcap$study_id) %>%
  # clean up ffq.date name
  mutate(ffq.date = as.Date(substr(monash.ffq.date, 
                           nchar(monash.ffq.date)-6,
                           nchar(monash.ffq.date)), "%d%b%y")) %>% 
  # remove NA entries
  subset(!is.na(ffq.date))

# manually add HM1053 withdrawal date: 07-05-2024
lsarp.uc.redcap.data.loaded$rs.trial.withdraw.date = ifelse(lsarp.uc.redcap.data.loaded$study_id == "HM1053.00", "2024-05-07", as.character(lsarp.uc.redcap.data.loaded$rs.trial.withdraw.date))
lsarp.uc.redcap.data.loaded$rs.trial.withdraw.date = as.Date(lsarp.uc.redcap.data.loaded$rs.trial.withdraw.date, "%Y-%m-%d")

# since these are new onset, we can simply subset to before withdrawn from study; diagnosis and any others will be ~within study period
lsarp.uc.redcap.data.loaded = subset(lsarp.uc.redcap.data.loaded, ffq.date <= rs.trial.withdraw.date)

# clean up to extract out necessary diet data
lsarp.uc.redcap.data.loaded = subset(lsarp.uc.redcap.data.loaded, !is.na(dose))[,c("study_id", "fecalcal_res","ffq.date",
                                                                             "energy", "resistant_starch","dietary_fibre", "dose", "dose_2", "dose_3",
                                                                             "percent_product", "percent_product_2", "percent_product_3",
                                                                             "lot_number", "lot_number_2", "lot_number_3", "m2")] %>% distinct()
unique(lsarp.uc.redcap.data.loaded$study_id) # all 15

## calculate adj.fibre and adj.rs
lsarp.uc.redcap.data.loaded$fibre_adj = lsarp.uc.redcap.data.loaded$dietary_fibre / lsarp.uc.redcap.data.loaded$energy * 4.184 * 1000
lsarp.uc.redcap.data.loaded$rs_adj = lsarp.uc.redcap.data.loaded$resistant_starch / lsarp.uc.redcap.data.loaded$energy * 4.184 * 1000

# take average of lots
lsarp.uc.redcap.data.loaded = lsarp.uc.redcap.data.loaded %>%
  group_by(study_id, ffq.date) %>%
  mutate(rs_dose = mean(na.omit(c(dose,dose_2,dose_3)))) %>%
  mutate(rs_perc = mean(na.omit(c(percent_product,percent_product_2,percent_product_3)))) %>%
  mutate(rs_intake = rs_dose * rs_perc/100) %>% data.frame()

lsarp.uc.redcap.data.loaded$HM = substr(lsarp.uc.redcap.data.loaded$study_id, 1, 6)
lsarp.uc.redcap.data.loaded = merge(lsarp.uc.redcap.data.loaded,
                                    lsarp.uc[,c("HM", "group")], by="HM")

# add calories from starch (more complex than assumed)
lsarp.uc.redcap.data.loaded$kcal = lsarp.uc.redcap.data.loaded$energy / 4.184
lsarp.uc.redcap.data.loaded$rs.kcal = 
  ifelse(lsarp.uc.redcap.data.loaded$group == "RS",
         # add calories from RS fraction
         (lsarp.uc.redcap.data.loaded$rs_dose * lsarp.uc.redcap.data.loaded$rs_perc/100 * 2.7) + # cite: Miketinas, 2020 https://pubmed.ncbi.nlm.nih.gov/32840627/
           # plus calories from non-RS fraction
           (lsarp.uc.redcap.data.loaded$rs_dose * (1-lsarp.uc.redcap.data.loaded$rs_perc/100) * 4), # 4 kcal/g carbohydrates
         ## FOR PLACEBO, back-calculate Amioca dose
         # note: RS dose = 7.5 g RS / m2 (BSA)
         # therefore: 7.5 / RS Dose = BSA
         # 4 / BSA = Amioca dose
         # Amioca dose * 4 = Amioca calories
         4 * 4 * lsarp.uc.redcap.data.loaded$m2) # 4 kcal/g carbohydrates for amioca
lsarp.uc.redcap.data.loaded = merge(lsarp.uc.redcap.data.loaded,
                                    lsarp.uc.redcap[,c("HM", "rs_start_date", "product_end")]%>%distinct(),
                                 by="HM")

# calculate days since starting RS
lsarp.uc.redcap.data.loaded$lsarp.uc.days = as.Date(lsarp.uc.redcap.data.loaded$ffq.date) - as.Date(lsarp.uc.redcap.data.loaded$rs_start_date)

# check energy values
hist(lsarp.uc.redcap.data.loaded$kcal)

# replace RS value with 0 for days outside of treatment period
lsarp.uc.redcap.data.loaded$rs_intake = ifelse(lsarp.uc.redcap.data.loaded$group == "RS" & lsarp.uc.redcap.data.loaded$ffq.date <= lsarp.uc.redcap.data.loaded$product_end &
                                              lsarp.uc.redcap.data.loaded$group == "RS" & lsarp.uc.redcap.data.loaded$ffq.date >= lsarp.uc.redcap.data.loaded$rs_start_date,
                                            lsarp.uc.redcap.data.loaded$rs_intake, 0)

# add dietary energy + dietary supplement
lsarp.uc.redcap.data.loaded$energy = lsarp.uc.redcap.data.loaded$kcal + lsarp.uc.redcap.data.loaded$rs.kcal
# total rs from diet and supplement
lsarp.uc.redcap.data.loaded$total_rs = lsarp.uc.redcap.data.loaded$resistant_starch + lsarp.uc.redcap.data.loaded$rs_intake
# re-normalize to energy
lsarp.uc.redcap.data.loaded$total_rs_adj = lsarp.uc.redcap.data.loaded$total_rs / lsarp.uc.redcap.data.loaded$energy * 1000


# :: ``` longitudinal fiber plot -----------------------------------------------------------

metadata.lsarp.uc.fiber.plot = ggplot(lsarp.uc.redcap.data.loaded %>% subset(!is.na(energy)),
                                   aes(x=as.numeric(lsarp.uc.days), 
                                       y=fibre_adj))+
  annotate("rect", xmin=0, xmax=max(na.omit(subset(lsarp.uc.redcap, phase == "onRS")$lsarp.days)),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="black",  se=F, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.uc.days > 0 & lsarp.uc.days < max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.uc.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Fiber Intake (g per day / 1000 kcal)")
metadata.lsarp.uc.fiber.plot


# :: ``` average fiber plot -----------------------------------------------------------

metadata.lsarp.uc.average.fiber.plot = ggplot()+
  geom_segment(data=lsarp.uc.redcap.data.loaded %>%
                 group_by(HM) %>%
                 mutate(min.fiber = min(na.omit(fibre_adj)),
                        max.fiber = max(na.omit(fibre_adj)),
                        mean.fiber = mean(na.omit(fibre_adj))), 
               aes(x=min.fiber, y=reorder(HM, mean.fiber),
                   xend=max.fiber, yend=reorder(HM, mean.fiber)),
               linewidth=0.3)+
  geom_point(data=lsarp.uc.redcap.data.loaded %>% group_by(HM) %>% 
               mutate(mean.fib.adj = mean(na.omit(fibre_adj))) %>%
               dplyr::select(HM, mean.fib.adj) %>% distinct(),
             aes(x=mean.fib.adj,
                 y=reorder(HM, mean.fib.adj), 
                 fill=scale(mean.fib.adj)), shape=21, size=3)+
  coord_flip()+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        axis.text.x=element_blank())+
  labs(x="Average Fiber Intake (g per day / 1000 kcal)", y="Participant", fill = "")
metadata.lsarp.uc.average.fiber.plot


# >> Fiber vs RS ----------------------------------------------------------

lsarp.uc.nutrient.data.fiber.rs.plot = ggplot(lsarp.uc.redcap.data.loaded,
                                           aes(x=fibre_adj, rs_adj))+
  geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(shape=21, color="white", fill="black", size=3)+
  geom_smooth(method="lm", color="black", se=F)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+
  labs(x="Fiber Intake g / 1000 kcal per day",
       y="RS Intake g / 1000 kcal per day")
lsarp.uc.nutrient.data.fiber.rs.plot


# :: ``` longitudinal rs plot ------------------------------------------------------------

metadata.lsarp.uc.rs.plot = ggplot(lsarp.uc.redcap.data.loaded,
                                aes(x=lsarp.uc.days, 
                                    y=resistant_starch))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="white",  se=T, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.uc.days > 0 & lsarp.uc.days < max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.uc.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Dietary Resistant Starch Intake (g per day)")
metadata.lsarp.uc.rs.plot


metadata.lsarp.uc.rs.supp.plot = ggplot(lsarp.uc.redcap.data.loaded,
                                     aes(x=lsarp.uc.days, 
                                         y=total_rs))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_smooth(color="black",  se=F, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.uc.days > 0 & lsarp.uc.days < max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.uc.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2, "grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none",
                        strip.background = element_rect(color="black"),
                        strip.text=element_text(size=10))+
  facet_wrap(~group)+
  labs(x="Days since starting Product", y="Total Resistant Starch Intake (g per day)")
metadata.lsarp.uc.rs.supp.plot

metadata.lsarp.uc.rs.plot+metadata.lsarp.uc.rs.supp.plot


# :: impute average fiber -------------------------------------------------

# calculate average fiber
lsarp.uc.redcap.data.loaded = lsarp.uc.redcap.data.loaded %>%
  group_by(HM) %>%
  mutate(ave.fiber = mean(na.omit(fibre_adj)))
# double check

# use average fiber intake to impute NA values
lsarp.uc.redcap.data.loaded = lsarp.uc.redcap.data.loaded %>%
  group_by(HM) %>%
  mutate(adj.fiber = ifelse(is.na(fibre_adj), ave.fiber, fibre_adj))
# double check



ggplot(lsarp.uc.redcap.data.loaded,
       aes(x=lsarp.uc.days, 
           y=adj.fiber))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white",  se=T, linewidth=0.5)+
  geom_point(aes(fill=ifelse(lsarp.uc.days > 0 & lsarp.uc.days < max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days), "1", "2"), 
                 shape=ifelse(lsarp.uc.days <= 0, "1", "2")), color="white", size=3)+
  scale_fill_manual(values=c(2,"grey","grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days since starting RS", 
       y="Fiber Intake (g per day)")


# add average fiber intake to main mapping
metadata.lsarp.uc.stool = merge(lsarp.uc.redcap, 
                             lsarp.uc.redcap.data.loaded[,c("HM", "ave.fiber")] %>% distinct(), 
                             by="HM", all.x = T)

lsarp.uc.fiber.mean.pval = 
  wilcox.test(distinct(subset(metadata.lsarp.uc.stool, group == "RS")[,c("HM", "ave.fiber")])$ave.fiber,
              distinct(subset(metadata.lsarp.uc.stool, group == "Placebo")[,c("HM", "ave.fiber")])$ave.fiber)$p.value

lsarp.uc.fiber.comparison.plot = ggplot(metadata.lsarp.uc.stool[,c("HM", "ave.fiber", "group")] %>% distinct(),
                                     aes(x=group, y=ave.fiber))+
  geom_boxplot(width=0.5, alpha=0.5)+
  #geom_boxplot(width=0.5, alpha=0.5, outlier.shape = NA, coef = 0)+
  ggbeeswarm::geom_beeswarm(aes(fill=group), shape=21, color="white", size=3)+
  ylim(c(5,20))+
  annotate(geom="text", y=Inf, x=1.5, vjust=1.5, size=4,
           label=paste("p: ", round(lsarp.uc.fiber.mean.pval, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Average Fiber Intake")+
  labs(x="", y="Average Fiber Intake g / 1000 kcal per day")
lsarp.uc.fiber.comparison.plot


# Brief methods write up:
# Fiber intakes were calculated from Monash FFQ and adjusted for calculated energy intake. 
# Values were calculated per stool. Missing values were imputed using the average of 3+ FFQs.

# Resistant starch doses & percentages were averaged per RS phase
# Total energy was adjusted for the RS (2.7 kcal/g) and digestible starch (4 kcal/g)

###
# :: load 16S data --------------------------------------------------------

lsarp.uc.asv = readRDS("./2024_11_20_rs_trial_16s_data_rarefied_gg2.Rds")

# reduce to LSARP-UC list
lsarp.uc.asv = phyloseq::subset_samples(lsarp.uc.asv, standard.name %in% lsarp.uc.redcap$standard.name)

# eliminate ASVs not present
lsarp.uc.asv <- phyloseq::prune_taxa(phyloseq::taxa_sums(lsarp.uc.asv) > 0, lsarp.uc.asv)

# visualize completion
phyloseq::sample_data(lsarp.uc.asv)$standard.name

# visualize stools + timings (old version)
ggplot(lsarp.uc.redcap,
       aes(x=lsarp.days, y=reorder(HM, lsarp.off)))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase=="onRS")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(color="black", size=3)+
  geom_point(color="white", size=2)+
  geom_point(aes(color=phase))+
  scale_color_manual(values=c("black",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  geom_point(data=lsarp.uc.redcap[,c("flare_day", "HM")] %>% distinct(), aes(x=flare_day, y=HM), shape=4, size=5)+
  geom_text(data=lsarp.uc.redcap[,c("lsarp.days","standard.name", "HM")] %>% mutate(asv = ifelse(!standard.name %in% phyloseq::sample_data(lsarp.uc.asv)$standard.name, "*", NA)),
            aes(x=lsarp.days, y=HM, label=asv), vjust=0.76, size=6)+
  theme_bw()+
  labs(x="Day since starting product", y="", color="Phase", alpha="")

# convert to median of extraction replicates
lsarp.uc.asv.median = speedyseq::psmelt(lsarp.uc.asv)

lsarp.uc.asv.median = reshape2::acast(lsarp.uc.asv.median,
                                      standard.name ~ OTU, value.var="Abundance", fun.aggregate=median)

# replace taxa names
lsarp.uc.asv.median.tax = lsarp.uc.asv.median
lsarp.uc.asv.tax = phyloseq::tax_table(lsarp.uc.asv) %>% data.frame()
colnames(lsarp.uc.asv.median.tax) = lsarp.uc.asv.tax$Taxa[match(colnames(lsarp.uc.asv.median.tax), rownames(lsarp.uc.asv.tax))]

# and save glom'ed
lsarp.uc.asv.median.glom = reshape2::acast(speedyseq::psmelt(lsarp.uc.asv),
                                      standard.name ~ LCA, value.var="Abundance", fun.aggregate=median)

# save these
lsarp.uc.asv.median.tax
lsarp.uc.asv.median.glom

# :: :: alpha -------------------------------------------------------------

## alpha diversity
lsarp.uc.asv.median.richness = vegan::specnumber(lsarp.uc.asv.median)
lsarp.uc.asv.median.shannon = vegan::diversity(lsarp.uc.asv.median, index="shannon")

lsarp.uc.asv.median.alpha = data.frame(
  standard.name = rownames(lsarp.uc.asv.median),
  richness = lsarp.uc.asv.median.richness,
  shannon = lsarp.uc.asv.median.shannon
)

# :: :: but i -------------------------------------------------------------

## butyrogens I (classic definition)

# load taxa
amplicon.data.gg138.tax.df = readRDS("./2024_11_20_rs_trial_16s_data_tax_table_gg2.Rds") # note: misnamed "gg2", should be "gg138"

# subset to butyrogens
lsarp.uc.asv.gg138.tax.df = subset(data.frame(amplicon.data.gg138.tax.df), Family=="f__Lachnospiraceae" | Genus=="g__Blautia" | Genus=="g__Roseburia" | Genus=="g__Eubacterium" | Genus=="g__Ruminococcus" | Genus=="g__Clostridium" | Genus=="g__Faecalibacterium") # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
# apply this filter to phyloseq data
lsarp.uc.butyrogens.i = speedyseq::psmelt(lsarp.uc.asv)
lsarp.uc.butyrogens.i = lsarp.uc.butyrogens.i %>%
  subset(OTU %in% lsarp.uc.asv.gg138.tax.df$ASV) %>%
  group_by(Sample) %>%
  mutate(but = sum(Abundance/50000)) %>%
  group_by(standard.name) %>%
  # take median of extraction replicates
  mutate(but = median(but)) %>%
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but) %>% distinct() %>% data.frame()
lsarp.uc.butyrogens.i %>% arrange(standard.name)

lsarp.uc.butyrogens.i = lsarp.uc.butyrogens.i[,c("standard.name", "but")]
colnames(lsarp.uc.butyrogens.i) = c("standard.name", "but.i")

# :: :: picrust -------------------------------------------------------------

## butyrogens II (Vital definition)
# prepare seq.table for picrust2
fasta_seqs <- Biostrings::DNAStringSet(data.frame(phyloseq::refseq(lsarp.uc.asv))[,1])
names(fasta_seqs) <- rownames(data.frame(phyloseq::refseq(lsarp.uc.asv)))  # assign ASV IDs as sequence names
# prepare asv.table for picrust2
lsarp.uc.picrust.abuntable = data.frame(phyloseq::otu_table(lsarp.uc.asv))
# export
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/lsarp.uc.picrust2.seqtable.fasta")
write.table(t(lsarp.uc.picrust.abuntable), "oars_picrust2/lsarp.uc.picrust2.abuntable.tsv", 
            sep="\t", quote = F, col.names = NA)

terminal=F
if(terminal==T){
# STEP 1: run through default picrust2 (to ensure picrust2 works)
cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2
conda activate oars_picrust2
# In R, may need to install Rcpp, jsonlite, lattice, Matrix, RSpectra, castor
rm -r lsarp_uc_picrust2_out

picrust2_pipeline.py \
-s lsarp.uc.picrust2.seqtable.fasta \
-i lsarp.uc.picrust2.abuntable.tsv \
-o lsarp_uc_picrust2_out \
-p 1

# optionally perform stratification to source functions to taxa
# use EC's (not preferable; 32 million rows)
metagenome_pipeline.py -i lsarp.uc.picrust2.abuntable.tsv -m lsarp_uc_picrust2_out/marker_predicted_and_nsti.tsv.gz -f lsarp_uc_picrust2_out/EC_predicted.tsv.gz \
-o lsarp_uc_picrust2_out/EC_metagenome_out --strat_out


# ::  :: but ii-------------------------------------------------------

cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2/predict_SCFA_producers
rm -r placement_working
place_seqs.py -s ../lsarp.uc.picrust2.seqtable.fasta -o placed_seqs.tre -p 1 --intermediate placement_working --ref_dir SCFA
# 18 sequences failed to align
hsp.py -t placed_seqs.tre --observed_trait_table SCFA/SCFA_pathwaydata.txt -o lsarp_uc_SCFA_predicted.tsv -p 1 -m emp_prob -n

# save files to new folder
cp ./lsarp_uc_SCFA_predicted.tsv ../picrust2_saved/lsarp_uc_SCFA_predicted.tsv
cp ../lsarp_uc_picrust2_out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz ../picrust2_saved/lsarp_uc_pred_metagenome_contrib.tsv

}

# import
lsarp.uc.butyrogens.ii = read.csv("./oars_picrust2/picrust2_saved/lsarp_uc_SCFA_predicted.tsv", sep="\t")
# subset to acetylcoa+but or buk alone
lsarp.uc.butyrogens.ii.acetylcoa = subset(lsarp.uc.butyrogens.ii, acetylcoa == 1)
lsarp.uc.butyrogens.ii.acetylcoa = subset(lsarp.uc.butyrogens.ii.acetylcoa, but == 1 | buk == 1)

# return to phyloseq object
lsarp.uc.butyrogens.ii <- phyloseq::otu_table(lsarp.uc.asv)[,lsarp.uc.butyrogens.ii.acetylcoa$sequence]
lsarp.uc.butyrogens.ii <- phyloseq::merge_phyloseq(lsarp.uc.butyrogens.ii, 
                                                     phyloseq::tax_table(lsarp.uc.asv), phyloseq::sample_data(lsarp.uc.asv))
# format to dataframe
lsarp.uc.butyrogens.ii = speedyseq::psmelt(lsarp.uc.butyrogens.ii)
lsarp.uc.butyrogens.ii$Abundance = lsarp.uc.butyrogens.ii$Abundance / 50000
# visualize
lsarp.uc.butyrogens.ii %>%
  group_by(LCA) %>%
  mutate(sum.abun = sum(Abundance)) %>% dplyr::select(LCA, sum.abun) %>% distinct() %>% data.frame() %>%
  slice_max(n=20, sum.abun) %>%
  ggplot(aes(x=reorder(LCA, sum.abun), y=sum.abun)) +
  coord_flip()+scale_y_log10()+
  geom_point(shape=21, aes(fill=sum.abun), size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none")+labs(x="", y="Total Abundance (log10)")
# top are Faecalibacterium, Gemmiger, Agathobacter, Anaerostipes, Anaerobutyricum, Roseburia

lsarp.uc.butyrogens.ii = lsarp.uc.butyrogens.ii %>%
  group_by(Sample) %>%
  mutate(but = sum(Abundance)) %>%
  group_by(standard.name) %>%
  mutate(but.ii = median(but)) %>%
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but.ii) %>% distinct() %>% data.frame()
# take median
lsarp.uc.butyrogens.ii %>% arrange(standard.name)

lsarp.uc.butyrogens.ii = lsarp.uc.butyrogens.ii[,c("standard.name", "but.ii")]
colnames(lsarp.uc.butyrogens.ii) = c("standard.name", "but.ii")

# good

# :: :: Functional Diversity ---------------------------------------------

run.fd = T
if(run.fd == T){
  # using PICRUSt2
  lsarp.uc.picrust <- read.csv("./oars_picrust2/picrust2_saved/lsarp_uc_pred_metagenome_contrib.tsv", sep="\t")
  
  # calculate redundancy per function
  picrust_data_redun = lsarp.uc.picrust %>%
    group_by(sample, function.) %>%
    mutate(ntaxa = length(unique(taxon))) %>%
    dplyr::select(sample, function., ntaxa) %>% distinct()
  
  # average functional redundancy (across functions) per sample
  picrust_data_redun_mean <- picrust_data_redun %>%
    group_by(sample) %>%
    # count n unique functions per sample
    mutate(nfunctions = length(unique(function.))) %>%
    mutate(meantaxa = mean(ntaxa)) %>%
    dplyr::select(sample, meantaxa, nfunctions) %>% distinct() %>%
    # functional diversity (richness)
    mutate(fd = nfunctions) %>%
    # functional redundancy
    mutate(fr = meantaxa) %>%
    dplyr::select(sample, fd, fr)
  # on average, each function has ~meantaxa contributing to it
  
  colnames(picrust_data_redun_mean)[1] = "dada2.sampleNames"
  
  # link to STL and take median
  lsarp.uc.asv.functionalredundancy = merge(picrust_data_redun_mean,
                                        data.frame(phyloseq::sample_data(lsarp.uc.asv))[,c("standard.name", "dada2.sampleNames")],
                                        by="dada2.sampleNames")
  lsarp.uc.asv.functionalredundancy = lsarp.uc.asv.functionalredundancy %>%
    group_by(standard.name) %>%
    mutate(fd = median(fd)) %>%
    mutate(fr = median(fr)) %>%
    dplyr::select(standard.name, fd, fr) %>% distinct() %>% data.frame()
  
  saveRDS(lsarp.uc.asv.functionalredundancy, "./2025_07_10_lsarp_uc_functional_redundancy.Rds")
  
}

lsarp.uc.asv.functionalredundancy = readRDS("./2025_07_10_lsarp_uc_functional_redundancy.Rds")

lsarp.uc.asv.functionalredundancy %>% arrange(standard.name)


# :: :: MLP ---------------------------------------------------------------


# load RDP data
amplicon.data.rdp.tax.df = readRDS("./2024_11_20_rs_trial_16s_data_tax_table_rdp.Rds")
# use RDP tax to apply MLP (check internal validation performance)
amplicon.data.rdp.tax.df$OTU = amplicon.data.rdp.tax.df$ASV

# Apply (but also, validate) Bork's Microbial Load Predictor
# collapse to median
lsarp.uc.asv.data.for.mlp = speedyseq::psmelt(phyloseq::subset_samples(lsarp.uc.asv, standard.name %in% metadata.lsarp.uc.stool$standard.name))
lsarp.uc.asv.data.for.mlp = lsarp.uc.asv.data.for.mlp[,c("OTU", "Abundance", "standard.name")]
lsarp.uc.asv.data.for.mlp$Abundance = lsarp.uc.asv.data.for.mlp$Abundance / 50000
# take median % of ASVs
lsarp.uc.asv.data.for.mlp = lsarp.uc.asv.data.for.mlp %>%
  group_by(standard.name, OTU) %>%
  mutate(Abundance = median(Abundance)) %>% distinct()
lsarp.uc.asv.data.for.mlp = merge(lsarp.uc.asv.data.for.mlp[,c("OTU", "standard.name", "Abundance")],
                                  amplicon.data.rdp.tax.df, by="OTU") #

# format lsarp.uc taxa (e.g. all genera, or uc_o/c/f) to match MLP
lsarp.uc.asv.data.for.mlp = lsarp.uc.asv.data.for.mlp %>%
  mutate(RDP_taxa = ifelse(!is.na(Genus), paste("g_", Genus,sep=""),
                           ifelse(!is.na(Family), paste("uc_f_", Family, sep=""),
                                  ifelse(!is.na(Order), paste("uc_o_", Order, sep=""),
                                         ifelse(!is.na(Class), paste("uc_c_", Class, sep=""),
                                                ifelse(!is.na(Phylum), paste("uc_p_", Phylum, sep=""),
                                                       "Bacteria"))))))

# make matrix; sum up glom'ed taxa using RDP annotation (range = 0-67.2%)
lsarp.uc.asv.data.for.mlp = reshape2::acast(lsarp.uc.asv.data.for.mlp,
                                            standard.name ~ RDP_taxa, value.var="Abundance",
                                            fun.aggregate=sum) %>% as.data.frame()
# Note: I'm doing this in a different order than other 16S processing
# e.g. here: median --> sum (if median is 0, lose that taxa)
# elsewhere: sum --> median (if median is 0, sum may rescue)

# add Shannon
lsarp.uc.asv.data.for.mlp.shannon = vegan::diversity(lsarp.uc.asv.data.for.mlp)
# pseudocount
asv.mlp.min = min(lsarp.uc.asv.data.for.mlp[lsarp.uc.asv.data.for.mlp!=0])/2

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
# add Shannon
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
# briefly return to lsarp.uc data
lsarp.uc.asv.data.for.mlp = log2(lsarp.uc.asv.data.for.mlp+asv.mlp.min.to.use)
lsarp.uc.asv.data.for.mlp$Shannon = lsarp.uc.asv.data.for.mlp.shannon
# --
# now back to MLP data
mlp.rdp.retrain.data = log2(mlp.rdp.retrain.data+asv.mlp.min.to.use)
mlp.rdp.retrain.data$Shannon = mlp.rdp.retrain.data.shannon
# 4. add load data
mlp.rdp.retrain.data = merge(mlp.rdp.retrain.data,
                             mlp.rdp.retrain.data.load[,c("Cell_count_per_gram", "X")], by="row.names")
mlp.rdp.retrain.data$X.x = NULL
mlp.rdp.retrain.data$X.y = NULL
rownames(mlp.rdp.retrain.data) = mlp.rdp.retrain.data$Row.names
mlp.rdp.retrain.data$Row.names = NULL
mlp.rdp.retrain.data.load = log10(mlp.rdp.retrain.data$Cell_count_per_gram) # Note: Log10
mlp.rdp.retrain.data$Cell_count_per_gram = NULL
# 5. take overlapping features
mlp.rdp.retrain.data = mlp.rdp.retrain.data[,colnames(mlp.rdp.retrain.data) %in% colnames(lsarp.uc.asv.data.for.mlp)]
mlp.rdp.retrain.data$load = mlp.rdp.retrain.data.load
dim(mlp.rdp.retrain.data) # 108 taxa
mlp.rdp.retrain.data = subset(mlp.rdp.retrain.data, !is.na(load))

# 6. build/validate model

# note: XGBoost predictions are bimodal
# I will need to re-validate their results using their validation data.
# Perhaps try different models (glmnet seems to work better, for instance)

t1 = Sys.time()
set.seed(25)
mlp.rdp.model.xgb = caret::train(load ~.,
                                 mlp.rdp.retrain.data,
                                 method = "xgbTree",
                                 trControl = caret::trainControl(method = "cv", number = 10,
                                                                 savePredictions=TRUE),
                                 verbose = FALSE
)
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
# Pearson Cor = 0.79

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
mlp.rdp.original.plot = ggplot(mlp.rdp.original.preds,
                               aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Original")+labs(x="Prediction", y="Observed")
mlp.rdp.original.plot
# Cor = 0.79
# Retrained is the same! Let's use it


# apply to lsarp.uc
lsarp.uc.asv.mlp.preds.xgb = predict(mlp.rdp.model.xgb, lsarp.uc.asv.data.for.mlp)
lsarp.uc.asv.mlp.preds.xgb = data.frame(standard.name = gsub("\\.", "-", rownames(lsarp.uc.asv.data.for.mlp)),
                                        load = lsarp.uc.asv.mlp.preds.xgb)

lsarp.uc.asv.mlp.preds.xgb.plot = ggplot(lsarp.uc.asv.mlp.preds.xgb,
                                         aes(x=load, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained XGBoost")+
  labs(x="Predicted Load", y="Density")
lsarp.uc.asv.mlp.preds.xgb.plot

# :: :: MLP Validation ---------------------------------------------------

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

# apply to lsarp.uc
lsarp.uc.asv.mlp.preds = predict(mlp.rdp.model.final, lsarp.uc.asv.data.for.mlp)
lsarp.uc.asv.mlp.preds = data.frame(standard.name = gsub("\\.", "-", rownames(lsarp.uc.asv.data.for.mlp)),
                                    load.asv = lsarp.uc.asv.mlp.preds)

lsarp.uc.asv.mlp.preds.rf.plot = ggplot(lsarp.uc.asv.mlp.preds,
                                        aes(x=load.asv, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained Ranger")+
  labs(x="Predicted Load", y="Density")
lsarp.uc.asv.mlp.preds.rf.plot
# Much smoother than xgBoost
# Good


# plus
(mlp.rdp.retrain.models.plots/
    (lsarp.uc.asv.mlp.preds.xgb.plot|
       lsarp.uc.asv.mlp.preds.rf.plot))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_mlp_plots.pdf",
         width=4, height=6,device = cairo_pdf)

# Correlate predictions
Hmisc::rcorr(reshape2::acast(mlp.rdp.retrain.models,
                             model ~ rowIndex, value.var="pred")%>%t() ,
             type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))


# :: append to mapping ----------------------------------------------------
# add final variables
lsarp.uc.metadata = metadata.lsarp.uc.stool

unique(lsarp.uc.metadata$HM) %in% excluded.lsarp.uc
# good

# add timing
lsarp.uc.metadata$timing = ifelse(lsarp.uc.metadata$trial.stool.timing == "baseline at time of RapidAIM collection", "0M",
                                  ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 1 post start of product", "1M",
                                         ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 2 post start of product", "2M",
                                                ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 3 post start of product", "3M",
                                                       ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 4 post start of product", "4M",
                                                              ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 5 post start of product", "5M",
                                                                     ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 6 post start of product", "6M",
                                                                            ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 8 post start of product", "8M",
                                                                                   ifelse(lsarp.uc.metadata$trial.stool.timing == "Month 10 post start of product", "10M",
                                                                                          ifelse(lsarp.uc.metadata$trial.stool.timing %in% c("Month 12 end of study", "Month 12 end of study / OARS baseline"), "12M", NA))))))))))
lsarp.uc.metadata$timing = factor(lsarp.uc.metadata$timing, levels=c("0M", "1M", "2M", "3M", "4M", "5M", "6M", "8M", "10M", "12M"))                                                                                   

# append RS selection
lsarp.uc.metadata = merge(lsarp.uc.metadata,
                          data.frame(lsarp.uc.rs.selections[,c("HM", "Selected")] %>% distinct() %>%
                                       mutate(rs_selected = Selected)), by="HM")

# add rs.col to grey out non-RS
lsarp.uc.metadata$rs.col = ifelse(lsarp.uc.metadata$phase == "onRS", as.character(lsarp.uc.metadata$rs_selected), "1")
lsarp.uc.metadata$rs.col = factor(lsarp.uc.metadata$rs.col, levels=c(1, rs.names))

metadata.lsarp.uc.stool = lsarp.uc.metadata

# now prepare ASV mapping
metadata.lsarp.uc.stool.asv = merge(lsarp.uc.metadata,
                                    lsarp.uc.asv.median.alpha, by="standard.name")
metadata.lsarp.uc.stool.asv = merge(metadata.lsarp.uc.stool.asv,
                                    lsarp.uc.butyrogens.i, by="standard.name")
metadata.lsarp.uc.stool.asv = merge(metadata.lsarp.uc.stool.asv,
                                    lsarp.uc.butyrogens.ii, by="standard.name")
metadata.lsarp.uc.stool.asv = merge(metadata.lsarp.uc.stool.asv,
                                    lsarp.uc.asv.functionalredundancy, by="standard.name")
metadata.lsarp.uc.stool.asv = merge(metadata.lsarp.uc.stool.asv,
                                    lsarp.uc.asv.mlp.preds, by="standard.name")
rownames(metadata.lsarp.uc.stool.asv) = metadata.lsarp.uc.stool.asv$standard.name



# :: save data -------------------------------------------------------------


# save these, load into analysis script
save(
  # colors and vectors
  labelcolors,
  rs.names,
  rs.names.pbs,
  # patient list
  lsarp.uc,
  # RS selections
  lsarp.uc.rs.selections,
  # Mapping
  metadata.lsarp.uc.stool,
  # Mapping + ASV data
  metadata.lsarp.uc.stool.asv,
  # ASV Tables
  lsarp.uc.asv.median.tax,
  lsarp.uc.asv.median.glom,
  # destination
  file = "./2026_01_08_lsarp_uc_data.Renv")


# :: ----------------------------------------------------------------------


# >> Analysis -------------------------------------------------------------

load(file = "./2026_01_08_lsarp_uc_data.Renv")


# ``` plot: stool timings -------------------------------------------------

# same as above in processing

lsarp.uc.stool.plot = ggplot(lsarp.uc.redcap %>% subset(!HM %in% excluded.lsarp.uc) %>%
                               # remove post-flare
                               subset(lsarp.days <= flare_day | is.na(flare_day)),
                             aes(x=lsarp.days, y=reorder(HM, flare_day)))+
  annotate("rect", xmin=0, xmax=max(subset(lsarp.uc.redcap, phase%in% c("preRS", "onRS"))$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  # background point
  geom_point(fill="white", shape=21, size=3)+
  # actual point
  geom_point(aes(fill=phase), 
             #alpha = ifelse(lsarp.days <= flare_day, "1", "0")),
             color="white",
             shape=21, size=5)+
  # add flare date
  geom_point(data=lsarp.uc.redcap[,c("flare_day", "lsarp.days","HM", "group")] %>% 
               subset(!HM %in% excluded.lsarp.uc) %>% distinct(),
             aes(x=flare_day, y=reorder(HM, (lsarp.days))), shape=4, size=5)+
  scale_fill_manual(values=c("black",2, "grey"))+
  scale_alpha_manual(values=c(0.25,1))+
  #geom_text(aes(label=missing.16s) ,color="white", nudge_y=-0.1, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-2, nchar(standard.name))), vjust=-1, size=2)+
  theme_classic()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group,nrow=1, scales="free")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.major.y=element_line(color="grey", linewidth=0.2),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.uc.stool.plot

# 5 samples missing


# :: Add RS end ---------------------------------------------------------------


# note: need to add custom "RS END BASELINE" for washout analyses
# because 5M is not the last stool for a few HMs
metadata.lsarp.uc.stool$rs.end = ifelse(
  # for HMs with end of RS at 4M:
  metadata.lsarp.uc.stool$HM %in% c("HM0871", "HM0941"), "4M",
  # for HMs with end of RS at 5M:
  ifelse(metadata.lsarp.uc.stool$HM %in% c("HM0869", "HM0897", "HM0898", "HM0914", "HM0947", "HM1033", "HM1053" ), "5M",
         NA)) # NA for those who did not make it to 5M

metadata.lsarp.uc.stool = metadata.lsarp.uc.stool %>%
  mutate(rs.end = ifelse(rs.end == timing, "rs.end", NA))
subset(metadata.lsarp.uc.stool, rs.end == "rs.end")[,c("HM", "timing", "rs.end", "phase")]

# went through 1 by 1 and recorded when the last stool on treatment (or + 3 days) was 4M or 5M
unique(metadata.lsarp.uc.stool$HM)
subset(metadata.lsarp.uc.stool, HM %in% unique(metadata.lsarp.uc.stool$HM)[1:15] & timing %in% c("3M", "4M", "5M", "6M"))[,c("HM", "lsarp.off", "timing")]

# repeat for .asv
metadata.lsarp.uc.stool.asv$rs.end = ifelse(
  # for HMs with end of RS at 4M:
  metadata.lsarp.uc.stool.asv$HM %in% c("HM0871", "HM0941"), "4M",
  # for HMs with end of RS at 5M:
  ifelse(metadata.lsarp.uc.stool.asv$HM %in% c("HM0869", "HM0897", "HM0898", "HM0914", "HM0947", "HM1033", "HM1053" ), "5M",
         NA)) # NA for those who did not make it to 5M

metadata.lsarp.uc.stool.asv = metadata.lsarp.uc.stool.asv %>%
  mutate(rs.end = ifelse(rs.end == timing, "rs.end", NA))
subset(metadata.lsarp.uc.stool.asv, rs.end == "rs.end")[,c("HM", "timing", "rs.end", "phase")]



# :: RS selections --------------------------------------------------------

lsarp.uc.rs.selections.rs.plot <- ggplot(subset(lsarp.uc.rs.selections,group=="RS"),
                                         aes(x=1, y=Z_score))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add points
  geom_point(aes(fill=RS_Name), shape=21,color="white", size=3, alpha=1)+
  # overlay selected RS (for clarity)
  geom_point(data = subset(lsarp.uc.rs.selections, group=="RS" & RS == Selected), 
             aes(fill=RS_Name), color="white", size=3, shape=21, alpha=1)+
  # add lines
  geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6)+
  # label selected RS
  geom_label(data = subset(lsarp.uc.rs.selections,  group=="RS" & RS_Name == Selected), 
             aes(label=RS_Name, color=RS_Name, x=1, y=3, vjust=0.55),
             size=2.5)+
  scale_y_continuous(limits=c(-2.2,3.2))+
  geom_hline(yintercept=1, linetype=2, color="red", alpha=0.5)+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Butyrogen Z-Score", title="RS")+
  facet_wrap(~HM, nrow=3)+
  theme_minimal()+theme(legend.position="none",
                        axis.text.x=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.minor.y=element_blank(),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.uc.rs.selections.rs.plot

lsarp.uc.rs.selections.placebo.plot <- ggplot(subset(lsarp.uc.rs.selections,group=="Plac") %>%
                                                mutate(Group = ifelse(group == "Plac", "Placebo", "RS")),
                                              aes(x=1, y=Z_score))+
  
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add points
  geom_point(aes(fill=RS_Name), shape=21, color="white", size=3, alpha=1)+
  # overlay selected RS (for clarity)
  geom_point(data = subset(lsarp.uc.rs.selections, group=="Plac" & RS == Selected), 
             aes(fill=RS_Name), color="white", size=3,  shape=21, alpha=1)+
  # add lines
  geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6)+
  # label selected RS
  geom_label(data = subset(lsarp.uc.rs.selections,  group=="Plac" & RS_Name == Selected), 
             aes(label=RS_Name, color=RS_Name, x=1, y=3, vjust=0.55),
             size=2.5)+
  scale_y_continuous(limits=c(-2.2,3.2))+
  geom_hline(yintercept=1, linetype=2, color="red", alpha=0.5)+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Butyrogen Z-Score", title="Placebo")+
  facet_wrap(~HM, nrow=3)+
  theme_minimal()+theme(legend.position="none",
                        axis.text.x=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.minor.y=element_blank(),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.uc.rs.selections.placebo.plot

# break down frequency of RS's being selected, per timepoint (Stacked barplot)
lsarp.uc.rs.selections.frequences = lsarp.uc.rs.selections %>%
  dplyr::select(HM, group, Selected) %>% distinct() %>% data.frame() %>%
  dplyr::select(Selected, group) %>%
  table() %>% data.frame() %>%
  group_by(group) %>%
  mutate(perc = Freq / sum(Freq)) %>% data.frame() 

lsarp.uc.rs.selections.frequences$Selected = factor(lsarp.uc.rs.selections.frequences$Selected, levels=rs.names)

lsarp.uc.rs.selections.frequences.plot = ggplot(lsarp.uc.rs.selections.frequences %>%
                                                  mutate(Group = ifelse(group == "Plac", "Placebo", "RS")) %>%
                                                  mutate(Group = factor(Group, levels=c("RS", "Placebo"))),
                                                aes(x=Selected, y=Freq))+
  geom_bar(stat="identity", fill="white", alpha=1,linewidth=0)+
  geom_bar(stat="identity", position="stack", 
           aes(fill=Selected), alpha=1,
           color="white", linewidth=0)+
  scale_fill_manual(values= labelcolors$cols[c(1:9)])+
  scale_y_continuous(breaks=seq(from=0, to=13, by=1))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y=element_line(color="grey", linewidth=0.25),
        strip.text = element_text(size=12),
        strip.background = element_rect(
          color="black"),
        legend.position="none")+
  labs(x="", y="Frequency Selected", fill="")+
  facet_wrap(~Group,  nrow=2)
lsarp.uc.rs.selections.frequences.plot


# ``` plot: rs selections -------------------------------------------------

lsarp.uc.rs.selections.placebo.plot+
  lsarp.uc.rs.selections.rs.plot+ 
  lsarp.uc.rs.selections.frequences.plot+
  patchwork::plot_layout(widths=c(4,2,2))


# :: Survival Curve -------------------------------------------------------

# only take treatment period
lsarp.uc.survival.data = subset(metadata.lsarp.uc.stool, phase != "postRS")[,c("HM","group","lsarp.days", "flare_day")] %>% unique()
lsarp.uc.survival.data$event = ifelse(is.na(lsarp.uc.survival.data$flare_day), 0, 1)
lsarp.uc.survival.data = lsarp.uc.survival.data %>%
  group_by(HM) %>% 
  mutate(time = ifelse(event == 0, max(lsarp.days), flare_day)) %>%
  slice_max(lsarp.days, n = 1)

# first, see if sig better prior to scope
lsarp.uc.survival <- survival::survfit(Surv(time, event) ~ group, 
                                           data = lsarp.uc.survival.data)
lsarp.uc.logrank <- survival::survdiff(Surv(time, event) ~ group, 
                                                 data = lsarp.uc.survival.data)
lsarp.uc.logrank


# ```plot: survival curve -------------------------------------------------

lsarp.uc.survival.plot = survminer::ggsurvplot(lsarp.uc.survival, 
                                                   data = lsarp.uc.survival.data,
                                                   pval = TRUE,             # show log-rank p-value
                                                   conf.int = TRUE,         # show confidence intervals
                                                   #risk.table = TRUE,       # show number at risk
                                                   conf.int.style="ribbon",
                                                   xlab = "Days since starting Product",
                                                   ylab = "Remission probability",
                                                   pval.size=4,
                                                   censor.shape=124,
                                                   legend="right",
                                                   legend.labs=c("Placebo", "RS"),
                                                   legend.title="Group",
                                                   #palette = c("#E7B800", "#2E9FDF"),
                                                   ggtheme = theme_classic())
lsarp.uc.survival.plot


# :: Fecal calprotectin ---------------------------------------------------


# does fecal calprotectin increase or decrease over RS treatment
# stats
lsarp.uc.stats.fcal.treatment = lmerTest::lmer(scale(log10(fecalcal_res)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                            subset(metadata.lsarp.uc.stool, 
                                                   phase %in% c("preRS", "onRS"))) %>%
  summary() %>% coef()

lsarp.uc.stats.fcal.washout = lmerTest::lmer(scale(log10(fecalcal_res)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                          subset(metadata.lsarp.uc.stool, 
                                                 phase == "postRS" | rs.end == "rs.end")) %>%
  summary() %>% coef()


# plot
lsarp.uc.fcal.plot = ggplot(subset(metadata.lsarp.uc.stool, !is.na(fecalcal_res)),
                                  aes(x=lsarp.days, 
                                      y=fecalcal_res))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                           phase%in%c("preRS", "onRS"))$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(span=1, color="white")+
  geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
  scale_shape_manual(values=c(23,21))+
  scale_alpha_manual(values=c(0.5, 1))+
  geom_hline(yintercept=250, color="red")+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Fecal Calprotectin (μg/g)",
       title=paste(paste("Treatment p:", round(lsarp.uc.stats.fcal.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.uc.stats.fcal.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.uc.fcal.plot
# no sig impact


# :: ASV Analysis ---------------------------------------------------------


# :: :: Richness ----------------------------------------------------------

  
  # does fecal calprotectin increase or decrease over RS treatment
  # stats
lsarp.uc.stats.richness.treatment = lmerTest::lmer(scale(log10(richness)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                 subset(metadata.lsarp.uc.stool.asv, 
                                                        phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
lsarp.uc.stats.richness.washout = lmerTest::lmer(scale(log10(richness)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                               subset(metadata.lsarp.uc.stool.asv, 
                                                      phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.richness.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(richness)),
                              aes(x=lsarp.days, 
                                  y=richness))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Richness",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.richness.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.richness.washout[5,5], 3))))+
    facet_wrap(~group)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.richness.plot
  # no sig impact
  

# :: :: Shannon -----------------------------------------------------------

  
  # does fecal calprotectin increase or decrease over RS treatment
  # stats
  lsarp.uc.stats.shannon.treatment = lmerTest::lmer(scale((shannon)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                     subset(metadata.lsarp.uc.stool.asv, 
                                                            phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.shannon.washout = lmerTest::lmer(scale((shannon)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                                   subset(metadata.lsarp.uc.stool.asv, 
                                                          phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.shannon.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(shannon)),
                                  aes(x=lsarp.days, 
                                      y=shannon))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    #geom_text(aes(label=HM))+
    labs(x="Days since starting Product", y="Shannon Diversity",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.shannon.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.shannon.washout[5,5], 3))))+
    facet_wrap(~group)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.shannon.plot
  # no sig impact
  
  # :: :: Butyrogens I -----------------------------------------------------------
  
  
  # does fecal calprotectin increase or decrease over RS treatment
  # stats
  lsarp.uc.stats.but.i.treatment = lmerTest::lmer(scale(log10(but.i)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                    subset(metadata.lsarp.uc.stool.asv, 
                                                           phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.but.i.washout = lmerTest::lmer(scale(log10(but.i)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                                  subset(metadata.lsarp.uc.stool.asv, 
                                                         phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.but.i.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(but.i)),
                                 aes(x=lsarp.days, 
                                     y=but.i*100))+
    scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Butyrogens (%)",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.but.i.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.but.i.washout[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=12),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.but.i.plot
  # no sig impact
  
  # :: :: Butyrogens II -----------------------------------------------------------
  
  
  # stats
  lsarp.uc.stats.but.ii.treatment = lmerTest::lmer(scale(log10(but.ii)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                  subset(metadata.lsarp.uc.stool.asv, 
                                                         phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.but.ii.washout = lmerTest::lmer(scale(log10(but.ii)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                                subset(metadata.lsarp.uc.stool.asv, 
                                                       phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.but.ii.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(but.ii)),
                               aes(x=lsarp.days, 
                                   y=but.ii*100))+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
    scale_y_log10()+
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Kircher Butyrogens (%)",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.but.ii.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.but.ii.washout[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.but.ii.plot
  # no sig impact
  
  
  # :: :: Functional Redundancy -----------------------------------------------------------
  
  
  # stats
  lsarp.uc.stats.fr.treatment = lmerTest::lmer(scale((fr)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                               subset(metadata.lsarp.uc.stool.asv, 
                                                      phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.fr.washout = lmerTest::lmer(scale((fr)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                             subset(metadata.lsarp.uc.stool.asv, 
                                                    phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.fr.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(fr)),
                            aes(x=lsarp.days, 
                                y=fr))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Functional Redundancy",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.fr.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.fr.washout[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.fr.plot
  # no sig impact
  
  # :: :: Functional Richness -----------------------------------------------------------
  
  
  # stats
  lsarp.uc.stats.fd.treatment = lmerTest::lmer(scale((fd)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                   subset(metadata.lsarp.uc.stool.asv, 
                                                          phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.fd.washout = lmerTest::lmer(scale((fd)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                                 subset(metadata.lsarp.uc.stool.asv, 
                                                        phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.fd.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(fd)),
                                aes(x=lsarp.days, 
                                    y=fd))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Functional Richness",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.fd.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.fd.washout[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.fd.plot
  # no sig impact
  
  # :: :: MLP -----------------------------------------------------------
  
  
  # stats
  lsarp.uc.stats.load.treatment = lmerTest::lmer(scale((load.asv)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                               subset(metadata.lsarp.uc.stool.asv, 
                                                      phase %in% c("preRS", "onRS"))) %>%
    summary() %>% coef()
  
  lsarp.uc.stats.load.washout = lmerTest::lmer(scale((load.asv)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                             subset(metadata.lsarp.uc.stool.asv, 
                                                    phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  
  # plot
  lsarp.uc.load.plot = ggplot(subset(metadata.lsarp.uc.stool.asv, !is.na(load.asv)),
                            aes(x=lsarp.days, 
                                y=load.asv))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool, 
                                             phase%in%c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    scale_alpha_manual(values=c(0.5, 1))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Microbial Load",
         title=paste(paste("Treatment p:", round(lsarp.uc.stats.fd.treatment[5,5], 3)), 
                     paste("  Washout p:", round(lsarp.uc.stats.fd.washout[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.load.plot
  # no sig impact
  
  

# :: :: Bray-Curtis -------------------------------------------------------


  # remove non-compliant (not necessary)
  
  # filter!
  ncol(lsarp.uc.asv.median.glom) # 549 taxa
  lsarp.uc.asv.median.glom.pa = lsarp.uc.asv.median.glom
  lsarp.uc.asv.median.glom.pa[lsarp.uc.asv.median.glom.pa > 0] = 1
  lsarp.uc.asv.median.glom = lsarp.uc.asv.median.glom.pa[,colSums(lsarp.uc.asv.median.glom.pa) > nrow(lsarp.uc.asv.median.glom.pa)*0.2]
  ncol(lsarp.uc.asv.median.glom) # 58 taxa
  
  # calculate Bray-Curtis dissimilarities
  lsarp.uc.asv.bray = vegan::vegdist(lsarp.uc.asv.median.glom[rownames(lsarp.uc.asv.median.glom) %in% 
                                                               metadata.lsarp.uc.stool.asv$standard.name,], method="bray") 
  # perform PCoA
  lsarp.uc.asv.pcoa = ape::pcoa(lsarp.uc.asv.bray)
  # extract data from pcoa
  lsarp.uc.asv.pcoa.df = data.frame(lsarp.uc.asv.pcoa$vectors[,c(1:2)])
  lsarp.uc.asv.pcoa.df$standard.name = rownames(lsarp.uc.asv.pcoa.df)
  # add metadata
  lsarp.uc.asv.pcoa.df = merge(lsarp.uc.asv.pcoa.df,
                               metadata.lsarp.uc.stool.asv, by="standard.name")
  # extract variance explained
  lsarp.uc.asv.pcoa.var_exp = lsarp.uc.asv.pcoa$values[c(1:2),2]
  lsarp.uc.asv.pcoa.df$var1 = round(lsarp.uc.asv.pcoa.var_exp[1]*100, digits=2)
  lsarp.uc.asv.pcoa.df$var2 = round(lsarp.uc.asv.pcoa.var_exp[2]*100, digits=2)
  
  # clean up "lsarp.ucon.rs" variable
  
  set.seed(25)
  t1 = Sys.time()
  lsarp.uc.asv.permanova = vegan::adonis2(lsarp.uc.asv.bray ~ group*on.rs + ave.fiber,
                                       lsarp.uc.asv.pcoa.df %>% mutate(on.rs = ifelse(phase == "onRS", "onRS", "offRS")),
                                       strata = lsarp.uc.asv.pcoa.df$HM,
                                       by="margin")
  t2 = Sys.time()
  t2 - t1
  lsarp.uc.asv.permanova # not-sig, but very sig association with fiber
  
  lsarp.uc.asv.pcoa.plot <- ggplot(
    data=lsarp.uc.asv.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
    aes(x=Axis.1, y=Axis.2))+
    geom_path(aes(group=HM), alpha=0.5, linetype=2, linewidth=0.3) + 
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1", "2")), color="white", size=2.5)+
    scale_shape_manual(values=c(23,21))+
    #scale_alpha_manual(values=c(1,0.5))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    #geom_text(aes(label=HM))+
    labs(x=paste("Axis 1: ", round(unique(lsarp.uc.asv.pcoa.df$var1), digits=2), "%", sep=""), 
         y=paste("Axis 2: ", round(unique(lsarp.uc.asv.pcoa.df$var2), digits=2), "%", sep=""),
         title=paste(paste("Treatment R²: ", round(data.frame(lsarp.uc.asv.permanova)[2,3], 3)*100, "%",
                           "  p: ", round(data.frame(lsarp.uc.asv.permanova)[2,5], 3), sep="")), sep="") + 
    facet_wrap(~group)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.asv.pcoa.plot
  

# :: :: Beta-Diversity -------------------------------------------------------

  
  beta.trajectory.data = beta.trajectory(lsarp.uc.asv.bray)
  
  beta.trajectory.data =   merge(metadata.lsarp.uc.stool.asv, 
                                 beta.trajectory.data, by="standard.name")
  
  
  # stats
  stats.beta.between.treatment = lmerTest::lmer(scale((between.beta)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                subset(beta.trajectory.data, phase %in% c("preRS", "postRS"))) %>%
    summary() %>% coef()
  
  stats.beta.between.washout = lmerTest::lmer(scale((between.beta)) ~ group*scale(lsarp.days)+  ave.fiber + (1|HM),
                                              subset(beta.trajectory.data,  phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  # plot
  metadata.lsarp.uc.stool.asv.between.plot = ggplot(beta.trajectory.data %>%
                                                      mutate(rs.col = factor(rs.col, levels=c(1, rs.names))),
                                                 aes(x=lsarp.days, 
                                                     y=between.beta))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool.asv,  phase %in% c("preRS", "onRS"))$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    #scale_alpha_manual(values=c(1,0.5))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Bray-Curtis Dissimilarity",
         title=paste(paste("Treatment p:", round(stats.beta.between.treatment[5,5], 3)), 
                     paste("  Washout p:", round(stats.beta.between.washout[5,5], 3))))+
    facet_wrap(~group)+
    #geom_text(aes(label=rs.col), size=2)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  metadata.lsarp.uc.stool.asv.between.plot
  
  
  
  # stats
  stats.beta.within.treatment = lmerTest::lmer(scale((within.beta)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                               subset(beta.trajectory.data, phase %in% c("preRS", "onRS")))  %>%
    summary() %>% coef()
  
  stats.beta.within.washout = lmerTest::lmer(scale((within.beta)) ~ group*scale(lsarp.days)+  ave.fiber + (1|HM),
                                             subset(beta.trajectory.data,  phase == "postRS" | rs.end == "rs.end")) %>%
    summary() %>% coef()
  
  # plot
  metadata.lsarp.uc.stool.asv.within.plot = ggplot(beta.trajectory.data %>%
                                                     mutate(rs.col = factor(rs.col, levels=c(1, rs.names))),
                                                aes(x=lsarp.days, 
                                                    y=within.beta))+
    #scale_y_log10()+
    annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.uc.stool.asv, phase=="onRS")$lsarp.days),
             ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_smooth(span=1, color="white")+
    geom_point(aes(fill=rs.col, shape=ifelse(timing == "0M", "1","2")),color="white", size=3)+
    scale_shape_manual(values=c(23,21))+
    #scale_alpha_manual(values=c(1,0.5))+
    scale_fill_manual(values=rs.colors, na.value="grey")+
    labs(x="Days since starting Product", y="Mean Bray-Curtis Dissimilarity",
         title=paste(paste("Treatment p:", round(stats.beta.within.treatment[5,5], 3)), 
                     paste("  Washout p:", round(stats.beta.within.washout[5,5], 3))))+
    facet_wrap(~group)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  metadata.lsarp.uc.stool.asv.within.plot
  
  

# :: :: Linear Model -------------------------------------------------------------

  # instead of using a function, I'll just write out the code here
  
  # first, filter
  dim(lsarp.uc.asv.median.glom)
  
  lsarp.uc.asv.median.glom.pa = lsarp.uc.asv.median.glom
  lsarp.uc.asv.median.glom.pa[lsarp.uc.asv.median.glom.pa > 0] = 1
  # detected in 20% of samples
  lsarp.uc.asv.median.glom.filtered = lsarp.uc.asv.median.glom[,colSums(lsarp.uc.asv.median.glom.pa)> 
                                                                 nrow(lsarp.uc.asv.median.glom.pa)*0.20]
  lsarp.uc.pseudo = min(lsarp.uc.asv.median.glom.filtered[lsarp.uc.asv.median.glom.filtered!=0])/2
  lsarp.uc.pseudo
  
  # loop through taxa
  lsarp.uc.asv.lmer = do.call(rbind, lapply(colnames(lsarp.uc.asv.median.glom.filtered), function(taxa){
    print(taxa)
    
    # identify min count
    min.val = min(lsarp.uc.asv.median.glom.filtered[lsarp.uc.asv.median.glom.filtered!=0])
    
    # subset taxa data
    data.subset = lsarp.uc.asv.median.glom.filtered[,taxa] %>% data.frame()
    data.subset$standard.name = rownames(lsarp.uc.asv.median.glom.filtered)
    # log2 transform w pseudo
    colnames(data.subset)[1] = "feature"
    data.subset$feature = log2(data.subset$feature + min.val/2)
    # merge with meta
    data.subset = merge(data.subset,
                        metadata.lsarp.uc.stool.asv, by="standard.name")
    # subset through time points
    do.call(rbind, lapply(c("treatment", "washout"),function(time){
      if(time == "washout"){
        data.subset = subset(data.subset, phase == "postRS" | timing == "5M")
      } else{
        data.subset = subset(data.subset, phase != "postRS")
      }
      
      # apply catch for low variance taxa
      if(sum(data.subset$feature>log2(min.val))>=(nrow(data.subset)*0.2)){
        
      
      lmer.output = lmerTest::lmer(feature ~ group*lsarp.days + (1|HM), data.subset) %>% summary() %>% coef()
      lmer.output = data.frame(lmer.output)[4,c(1,5)]
      lmer.output$taxa = taxa
      colnames(lmer.output)[1] = "estimate"
      colnames(lmer.output)[2] = "pval"
      lmer.output$phase = time
      lmer.output
      } else {lmer.output = data.frame(estimate = NA, pval = NA, taxa = taxa, phase = time)}
    }))
  }))
    
  lsarp.uc.asv.lmer$padj = p.adjust(lsarp.uc.asv.lmer$pval, method="BH")
  lsarp.uc.asv.lmer = subset(lsarp.uc.asv.lmer, !is.na(pval))
  lsarp.uc.asv.lmer = lsarp.uc.asv.lmer %>% arrange(pval)

  lsarp.uc.asv.lmer
  # B. sanguini is sig enriched in RS
  
# fix names
  lsarp.uc.asv.lmer$taxa =  ifelse(grepl("s__", (lsarp.uc.asv.lmer$taxa)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.uc.asv.lmer$taxa)))), sep=""),
                                                            paste("(", substr((lsarp.uc.asv.lmer$taxa), 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.uc.asv.lmer$taxa)))), sep=""))
  
  
  
  lsarp.uc.asv.lmer.interactions.plot = ggplot(lsarp.uc.asv.lmer %>%
                                              mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                            aes(x=estimate, y=padj))+
    geom_point(shape=21, aes(fill=estimate))+
    geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
    scale_fill_gradient2(low="blue", high="red")+
    ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                             size=2.5)+
    ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                             size=2.5, color="grey")+
    scale_y_continuous(transform=neg_log10_trans,
                       breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
    labs(x="Interaction Coefficient",
         y="FDR",
         title="ASV Interaction")+
    facet_wrap(~phase)+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12),
                          strip.text = element_text(size=10),
                          strip.background = element_rect(
                            color="black"))
  lsarp.uc.asv.lmer.interactions.plot


# :: :: Heatmaps  -------------------------------------------------------

  
  # subset to nominally sig taxa
  lsarp.uc.asv.lmer.heatmap = lsarp.uc.asv.median.tax[,colnames(lsarp.uc.asv.median.tax) %in% 
                                                           subset(lsarp.uc.asv.lmer, pval<0.05 & phase == "treatment")$taxa]
  # processing: TSS, pseudo, log2
  lsarp.uc.asv.lmer.heatmap = lsarp.uc.asv.lmer.heatmap / 50000
  lsarp.uc.asv.lmer.heatmap = log2(lsarp.uc.asv.lmer.heatmap+min(lsarp.uc.asv.lmer.heatmap[lsarp.uc.asv.lmer.heatmap!=0])/2)
  # melt
  lsarp.uc.asv.lmer.heatmap = lsarp.uc.asv.lmer.heatmap %>% reshape2::melt()
  colnames(lsarp.uc.asv.lmer.heatmap) = c("standard.name", "taxa", "value")
  # merge
  lsarp.uc.asv.lmer.heatmap = merge(lsarp.uc.asv.lmer.heatmap,
                                    metadata.lsarp.uc.stool.asv, by="standard.name")
  # Log2FC earliest vs latest on RS
  lsarp.uc.asv.lmer.heatmap = lsarp.uc.asv.lmer.heatmap %>%
    subset(phase != "postRS") %>%
    group_by(HM, taxa) %>%
    filter(lsarp.days == min(lsarp.days) | lsarp.days == max(lsarp.days)) %>%
    group_by(HM, taxa) %>%
    mutate(lfc.rs = value - value[lsarp.days == min(lsarp.days)]) %>%
    # keep only 6M sample
    slice_max(order_by=lsarp.days, n=1) %>% data.frame()
  # cast
  lsarp.uc.asv.lmer.heatmap = reshape2::acast(lsarp.uc.asv.lmer.heatmap,
                                                     HM ~ taxa, value.var="lfc.rs")
  pheatmap::pheatmap(scale(lsarp.uc.asv.lmer.heatmap),
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     color=colorRampPalette(c("blue","white", "red"))(100),
                     breaks=c(seq(min(lsarp.uc.asv.lmer.heatmap), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(lsarp.uc.asv.lmer.heatmap)/100, max(lsarp.uc.asv.lmer.heatmap), length.out=floor(100/2))))
  
  # HM annotations
  lsarp.uc.asv.lmer.heatmap.hm.map = data.frame(HM = rownames(lsarp.uc.asv.lmer.heatmap)) %>%
    merge(metadata.lsarp.uc.stool.asv[,c("HM", "group", "ave.fiber")] %>% distinct(),
          by="HM")
  
  rownames(lsarp.uc.asv.lmer.heatmap.hm.map) = lsarp.uc.asv.lmer.heatmap.hm.map$HM
  lsarp.uc.asv.lmer.heatmap.hm.map$HM = NULL
  colnames(lsarp.uc.asv.lmer.heatmap.hm.map) = c("Group", "Fiber intake")
  
  # Feature annotations
  lsarp.uc.asv.lmer.heatmap.feature.map = subset(lsarp.uc.asv.lmer, pval<0.05 & phase == "treatment")
  lsarp.uc.asv.lmer.heatmap.feature.map = lsarp.uc.asv.lmer.heatmap.feature.map[,c("taxa", "estimate")]
  rownames(lsarp.uc.asv.lmer.heatmap.feature.map) = lsarp.uc.asv.lmer.heatmap.feature.map$taxa
  lsarp.uc.asv.lmer.heatmap.feature.map$taxa = NULL
  colnames(lsarp.uc.asv.lmer.heatmap.feature.map) = c("Interaction")
  
  lsarp.uc.asv.lmer.heatmap
  
  # lastly, scale Log2FC by feature
  lsarp.uc.asv.lmer.heatmap = reshape2::melt(as.matrix(lsarp.uc.asv.lmer.heatmap)) %>%
    group_by(Var2) %>%
    mutate(slfc = scale(value)) %>%
    reshape2::acast(Var2 ~ Var1, value.var="slfc")
  
  # fix names
  rownames(lsarp.uc.asv.lmer.heatmap) =  ifelse(grepl("s__", rownames(lsarp.uc.asv.lmer.heatmap)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", rownames(lsarp.uc.asv.lmer.heatmap)))), sep=""),
                               paste("(", substr(rownames(lsarp.uc.asv.lmer.heatmap), 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", rownames(lsarp.uc.asv.lmer.heatmap)))), sep=""))
  rownames(lsarp.uc.asv.lmer.heatmap.feature.map) =  ifelse(grepl("s__", rownames(lsarp.uc.asv.lmer.heatmap.feature.map)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", rownames(lsarp.uc.asv.lmer.heatmap.feature.map)))), sep=""),
                                                paste("(", substr(rownames(lsarp.uc.asv.lmer.heatmap.feature.map), 1, 1), ") ", gsub("g__", "", gsub("f__", "", gsub("s__", "", rownames(lsarp.uc.asv.lmer.heatmap.feature.map)))), sep=""))
  
  
  
  # lfc heatmap
  pheatmap::pheatmap(lsarp.uc.asv.lmer.heatmap,
                     #scale_rows=T, # 
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     fontsize_row = 7,
                     fontsize_col = 7,
                     color=colorRampPalette(c("blue","white", "red"))(100),
                     breaks=c(seq(min(lsarp.uc.asv.lmer.heatmap), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(lsarp.uc.asv.lmer.heatmap)/100, max(lsarp.uc.asv.lmer.heatmap), length.out=floor(100/2))),
                     annotation_col=lsarp.uc.asv.lmer.heatmap.hm.map,
                     annotation_row = lsarp.uc.asv.lmer.heatmap.feature.map,
                     annotation_colors = list(`Group` = c(Placebo = gg_color_hue(2)[1],
                                                          RS = gg_color_hue(2)[2]),
                                              `Interaction` = colorRampPalette(c("blue","white", "red"))(100),
                                              `Fiber intake`= colorRampPalette(c("blue","white", "red"))(100)))
  

# >> Plots ----------------------------------------------------------------

  
  # final version
  (metadata.lsarp.uc.rs.plot+metadata.lsarp.uc.rs.supp.plot+
     metadata.lsarp.uc.fiber.plot+lsarp.uc.fiber.comparison.plot) %>%
    ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_ffq_fiber_rs.pdf",
           width=10, height=6.5)

  # Stools
  ((lsarp.uc.stool.plot + lsarp.uc.survival.plot$plot)+
    patchwork::plot_layout(nrow=1 ,widths=c(1.5,1)))%>%
    ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_stool_km_curve_plot_1.pdf",
           width=12, height=4, device = cairo_pdf)
  
  # RS Selections
  (lsarp.uc.rs.selections.placebo.plot+
      lsarp.uc.rs.selections.rs.plot+ 
    lsarp.uc.rs.selections.frequences.plot+
    patchwork::plot_layout(widths=c(4,2,2)))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_rs_selections_2.pdf",
         width=12, height=6, device = cairo_pdf)
  
  # Main plots + others below
  ((lsarp.uc.but.i.plot + lsarp.uc.fcal.plot)/
    (lsarp.uc.richness.plot+
       lsarp.uc.shannon.plot+
       lsarp.uc.fd.plot+lsarp.uc.fr.plot+
       metadata.lsarp.uc.stool.asv.between.plot+
       lsarp.uc.load.plot)+
      patchwork::plot_layout(heights=c(1,2)))%>%
    ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_variables_3.pdf",
           width=14, height=8, device = cairo_pdf)
  
  # Bray-Curtis + Volcano
    (lsarp.uc.asv.pcoa.plot/
      lsarp.uc.asv.lmer.interactions.plot) %>%
      ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_pcoa_omics_4.pdf",
             width=7, height=8, device = cairo_pdf)
  
    # Heatmap
    
  pheatmap::pheatmap(lsarp.uc.asv.lmer.heatmap,
                     scale_rows=T,
                     fontsize_row = 7,
                     fontsize_col = 7,
                     color=colorRampPalette(c("blue","white", "red"))(100),
                     breaks=c(seq(min(lsarp.uc.asv.lmer.heatmap), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(lsarp.uc.asv.lmer.heatmap)/100, max(lsarp.uc.asv.lmer.heatmap), length.out=floor(100/2))),
                     annotation_col=lsarp.uc.asv.lmer.heatmap.hm.map,
                     annotation_row = lsarp.uc.asv.lmer.heatmap.feature.map,
                     annotation_colors = list(`Group` = c(Placebo = gg_color_hue(2)[1],
                                                          RS = gg_color_hue(2)[2]),
                                              `Interaction` = colorRampPalette(c("blue","white", "red"))(100))) %>% 
    #`Fiber intake`= colorRampPalette(c("white","pink", "purple"))(100)))%>%
    ggsave(filename="./lsarp_plots/2026_01_08_lsarp_uc_heatmap_5.pdf",
           width=8, height=4, device = cairo_pdf)
  
  
  ## numbers
  metadata.lsarp.uc.stool.asv[,c("HM", "group")] %>% table()
  metadata.lsarp.uc.stool.asv[,c("HM", "group")] %>% distinct() %>% dplyr::select(group) %>% table()
  
  