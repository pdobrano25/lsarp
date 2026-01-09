### 2025_06_15  LSARP Stool Analysis

# Conducts major stool multi-omic analyses, including group-level, response group-level, machine learning.


# :: save -----------------------------------------------------------------

v.date = "2026_01_06" # saved with all projects included (some processing, analysis, proteomics, UC)
# save.image("./2026_01_08_lsarp_analysis.Renv")

# load("./2026_01_08_lsarp_analysis.Renv")

# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse")

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
omics.colors = c("ASV" = "#8DD3C7",
                 "Species" = "#FFFFB3",
                 "Pathway" =  "#BEBADA",
                 "COG" = "#FB8072",
                 "CAZy" = "#80B1D3",
                 "Metabolite" = "#FDB462")
# 
omics.shapes = c("ASV" = 21,
                 "Species" = 22,
                 "Pathway" =  23,
                 "COG" = 23,
                 "CAZy" = 23,
                 "Metabolite" = 24)


# >> Exclusions -----------------------------------------------------------

# after critically evaluating HMs to keep (and basically relying on Dave's judgement)
# remove these:

excluded.by.dave = c("HM1045") # "HM0998"? No, withdrew after flaring

# for reference, see: 2025_09_11_lsarp_cd_patient_list_exclusions

metadata.lsarp.stool = subset(metadata.lsarp.stool, !HM %in% excluded.by.dave)
metadata.lsarp.stool.asv = subset(metadata.lsarp.stool.asv, !HM %in% excluded.by.dave)

# metadata has now been filtered to included participants (n=26, good)
unique(metadata.lsarp.stool$HM)


# >> Add RS end --------------------------------------------------

# note: need to add custom "RS END BASELINE" for washout analyses
# because 5M is not the last stool for a few HMs
metadata.lsarp.stool$rs.end = ifelse(
  # for HMs with end of RS at 4M:
  metadata.lsarp.stool$HM %in% c("HM0860","HM0865"), "4M",
  # for HMs with end of RS at 5M:
  ifelse(metadata.lsarp.stool$HM %in% c("HM0878","HM0912", "HM0918", "HM0921", "HM0924",
                                        "HM0926", "HM0933", "HM0938", "HM0960", "HM0962",
                                        "HM0978", "HM0999", "HM1007", "HM1022", "HM1023",
                                        "HM1035", "HM1051"), "5M",
  NA)) # NA for those who did not make it to 5M

metadata.lsarp.stool = metadata.lsarp.stool %>%
  mutate(rs.end = ifelse(rs.end == timing, "rs.end", NA))
subset(metadata.lsarp.stool, rs.end == "rs.end")[,c("HM", "timing", "rs.end", "phase")]

# went through 1 by 1 and recorded when the last stool on treatment (or + 3 days) was 4M or 5M
unique(metadata.lsarp.stool$HM)
subset(metadata.lsarp.stool, HM %in% unique(metadata.lsarp.stool$HM)[15:26] & timing %in% c("3M", "4M", "5M", "6M"))[,c("HM", "lsarp.off", "timing")]

# repeat for .asv
metadata.lsarp.stool.asv$rs.end = ifelse(
  # for HMs with end of RS at 4M:
  metadata.lsarp.stool.asv$HM %in% c("HM0860","HM0865"), "4M",
  # for HMs with end of RS at 5M:
  ifelse(metadata.lsarp.stool.asv$HM %in% c("HM0878","HM0912", "HM0918", "HM0921", "HM0924",
                                        "HM0926", "HM0933", "HM0938", "HM0960", "HM0962",
                                        "HM0978", "HM0999", "HM1007", "HM1022", "HM1023",
                                        "HM1035", "HM1051"), "5M",
         NA)) # NA for those who did not make it to 5M

metadata.lsarp.stool.asv = metadata.lsarp.stool.asv %>%
  mutate(rs.end = ifelse(rs.end == timing, "rs.end", NA))
subset(metadata.lsarp.stool.asv, rs.end == "rs.end")[,c("HM", "timing", "rs.end", "phase")]



# :: ``` stool collections plot -----------------------------------------------

# visualize stools + timings
metadata.lsarp.stool.plot = ggplot(metadata.lsarp.stool %>% 
                                     arrange(lsarp.days) %>%
                                     # remove post-flare samples
                                     subset(lsarp.days <= exit.day | is.na(exit.day)),
                                   aes(x=lsarp.days, y=reorder(HM, (exit.day))))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  # add background point to look neater
  geom_point(fill = "white", shape=21, size=5)+
  # then plot real points (they overlap in the window, but look good when saved at full size)
  geom_point(aes(fill=lsarp.on.rs), color="white",
           #      alpha = ifelse(lsarp.days <= exit.day, "1", "0")),
             shape=21, size=5)+
  # now plot endoscopy
  geom_point(aes(x=scope.day),
             #      alpha = ifelse(lsarp.days <= exit.day, "1", "0")),
             shape=1, size=3)+
  # now plot recorded flare date
  geom_point(data=metadata.lsarp.stool[,c("exit.day","exit.reason", "lsarp.days","HM", "group")] %>% 
               distinct() %>% 
               group_by(HM)%>%
               filter(lsarp.days == max(lsarp.days)),
             aes(x=exit.day, y=reorder(HM, (lsarp.days)), shape=exit.reason), 
             size=4)+
  scale_fill_manual(values=c(2,"grey","black"))+
  #scale_alpha_manual(values=c(0.25,1))+
  scale_shape_manual(values=c(4))+
  # geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), color="white", nudge_y=-0.1, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-1, nchar(standard.name))), size=2)+
  theme_minimal()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group,nrow=1, scales="free")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        panel.grid.major.y=element_line(color="grey", linewidth=0.2),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.plot


## Visualize dataset completeness

# load mapping file that was sent to Sara/Figeys
mpx.original.list = read.csv("~/Downloads/2024_11_30_mpx_stool_list_r.csv")

# visualize stools + timings
ggplot(metadata.lsarp.stool %>%
         mutate(
           fcal.complete = ifelse(!is.na(fcal), "fcal.complete", "fcal.missing"),
           water.complete = ifelse(!is.na(stool_water_perc), "water.complete", "water.missing"),
           asv.complete = ifelse(standard.name %in% rownames(lsarp.asv.data.glom), "asv.complete", "asv.missing"),
           mgx.complete = ifelse(standard.name %in% rownames(lsarp.mgx.taxa), "mgx.complete", "mgx.missing"),
           mpx.complete = ifelse(standard.name %in% rownames(lsarp.cd.mpx.cazy.mat), "mpx.complete", "mpx.missing"),
           mpx.requested = ifelse(standard.name %in% subset(mpx.original.list, DEPLETED != "DEPLETED")$standard.name, "mpx.requested", "mpx.not.requested")),
       aes(x=lsarp.days, y=reorder(HM, (flare.day))))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(aes(fill=mpx.complete, shape=mpx.requested), 
             size=3)+
  scale_fill_manual(values=c("grey","2"))+
  scale_alpha_manual(values=c(0.25,1))+
  scale_shape_manual(values=c(21,22))+
  #geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), nudge_y=-0.15, size=6)+
  geom_text(aes(label=substr(standard.name, nchar(standard.name)-1, nchar(standard.name))), size=2)+
  theme_minimal()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group, nrow=1, scales="free")+
  theme_minimal()+theme(#legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))

## Plot statistic test groupings samples
ggplot()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(data = subset(metadata.lsarp.stool, lsarp.on.rs %in% c("pre.rs", "on.rs")),
             aes(x=lsarp.days, y=reorder(HM, (flare.day)), fill=1, alpha = ifelse(lsarp.days <= flare.day, "1", "0")),
             shape=21, size=5)+
  geom_point(data = subset(metadata.lsarp.stool, timing == "5M" | lsarp.on.rs == "post.rs"),
             aes(x=lsarp.days, y=reorder(HM, (flare.day)), fill=2, 
                 alpha = ifelse(lsarp.days <= flare.day, "1", "0")),
             shape=21, size=5)+
  #scale_fill_manual(values=c(2,"grey","black"))+
  scale_alpha_manual(values=c(0.25,1))+
  # geom_text(aes(label=ifelse(missing.16s == "missing", "*", "")), color="white", nudge_y=-0.1, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-1, nchar(standard.name))), size=2)+
  theme_minimal()+
  labs(x="Day since starting Product", y="", color="")+
  facet_wrap(~group,nrow=1, scales="free")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        panel.grid.major.y=element_line(),
                        strip.background = element_rect(
                          color="black"))

# >> Functions ------------------------------------------------------------

# this function applies a prevalence filter to omics data (using only samples of interest)
# applies a log(+pseudo) transform
# and optionally scales to 100%
lsarp.delta.omic.prepare = function(data = lsarp.mpx.cazy.mat,
                                    filter = T,
                              prev = 0.20,
                              #min.abun = 1,
                              pseudo = T,
                              normalize = F){

  # loop through taxa (present in >20% samples)
  n.samples = nrow(data)*prev
  # replace NA with 0
  data[is.na(data)] = 0
  # min value
  min.abun = min(data[data!=0])  # grants flexibility if data are counts or %
  if(filter == T){
  # filter
  data.filtered = data[,rowSums(apply(data, 1, function(x) {x >= min.abun}))>n.samples] %>% as.data.frame()
  }else{
    data.filtered = data
  }
  # scale to 100%
  if(normalize == T){
    data.filtered = as.data.frame(((data.filtered) / rowSums(data.filtered)) * 100)
  }
  # add pseudo (note: this order (scale before +pseudo), means total % will exceed 100)
  pseudo = min(data.filtered[data.filtered!=0])/2
  data.filtered = data.filtered+pseudo
  data.filtered$standard.name = rownames(data.filtered)
  #
  return(data.filtered)
}

# lmer wrapper that either tests (split = F) INTERACTION; or (split=T) SEPARATE LFC PER GROUP
# Note: split=T for response groups is not actually relevant because LSARP-RS had only 1 responder

lsarp.delta.omic.lmer = function(data = lsarp.lmer.asv.data,
                                 prev = 0.20, # prevalence filter, adjust to 0.80 for MBX
                                 split = FALSE) {
  
  meta.cols = c("HM","group","phase","timing", "standard.name", "lsarp.days","rs.end", "ave.fiber")
  
  if(split == FALSE){
    data = merge(data, 
                 metadata.lsarp.stool[,meta.cols], 
                 by="standard.name")
  }
  
  
  lm.output = do.call(rbind, lapply((colnames(data)[!colnames(data) %in% meta.cols]), function(x){
    print(x)
    # x = (colnames(data)[!colnames(data) %in% meta.cols])[1]
    # x = "c__Bacilli_A"
    # isolate per taxa and compare to lfc or response
    data.subset = data[,colnames(data) %in% c(x, meta.cols)] %>% as.data.frame()
    colnames(data.subset)[2] = "feature"

    # identify min val
    min.val = min((data[,!colnames(data) %in% meta.cols]))
    
    # log scale
    data.subset$feature = log2(data.subset$feature)
    
    # subset based on response level
    if(split == FALSE){
      lmer.output = do.call(rbind, lapply(c("treatment", "washout"), function(p){
       if(p == "washout"){
        data.response = subset(data.subset, phase == p | rs.end == "rs.end") %>%
          # and remove samples with only 1 time point in washout (5M, prior to flare)
          subset(!HM %in% excluded.by.dave)
       } else{
         data.response = subset(data.subset, phase == p)
       }
        
        # skip feature without 20% prev (or 80% metabolite) PER GROUP
        #prev.placebo = sum(subset(data.response, group == "Placebo")$feature>log2(min.val))>=(nrow(subset(data.response, group == "Placebo"))*prev)
        #prev.rs =  sum(subset(data.response, group == "RS")$feature>log2(min.val))>=(nrow(subset(data.response, group == "RS"))*prev)
        
        # skip feature without 20% prev (or 80% metabolite) IN TOTAL
        prev.total = sum(data.response$feature>log2(min.val))>=(nrow(data.response)*prev)

        if(sum(prev.total)>0){
          #
          lmer.output = lmerTest::lmer((feature) ~ lsarp.days*group +(1|HM), data.response) %>% summary() %>% coef()
          #
          lmer.output = data.frame(lmer.output)[4,c(1,5)]
          lmer.output$taxa = x
          colnames(lmer.output)[1] = "estimate"
          colnames(lmer.output)[2] = "pval"
          lmer.output$phase = p
          lmer.output
        } else {lmer.output = data.frame(estimate = NA, pval = NA, taxa = x, phase = p)}
    }))
    }
    return(lmer.output)
  }))
  lm.output$padj = p.adjust(lm.output$pval, method="BH")
  lm.output = lm.output %>% arrange(pval)
  return(lm.output)
}

# function for responders vs non-responders

lsarp.delta.omic.resp.lmer = function(data = lsarp.lmer.asv.data,
                                      prev = 0.20, # prevalence filter
                                      split = F){
  
  meta.cols = c("HM","group","phase","timing", "standard.name", "lsarp.days","rs.end", "ave.fiber", "flare.group")
  
  
  #if(split == FALSE){
    data = merge(data, lsarp.metadata.responders[,meta.cols], by="standard.name")
  
  lm.output = do.call(rbind, lapply((colnames(data)[!colnames(data) %in% meta.cols]), function(x){
    print(x)
    # x = colnames(data)[6]
    # isolate per taxa and compare to lfc or response
    data.subset = data[,colnames(data) %in% c(x, meta.cols)] %>% as.data.frame()
    colnames(data.subset)[2] = "feature"
    
    # identify min value
    min.val = min(data.subset$feature)
    
    # log scale
    data.subset$feature = log2(data.subset$feature)
    
    # subset based on flare status
    if(split == FALSE){
      lmer.output = do.call(rbind, lapply(unique(lsarp.metadata.responders$group), function(p){
        
        data.response = subset(data.subset, group == p)
        
        #  skip features with < 20% prev (or 80% metabolite) PER GROUP
        # prev.placebo = sum(subset(data.response, flare.group == "Remit")$feature>log2(min.val))>=(nrow(subset(data.response, flare.group == "Remit"))*prev)
        # prev.rs =  sum(subset(data.response, flare.group == "Relapse")$feature>log2(min.val))>=(nrow(subset(data.response, flare.group == "Relapse"))*prev)
        
        # skip features with < 20% prev (or 80% metabolite) ACROSS BOTH GROUPS
        prev.total = sum(data.response$feature > log2(min.val))>=(nrow(data.response)*prev)
        # I think it is more appropriate to do it this way, since we're building a single model where both groups are considered
        
        if(sum(prev.total)>0){
          #
          lmer.output = lmerTest::lmer((feature) ~ lsarp.days*flare.group + (1|HM), data.response) %>% summary() %>% coef()
          #
          lmer.output = data.frame(lmer.output)[4,c(1,5)]
          lmer.output$taxa = x
          colnames(lmer.output)[1] = "estimate"
          colnames(lmer.output)[2] = "pval"
          lmer.output$group = p
          lmer.output
        } else {lmer.output = data.frame(estimate = NA, pval = NA, taxa = x, group = p)}
      }))
    }
    return(lmer.output)
  }))
  lm.output$padj = p.adjust(lm.output$pval, method="BH")
  lm.output = lm.output %>% arrange(pval)
  return(lm.output)
}
  

# this function calculates within/between sample beta-diversity
# into a format compatible with standard variable analyses
beta.trajectory = function(data=lsarp.asv.bray){
  # melt
  data.subset = reshape2::melt(as.matrix(data))
  # note, I'm keeping the duplicate values (eg. X-Y and Y-X)
  # because attempting to remove them tends to introduce new errors (hard to fix with code)
  # and the average will be the same if all values are duplicated
  
  # but we will remove self comparisons
  data.subset = subset(data.subset, Var1 != Var2)
  
  # now calculate between sample beta
  
  # unpack sample names
  data.subset = tidyr::separate(data.subset, col=Var1, into=c("HM.A", "STL.A", "Number.A"), sep="-", remove=F)
  data.subset = tidyr::separate(data.subset, col=Var2, into=c("HM.B", "STL.B", "Number.B"), sep="-", remove=F)
  # calculate between sample diversity (between ALL samples)
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


# function to perform RF and associated analyses
rf.function = function(data.types = "All",
                       iters = 15,
                       balance = T,
                       reduce_features = F, # ttest
                       p.threshold = NA, # wilcox test p value threshold
                       n.threshold = NA, # top 10 features by pval
                       output = "AUC"){ # | ROC | importances | interactions
  
  # start timer
  t1 = Sys.time()
  
  rf.function.output = 
    
    # select data.type
    do.call(rbind, lapply(data.types, function(data.type){
      
      # select iteration (seed)
      do.call(rbind, lapply(1:iters, function(iter){
        
        # track progress
        print(paste(data.type, iter))
        
        # select data
        if(data.type == "ASV"){
          lsarp.data = lsarp.asv.baseline
        }
        if(data.type == "Species"){
          lsarp.data = lsarp.mgx.baseline
        }
        if(data.type == "Pathway"){
          lsarp.data = lsarp.kegg.baseline
        }
        if(data.type == "COG"){
          lsarp.data = lsarp.cog.baseline
        }
        if(data.type == "CAZy"){
          lsarp.data = lsarp.cazy.baseline
        }
        if(data.type == "Metabolite"){
          lsarp.data = lsarp.mbx.baseline
        }
        if(data.type == "Multi-Omic"){
          lsarp.data = lsarp.omics.baseline
        }
        
        # subset to matching data across omics (n=13 responders + non-responders, not n=14)
        lsarp.data = lsarp.data[lsarp.data$standard.name %in% lsarp.omics.baseline$standard.name,]
        
        # add potentially predictive features
        lsarp.data = merge(lsarp.data, 
                           metadata.lsarp.stool.asv[,c("standard.name","richness","shannon","but.i","but.ii","fd","fr", "load.asv","stool_water_perc")],
                           by="standard.name")
        # add RS + survival data
        lsarp.data = merge(lsarp.data, # Group so we can subset to RS only
                           metadata.lsarp.stool.asv[,c("Group", "RS_Name", "standard.name")],
                           by="standard.name")
        lsarp.data$HM = substr(lsarp.data$standard.name, 1, 6)
        lsarp.data$standard.name = NULL
        lsarp.data = merge(lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                           lsarp.data,
                           by="HM")
        # subset to RS group (n = 13; 8 relapse + 5 remit)
        lsarp.data = subset(lsarp.data, Group == "RS")
        # ready to run
        lsarp.data = lsarp.data %>%
          mutate(HM = NULL,
                 time = NULL,
                 Group = as.factor(Group),
                 RS_Name = as.factor(RS_Name),
                 flare.group = as.factor(flare.group))
        
        # for "output == interactions", subset features to important; somewhat // defunct
        if(output == "interactions"){
          lsarp.data = lsarp.data[,colnames(lsarp.data) %in% c("flare.group", lsarp.rf.models.importances.df$feature)]
        }
        
        # RUN LOOCV ANALYSIS
        # note: this do.call gives output of LOO predictions
        lsarp.loocv = do.call(rbind,lapply(1:nrow(lsarp.data),function(x){
          print(x)
          # first, we remove one of the samples
          data.subset = lsarp.data[-x,]
          # optional: then balance to smaller class size (either 5 or 4; 4 if a remiter is leaved out)
          if(balance == T){
            # balance data (undersample)
            set.seed(iter)
            data.subset  = data.subset %>%
              group_by(flare.group) %>%
              sample_n(size=min(table(data.subset$flare.group)), replace=F) %>% # n = min number in minority class
              as.data.frame() %>%
              mutate(flare.group = factor(flare.group, levels=c("Remit", "Relapse")))
          }
          
          # remove features with no variance (now that we've subsampled)
          no_var_cols = apply(data.subset[,!colnames(data.subset) %in% c("flare.group")], 2, sd) %>% data.frame() %>% subset(. == 0)
          data.subset = data.subset[,!colnames(data.subset) %in% rownames(no_var_cols)]
          
          # additionally, we can identify features with low (non-0) variance using caret::nearZeroVar
          low_var_cols <- caret::nearZeroVar(data.subset[,!colnames(data.subset) %in% c("flare.group")], 
                                             saveMetrics = TRUE)
          
          # select features with >= 20% unique values to keep consistent with other filtering
          data.subset = data.subset[,!colnames(data.subset) %in% 
                                      rownames(subset(low_var_cols, percentUnique <20))]
          
          # feature reduction // defunct
          if(reduce_features == T){
            # run wilcox and select sig features
            features = colnames(data.subset)
            # remove group + RS_Name
            features = features[!features%in%c("flare.group", "RS_Name")]
            
            test.results = do.call(rbind, lapply(features, function(feature){
              #feature = features[1]
              data.subset.feature = data.subset[,c(feature, "flare.group")]
              # run wilcox
              test.result = wilcox.test(subset(data.subset.feature, flare.group=="Relapse")[,1],
                                        subset(data.subset.feature, flare.group=="Remit")[,1])
              if(!is.na(p.threshold)){
                test.results = data.frame(feature = feature,
                                          pval = test.result$p.value) %>% arrange(pval) %>% subset(pval < p.threshold)
              }
              if(!is.na(n.threshold)){
                test.results = data.frame(feature = feature,
                                          pval = test.result$p.value) %>% arrange(pval) %>% slice_min(order_by=pval, n=n.threshold)
              }
            }))
            # apply filter
            data.subset = data.subset[,colnames(data.subset) %in% c(test.results$feature, "flare.group")]
            # and to full data (note: the "test" sample is here as df[x,])
            lsarp.data = lsarp.data[,colnames(lsarp.data) %in% c(test.results$feature, "flare.group")]
            
          }
          
          # build the RF model
          # use rfsrc so we can extract interaction scores later
          # use probability so we can calculate AUC/ROC
          set.seed(iter)
          model.rf = randomForestSRC::rfsrc(flare.group~ ., data.subset %>%
                                                   as.data.frame(), 
                                                 probability = TRUE)
          
          ## OUTPUT 1 = AUC | ROC (performances of model(s)); continues outside of loop
          if(output %in% c("AUC", "ROC")){
            # now we bring back the leaved out test sample
            pred = predict(model.rf, 
                           lsarp.data[x, ])
            
            output.df = data.frame(relapse = pred$predicted[1], # extract prob of "Relapse"
                                   true = factor(lsarp.data[x, ]$flare.group, levels=c("Relapse", "Remit")), # if "high", true
                                   index = x,
                                   iter = iter)
            return(output.df)
          }
          
          ## OUTPUT 2 = Importances
          if(output == "importances"){
            
            # calculate importances
            set.seed(iter)
            ints = randomForestSRC::vimp(model.rf, method = "permute", joint=F) 
            ints = data.frame(imp = ints$importance[,1]) %>%
              mutate(iter = iter,
                     index = x) %>%
              mutate(feature = rownames(ints$importance)) %>%
              arrange(imp)
            rownames(ints) = NULL
            return(ints)
          }
          
          ## OUTPUT 3 = Interactions
          if(output == "interactions"){
            # calculate interactions
            set.seed(iter)
            ints = randomForestSRC::find.interaction(model.rf,
                                                     method = "vimp", verbose=F)
            ints.df = ints %>% 
              as.data.frame() %>%
              mutate(features = rownames(.))%>%
              tidyr::separate(features, into=c("var1", "var2"), sep=":", remove=F)%>%
              # save df identifies
              mutate(
                iter = iter,
                index = x)
            return(ints.df)
          }
          
        })) # end loop for a single leaved out sample (still within loop for a specific iteration, though)
        
        # for "output==AUC", calculate AUC of that iteration
        if(output == "AUC"){
          # calculate LOOCV Accuracy
          # why? there are cases where the AUC is much higher than accuracy
          # this is because if you change your "confidence" threshold from 0.5 to 0.4, accuracy increases
          # but that (and those) outcome(s) is/are not considered if you use only 1 threshold
          lsarp.loocv = lsarp.loocv %>% 
            mutate(guess = ifelse(relapse >= 0.5, "Relapse", "Remit"),
                   correct = ifelse(true == guess, TRUE, FALSE),
                   accuracy = mean(correct))
          # calculate LOOCV AUC (hard code levels!)
          lsarp.loocv.auc = lsarp.loocv %>%
            mutate(auc = pROC::auc(true, relapse, 
                                   levels=c("Remit", "Relapse"),  # define case = "Relapse"
                                   direction="<")[1], # if probability > 0.5 (e.g.), then it is a case
                   iter = iter,
                   data.type = data.type) 
          return(lsarp.loocv.auc)
        }else{
          return(lsarp.loocv)
        }
      }))
    }))
  
  print(Sys.time() - t1)
  
  return(rf.function.output)
}
# stop timer outside

# function to calculate trajectory (lm coefficient) or Log2FC (baseline vs last recorded entry)

lsarp.delta = function(data = lsarp.fecalcal.omics.data[,colnames(lsarp.fecalcal.omics.data)%in% c(subset(metadata.lsarp.stool.omics.cor.bh, padj < 0.20)$feature)],
                       time.phase = "treatment",
                       responders = F,
                       type = "trajectory"){ # trajectory | log2fc
  if(type == "trajectory"){
    # do treatment first:
    data.output = do.call(rbind, lapply(1:(ncol(data)), function(col){
      print(col)
      data.subset = data[,c(col, ncol(data))]
      # subset through HM
      data.subset$standard.name = rownames(data.subset)
      feature.name = colnames(data.subset)[1]
      # change colname
      colnames(data.subset)[1] = "feature"
      # replace NA with 0
      data.subset$feature = ifelse(is.na(data.subset$feature), 0, data.subset$feature)
      # pseudo per feature
      feature.pseudo = min(data.subset$feature[data.subset$feature != 0])/2
      # log2 scale
      data.subset$feature = log2(data.subset$feature + feature.pseudo)
      
      data.subset = merge(data.subset, metadata.lsarp.stool.asv[,c("HM", "standard.name", "lsarp.days", "phase", "rs.end")], by="standard.name")
      
      data.output = do.call(rbind, parallel::mclapply(unique(data.subset$HM), function(hm){
        do.call(rbind, lapply(c("treatment", "washout"), function(time){
          # subset to treatment or washout + last treatment
          if(time == "treatment"){
            data.subset.hm = subset(data.subset, HM == hm & phase == "treatment")
          } else {
            data.subset.hm = subset(data.subset, HM == hm & phase == "washout")
            data.subset.hm = rbind(data.subset.hm,
                                   subset(data.subset, HM == hm & rs.end == "rs.end")) %>% as.data.frame() %>% distinct()
          }
        if(length(unique(data.subset.hm$feature)) >= 3){
        # if sufficient samples, run LM
        lm.output = lm(feature ~ lsarp.days, data.subset.hm)
        data.frame(coef = coef(lm.output)[2],
                   HM = hm,
                   feature = feature.name,
                   phase = time)
        } else{
          data.frame(coef = NA,
                     HM = hm,
                     feature = feature.name,
                     phase = time)
        }
        }))
      }))
    }))
      return(data.output)
    }
  if(type == "log2fc"){
  # calculate all Log2FC to baseline
    
    # merge all data
    lsarp.lmer.omics = merge(lsarp.lmer.asv.data,
                             lsarp.lmer.mgx.data, by="standard.name")
    lsarp.lmer.omics = merge(lsarp.lmer.omics,
                             lsarp.lmer.mpx.kegg.data%>%as.data.frame()%>%mutate(standard.name = rownames(.)), by="standard.name")
    lsarp.lmer.omics = merge(lsarp.lmer.omics,
                             lsarp.lmer.mpx.cog.data%>%as.data.frame()%>%mutate(standard.name = rownames(.)), by="standard.name")
    lsarp.lmer.omics = merge(lsarp.lmer.omics,
                             lsarp.lmer.mpx.cazy.data%>%as.data.frame()%>%mutate(standard.name = rownames(.)), by="standard.name")
    lsarp.lmer.omics = merge(lsarp.lmer.omics,
                             lsarp.lmer.mbx.data, by="standard.name")
    
    # subset to sig features
    if(responders == F){
    lmer.sig.features = 
      c(subset(lsarp.lmer.asv.interactions, padj < 0.20 & phase == time.phase)$taxa,
        subset(lsarp.lmer.mgx.interactions, padj < 0.20 & phase == time.phase)$taxa,
        subset(lsarp.lmer.mpx.kegg.interactions, padj < 0.20 & phase == time.phase)$taxa,
        subset(lsarp.lmer.mpx.cog.interactions, padj < 0.20 & phase == time.phase)$taxa,
        subset(lsarp.lmer.mpx.cazy.interactions, padj < 0.20 & phase == time.phase)$taxa,
        subset(lsarp.lmer.mbx.interactions, padj < 0.20 & phase == time.phase)$taxa)
    lmer.meta = metadata.lsarp.stool.asv[,c("standard.name","HM","RS_Name", "Group", "lsarp.days", "baseline", "phase")]
    }
    if(responders == T){
      lmer.sig.features = 
        c(subset(lsarp.lmer.asv.resp.interactions, padj < 0.20 & group == "RS")$taxa,
          subset(lsarp.lmer.mgx.resp.interactions, padj < 0.20 & group == "RS")$taxa,
          subset(lsarp.lmer.mpx.kegg.resp.interactions, padj < 0.20 & group == "RS")$taxa,
          subset(lsarp.lmer.mpx.cog.resp.interactions, padj < 0.20 & group == "RS")$taxa,
          subset(lsarp.lmer.mpx.cazy.resp.interactions, padj < 0.20 & group == "RS")$taxa,
          subset(lsarp.lmer.mbx.resp.interactions, padj < 0.20 & group == "RS")$taxa)
      lmer.meta = lsarp.metadata.responders[,c("standard.name","HM","RS_Name","phase", "Group", "lsarp.days", "baseline", "flare.group")]
    }
    # apply subset
    lsarp.lmer.omics = lsarp.lmer.omics[,colnames(lsarp.lmer.omics) %in% c("standard.name",
                                                                           lmer.sig.features)]
    # merge
    lsarp.lmer.omics = merge(lsarp.lmer.omics,lmer.meta, by="standard.name")
    # subset to phase
    lsarp.lmer.omics = subset(lsarp.lmer.omics, phase == time.phase)
    
    # loop through taxa
    lfc.output = do.call(rbind, lapply(lmer.sig.features, function(feature){
      if(responders == F){
        data.subset = lsarp.lmer.omics[,colnames(lsarp.lmer.omics) %in% c(feature, "standard.name","HM","RS_Name", "Group", "lsarp.days", "baseline")]
      }
      if(responders == T){
        data.subset = lsarp.lmer.omics[,colnames(lsarp.lmer.omics) %in% c(feature, "standard.name","HM","RS_Name", "Group", "lsarp.days", "baseline", "flare.group")]
      }
      
      colnames(data.subset)[2] = "feature"
       # add pseudo (if needed)
      
      # calculate Log2FC
      data.subset = data.subset %>% 
        group_by(HM) %>%
        mutate(lfc = log2(feature)-log2(feature[lsarp.days == min(lsarp.days)])) %>%
        # remove baseline timepoint
        filter(lsarp.days != min(lsarp.days)) %>%
        as.data.frame()
      # add month
      data.subset$month = data.subset$lsarp.days / 30
      data.subset
      # add feature
      data.subset$feature.name = feature
      data.subset
      
    }))
  return(lfc.output)
  }
}



# >> RS Selections --------------------------------------------------------

lsarp.cd.rs.selections = readRDS("./lsarp_rs_scores/2025_06_29_lsarp.cd.rs.selections.RDS")

lsarp.cd.rs.selections$RS_Name = factor(lsarp.cd.rs.selections$Selected, levels=rs.names)

lsarp.cd.rs.selections$HM = lsarp.cd.rs.selections$sample

lsarp.cd.rs.selections = merge(lsarp.cd.rs.selections,
                               distinct(metadata.lsarp.stool[,c("HM", "group")]), by="HM")

lsarp.cd.rs.selections.rs.plot <- ggplot(subset(lsarp.cd.rs.selections,group=="RS"),
                                              aes(x=1, y=Z_score))+
  
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add points
  geom_point(aes(fill=RS), color="white", shape=21, size=2.5, alpha=1)+
  # overlay selected RS (for clarity)
  geom_point(data = subset(lsarp.cd.rs.selections, group=="RS" & RS == Selected), 
             aes(fill=RS), color="white", size=2.5, shape=21, alpha=1)+
  # add lines
  #geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6)+
  # label selected RS
  geom_label(data = subset(lsarp.cd.rs.selections,  group=="RS" & RS == Selected), 
             aes(label=RS_Name, color=RS_Name, x=1, y=3, vjust=0.55),
             size=2.5)+
  scale_y_continuous(limits=c(-2.2,3.2))+
  geom_hline(yintercept=1, linetype=2, color="red", alpha=0.5)+
  scale_color_manual(values = rs.colors)+
  scale_fill_manual(values = rs.colors)+
  theme_minimal()+theme(legend.position="none")+
  labs(x="", y=" ", title="RS")+
  facet_wrap(~HM)+
  theme_minimal()+theme(legend.position="none",
                        axis.text.x=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.x = element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.cd.rs.selections.rs.plot

lsarp.cd.rs.selections.placebo.plot <- ggplot(subset(lsarp.cd.rs.selections,group=="Placebo"),
                                              aes(x=1, y=Z_score))+
  
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add points
  geom_point(aes(fill=RS), color="white", shape=21, size=2.5, alpha=1)+
  # overlay selected RS (for clarity)
  geom_point(data = subset(lsarp.cd.rs.selections, group=="Placebo" & RS == Selected), 
             aes(fill=RS), color="white", size=2.5, shape=21, alpha=1)+
  # add lines
  #geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6)+
  # label selected RS
  geom_label(data = subset(lsarp.cd.rs.selections,  group=="Placebo" & RS == Selected), 
             aes(label=RS_Name, color=RS_Name, x=1, y=3, vjust=0.55),
             size=2.5)+
  scale_y_continuous(limits=c(-2.2,3.2))+
  geom_hline(yintercept=1, linetype=2, color="red", alpha=0.5)+
  scale_color_manual(values = rs.colors)+
  scale_fill_manual(values = rs.colors)+
  theme_minimal()+theme(legend.position="none")+
  labs(x="", y="Butyrogen Z-Score", title="Placebo")+
  facet_wrap(~HM, nrow=4)+
  theme_minimal()+theme(legend.position="none",
                        axis.text.x=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.x = element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.cd.rs.selections.placebo.plot

# break down frequency of RS's being selected, per timepoint (Stacked barplot)
lsarp.cd.rs.selections.frequences = lsarp.cd.rs.selections %>%
  dplyr::select(HM, group, Selected) %>% distinct() %>% as.data.frame() %>%
  dplyr::select(Selected, group) %>%
  table() %>% as.data.frame() %>%
  group_by(group) %>%
  mutate(perc = Freq / sum(Freq)) %>% as.data.frame() 

lsarp.cd.rs.selections.frequences$Selected = factor(lsarp.cd.rs.selections.frequences$Selected, levels=rs.names)

lsarp.cd.rs.selections.frequences.plot = ggplot(lsarp.cd.rs.selections.frequences %>%
                                                  mutate(Group = factor(group, levels=c("RS", "Placebo")))%>%
                                                  rbind(data.frame(Selected = "BobsRedMill", group="RS", Group = "RS", Freq = 0, perc = 0)),
                                                  aes(x=Selected, y=Freq))+
  geom_bar(stat="identity", fill="white", alpha=1)+
  geom_bar(stat="identity", position="stack", 
           aes(fill=Selected), alpha=1,linewidth=0,
           color="white")+
  scale_fill_manual(values= rs.colors)+
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
  facet_wrap(~Group,  nrow=4)
lsarp.cd.rs.selections.frequences.plot

library("patchwork")
lsarp.cd.rs.selections.placebo.plot+
  lsarp.cd.rs.selections.rs.plot+ 
  lsarp.cd.rs.selections.frequences.plot+
  patchwork::plot_layout(widths=c(2,3,2))



# :: Fermentation Responses -----------------------------------------------


lsarp.cd.ph.rs.boxplot = ggplot(metadata.lsarp.stool %>% dplyr::select(HM, Group, rs.col, delta.ph) %>% distinct() %>%
                              # add blanks
                              rbind(data.frame(HM = NA, Group = "RS", rs.col = "FibersymRW", delta.ph=NA),
                                    data.frame(HM = NA, Group = "RS", rs.col = "BobsRedMill", delta.ph=NA)) %>%
                                mutate(Group = factor(Group, levels=c("Placebo", "RS"))) %>%
                                mutate(rs.col = factor(rs.col, levels=rs.names)),
                            aes(x=rs.col, y=delta.ph))+
  geom_boxplot(width=0.5, alpha=0.5, color="black")+
  ggbeeswarm::geom_beeswarm(shape=21, 
             aes(fill=rs.col), color="white", size=3)+
  geom_hline(yintercept=-1.27,
             linetype=2, color="red")+
  #ggrepel::geom_text_repel(aes(label=HM_time),size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(),
                        panel.grid.major.y=element_line(color="grey", linewidth=0.25),
                        axis.text.x=element_text(angle=45, hjust=1))+
  facet_wrap(~Group)+
  labs(y="Δ pH of Selected RS",
       x="")
lsarp.cd.ph.rs.boxplot

# Note: I doubled check these values by going back to original data (from Patient Reports)
subset(metadata.lsarp.stool, HM %in% unique(metadata.lsarp.stool$HM))[,c("HM", "delta.ph", "rs.col")] %>% distinct()
# note: strong Novelose330; pH measurement was performed after freeze-thaw (when device broke)

# fisher test on likelihood of being fermenter/non-fermenter per group:
metadata.lsarp.stool %>% 
  dplyr::select(HM, Group, rs.col, delta.ph) %>% 
  distinct() %>%
  mutate(response = ifelse(delta.ph < -1.27, "strong", "weak")) %>%
  dplyr::select(Group, response) %>%
  table() %>%
  fisher.test()
# threshold must be -1.28 to -1.4 to be sig, not -1.27

# ttest plot
lsarp.ph.ttest = wilcox.test(subset(metadata.lsarp.stool, baseline == "baseline" & group == "RS")$delta.ph,
            subset(metadata.lsarp.stool, baseline == "baseline" & group == "Placebo")$delta.ph)

lsarp.cd.ph.rs.ttest.plot = ggplot(metadata.lsarp.stool %>% 
                                     dplyr::select(HM, Group, rs.col, delta.ph) %>% 
                                     distinct() %>%
                                     mutate(Group = factor(Group, levels=c("Placebo", "RS"))),
       aes(x=Group, y=delta.ph))+
  geom_boxplot(width=0.5, alpha=0.5, color="black")+
  ggbeeswarm::geom_beeswarm(shape=21, 
                            aes(fill=Group), color="white", size=3)+
  geom_hline(yintercept=-1.27,
             linetype=2, color="red")+
  #ggrepel::geom_text_repel(aes(label=HM_time),size=3)+
  #scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.major.y=element_line(color="grey", linewidth=0.25))+
  labs(y="Δ pH of Selected RS",
       x="")
lsarp.cd.ph.rs.ttest.plot

# :: WHY -----------------------------------------------------------

# conduct analyses to understand HOW this happened

metadata.lsarp.stool.ph = subset(metadata.lsarp.stool, timing == "0M")
nrow(metadata.lsarp.stool.ph)
# n = 26


# :: :: baseline microbiome --------------------------------------------------


lsarp.cd.baseline.ph.samples = metadata.lsarp.stool.ph$standard.name
# ASV
lsarp.asv.baseline.ph = lsarp.asv.data.glom[rownames(lsarp.asv.data.glom) %in% lsarp.cd.baseline.ph.samples,] / 50000
lsarp.asv.baseline.ph.pa = lsarp.asv.baseline.ph
lsarp.asv.baseline.ph.pa[lsarp.asv.baseline.ph.pa>0] = 1
lsarp.asv.baseline.ph = lsarp.asv.baseline.ph[,colSums(lsarp.asv.baseline.ph.pa) >= nrow(lsarp.asv.baseline.ph.pa)*0.2]
lsarp.asv.baseline.ph = lsarp.asv.baseline.ph %>% as.data.frame()
# filtered
# lsarp.asv.baseline.ph = log2(lsarp.asv.baseline.ph+min(lsarp.asv.baseline.ph[lsarp.asv.baseline.ph!=0])/2) %>% as.data.frame()
lsarp.asv.baseline.ph$standard.name = rownames(lsarp.asv.baseline.ph)
dim(lsarp.asv.baseline.ph) # 26 samples * 280 features

# calculate Bray-Curtis dissimilarities
lsarp.rapidaim.raw.asv.stool.bray = vegan::vegdist(lsarp.asv.baseline.ph[,colnames(lsarp.asv.baseline.ph) != "standard.name"], method="bray") 
lsarp.rapidaim.raw.asv.stool.pcoa = ape::pcoa(lsarp.rapidaim.raw.asv.stool.bray)
lsarp.rapidaim.raw.asv.stool.bray.df = data.frame(lsarp.rapidaim.raw.asv.stool.pcoa$vectors[,c(1:2)])
lsarp.rapidaim.raw.asv.stool.bray.df$standard.name = rownames(lsarp.rapidaim.raw.asv.stool.bray.df)
lsarp.rapidaim.raw.asv.stool.bray.df = merge(lsarp.rapidaim.raw.asv.stool.bray.df,
                                         metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.rapidaim.raw.asv.stool.bray.var_exp = lsarp.rapidaim.raw.asv.stool.pcoa$values[c(1:2),2]
lsarp.rapidaim.raw.asv.stool.bray.df$var1 = round(lsarp.rapidaim.raw.asv.stool.bray.var_exp[1]*100, digits=2)
lsarp.rapidaim.raw.asv.stool.bray.df$var2 = round(lsarp.rapidaim.raw.asv.stool.bray.var_exp[2]*100, digits=2)
set.seed(25)
lsarp.rapidaim.raw.asv.stool.permanova = vegan::adonis2(lsarp.rapidaim.raw.asv.stool.bray ~ group,
                                                    data = lsarp.rapidaim.raw.asv.stool.bray.df,
                                                    by="margin")
lsarp.rapidaim.raw.asv.stool.permanova # not sig


# MGX
lsarp.mgx.baseline.ph = lsarp.mgx.taxa[rownames(lsarp.mgx.taxa) %in% lsarp.cd.baseline.ph.samples,] / 50000
lsarp.mgx.baseline.ph.pa = lsarp.mgx.baseline.ph
lsarp.mgx.baseline.ph.pa[lsarp.mgx.baseline.ph.pa>0] = 1
lsarp.mgx.baseline.ph = lsarp.mgx.baseline.ph[,colSums(lsarp.mgx.baseline.ph.pa) >= nrow(lsarp.mgx.baseline.ph.pa)*0.2]
lsarp.mgx.baseline.ph = lsarp.mgx.baseline.ph %>% as.data.frame()
# filtered
# lsarp.mgx.baseline.ph = log2(lsarp.mgx.baseline.ph+min(lsarp.mgx.baseline.ph[lsarp.mgx.baseline.ph!=0])/2) %>% as.data.frame()
lsarp.mgx.baseline.ph$standard.name = rownames(lsarp.mgx.baseline.ph)
dim(lsarp.mgx.baseline.ph) # 26 samples * 302 features

# calculate Bray-Curtis dissimilarities
lsarp.rapidaim.raw.mgx.stool.bray = vegan::vegdist(lsarp.mgx.baseline.ph[,colnames(lsarp.mgx.baseline.ph) != "standard.name"], method="bray") 
lsarp.rapidaim.raw.mgx.stool.pcoa = ape::pcoa(lsarp.rapidaim.raw.mgx.stool.bray)
lsarp.rapidaim.raw.mgx.stool.bray.df = data.frame(lsarp.rapidaim.raw.mgx.stool.pcoa$vectors[,c(1:2)])
lsarp.rapidaim.raw.mgx.stool.bray.df$standard.name = rownames(lsarp.rapidaim.raw.mgx.stool.bray.df)
lsarp.rapidaim.raw.mgx.stool.bray.df = merge(lsarp.rapidaim.raw.mgx.stool.bray.df,
                                             metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.rapidaim.raw.mgx.stool.bray.var_exp = lsarp.rapidaim.raw.mgx.stool.pcoa$values[c(1:2),2]
lsarp.rapidaim.raw.mgx.stool.bray.df$var1 = round(lsarp.rapidaim.raw.mgx.stool.bray.var_exp[1]*100, digits=2)
lsarp.rapidaim.raw.mgx.stool.bray.df$var2 = round(lsarp.rapidaim.raw.mgx.stool.bray.var_exp[2]*100, digits=2)
set.seed(25)
lsarp.rapidaim.raw.mgx.stool.permanova = vegan::adonis2(lsarp.rapidaim.raw.mgx.stool.bray ~ group,
                                                        data = lsarp.rapidaim.raw.mgx.stool.bray.df,
                                                        by="margin")
lsarp.rapidaim.raw.mgx.stool.permanova # not sig


# Pathways
lsarp.kegg.baseline.ph = lsarp.cd.mpx.kegg.mat[rownames(lsarp.cd.mpx.kegg.mat) %in% lsarp.cd.baseline.ph.samples,]
lsarp.kegg.baseline.ph[is.na(lsarp.kegg.baseline.ph)] = 0
lsarp.kegg.baseline.ph.pa = lsarp.kegg.baseline.ph
lsarp.kegg.baseline.ph.pa[lsarp.kegg.baseline.ph.pa>0] = 1
lsarp.kegg.baseline.ph = lsarp.kegg.baseline.ph[,colSums(lsarp.kegg.baseline.ph.pa) >= nrow(lsarp.kegg.baseline.ph.pa)*0.2]
lsarp.kegg.baseline.ph = lsarp.kegg.baseline.ph %>% as.data.frame()
# filtered
# lsarp.kegg.baseline.ph = log2(lsarp.kegg.baseline.ph+min(lsarp.kegg.baseline.ph[lsarp.kegg.baseline.ph!=0])/2) %>% as.data.frame()
lsarp.kegg.baseline.ph$standard.name = rownames(lsarp.kegg.baseline.ph)
dim(lsarp.kegg.baseline.ph) # 25 samples * 179 features
# PCA
lsarp.kegg.baseline.ph[is.na(lsarp.kegg.baseline.ph)] = 0
lsarp.baseline.kegg.pca.pseudo = min(lsarp.kegg.baseline.ph[,colnames(lsarp.kegg.baseline.ph) != "standard.name"][lsarp.kegg.baseline.ph[,colnames(lsarp.kegg.baseline.ph) != "standard.name"]!=0])/2
lsarp.baseline.kegg.pca = prcomp(log2(lsarp.kegg.baseline.ph[,colnames(lsarp.kegg.baseline.ph) != "standard.name"]+lsarp.baseline.kegg.pca.pseudo))
lsarp.baseline.kegg.pca.df = data.frame(lsarp.baseline.kegg.pca$x[,c(1:2)])
lsarp.baseline.kegg.pca.df$standard.name = rownames(lsarp.baseline.kegg.pca.df)
lsarp.baseline.kegg.pca.df = merge(lsarp.baseline.kegg.pca.df,
                                   metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.baseline.kegg.pca.var_exp = (lsarp.baseline.kegg.pca$sdev)^2 / sum(lsarp.baseline.kegg.pca$sdev^2) * 100
rownames(lsarp.baseline.kegg.pca.df) = lsarp.baseline.kegg.pca.df$standard.name
set.seed(25)
lsarp.baseline.kegg.pca.permanova = vegan::adonis2(dist(lsarp.baseline.kegg.pca$x) ~ group,
                                               lsarp.baseline.kegg.pca.df,
                                               by="margin")
lsarp.baseline.kegg.pca.permanova

# COGs
lsarp.cog.baseline.ph = lsarp.cd.mpx.cog.mat[rownames(lsarp.cd.mpx.cog.mat) %in% lsarp.cd.baseline.ph.samples,]
lsarp.cog.baseline.ph[is.na(lsarp.cog.baseline.ph)] = 0
lsarp.cog.baseline.ph.pa = lsarp.cog.baseline.ph
lsarp.cog.baseline.ph.pa[lsarp.cog.baseline.ph.pa>0] = 1
lsarp.cog.baseline.ph = lsarp.cog.baseline.ph[,colSums(lsarp.cog.baseline.ph.pa) >= nrow(lsarp.cog.baseline.ph.pa)*0.2]
lsarp.cog.baseline.ph = lsarp.cog.baseline.ph %>% as.data.frame()
# filtered
# lsarp.cog.baseline.ph = log2(lsarp.cog.baseline.ph+min(lsarp.cog.baseline.ph[lsarp.cog.baseline.ph!=0])/2) %>% as.data.frame()
lsarp.cog.baseline.ph$standard.name = rownames(lsarp.cog.baseline.ph)
dim(lsarp.cog.baseline.ph) # 25 samples * 2369 features
# PCA
lsarp.cog.baseline.ph[is.na(lsarp.cog.baseline.ph)] = 0
lsarp.baseline.cog.pca.pseudo = min(lsarp.cog.baseline.ph[,colnames(lsarp.cog.baseline.ph) != "standard.name"][lsarp.cog.baseline.ph[,colnames(lsarp.cog.baseline.ph) != "standard.name"]!=0])/2
lsarp.baseline.cog.pca = prcomp(log2(lsarp.cog.baseline.ph[,colnames(lsarp.cog.baseline.ph) != "standard.name"]+lsarp.baseline.cog.pca.pseudo))
lsarp.baseline.cog.pca.df = data.frame(lsarp.baseline.cog.pca$x[,c(1:2)])
lsarp.baseline.cog.pca.df$standard.name = rownames(lsarp.baseline.cog.pca.df)
lsarp.baseline.cog.pca.df = merge(lsarp.baseline.cog.pca.df,
                                  metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.baseline.cog.pca.var_exp = (lsarp.baseline.cog.pca$sdev)^2 / sum(lsarp.baseline.cog.pca$sdev^2) * 100
rownames(lsarp.baseline.cog.pca.df) = lsarp.baseline.cog.pca.df$standard.name
set.seed(25)
lsarp.baseline.cog.pca.permanova = vegan::adonis2(dist(lsarp.baseline.cog.pca$x) ~ group,
                                                   lsarp.baseline.cog.pca.df,
                                                   by="margin")
lsarp.baseline.cog.pca.permanova
# not sig

# CAZy
lsarp.cazy.baseline.ph = lsarp.cd.mpx.cazy.mat[rownames(lsarp.cd.mpx.cazy.mat) %in% lsarp.cd.baseline.ph.samples,]
lsarp.cazy.baseline.ph[is.na(lsarp.cazy.baseline.ph)] = 0
lsarp.cazy.baseline.ph.pa = lsarp.cazy.baseline.ph
lsarp.cazy.baseline.ph.pa[lsarp.cazy.baseline.ph.pa>0] = 1
lsarp.cazy.baseline.ph = lsarp.cazy.baseline.ph[,colSums(lsarp.cazy.baseline.ph.pa) >= nrow(lsarp.cazy.baseline.ph.pa)*0.2]
lsarp.cazy.baseline.ph = lsarp.cazy.baseline.ph %>% as.data.frame()
# filtered
# lsarp.cazy.baseline.ph = log2(lsarp.cazy.baseline.ph+min(lsarp.cazy.baseline.ph[lsarp.cazy.baseline.ph!=0])/2) %>% as.data.frame()
lsarp.cazy.baseline.ph$standard.name = rownames(lsarp.cazy.baseline.ph)
dim(lsarp.cazy.baseline.ph) # 25 samples * 64 features
# PCA
lsarp.cazy.baseline.ph[is.na(lsarp.cazy.baseline.ph)] = 0
lsarp.baseline.cazy.pca.pseudo = min(lsarp.cazy.baseline.ph[,colnames(lsarp.cazy.baseline.ph) != "standard.name"][lsarp.cazy.baseline.ph[,colnames(lsarp.cazy.baseline.ph) != "standard.name"]!=0])/2
lsarp.baseline.cazy.pca = prcomp(log2(lsarp.cazy.baseline.ph[,colnames(lsarp.cazy.baseline.ph) != "standard.name"]+lsarp.baseline.cazy.pca.pseudo))
lsarp.baseline.cazy.pca.df = data.frame(lsarp.baseline.cazy.pca$x[,c(1:2)])
lsarp.baseline.cazy.pca.df$standard.name = rownames(lsarp.baseline.cazy.pca.df)
lsarp.baseline.cazy.pca.df = merge(lsarp.baseline.cazy.pca.df,
                                   metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.baseline.cazy.pca.var_exp = (lsarp.baseline.cazy.pca$sdev)^2 / sum(lsarp.baseline.cazy.pca$sdev^2) * 100
rownames(lsarp.baseline.cazy.pca.df) = lsarp.baseline.cazy.pca.df$standard.name
set.seed(25)
lsarp.baseline.cazy.pca.permanova = vegan::adonis2(dist(lsarp.baseline.cazy.pca$x) ~ group,
                                                  lsarp.baseline.cazy.pca.df,
                                                  by="margin") 
lsarp.baseline.cazy.pca.permanova
# not sig

# MBX by 80% per response group
lsarp.mbx.baseline.ph = lsarp.mbx.annotated.mat[rownames(lsarp.mbx.annotated.mat) %in% lsarp.cd.baseline.ph.samples,]
lsarp.mbx.baseline.ph[is.na(lsarp.mbx.baseline.ph)] = 0
lsarp.mbx.baseline.ph.pa = lsarp.mbx.baseline.ph
lsarp.mbx.baseline.ph.pa[lsarp.mbx.baseline.ph.pa>0] = 1
lsarp.mbx.baseline.ph = lsarp.mbx.baseline.ph[,colSums(lsarp.mbx.baseline.ph.pa) >= nrow(lsarp.mbx.baseline.ph.pa)*0.8]
lsarp.mbx.baseline.ph = lsarp.mbx.baseline.ph %>% as.data.frame()
# filtered
# lsarp.mbx.baseline = log2(lsarp.mbx.baseline+min(lsarp.mbx.baseline[lsarp.mbx.baseline!=0])/2) %>% as.data.frame()
lsarp.mbx.baseline.ph$standard.name = rownames(lsarp.mbx.baseline.ph)
dim(lsarp.mbx.baseline.ph) # 26 samples * 168 features
# PCA
lsarp.mbx.baseline.ph[is.na(lsarp.mbx.baseline.ph)] = 0
lsarp.baseline.mbx.pca.pseudo = min(lsarp.mbx.baseline.ph[,colnames(lsarp.mbx.baseline.ph) != "standard.name"][lsarp.mbx.baseline.ph[,colnames(lsarp.mbx.baseline.ph) != "standard.name"]!=0])/2
lsarp.baseline.mbx.pca = prcomp(log2(lsarp.mbx.baseline.ph[,colnames(lsarp.mbx.baseline.ph) != "standard.name"]+lsarp.baseline.mbx.pca.pseudo))
lsarp.baseline.mbx.pca.df = data.frame(lsarp.baseline.mbx.pca$x[,c(1:2)])
lsarp.baseline.mbx.pca.df$standard.name = rownames(lsarp.baseline.mbx.pca.df)
lsarp.baseline.mbx.pca.df = merge(lsarp.baseline.mbx.pca.df,
                                  metadata.lsarp.stool.ph[,c("standard.name","group", "delta.ph")], by="standard.name")
lsarp.baseline.mbx.pca.var_exp = (lsarp.baseline.mbx.pca$sdev)^2 / sum(lsarp.baseline.mbx.pca$sdev^2) * 100
rownames(lsarp.baseline.mbx.pca.df) = lsarp.baseline.mbx.pca.df$standard.name
set.seed(25)
lsarp.baseline.mbx.pca.permanova = vegan::adonis2(dist(lsarp.baseline.mbx.pca$x) ~ group,
                                                  lsarp.baseline.mbx.pca.df,
                                                  by="margin")
lsarp.baseline.mbx.pca.permanova
# not sig

# plot

lsarp.baseline.multiome.permanova = 
  rbind(
    as.data.frame(lsarp.rapidaim.raw.asv.stool.permanova)[1,] %>% mutate(type = "ASV"),
    as.data.frame(lsarp.rapidaim.raw.mgx.stool.permanova)[1,] %>% mutate(type = "Species"),
    as.data.frame(lsarp.baseline.kegg.pca.permanova)[1,] %>% mutate(type = "Pathway"),
    as.data.frame(lsarp.baseline.cog.pca.permanova)[1,] %>% mutate(type = "COG"),
    as.data.frame(lsarp.baseline.cazy.pca.permanova)[1,] %>% mutate(type = "CAZy"),
    as.data.frame(lsarp.baseline.mbx.pca.permanova)[1,] %>% mutate(type = "Metabolite"))

ggplot(lsarp.baseline.multiome.permanova %>%
         mutate(type = factor(type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))),
       aes(x=R2, y=reorder(type, R2)))+
  geom_bar(stat="identity", aes(fill=type), linewidth=0.5, color="black")+
  scale_fill_manual(values=omics.colors)+
  theme_classic()+theme(legend.position="none")+
  labs(x="R2", y=NULL)
# none are significant

# :: :: baseline multi-ome -----------------------------------------------

# combine
lsarp.omics.baseline.ph = merge(lsarp.asv.baseline.ph,
                             lsarp.mgx.baseline.ph, by="standard.name")
lsarp.omics.baseline.ph = merge(lsarp.omics.baseline.ph,
                             lsarp.kegg.baseline.ph, by="standard.name")
lsarp.omics.baseline.ph = merge(lsarp.omics.baseline.ph,
                             lsarp.cog.baseline.ph, by="standard.name")
lsarp.omics.baseline.ph = merge(lsarp.omics.baseline.ph,
                             lsarp.cazy.baseline.ph, by="standard.name")
lsarp.omics.baseline.ph = merge(lsarp.omics.baseline.ph,
                             lsarp.mbx.baseline.ph, by="standard.name")

# add ph data
lsarp.omics.baseline.ph = merge(lsarp.omics.baseline.ph,
                                metadata.lsarp.stool.ph[,c("standard.name", "delta.ph")],
                                by="standard.name")
rownames(lsarp.omics.baseline.ph) = lsarp.omics.baseline.ph$standard.name
lsarp.omics.baseline.ph$standard.name = NULL

# RandomForest
set.seed(25)
lsarp.omics.baseline.ph.rf = ranger::ranger(delta.ph ~ .,
               data.frame(lsarp.omics.baseline.ph),
               importance = "permutation")
lsarp.omics.baseline.ph.rf
# no predictive performance
data.frame(imp = lsarp.omics.baseline.ph.rf$variable.importance) %>%
  arrange(-imp) %>% 
  head(n=20)
# mostly nonsense features

# systematic correlation
lsarp.omics.baseline.ph.cor = Hmisc::rcorr(as.matrix(lsarp.omics.baseline.ph), type="spearman")
lsarp.omics.baseline.ph.cor = cbind(reshape2::melt(lsarp.omics.baseline.ph.cor$r),
                                    reshape2::melt(lsarp.omics.baseline.ph.cor$P)[,3]) %>%
  as.data.frame()
colnames(lsarp.omics.baseline.ph.cor) = c("Var1", "Var2", "r", "P")
lsarp.omics.baseline.ph.cor = subset(lsarp.omics.baseline.ph.cor, Var1 == "delta.ph")
lsarp.omics.baseline.ph.cor$padj = p.adjust(lsarp.omics.baseline.ph.cor$P, method="BH")
lsarp.omics.baseline.ph.cor %>% 
  arrange(P) %>%
  head()
# none pass FDR < 0.20


# :: :: baseline 16S cor -----------------------------------------------

# read from lsarp_rapidaim_test
lsarp.rapidaim.asv.mat = readRDS("./lsarp_rapidaim_otu_mat.Rds")

lsarp.rapidaim.asv.stool = lsarp.rapidaim.asv.mat[grepl("Stool", rownames(lsarp.rapidaim.asv.mat)),]
lsarp.rapidaim.asv.stool = lsarp.rapidaim.asv.stool[!grepl("RawStool", rownames(lsarp.rapidaim.asv.stool)),]
# good
lsarp.rapidaim.asv.stool = as.data.frame(lsarp.rapidaim.asv.stool)

lsarp.rapidaim.asv.stool$HM = substr(rownames(lsarp.rapidaim.asv.stool), 1, 6)

lsarp.rapidaim.asv.stool = merge(lsarp.rapidaim.asv.stool,
                                 metadata.lsarp.stool.ph[,c("HM", "delta.ph")])
rownames(lsarp.rapidaim.asv.stool) = lsarp.rapidaim.asv.stool$HM
lsarp.rapidaim.asv.stool$HM = NULL

# RandomForest
set.seed(25)
lsarp.rapidaim.asv.stool.rf = ranger::ranger(delta.ph ~ .,
                                            data.frame(lsarp.rapidaim.asv.stool),
                                            importance = "permutation")
lsarp.rapidaim.asv.stool.rf
# weak predictive performance
data.frame(imp = lsarp.rapidaim.asv.stool.rf$variable.importance) %>%
  arrange(-imp) %>% 
  head(n=20)

# systematic correlation
lsarp.rapidaim.asv.stool.cor = Hmisc::rcorr(as.matrix(lsarp.rapidaim.asv.stool), type="spearman")
lsarp.rapidaim.asv.stool.cor = cbind(reshape2::melt(lsarp.rapidaim.asv.stool.cor$r),
                                    reshape2::melt(lsarp.rapidaim.asv.stool.cor$P)[,3]) %>%
  as.data.frame()
colnames(lsarp.rapidaim.asv.stool.cor) = c("Var1", "Var2", "r", "P")
lsarp.rapidaim.asv.stool.cor = subset(lsarp.rapidaim.asv.stool.cor, Var1 == "delta.ph")
lsarp.rapidaim.asv.stool.cor$padj = p.adjust(lsarp.rapidaim.asv.stool.cor$P, method="BH")
lsarp.rapidaim.asv.stool.cor %>% 
  arrange(padj) %>%
  head()
# F. prausnitzii and Gemellaceae pass FDR < 0.20


# :: :: loop through RS ---------------------------------------------------

# need to regenerate source data!

# :: :: F. prausnitzii ----------------------------------------------------

lsarp.rapidaim.asv.stool.fprau.plot = ggplot(lsarp.rapidaim.asv.stool %>%
         mutate(HM = rownames(.)) %>%
         merge(metadata.lsarp.stool.ph[,c("HM","group")],
                      by="HM"),
       aes(x=Faecalibacterium_prausnitzii*100, y=delta.ph))+
  ggpubr::stat_cor(method="spearman",
                   label.y.npc = "bottom",label.x.npc = "left")+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  #scale_x_log10()+
  theme_classic()+
  labs(x="Baseline F. prausnitzii (%)",
       y="Δ pH",
       fill="Group")
lsarp.rapidaim.asv.stool.fprau.plot

# Lachnospira and F. prausnitzii predict delta-pH (sig)
# F. prausnitzii is strongly and sig correlated

# are they different between groups?
lsarp.rapidaim.asv.stool.meta = lsarp.rapidaim.asv.stool %>%
  mutate(HM = rownames(.)) %>%
  merge(metadata.lsarp.stool.ph[,c("HM","group")],
        by="HM")
  
ggplot(lsarp.rapidaim.asv.stool.meta,
       aes(x=group, y=Faecalibacterium_prausnitzii))+
  geom_boxplot()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL, fill="Group", y="Baseline F. prausnitzii (%)")


# :: :: baseline bray -----------------------------------------------------

# calculate Bray-Curtis dissimilarities
lsarp.rapidaim.asv.stool.bray = vegan::vegdist(lsarp.rapidaim.asv.stool[,colnames(lsarp.rapidaim.asv.stool) != "delta.ph"], method="bray") 
# perform PCoA
lsarp.rapidaim.asv.stool.pcoa = ape::pcoa(lsarp.rapidaim.asv.stool.bray)
# extract data from pcoa
lsarp.rapidaim.asv.stool.pcoa.df = data.frame(lsarp.rapidaim.asv.stool.pcoa$vectors[,c(1:2)])
lsarp.rapidaim.asv.stool.pcoa.df$HM = rownames(lsarp.rapidaim.asv.stool.pcoa.df)
# add metadata
lsarp.rapidaim.asv.stool.pcoa.df = merge(lsarp.rapidaim.asv.stool.pcoa.df,
                                         metadata.lsarp.stool.ph[,c("HM","group", "delta.ph")], by="HM")
# extract variance explained
lsarp.rapidaim.asv.stool.pcoa.var_exp = lsarp.rapidaim.asv.stool.pcoa$values[c(1:2),2]
lsarp.rapidaim.asv.stool.pcoa.df$var1 = round(lsarp.rapidaim.asv.stool.pcoa.var_exp[1]*100, digits=2)
lsarp.rapidaim.asv.stool.pcoa.df$var2 = round(lsarp.rapidaim.asv.stool.pcoa.var_exp[2]*100, digits=2)

# clean up "lsarp.on.rs" variable

set.seed(25)
t1 = Sys.time()
lsarp.rapidaim.asv.stool.permanova = vegan::adonis2(lsarp.rapidaim.asv.stool.bray ~ delta.ph,
                                     data = lsarp.rapidaim.asv.stool.pcoa.df,
                                     by="margin")
t2 = Sys.time()
t2 - t1
lsarp.rapidaim.asv.stool.permanova # not sig

lsarp.rapidaim.asv.stool.pcoa.plot <- ggplot(
  data=lsarp.rapidaim.asv.stool.pcoa.df %>%
    mutate(med1 = median(Axis.1),
           med2 = median(Axis.2)),
  aes(x=Axis.1, y=Axis.2))+
  geom_segment(aes(x=med1, xend=Axis.1,
                   y=med2, yend=Axis.2), linewidth=0.2, linetype=2)+
  geom_point(shape=21, aes(fill=scale(delta.ph)), color="black", size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  #geom_text(aes(label=standard.name), size=3, vjust=1.5)+
  labs(x=paste("Axis 1: ", round(unique(lsarp.rapidaim.asv.stool.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.rapidaim.asv.stool.pcoa.df$var2), digits=2), "%", sep=""),
       title=paste(paste("Δ pH R²: ", round(data.frame(lsarp.rapidaim.asv.stool.permanova)[1,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.rapidaim.asv.stool.permanova)[1,5], 3), sep="")),
       fill="Scaled Δ pH") + 
  theme_classic()+theme(#legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.rapidaim.asv.stool.pcoa.plot

# does it correlate with distance from median?
lsarp.rapidaim.asv.stool.pcoa.df = lsarp.rapidaim.asv.stool.pcoa.df %>%
  mutate(dis1 = Axis.1 - median(Axis.1),
         dis2 = Axis.2 - median(Axis.2)) %>%
  mutate(dis = sqrt(dis1^2 + dis2^2))

lsarp.rapidaim.asv.stool.pcoa.dis.plot = ggplot(lsarp.rapidaim.asv.stool.pcoa.df,
       aes(x=dis, y = delta.ph))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group),color="white", size=3)+
  theme_classic()+#theme(legend.position="none")+
  labs(x="Distance to Medioid",
       y="Δ pH",
       fill="Group")  
lsarp.rapidaim.asv.stool.pcoa.dis.plot

lsarp.rapidaim.asv.stool.dis.ttest = wilcox.test(subset(lsarp.rapidaim.asv.stool.pcoa.df, group == "Placebo")$dis,
            subset(lsarp.rapidaim.asv.stool.pcoa.df, group == "RS")$dis)
# p = 0.085
lsarp.rapidaim.asv.stool.pcoa.dis.ttest.plot = ggplot(lsarp.rapidaim.asv.stool.pcoa.df,
       aes(x=group, y = dis))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  annotate(geom="text", x=1.5, y=Inf, vjust=2,
            label=paste0("p: ", round(lsarp.rapidaim.asv.stool.dis.ttest$p.value, digits=2), sep=""))+
  labs(x=NULL,
       y="Distance to Medioid",
       fill="Group")  
lsarp.rapidaim.asv.stool.pcoa.dis.ttest.plot

# correlate with F. prausnitzii
lsarp.rapidaim.asv.stool.pcoa.df = merge(lsarp.rapidaim.asv.stool.pcoa.df,
                                         lsarp.rapidaim.asv.stool.meta[,c("HM", "Faecalibacterium_prausnitzii")], by="HM")

lsarp.rapidaim.asv.stool.pcoa.dis.fp.plot = ggplot(lsarp.rapidaim.asv.stool.pcoa.df,
                                                aes(x=dis, y = Faecalibacterium_prausnitzii))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",
                   label.y.npc = "top",
                   label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+#theme(legend.position="none")+
  labs(x="Distance to Medioid",
       y="Baseline F. prausnitzii (%)",
       fill="Group")  
lsarp.rapidaim.asv.stool.pcoa.dis.fp.plot
# strongly correlate


lsarp.rapidaim.asv.stool.fprau.plot+
lsarp.rapidaim.asv.stool.pcoa.plot+
lsarp.rapidaim.asv.stool.pcoa.dis.plot+
lsarp.rapidaim.asv.stool.pcoa.dis.ttest.plot


# :: :: RS vs Placebo --------------------------------------------------------

set.seed(25)
lsarp.rapidaim.asv.stool.group.permanova = vegan::adonis2(lsarp.rapidaim.asv.stool.bray ~ group,
               data = lsarp.rapidaim.asv.stool.pcoa.df,
               by="margin")

lsarp.rapidaim.asv.stool.group.pcoa.plot <- ggplot(
  data=lsarp.rapidaim.asv.stool.pcoa.df %>%
    mutate(med1 = median(Axis.1),
           med2 = median(Axis.2)),
  aes(x=Axis.1, y=Axis.2))+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  stat_ellipse(aes(color=group))+
  #scale_fill_gradient2(low="blue", high="red")+
  #geom_text(aes(label=standard.name), size=3, vjust=1.5)+
  labs(x=paste("Axis 1: ", round(unique(lsarp.rapidaim.asv.stool.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.rapidaim.asv.stool.pcoa.df$var2), digits=2), "%", sep=""),
       title=paste(paste("Group R²: ", round(data.frame(lsarp.rapidaim.asv.stool.group.permanova)[1,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.rapidaim.asv.stool.group.permanova)[1,5], 3), sep="")),
       fill="Group", color="Group") + 
  theme_classic()+theme(#legend.position="none",
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))
lsarp.rapidaim.asv.stool.group.pcoa.plot


# maaslin2
rownames(lsarp.rapidaim.asv.stool.pcoa.df) = lsarp.rapidaim.asv.stool.pcoa.df$HM
lsarp.rapidaim.asv.stool.group.maaslin = Maaslin2::Maaslin2(input_data = lsarp.rapidaim.asv.stool[,colnames(lsarp.rapidaim.asv.stool) != "delta.ph"],
                                                    input_metadata = lsarp.rapidaim.asv.stool.pcoa.df,
                                                    output = "~/Downloads",
                                                    fixed_effects = c("group"),  # Example fixed effects
                                                    #random_effects = c("HM"),       # Example random effects
                                                    normalization = "TSS",                       # Total Sum Scaling normalization
                                                    transform = "LOG",                           # Log transformation
                                                    analysis_method = "LM",                      # Linear model
                                                    plot_scatter = FALSE,                        # Disable scatterplot generation
                                                    plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                    max_significance = 0.05,                     # Significance threshold for q-values
                                                    standardize = TRUE                           # Disable standardization (optional)
)
lsarp.rapidaim.asv.stool.group.maaslin$results %>% arrange(pval)

# are RS vs Placebo different?
set.seed(25)
ranger::ranger(as.factor(group) ~ .,
               data.frame(lsarp.rapidaim.asv.stool.meta)[,colnames(lsarp.rapidaim.asv.stool.meta) != "delta.ph"],
               importance = "permutation")
ranger::ranger(as.factor(group) ~ .,
               data.frame(lsarp.rapidaim.asv.stool.meta)[,colnames(lsarp.rapidaim.asv.stool.meta) != "delta.ph"],
               importance = "permutation")$variable.importance %>% sort() %>% tail()
# not predictable

# Prevotella_copri
ggplot(lsarp.rapidaim.asv.stool.meta,
       aes(x=group, y=Prevotella_copri))+
  scale_y_log10()+
  geom_boxplot()+
  geom_point(aes(fill=group), color="white", shape=21, size=3)+
  theme_classic()+
  labs(fill="Group", x= NULL)
# not sig




# :: :: Technical ---------------------------------------------------------

lsarp.rapidaim.asv.stool[,c("delta.ph")]
# randomization scheme
# culture batch
# culture technician
# stool processing notes

lsarp.rapidaims.data = read.csv("./2025_12_08_lsarp_cd_rapidaims_data.csv")

# append 16S data
lsarp.rapidaims.data = merge(lsarp.rapidaims.data %>% dplyr::select(-group),
                             lsarp.rapidaim.asv.stool.pcoa.df, by="HM")

# Differences in Time:
lsarp.rapidaims.ttest.process = wilcox.test(subset(lsarp.rapidaims.data, group=="RS")$days_to_process,
            subset(lsarp.rapidaims.data, group=="Placebo")$days_to_process)
# p = 0.127
ggplot(lsarp.rapidaims.data,
       aes(x=group, y=days_to_process))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0.1,
              aes(fill=group), color="white", shape=21, size=3)+
  annotate(geom="text", x=1.5, y=Inf, vjust=2,
           label=paste0("p: ", round(lsarp.rapidaims.ttest.process$p.value, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL,
       y="Days to process")

lsarp.rapidaims.ttest.culture = wilcox.test(subset(lsarp.rapidaims.data, group=="RS")$time_to_culture,
            subset(lsarp.rapidaims.data, group=="Placebo")$time_to_culture)
# p = 0.017
ggplot(lsarp.rapidaims.data,
       aes(x=group, y=time_to_culture))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0.1,
              aes(fill=group), color="white", shape=21, size=3)+
  annotate(geom="text", x=1.5, y=Inf, vjust=2,
           label=paste0("p: ", round(lsarp.rapidaims.ttest.culture$p.value, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL,
       y="Days to culture (- processing)")

lsarp.rapidaims.ttest.both = wilcox.test(subset(lsarp.rapidaims.data, group=="RS")$days_to_process+subset(lsarp.rapidaims.data, group=="RS")$time_to_culture,
            subset(lsarp.rapidaims.data, group=="Placebo")$days_to_process+subset(lsarp.rapidaims.data, group=="Placebo")$time_to_culture)
# p = 0.009
lsarp.rapidaims.data.plot.1 = ggplot(lsarp.rapidaims.data,
       aes(x=group, y=days_to_process+time_to_culture))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0.1,
              aes(fill=group), color="white", shape=21, size=3)+
  annotate(geom="text", x=1.5, y=Inf, vjust=2,
           label=paste0("p: ", round(lsarp.rapidaims.ttest.both$p.value, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL,
   y="Days to culture (+ processing)")
lsarp.rapidaims.data.plot.1

set.seed(25)
lsarp.rapidaim.asv.stool.group.permanova = vegan::adonis2(lsarp.rapidaim.asv.stool.bray ~ total,
                                                          data = lsarp.rapidaim.asv.stool.pcoa.df %>%
                                                            merge(lsarp.rapidaims.data[,c("HM", "days_to_process", "time_to_culture")]%>%
                                                                    mutate(total = days_to_process + time_to_culture)) %>%
                                                            mutate(total = ifelse(is.na(total), mean(na.omit(total)), total)),
                                                          by="margin")
# not sig

# random forest

# Correlate with fermentation:
ggplot(lsarp.rapidaims.data,
       aes(x=days_to_process, y = delta.ph.y))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days to process",
       y="Δ pH",
       fill = "Group")

ggplot(lsarp.rapidaims.data,
       aes(x=time_to_culture, y = delta.ph.y))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days to culture (after processing)",
      y="Δ pH",
      fill = "Group")

lsarp.rapidaims.data.plot.2 = ggplot(lsarp.rapidaims.data,
       aes(x=time_to_culture+days_to_process, y = delta.ph.y))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days to culture (+ processing)",
       y="Δ pH",
       fill = "Group")
lsarp.rapidaims.data.plot.2

ggplot(lsarp.rapidaims.data,
       aes(x=time_to_culture+days_to_process, y = Faecalibacterium_prausnitzii))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days to culture (+ processing)",
       y="Baseline F. prausnitzii (%)",
       fill = "Group")

ggplot(lsarp.rapidaims.data,
       aes(x=time_to_culture+days_to_process, y = dis))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days to culture (+ processing)",
       y="Distance to Medioid",
       fill = "Group")


# random forest
lsarp.rapidaim.asv.stool.time = lsarp.rapidaim.asv.stool %>%
  mutate(HM = rownames(.)) %>%
  merge(lsarp.rapidaims.data[,c("HM", "days_to_process", "time_to_culture")] %>%
          mutate(total = days_to_process+time_to_culture) %>%
          subset(!is.na(total)), by="HM") %>%
  dplyr::select(-delta.ph, -days_to_process, -time_to_culture)
# RandomForest
set.seed(25)
lsarp.rapidaim.asv.stool.time.rf = ranger::ranger(total ~ .,
                                            data.frame(lsarp.rapidaim.asv.stool.time),
                                            importance = "permutation")
lsarp.rapidaim.asv.stool.time.rf
# no predictive performance
data.frame(imp = lsarp.rapidaim.asv.stool.time.rf$variable.importance) %>%
  arrange(-imp) %>% 
  head(n=20)
# Prevotella_copri
lsarp.rapidaim.asv.stool.time.pcopri = ggplot(lsarp.rapidaim.asv.stool.time %>%
         merge(lsarp.rapidaims.data[,c("HM","days_to_process", "group")], by="HM"),
       aes(x=`Prevotella_copri`, y=total))+
  geom_smooth(method="lm", color="black")+
  #scale_x_log10()+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Baseline P. copri (%)",
       y="Days to culture (+ processing)",
       fill="Group")  
lsarp.rapidaim.asv.stool.time.pcopri

# [Eubacterium]_dolichum (not sig when log10 is not used)
ggplot(lsarp.rapidaim.asv.stool.time %>%
         merge(lsarp.rapidaims.data[,c("HM","days_to_process", "group")], by="HM"),
       aes(x=`[Eubacterium]_dolichum`, y=total))+
  geom_smooth(method="lm", color="black")+
  #scale_x_log10()+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Baseline E. dolichum (%)",
       y="Days to culture (+ processing)",
       fill="Group")  

# Order received:
lsarp.rapidaims.data.order.plot.1 = ggplot(lsarp.rapidaims.data %>%
         mutate(rank.hm = as.numeric(rownames(.))),
       aes(x=rank.hm, y=delta.ph.y))+
  geom_line(linetype=2, linewidth=0.2)+
  geom_smooth(method="lm", color="black")+
  #scale_x_log10()+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(y="Δ pH",
       x="Order received") 
lsarp.rapidaims.data.order.plot.1

lsarp.rapidaims.data = lsarp.rapidaims.data %>%
  mutate(rank.hm = as.numeric(rownames(.)))

lsarp.rapidaims.data.order.ttest = wilcox.test(subset(lsarp.rapidaims.data, group == "RS")$rank.hm,
                                               subset(lsarp.rapidaims.data, group == "Placebo")$rank.hm)
lsarp.rapidaims.data.order.plot.2 = ggplot(lsarp.rapidaims.data,
       aes(x=group, y=rank.hm))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0.1,
              shape=21, aes(fill=group), color="white", size=3)+
  annotate(geom="text", x=1.5, y=Inf, vjust=2,
           label=paste0("p: ", round(lsarp.rapidaims.data.order.ttest$p.value, digits=3), sep=""))+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL,
       y="Order received") 
lsarp.rapidaims.data.order.plot.2

# Order received:
lsarp.rapidaims.data.order.plot.3 = ggplot(lsarp.rapidaims.data %>%
                                             mutate(rank.hm = as.numeric(rownames(.)),
                                                    total = days_to_process+time_to_culture),
                                           aes(x=rank.hm, y=total))+
  geom_line(linetype=2, linewidth=0.2)+
  geom_smooth(method="lm", color="black")+
  #scale_x_log10()+
  ggpubr::stat_cor(method="spearman",label.y.npc = "bottom",label.x.npc = "middle")+
  geom_point(shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(y="Days to culture (+ processing)",
       x="Order received") 
lsarp.rapidaims.data.order.plot.3

# Order received:
lsarp.rapidaims.data.order.plot.4 = ggplot(lsarp.rapidaims.data %>%
                                             mutate(rank.hm = as.numeric(rownames(.)),
                                                    total = days_to_process+time_to_culture),
                                           aes(x=rank.hm, y=total))+
  geom_tile(aes(fill=group, x=rank.hm, y=group), color="black")+
  theme_classic()+theme(legend.position="none")+
  coord_flip()+
  labs(y=NULL,
       x="Order received") 
lsarp.rapidaims.data.order.plot.4


# Person processing
lsarp.rapidaims.data.personelle.ph.plot = ggplot(lsarp.rapidaims.data,
       aes(x=person_processing, y=delta.ph.y))+
  geom_boxplot(outlier.shape=NA, width=0.5)+
  geom_jitter(width=0.1, height=0.1,
             shape=21, aes(fill=group), color="white", size=3)+
  theme_classic()+
  labs(x="Personelle Processing",
       y="Δ pH",
       fill="Group")
lsarp.rapidaims.data.personelle.ph.plot
kruskal.test(delta.ph.y ~ person_processing, lsarp.rapidaims.data)
# kruskal-wallis test is not sig
lm(delta.ph.y ~ person_processing, lsarp.rapidaims.data) %>% summary()
# lm is not sig

lsarp.rapidaims.data.personelle.process.plot = ggplot(lsarp.rapidaims.data,
       aes(x=person_processing))+
  geom_bar(position="stack", stat="count",
           aes(fill=group), color="white")+
  theme_classic()+
  labs(x="Personel Processing",
       y="Count",
       fill="Group")
lsarp.rapidaims.data.personelle.process.plot
# no difference

# :: Plots ----------------------------------------------------------------

lsarp.cd.ph.rs.boxplot+
  lsarp.cd.ph.rs.ttest.plot+
  patchwork::plot_layout(nrow=1, widths=c(2,1))
# difference in fermentation b/w groups

lsarp.rapidaim.asv.stool.group.pcoa.plot
# sig difference in baseline microbiome b/w groups


lsarp.rapidaim.asv.stool.fprau.plot+
  lsarp.rapidaim.asv.stool.pcoa.plot+
  lsarp.rapidaim.asv.stool.pcoa.dis.plot+
  lsarp.rapidaim.asv.stool.pcoa.dis.ttest.plot
# fermentation correlates w baseline F. prausnitzii + distance from medioid


lsarp.rapidaims.data.order.plot.4+
  lsarp.rapidaims.data.order.plot.2+
lsarp.rapidaims.data.order.plot.1+
  lsarp.rapidaims.data.order.plot.3+
  patchwork::plot_layout(nrow=1, widths=c(1,1,2,2))
# no association w order of collection

lsarp.rapidaims.data.personelle.process.plot+
  lsarp.rapidaims.data.personelle.ph.plot

lsarp.rapidaims.data.plot.1+
lsarp.rapidaims.data.plot.2+
  lsarp.rapidaim.asv.stool.time.pcopri+
  patchwork::plot_layout(nrow=1, widths=c(1,2,2))
# sig difference in time to culture b/w groups (SOMEHOW)




# >>> 1. STANDARD ANALYSES -----------------------------------

# :: 6M Survival Curve --------------------------------------------------

# add flare day
redcap.lsarp.wpcdai = read.csv("./2024_07_31_unblinding_DM.csv")   
# ensure we only use compliant
redcap.lsarp.wpcdai$HM = substr(redcap.lsarp.wpcdai$study_id, 1, 6)
redcap.lsarp.wpcdai = subset(redcap.lsarp.wpcdai, HM %in% metadata.lsarp.stool$HM)
nrow(redcap.lsarp.wpcdai)
# n = 26, good
# clean dataframe
redcap.lsarp.wpcdai = redcap.lsarp.wpcdai[,c("HM", "group","days.taken", "flare.day", "flare.call","flare.call_verified", "scope.day")]
# add status "event" (1) if flare occurs before scope
redcap.lsarp.wpcdai$status.pre = ifelse(redcap.lsarp.wpcdai$flare.call == "ifx.by.visit2", 1,0)
# add status event (1) if flare occurs at any time in study
redcap.lsarp.wpcdai$status.all = ifelse(is.na(redcap.lsarp.wpcdai$flare.day), 0, 1)

subset(redcap.lsarp.wpcdai, group == "RS") %>% arrange(flare.day)

# add ave.fiber
redcap.lsarp.wpcdai = merge(redcap.lsarp.wpcdai,
                            metadata.lsarp.stool[,c("HM", "ave.fiber")]%>%distinct(), by="HM")

# make time-to-event the last observation before re-scope
redcap.lsarp.wpcdai$time.to.event = ifelse(!is.na(redcap.lsarp.wpcdai$scope.day),
                                           redcap.lsarp.wpcdai$scope.day, 
                                           redcap.lsarp.wpcdai$flare.day)

redcap.lsarp.wpcdai[,c("HM", "group", "days.taken", "flare.day", "scope.day", "time.to.event")]

# numbers
redcap.lsarp.wpcdai[,c("HM", "group")] %>% distinct() %>% dplyr::select(group) %>% table()

# first question is: "is remission more likely in patients on RS up to scope"
redcap.lsarp.survival <- survival::survfit(Surv(time.to.event, status.pre) ~ group, 
                                 data = redcap.lsarp.wpcdai)
redcap.lsarp.log.rank.test <- survival::survdiff(Surv(time.to.event, status.pre) ~ group, 
                       data = redcap.lsarp.wpcdai)
redcap.lsarp.log.rank.test
# not significant

redcap.lsarp.survival.plot = survminer::ggsurvplot(redcap.lsarp.survival, 
           data = redcap.lsarp.wpcdai,
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
# Right-censored at 2nd scope, flare rate

redcap.lsarp.survival.plot
# cross hatch = Standard of care scope date

# Completed 12 m without flare
# RS: 938, 960, 1035
# Plac: 978, 1022

# 878? flared near very end of 12M

# use fisher test including endoscope (64% vs 83%)
data.frame(
  #condition = c("placebo", "rs"),
  total = c(12, 14),
  flare = c(10, 9)) %>%
  fisher.test()

# use fisher test up to endoscope ( 79% vs 58%)
data.frame(
  #condition = c("placebo", "rs"),
  total = c(12, 14),
  flare = c(7, 11)) %>%
  fisher.test()
# not sig

# :: Fecal Calprotectin ----------------------------------------------------------

# numbers
metadata.lsarp.stool[,c("HM", "group", "flare")] %>% subset(flare != "flare" & !HM %in% excluded.by.dave) %>% table()
metadata.lsarp.stool[,c("HM", "group", "flare")] %>% subset(flare != "flare" & !HM %in% excluded.by.dave) %>% distinct() %>% dplyr::select(group) %>% table()

# does fecal calprotectin increase or decrease over RS treatment
# stats
lsarp.stats.fcal.treatment = lmerTest::lmer(scale(log10(fcal)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                subset(metadata.lsarp.stool, 
                                       phase %in% c("treatment"))%>%
                                  # remove samples collected after flare
                                  subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.fcal.washout = lmerTest::lmer(scale(log10(fcal)) ~ group*scale(lsarp.days) + ave.fiber +  (1|HM),
                                          subset(metadata.lsarp.stool, 
                                                 phase == "washout" | rs.end == "rs.end") %>%
                                            # remove samples collected after flare
                                            subset(flare != "flare") %>%
                                            # for washout
                                            subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.fcal.plot = ggplot(metadata.lsarp.stool %>%
                                    subset(!is.na(fcal))%>%
                                    # remove samples collected after flare
                                    subset(flare != "flare"& !HM %in% excluded.by.dave),
                                 aes(x=lsarp.days, y=fcal))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, 
                                           phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_alpha_manual(values=c(0.5, 1))+
  geom_hline(yintercept=250, color="red")+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  labs(x="Days since starting Product", y="Fecal Calprotectin (μg/g)",
       title=paste(paste("Treatment p:", round(lsarp.stats.fcal.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.fcal.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.fcal.plot


# :: Stool Water ----------------------------------------------------------

# Note: change to "moisture" because dry weight excludes more than water (i.e. other volatiles)
# does stool water increase or decrease over RS treatment

# stats
lsarp.stats.water.treatment = lmerTest::lmer(scale((stool_water_perc)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                             subset(metadata.lsarp.stool, 
                                                    phase %in% c("treatment"))%>%
                                               # remove samples collected after flare
                                               subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.water.washout = lmerTest::lmer(scale((stool_water_perc)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                           subset(metadata.lsarp.stool, phase %in% "washout" | rs.end == "rs.end")%>%
                                             # remove samples collected after flare
                                             subset(flare != "flare") %>%
                                             # for washout
                                             subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.water.plot = ggplot(metadata.lsarp.stool %>%
                                     subset(!is.na(stool_water_perc))%>%
                                     # remove samples collected after flare
                                     subset(flare != "flare"& !HM %in% excluded.by.dave),
                                  aes(x=lsarp.days, 
                                      y=stool_water_perc*100))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  labs(x="Days since starting Product", y="Stool Moisture (%)",
       title=paste(paste("Treatment p:", round(lsarp.stats.water.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.water.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.water.plot


# :: ASV Microbial Load -------------------------------------------------------

# stats
lsarp.stats.load.treatment = lmerTest::lmer(scale(load.asv) ~ group*scale(lsarp.days)+ ave.fiber + (1|HM),
                                subset(metadata.lsarp.stool.asv, phase %in% "treatment") %>%
                                  # remove samples collected after flare
                                  subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.load.washout = lmerTest::lmer(scale(load.asv) ~ group*scale(lsarp.days)+ ave.fiber + (1|HM),
                                 subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end") %>%
                                   # remove samples collected after flare
                                   subset(flare != "flare") %>%
                                   # for washout
                                   subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.load.plot = ggplot(metadata.lsarp.stool.asv %>%
                                    subset(!is.na(load.asv))%>%
                                    # remove samples collected after flare
                                    subset(flare != "flare"& !HM %in% excluded.by.dave),
                                 aes(x=lsarp.days,
                                     y=load.asv))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white",size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  labs(x="Days since starting Product", y="Microbial Load",
       title=paste(paste("Treatment p:", round(lsarp.stats.load.treatment[5,5], 3)),
                   paste("  Washout p:", round(lsarp.stats.load.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.load.plot


# :: ASV Richness ------------------------------------------------------

# stats
lsarp.stats.richness.treatment = lmerTest::lmer(scale(log10(richness)) ~ group*scale(lsarp.days)+ ave.fiber + (1|HM),
                                    subset(metadata.lsarp.stool.asv, phase == "treatment") %>%
                                      # remove samples collected after flare
                                      subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.richness.washout = lmerTest::lmer(scale(log10(richness)) ~ group*scale(lsarp.days)+ave.fiber +  (1|HM),
                                     subset(metadata.lsarp.stool.asv,  phase == "washout" | rs.end == "rs.end") %>%
                                       # remove samples collected after flare
                                       subset(flare != "flare") %>%
                                       # for washout
                                       subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.richness.plot = ggplot(metadata.lsarp.stool.asv %>%
                                                  subset(!is.na(richness))%>%
                                                  # remove samples collected after flare
                                                  subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=richness))+
 # scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="ASV Richness",
       title=paste(paste("Treatment p:", round(lsarp.stats.richness.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.richness.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.richness.plot


# :: ASV Shannon ------------------------------------------------------


# stats
lsarp.stats.shannon.treatment = lmerTest::lmer(scale((shannon)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                                 # remove samples collected after flare
                                                 subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.shannon.washout = lmerTest::lmer(scale((shannon)) ~ group*scale(lsarp.days)+  ave.fiber + (1|HM),
                                              subset(metadata.lsarp.stool.asv,  phase == "washout" | rs.end == "rs.end") %>%
                                               # remove samples collected after flare
                                               subset(flare != "flare") %>%
                                               # for washout
                                               subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.shannon.plot = ggplot(metadata.lsarp.stool.asv  %>%
                                                 subset(!is.na(shannon))%>%
                                                 # remove samples collected after flare
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                                aes(x=lsarp.days, 
                                                    y=shannon))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  labs(x="Days since starting Product", y="Shannon Diversity",
       title=paste(paste("Treatment p:", round(lsarp.stats.shannon.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.shannon.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.shannon.plot


# :: ASV Functional Redundancy ------------------------------------------------------

# check Functional Richness vs Taxa Richness
ggplot(metadata.lsarp.stool.asv,
       aes(x=fd, y=richness))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  geom_point(fill="black", color="white", shape=21, size=3)+
  theme_classic()+
  labs(x="Functional Richness (Number of unique ECs)",
       y="Taxa Richness (Number of unique ASVs)")
# Functional Richness correlates with Taxa richness (R = 0.46)

# check FR vs FD
ggplot(metadata.lsarp.stool.asv,
       aes(x=fd, y=fr))+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  geom_point(fill="black", color="white", shape=21, size=3)+
  theme_classic()+
  labs(x="Functional Richness",
       y="Functional Redundancy")
# Functional Richness correlates with Functional Redundancy (R = 0.40)

# stats
lsarp.stats.fr.treatment = lmerTest::lmer(scale((fr)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                               subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                            # remove samples collected after flare
                                            subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.fr.washout = lmerTest::lmer(scale((fr)) ~ group*scale(lsarp.days) +  ave.fiber + (1|HM),
                                             subset(metadata.lsarp.stool.asv,  phase == "washout" | rs.end == "rs.end")  %>%
                                          # remove samples collected after flare
                                          subset(flare != "flare") %>%
                                          # for washout
                                          subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.fr.plot = ggplot(metadata.lsarp.stool.asv %>%
                                            subset(!is.na(fr))%>%
                                            # remove samples collected after flare
                                            subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=fr))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Functional Redundancy",
       title=paste(paste("Treatment p:", round(lsarp.stats.fr.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.fr.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.fr.plot


# :: ASV Functional Richness ------------------------------------------------------

# stats
lsarp.stats.fd.treatment = lmerTest::lmer(scale((fd)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                          subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                            # remove samples collected after flare
                                            subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.fd.washout = lmerTest::lmer(scale((fd)) ~ group*scale(lsarp.days) +  ave.fiber + (1|HM),
                                        subset(metadata.lsarp.stool.asv,  phase == "washout" | rs.end == "rs.end")  %>%
                                          # remove samples collected after flare
                                          subset(flare != "flare") %>%
                                          # for washout
                                          subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.fd.plot = ggplot(metadata.lsarp.stool.asv %>%
                                            subset(!is.na(fd))%>%
                                            # remove samples collected after flare
                                            subset(flare != "flare"& !HM %in% excluded.by.dave),
                                          aes(x=lsarp.days, 
                                              y=fd))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Functional Richness",
       title=paste(paste("Treatment p:", round(lsarp.stats.fd.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.fd.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.fd.plot


# :: ASV Butyrogens I -------------------------------------------------------

# stats
lsarp.stats.but.i.treatment = lmerTest::lmer(scale(log10(but.i)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                             subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                               # remove samples collected after flare
                                               subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.but.i.washout = lmerTest::lmer(scale(log10(but.i)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                           subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end")  %>%
                                             # remove samples collected after flare
                                             subset(flare != "flare") %>%
                                             # for washout
                                             subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.but.i.plot = ggplot(metadata.lsarp.stool.asv %>%
                                     subset(!is.na(but.i))%>%
                                     # remove samples collected after flare
                                     subset(flare != "flare"& !HM %in% excluded.by.dave),
                                   aes(x=lsarp.days, 
                                       y=but.i*100))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Butyrogens (%)",
       title=paste(paste("Treatment p:", round(lsarp.stats.but.i.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.but.i.washout[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.but.i.plot


# :: ASV Butyrogens II (Vital) -------------------------------------------------------

# stats
lsarp.stats.but.ii.treatment = lmerTest::lmer(scale(log10(but.ii)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                              subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                                # remove samples collected after flare
                                                subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.but.ii.washout = lmerTest::lmer(scale(log10(but.ii)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                            subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end")  %>%
                                              # remove samples collected after flare
                                              subset(flare != "flare") %>%
                                              # for washout
                                              subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.but.ii.plot = ggplot(metadata.lsarp.stool.asv %>%
                                      subset(!is.na(but.ii))%>%
                                      # remove samples collected after flare
                                      subset(flare != "flare"& !HM %in% excluded.by.dave),
                                    aes(x=lsarp.days, 
                                        y=but.ii*100))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value= "grey")+
  labs(x="Days since starting Product", y="Kircher Butyrogens (%)",
       title=paste(paste("Treatment p:", round(lsarp.stats.but.ii.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.but.ii.washout[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.but.ii.plot


# :: ASV Bray-Curtis ------------------------------------------------------

# note: HM0860-STL-07 ONLY contains ASVs not detected in >10% of other samples
# therefore, use glommed taxa to enable proper comparisons

# apply 20% prevalence filter (first time doing this in this analysis)
# keep important samples (not flare and not excluded)
lsarp.asv.data.glom.filt = lsarp.asv.data.glom[rownames(lsarp.asv.data.glom) %in% 
                             subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,]

dim(lsarp.asv.data.glom.filt) # 1554 glommed taxa
lsarp.asv.data.median.pa = lsarp.asv.data.glom.filt
lsarp.asv.data.median.pa[lsarp.asv.data.median.pa > 0] = 1
# 10%
lsarp.asv.data.median.filt = lsarp.asv.data.glom.filt[,colSums(lsarp.asv.data.median.pa) >= nrow(lsarp.asv.data.median.pa)*0.1]
dim(lsarp.asv.data.median.filt) # 431 taxa
# 20%
lsarp.asv.data.median.filt = lsarp.asv.data.glom.filt[,colSums(lsarp.asv.data.median.pa) >= nrow(lsarp.asv.data.median.pa)*0.2]
ncol(lsarp.asv.data.median.filt) # 313 taxa


# calculate Bray-Curtis dissimilarities
lsarp.asv.bray = vegan::vegdist(lsarp.asv.data.median.filt, method="bray") 
# perform PCoA
lsarp.asv.pcoa = ape::pcoa(lsarp.asv.bray)
# extract data from pcoa
lsarp.asv.pcoa.df = data.frame(lsarp.asv.pcoa$vectors[,c(1:2)])
lsarp.asv.pcoa.df$standard.name = rownames(lsarp.asv.pcoa.df)
# add metadata
lsarp.asv.pcoa.df = merge(lsarp.asv.pcoa.df,
                         metadata.lsarp.stool.asv, by="standard.name")
# extract variance explained
lsarp.asv.pcoa.var_exp = lsarp.asv.pcoa$values[c(1:2),2]
lsarp.asv.pcoa.df$var1 = round(lsarp.asv.pcoa.var_exp[1]*100, digits=2)
lsarp.asv.pcoa.df$var2 = round(lsarp.asv.pcoa.var_exp[2]*100, digits=2)

# clean up "lsarp.on.rs" variable

set.seed(25)
t1 = Sys.time()
lsarp.asv.permanova = vegan::adonis2(lsarp.asv.bray ~ group*on.rs + ave.fiber,
                                    lsarp.asv.pcoa.df %>% mutate(on.rs = ifelse(lsarp.on.rs == "on.rs", "on.rs", "off.rs")),
                                    strata = lsarp.asv.pcoa.df$HM,
                                    by="margin")
t2 = Sys.time()
t2 - t1
lsarp.asv.permanova # sig

lsarp.asv.pcoa.plot = ggplot(
  data=lsarp.asv.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), alpha=0.5, linetype=2, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value= "grey")+
  #geom_text(aes(label=standard.name), size=3, vjust=1.5)+
  labs(x=paste("Axis 1: ", round(unique(lsarp.asv.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.asv.pcoa.df$var2), digits=2), "%", sep=""),
       title=paste(paste("ASV\nTreatment R²: ", round(data.frame(lsarp.asv.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.asv.permanova)[2,5], 3), sep="")), sep="") + 
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.asv.pcoa.plot



# :: ASV Beta Diversity ---------------------------------------------------

beta.trajectory.data = beta.trajectory(lsarp.asv.bray)

beta.trajectory.data =   merge(metadata.lsarp.stool.asv, 
                               beta.trajectory.data, by="standard.name")

# stats
stats.beta.between.treatment = lmerTest::lmer(scale((between.beta)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                               subset(beta.trajectory.data, phase == "treatment")  %>%
                                                 # remove samples collected after flare
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

stats.beta.between.washout = lmerTest::lmer(scale((between.beta)) ~ group*scale(lsarp.days)+  ave.fiber + (1|HM),
                                             subset(beta.trajectory.data,  phase == "washout" | rs.end == "rs.end") %>%
                                              # remove samples collected after flare
                                              subset(flare != "flare") %>%
                                              # for washout
                                              subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.between.plot = ggplot(beta.trajectory.data  %>%
                                                 # remove samples collected after flare
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=between.beta))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Bray-Curtis Dissimilarity",
       title=paste(paste("Treatment p:", round(stats.beta.between.treatment[5,5], 3)), 
                   paste("  Washout p:", round(stats.beta.between.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.between.plot



# stats
stats.beta.within.treatment = lmerTest::lmer(scale((within.beta)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                              subset(beta.trajectory.data, phase == "treatment")  %>%
                                                # remove samples collected after flare
                                                subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

stats.beta.within.washout = lmerTest::lmer(scale((within.beta)) ~ group*scale(lsarp.days)+  ave.fiber + (1|HM),
                                            subset(beta.trajectory.data,  phase == "washout" | rs.end == "rs.end") %>%
                                             # remove samples collected after flare
                                             subset(flare != "flare") %>%
                                             # for washout
                                             subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.within.plot = ggplot(beta.trajectory.data  %>%
                                                 # remove samples collected after flare
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=within.beta))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  labs(x="Days since starting Product", y="Bray-Curtis Dissimilarity",
       title=paste(paste("Treatment p:", round(stats.beta.within.treatment[5,5], 3)), 
                   paste("  Washout p:", round(stats.beta.within.washout[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.within.plot


# :: ASV Interaction ------------------------------------------------------

# these are the stools intended for group level analyses; same as for the individual microbiome variables:

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()
# n = 138 total eligible
lsarp.lmer.asv.data = lsarp.delta.omic.prepare(lsarp.asv.data.glom[rownames(lsarp.asv.data.glom) %in% stools.not.flare,],
                                               normalize = T)
nrow(lsarp.lmer.asv.data)
# n = 135 with available data

lsarp.lmer.asv.interactions = lsarp.delta.omic.lmer(lsarp.lmer.asv.data,
                                                    split=F)
# fix names
lsarp.lmer.asv.interactions$taxa =  ifelse(grepl("s__", (lsarp.lmer.asv.interactions$taxa)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.lmer.asv.interactions$taxa)))), sep=""),
                                 paste("(", substr((lsarp.lmer.asv.interactions$taxa), 1, 1), ") ", gsub("c__", "", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.lmer.asv.interactions$taxa))))), sep=""))

lsarp.lmer.asv.interactions.plot = ggplot(lsarp.lmer.asv.interactions %>%
                                            mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                          aes(x=estimate, y=padj))+
  geom_point(shape=21, color="black", aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Group Interaction (ASV)",
       y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.asv.interactions.plot
lsarp.lmer.asv.interactions
# no padj < 0.20 during treatment

subset(lsarp.lmer.asv.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.asv.interactions, padj < 0.20 & phase == "washout")


# >>> 2. MULTI-OMIC ANALYSES -------------------------------------------------------------

# :: MGX Bray-Curtis ------------------------------------------------------

# apply 20% prevalence filter
dim(lsarp.mgx.taxa) # 1224 taxa
lsarp.mgx.taxa.pa = lsarp.mgx.taxa
lsarp.mgx.taxa.pa[lsarp.mgx.taxa.pa > 0] = 1
# 10%
lsarp.mgx.taxa.filt = lsarp.mgx.taxa[,colSums(lsarp.mgx.taxa.pa) >= nrow(lsarp.mgx.taxa.pa)*0.1]
dim(lsarp.mgx.taxa.filt) # 475 taxa
# 20%
lsarp.mgx.taxa.filt = lsarp.mgx.taxa[,colSums(lsarp.mgx.taxa.pa) >= nrow(lsarp.mgx.taxa.pa)*0.2]
ncol(lsarp.mgx.taxa.filt) # 311 taxa


# calculate Bray-Curtis dissimilarities
lsarp.mgx.bray = vegan::vegdist(lsarp.mgx.taxa.filt[rownames(lsarp.mgx.taxa.filt) %in% 
                                                        subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,], method="bray") 
# perform PCoA
lsarp.mgx.pcoa = ape::pcoa(lsarp.mgx.bray)
# extract data from pcoa
lsarp.mgx.pcoa.df = data.frame(lsarp.mgx.pcoa$vectors[,c(1:2)])
lsarp.mgx.pcoa.df$standard.name = rownames(lsarp.mgx.pcoa.df)
# add metadata
lsarp.mgx.pcoa.df = merge(lsarp.mgx.pcoa.df,
                          metadata.lsarp.stool.asv, by="standard.name")
# extract variance explained
lsarp.mgx.pcoa.var_exp = lsarp.mgx.pcoa$values[c(1:2),2]
lsarp.mgx.pcoa.df$var1 = round(lsarp.mgx.pcoa.var_exp[1]*100, digits=2)
lsarp.mgx.pcoa.df$var2 = round(lsarp.mgx.pcoa.var_exp[2]*100, digits=2)

# clean up "lsarp.on.rs" variable

set.seed(25)
t1 = Sys.time()
lsarp.mgx.permanova = vegan::adonis2(lsarp.mgx.bray ~ group*on.rs + ave.fiber,
                                     lsarp.mgx.pcoa.df %>% mutate(on.rs = ifelse(lsarp.on.rs == "on.rs", "on.rs", "off.rs")),
                                     strata = lsarp.mgx.pcoa.df$HM,
                                     by="margin")
t2 = Sys.time()
t2 - t1
lsarp.mgx.permanova # sig

lsarp.mgx.pcoa.plot <- ggplot(
  data=lsarp.mgx.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), alpha=0.5, linetype=2, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  #geom_text(aes(label=HM))+
  labs(x=paste("Axis 1: ", round(unique(lsarp.mgx.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.mgx.pcoa.df$var2), digits=2), "%", sep=""),
       title=paste(paste("Species\nTreatment R²: ", round(data.frame(lsarp.mgx.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.mgx.permanova)[2,5], 3), sep="")), sep="") + 
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.mgx.pcoa.plot



# :: MGX Interaction ------------------------------------------------------

# Goal: Perform linear regression to test interaction between Group*DaysonRS
# this preserves power and shouldn't suffer as much from imbalanced groups

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()
# n = 138 total eligible

lsarp.lmer.mgx.data = lsarp.delta.omic.prepare(lsarp.mgx.taxa[rownames(lsarp.mgx.taxa) %in% stools.not.flare,],
                                               normalize=T)
nrow(lsarp.lmer.mgx.data)
# n = 128 with available data

lsarp.lmer.mgx.interactions = lsarp.delta.omic.lmer(lsarp.lmer.mgx.data,
                                                    split=F)

lsarp.lmer.mgx.interactions.plot = ggplot(lsarp.lmer.mgx.interactions%>%
                                            mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                          aes(x=estimate, y=padj))+
                                          geom_point(shape=21, aes(fill=estimate))+
                                            geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
                                            scale_fill_gradient2(low="blue", high="red")+
                                            scale_y_continuous(transform=neg_log10_trans,
                                                               breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
                                            ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                                                                     size=2.5)+
                                            ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                                                                     size=2.5, color="grey")+
                                            labs(x="Interaction Coefficient",
                                                 title="Group Interaction (Species)",
                                                 y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))

lsarp.lmer.mgx.interactions.plot

lsarp.lmer.mgx.interactions
# Washout + RS seems to increase many important taxa; perhaps because of treatment

# Treatment decreases B. adolescentis, increases 3 other taxa

subset(lsarp.lmer.mgx.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.mgx.interactions, padj < 0.20 & phase == "washout")


# :: Pathway PCA ------------------------------------------------------

# fix rownames
rownames(lsarp.cd.mpx.kegg.mat) = gsub("\\_", "\\-", rownames(lsarp.cd.mpx.kegg.mat))
# remove outlier
lsarp.kegg.lsarp.pca = lsarp.cd.mpx.kegg.mat[rownames(lsarp.cd.mpx.kegg.mat) != "HM0999-STL-06",]
# subset to no-flare
lsarp.kegg.lsarp.pca = lsarp.kegg.lsarp.pca[rownames(lsarp.kegg.lsarp.pca) %in% subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,]

lsarp.kegg.lsarp.pca = lsarp.kegg.lsarp.pca[,colnames(lsarp.kegg.lsarp.pca) != "X."]
lsarp.kegg.lsarp.pca[is.na(lsarp.kegg.lsarp.pca)] = 0

# apply 20% prevalence filter
dim(lsarp.kegg.lsarp.pca) # 181 Pathways (84 samples)
lsarp.kegg.lsarp.pca.pa = lsarp.kegg.lsarp.pca
lsarp.kegg.lsarp.pca.pa[lsarp.kegg.lsarp.pca.pa > 0] = 1
# 10%
lsarp.kegg.lsarp.pca.filt = lsarp.kegg.lsarp.pca[,colSums(lsarp.kegg.lsarp.pca.pa) >= nrow(lsarp.kegg.lsarp.pca.pa)*0.1]
dim(lsarp.kegg.lsarp.pca.filt) # 179 Pathways
# 20%
lsarp.kegg.lsarp.pca.filt = lsarp.kegg.lsarp.pca[,colSums(lsarp.kegg.lsarp.pca.pa) >= nrow(lsarp.kegg.lsarp.pca.pa)*0.2]
dim(lsarp.kegg.lsarp.pca.filt) # 177 Pathways

# log transform
lsarp.kegg.lsarp.pca.filt = log2(lsarp.kegg.lsarp.pca.filt+(min(lsarp.kegg.lsarp.pca.filt[lsarp.kegg.lsarp.pca.filt!=0])/2))
lsarp.kegg.lsarp.pca = prcomp((lsarp.kegg.lsarp.pca.filt), scale=T)
lsarp.kegg.lsarp.pca.df = lsarp.kegg.lsarp.pca$x[,c(1,2)] %>% as.data.frame() %>%
  rownames_to_column("standard.name")
lsarp.kegg.lsarp.pca.df = merge(lsarp.kegg.lsarp.pca.df,
                               metadata.lsarp.stool, by="standard.name")
lsarp.kegg.lsarp.pca.var <- (lsarp.kegg.lsarp.pca$sdev)^2 / sum(lsarp.kegg.lsarp.pca$sdev^2) * 100

rownames(lsarp.cd.mpx.kegg.mat) %in% lsarp.kegg.lsarp.pca.df$standard.name

# permanova (euclidean distance of PC1 and 2)
set.seed(25)
t1 = Sys.time()
lsarp.kegg.lsarp.pca.permanova = vegan::adonis2(dist(lsarp.kegg.lsarp.pca$x) ~ group*on.rs + ave.fiber,
                                               lsarp.kegg.lsarp.pca.df %>% mutate(on.rs = as.factor(ifelse(lsarp.on.rs=="on.rs", "on.rs", "off.rs"))),
                                               strata = lsarp.kegg.lsarp.pca.df$HM,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.kegg.lsarp.pca.permanova
# sig

lsarp.kegg.lsarp.pca.plot <- ggplot(
  data=lsarp.kegg.lsarp.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2) %>%
    mutate(Group = group), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "grey")+
  facet_wrap(~Group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round((lsarp.kegg.lsarp.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round((lsarp.kegg.lsarp.pca.var)[2], digits=2), "%", sep=""),
       title=paste(paste("Pathway\nTreatment R²: ", round(data.frame(lsarp.kegg.lsarp.pca.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.kegg.lsarp.pca.permanova)[2,5], 3), sep="")), sep="")
lsarp.kegg.lsarp.pca.plot

# outlier: HM0999-STL-06

# :: Pathway Interaction -----------------------------------------------

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()

lsarp.lmer.mpx.kegg.data = lsarp.delta.omic.prepare(lsarp.cd.mpx.kegg.mat[rownames(lsarp.cd.mpx.kegg.mat) %in% stools.not.flare,],
                                                    normalize=F)

lsarp.lmer.mpx.kegg.data = lsarp.lmer.mpx.kegg.data[,colnames(lsarp.lmer.mpx.kegg.data) != "X."]

lsarp.lmer.mpx.kegg.interactions = lsarp.delta.omic.lmer(lsarp.lmer.mpx.kegg.data,
                                                         split=F)

# fix feature names
lsarp.kegg.feature.map = data.frame(clean.name = colnames(lsarp.cd.mpx.kegg.mat),
                                    taxa = make.names(colnames(lsarp.cd.mpx.kegg.mat)))
lsarp.lmer.mpx.kegg.interactions$clean.name = lsarp.kegg.feature.map$clean.name[match(lsarp.lmer.mpx.kegg.interactions$taxa,
                                                                                      lsarp.kegg.feature.map$taxa)]

lsarp.lmer.mpx.kegg.interactions.plot = ggplot(lsarp.lmer.mpx.kegg.interactions%>%
                                                 mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, taxa, NA)),
                           size=2.5)+
  ggnetwork::geom_nodetext_repel(aes(label=ifelse(pval < 0.05 & padj > 0.2, taxa, NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Group Interaction (Pathway)",
       y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.kegg.interactions.plot
lsarp.lmer.mpx.kegg.interactions

# nothing sig during treatment
subset(lsarp.lmer.mpx.kegg.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.mpx.kegg.interactions, padj < 0.20 & phase == "washout")


# :: COG PCA ------------------------------------------------------

# fix rownames
rownames(lsarp.cd.mpx.cog.mat) = gsub("\\_", "\\-", rownames(lsarp.cd.mpx.cog.mat))
# remove outlier
lsarp.cog.lsarp.pca = lsarp.cd.mpx.cog.mat[rownames(lsarp.cd.mpx.cog.mat) != "HM0999-STL-06",]
# subset to no-flare
lsarp.cog.lsarp.pca = lsarp.cog.lsarp.pca[rownames(lsarp.cog.lsarp.pca) %in% subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,]

lsarp.cog.lsarp.pca = lsarp.cog.lsarp.pca[,colnames(lsarp.cog.lsarp.pca) != "X."]
lsarp.cog.lsarp.pca[is.na(lsarp.cog.lsarp.pca)] = 0

# apply 10% prevalence filter
dim(lsarp.cog.lsarp.pca) # 2696 COGs (122 samples)
lsarp.cog.lsarp.pca.pa = lsarp.cog.lsarp.pca
lsarp.cog.lsarp.pca.pa[lsarp.cog.lsarp.pca.pa > 0] = 1
# 10%
lsarp.cog.lsarp.pca.filt = lsarp.cog.lsarp.pca[,colSums(lsarp.cog.lsarp.pca.pa) >= nrow(lsarp.cog.lsarp.pca.pa)*0.1]
dim(lsarp.cog.lsarp.pca.filt) # 2491 COGs
# 20%
lsarp.cog.lsarp.pca.filt = lsarp.cog.lsarp.pca[,colSums(lsarp.cog.lsarp.pca.pa) >= nrow(lsarp.cog.lsarp.pca.pa)*0.2]
dim(lsarp.cog.lsarp.pca.filt) # 2342 COGs

# log transform
lsarp.cog.lsarp.pca.filt = log2(lsarp.cog.lsarp.pca.filt+(min(lsarp.cog.lsarp.pca.filt[lsarp.cog.lsarp.pca.filt!=0])/2))

lsarp.cog.lsarp.pca = prcomp((lsarp.cog.lsarp.pca.filt), scale=T)
lsarp.cog.lsarp.pca.df = lsarp.cog.lsarp.pca$x[,c(1,2)] %>% as.data.frame() %>%
  rownames_to_column("standard.name")
lsarp.cog.lsarp.pca.df = merge(lsarp.cog.lsarp.pca.df,
                                metadata.lsarp.stool, by="standard.name")
lsarp.cog.lsarp.pca.var <- (lsarp.cog.lsarp.pca$sdev)^2 / sum(lsarp.cog.lsarp.pca$sdev^2) * 100

rownames(lsarp.cd.mpx.cog.mat) %in% lsarp.cog.lsarp.pca.df$standard.name
# many samples excluded because of after flare

# permanova
set.seed(25)
t1 = Sys.time()
lsarp.cog.lsarp.pca.permanova = vegan::adonis2(dist(lsarp.cog.lsarp.pca$x) ~ group*on.rs +ave.fiber,
                                                lsarp.cog.lsarp.pca.df %>% mutate(on.rs = as.factor(ifelse(lsarp.on.rs=="on.rs", "on.rs", "off.rs"))),
                                                strata = lsarp.cog.lsarp.pca.df$HM,
                                                by="margin")
t2 = Sys.time()
t2 - t1
lsarp.cog.lsarp.pca.permanova
# sig

lsarp.cog.lsarp.pca.plot <- ggplot(
  data=lsarp.cog.lsarp.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2) %>%
    mutate(Group = group), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  facet_wrap(~Group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round((lsarp.cog.lsarp.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round((lsarp.cog.lsarp.pca.var)[2], digits=2), "%", sep=""),
       title=paste(paste("COG\nTreatment R²: ", round(data.frame(lsarp.cog.lsarp.pca.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.cog.lsarp.pca.permanova)[2,5], 3), sep="")), sep="")
lsarp.cog.lsarp.pca.plot

# :: COG Interaction -----------------------------------------------

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()

lsarp.lmer.mpx.cog.data = lsarp.delta.omic.prepare(lsarp.cd.mpx.cog.mat[rownames(lsarp.cd.mpx.cog.mat) %in% stools.not.flare,],
                                                    normalize=F)
lsarp.lmer.mpx.cog.data = lsarp.lmer.mpx.cog.data[,colnames(lsarp.lmer.mpx.cog.data) != "X."]

lsarp.lmer.mpx.cog.interactions = lsarp.delta.omic.lmer(lsarp.lmer.mpx.cog.data,
                                                         split=F)

# fix feature names
lsarp.cog.feature.map = data.frame(clean.name = colnames(lsarp.cd.mpx.cog.mat),
                                    taxa = make.names(colnames(lsarp.cd.mpx.cog.mat)))
lsarp.lmer.mpx.cog.interactions$clean.name = lsarp.cog.feature.map$clean.name[match(lsarp.lmer.mpx.cog.interactions$taxa,
                                                                                      lsarp.cog.feature.map$taxa)]

lsarp.lmer.mpx.cog.interactions.plot = ggplot(lsarp.lmer.mpx.cog.interactions%>%
                                                 mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, taxa, NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Group Interaction (COG)",
       y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.cog.interactions.plot
lsarp.lmer.mpx.cog.interactions

# nothing sig in Treatment
subset(lsarp.lmer.mpx.cog.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.mpx.cog.interactions, padj < 0.20 & phase == "washout")

# :: CAZy PCA ------------------------------------------------------

# fix rownames
rownames(lsarp.cd.mpx.cazy.mat) = gsub("\\_", "\\-", rownames(lsarp.cd.mpx.cazy.mat))
# remove outlier
lsarp.cazy.lsarp.pca = lsarp.cd.mpx.cazy.mat[rownames(lsarp.cd.mpx.cazy.mat) != "HM0999-STL-06",]
# subset to no-flare
lsarp.cazy.lsarp.pca = lsarp.cazy.lsarp.pca[rownames(lsarp.cazy.lsarp.pca) %in% subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,]

lsarp.cazy.lsarp.pca = lsarp.cazy.lsarp.pca[,colnames(lsarp.cazy.lsarp.pca) != "X."]
lsarp.cazy.lsarp.pca[is.na(lsarp.cazy.lsarp.pca)] = 0

# apply 20% prevalence filter
dim(lsarp.cazy.lsarp.pca) # 66 CAZy (125 samples)
lsarp.cazy.lsarp.pca.pa = lsarp.cazy.lsarp.pca
lsarp.cazy.lsarp.pca.pa[lsarp.cazy.lsarp.pca.pa > 0] = 1
# 10%
lsarp.cazy.lsarp.pca.filt = lsarp.cazy.lsarp.pca[,colSums(lsarp.cazy.lsarp.pca.pa) >= nrow(lsarp.cazy.lsarp.pca.pa)*0.1]
dim(lsarp.cazy.lsarp.pca.filt) # 65 CAZy
# 20%
lsarp.cazy.lsarp.pca.filt = lsarp.cazy.lsarp.pca[,colSums(lsarp.cazy.lsarp.pca.pa) >= nrow(lsarp.cazy.lsarp.pca.pa)*0.2]
dim(lsarp.cazy.lsarp.pca.filt) # 63 CAZy

# log transform
lsarp.cazy.lsarp.pca.filt = log2(lsarp.cazy.lsarp.pca.filt+(min(lsarp.cazy.lsarp.pca.filt[lsarp.cazy.lsarp.pca.filt!=0])/2))
lsarp.cazy.lsarp.pca = prcomp((lsarp.cazy.lsarp.pca.filt), scale=T)
lsarp.cazy.lsarp.pca.df = lsarp.cazy.lsarp.pca$x[,c(1,2)] %>% as.data.frame() %>%
  rownames_to_column("standard.name")
lsarp.cazy.lsarp.pca.df = merge(lsarp.cazy.lsarp.pca.df,
                               metadata.lsarp.stool, by="standard.name")
lsarp.cazy.lsarp.pca.var <- (lsarp.cazy.lsarp.pca$sdev)^2 / sum(lsarp.cazy.lsarp.pca$sdev^2) * 100

rownames(lsarp.cd.mpx.cazy.mat) %in% lsarp.cazy.lsarp.pca.df$standard.name

# permanova
set.seed(25)
t1 = Sys.time()
lsarp.cazy.lsarp.pca.permanova = vegan::adonis2(dist(lsarp.cazy.lsarp.pca$x) ~ group*on.rs +ave.fiber,
                                               lsarp.cazy.lsarp.pca.df %>% mutate(on.rs = as.factor(ifelse(lsarp.on.rs=="on.rs", "on.rs", "off.rs"))),
                                               strata = lsarp.cazy.lsarp.pca.df$HM,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.cazy.lsarp.pca.permanova
# not sig

lsarp.cazy.lsarp.pca.plot <- ggplot(
  data=lsarp.cazy.lsarp.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2) %>%
    mutate(Group = group), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline),color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value= "grey")+
  facet_wrap(~Group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round((lsarp.cazy.lsarp.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round((lsarp.cazy.lsarp.pca.var)[2], digits=2), "%", sep=""),
       title=paste(paste("CAZy\nTreatment R²: ", round(data.frame(lsarp.cazy.lsarp.pca.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.cazy.lsarp.pca.permanova)[2,5], 3), sep="")), sep="")
lsarp.cazy.lsarp.pca.plot

# outlier: HM0999-STL-06

# numbers
lsarp.cazy.lsarp.pca.df[,c("HM", "Group")] %>% table()
lsarp.cazy.lsarp.pca.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()

# :: CAZy Interaction -----------------------------------------------

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()

lsarp.lmer.mpx.cazy.data = lsarp.delta.omic.prepare(lsarp.cd.mpx.cazy.mat[rownames(lsarp.cd.mpx.cazy.mat) %in% stools.not.flare,],
                                                    normalize=F)
lsarp.lmer.mpx.cazy.data = lsarp.lmer.mpx.cazy.data[,colnames(lsarp.lmer.mpx.cazy.data) != "X."]

lsarp.lmer.mpx.cazy.interactions = lsarp.delta.omic.lmer(lsarp.lmer.mpx.cazy.data,
                                                         split=F)

lsarp.lmer.mpx.cazy.interactions.plot = ggplot(lsarp.lmer.mpx.cazy.interactions%>%
                                                 mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Group Interaction (CAZy)",
       y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.cazy.interactions.plot
lsarp.lmer.mpx.cazy.interactions

# nothing is padj < 0.20 for treatment
subset(lsarp.lmer.mpx.cazy.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.mpx.cazy.interactions, padj < 0.20 & phase == "washout")

# :: CAZy Starch ----------------------------------------------------

# stats
lsarp.stats.starch.treatment = lmerTest::lmer(scale(log10(starch)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                    subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                                      # remove samples collected after flare
                                                      subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.starch.washout = lmerTest::lmer(scale(log10(starch)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                  subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end")  %>%
                                                    # remove samples collected after flare
                                                    subset(flare != "flare") %>%
                                                    # for washout
                                                    subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.starch.plot = ggplot(metadata.lsarp.stool.asv %>%
                                      subset(!is.na(starch))%>%
                                            # remove samples collected after flare
                                            subset(flare != "flare"& !HM %in% excluded.by.dave),
                                          aes(x=lsarp.days, 
                                              y=starch))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  #geom_text(aes(label=HM), size=3, vjust=-1)+
  labs(x="Days since starting Product", y="Starch CAZy Intensity",
       title=paste(paste("Treatment p:", round(lsarp.stats.starch.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.starch.washout[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.starch.plot


# :: CAZy Mucin ----------------------------------------------------

# stats
lsarp.stats.mucin.treatment = lmerTest::lmer(scale(log10(mucin)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                    subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                                      # remove samples collected after flare
                                                      subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.mucin.washout = lmerTest::lmer(scale(log10(mucin)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                                  subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end")  %>%
                                                    # remove samples collected after flare
                                                    subset(flare != "flare") %>%
                                                    # for washout
                                                    subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.mucin.plot = ggplot(metadata.lsarp.stool.asv %>%
                                     subset(!is.na(starch))%>%
                                            # remove samples collected after flare
                                            subset(flare != "flare"& !HM %in% excluded.by.dave),
                                          aes(x=lsarp.days, 
                                              y=mucin))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  # geom_text(aes(label=HM), size=3, vjust=-1)+
  labs(x="Days since starting Product", y="Mucin CAZy Intensity",
       title=paste(paste("Treatment p:", round(lsarp.stats.mucin.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.mucin.washout[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.mucin.plot

# :: CAZy Starch:Mucin ----------------------------------------------------

# stats
lsarp.stats.starch.mucin.treatment = lmerTest::lmer(scale((starch.mucin)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                             subset(metadata.lsarp.stool.asv, phase == "treatment")  %>%
                                               # remove samples collected after flare
                                               subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.stats.starch.mucin.washout = lmerTest::lmer(scale((starch.mucin)) ~ group*scale(lsarp.days) + ave.fiber + (1|HM),
                                           subset(metadata.lsarp.stool.asv, phase == "washout" | rs.end == "rs.end")  %>%
                                             # remove samples collected after flare
                                             subset(flare != "flare") %>%
                                             # for washout
                                             subset(!HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.starch.mucin.plot = ggplot(metadata.lsarp.stool.asv %>%
                                            subset(!is.na(starch))%>%
                                     # remove samples collected after flare
                                     subset(flare != "flare"& !HM %in% excluded.by.dave),
                                   aes(x=lsarp.days, 
                                       y=starch.mucin))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white", span=1)+
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  labs(x="Days since starting Product", y="Starch:Mucin Ratio",
       title=paste(paste("Treatment p:", round(lsarp.stats.starch.mucin.treatment[5,5], 3)), 
                   paste("  Washout p:", round(lsarp.stats.starch.mucin.washout[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.starch.mucin.plot



# :: MBX Filter -----------------------------------------------------------

# lsarp.mbx.raw.mat = unannotated --> no sig results
# lsarp.mbx.annotated.mat = only annotated

# 80% prevalence; PER TIME PHASE (preRS, onRS, postRS)
#dim(lsarp.mbx.raw.mat)
lsarp.mbx.raw.mat.filt.1 = lsarp.mbx.annotated.mat %>% as.matrix() %>%
  reshape2::melt()
colnames(lsarp.mbx.raw.mat.filt.1)[1] = "standard.name"
lsarp.mbx.raw.mat.filt.1 = merge(lsarp.mbx.raw.mat.filt.1,
                                metadata.lsarp.stool[,c("standard.name", "group")], by="standard.name")
lsarp.mbx.raw.mat.filt.1$value = ifelse(lsarp.mbx.raw.mat.filt.1$value >0, 1, 0)
# calc prevalence per group
lsarp.mbx.raw.mat.filt.1 = lsarp.mbx.raw.mat.filt.1 %>%
  group_by(Var2, group) %>%
  summarise(prevalence = mean(value, na.rm = TRUE)) %>% # cool, use mean of presence to calculate prevalence
  ungroup() %>%
  subset(prevalence >= .80)
# keep features with 80% prevalence in at least one of the groups
length(unique(lsarp.mbx.raw.mat.filt.1$Var2))
# n = 231 features (unannotated)
# n = 149 features (annotated)
# apply filter
lsarp.mbx.raw.mat.filt.1 = lsarp.mbx.annotated.mat[,colnames(lsarp.mbx.annotated.mat) %in% unique(lsarp.mbx.raw.mat.filt.1$Var2)]
dim(lsarp.mbx.raw.mat.filt.1)



# :: MBX PCA --------------------------------------------------------------


# subset to no-flare
lsarp.mbx.lsarp.pca = lsarp.mbx.raw.mat.filt.1[rownames(lsarp.mbx.raw.mat.filt.1) %in% subset(metadata.lsarp.stool.asv, flare != "flare" & !HM %in% excluded.by.dave)$standard.name,]

# already 80% filtered

# log transform
lsarp.mbx.lsarp.pca = log2(lsarp.mbx.lsarp.pca+(min(lsarp.mbx.lsarp.pca[lsarp.mbx.lsarp.pca!=0])/2))

lsarp.mbx.lsarp.pca = prcomp((lsarp.mbx.lsarp.pca), scale=T)
lsarp.mbx.lsarp.pca.df = lsarp.mbx.lsarp.pca$x[,c(1,2)] %>% as.data.frame() %>%
  rownames_to_column("standard.name")
lsarp.mbx.lsarp.pca.df = merge(lsarp.mbx.lsarp.pca.df,
                               metadata.lsarp.stool, by="standard.name")
lsarp.mbx.lsarp.pca.var <- (lsarp.mbx.lsarp.pca$sdev)^2 / sum(lsarp.mbx.lsarp.pca$sdev^2) * 100

# permanova
set.seed(25)
t1 = Sys.time()
lsarp.mbx.lsarp.pca.permanova = vegan::adonis2(dist(lsarp.mbx.lsarp.pca$x) ~ group*on.rs + ave.fiber,
                                               lsarp.mbx.lsarp.pca.df %>% mutate(on.rs = as.factor(ifelse(lsarp.on.rs=="on.rs", "on.rs", "off.rs"))),
                                               strata = lsarp.mbx.lsarp.pca.df$HM,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.mbx.lsarp.pca.permanova
# not sig

lsarp.mbx.lsarp.pca.plot <- ggplot(
  data=lsarp.mbx.lsarp.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2) %>%
    mutate(Group = group), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=ifelse(lsarp.on.rs=="on.rs", RS_Name, NA), shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value="grey")+
  facet_wrap(~Group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round((lsarp.mbx.lsarp.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round((lsarp.mbx.lsarp.pca.var)[2], digits=2), "%", sep=""),
       title=paste(paste("Metabolite\nTreatment R²: ", round(data.frame(lsarp.mbx.lsarp.pca.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.mbx.lsarp.pca.permanova)[2,5], 3), sep="")), sep="")
lsarp.mbx.lsarp.pca.plot


# :: MBX Interaction -----------------------------------------------

stools.not.flare = subset(metadata.lsarp.stool, 
                          flare != "flare" & !HM %in% excluded.by.dave)$standard.name %>% unique()
# add pseudocount, already filtered
lsarp.lmer.mbx.data = lsarp.delta.omic.prepare(lsarp.mbx.raw.mat.filt.1[rownames(lsarp.mbx.raw.mat.filt.1) %in% stools.not.flare,],
                                               filter = F,
                                               normalize=F)

lsarp.lmer.mbx.interactions = lsarp.delta.omic.lmer(lsarp.lmer.mbx.data,
                                                         split=F)

lsarp.lmer.mbx.interactions.plot = ggplot(lsarp.lmer.mbx.interactions%>%
                                                 mutate(phase = ifelse(phase == "treatment", "Treatment", "Washout")),
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub(" \\|.*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub(" \\|.*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Group Interaction (Metabolite)",
       y="FDR")+
  facet_wrap(~phase)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mbx.interactions.plot
# using annotated (closer to significance than unannotated)

subset(lsarp.lmer.mbx.interactions, padj < 0.20 & phase == "treatment")
subset(lsarp.lmer.mbx.interactions, padj < 0.20 & phase == "washout")

# :: Plots ----------------------------------------------------------------


lsarp.lmer.asv.interactions.plot+
  lsarp.lmer.mgx.interactions.plot+
  lsarp.lmer.mpx.kegg.interactions.plot+
  lsarp.lmer.mpx.cog.interactions.plot+
  lsarp.lmer.mpx.cazy.interactions.plot+
  lsarp.lmer.mbx.interactions.plot



# :: Group-Level Heatmaps -------------------------------------------------------------

# Take statistically significant features and plot Log2FC(?) heatmaps
# For Treatment and Washout separately
# Although, Washout is severely imbalanced and underpowered

# Log2FC from baseline

lsarp.lmer.omics.lfc.treatment = lsarp.delta(time.phase = "treatment",
                                             type = "log2fc")

# make matrix
lsarp.lmer.omics.lfc.treatment.mat = reshape2::acast(lsarp.lmer.omics.lfc.treatment,
                                                     standard.name ~ feature.name, value.var="lfc")
# make annotation maps
lsarp.lmer.omics.lfc.treatment.sample.map = lsarp.lmer.omics.lfc.treatment[,c("standard.name", "HM", "RS_Name", "Group", "month")] %>% distinct()
rownames(lsarp.lmer.omics.lfc.treatment.sample.map) = lsarp.lmer.omics.lfc.treatment.sample.map$standard.name
lsarp.lmer.omics.lfc.treatment.sample.map$standard.name = NULL
colnames(lsarp.lmer.omics.lfc.treatment.sample.map) = c("HM", "RS_Name", "Group", "Months on Product")

lsarp.lmer.omics.lfc.treatment.feature.map = rbind(subset(lsarp.lmer.asv.interactions, padj < 0.20 & phase == "treatment") %>% mutate(data.type = "ASV"),
                                                   subset(lsarp.lmer.mgx.interactions, padj < 0.20& phase == "treatment")%>% mutate(data.type = "Species"),
                                                   subset(lsarp.lmer.mpx.kegg.interactions %>% dplyr::select(-clean.name)%>% mutate(data.type = "Pathway"), padj < 0.20& phase == "treatment"),
                                                   subset(lsarp.lmer.mpx.cog.interactions %>% dplyr::select(-clean.name)%>% mutate(data.type = "COG"), padj < 0.20& phase == "treatment"),
                                                   subset(lsarp.lmer.mpx.cazy.interactions, padj < 0.20 & phase == "treatment") %>% mutate(data.type = "CAZy"),
                                                   subset(lsarp.lmer.mbx.interactions, padj < 0.20 & phase == "treatment") %>% mutate(data.type = "Metabolite"))
rownames(lsarp.lmer.omics.lfc.treatment.feature.map) = lsarp.lmer.omics.lfc.treatment.feature.map$taxa
lsarp.lmer.omics.lfc.treatment.feature.map = lsarp.lmer.omics.lfc.treatment.feature.map %>% dplyr::select(estimate, data.type)
colnames(lsarp.lmer.omics.lfc.treatment.feature.map) = c("Interaction", "Data type")

# make heatmap

# note: remove correlated features
#lsarp.lmer.omics.lfc.treatment.mat = lsarp.lmer.omics.lfc.treatment.mat %>% as.data.frame() %>% dplyr::select(-c(`Vibrio cholerae infection`, `Influenza A`, `Pathogenic Escherichia coli infection`, `Shigellosis`))

# scale by sample to find patient-level trends (ignore trends across features)
lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap = t(((lsarp.lmer.omics.lfc.treatment.mat )))
colMeans(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap)
# approx 0, means features were scaled per sample

pheatmap::pheatmap(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = lsarp.lmer.omics.lfc.treatment.sample.map %>% dplyr::select(-HM),
                   annotation_colors = list(Group = c(`RS` = gg_color_hue(2)[1],
                                                            `Placebo` = gg_color_hue(2)[2]),
                                            `Interaction` =  colorRampPalette(c("blue","white", "red"))(100),
                                            `Data type` = omics.colors,
                                            `RS_Name` = rs.colors),
                   annotation_row = lsarp.lmer.omics.lfc.treatment.feature.map,
                   breaks=c(seq(min(na.omit(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap))/100, max(na.omit(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap)), length.out=floor(100/2))))
# correlation clustering
# scaled Log2FC

# >>> 3. CLINICAL RESPONDERS --------------------------------------------------

lsarp.patient.flare = metadata.lsarp.stool.asv[,c("HM", "flare.call")] %>% distinct() %>%
  mutate(flare = ifelse(flare.call %in% c("No Flare") | flare.call == "ifx.between.visit2.3", 
                        "Remit", "Relapse"))
# code "no flare by rescope INCLUDING rescope" as "Remit"
# code "flare prior to rescope, or at rescope" as "Relapse"
# these have been double checked (based on records available to me) ad nauseum

# subset to treatment phase (plus baseline)
lsarp.metadata.responders = merge(metadata.lsarp.stool, lsarp.patient.flare[,c("HM", "flare")]%>%distinct(), by="HM") %>% subset(phase %in% c("treatment"))
lsarp.metadata.responders.asv = merge(metadata.lsarp.stool.asv, lsarp.patient.flare[,c("HM", "flare")]%>%distinct(), by="HM") %>% subset(phase %in% c("treatment"))
# fix names
colnames(lsarp.metadata.responders)[colnames(lsarp.metadata.responders)=="flare.x"] = "flare" # relapsed before endoscope
colnames(lsarp.metadata.responders)[colnames(lsarp.metadata.responders)=="flare.y"] = "flare.group" # group
colnames(lsarp.metadata.responders.asv)[colnames(lsarp.metadata.responders.asv)=="flare.x"] = "flare"
colnames(lsarp.metadata.responders.asv)[colnames(lsarp.metadata.responders.asv)=="flare.y"] = "flare.group"

lsarp.metadata.responders[,c("group", "flare.group", "HM")] %>% distinct() %>%
  dplyr::select(flare.group, group) %>% table()

saveRDS(lsarp.metadata.responders, "./2025_09_10_lsarp_cd_response_list.Rds")


# note: HM0865 and HM0938 have A in TI
# need precise definition of "flare"

# lsarp.metadata.responders$flare = ifelse(lsarp.metadata.responders$HM %in% c("HM0865", "HM0938"), "Relapse", lsarp.metadata.responders$flare)
# lsarp.metadata.responders.asv$flare = ifelse(lsarp.metadata.responders.asv$HM %in% c("HM0865", "HM0938"), "Relapse", lsarp.metadata.responders.asv$flare)
# Note: if you do this, the Butyrogen p value becomes > 0.3

# :: 6M Scope Outcome -----------------------------------------------------

# remove excluded patients
redcap.lsarp.wpcdai.filtered = redcap.lsarp.wpcdai %>%
  subset(!HM %in% excluded.by.dave)

redcap.lsarp.wpcdai.filtered = redcap.lsarp.wpcdai.filtered[,c("group","flare.call", "HM")] %>%
  as.data.frame() %>%
  mutate(condition = ifelse(flare.call %in% c("No Flare", "ifx.between.visit2.3"), "Remission", "Relapse"))

# tally up patients who did not flare by scope (6M)
redcap.lsarp.wpcdai.6m = redcap.lsarp.wpcdai.filtered %>% dplyr::select(group, condition) %>% table() %>% as.data.frame() %>%
  group_by(condition, group) %>% 
  mutate(total = sum(Freq)) %>% dplyr::select(condition, group, total) %>% distinct()

redcap.lsarp.wpcdai.6m.plot = ggplot(redcap.lsarp.wpcdai.6m %>% mutate(Group = ifelse(group == "Plac", "Placebo", "RS")),
       aes(x=Group, y=condition))+
  geom_tile(fill="white")+
  geom_tile(aes(fill=condition, alpha=total), color="white")+
  geom_text(aes(label=total), size=7)+
  #scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.major = element_blank())+
  labs(x="", y="")
redcap.lsarp.wpcdai.6m.plot

# add this to meta, too
lsarp.metadata.responders = merge(lsarp.metadata.responders,
                                  redcap.lsarp.wpcdai.filtered[,c("HM", "condition")],
                                  by="HM")
# flare.group = "did patient make it to scope without clinical flare"
# condition = "did patient make it past scope without flare" (e.g. no flare at )

table(paste(lsarp.metadata.responders$flare.group, lsarp.metadata.responders$condition))
# they're the same!


# >>> 4. RESPONSE ---------------------------------------------------------

# numbers
subset(lsarp.metadata.responders, group %in% c("RS")) %>%
  subset(flare != "flare"& !HM %in% excluded.by.dave) %>%
  dplyr::select(HM, flare.group) %>% table() # 1 to 6
subset(lsarp.metadata.responders, group %in% c("RS")) %>%
  subset(flare != "flare"& !HM %in% excluded.by.dave) %>%
  dplyr::select(HM, flare.group) %>% distinct() %>%
  dplyr::select(flare.group) %>% table() #  5 vs 9

# :: Fecal Calprotectin ----------------------------------------------------------

# does fecal calprotectin increase or decrease over RS treatment in Responders vs Non-Responders
# stats
lsarp.resp.stats.fcal.treatment = lmerTest::lmer(scale(log10(fcal)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                            subset(lsarp.metadata.responders, group %in% c("RS")) %>%
                                              subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.fcal.placebo = lmerTest::lmer(scale(log10(fcal)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                 subset(lsarp.metadata.responders, group %in% c("Placebo")) %>%
                                                 subset(flare != "flare" & !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.fcal.resp.plot = ggplot(lsarp.metadata.responders  %>%
                                         subset(!is.na(fcal))%>%
                                         subset(flare != "flare"& !HM %in% excluded.by.dave),
                                  aes(x=lsarp.days, 
                                      y=fcal))+
  scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders, 
  #                                         phase=="treatment")$lsarp.days),
  #         ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  geom_hline(yintercept=250, color="red")+
  labs(x="Days since starting Product", y="Fecal Calprotectin (μg/g)",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.fcal.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.fcal.treatment[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.fcal.resp.plot


# :: Stool Water ----------------------------------------------------------

# Note: change to "moisture" because dry weight excludes more than water (i.e. other volatiles)
# does stool water increase or decrease over RS treatment in Responders vs Non-Responders

# stats
lsarp.resp.stats.water.treatment = lmerTest::lmer(scale(stool_water_perc) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                 subset(lsarp.metadata.responders, group %in% c("RS")) %>%
                                                   subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.water.placebo = lmerTest::lmer(scale(stool_water_perc) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                               subset(lsarp.metadata.responders, group %in% c("Placebo")) %>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.water.resp.plot = ggplot(lsarp.metadata.responders %>%
                                         subset(!is.na(stool_water_perc))%>%
                                     subset(flare != "flare" & !HM %in% excluded.by.dave),
                                   aes(x=lsarp.days, 
                                       y=stool_water_perc*100))+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Stool Moisture (%)",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.water.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.water.treatment[5,5], 3))))+
   facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.water.resp.plot


# :: Microbial Load -------------------------------------------------------

# stats
lsarp.resp.stats.load.treatment = lmerTest::lmer(scale(load.asv) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                 subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                   subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.load.placebo = lmerTest::lmer(scale(load.asv) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                               subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                 subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.load.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                         subset(!is.na(load.asv))%>%
                                         subset(flare!="flare"& !HM %in% excluded.by.dave),
                                  aes(x=lsarp.days,
                                      y=load.asv))+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Microbial Load",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.load.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.load.treatment[5,5], 3))))+
    facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.load.resp.plot


# :: ASV Richness ------------------------------------------------------

# stats
lsarp.resp.stats.richness.resp.treatment = lmerTest::lmer(scale(log10(richness)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                 subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                   subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.richness.resp.placebo = lmerTest::lmer(scale(log10(richness)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                               subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                 subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.stool.asv.richness.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                                       subset(!is.na(richness))%>%
                                                       subset(flare!="flare"& !HM %in% excluded.by.dave),
                                                aes(x=lsarp.days, 
                                                    y=richness))+
  # scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="ASV Richness",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.richness.resp.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.richness.resp.treatment[5,5], 3))))+
    facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.richness.resp.plot


# :: ASV Shannon ------------------------------------------------------


# stats
lsarp.resp.stats.shannon.treatment = lmerTest::lmer(scale(shannon) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                     subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                      subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.shannon.placebo = lmerTest::lmer(scale(shannon) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                   subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                    subset(flare!="flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

# plot
metadata.lsarp.stool.asv.shannon.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                                      subset(!is.na(shannon))%>%
                                                      subset(flare!="flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=shannon))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Shannon Diversity",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.shannon.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.shannon.treatment[5,5], 3))))+
    facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.shannon.resp.plot


# :: ASV Functional Redundancy ------------------------------------------------------

# stats
lsarp.resp.stats.fd.treatment = lmerTest::lmer(scale(fd) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                    subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.fd.placebo = lmerTest::lmer(scale(fd) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                  subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                               subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.stool.asv.fd.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                                 subset(!is.na(fd))%>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                          aes(x=lsarp.days, 
                                              y=fd))+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Functional Redundancy",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.fd.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.fd.treatment[5,5], 3))))+
    facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.fd.resp.plot


# :: ASV Functional Richness ------------------------------------------------------

# stats
lsarp.resp.stats.fr.treatment = lmerTest::lmer(scale(fr) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                               subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.fr.placebo = lmerTest::lmer(scale(fr) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                             subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                               subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.stool.asv.fr.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                                 subset(!is.na(fr))%>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=fr))+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Functional Richness",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.fr.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.fr.treatment[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.stool.asv.fr.resp.plot

# :: ASV Butyrogens I -------------------------------------------------------

# stats
lsarp.resp.stats.but.i.treatment = lmerTest::lmer(scale(log10(but.i)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                               subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.but.i.placebo = lmerTest::lmer(scale(log10(but.i)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                             subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                               subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.but.i.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                          subset(!is.na(but.i))%>%
                                          subset(flare !="flare"& !HM %in% excluded.by.dave),
                                   aes(x=lsarp.days, 
                                       y=but.i*100))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Butyrogens %",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.but.i.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.but.i.treatment[5,5], 3))))+
    facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.but.i.resp.plot


# :: ASV Butyrogens II (Vital) -------------------------------------------------------

# stats
lsarp.resp.stats.but.ii.treatment = lmerTest::lmer(scale(log10(but.ii)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                  subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                    subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.but.ii.placebo = lmerTest::lmer(scale(log10(but.ii)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                  subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.but.ii.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                           subset(!is.na(but.ii))%>%
                                           subset(flare != "flare"& !HM %in% excluded.by.dave),
                                        aes(x=lsarp.days, 
                                            y=but.ii*100))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Butyrogens %",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.but.ii.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.but.ii.treatment[5,5], 3))))+
    facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.but.ii.resp.plot


# :: ASV Bray-Curtis ------------------------------------------------------

# perform PCoA

# to make results comparable to group-level, must glom and prevalence filter

# apply 20% prevalence filter (first time doing this in this analysis)
# keep important samples (not flare and not excluded)
lsarp.asv.data.resp.median = lsarp.asv.data.glom[rownames(lsarp.asv.data.glom) %in% 
                                                     subset(lsarp.metadata.responders.asv, flare !="flare"& !HM %in% excluded.by.dave)$standard.name,]

dim(lsarp.asv.data.resp.median) # 1554 glommed taxa
lsarp.asv.data.resp.median.pa = lsarp.asv.data.resp.median
lsarp.asv.data.resp.median.pa[lsarp.asv.data.resp.median.pa > 0] = 1
# 10%
lsarp.asv.data.resp.median.filt = lsarp.asv.data.resp.median[,colSums(lsarp.asv.data.resp.median.pa) >= nrow(lsarp.asv.data.resp.median.pa)*0.1]
dim(lsarp.asv.data.resp.median.filt) # 426 taxa
# 20%
lsarp.asv.data.resp.median.filt = lsarp.asv.data.resp.median[,colSums(lsarp.asv.data.resp.median.pa) >= nrow(lsarp.asv.data.resp.median.pa)*0.2]
ncol(lsarp.asv.data.resp.median.filt) # 304 taxa


# calculate Bray-Curtis dissimilarities [ASVs are not prevalence filtered!]
lsarp.asv.resp.bray = vegan::vegdist(lsarp.asv.data.resp.median.filt, method="bray") 
# perform PCoA
lsarp.asv.resp.pcoa = ape::pcoa(lsarp.asv.resp.bray)
# extract data from pcoa
lsarp.asv.resp.pcoa.df = data.frame(lsarp.asv.resp.pcoa$vectors[,c(1:2)])
lsarp.asv.resp.pcoa.df$standard.name = rownames(lsarp.asv.resp.pcoa.df)
# add metadata
lsarp.asv.resp.pcoa.df = merge(lsarp.asv.resp.pcoa.df,
                               lsarp.metadata.responders.asv, by="standard.name")
# extract variance explained
lsarp.asv.resp.pcoa.var_exp = lsarp.asv.resp.pcoa$values[c(1:2),2]
lsarp.asv.resp.pcoa.df$var1 = round(lsarp.asv.resp.pcoa.var_exp[1]*100, digits=2)
lsarp.asv.resp.pcoa.df$var2 = round(lsarp.asv.resp.pcoa.var_exp[2]*100, digits=2)

# clean up "lsarp.on.rs" variable

set.seed(25)
t1 = Sys.time()
lsarp.asv.resp.permanova = vegan::adonis2(lsarp.asv.resp.bray ~ flare.group*group*lsarp.days + ave.fiber,
                                     lsarp.asv.resp.pcoa.df %>% mutate(on.rs = ifelse(lsarp.on.rs == "on.rs", "on.rs", "off.rs")),
                                     strata = lsarp.asv.resp.pcoa.df$HM,
                                     by="margin")
t2 = Sys.time()
t2 - t1
lsarp.asv.resp.permanova # not-sig

lsarp.asv.resp.pcoa.plot <- ggplot(
  data=lsarp.asv.resp.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM, color=flare.group), alpha=0.5, linetype=2, linewidth=0.3) + 
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2.5)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  #geom_text(aes(label=HM))+
  labs(x=paste("Axis 1: ", round(unique(lsarp.asv.resp.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.asv.resp.pcoa.df$var2), digits=2), "%", sep=""),
       title=paste(paste("Response:Group:Time R²: ", round(data.frame(lsarp.asv.resp.permanova)[2,3], 3)*100, "%",
                         "  p: ", round(data.frame(lsarp.asv.resp.permanova)[2,5], 3), sep="")), sep="") + 
  facet_wrap(~group,scales="free")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.asv.resp.pcoa.plot


# :: ASV Beta-Diversity ------------------------------------------------------

beta.trajectory.resp.data = beta.trajectory(lsarp.asv.resp.bray)

# Between (using new Bray-Curtis)
beta.trajectory.resp.data = merge(beta.trajectory.resp.data,
                             lsarp.metadata.responders[,c("standard.name","group","ave.fiber", "HM","lsarp.days","flare","baseline", "flare.group")],
                             by="standard.name")

# stats
lsarp.resp.stats.between.beta.treatment = lmerTest::lmer(scale(log10(between.beta)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                  subset(beta.trajectory.resp.data, group %in% c("RS")) %>%
                                                    subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.between.beta.placebo = lmerTest::lmer(scale(log10(between.beta)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                subset(beta.trajectory.resp.data, group %in% c("Placebo")) %>%
                                                  subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.between.beta.resp.plot = ggplot(beta.trajectory.resp.data %>%
                                                 subset(!is.na(between.beta))%>%
                                                 subset(flare !="flare"& !HM %in% excluded.by.dave),
                                        aes(x=lsarp.days, 
                                            y=between.beta*100))+
  scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Bray-Curtis Dissimilarity",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.between.beta.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.between.beta.treatment[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.between.beta.resp.plot

# Within

# stats
lsarp.resp.stats.within.beta.treatment = lmerTest::lmer(scale(log10(within.beta)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                         subset(beta.trajectory.resp.data, group %in% c("RS")) %>%
                                                           subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.within.beta.placebo = lmerTest::lmer(scale(log10(within.beta)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                       subset(beta.trajectory.resp.data, group %in% c("Placebo")) %>%
                                                         subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.within.beta.resp.plot = ggplot(beta.trajectory.resp.data %>%
                                                subset(!is.na(within.beta))%>%
                                                subset(flare !="flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=within.beta*100))+
  scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=3)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Bray-Curtis Dissimilarity",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.within.beta.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.within.beta.treatment[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.within.beta.resp.plot


# :: ASV Interaction ------------------------------------------------------


lsarp.lmer.asv.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.lmer.asv.data[rownames(lsarp.lmer.asv.data) %in% stools.not.flare,],
                                                              split=F)

# fix names
lsarp.lmer.asv.resp.interactions$taxa =  ifelse(grepl("s__", (lsarp.lmer.asv.resp.interactions$taxa)), paste("(s) ", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.lmer.asv.resp.interactions$taxa)))), sep=""),
                                 paste("(", substr((lsarp.lmer.asv.resp.interactions$taxa), 1, 1), ") ", gsub("o__", "", gsub("c__", "", gsub("g__", "", gsub("f__", "", gsub("s__", "", (lsarp.lmer.asv.resp.interactions$taxa)))))), sep=""))


lsarp.lmer.asv.resp.interactions.plot = ggplot(lsarp.lmer.asv.resp.interactions,
                                          aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (ASV)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.asv.resp.interactions.plot
lsarp.lmer.asv.resp.interactions
# several padj < 0.20, including RS responders enriched w F. prausnitzii

subset(lsarp.lmer.asv.resp.interactions, padj < 0.20 & group == "Placebo")
subset(lsarp.lmer.asv.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.asv.resp.interactions, padj < 0.20)[,c("taxa")] %>% table()
# Lachnospiraceae was impacted in both RS and placebo groups, so exclude from interpretation

# :: MGX Interaction ------------------------------------------------------

lsarp.lmer.mgx.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.lmer.mgx.data[rownames(lsarp.lmer.mgx.data) %in% stools.not.flare,])

lsarp.lmer.mgx.resp.interactions.plot = ggplot(lsarp.lmer.mgx.resp.interactions,
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj > 0.2, gsub("\\..*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (Species)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mgx.resp.interactions.plot
lsarp.lmer.mgx.resp.interactions
# few padj < 0.20

subset(lsarp.lmer.mgx.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.mgx.resp.interactions, padj < 0.20)[,c("taxa")] %>% table()


# :: Pathway Interaction ------------------------------------------------------


lsarp.lmer.mpx.kegg.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.lmer.mpx.kegg.data[rownames(lsarp.lmer.mpx.kegg.data) %in% stools.not.flare,],
                                                                   split=F)

lsarp.lmer.mpx.kegg.resp.interactions.plot = ggplot(lsarp.lmer.mpx.kegg.resp.interactions,
                                                    aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\.", " ", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\.", " ", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (Pathway)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.kegg.resp.interactions.plot
lsarp.lmer.mpx.kegg.resp.interactions
# no padj < 0.20 in RS
# but Butanoate metabolism increased in Placebo responders

subset(lsarp.lmer.mpx.kegg.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.mpx.kegg.resp.interactions, padj < 0.20)[,c("taxa")] %>% table()


# :: COG Interaction ------------------------------------------------------


lsarp.lmer.mpx.cog.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.lmer.mpx.cog.data[rownames(lsarp.lmer.mpx.cog.data) %in% stools.not.flare,],
                                                                   split=F)

lsarp.lmer.mpx.cog.resp.interactions.plot = ggplot(lsarp.lmer.mpx.cog.resp.interactions,
                                                    aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, taxa, NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (COG)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.cog.resp.interactions.plot
lsarp.lmer.mpx.cog.resp.interactions

subset(lsarp.lmer.mpx.cog.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.mpx.cog.resp.interactions, padj < 0.20)[,c("taxa")] %>% table() %>% range()


# :: CAZy Interaction ------------------------------------------------------


lsarp.lmer.mpx.cazy.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.lmer.mpx.cazy.data[rownames(lsarp.lmer.mpx.cazy.data) %in% stools.not.flare,],
                                                                   split=F)

lsarp.lmer.mpx.cazy.resp.interactions.plot = ggplot(lsarp.lmer.mpx.cazy.resp.interactions,
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub("\\..*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (CAZy)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mpx.cazy.resp.interactions.plot
lsarp.lmer.mpx.cazy.resp.interactions
# enriched GT84 in RS responders

subset(lsarp.lmer.mpx.cazy.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.mpx.cazy.resp.interactions, padj < 0.20)[,c("taxa")] %>% table()


# :: CAZy Starch Interaction -------------------------------------------------------

# stats
lsarp.resp.stats.starch.treatment = lmerTest::lmer(scale(log10(starch)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                         subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                           subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.starch.placebo = lmerTest::lmer(scale(log10(starch)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                       subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                         subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.starch.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                           subset(!is.na(starch))%>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=starch))+
  scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Starch CAZy Intensity",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.starch.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.starch.treatment[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.starch.resp.plot


# :: CAZy Mucin Interaction -------------------------------------------------------

# stats
lsarp.resp.stats.mucin.treatment = lmerTest::lmer(scale(log10(mucin)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                         subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                           subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.mucin.placebo = lmerTest::lmer(scale(log10(mucin)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                       subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                         subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.mucin.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                          subset(!is.na(mucin))%>%
                                          subset(flare != "flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=mucin))+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Mucin CAZy Intensity",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.mucin.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.mucin.treatment[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.mucin.resp.plot

# :: CAZy Starch:Mucin Interaction -------------------------------------------------------

# stats
lsarp.resp.stats.starch.mucin.treatment = lmerTest::lmer(scale((starch.mucin)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                   subset(lsarp.metadata.responders.asv, group %in% c("RS")) %>%
                                                     subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.starch.mucin.placebo = lmerTest::lmer(scale((starch.mucin)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                 subset(lsarp.metadata.responders.asv, group %in% c("Placebo")) %>%
                                                   subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()


# plot
metadata.lsarp.starch.mucin.resp.plot = ggplot(lsarp.metadata.responders.asv %>%
                                                 subset(!is.na(starch.mucin))%>%
                                                 subset(flare != "flare"& !HM %in% excluded.by.dave),
                                         aes(x=lsarp.days, 
                                             y=starch.mucin))+
  #annotate("rect", xmin=0, xmax=max(subset(lsarp.metadata.responders.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=2)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Starch:Mucin Ratio",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.starch.mucin.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.starch.mucin.treatment[5,5], 3))))+
  facet_wrap(~group,scales="free_x")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.starch.mucin.resp.plot


# :: MBX Filter ------------------------------------------------------------------

# filter again based on responder vs non-responder

# 80% prevalence; in at least one TREATMENT GROUP * RESPONSE GROUP (4 possible groups)
dim(lsarp.mbx.raw.mat)
lsarp.mbx.raw.mat.filt.2 = lsarp.mbx.annotated.mat %>% as.matrix() %>%
  reshape2::melt()
colnames(lsarp.mbx.raw.mat.filt.2)[1] = "standard.name"
lsarp.mbx.raw.mat.filt.2 = merge(lsarp.mbx.raw.mat.filt.2,
                                 metadata.lsarp.stool[,c("standard.name", "Group", "flare")], by="standard.name")
lsarp.mbx.raw.mat.filt.2$value = ifelse(lsarp.mbx.raw.mat.filt.2$value >0, 1, 0)
# calc prevalence per group
lsarp.mbx.raw.mat.filt.2 = lsarp.mbx.raw.mat.filt.2 %>%
  group_by(Var2, flare, Group) %>%
  summarise(prevalence = mean(value, na.rm = TRUE)) %>% # cool, use mean of presence to calculate prevalence
  ungroup() %>%
  subset(prevalence >= .80)
# keep features with 80% prevalence in at least one of the groups
length(unique(lsarp.mbx.raw.mat.filt.2$Var2))
# n = 187 features (annotated)
# apply filter
lsarp.mbx.raw.mat.filt.2 = lsarp.mbx.annotated.mat[,colnames(lsarp.mbx.annotated.mat) %in% unique(lsarp.mbx.raw.mat.filt.2$Var2)]
dim(lsarp.mbx.raw.mat.filt.2)


# :: MBX Interaction -----------------------------------------------

# prep
lsarp.mbx.raw.mat.filt.2

lsarp.mbx.resp.prep = lsarp.delta.omic.prepare(lsarp.mbx.raw.mat.filt.2,
                                               filter = F)

lsarp.lmer.mbx.resp.interactions = lsarp.delta.omic.resp.lmer(lsarp.mbx.resp.prep[rownames(lsarp.mbx.resp.prep) %in% stools.not.flare,])

lsarp.lmer.mbx.resp.interactions.plot = ggplot(lsarp.lmer.mbx.resp.interactions,
                                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub(" \\|.*", "", taxa), NA)),
                           size=2.5)+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05 & padj >0.2, gsub(" \\|.*", "", taxa), NA)),
                           size=2.5, color="grey")+
  labs(x="Interaction Coefficient",
       title="Response Interaction (Metabolite)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
lsarp.lmer.mbx.resp.interactions.plot

# tetraethylene glycol enriched in RS responders

subset(lsarp.lmer.mbx.resp.interactions, padj < 0.20 & group == "RS")
subset(lsarp.lmer.mbx.resp.interactions, padj < 0.20)[,c("taxa")] %>% table()


# :: Responders Heatmaps -------------------------------------------------------------


lsarp.lmer.omics.lfc.resp.treatment = lsarp.delta(time.phase = "treatment",
                                                  responders=T,
                                             type = "log2fc")

# make matrix
lsarp.lmer.omics.lfc.resp.treatment.mat = reshape2::acast(lsarp.lmer.omics.lfc.resp.treatment,
                                                     standard.name ~ feature.name, value.var="lfc")
# make annotation maps
lsarp.lmer.omics.lfc.resp.treatment.sample.map = lsarp.lmer.omics.lfc.resp.treatment[,c("standard.name", "HM", "RS_Name", "Group", "month", "flare.group")] %>% distinct()
rownames(lsarp.lmer.omics.lfc.resp.treatment.sample.map) = lsarp.lmer.omics.lfc.resp.treatment.sample.map$standard.name
lsarp.lmer.omics.lfc.resp.treatment.sample.map$standard.name = NULL
colnames(lsarp.lmer.omics.lfc.resp.treatment.sample.map) = c("HM", "RS_Name", "Group", "Months on Product", "Response")

lsarp.lmer.omics.lfc.resp.treatment.feature.map = rbind(subset(lsarp.lmer.asv.resp.interactions, padj < 0.20 & group == "RS") %>% mutate(data.type = "ASV"),
                                                   subset(lsarp.lmer.mgx.resp.interactions, padj < 0.20 & group == "RS")%>% mutate(data.type = "Species"),
                                                   subset(lsarp.lmer.mpx.kegg.resp.interactions %>% mutate(data.type = "Pathway"), padj < 0.20 & group == "RS"),
                                                   subset(lsarp.lmer.mpx.cog.resp.interactions %>% mutate(data.type = "COG"), padj < 0.20 & group == "RS"),
                                                   subset(lsarp.lmer.mpx.cazy.resp.interactions, padj < 0.20 & group == "RS") %>% mutate(data.type = "CAZy"),
                                                   subset(lsarp.lmer.mbx.resp.interactions, padj < 0.20 & group == "RS") %>% mutate(data.type = "Metabolite"))
rownames(lsarp.lmer.omics.lfc.resp.treatment.feature.map) = lsarp.lmer.omics.lfc.resp.treatment.feature.map$taxa
lsarp.lmer.omics.lfc.resp.treatment.feature.map = lsarp.lmer.omics.lfc.resp.treatment.feature.map %>% dplyr::select(estimate, data.type)
colnames(lsarp.lmer.omics.lfc.resp.treatment.feature.map) = c("Interaction", "Data type")

# make heatmap

# note: remove correlated features
# lsarp.lmer.omics.lfc.resp.treatment.mat = lsarp.lmer.omics.lfc.resp.treatment.mat %>% as.data.frame() %>% dplyr::select(-c(`Vibrio cholerae infection`, `Influenza A`, `Pathogenic Escherichia coli infection`, `Shigellosis`))

# scale by sample to find patient-level trends (ignore trends across features)
lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap = ((t(lsarp.lmer.omics.lfc.resp.treatment.mat )))
colMeans(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap)
# approx 0, means features were scaled per sample

pheatmap::pheatmap(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap,
                   #scale = "row",
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = lsarp.lmer.omics.lfc.resp.treatment.sample.map %>% dplyr::select(-HM),
                   annotation_colors = list(Response = c(`Relapse` = gg_color_hue(2)[1],
                                                      `Remit` = gg_color_hue(2)[2]),
                                            Group = c(`RS` = gg_color_hue(2)[1],
                                                         `Placebo` = gg_color_hue(2)[2]),
                                            `Interaction` =  colorRampPalette(c("blue","white", "red"))(100),
                                            `Data type` = omics.colors,
                                            `RS_Name` = rs.colors),
                   annotation_row = lsarp.lmer.omics.lfc.resp.treatment.feature.map,
                   breaks=c(seq(min(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap))/100, max(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap)), length.out=floor(100/2))))
# correlation clustering



# >>> 5. FCAL-OMICS ---------------------------------------------------

# merge all multi-omic data
lsarp.fecalcal.omics.data = merge(lsarp.asv.data.glom/50000,
                                  lsarp.mgx.taxa,
                                  by="row.names", all=T)
lsarp.fecalcal.omics.data = merge(lsarp.fecalcal.omics.data,
                                  lsarp.cd.mpx.cog.mat %>%as.data.frame()%>%
                                    mutate(`Row.names` = rownames(.)),
                                  by="Row.names", all=T)
lsarp.fecalcal.omics.data = merge(lsarp.fecalcal.omics.data,
                                  lsarp.cd.mpx.kegg.mat%>%as.data.frame()%>%
                                    mutate(`Row.names` = rownames(.)),
                                  by="Row.names", all=T)
lsarp.fecalcal.omics.data = merge(lsarp.fecalcal.omics.data,
                                  lsarp.cd.mpx.cazy.mat%>%as.data.frame()%>%
                                    mutate(`Row.names` = rownames(.)),
                                  by="Row.names", all=T)
lsarp.fecalcal.omics.data = merge(lsarp.fecalcal.omics.data,
                                  lsarp.mbx.raw.mat.filt.1%>%as.data.frame()%>%
                                    mutate(`Row.names` = rownames(.)),
                                  by="Row.names", all=T)


lsarp.fecalcal.omics.data$standard.name = lsarp.fecalcal.omics.data$`Row.names`
# add fecal cal
lsarp.fecalcal.omics.data = merge(lsarp.fecalcal.omics.data,
                                  metadata.lsarp.stool[,c("fcal", "standard.name")],
                                  by="standard.name")
rownames(lsarp.fecalcal.omics.data) = lsarp.fecalcal.omics.data$Row.names
lsarp.fecalcal.omics.data$Row.names=NULL
# ready

dim(lsarp.fecalcal.omics.data) # 5872 features x 192 samples
rownames(lsarp.fecalcal.omics.data) = lsarp.fecalcal.omics.data$standard.name
lsarp.fecalcal.omics.data$standard.name = NULL
# replace NA with 0
lsarp.fecalcal.omics.data[is.na(lsarp.fecalcal.omics.data)] = 0

# remove low prevalent (< 50%)
lsarp.fecalcal.omics.data.pa = lsarp.fecalcal.omics.data
lsarp.fecalcal.omics.data.pa[lsarp.fecalcal.omics.data.pa!=0] = 1
cols.to.keep = data.frame(keep = colSums(lsarp.fecalcal.omics.data.pa)>
                            (nrow(lsarp.fecalcal.omics.data)*0.50))

lsarp.fecalcal.omics.data = lsarp.fecalcal.omics.data[,rownames(subset(cols.to.keep, keep == T))]
dim(lsarp.fecalcal.omics.data) # 2571 features (192 samples)
colnames(lsarp.fecalcal.omics.data) %>% tail()

# convert to lfc (one feature at a time)
lsarp.fecalcal.omics.data.lfc = do.call(cbind, lapply(1:ncol(lsarp.fecalcal.omics.data), function(col){
  print(col)
  data.subset = data.frame(sample = rownames(lsarp.fecalcal.omics.data),
                           feature = lsarp.fecalcal.omics.data[,col])
  data.subset = tidyr::separate(data.subset, col="sample", into=c("HM", "stl", "no"), sep="-", remove=F)
  data.subset = data.subset %>%
    group_by(HM) %>%
    mutate(feature = feature + (min(feature[feature!=0])/2)) %>%
    arrange(no) %>% # ensure stools are ordered by collection time
    mutate(lfc = log2(feature / lag(feature))) # lag takes the preceding value
  new.data = data.frame(lfc = data.subset$lfc)
  colnames(new.data)[1] = colnames(lsarp.fecalcal.omics.data)[col]
  rownames(new.data) = data.subset$sample
  return(new.data)
})) %>% as.data.frame()

# save
lsarp.fecalcal.omics.data.saved = lsarp.fecalcal.omics.data.lfc
lsarp.fecalcal.omics.data.lfc = lsarp.fecalcal.omics.data.saved

rownames(lsarp.fecalcal.omics.data.lfc) # = lsarp.fecalcal.omics.data$standard.name
lsarp.fecalcal.omics.data.lfc$standard.name = NULL


# :: Correlation ----------------------------------------------------------

# spearman correlation (takes a min)
metadata.lsarp.stool.omics.cor = Hmisc::rcorr(as.matrix(lsarp.fecalcal.omics.data.lfc[,colnames(lsarp.fecalcal.omics.data.lfc)!="standard.name"]+1), type="spearman")
dim(lsarp.fecalcal.omics.data.lfc) # 192 samples

metadata.lsarp.stool.omics.cor.df = 
  data.frame(reshape2::melt(metadata.lsarp.stool.omics.cor$r)) %>%
  mutate(pval = reshape2::melt(metadata.lsarp.stool.omics.cor$P)$value) %>%
  subset(!is.na(pval)) %>%
  arrange(pval)

metadata.lsarp.stool.omics.cor.bh = subset(metadata.lsarp.stool.omics.cor.df, Var1 == "fcal") %>%
  mutate(feature = as.character(Var2))%>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

lsarp.fecalcal.omics.cor.plot = ggplot(metadata.lsarp.stool.omics.cor.bh,
                                       aes(x=value, y=(padj)))+
  geom_point(shape=21, size=2.5, aes(fill = ifelse(padj < 0.20, "sig", "notsig")))+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)),
                           size=3)+  
  scale_fill_manual(values=c("black", "white"))+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  facet_wrap(~"Fecal Calprotectin Correlation")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="FDR")
lsarp.fecalcal.omics.cor.plot

lsarp.fecalcal.omics.cor.df = reshape2::melt(lsarp.fecalcal.omics.data.lfc[,colnames(lsarp.fecalcal.omics.data.lfc)%in%c("fcal", subset(metadata.lsarp.stool.omics.cor.bh, padj < 0.20)$feature)],
                                             id.vars="fcal")

lsarp.fecalcal.omics.cor.df$clean.feature = ifelse(nchar(as.character(lsarp.fecalcal.omics.cor.df$variable))>30, paste(substr(lsarp.fecalcal.omics.cor.df$variable, 1, 30), "...", sep=""), as.character(lsarp.fecalcal.omics.cor.df$variable))

# add omics shape

lsarp.fecalcal.omics.cor.df = lsarp.fecalcal.omics.cor.df %>%
  mutate(data.type = ifelse(variable %in% colnames(lsarp.mgx.taxa), "Species",
                            ifelse(variable %in% colnames(lsarp.cd.mpx.kegg.mat), "Pathway", 
                                   ifelse(variable %in% colnames(lsarp.cd.mpx.cog.mat), "COG", 
                                          ifelse(variable %in% colnames(lsarp.cd.mpx.cazy.mat), "CAZy", 
                                                 ifelse(variable %in% colnames(lsarp.asv.data.glom), "ASV",
                                                        ifelse(variable %in% colnames(lsarp.mbx.raw.mat.filt.1), "Metabolite", "Fecal calprotectin")))))))

lsarp.fecalcal.omics.cor.df$data.type = factor(lsarp.fecalcal.omics.cor.df$data.type, levels=c("Fecal calprotectin", "ASV", "Species", "COG", "Pathway", "CAZy", "Metabolite"))

lsarp.fecalcal.omics.cor.data = lsarp.fecalcal.omics.cor.df %>%
  subset(variable %in% slice_min(metadata.lsarp.stool.omics.cor.bh, order_by=abs(pval), n=8)$feature) %>%
  arrange(order_by=value)  %>%
  subset(!is.na(value)) %>%
  mutate(Var2 = variable)

lsarp.fecalcal.omics.cor.data.order = merge(lsarp.fecalcal.omics.cor.data[,c("Var2", "clean.feature")],
                                            metadata.lsarp.stool.omics.cor.bh[,c("Var2", "value")], by="Var2") %>%
  dplyr::select(clean.feature, value) %>% distinct() %>%
  arrange(value)
lsarp.fecalcal.omics.cor.data = lsarp.fecalcal.omics.cor.data %>%
  mutate(clean.feature = factor(clean.feature, levels=lsarp.fecalcal.omics.cor.data.order$clean.feature))

lsarp.fecalcal.omics.cor.plots = ggplot(lsarp.fecalcal.omics.cor.data,
                                        aes(x=value,
                                            y=fcal))+
  geom_point(size = 2, aes(shape=data.type, fill=data.type), color="black", alpha=0.9) +
  scale_shape_manual(values=omics.shapes)+
  guides(fill=FALSE, shape=FALSE, size = FALSE)+
  scale_fill_manual(values= omics.colors)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(aes(label = ..r.label..), method="spearman", size=3, 
                   label.y=Inf, vjust=1.5,
                   label.x=Inf,hjust=1.5)+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_classic()+theme(strip.text=element_text(size=8))+
  facet_wrap(~clean.feature, scales="free", nrow=2)+
  labs(y=expression(Calprotectin~Log[2]*FC),
       x=expression(Feature~Log[2]*FC))
lsarp.fecalcal.omics.cor.plots


lsarp.fecalcal.omics.cor.plot+
lsarp.fecalcal.omics.cor.plots
# save these

# >>> 6. PREDICTORS -------------------------------------------------------

# Note: ML may be inappropriate to use, given the small sample size (minority class: n = 5)

# Just do PC/oA and Wilcoxon


# :: Data -----------------------------------------------------------------

# omics data:
lsarp.asv.data.glom
lsarp.mgx.taxa
lsarp.cd.mpx.kegg.mat
lsarp.cd.mpx.cog.mat
lsarp.cd.mpx.cazy.mat

# baseline RS-group samples
lsarp.cd.baseline.samples = subset(metadata.lsarp.stool.asv, 
                                   # Subset to baseline
                                   baseline=="baseline" & 
                                     # AND RS group only
                                     Group == "RS")
# and remove excluded participants
lsarp.cd.baseline.samples = subset(lsarp.cd.baseline.samples, 
                                   !HM %in% excluded.by.dave)$standard.name

# process + 20% filter each dataset:

# ASV
lsarp.asv.baseline = lsarp.asv.data.glom[rownames(lsarp.asv.data.glom) %in% lsarp.cd.baseline.samples,] / 50000
lsarp.asv.baseline.pa = lsarp.asv.baseline
lsarp.asv.baseline.pa[lsarp.asv.baseline.pa>0] = 1
lsarp.asv.baseline = lsarp.asv.baseline[,colSums(lsarp.asv.baseline.pa) >= nrow(lsarp.asv.baseline.pa)*0.2]
lsarp.asv.baseline = lsarp.asv.baseline %>% as.data.frame()
# filtered
# lsarp.asv.baseline = log2(lsarp.asv.baseline+min(lsarp.asv.baseline[lsarp.asv.baseline!=0])/2) %>% as.data.frame()
lsarp.asv.baseline$standard.name = rownames(lsarp.asv.baseline)
dim(lsarp.asv.baseline) # 14 samples * 294 features

# MGX
lsarp.mgx.baseline = lsarp.mgx.taxa[rownames(lsarp.mgx.taxa) %in% lsarp.cd.baseline.samples,] / 50000
lsarp.mgx.baseline.pa = lsarp.mgx.baseline
lsarp.mgx.baseline.pa[lsarp.mgx.baseline.pa>0] = 1
lsarp.mgx.baseline = lsarp.mgx.baseline[,colSums(lsarp.mgx.baseline.pa) >= nrow(lsarp.mgx.baseline.pa)*0.2]
lsarp.mgx.baseline = lsarp.mgx.baseline %>% as.data.frame()
# filtered
# lsarp.mgx.baseline = log2(lsarp.mgx.baseline+min(lsarp.mgx.baseline[lsarp.mgx.baseline!=0])/2) %>% as.data.frame()
lsarp.mgx.baseline$standard.name = rownames(lsarp.mgx.baseline)
dim(lsarp.mgx.baseline) # 14 samples * 337 features

# Pathways
lsarp.kegg.baseline = lsarp.cd.mpx.kegg.mat[rownames(lsarp.cd.mpx.kegg.mat) %in% lsarp.cd.baseline.samples,]
lsarp.kegg.baseline[is.na(lsarp.kegg.baseline)] = 0
lsarp.kegg.baseline.pa = lsarp.kegg.baseline
lsarp.kegg.baseline.pa[lsarp.kegg.baseline.pa>0] = 1
lsarp.kegg.baseline = lsarp.kegg.baseline[,colSums(lsarp.kegg.baseline.pa) >= nrow(lsarp.kegg.baseline.pa)*0.2]
lsarp.kegg.baseline = lsarp.kegg.baseline %>% as.data.frame()
# filtered
# lsarp.kegg.baseline = log2(lsarp.kegg.baseline+min(lsarp.kegg.baseline[lsarp.kegg.baseline!=0])/2) %>% as.data.frame()
lsarp.kegg.baseline$standard.name = rownames(lsarp.kegg.baseline)
dim(lsarp.kegg.baseline) # 13 samples * 180 features

# COGs
lsarp.cog.baseline = lsarp.cd.mpx.cog.mat[rownames(lsarp.cd.mpx.cog.mat) %in% lsarp.cd.baseline.samples,]
lsarp.cog.baseline[is.na(lsarp.cog.baseline)] = 0
lsarp.cog.baseline.pa = lsarp.cog.baseline
lsarp.cog.baseline.pa[lsarp.cog.baseline.pa>0] = 1
lsarp.cog.baseline = lsarp.cog.baseline[,colSums(lsarp.cog.baseline.pa) >= nrow(lsarp.cog.baseline.pa)*0.2]
lsarp.cog.baseline = lsarp.cog.baseline %>% as.data.frame()
# filtered
# lsarp.cog.baseline = log2(lsarp.cog.baseline+min(lsarp.cog.baseline[lsarp.cog.baseline!=0])/2) %>% as.data.frame()
lsarp.cog.baseline$standard.name = rownames(lsarp.cog.baseline)
dim(lsarp.cog.baseline) # 13 samples * 2355 features

# CAZy
lsarp.cazy.baseline = lsarp.cd.mpx.cazy.mat[rownames(lsarp.cd.mpx.cazy.mat) %in% lsarp.cd.baseline.samples,]
lsarp.cazy.baseline[is.na(lsarp.cazy.baseline)] = 0
lsarp.cazy.baseline.pa = lsarp.cazy.baseline
lsarp.cazy.baseline.pa[lsarp.cazy.baseline.pa>0] = 1
lsarp.cazy.baseline = lsarp.cazy.baseline[,colSums(lsarp.cazy.baseline.pa) >= nrow(lsarp.cazy.baseline.pa)*0.2]
lsarp.cazy.baseline = lsarp.cazy.baseline %>% as.data.frame()
# filtered
# lsarp.cazy.baseline = log2(lsarp.cazy.baseline+min(lsarp.cazy.baseline[lsarp.cazy.baseline!=0])/2) %>% as.data.frame()
lsarp.cazy.baseline$standard.name = rownames(lsarp.cazy.baseline)
dim(lsarp.cazy.baseline) # 13 samples * 64 features

# MBX by 80% per response group
lsarp.mbx.baseline = lsarp.mbx.annotated.mat[rownames(lsarp.mbx.annotated.mat) %in% lsarp.cd.baseline.samples,]
lsarp.mbx.baseline[is.na(lsarp.mbx.baseline)] = 0
lsarp.mbx.baseline.pa = lsarp.mbx.baseline
lsarp.mbx.baseline.pa[lsarp.mbx.baseline.pa>0] = 1
lsarp.mbx.baseline = lsarp.mbx.baseline[,colSums(lsarp.mbx.baseline.pa) >= nrow(lsarp.mbx.baseline.pa)*0.8]
lsarp.mbx.baseline = lsarp.mbx.baseline %>% as.data.frame()
# filtered
# lsarp.mbx.baseline = log2(lsarp.mbx.baseline+min(lsarp.mbx.baseline[lsarp.mbx.baseline!=0])/2) %>% as.data.frame()
lsarp.mbx.baseline$standard.name = rownames(lsarp.mbx.baseline)
dim(lsarp.mbx.baseline) # 14 samples * 148 features


# No FFQ because most participants were missing baseline

# combine
lsarp.omics.baseline = merge(lsarp.asv.baseline,
                             lsarp.mgx.baseline, by="standard.name")
lsarp.omics.baseline = merge(lsarp.omics.baseline,
                             lsarp.kegg.baseline, by="standard.name")
lsarp.omics.baseline = merge(lsarp.omics.baseline,
                             lsarp.cog.baseline, by="standard.name")
lsarp.omics.baseline = merge(lsarp.omics.baseline,
                             lsarp.cazy.baseline, by="standard.name")
lsarp.omics.baseline = merge(lsarp.omics.baseline,
                             lsarp.mbx.baseline, by="standard.name")

dim(lsarp.omics.baseline) # 3373 features

pheatmap::pheatmap(lsarp.omics.baseline[c(1:10),c(2:10)])

# n = 23 samples in All; but more in others

data.sets = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "Multi-Omic")


# :: Wilcoxon -------------------------------------------------------------

# calculate wilcox
lsarp.rf.models.importances.wilcox = #do.call(rbind, lapply(lsarp.rf.models.importances.df$feature, function(x){
  do.call(rbind, lapply(colnames(lsarp.omics.baseline[,!colnames(lsarp.omics.baseline) %in% c("RS_Name", "flare.group", "standard.name")]), function(x){
    
    print(x)
    #x = lsarp.rf.models.importances.df$feature[1]
    # add meta
    new.data = lsarp.omics.baseline
    new.data = merge(new.data,
                     metadata.lsarp.stool.asv[,c("Group", "RS_Name", "standard.name","ave.fiber", "richness", "shannon", "but.i", "but.ii", "fd","load.asv","load.mgx","starch","mucin","starch.mucin")],
                     by="standard.name")
    new.data$HM = substr(new.data$standard.name, 1, 6)
    new.data$standard.name = NULL
    # subset to RS group
    new.data = subset(new.data, Group == "RS")
    # add flare.group
    new.data = merge(lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                     new.data,
                     by="HM")
    
    data.subset = new.data[,colnames(new.data) %in% c(x, "flare.group")]

    # for omics, log2-transform
    if(x %in% colnames(lsarp.omics.baseline)){
      omics.values = data.subset[,2]
      data.subset[,2] = log2(omics.values + min(omics.values[omics.values!=0])/2)
    }
    # run wilcox
    wilcox.p = wilcox.test(subset(data.subset, flare.group == "Relapse")[,2],
                           subset(data.subset, flare.group == "Remit")[,2])
    # run t-test
    ttest.p = t.test(subset(data.subset, flare.group == "Relapse")[,2],
                           subset(data.subset, flare.group == "Remit")[,2])
    
    coef = mean((subset(data.subset, flare.group == "Remit")[,2])) - mean((subset(data.subset, flare.group == "Relapse")[,2]))
    
    data.frame(feature = x,
               wilcox.pval = wilcox.p$p.value,
               ttest.pval = ttest.p$p.value,
               coef = coef)
  }))

# append feature type
lsarp.rf.models.importances.wilcox = lsarp.rf.models.importances.wilcox %>%
  mutate(data.type = ifelse(feature %in% colnames(lsarp.asv.baseline), "ASV",
                            ifelse(feature %in% colnames(lsarp.mgx.baseline), "Species",
                                   ifelse(feature %in% colnames(lsarp.kegg.baseline), "Pathway",
                                          ifelse(feature %in% colnames(lsarp.cog.baseline), "COG",
                                                 ifelse(feature %in% colnames(lsarp.cazy.baseline), "CAZy", 
                                                        ifelse(feature %in% colnames(lsarp.mbx.baseline), "Metabolite", "Other")))))))
lsarp.rf.models.importances.wilcox$wilcox.padj = p.adjust(lsarp.rf.models.importances.wilcox$wilcox.pval, method="BH")
lsarp.rf.models.importances.wilcox$ttest.padj = p.adjust(lsarp.rf.models.importances.wilcox$ttest.pval, method="BH")

lsarp.rf.models.importances.wilcox$data.type = factor(lsarp.rf.models.importances.wilcox$data.type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))

# select specific features to show
ml.feature.selected = c("Faecalibacterium_SGB15346",
                        "Uridine kinase",
                        "Ethanolamine ammonia-lyase, large subunit",
                        "g__Faecalibacterium_s__duncaniae",
                        "Propanediol dehydratase, large subunit",
                        "Sulfur metabolism")

# add FDR per omic
lsarp.rf.models.importances.wilcox = lsarp.rf.models.importances.wilcox %>%
  group_by(data.type) %>%
  mutate(wilcox.padj.per.omic = p.adjust(wilcox.pval, method="BH"))%>%
  mutate(ttest.padj.per.omic = p.adjust(ttest.pval, method="BH"))
  

lsarp.rf.models.importances.wilcox %>%
  arrange(ttest.pval)

# volcano plot
lsarp.rf.models.importances.wilcox.plot = ggplot(lsarp.rf.models.importances.wilcox,
        aes(x=coef, y=wilcox.pval))+
  geom_point(aes(alpha= ifelse(wilcox.pval < 0.05, "A", "B"), fill=data.type, shape=data.type))+
  geom_hline(yintercept=-log10(0.05), linetype=0.2, alpha=1)+
  scale_shape_manual(values=omics.shapes)+
  scale_fill_manual(values=omics.colors)+
  scale_alpha_manual(values=c(1,0.1), guide="none")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.005, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggnetwork::geom_nodetext_repel(
    data = subset(lsarp.rf.models.importances.wilcox,
                  feature %in% ml.feature.selected),
    aes(label= gsub("\\..*", "", feature)),
    size=2.5)+
  facet_wrap(~"Log2FC ~ Wilcoxon p value")+
  theme_classic()+theme(legend.position="bottom",
                        legend.direction="horizontal")+
  labs(x=expression(Log[2]*FC),
       y="p value",
       fill = "Data type", shape = "Data type",
       alpha=NULL)
lsarp.rf.models.importances.wilcox.plot

subset(lsarp.rf.models.importances.wilcox, wilcox.pval < 0.05)  %>% nrow()
subset(lsarp.rf.models.importances.wilcox, wilcox.pval < 0.05) %>%
  mutate(coef.sign = sign(coef)) %>% dplyr::select(coef.sign) %>% table() %>% colSums()
# 74 / 91 = 81%

# :: feature plots --------------------------------------------------------

features.to.keep = lsarp.rf.models.importances.wilcox %>% slice_min(order_by = wilcox.pval, n=20) %>% arrange(wilcox.pval)

# loop through top features and plot; extract microbiome values
lsarp.rf.models.importances.ttest = 
  do.call(rbind, lapply(features.to.keep$feature, function(x){
    print(x)
    # add meta
    new.data = lsarp.omics.baseline
    new.data = merge(new.data,
                     metadata.lsarp.stool.asv[,c("Group", "standard.name","ave.fiber", "richness", "shannon", "but.i", "but.ii", "fd","load.asv","load.mgx","starch","mucin","starch.mucin")],
                     by="standard.name")
    new.data = subset(new.data, Group == "RS")
    new.data$HM = substr(new.data$standard.name, 1, 6)
    new.data$standard.name = NULL
    new.data = merge(lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                     new.data,
                     by="HM")
    
    # extract values 
    data.subset = new.data[,c(x, "flare.group")]
    data.subset$feature = colnames(data.subset)[1]
    data.subset$sample = rownames(data.subset)
    colnames(data.subset)[1] = "value"
    data.subset
  }))
# add stats
lsarp.rf.models.importances.ttest = merge(lsarp.rf.models.importances.ttest,
                                          lsarp.rf.models.importances.wilcox, by="feature")
lsarp.rf.models.importances.ttest.p = lsarp.rf.models.importances.ttest[,c("feature", "wilcox.pval")] %>% distinct()

# now clip names to 30 char
lsarp.rf.models.importances.ttest$short.feature = ifelse(nchar(lsarp.rf.models.importances.ttest$feature)>30, paste(substr(lsarp.rf.models.importances.ttest$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.ttest$feature)
lsarp.rf.models.importances.ttest.p$short.feature = ifelse(nchar(lsarp.rf.models.importances.ttest.p$feature)>30, paste(substr(lsarp.rf.models.importances.ttest.p$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.ttest.p$feature)
#lsarp.rf.models.importances.df$short.feature = ifelse(nchar(lsarp.rf.models.importances.df$feature)>30, paste(substr(lsarp.rf.models.importances.df$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.df$feature)

lsarp.rf.models.importances.ttest$short.feature = gsub("shannon", "Shannon Diversity", lsarp.rf.models.importances.ttest$short.feature)
lsarp.rf.models.importances.ttest$short.feature = gsub("but.ii", "Kircher Butyrogens", lsarp.rf.models.importances.ttest$short.feature)
lsarp.rf.models.importances.ttest$short.feature = gsub("but.i", "Butyrogens", lsarp.rf.models.importances.ttest$short.feature)

lsarp.rf.models.importances.ttest.p$short.feature = gsub("shannon", "Shannon Diversity", lsarp.rf.models.importances.ttest.p$short.feature)
lsarp.rf.models.importances.ttest.p$short.feature = gsub("but.ii", "Kircher Butyrogens", lsarp.rf.models.importances.ttest.p$short.feature)
lsarp.rf.models.importances.ttest.p$short.feature = gsub("but.i", "Butyrogens", lsarp.rf.models.importances.ttest.p$short.feature)

unique(lsarp.rf.models.importances.ttest$feature)


lsarp.rf.models.importances.ttest.plot = ggplot(lsarp.rf.models.importances.ttest %>%
                                                  # subset to specific features
                                                  subset(feature %in% ml.feature.selected) %>% 
                                                  # clean names
                                                  mutate(Group = flare.group) %>%
                                                  # add indicator of not present
                                                  group_by(feature) %>%
                                                  mutate(pseudo = ifelse(value == min(value), "pseudo", "real")) %>%
                                                  # reorder taxa based on imp
                                                  mutate(short.feature = factor(short.feature, levels=arrange(distinct(subset(lsarp.rf.models.importances.ttest.p,
                                                                                                                              feature %in% ml.feature.selected)[,c("short.feature", "wilcox.pval")]),wilcox.pval)$short.feature)),
                                                aes(x=Group, y=value))+
  geom_boxplot(width=0.3, outlier.shape=NA)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  #scale_y_log10()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=Group, alpha=pseudo), size=2)+
  scale_alpha_manual(values=c(0.2, 1), guide="none")+
  theme_classic()+theme(legend.position="none",
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=8),
                        strip.background = element_rect(color="black"))+
  geom_text(data=lsarp.rf.models.importances.ttest.p %>%
              # subset to specific features
              subset(feature %in% ml.feature.selected) %>%
              # clean names
              mutate(short.feature = factor(short.feature, levels=arrange(distinct(subset(lsarp.rf.models.importances.ttest.p,
                                                                                          feature %in% ml.feature.selected)[,c("short.feature", "wilcox.pval")]),wilcox.pval)$short.feature)),
            x=1.5, y=Inf, vjust=1.2, 
            aes(label = paste("p =", round(wilcox.pval, digits=3))),
            size=2.5)+
  labs(x="", y="Feature Abundance")+
  facet_wrap(~short.feature, ncol=2, scales="free")
lsarp.rf.models.importances.ttest.plot



# :: RF -------------------------------------------------------------------

# Goal: use RandomSurvivalForest to predict relapse prior to, and including SOC scope

# Note: Only n=11 participants had baseline FFQ

data.sets = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "Multi-Omic")

lsarp.omics.baseline$standard.name %>% length()

# n = 13

# starts with baseline data, processed above
lsarp.rf.validation.output = rf.function(data.types = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "Multi-Omic"),
                                         iters = 15,
                                         reduce_features = F,
                                         # p.threshold = 0.10,
                                         output = "AUC")
# 25 min (or 1.5 h with reduce_features)                              

# summarize
lsarp.rf.validation.output.df = lsarp.rf.validation.output %>%
  group_by(data.type) %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         mean.acc = mean(accuracy),
         median.acc = median(accuracy),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high, mean.acc, median.acc, data.type) %>% distinct()
lsarp.rf.validation.output.df

lsarp.rf.validation.output.plot = lsarp.rf.validation.output.df %>%
  mutate(data.set = factor(data.type, levels=arrange(lsarp.rf.validation.output.df, mean.auc)$data.type)) %>%
  ggplot(aes(y=data.set, x=mean.auc))+
  geom_bar(stat="identity", width=0.75,
           aes(fill=data.type),color="black")+
  geom_segment(aes(x=auc.low, xend=auc.high))+
  scale_x_continuous(breaks=seq(0,1, by=0.1),
                     limits=c(0,1))+
  geom_text(aes(y=data.type, x=auc.high,  
                label=paste(" ", round(mean.auc, digits=2))), hjust=0, size=3)+
  scale_fill_manual(values=c(omics.colors, "Multi-Omic" = "black"))+
  #geom_segment()+
  facet_wrap(~"15x LOOCV RandomForest")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(y="", x="LOOCV AUC")
lsarp.rf.validation.output.plot


# :: RF ROC ---------------------------------------------

# select which RF to analyze

lsarp.rf.selected.output = rf.function(data.types = c("Multi-Omic"),
                                       iters = 15,
                                       balance = T, # need to balance, otherwise biased towards majority
                                       reduce_features = F,
                                       output = "ROC")

# calculate values and extract best
lsarp.rf.selected.output.df = lsarp.rf.selected.output %>%
  group_by(iter) %>%
  summarize(auc = pROC::auc(true, relapse, 
                            levels=c("Remit", "Relapse"),  # define case = "Relapse"
                            direction="<")[1])%>% # if probability > 0.5 (e.g.), then it is a case
            
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high) %>% distinct()

# extract best
lsarp.rf.selected.output.stats = subset(lsarp.rf.selected.output.df, mean.auc == max(lsarp.rf.selected.output.df$mean.auc))

# plot
lsarp.rf.selected.output.roc = do.call(rbind, lapply(1:15, function(seed){
  data.subset = subset(lsarp.rf.selected.output, iter == seed)
  # calibrate pred; not necessary
  #data.subset$pred = 1 / (1 + exp(-data.subset$pred))
  data.frame(sens = pROC::roc(response = data.subset$true, predictor = data.subset$relapse, levels=c("Remit", "Relapse"),  # define case = "Relapse"
                              direction="<")$sensitivities,
             spec = pROC::roc(response = data.subset$true, predictor = data.subset$relapse,levels=c("Remit", "Relapse"),  # define case = "Relapse"
                              direction="<")$specificities,
             iter = seed)
}))

# Step 2: Compute mean and standard error for sensitivities across iterations
lsarp.rf.selected.output.roc_summary <- lsarp.rf.selected.output.roc %>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


lsarp.rf.selected.output.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_smooth(data = lsarp.rf.selected.output.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter),
  #          se=T, linewidth=0.3, color="grey")+
  # Error ribbon (95% CI)
  geom_ribbon(data = lsarp.rf.selected.output.roc_summary, 
              aes(x = fpr, ymin = lower, ymax = upper, fill = "Multi-Omic"),
              alpha=0.2)+
  scale_fill_manual(values = omics.colors)+
  # Mean ROC curve
  geom_path(data = lsarp.rf.selected.output.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black", linewidth=1.5)+
  geom_path(data = lsarp.rf.selected.output.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="white", linewidth=1)+
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("LOOCV\nAUC: ", round(lsarp.rf.selected.output.stats$mean.auc, digits=2),
                       "\n(", round(lsarp.rf.selected.output.stats$auc.low, digits=2), 
                       ", ", round(lsarp.rf.selected.output.stats$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position="none")+
  facet_wrap(~"Multi-Omic RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
lsarp.rf.selected.output.roc.plot


# :: RF Feature Importances ----------------------------------------------------------

lsarp.rf.selected.importances.output = rf.function(data.types = c("Multi-Omic"),
                                       iters = 15,
                                       balance = T, # need to balance, otherwise biased towards majority
                                       reduce_features = F,
                                       output = "importances")

lsarp.rf.selected.importances.output.df = lsarp.rf.selected.importances.output %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(mean.imp = mean(na.omit(imp))) %>%
  mutate(imp.low = mean(imp)- (sd(imp)/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(imp)+ (sd(imp)/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high) %>% distinct() %>%
  arrange(-mean.imp) %>% as.data.frame()
  
# evaluate sign of these features
lsarp.rf.selected.importances.output.plot =  lsarp.omics.baseline[,colnames(lsarp.omics.baseline) %in% c("standard.name", (lsarp.rf.selected.importances.output.df %>%
                       slice_max(order_by=mean.imp, n = 20) %>% dplyr::select(feature) %>% unlist() %>% as.vector()))] %>%
  merge(lsarp.metadata.responders[,c("standard.name", "flare.group")],
        by="standard.name") %>%
  reshape2::melt() %>%
  group_by(variable) %>%
  # calculate lfc ( // defunct)
  mutate(coef = mean(value[flare.group == "Remit"]) / mean(value[flare.group == "Relapse"])) %>%
  mutate(coef = log2(coef))%>%
  dplyr::select(variable, coef) %>% distinct() %>%
  mutate(feature = variable) %>% dplyr::select(-variable)%>%
  merge(lsarp.rf.selected.importances.output.df,
        by="feature") %>%
  # add omic type
  mutate(data.type = ifelse(feature %in% colnames(lsarp.asv.baseline), "ASV",
                            ifelse(feature %in% colnames(lsarp.mgx.baseline), "Species",
                                   ifelse(feature %in% colnames(lsarp.kegg.baseline), "Pathway",
                                          ifelse(feature %in% colnames(lsarp.cog.baseline), "COG",
                                                 ifelse(feature %in% colnames(lsarp.cazy.baseline), "CAZy", 
                                                        ifelse(feature %in% colnames(lsarp.mbx.baseline), "Metabolite", "Other"))))))) %>%
  mutate(data.type = factor(data.type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))) %>%
  # bolden axis labels
  mutate(bolded = ifelse(feature %in% ml.feature.selected, "bold", "plain")) %>%
  # ugh
  mutate(bolded_feature = ifelse(bolded == "bold", paste0("<b>", feature, "</b>"), as.character(feature)))%>%
  mutate(bolded_feature = ifelse(bolded_feature == "<b>g__Faecalibacterium_s__duncaniae</b>",
                                 "<b>(s) Faecalibacterium_duncaniae</b>", bolded_feature)) %>%
  mutate(feature = ifelse(feature == "g__Faecalibacterium_s__duncaniae",
                                 "(s) Faecalibacterium_duncaniae", feature)) %>%
  # plot
  ggplot(aes(x=reorder(bolded_feature, mean.imp), y=mean.imp))+
  geom_segment(aes(x=reorder(bolded_feature, mean.imp), y=imp.low, yend=imp.high))+
  geom_point(aes(shape = data.type, fill = data.type), size=3)+
  scale_shape_manual(values=omics.shapes)+
  scale_fill_manual(values=omics.colors)+
  #scale_fill_gradient2(low="blue", high="red")+
  coord_flip()+
  theme_classic()+
  theme(legend.position=c(0.8,0.35),
        strip.text = element_text(size=10),
        axis.text.y = ggtext::element_markdown()) + 
  facet_wrap(~"Feature Importance")+
  labs(y="Mean Decrease in Accuracy", x=NULL, 
       fill="Omic type", shape = "Omic type")
lsarp.rf.selected.importances.output.plot


# :: PCA ----------------------------------------------------

# Use 20% prevalence filtered

# baseline samples
lsarp.baseline.samples = subset(lsarp.metadata.responders,
                                baseline == "baseline" & Group == "RS")$standard.name

## ASV (glommed)

lsarp.lmer.asv.bray = vegan::vegdist(lsarp.asv.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples)%>%dplyr::select(-standard.name), method="bray") 
# perform PCoA
lsarp.lmer.asv.pcoa = ape::pcoa(lsarp.lmer.asv.bray)
# extract data from pcoa
lsarp.lmer.asv.pcoa.df = data.frame(lsarp.lmer.asv.pcoa$vectors[,c(1:2)])
lsarp.lmer.asv.pcoa.df$standard.name = rownames(lsarp.lmer.asv.pcoa.df)
# add metadata
lsarp.lmer.asv.pcoa.df = merge(lsarp.lmer.asv.pcoa.df,
                               lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.asv.pcoa.var_exp = lsarp.lmer.asv.pcoa$values[c(1:2),2]
lsarp.lmer.asv.pcoa.df$var1 = round(lsarp.lmer.asv.pcoa.var_exp[1]*100, digits=2)
lsarp.lmer.asv.pcoa.df$var2 = round(lsarp.lmer.asv.pcoa.var_exp[2]*100, digits=2)
rownames(lsarp.lmer.asv.pcoa.df) = lsarp.lmer.asv.pcoa.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.asv.pcoa.permanova = vegan::adonis2(lsarp.lmer.asv.bray ~ flare.group,
                                          lsarp.lmer.asv.pcoa.df,
                                          by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.asv.pcoa.permanova # not-sig

lsarp.lmer.asv.pcoa.plot <- ggplot(
  data=lsarp.lmer.asv.pcoa.df,
  aes(x=Axis.1, y=Axis.2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.asv.pcoa.permanova)[1,3], 3)*100, "%",
                       "  p: ", round(data.frame(lsarp.lmer.asv.pcoa.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("Axis 1: ", round(unique(lsarp.lmer.asv.pcoa.var_exp[1])*100, digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.lmer.asv.pcoa.var_exp[2])*100, digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"ASV")
lsarp.lmer.asv.pcoa.plot


### MGX

lsarp.lmer.mgx.bray = vegan::vegdist(lsarp.mgx.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples)%>%dplyr::select(-standard.name), method="bray") 
# perform PCoA
lsarp.lmer.mgx.pcoa = ape::pcoa(lsarp.lmer.mgx.bray)
# extract data from pcoa
lsarp.lmer.mgx.pcoa.df = data.frame(lsarp.lmer.mgx.pcoa$vectors[,c(1:2)])
lsarp.lmer.mgx.pcoa.df$standard.name = rownames(lsarp.lmer.mgx.pcoa.df)
# add metadata
lsarp.lmer.mgx.pcoa.df = merge(lsarp.lmer.mgx.pcoa.df,
                               lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.mgx.pcoa.var_exp = lsarp.lmer.mgx.pcoa$values[c(1:2),2]
lsarp.lmer.mgx.pcoa.df$var1 = round(lsarp.lmer.mgx.pcoa.var_exp[1]*100, digits=2)
lsarp.lmer.mgx.pcoa.df$var2 = round(lsarp.lmer.mgx.pcoa.var_exp[2]*100, digits=2)
rownames(lsarp.lmer.mgx.pcoa.df) = lsarp.lmer.mgx.pcoa.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.mgx.pcoa.permanova = vegan::adonis2((lsarp.lmer.mgx.pcoa$vectors[,c(1,2)]) ~ flare.group,
                                               lsarp.lmer.mgx.pcoa.df,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.mgx.pcoa.permanova # not-sig

lsarp.lmer.mgx.pcoa.plot <- ggplot(
  data=lsarp.lmer.mgx.pcoa.df,
  aes(x=Axis.1, y=Axis.2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.mgx.pcoa.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(lsarp.lmer.mgx.pcoa.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("Axis 1: ", round(unique(lsarp.lmer.mgx.pcoa.var_exp[1])*100, digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(lsarp.lmer.mgx.pcoa.var_exp[2])*100, digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Species")
lsarp.lmer.mgx.pcoa.plot


### Pathway

# Perform PCA
lsarp.lmer.kegg.pca = lsarp.kegg.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples)%>%dplyr::select(-standard.name)
lsarp.lmer.kegg.pca.pseudo = min(lsarp.lmer.kegg.pca[lsarp.lmer.kegg.pca!=0])/2
lsarp.lmer.kegg.pca = prcomp(log2(lsarp.lmer.kegg.pca+lsarp.lmer.kegg.pca.pseudo))

# extract data from pcoa
lsarp.lmer.kegg.pca.df = data.frame(lsarp.lmer.kegg.pca$x[,c(1:2)])
lsarp.lmer.kegg.pca.df$standard.name = rownames(lsarp.lmer.kegg.pca.df)
# add metadata
lsarp.lmer.kegg.pca.df = merge(lsarp.lmer.kegg.pca.df,
                               lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.kegg.pca.var_exp = (lsarp.lmer.kegg.pca$sdev)^2 / sum(lsarp.lmer.kegg.pca$sdev^2) * 100

rownames(lsarp.lmer.kegg.pca.df) = lsarp.lmer.kegg.pca.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.kegg.pca.permanova = vegan::adonis2(dist(lsarp.lmer.kegg.pca$x) ~ flare.group,
                                               lsarp.lmer.kegg.pca.df,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.kegg.pca.permanova # not-sig

lsarp.lmer.kegg.pca.plot <- ggplot(
  data=lsarp.lmer.kegg.pca.df,
  aes(x=PC1, y=PC2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.kegg.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(lsarp.lmer.kegg.pca.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("PC 1: ", round(unique(lsarp.lmer.kegg.pca.var_exp[1]), digits=2), "%", sep=""), 
       y=paste("PC 2: ", round(unique(lsarp.lmer.kegg.pca.var_exp[2]), digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Pathway")
lsarp.lmer.kegg.pca.plot


### COG

# 
lsarp.lmer.cog.pca = lsarp.cog.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples) %>% dplyr::select(-standard.name)
lsarp.lmer.cog.pca.pseudo = min(lsarp.lmer.cog.pca[lsarp.lmer.cog.pca!=0])/2
lsarp.lmer.cog.pca = prcomp(log2(lsarp.lmer.cog.pca+lsarp.lmer.cog.pca.pseudo))

# extract data from pcoa
lsarp.lmer.cog.pca.df = data.frame(lsarp.lmer.cog.pca$x[,c(1:2)])
lsarp.lmer.cog.pca.df$standard.name = rownames(lsarp.lmer.cog.pca.df)
# add metadata
lsarp.lmer.cog.pca.df = merge(lsarp.lmer.cog.pca.df,
                              lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.cog.pca.var_exp = (lsarp.lmer.cog.pca$sdev)^2 / sum(lsarp.lmer.cog.pca$sdev^2) * 100

rownames(lsarp.lmer.cog.pca.df) = lsarp.lmer.cog.pca.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.cog.pca.permanova = vegan::adonis2(dist(lsarp.lmer.cog.pca$x) ~ flare.group,
                                               lsarp.lmer.cog.pca.df,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.cog.pca.permanova # not-sig

lsarp.lmer.cog.pca.plot <- ggplot(
  data=lsarp.lmer.cog.pca.df,
  aes(x=PC1, y=PC2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.cog.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(lsarp.lmer.cog.pca.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("PC 1: ", round(unique(lsarp.lmer.cog.pca.var_exp[1]), digits=2), "%", sep=""), 
       y=paste("PC 2: ", round(unique(lsarp.lmer.cog.pca.var_exp[2]), digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"COG")
lsarp.lmer.cog.pca.plot

### CAZy

# 
lsarp.lmer.cazy.pca = lsarp.cazy.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples) %>% dplyr::select(-standard.name)
lsarp.lmer.cazy.pca.pseudo = min(lsarp.lmer.cazy.pca[lsarp.lmer.cazy.pca!=0])/2
lsarp.lmer.cazy.pca = prcomp(log2(lsarp.lmer.cazy.pca+lsarp.lmer.cazy.pca.pseudo))

# extract data from pcoa
lsarp.lmer.cazy.pca.df = data.frame(lsarp.lmer.cazy.pca$x[,c(1:2)])
lsarp.lmer.cazy.pca.df$standard.name = rownames(lsarp.lmer.cazy.pca.df)
# add metadata
lsarp.lmer.cazy.pca.df = merge(lsarp.lmer.cazy.pca.df,
                               lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.cazy.pca.var_exp = (lsarp.lmer.cazy.pca$sdev)^2 / sum(lsarp.lmer.cazy.pca$sdev^2) * 100

rownames(lsarp.lmer.cazy.pca.df) = lsarp.lmer.cazy.pca.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.cazy.pca.permanova = vegan::adonis2(dist(lsarp.lmer.cazy.pca$x) ~ flare.group,
                                               lsarp.lmer.cazy.pca.df,
                                               by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.cazy.pca.permanova # not-sig

lsarp.lmer.cazy.pca.plot <- ggplot(
  data=lsarp.lmer.cazy.pca.df,
  aes(x=PC1, y=PC2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.cazy.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(lsarp.lmer.cazy.pca.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("PC 1: ", round(unique(lsarp.lmer.cazy.pca.var_exp[1]), digits=2), "%", sep=""), 
       y=paste("PC 2: ", round(unique(lsarp.lmer.cazy.pca.var_exp[2]), digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"CAZy")
lsarp.lmer.cazy.pca.plot

### MBX

# 
lsarp.lmer.mbx.pca = lsarp.mbx.baseline%>%subset(standard.name %in% lsarp.cd.baseline.samples) %>% dplyr::select(-standard.name)
lsarp.lmer.mbx.pca.pseudo = min(lsarp.lmer.mbx.pca[lsarp.lmer.mbx.pca!=0])/2
lsarp.lmer.mbx.pca = prcomp(log2(lsarp.lmer.mbx.pca+lsarp.lmer.mbx.pca.pseudo))

# extract data from pcoa
lsarp.lmer.mbx.pca.df = data.frame(lsarp.lmer.mbx.pca$x[,c(1:2)])
lsarp.lmer.mbx.pca.df$standard.name = rownames(lsarp.lmer.mbx.pca.df)
# add metadata
lsarp.lmer.mbx.pca.df = merge(lsarp.lmer.mbx.pca.df,
                              lsarp.metadata.responders[,c("standard.name","flare.group")], by="standard.name")
# extract variance explained
lsarp.lmer.mbx.pca.var_exp = (lsarp.lmer.mbx.pca$sdev)^2 / sum(lsarp.lmer.mbx.pca$sdev^2) * 100

rownames(lsarp.lmer.mbx.pca.df) = lsarp.lmer.mbx.pca.df$standard.name

set.seed(25)
t1 = Sys.time()
lsarp.lmer.mbx.pca.permanova = vegan::adonis2(dist(lsarp.lmer.mbx.pca$x) ~ flare.group,
                                              lsarp.lmer.mbx.pca.df,
                                              by="margin")
t2 = Sys.time()
t2 - t1
lsarp.lmer.mbx.pca.permanova # not-sig

lsarp.lmer.mbx.pca.plot <- ggplot(
  data=lsarp.lmer.mbx.pca.df,
  aes(x=PC1, y=PC2))+
  stat_ellipse(aes(group=flare.group, color=flare.group), alpha=0.5)+
  geom_point(aes(fill=flare.group),shape=21, color="white", size=2.5)+
  annotate(geom="text", x=0, y=Inf, hjust=0.5, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(lsarp.lmer.mbx.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(lsarp.lmer.mbx.pca.permanova)[1,5], 3), sep=""), sep=""))+
  labs(x=paste("PC 1: ", round(unique(lsarp.lmer.mbx.pca.var_exp[1]), digits=2), "%", sep=""), 
       y=paste("PC 2: ", round(unique(lsarp.lmer.mbx.pca.var_exp[2]), digits=2), "%", sep=""))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Metabolite")
lsarp.lmer.mbx.pca.plot


rf.pca.p1 = ((lsarp.lmer.asv.pcoa.plot+
         lsarp.lmer.mgx.pcoa.plot+
         lsarp.lmer.kegg.pca.plot+
         lsarp.lmer.cog.pca.plot+
         lsarp.lmer.cazy.pca.plot+
           lsarp.lmer.mbx.pca.plot)+patchwork::plot_layout(nrow=1))

rf.p1 = (rf.pca.p1/
           (lsarp.rf.models.importances.wilcox.plot+lsarp.rf.models.importances.ttest.plot+facet_wrap(~short.feature, nrow=2, scales="free"))+
           patchwork::plot_layout(heights=c(1,3)))

rf.p2 = lsarp.rf.validation.output.plot/patchwork::free(lsarp.rf.selected.output.roc.plot, type="label")
rf.p2 = lsarp.rf.validation.output.plot/(lsarp.rf.selected.output.roc.plot)

(rf.p1 | rf.p2) + patchwork::plot_layout(widths=c(2,1))


# :: Change over time -----------------------------------------------------

# Are the important features impacted by RS?

lsarp.rf.models.importances.ttest.p$feature

# combine all stats results for RS responder vs non-responder
lsarp.rf.feature.change.with.rs = rbind(
      lsarp.lmer.asv.resp.interactions,
      lsarp.lmer.mgx.resp.interactions,
      lsarp.lmer.mpx.kegg.resp.interactions,
      lsarp.lmer.mpx.cog.resp.interactions,
      lsarp.lmer.mpx.cazy.resp.interactions,
      lsarp.lmer.mbx.resp.interactions) %>%
  as.data.frame() %>%
  subset(group == "RS") %>%
  subset(taxa %in% lsarp.rf.models.importances.ttest.p$feature) %>%
  arrange(padj)
lsarp.rf.feature.change.with.rs
# ATP-dependent protease decreased
# Lachnospiraceae increased
# F. prausnitzii increased
# Propanediol dehydratase increasede

# add omic type
lsarp.rf.feature.change.with.rs = lsarp.rf.feature.change.with.rs %>%
  mutate(data.type = ifelse(taxa %in% colnames(lsarp.asv.baseline), "ASV",
                            ifelse(taxa %in% colnames(lsarp.mgx.baseline), "Species",
                                   ifelse(taxa %in% colnames(lsarp.kegg.baseline), "Pathway",
                                          ifelse(taxa %in% colnames(lsarp.cog.baseline), "COG",
                                                 ifelse(taxa %in% colnames(lsarp.cazy.baseline), "CAZy", 
                                                        ifelse(taxa %in% colnames(lsarp.mbx.baseline), "Metabolite", "Other")))))))
# plot
lsarp.rf.feature.change.with.rs.volcano = ggplot(lsarp.rf.feature.change.with.rs %>%
                                                   mutate(taxa = ifelse(taxa == "g__Faecalibacterium_s__duncaniae",
                                                                           "(s) Faecalibacterium_duncaniae", taxa)),
       aes(x=estimate, y=padj))+
  geom_point(aes(fill=data.type, shape=data.type))+
  scale_shape_manual(values=omics.shapes)+
  scale_fill_manual(values=omics.colors)+
  geom_hline(yintercept = 0.20, linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", taxa), NA)),
                           size=2.5)+
  facet_wrap(~"Interaction p value (Responders vs Non-Responders)")+
  theme_classic()+theme(legend.position="bottom",
                        legend.direction="horizontal")+
  labs(x="Interaction",
       y="FDR",
       fill = "Data type", shape = "Data type")
lsarp.rf.feature.change.with.rs.volcano


# plot the actual change over time

lsarp.ml.predictors.over.time = lsarp.asv.data.glom[,colnames(lsarp.asv.data.glom) %in% ml.feature.selected]

# merge
lsarp.ml.predictors.over.time = 
  rownames_to_column(as.data.frame(lsarp.asv.data.glom/50000), var= "standard.name") %>%
  full_join(., rownames_to_column(as.data.frame(lsarp.mgx.taxa/100), var= "standard.name")) %>%
  full_join(., rownames_to_column(as.data.frame(lsarp.cd.mpx.kegg.mat), var= "standard.name")) %>%
  full_join(., rownames_to_column(as.data.frame(lsarp.cd.mpx.cog.mat), var= "standard.name")) %>%
  full_join(., rownames_to_column(as.data.frame(lsarp.cd.mpx.cazy.mat), var= "standard.name")) %>%
  full_join(., rownames_to_column(as.data.frame(lsarp.cd.mpx.cazy.mat), var= "standard.name")) %>%
  dplyr::select(standard.name, ml.feature.selected) %>%
  reshape2::melt()
# add patient data
lsarp.ml.predictors.over.time = 
  merge(lsarp.ml.predictors.over.time,
        lsarp.metadata.responders.asv, by="standard.name")
# add pseudocount and log scale
lsarp.ml.predictors.over.time = lsarp.ml.predictors.over.time %>%
  group_by(variable) %>%
  mutate(value = ifelse(is.na(value), 0, value))%>%
  mutate(value = (value + min(value[value!=0]/2)))
# now plot trends

# get p values from interaction
lsarp.rf.feature.change.with.rs %>%
  subset(taxa %in% lsarp.ml.predictors.over.time$variable)

lsarp.ml.predictors.over.time.plot = ggplot(lsarp.ml.predictors.over.time%>%
         subset(flare!="flare"& !HM %in% excluded.by.dave) %>%
         subset(Group == "RS") %>%
         # remove sample with missing data (i.e. keep only those used in RF)
         subset(HM %in% substr(lsarp.omics.baseline$standard.name, 1, 6)) %>%
           mutate(variable = ifelse(as.character(variable) == "g__Faecalibacterium_s__duncaniae",
                                    "(s) Faecalibacterium_duncaniae", as.character(variable))),
       aes(x=lsarp.days, y=(na.omit(value))))+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  scale_y_log10()+
  #geom_point(shape=21, aes(fill = flare.group))+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(data = lsarp.ml.predictors.over.time%>%
               subset(flare!="flare"& !HM %in% excluded.by.dave) %>%
               subset(Group == "RS") %>%
               # remove sample with missing data (i.e. keep only those used in RF)
               subset(HM %in% substr(lsarp.omics.baseline$standard.name, 1, 6)) %>%
               subset(lsarp.days < 0) %>%
               mutate(variable = ifelse(as.character(variable) == "g__Faecalibacterium_s__duncaniae",
                                        "(s) Faecalibacterium_duncaniae", as.character(variable))),
             shape=21, size=2.5,
             aes(fill = flare.group, x = lsarp.days, y=value))+
  # add lm p value
  geom_text(
    data = lsarp.rf.feature.change.with.rs %>%
      subset(taxa %in% lsarp.ml.predictors.over.time$variable) %>%
      mutate(variable = taxa) %>%
      mutate(variable = ifelse(as.character(variable) == "g__Faecalibacterium_s__duncaniae",
                               "(s) Faecalibacterium_duncaniae", as.character(variable))),
    aes(y=Inf, label = paste("p:", round(pval, digits=3))), 
    x=55, vjust=2, hjust=0)+
  # add baseline wilcox p value
  geom_text(
    data = lsarp.rf.models.importances.wilcox %>%
      subset(feature %in% lsarp.ml.predictors.over.time$variable) %>%
      mutate(variable = feature) %>%
      mutate(variable = ifelse(variable == "g__Faecalibacterium_s__duncaniae",
                               "(s) Faecalibacterium_duncaniae", variable)),
    aes(y=Inf, label = paste("p:", round(wilcox.pval, digits=3))), 
    x=0, vjust=2, hjust=0.5)+
  facet_wrap(~variable, scales="free", nrow=2)+
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Days on product", y="Feature Abundance")
lsarp.ml.predictors.over.time.plot

# Generally, responders start with lower, and increase with RS

# >>> 7. RapidAIM Predictors ----------------------------------------------

# Preliminary ex vivo data analysis is performed in 
# 2025_09_20_lsarp_rapidaim_test.R
# In there, we observe that responders had a significantly higher increase
# in Phasco/Dialister (using OTU data) versus non-responders

# Let's see if Phasco/Dialister were positively impacted in vivo


# :: ASV Phasco/Dialister -------------------------------------------------------

# calculate log2 ratio of Phasco / Dialister
lsarp.metadata.responders.asv.phasdial = data.frame(phasdial = log2((lsarp.lmer.asv.data[,grep("Phasco", colnames(lsarp.lmer.asv.data))])) - 
  log2(rowSums(lsarp.lmer.asv.data[,grep("Dialister", colnames(lsarp.lmer.asv.data))])),
  standard.name = rownames(lsarp.lmer.asv.data))
lsarp.metadata.responders.asv.phasdial = merge(lsarp.metadata.responders.asv.phasdial,
                                               lsarp.metadata.responders.asv, by="standard.name")

# stats
lsarp.resp.stats.phasdial.treatment = lmerTest::lmer(scale((phasdial)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                  subset(lsarp.metadata.responders.asv.phasdial, group %in% c("RS")) %>%
                                                    subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.phasdial.placebo = lmerTest::lmer(scale(phasdial) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                subset(lsarp.metadata.responders.asv.phasdial, group %in% c("Placebo")) %>%
                                                  subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.phasdial.resp.plot = ggplot(lsarp.metadata.responders.asv.phasdial %>%
                                          subset(!is.na(phasdial))%>%
                                          subset(flare !="flare"& !HM %in% excluded.by.dave),
                                        aes(x=lsarp.days, 
                                            y=phasdial))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=3)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Phascolarctobacterium/Dialister (ASV)",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.phasdial.treatment[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.phasdial.placebo[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.phasdial.resp.plot



# :: MGX Phasco/Dialister -------------------------------------------------------

lsarp.metadata.responders.mgx.phasdial = data.frame(phasdial = log2((lsarp.lmer.mgx.data[,grep("Phasco", colnames(lsarp.lmer.mgx.data))])) - 
                                                      log2((lsarp.lmer.mgx.data[,grep("Dialister", colnames(lsarp.lmer.mgx.data))])),
                                                    standard.name = rownames(lsarp.lmer.mgx.data))
lsarp.metadata.responders.mgx.phasdial = merge(lsarp.metadata.responders.mgx.phasdial,
                                               lsarp.metadata.responders.asv, by="standard.name")

# stats
lsarp.resp.stats.mgx.phasdial.treatment = lmerTest::lmer(scale((phasdial)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                     subset(lsarp.metadata.responders.mgx.phasdial, group %in% c("RS")) %>%
                                                       subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.mgx.phasdial.placebo = lmerTest::lmer(scale(phasdial) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                   subset(lsarp.metadata.responders.mgx.phasdial, group %in% c("Placebo")) %>%
                                                     subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.mgx.phasdial.resp.plot = ggplot(lsarp.metadata.responders.mgx.phasdial %>%
                                             subset(!is.na(phasdial))%>%
                                             subset(flare !="flare"& !HM %in% excluded.by.dave),
                                           aes(x=lsarp.days, 
                                               y=phasdial))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=3)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Phascolarctobacterium/Dialister (Metagenomics)",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.mgx.phasdial.treatment[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.mgx.phasdial.placebo[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.mgx.phasdial.resp.plot

# not sig

# :: MPX Succinate -------------------------------------------------------

# which proteins (COGs) are involved in succinate metabolism?
# Propionyl CoA:succinate CoA transferase
# Succinate dehydrogenase/fumarate reductase, cytochrome b subunit
# Succinate dehydrogenase/fumarate reductase, Fe-S protein subunit
# Succinate dehydrogenase/fumarate reductase, flavoprotein subunit
# Succinyl-CoA synthetase, alpha subunit
# Succinyl-CoA synthetase, beta subunit

lsarp.metadata.responders.mpx.succinate = data.frame(succ = log10(rowSums(lsarp.lmer.mpx.cog.data[,colnames(lsarp.lmer.mpx.cog.data) %in%
                                                                                                   c("Propionyl CoA:succinate CoA transferase",
                                                                                                     "Succinate dehydrogenase/fumarate reductase, cytochrome b subunit",
                                                                                                     "Succinate dehydrogenase/fumarate reductase, Fe-S protein subunit",
                                                                                                     "Succinate dehydrogenase/fumarate reductase, flavoprotein subunit",
                                                                                                     "Succinyl-CoA synthetase, alpha subunit",
                                                                                                     "Succinyl-CoA synthetase, beta subunit")])),
                                                    standard.name = rownames(lsarp.lmer.mpx.cog.data))
lsarp.metadata.responders.mpx.succinate = merge(lsarp.metadata.responders.mpx.succinate,
                                               lsarp.metadata.responders.asv, by="standard.name")

# stats
lsarp.resp.stats.mpx.succinate.treatment = lmerTest::lmer(scale((succ)) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                         subset(lsarp.metadata.responders.mpx.succinate, group %in% c("RS")) %>%
                                                           subset(flare != "flare"& !HM %in% excluded.by.dave)) %>%
  summary() %>% coef()

lsarp.resp.stats.mpx.succinate.placebo = lmerTest::lmer(scale(succ) ~ scale(lsarp.days)*flare.group + ave.fiber + (1|HM),
                                                       subset(lsarp.metadata.responders.mpx.succinate, group %in% c("Placebo")) %>%
                                                         subset(flare != "flare"& !HM %in% excluded.by.dave))%>%
  summary() %>% coef()


# plot
metadata.lsarp.mpx.succinate.resp.plot = ggplot(lsarp.metadata.responders.mpx.succinate %>%
                                                 subset(!is.na(succ))%>%
                                                 subset(flare !="flare"& !HM %in% excluded.by.dave),
                                               aes(x=lsarp.days, 
                                                   y=succ))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(method="lm", se=T,aes(fill=flare.group, color=flare.group))+
  geom_smooth(method="lm", se=F, aes(fill=flare.group), color="white", size=2)+
  geom_smooth(method="lm", se=F, aes(group=flare.group, color=flare.group), size=1)+
  geom_point(aes(fill=flare.group, shape=baseline), color="white", size=3)+
  scale_shape_manual(values=c(23,21))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Sum(Succinate-related COGs)",
       title=paste(paste("Placebo p:", round(lsarp.resp.stats.mpx.succinate.placebo[5,5], 3)),
                   paste("   RS p:", round(lsarp.resp.stats.mpx.succinate.treatment[5,5], 3))))+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
metadata.lsarp.mpx.succinate.resp.plot

# not sig

(metadata.lsarp.phasdial.resp.plot+
    metadata.lsarp.mgx.phasdial.resp.plot+
    metadata.lsarp.mpx.succinate.resp.plot) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_phasco_invivo.pdf",
         width=18, height=4,device = cairo_pdf)


# >>> FIGURES -------------------------------------------------------------



# :: Main Figures ---------------------------------------------------------


# statistics are adjusted for average fiber intake

(metadata.lsarp.stool.plot +
  redcap.lsarp.survival.plot$plot +
   patchwork::plot_layout(widths=c(2,1.2))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_stools_survival.pdf",
         width=14, height=5,device = cairo_pdf)

(lsarp.cd.rs.selections.placebo.plot+
    lsarp.cd.rs.selections.rs.plot+ 
    lsarp.cd.rs.selections.frequences.plot+
    patchwork::plot_layout(widths=c(2.3,3,2))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_rs_selections.pdf",
         width=16, height=8, device = cairo_pdf)


((metadata.lsarp.but.i.plot+metadata.lsarp.fcal.plot)/
    (metadata.lsarp.stool.asv.richness.plot+
    metadata.lsarp.stool.asv.shannon.plot+
    metadata.lsarp.stool.asv.fd.plot+
      metadata.lsarp.stool.asv.fr.plot+
      metadata.lsarp.stool.asv.between.plot+
      metadata.lsarp.water.plot+
      metadata.lsarp.load.plot+
      metadata.lsarp.starch.plot+
      metadata.lsarp.mucin.plot+
      metadata.lsarp.starch.mucin.plot+patchwork::plot_layout(nrow=4))+
  patchwork::plot_layout(heights=c(1,3.5))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_group_standard.pdf",
         width=14, height=15, device = cairo_pdf)

# MULTIOMICS
(lsarp.asv.pcoa.plot+
    lsarp.mgx.pcoa.plot+
    lsarp.kegg.lsarp.pca.plot+
    lsarp.cog.lsarp.pca.plot+
    lsarp.cazy.lsarp.pca.plot+
    lsarp.mbx.lsarp.pca.plot+
    patchwork::plot_layout(nrow=3))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_group_pca.pdf",
         width=12, height=10, device = cairo_pdf)


# heatmaps

pheatmap::pheatmap(lsarp.lmer.omics.lfc.treatment.mat.for.pheatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = lsarp.lmer.omics.lfc.treatment.sample.map %>% dplyr::select(-HM),
                   annotation_colors = list(Group = c(`RS` = gg_color_hue(2)[1],
                                                      `Placebo` = gg_color_hue(2)[2]),
                                            `interaction` =  colorRampPalette(c("blue","white", "red"))(100),
                                            `data.type` = omics.colors,
                                            `RS_Name` = rs.colors),
                   annotation_row = lsarp.lmer.omics.lfc.treatment.feature.map,
                   breaks=c(seq(min(na.omit(lsarp.lmer.omics.lfc.treatment.mat)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(lsarp.lmer.omics.lfc.treatment.mat))/100, max(na.omit(lsarp.lmer.omics.lfc.treatment.mat)), length.out=floor(100/2))))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_1_group_heatmap.pdf",
         width=12, height=5, device = cairo_pdf)


redcap.lsarp.wpcdai.6m.plot %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_2_response_rates.pdf",
         width=2, height=2, device = cairo_pdf)

((metadata.lsarp.but.i.resp.plot+metadata.lsarp.fcal.resp.plot)/
   (metadata.lsarp.stool.asv.richness.resp.plot+
   metadata.lsarp.stool.asv.shannon.resp.plot+
     metadata.lsarp.stool.asv.fr.resp.plot+
     metadata.lsarp.stool.asv.fd.resp.plot+
    metadata.lsarp.between.beta.resp.plot+
     metadata.lsarp.water.resp.plot+
     metadata.lsarp.load.resp.plot+
     metadata.lsarp.starch.resp.plot+
     metadata.lsarp.mucin.resp.plot+
     metadata.lsarp.starch.mucin.resp.plot+patchwork::plot_layout(nrow=4))+
  patchwork::plot_layout(heights=c(1,3.5))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_2_response_standard.pdf",
         width=14, height=15, device = cairo_pdf)


# resp.heatmaps

pheatmap::pheatmap(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap,
                   #scale = "row",
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = lsarp.lmer.omics.lfc.resp.treatment.sample.map %>% dplyr::select(-HM),
                   annotation_colors = list(Response = c(`Relapse` = gg_color_hue(2)[1],
                                                         `Remit` = gg_color_hue(2)[2]),
                                            Group = c(`RS` = gg_color_hue(2)[1],
                                                      `Placebo` = gg_color_hue(2)[2]),
                                            `Interaction` =  colorRampPalette(c("blue","white", "red"))(100),
                                            `Data type` = omics.colors,
                                            `RS_Name` = rs.colors),
                   annotation_row = lsarp.lmer.omics.lfc.resp.treatment.feature.map,
                   breaks=c(seq(min(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap))/100, max(na.omit(lsarp.lmer.omics.lfc.resp.treatment.mat.for.pheatmap)), length.out=floor(100/2))))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_2_response_omics_heatmap.pdf",
         width=18, height=9, device = cairo_pdf)


# v2
lsarp.rf.selected.importances.output.plot

rf.pca.p1 = ((lsarp.lmer.asv.pcoa.plot+
                lsarp.lmer.mgx.pcoa.plot+
                lsarp.lmer.kegg.pca.plot+
                lsarp.lmer.cog.pca.plot+
                lsarp.lmer.cazy.pca.plot+
                lsarp.lmer.mbx.pca.plot)+patchwork::plot_layout(nrow=1))


rf.p2 = lsarp.rf.validation.output.plot/patchwork::free(lsarp.rf.selected.output.roc.plot, type="label")

lsarp.rf.models.importances.ttest.plot


(((patchwork::free(rf.pca.p1, type="label")/
  (patchwork::free(rf.p2, type="label") |
  lsarp.rf.selected.importances.output.plot))/
  (lsarp.ml.predictors.over.time.plot))+
  patchwork::plot_layout(heights=c(0.75,3,2))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_2_response_predictors.pdf",
         width=13, height=15, device = cairo_pdf)

# :: Supplementals --------------------------------------------------------


lsarp.cd.ph.rs.boxplot%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_fermentation_responses.pdf",
         width=8, height=4, device = cairo_pdf)


(lsarp.lmer.asv.interactions.plot+
    lsarp.lmer.mgx.interactions.plot+
    lsarp.lmer.mpx.kegg.interactions.plot+
    lsarp.lmer.mpx.cog.interactions.plot+
    lsarp.lmer.mpx.cazy.interactions.plot+
    lsarp.lmer.mbx.interactions.plot+
    patchwork::plot_layout(nrow=3))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_group_interaction_volcano.pdf",
         width=12, height=10, device = cairo_pdf)


(lsarp.lmer.asv.resp.interactions.plot+
    lsarp.lmer.mgx.resp.interactions.plot+
    lsarp.lmer.mpx.kegg.resp.interactions.plot+
    lsarp.lmer.mpx.cog.resp.interactions.plot+
    lsarp.lmer.mpx.cazy.resp.interactions.plot+
    lsarp.lmer.mbx.resp.interactions.plot+
    patchwork::plot_layout(nrow=3))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_responders_interaction_volcano.pdf",
         width=12, height=10, device = cairo_pdf)

# fecal cal correlates
(lsarp.fecalcal.omics.cor.plot+
  lsarp.fecalcal.omics.cor.plots+patchwork::plot_layout(widths=c(1,2))) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_responders_fecal_cor_plots.pdf",
         width=20, height=6, device = cairo_pdf)

# predictor change over time

(lsarp.rf.feature.change.with.rs.volcano) %>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_supp_predictors_change_over_RS.pdf",
         width=5, height=5, device = cairo_pdf)


# :: Graphical Abstract ---------------------------------------------------

# see other file

# :: ----------------------------------------------------------------------


# >> Extra ----------------------------------------------------------------


# :: check numbers --------------------------------------------------------

## placebo vs RS

# ASV numbers 12 v 14
lsarp.asv.pcoa.df[,c("HM", "Group")] %>% table()
lsarp.asv.pcoa.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()

# MGX numbers 12 v 14
lsarp.mgx.resp.pcoa.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% dplyr::select(HM, flare.group) %>% table()
lsarp.mgx.resp.pcoa.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% dplyr::select(HM, flare.group) %>% distinct() %>% dplyr::select(flare.group) %>% table()

# MPX numbers 12 v 14
lsarp.cazy.lsarp.pca.df[,c("HM", "Group")] %>% table()
lsarp.cazy.lsarp.pca.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()

# MBX numbers 12 v 14
lsarp.mbx.lsarp.pca.df[,c("HM", "Group")] %>% table()
lsarp.mbx.lsarp.pca.df[,c("HM", "Group")] %>% distinct() %>% dplyr::select(Group) %>% table()



## responders vs non-responders

# ASV numbers 5 v 9, up to 6
lsarp.asv.resp.pcoa.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% dplyr::select(HM, flare.group) %>% table()
lsarp.asv.resp.pcoa.df[,c("HM", "Group", "flare.group")] %>% subset(Group == "RS") %>% dplyr::select(HM, flare.group) %>% distinct() %>% dplyr::select(flare.group) %>% table()

# MGX numbers 5 v 9, up to 6
lsarp.lmer.mgx.data[rownames(lsarp.lmer.mgx.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% table()
lsarp.lmer.mgx.data[rownames(lsarp.lmer.mgx.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% distinct() %>% dplyr::select(flare.group) %>% table()

# MPX numbers 5 v 9, up to 6
lsarp.lmer.mpx.cog.data[rownames(lsarp.lmer.mpx.cog.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% table()
lsarp.lmer.mpx.cog.data[rownames(lsarp.lmer.mpx.cog.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% distinct() %>% dplyr::select(flare.group) %>% table()

# MBX numbers 5 v 9, up to 6
lsarp.lmer.mbx.data[rownames(lsarp.lmer.mbx.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% table()
lsarp.lmer.mbx.data[rownames(lsarp.lmer.mbx.data) %in% stools.not.flare,] %>%
  merge(lsarp.metadata.responders.asv, by="standard.name") %>% 
  subset(group == "RS") %>%
  dplyr::select(HM, flare.group) %>% distinct() %>% dplyr::select(flare.group) %>% table()


# :: Spearman Cluster ASV -----------------------------------------------------

# convert to clr
lsarp.cd.asv.data.glom.clr = lsarp.asv.data.glom
lsarp.cd.asv.data.glom.clr.pa = lsarp.asv.data.glom
lsarp.cd.asv.data.glom.clr.pa[lsarp.cd.asv.data.glom.clr.pa > 0] = 1
lsarp.cd.asv.data.glom.clr = compositions::clr(lsarp.cd.asv.data.glom.clr) %>% as.data.frame()
lsarp.cd.asv.data.glom.clr[lsarp.cd.asv.data.glom.clr.pa == 0] = min(lsarp.cd.asv.data.glom.clr[lsarp.cd.asv.data.glom.clr!=0])-1

lsarp.cd.asv.data.glom.clr.dist = as.dist(1-Hmisc::rcorr(t(lsarp.cd.asv.data.glom.clr), type="spearman")$r)
lsarp.cd.asv.data.glom.clr.dist.tree <- ape::as.phylo(hclust(lsarp.cd.asv.data.glom.clr.dist, method = "complete"))
lsarp.cd.asv.data.glom.clr.dist.tree$tip.label <- rownames(as.matrix(lsarp.cd.asv.data.glom.clr.dist))
library("ggtree")
lsarp.cd.asv.data.glom.clr.dist.tree.plot <- ggtree(lsarp.cd.asv.data.glom.clr.dist.tree) %<+% metadata.lsarp.stool +
  geom_tiplab(size = 2) +
  #coord_flip()+
  geom_tippoint(aes(fill = Group), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(tree$edge.length) *1.4)+  # widen x-axis beyond tree height
  facet_wrap(~"LSARP ASV")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
lsarp.cd.asv.data.glom.clr.dist.tree.plot


# :: Spearman Cluster MGX -----------------------------------------------------

# convert to clr
lsarp.cd.mgx.data.clr = lsarp.mgx.taxa
lsarp.cd.mgx.data.clr.pa = lsarp.mgx.taxa
lsarp.cd.mgx.data.clr.pa[lsarp.cd.mgx.data.clr.pa > 0] = 1
lsarp.cd.mgx.data.clr = compositions::clr(lsarp.cd.mgx.data.clr) %>% as.data.frame()
lsarp.cd.mgx.data.clr[lsarp.cd.mgx.data.clr.pa == 0] = min(lsarp.cd.mgx.data.clr[lsarp.cd.mgx.data.clr!=0])-1

lsarp.cd.mgx.data.clr.dist = as.dist(1-Hmisc::rcorr(t(lsarp.cd.mgx.data.clr), type="spearman")$r)
lsarp.cd.mgx.data.clr.dist.tree <- ape::as.phylo(hclust(lsarp.cd.mgx.data.clr.dist, method = "complete"))
lsarp.cd.mgx.data.clr.dist.tree$tip.label <- rownames(as.matrix(lsarp.cd.mgx.data.clr.dist))
library("ggtree")
lsarp.cd.mgx.data.clr.dist.tree.plot <- ggtree(lsarp.cd.mgx.data.clr.dist.tree) %<+% metadata.lsarp.stool +
  geom_tiplab(size = 2) +
  #coord_flip()+
  geom_tippoint(aes(fill = Group), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(lsarp.cd.mgx.data.clr.dist.tree$edge.length)*2.5)+  # widen x-axis beyond tree height
  facet_wrap(~"LSARP MGX")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
lsarp.cd.mgx.data.clr.dist.tree.plot


# :: Spearman Cluster MPX -----------------------------------------------------

lsarp.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_lsarp_mpx_cog.Rds")
#lsarp.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_lsarp_mpx_cazy.Rds")

# convert to log
lsarp.cd.mpx.data.log = lsarp.mpx.cazy.mat
lsarp.cd.mpx.data.log[is.na(lsarp.cd.mpx.data.log)] = 0
lsarp.cd.mpx.data.log.pseudo = min(lsarp.cd.mpx.data.log[lsarp.cd.mpx.data.log!=0])/2
lsarp.cd.mpx.data.log = log2(lsarp.cd.mpx.data.log+lsarp.cd.mpx.data.log.pseudo)
lsarp.cd.mpx.data.log[lsarp.cd.mpx.data.log == 0] = min(lsarp.cd.mpx.data.log[lsarp.cd.mpx.data.log!=0])-1

lsarp.cd.mpx.data.log.dist = as.dist(1-Hmisc::rcorr(t(lsarp.cd.mpx.data.log), type="spearman")$r)
lsarp.cd.mpx.data.log.dist.tree <- ape::as.phylo(hclust(lsarp.cd.mpx.data.log.dist, method = "complete"))
lsarp.cd.mpx.data.log.dist.tree$tip.label <- rownames(as.matrix(lsarp.cd.mpx.data.log.dist))
library("ggtree")
lsarp.cd.mpx.data.log.dist.tree.plot <- ggtree(lsarp.cd.mpx.data.log.dist.tree) %<+% metadata.lsarp.stool +
  geom_tiplab(size = 2) +
  #coord_flip()+
  geom_tippoint(aes(fill = Group), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(lsarp.cd.mpx.data.log.dist.tree$edge.length) *1.2)+  # widen x-axis beyond tree height
  facet_wrap(~"LSARP MGX")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
lsarp.cd.mpx.data.log.dist.tree.plot


lsarp.cd.asv.data.glom.log.plot+
  lsarp.cd.asv.mgx.log.plot+
  lsarp.cd.asv.mpx.log.plot+
  patchwork::plot_layout(guides = "collect")


# :: ----------------------------------------------------------------------


# :: ----------------------------------------------------------------------


# >>> Old -----------------------------------------------------------------



# :: ASV Maaslin2 // defunct ---------------------------------------------------------

# Defunct because Placebo and RS have different sample sizes
# which confounds the difference in p values
# hence, use interactions

# asv or glom: lsarp.asv.data.median | lsarp.asv.data.glom
rownames(metadata.lsarp.stool.asv) = metadata.lsarp.stool.asv$standard.name

lsarp.asv.maaslin.treatment.rs = Maaslin2::Maaslin2(input_data = lsarp.asv.data.glom,
                                                    input_metadata = subset(metadata.lsarp.stool.asv, 
                                                                            phase == "treatment" & group == "RS"),
                                                    output = "~/Downloads",
                                                    fixed_effects = c("lsarp.days"),  # Example fixed effects
                                                    random_effects = c("HM"),       # Example random effects
                                                    normalization = "TSS",                       # Total Sum Scaling normalization
                                                    transform = "LOG",                           # Log transformation
                                                    analysis_method = "LM",                      # Linear model
                                                    plot_scatter = FALSE,                        # Disable scatterplot generation
                                                    plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                    max_significance = 0.05,                     # Significance threshold for q-values
                                                    standardize = TRUE                           # Disable standardization (optional)
)

lsarp.asv.maaslin.treatment.plac = Maaslin2::Maaslin2(input_data = lsarp.asv.data.glom,
                                                      input_metadata = subset(metadata.lsarp.stool.asv, 
                                                                              phase == "treatment" & group == "Placebo"),
                                                      output = "~/Downloads",
                                                      fixed_effects = c("lsarp.days"),  # Example fixed effects
                                                      random_effects = c("HM"),       # Example random effects
                                                      normalization = "TSS",                       # Total Sum Scaling normalization
                                                      transform = "LOG",                           # Log transformation
                                                      analysis_method = "LM",                      # Linear model
                                                      plot_scatter = FALSE,                        # Disable scatterplot generation
                                                      plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                      max_significance = 0.05,                     # Significance threshold for q-values
                                                      standardize = TRUE                           # Disable standardization (optional)
)

lsarp.asv.maaslin.treatment.rs = lsarp.asv.maaslin.treatment.rs$results %>% as.data.frame() %>% arrange(pval)
lsarp.asv.maaslin.treatment.rs$group = "RS"
lsarp.asv.maaslin.treatment.plac = lsarp.asv.maaslin.treatment.plac$results %>% as.data.frame() %>% arrange(pval)
lsarp.asv.maaslin.treatment.plac$group = "Placebo"

lsarp.asv.maaslin.treatment = rbind(lsarp.asv.maaslin.treatment.rs,
                                    lsarp.asv.maaslin.treatment.plac) %>% as.data.frame()

# recalculate padj minus diagnosis fixed effect
lsarp.asv.maaslin.treatment = lsarp.asv.maaslin.treatment %>%
  subset(metadata == "lsarp.days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

lsarp.asv.maaslin.treatment.volcano = ggplot(lsarp.asv.maaslin.treatment,
                                             aes(x=coef, y=padj))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20 & abs(coef) > 0.7, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  labs(x="Adjusted Coefficient",
       title="Treatment Phase (ASV)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.asv.maaslin.treatment.volcano
# Veillonella is sig depleted over RS
# are these "more significant" (lower p val) because there is more power with RS?
(metadata.lsarp.stool[,c("HM", "group")])$group %>% table()
# most likely, which is probably why this section is defunct

lsarp.asv.maaslin.washout.rs = Maaslin2::Maaslin2(input_data = lsarp.asv.data.glom,
                                                  input_metadata = subset(subset(metadata.lsarp.stool.asv, group=="RS"), phase == "washout" | timing == "5M"),
                                                  output = "~/Downloads",
                                                  fixed_effects = c("lsarp.days"),  # Example fixed effects
                                                  random_effects = c("HM"),       # Example random effects
                                                  normalization = "TSS",                       # Total Sum Scaling normalization
                                                  transform = "LOG",                           # Log transformation
                                                  analysis_method = "LM",                      # Linear model
                                                  plot_scatter = FALSE,                        # Disable scatterplot generation
                                                  plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                  max_significance = 0.05,                     # Significance threshold for q-values
                                                  standardize = TRUE                           # Disable standardization (optional)
)
lsarp.asv.maaslin.washout.plac = Maaslin2::Maaslin2(input_data = lsarp.asv.data.glom,
                                                    input_metadata = subset(subset(metadata.lsarp.stool.asv, group=="Placebo"), phase == "washout" | timing == "5M"),
                                                    output = "~/Downloads",
                                                    fixed_effects = c("lsarp.days"),  # Example fixed effects
                                                    random_effects = c("HM"),       # Example random effects
                                                    normalization = "TSS",                       # Total Sum Scaling normalization
                                                    transform = "LOG",                           # Log transformation
                                                    analysis_method = "LM",                      # Linear model
                                                    plot_scatter = FALSE,                        # Disable scatterplot generation
                                                    plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                    max_significance = 0.05,                     # Significance threshold for q-values
                                                    standardize = TRUE                           # Disable standardization (optional)
)

lsarp.asv.maaslin.washout.rs = lsarp.asv.maaslin.washout.rs$results %>% as.data.frame() %>% arrange(pval)
lsarp.asv.maaslin.washout.rs$group = "RS"
lsarp.asv.maaslin.washout.plac = lsarp.asv.maaslin.washout.plac$results %>% as.data.frame() %>% arrange(pval)
lsarp.asv.maaslin.washout.plac$group = "Placebo"

lsarp.asv.maaslin.washout = rbind(lsarp.asv.maaslin.washout.rs,
                                  lsarp.asv.maaslin.washout.plac) %>% as.data.frame()

# recalculate padj minus diagnosis fixed effect
lsarp.asv.maaslin.washout = lsarp.asv.maaslin.washout %>% 
  subset(metadata == "lsarp.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

lsarp.asv.maaslin.washout.volcano = ggplot(lsarp.asv.maaslin.washout,
                                           aes(x=coef, y=padj))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20 & abs(coef) > 0.5, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  labs(x="Adjusted Coefficient",
       title="Washout Phase (ASV)",
       y="FDR")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.asv.maaslin.washout.volcano
# nothing sig over washout
# certainly because underpowered


# :: MGX Maaslin2 // defunct ---------------------------------------------------------

# 
rownames(metadata.lsarp.stool) = metadata.lsarp.stool$standard.name
lsarp.mgx.maaslin.treatment.rs = Maaslin2::Maaslin2(input_data = lsarp.mgx.taxa,
                                                    input_metadata = subset(metadata.lsarp.stool, phase == "treatment" & group == "RS"),
                                                    output = "~/Downloads",
                                                    fixed_effects = c("lsarp.days"),  # Example fixed effects
                                                    random_effects = c("HM"),       # Example random effects
                                                    normalization = "TSS",                       # Total Sum Scaling normalization
                                                    transform = "LOG",                           # Log transformation
                                                    analysis_method = "LM",                      # Linear model
                                                    plot_scatter = FALSE,                        # Disable scatterplot generation
                                                    plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                    max_significance = 0.05,                     # Significance threshold for q-values
                                                    standardize = TRUE                           # Disable standardization (optional)
)
lsarp.mgx.maaslin.treatment.plac = Maaslin2::Maaslin2(input_data = lsarp.mgx.taxa,
                                                      input_metadata = subset(metadata.lsarp.stool, phase == "treatment" & group == "Plac"),
                                                      output = "~/Downloads",
                                                      fixed_effects = c("lsarp.days", "ave.fiber"),  # Example fixed effects
                                                      random_effects = c("HM"),       # Example random effects
                                                      normalization = "TSS",                       # Total Sum Scaling normalization
                                                      transform = "LOG",                           # Log transformation
                                                      analysis_method = "LM",                      # Linear model
                                                      plot_scatter = FALSE,                        # Disable scatterplot generation
                                                      plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                      max_significance = 0.05,                     # Significance threshold for q-values
                                                      standardize = TRUE                           # Disable standardization (optional)
)

lsarp.mgx.maaslin.treatment.rs = lsarp.mgx.maaslin.treatment.rs$results %>% as.data.frame() %>% arrange(pval)
lsarp.mgx.maaslin.treatment.rs$group = "RS"
lsarp.mgx.maaslin.treatment.plac = lsarp.mgx.maaslin.treatment.plac$results %>% as.data.frame() %>% arrange(pval)
lsarp.mgx.maaslin.treatment.plac$group = "Plac"

lsarp.mgx.maaslin.treatment = rbind(lsarp.mgx.maaslin.treatment.rs,
                                    lsarp.mgx.maaslin.treatment.plac) %>% as.data.frame()

# recalculate padj minus diagnosis fixed effect
lsarp.mgx.maaslin.treatment = lsarp.mgx.maaslin.treatment %>%
  subset(metadata == "lsarp.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

lsarp.mgx.maaslin.treatment.volcano = ggplot(lsarp.mgx.maaslin.treatment,
                                             aes(x=coef, y=-log(pval)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.05 & abs(coef) > 0.5, gsub("\\..*", "", feature), NA)),
                           size=2.5)+  
  labs(x="Adjusted Coefficient",
       title="MGX Treatment")+
  facet_wrap(~group,labeller = label_both)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.mgx.maaslin.treatment.volcano
# 
lsarp.mgx.maaslin.treatment %>% arrange(padj)


# washout
lsarp.mgx.maaslin.washout.rs = Maaslin2::Maaslin2(input_data = lsarp.mgx.taxa,
                                                  input_metadata = subset(subset(metadata.lsarp.stool,group =="RS"), phase == "washout" | timing == "5M"),
                                                  output = "~/Downloads",
                                                  fixed_effects = c("lsarp.days", "ave.fiber"),  # Example fixed effects
                                                  random_effects = c("HM"),       # Example random effects
                                                  normalization = "TSS",                       # Total Sum Scaling normalization
                                                  transform = "LOG",                           # Log transformation
                                                  analysis_method = "LM",                      # Linear model
                                                  plot_scatter = FALSE,                        # Disable scatterplot generation
                                                  plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                  max_significance = 0.05,                     # Significance threshold for q-values
                                                  standardize = TRUE                           # Disable standardization (optional)
)
lsarp.mgx.maaslin.washout.plac = Maaslin2::Maaslin2(input_data = lsarp.mgx.taxa,
                                                    input_metadata = subset(subset(metadata.lsarp.stool,group =="Plac"), phase == "washout" | timing == "5M"),
                                                    output = "~/Downloads",
                                                    fixed_effects = c("lsarp.days", "ave.fiber"),  # Example fixed effects
                                                    random_effects = c("HM"),       # Example random effects
                                                    normalization = "TSS",                       # Total Sum Scaling normalization
                                                    transform = "LOG",                           # Log transformation
                                                    analysis_method = "LM",                      # Linear model
                                                    plot_scatter = FALSE,                        # Disable scatterplot generation
                                                    plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                    max_significance = 0.05,                     # Significance threshold for q-values
                                                    standardize = TRUE                           # Disable standardization (optional)
)

lsarp.mgx.maaslin.washout.rs = lsarp.mgx.maaslin.washout.rs$results %>% as.data.frame() %>% arrange(pval)
lsarp.mgx.maaslin.washout.rs$group = "RS"
lsarp.mgx.maaslin.washout.plac = lsarp.mgx.maaslin.washout.plac$results %>% as.data.frame() %>% arrange(pval)
lsarp.mgx.maaslin.washout.plac$group = "Plac"

lsarp.mgx.maaslin.washout = rbind(lsarp.mgx.maaslin.washout.rs,
                                  lsarp.mgx.maaslin.washout.plac) %>% as.data.frame()

# recalculate padj minus diagnosis fixed effect
lsarp.mgx.maaslin.washout = lsarp.mgx.maaslin.washout %>%
  subset(metadata == "lsarp.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

lsarp.mgx.maaslin.washout.volcano = ggplot(lsarp.mgx.maaslin.washout,
                                           aes(x=coef, y=-log(pval)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.05 & abs(coef) > 0.5, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  labs(x="Adjusted Coefficient",
       title="MGX Washout")+
  facet_wrap(~group,labeller = label_both)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.mgx.maaslin.washout.volcano
# 

lsarp.mgx.maaslin.washout %>% arrange(padj)
# nothing sig



# :: MPX CAZy Maaslin2 // defunct ----------------------------------------------------------

lsarp.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_lsarp_mpx_cazy.Rds")

lsarp.mpx.cazy.mat[is.na(lsarp.mpx.cazy.mat)] = 0

rownames(lsarp.mpx.cazy.mat) = gsub("_", "-", rownames(lsarp.mpx.cazy.mat))

rownames(metadata.lsarp.stool) = metadata.lsarp.stool$standard.name
# cazy or cog

lsarp.mpx.cazy.maaslin.treatment.rs = Maaslin2::Maaslin2(input_data = lsarp.mpx.cazy.mat,
                                                         input_metadata = subset(metadata.lsarp.stool, 
                                                                                 phase == "treatment" & group == "RS"),
                                                         output = "~/Downloads",
                                                         fixed_effects = c("lsarp.days","ph_response", "ave.fiber"),  # Example fixed effects
                                                         random_effects = c("HM"),       # Example random effects
                                                         normalization = "NONE",                       # Total Sum Scaling normalization
                                                         transform = "LOG",                           # Log transformation
                                                         analysis_method = "LM",                      # Linear model
                                                         plot_scatter = FALSE,                        # Disable scatterplot generation
                                                         plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                         max_significance = 0.05,                     # Significance threshold for q-values
                                                         standardize = TRUE                           # Disable standardization (optional)
)


lsarp.mpx.cazy.maaslin.treatment.placebo = Maaslin2::Maaslin2(input_data = lsarp.mpx.cazy.mat,
                                                              input_metadata = subset(metadata.lsarp.stool, 
                                                                                      phase == "treatment" & group == "Plac"),
                                                              output = "~/Downloads",
                                                              fixed_effects = c("lsarp.days","ph_response", "ave.fiber"),  # Example fixed effects
                                                              random_effects = c("HM"),       # Example random effects
                                                              normalization = "NONE",                       # Total Sum Scaling normalization
                                                              transform = "LOG",                           # Log transformation
                                                              analysis_method = "LM",                      # Linear model
                                                              plot_scatter = FALSE,                        # Disable scatterplot generation
                                                              plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                              max_significance = 0.05,                     # Significance threshold for q-values
                                                              standardize = TRUE                           # Disable standardization (optional)
)


lsarp.mpx.cazy.maaslin.treatment.rs = lsarp.mpx.cazy.maaslin.treatment.rs$results %>% as.data.frame() %>% arrange(pval)
lsarp.mpx.cazy.maaslin.treatment.rs$group = "RS"
lsarp.mpx.cazy.maaslin.treatment.placebo = lsarp.mpx.cazy.maaslin.treatment.placebo$results %>% as.data.frame() %>% arrange(pval)
lsarp.mpx.cazy.maaslin.treatment.placebo$group = "Plac"

lsarp.mpx.cazy.maaslin.treatment = rbind(lsarp.mpx.cazy.maaslin.treatment.rs,
                                         lsarp.mpx.cazy.maaslin.treatment.placebo) %>% as.data.frame()

# recalculate padj minus diagnosis fixed effect
lsarp.mpx.cazy.maaslin.treatment = lsarp.mpx.cazy.maaslin.treatment %>% 
  subset(metadata == "lsarp.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

lsarp.mpx.cazy.maaslin.treatment.volcano = ggplot(lsarp.mpx.cazy.maaslin.treatment,
                                                  aes(x=coef, y=-log(pval)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  labs(x="Adjusted Coefficient",
       title="MPX (CAZyme) Treatment")+
  facet_wrap(~group,labeller = label_both, scales="free")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
lsarp.mpx.cazy.maaslin.treatment.volcano
# nothing sig over treatment


# :: Heatmaps (Responders) ----------------------------------------------------

# GOAL: Make Log2FC heatmap of top pval < 0.05 features in Resp vs Non-resp

# Subset data to nominally sig taxa

# features to include
features.to.include = metadata.lsarp.stool.omics.cor.bh %>%
  subset(padj < 0.20) %>%
  arrange(abs(value)) # do not filter beyond this; filter later
features.to.include = features.to.include$feature

# ASV Taxa
lsarp.asv.lfc = merge(lsarp.lmer.asv.data,
                      lsarp.metadata.responders.asv[,c("standard.name","HM", "phase","timing", "flare.group")],
                      by="standard.name") %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(lsarp.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% features.to.include)
#%>% slice_min(order_by=abs(estimate), n=thresh) %>% dplyr::select(taxa) %>% as.vector() %>% unlist()
lsarp.asv.lfc.map = lsarp.asv.lfc[,c("sample", "flare.group")] %>% as.data.frame() %>%distinct()
rownames(lsarp.asv.lfc.map) = lsarp.asv.lfc.map$sample
lsarp.asv.lfc.map$sample = NULL

# MGX Taxa
lsarp.mgx.lfc = merge(lsarp.lmer.mgx.data,
                      lsarp.metadata.responders.asv[,c("standard.name","HM", "phase","timing", "flare.group")],
                      by="standard.name") %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(lsarp.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% features.to.include)
#%>% slice_min(order_by=abs(estimate), n=thresh) %>% dplyr::select(taxa) %>% as.vector() %>% unlist()
lsarp.mgx.lfc.map = lsarp.mgx.lfc[,c("sample", "flare.group")] %>% as.data.frame() %>%distinct()
rownames(lsarp.mgx.lfc.map) = lsarp.mgx.lfc.map$sample
lsarp.mgx.lfc.map$sample = NULL

# MPX KEGG
lsarp.mpx.kegg.lfc = merge(lsarp.lmer.mpx.kegg.data,
                           lsarp.metadata.responders.asv[,c("standard.name","HM", "phase","timing", "flare.group")],
                           by="standard.name") %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(lsarp.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% make.names(features.to.include))
#%>% slice_min(order_by=abs(estimate), n=thresh) %>% dplyr::select(taxa) %>% as.vector() %>% unlist()
# perfectly correlated; just take one
lsarp.mpx.kegg.lfc.map = lsarp.mpx.kegg.lfc[,c("sample", "flare.group")] %>% as.data.frame() %>%distinct()
rownames(lsarp.mpx.kegg.lfc.map) = lsarp.mpx.kegg.lfc.map$sample
lsarp.mpx.kegg.lfc.map$sample = NULL

# MPX COG
lsarp.mpx.cog.lfc = merge(lsarp.lmer.mpx.cog.data,
                          lsarp.metadata.responders.asv[,c("standard.name","HM", "phase","timing", "flare.group")],
                          by="standard.name") %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(lsarp.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% make.names(features.to.include))
#%>% slice_min(order_by=abs(estimate), n=thresh) %>% dplyr::select(taxa) %>% as.vector() %>% unlist()
# perfectly correlated; just take one
lsarp.mpx.cog.lfc.map = lsarp.mpx.cog.lfc[,c("sample", "flare.group")] %>% as.data.frame() %>%distinct()
rownames(lsarp.mpx.cog.lfc.map) = lsarp.mpx.cog.lfc.map$sample
lsarp.mpx.cog.lfc.map$sample = NULL

# MPX CAZy
lsarp.mpx.cazy.lfc = merge(lsarp.lmer.mpx.cazy.data,
                           lsarp.metadata.responders.asv[,c("standard.name","HM", "phase","timing", "flare.group")],
                           by="standard.name") %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(lsarp.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% make.names(features.to.include))
#%>% slice_min(order_by=abs(estimate), n=thresh) %>% dplyr::select(taxa) %>% as.vector() %>% unlist()
# perfectly correlated; just take one
lsarp.mpx.cazy.lfc.map = lsarp.mpx.cazy.lfc[,c("sample", "flare.group")] %>% as.data.frame() %>%distinct()
rownames(lsarp.mpx.cazy.lfc.map) = lsarp.mpx.cazy.lfc.map$sample
lsarp.mpx.cazy.lfc.map$sample = NULL


# merge
lsarp.asv.mgx.mpx.lfc = rbind(
  lsarp.asv.lfc,
  lsarp.mgx.lfc,
  lsarp.mpx.kegg.lfc,
  lsarp.mpx.cog.lfc,
  lsarp.mpx.cazy.lfc) %>% as.data.frame()

# scale per variable
lsarp.asv.mgx.mpx.lfc = lsarp.asv.mgx.mpx.lfc %>% group_by(variable) %>% mutate(slfc = scale(lfc))
# add feature map
lsarp.asv.mgx.mpx.lfc.map = lsarp.asv.mgx.mpx.lfc[,c("variable")] %>% distinct() %>%
  mutate(datatype = ifelse(variable %in% lsarp.asv.lfc$variable, "ASV",
                           mutate(datatype = ifelse(variable %in% lsarp.mgx.lfc$variable, "Species",
                                                    ifelse(variable %in% lsarp.mpx.kegg.lfc$variable, "Pathway", 
                                                           ifelse(variable %in% lsarp.mpx.cog.lfc$variable, "COG", 
                                                                  ifelse(variable %in% lsarp.mpx.cazy.lfc$variable, "CAZy", "other"))))))) %>% as.data.frame()
rownames(lsarp.asv.mgx.mpx.lfc.map) = lsarp.asv.mgx.mpx.lfc.map$variable
colnames(lsarp.asv.mgx.mpx.lfc.map)[1] = "feature"

rownames(lsarp.asv.mgx.mpx.lfc.map) = lsarp.asv.mgx.mpx.lfc.map$feature
lsarp.asv.mgx.mpx.lfc.map$feature = NULL
colnames(lsarp.asv.mgx.mpx.lfc.map)[1] = "Data type"

# refresh response names
colnames(lsarp.mpx.cog.lfc.map) = "Response"

# refresh variable names
lsarp.asv.mgx.mpx.lfc.feature.map = data.frame(variable = make.names(features.to.include),
                                               clean.name = features.to.include)
lsarp.asv.mgx.mpx.lfc$feature = lsarp.asv.mgx.mpx.lfc.feature.map[match(lsarp.asv.mgx.mpx.lfc$variable,
                                                                        lsarp.asv.mgx.mpx.lfc.feature.map$variable),]$clean.name
rownames(lsarp.asv.mgx.mpx.lfc.map) = lsarp.asv.mgx.mpx.lfc.feature.map[match(rownames(lsarp.asv.mgx.mpx.lfc.map),
                                                                              lsarp.asv.mgx.mpx.lfc.feature.map$variable),]$clean.name
lsarp.asv.mgx.mpx.lfc.map$`Fcal Correlation` = metadata.lsarp.stool.omics.cor.bh[match(rownames(lsarp.asv.mgx.mpx.lfc.map),
                                                                                       metadata.lsarp.stool.omics.cor.bh$Var2),]$value
# shorten
rownames(lsarp.asv.mgx.mpx.lfc.map) = ifelse(nchar(as.character(rownames(lsarp.asv.mgx.mpx.lfc.map)))>40, paste(substr(rownames(lsarp.asv.mgx.mpx.lfc.map), 1, 40), "...", sep=""), as.character(rownames(lsarp.asv.mgx.mpx.lfc.map)))
lsarp.asv.mgx.mpx.lfc$feature = ifelse(nchar(as.character(lsarp.asv.mgx.mpx.lfc$feature))>40, paste(substr(lsarp.asv.mgx.mpx.lfc$feature, 1, 40), "...", sep=""), as.character(lsarp.asv.mgx.mpx.lfc$feature))

# clean
lsarp.asv.mgx.mpx.lfc.map$feature = NULL

# :: PLS-DA -----------------------------------------------------------------


# Step 1: Prepare data
# Assume lfc_matrix is your data (rows = samples, columns = features)
# Example: lfc_matrix <- metadata.lsarp.stool.omics.cor[, numeric_cols]
# labels is a factor vector of "strong" vs "weak"
# Replace with your actual data

lfc_matrix <- reshape2::acast(lsarp.asv.mgx.mpx.lfc,
                              sample ~ feature, value.var="lfc") %>% as.matrix()
# keep only complete data (missing ~50% of proteomics!)
lfc_matrix = lfc_matrix[rownames(lsarp.mpx.cog.lfc.map),]
labels <- factor(lsarp.mpx.cog.lfc.map[rownames(lfc_matrix),])  # Replace with your "strong" vs "weak" column


# Check for NAs and impute if needed
lfc_matrix[is.na(lfc_matrix)] <- 0  # Simple imputation; consider mixOmics::nipals for better handling

# Step 2: Run PLS-DA
plsda_model <- mixOmics::plsda(X = lfc_matrix, 
                               Y = labels, 
                               ncomp = 5)  # Try 5 components

# Tune number of components using cross-validation
perf_plsda <- mixOmics::perf(plsda_model, 
                             #validation = "Mfold", 
                             nrepeat = 10,
                             folds = 5, progressBar = TRUE)

ggplot(perf_plsda$error.rate.class$mahalanobis.dist%>%
         data.frame() %>%
         summarize(total.err = colSums(.)) %>%
         mutate(ncomp = seq(1:5)))+
  aes(x=ncomp, y=total.err*100)+
  geom_line(color="red", linewidth=1)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Method: Mahalanobis Distance")+
  labs(x="Number of Components",
       y="Total Error (%)")

ncomp_opt = 2

# Refit with optimal components
plsda_model <- mixOmics::plsda(X = lfc_matrix, Y = labels, ncomp = ncomp_opt)

# Step 4: Extract discriminative features (VIP scores)
vip_scores <- mixOmics::vip(plsda_model)  # VIP for each component
vip_df <- data.frame(
  Feature = colnames(lfc_matrix),
  VIP = vip_scores[, ncomp_opt]  # Use last component
) %>%
  arrange(desc(VIP))

# Filter top discriminative features (e.g., VIP > 1)
top_features <- vip_df %>%
  filter(VIP > 1) %>%
  head(20)  # Top 20 features

# Step 5: Visualize
# Score plot (sample separation)
mixOmics::plotIndiv(plsda_model, comp = c(1, 2), group = labels, legend = TRUE,
                    title = "PLS-DA Score Plot (Strong vs Weak)", ellipse = TRUE)
# plot manually

plsda_model_plot = plsda_model$variates$X[,c(1,2)] %>% as.data.frame() %>%
  merge(lsarp.mpx.cog.lfc.map, by="row.names") %>%
  ggplot(aes(x=comp1, y=comp2))+
  stat_ellipse(aes(color=Response), linewidth=1.5, alpha=0.5)+
  scale_color_discrete(guide = "none")+
  geom_point(aes(fill=Response), shape=21, size=3)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position=c(0.8, 0.15),
                        legend.title = element_blank(),
                        legend.background = element_rect(color="black"))+
  facet_wrap(~"PLS-DA")+
  labs(x=paste("Comp 1: ", round(plsda_model$prop_expl_var$X[1]*100, digits=1), "%", sep=""),
       y=paste("Comp 2: ", round(plsda_model$prop_expl_var$X[2]*100, digits=1), "%", sep=""))
plsda_model_plot
# poor discimination

# select top 20 features and add VIMP score
top_features # done

# :: Heatmap --------------------------------------------------------------

# Note: cannot continue this analysis because of vast NA count
# And confounded by right-censoring (LFC at flare at 2M doesn't make sense to include)

# :: Omics Heatmaps // defunct-------------------------------------------------------------

# Made defunct because there are only 4 sig features across omics
# Perhaps the story will be with MBX?

# select significant features and plot a heatmap of Log2FC relative to baseline
# note: several samples are missing for metaproteomics
# which means we will need to use Log2FC relative to earliest and latest sample

# subset to nominally sig taxa
lsarp.asv.data.glom.heatmap = lsarp.asv.data.glom[,colnames(lsarp.asv.data.glom) %in% 
                                                    subset(lsarp.lmer.asv.interactions, padj<0.20 & phase == "treatment")$taxa]
# processing: TSS, pseudo, log2
lsarp.asv.data.glom.heatmap = lsarp.asv.data.glom.heatmap / 50000
lsarp.asv.data.glom.heatmap = log2(lsarp.asv.data.glom.heatmap+min(lsarp.asv.data.glom.heatmap[lsarp.asv.data.glom.heatmap!=0])/2)
# melt
lsarp.asv.data.glom.heatmap = lsarp.asv.data.glom.heatmap %>% reshape2::melt()
colnames(lsarp.asv.data.glom.heatmap) = c("standard.name", "taxa", "value")
# merge
lsarp.asv.data.glom.heatmap = merge(lsarp.asv.data.glom.heatmap,
                                    metadata.lsarp.stool, by="standard.name")
# Log2FC earliest vs latest on RS
lsarp.asv.data.glom.heatmap = lsarp.asv.data.glom.heatmap %>%
  subset(phase == "treatment") %>%
  group_by(HM, taxa) %>%
  filter(lsarp.days == min(lsarp.days) | lsarp.days == max(lsarp.days))%>%
  group_by(HM, taxa) %>%
  mutate(lfc.rs = value - value[lsarp.days == min(lsarp.days)]) %>%
  # keep only 6M sample
  slice_max(order_by=lsarp.days, n=1) %>% as.data.frame()
# cast
lsarp.asv.data.glom.heatmap = reshape2::acast(lsarp.asv.data.glom.heatmap,
                                              HM ~ taxa, value.var="lfc.rs")
pheatmap::pheatmap(lsarp.asv.data.glom.heatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(lsarp.asv.data.glom.heatmap), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(lsarp.asv.data.glom.heatmap)/100, max(lsarp.asv.data.glom.heatmap), length.out=floor(100/2))))

## Repeat for MGX

# subset to nominally sig taxa
lsarp.mgx.heatmap = lsarp.mgx.taxa[,colnames(lsarp.mgx.taxa) %in% 
                                     subset(lsarp.lmer.mgx.interactions, pval<0.05 & phase == "treatment")$taxa]
# processing: TSS, pseudo, log2
lsarp.mgx.heatmap = log2(lsarp.mgx.heatmap+min(lsarp.mgx.heatmap[lsarp.mgx.heatmap!=0])/2)
# melt
lsarp.mgx.heatmap = lsarp.mgx.heatmap %>% reshape2::melt()
colnames(lsarp.mgx.heatmap) = c("standard.name", "taxa", "value")
# merge
lsarp.mgx.heatmap = merge(lsarp.mgx.heatmap,
                          metadata.lsarp.stool, by="standard.name")
# Log2FC earliest vs latest on RS
lsarp.mgx.heatmap = lsarp.mgx.heatmap %>%
  subset(phase == "treatment") %>%
  group_by(HM, taxa) %>%
  filter(lsarp.days == min(lsarp.days) | lsarp.days == max(lsarp.days))%>%
  group_by(HM, taxa) %>%
  mutate(lfc.rs = value - value[lsarp.days == min(lsarp.days)]) %>%
  # keep only 6M sample
  slice_max(order_by=lsarp.days, n=1) %>% as.data.frame()
# cast
lsarp.mgx.heatmap = reshape2::acast(lsarp.mgx.heatmap,
                                    HM ~ taxa, value.var="lfc.rs")
pheatmap::pheatmap(lsarp.mgx.heatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(lsarp.mgx.heatmap), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(lsarp.mgx.heatmap)/100, max(lsarp.mgx.heatmap), length.out=floor(100/2))))

## Repeat for MPX

# subset to nominally sig taxa
lsarp.mpx.heatmap = lsarp.mpx.cazy.mat[,colnames(lsarp.mpx.cazy.mat) %in% 
                                         subset(lsarp.lmer.mpx.cazy.interactions, pval<0.05 & phase == "treatment")$taxa]
# processing: NA-->0, pseudo, log2
lsarp.mpx.heatmap[is.na(lsarp.mpx.heatmap)] = 0
lsarp.mpx.heatmap = log2(lsarp.mpx.heatmap+min(lsarp.mpx.heatmap[lsarp.mpx.heatmap!=0])/2)
# melt
lsarp.mpx.heatmap = as.matrix(lsarp.mpx.heatmap) %>% reshape2::melt()
colnames(lsarp.mpx.heatmap) = c("standard.name", "taxa", "value")
# merge
lsarp.mpx.heatmap = merge(lsarp.mpx.heatmap,
                          metadata.lsarp.stool, by="standard.name")
# Log2FC earliest vs latest on RS
lsarp.mpx.heatmap = lsarp.mpx.heatmap %>%
  subset(phase == "treatment") %>%
  group_by(HM, taxa) %>%
  filter(lsarp.days == min(lsarp.days) | lsarp.days == max(lsarp.days))%>%
  group_by(HM, taxa) %>%
  mutate(lfc.rs = value - value[lsarp.days == min(lsarp.days)]) %>%
  # keep only 6M sample
  slice_max(order_by=lsarp.days, n=1) %>% as.data.frame()
# cast
lsarp.mpx.heatmap = reshape2::acast(lsarp.mpx.heatmap,
                                    HM ~ taxa, value.var="lfc.rs")
pheatmap::pheatmap(lsarp.mpx.heatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(lsarp.mpx.heatmap), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(lsarp.mpx.heatmap)/100, max(lsarp.mpx.heatmap), length.out=floor(100/2))))

# combine datasets and create mapping/annotation file
lsarp.omic.heatmap = merge(lsarp.asv.data.glom.heatmap,
                           lsarp.mgx.heatmap, by="row.names")
rownames(lsarp.omic.heatmap) = lsarp.omic.heatmap$Row.names
lsarp.omic.heatmap$Row.names = NULL
lsarp.omic.heatmap = merge(lsarp.omic.heatmap,
                           lsarp.mpx.heatmap, by="row.names")
rownames(lsarp.omic.heatmap) = lsarp.omic.heatmap$Row.names
lsarp.omic.heatmap$Row.names = NULL
# HM annotations
lsarp.omic.heatmap.hm.map = data.frame(HM = rownames(lsarp.omic.heatmap)) %>%
  merge(metadata.lsarp.stool[,c("HM", "group", "ave.fiber")] %>% distinct(),
        by="HM")
rownames(lsarp.omic.heatmap.hm.map) = lsarp.omic.heatmap.hm.map$HM
lsarp.omic.heatmap.hm.map$HM = NULL
colnames(lsarp.omic.heatmap.hm.map) = c("group","fiber intake")


# Feature annotations
lsarp.omic.heatmap.feature.map = rbind(lsarp.lmer.asv.interactions %>% mutate(data.type = "ASV"),
                                       lsarp.lmer.mgx.interactions %>% mutate(data.type = "MGX"),
                                       lsarp.lmer.mpx.cazy.interactions %>% mutate(data.type = "MPX")) %>% subset(pval < 0.05 & phase == "treatment")
lsarp.omic.heatmap.feature.map = lsarp.omic.heatmap.feature.map[,c("taxa", "estimate", "data.type")]
rownames(lsarp.omic.heatmap.feature.map) = lsarp.omic.heatmap.feature.map$taxa
lsarp.omic.heatmap.feature.map$taxa = NULL
colnames(lsarp.omic.heatmap.feature.map) = c("interaction","data type")

# lastly, scale Log2FC by feature
lsarp.omic.heatmap = reshape2::melt(as.matrix(lsarp.omic.heatmap)) %>%
  group_by(Var2) %>%
  mutate(slfc = scale(value)) %>%
  reshape2::acast(Var2 ~ Var1, value.var="slfc")

# heatmaps
pheatmap::pheatmap(lsarp.omic.heatmap,
                   scale_rows=T,
                   fontsize_row = 7,
                   fontsize_col = 7,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(lsarp.omic.heatmap), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(lsarp.omic.heatmap)/100, max(lsarp.omic.heatmap), length.out=floor(100/2))),
                   annotation_col=lsarp.omic.heatmap.hm.map,
                   annotation_row = lsarp.omic.heatmap.feature.map,
                   annotation_colors = list(`group` = c(Placebo = gg_color_hue(2)[1],
                                                        RS = gg_color_hue(2)[2]),
                                            `interaction` = colorRampPalette(c("blue","white", "red"))(100),
                                            `fiber intake`= colorRampPalette(c("blue","white", "red"))(100),
                                            `data type` = c(ASV = RColorBrewer::brewer.pal(n = 3, name = "Set3")[1],
                                                            MGX = RColorBrewer::brewer.pal(n = 3, name = "Set3")[2],
                                                            MPX = RColorBrewer::brewer.pal(n = 3, name = "Set3")[3])))


# :: PLS-DA // defunct -----------------------------------------------------------------

# Defunct

# Goal: Find features correlated with fecal calprotectin
# that discriminate Responders from Non-Responders

# Step 1: Prepare data
# Assume lfc_matrix is your data (rows = samples, columns = features)
# Example: lfc_matrix <- metadata.lsarp.stool.omics.cor[, numeric_cols]
# labels is a factor vector of "strong" vs "weak"
# Replace with your actual data

# subset to RS group; features sig associated with fecal calprotectin (padj < 0.20)
lfc_matrix <- lsarp.delta(data = lsarp.fecalcal.omics.data[,colnames(lsarp.fecalcal.omics.data)%in% c(subset(metadata.lsarp.stool.omics.cor.bh, padj < 0.20)$feature)],
                          type = "trajectory")
# convert to matrix
lfc_matrix = reshape2::acast(subset(lfc_matrix, phase == "treatment"), HM ~ feature, value.var = "coef")
dim(lfc_matrix) # n = 26 x 33

lfc_matrix = lfc_matrix[rownames(lfc_matrix) %in% subset(lsarp.metadata.responders.asv, Group == "RS")$HM,]
dim(lfc_matrix) # n = 14 x 33

labels = subset(lsarp.metadata.responders.asv, Group == "RS" & HM %in% rownames(lfc_matrix)) %>% dplyr::select(HM, ave.fiber, Group, RS_Name, flare.group) %>% distinct()
rownames(labels) = labels$HM
labels$HM = NULL

# Check for NAs and impute if needed
lfc_matrix[is.na(lfc_matrix)] <- 0  # Simple imputation; consider mixOmics::nipals for better handling

# Step 2: Run PLS-DA
plsda_model <- mixOmics::plsda(X = scale(lfc_matrix), 
                               Y = labels$flare.group, 
                               ncomp = 10)  # Try 5 components

# Tune number of components using cross-validation
set.seed(25)
perf_plsda <- mixOmics::perf(plsda_model, 
                             #validation = "Mfold", 
                             nrepeat = 10,
                             folds = 5, progressBar = TRUE)

ggplot(perf_plsda$error.rate.class$mahalanobis.dist%>%
         data.frame() %>%
         summarize(total.err = colSums(.)) %>%
         mutate(ncomp = seq(1:10)))+
  aes(x=ncomp, y=total.err*100)+
  scale_x_continuous(breaks=seq(1:10))+
  geom_line(color="red", linewidth=1)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Method: Mahalanobis Distance")+
  labs(x="Number of Components",
       y="Total Error (%)")

ncomp_opt = 4

# Refit with optimal components
plsda_model <- mixOmics::plsda(X = scale(lfc_matrix), 
                               Y = labels$flare.group, 
                               ncomp = ncomp_opt)

# Step 4: Extract discriminative features (VIP scores)
vip_scores <- mixOmics::vip(plsda_model)  # VIP for each component
vip_df <- data.frame(
  Feature = colnames(lfc_matrix),
  VIP = vip_scores[, ncomp_opt]  # Use last component
) %>%
  arrange(desc(VIP))

ggplot(vip_df,
       aes(x=reorder(Feature, VIP), y=VIP))+
  coord_flip()+
  geom_hline(yintercept=1, linetype=2, alpha=0.5)+
  geom_point(shape=21, aes(fill=scale(VIP)))+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none")

# Filter top discriminative features (e.g., VIP > 1)
top_features <- vip_df %>%
  arrange(-VIP) %>%
  filter(VIP > 1)# %>%
#head(15)  # Top 20 features
nrow(top_features)

# Step 5: Visualize
# Score plot (sample separation)
mixOmics::plotIndiv(plsda_model, comp = c(1, 2), group = labels$flare.group, legend = TRUE,
                    title = "PLS-DA Score Plot (Strong vs Weak)", ellipse = TRUE)
# plot manually

plsda_model_plot = plsda_model$variates$X[,c(1,2)] %>% data.frame() %>%
  merge(labels, by="row.names") %>%
  ggplot(aes(x=comp1, y=comp2))+
  stat_ellipse(aes(color=flare.group), linewidth=1.5, alpha=0.5)+
  scale_color_discrete(guide = "none")+
  geom_point(aes(fill=flare.group), shape=21, size=3)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position=c(0.8, 0.15),
                        legend.title = element_blank(),
                        legend.background = element_rect(color="black"))+
  facet_wrap(~"PLS-DA")+
  labs(x=paste("Comp 1: ", round(plsda_model$prop_expl_var$X[1]*100, digits=1), "%", sep=""),
       y=paste("Comp 2: ", round(plsda_model$prop_expl_var$X[2]*100, digits=1), "%", sep=""))
plsda_model_plot

# select top features and add VIMP score
top_features # done

# :: Heatmap --------------------------------------------------------------
lsarp.omics.fcal.lfc.map = subset(metadata.lsarp.stool.omics.cor.bh, padj < 0.20)
rownames(lsarp.omics.fcal.lfc.map) = lsarp.omics.fcal.lfc.map$feature
lsarp.omics.fcal.lfc.map = merge(lsarp.omics.fcal.lfc.map,
                                 vip_df, by="row.names")
rownames(lsarp.omics.fcal.lfc.map) = lsarp.omics.fcal.lfc.map$Row.names
lsarp.omics.fcal.lfc.map$Row.names = NULL
lsarp.omics.fcal.lfc.map$Feature = NULL
# add datatype
lsarp.omics.fcal.lfc.map = lsarp.omics.fcal.lfc.map %>%
  mutate(data.type = ifelse(feature %in% colnames(lsarp.asv.data.glom), "ASV",
                            ifelse(feature %in% colnames(lsarp.mgx.taxa), "Species", 
                                   ifelse(feature %in% colnames(lsarp.cd.mpx.kegg.mat), "Pathway", 
                                          ifelse(feature %in% colnames(lsarp.cd.mpx.cog.mat), "COG", 
                                                 ifelse(feature %in% colnames(lsarp.cd.mpx.cazy.mat), "CAZy", "Fecal calprotectin"))))))
# remove HM with <3 samples
lfc_matrix_hm_keep = data.frame(table(subset(lsarp.metadata.responders.asv, HM %in% rownames(lfc_matrix) & phase == "treatment")$HM)) %>%
  subset(Freq > 2) %>% dplyr::select(Var1)

pheatmap::pheatmap(t(lfc_matrix[rownames(lfc_matrix) %in% lfc_matrix_hm_keep$Var1,]),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = labels[rownames(labels) %in% lfc_matrix_hm_keep$Var1,] %>% dplyr::select(-Group),
                   annotation_colors = list(flare.group = c(`Relapse` = gg_color_hue(2)[1],
                                                            `Remit` = gg_color_hue(2)[2]),
                                            `value` = colorRampPalette(c("blue","white", "red"))(100),
                                            `VIP` = colorRampPalette(c("white", "red"))(100),
                                            `data.type` = omics.colors,
                                            `RS_Name` = rs.colors),
                   annotation_row = lsarp.omics.fcal.lfc.map[,c("VIP","value", "data.type")],
                   breaks=c(seq(min(na.omit(lfc_matrix)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(lfc_matrix))/100, max(na.omit(lfc_matrix)), length.out=floor(100/2))))
# correlation clustering
# scaled Log2FC



# :: Random Survival Forest // defunct ------------------------------------

# Survival data:
redcap.lsarp.wpcdai$time = ifelse(is.na(redcap.lsarp.wpcdai$time.to.event),
                                  max(na.omit(redcap.lsarp.wpcdai$time.to.event)),
                                  (redcap.lsarp.wpcdai$time.to.event))
redcap.lsarp.wpcdai$event = ifelse(is.na(redcap.lsarp.wpcdai$time.to.event), 0, 1)

# note: cannot predict responder vs nonresponder among RS group

# loop through datasets (add clinical and survival data in loop)
t1 = Sys.time()
lsarp.rsf.validation.output = do.call(rbind, lapply(1:15, function(iter){
  do.call(rbind, lapply(data.sets, function(data.set){
    print(paste(iter, data.set))
    # select dataset
    if(data.set == "ASV"){
      lsarp.data = lsarp.asv.baseline
    }
    if(data.set == "MGX"){
      lsarp.data = lsarp.mgx.baseline
    }
    if(data.set == "KEGG"){
      lsarp.data = lsarp.kegg.baseline
    }
    if(data.set == "COG"){
      lsarp.data = lsarp.cog.baseline
    }
    if(data.set == "CAZy"){
      lsarp.data = lsarp.cazy.baseline
    }
    if(data.set == "ALL"){
      lsarp.data = lsarp.omics.baseline
    }
    # add RS + survival data
    lsarp.data = merge(lsarp.data,
                       metadata.lsarp.stool.asv[,c("Group", "RS_Name", "standard.name")],
                       by="standard.name")
    lsarp.data$HM = substr(lsarp.data$standard.name, 1, 6)
    lsarp.data$standard.name = NULL
    lsarp.data = merge(redcap.lsarp.wpcdai[,c("HM", "time", "event")],
                       lsarp.data,
                       by="HM")
    # subset to RS group
    lsarp.data = subset(lsarp.data, Group == "RS")
    # ready to run
    lsarp.data = lsarp.data %>%
      mutate(HM = NULL,
             Group = as.factor(Group),
             RS_Name = as.factor(RS_Name),
             time = as.numeric(time),
             event = as.numeric(event))
    
    set.seed(iter)
    rfsrc_output = randomForestSRC::rfsrc(Surv(time, event) ~ .,
                                          lsarp.data, ntree=500,
                                          nodesize = 3,
                                          importance=T)
    data.frame(imp = rfsrc_output$importance) %>% arrange(-imp)
    error.rate = rfsrc_output$err.rate[length(rfsrc_output$err.rate)]
    # output data
    data.frame(err = error.rate,
               iter = iter,
               data.set = data.set)
  }))
}))
t2 = Sys.time()
t2-t1 # 1 min

# tidy up for plot

ggplot(lsarp.rsf.validation.output%>%
         mutate(data.set = ifelse(data.set=="KEGG", "Pathway", data.set))%>%
         mutate(data.set = factor(data.set, levels=c("ASV", "MGX", "Pathway", "COG", "CAZy", "ALL"))),
       aes(x=data.set, y= err))+
  geom_boxplot()+
  geom_point(shape=21, aes(fill=data.set))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(n = 5, name = "Set3"), "black"))+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Random Survival Forest")+
  labs(x="Data type", y="OOB Prediction Error")

# note: try comparing with or without subsetting to RS group


# :: RF features ----------------------------------------------------------



# :: RF Importances -------------------------------------------------------

lsarp.rf.models.importances = rf.function(data.types = c("Pathway"),
                                          iters = 15,
                                          reduce_features = T,
                                          p.threshold = 0.10,
                                          output = "importances")

lsarp.rf.models.importances.df = lsarp.rf.models.importances %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(mean.imp = mean(na.omit(imp))) %>%
  mutate(imp.low = mean(imp)- (sd(imp)/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(imp)+ (sd(imp)/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high) %>% distinct() %>%
  arrange(-mean.imp) %>% as.data.frame()

# log2fc
lsarp.rf.models.importances.df = lsarp.rf.models.importances.df %>% 
  subset(!is.na(mean.imp))%>%
  arrange(mean.imp) %>% slice_max(mean.imp, n=20) 
# calculate wilcox
lsarp.rf.models.importances.wilcox = #do.call(rbind, lapply(lsarp.rf.models.importances.df$feature, function(x){
  do.call(rbind, lapply(colnames(lsarp.omics.baseline[,!colnames(lsarp.omics.baseline) %in% c("RS_Name", "flare.group", "standard.name")]), function(x){
    
    print(x)
    #x = lsarp.rf.models.importances.df$feature[1]
    # add meta
    new.data = lsarp.omics.baseline
    new.data = merge(new.data,
                     metadata.lsarp.stool.asv[,c("Group", "RS_Name", "standard.name","ave.fiber", "richness", "shannon", "but.i", "but.ii", "fd","load.asv","load.mgx","starch","mucin","starch.mucin")],
                     by="standard.name")
    new.data$HM = substr(new.data$standard.name, 1, 6)
    new.data$standard.name = NULL
    # subset to RS group
    new.data = subset(new.data, Group == "RS")
    # add flare.group
    new.data = merge(lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                     new.data,
                     by="HM")
    
    data.subset = new.data[,colnames(new.data) %in% c(x, "flare.group")]
    # run wilcox
    wilcox.p = wilcox.test(subset(data.subset, flare.group == "Relapse")[,2],
                           subset(data.subset, flare.group == "Remit")[,2])
    # for omics, log2-transform
    if(x %in% colnames(lsarp.omics.baseline)){
      omics.values = data.subset[,2]
      data.subset[,2] = log2(omics.values + min(omics.values[omics.values!=0])/2)
    }
    coef = mean((subset(data.subset, flare.group == "Remit")[,2])) - mean((subset(data.subset, flare.group == "Relapse")[,2]))
    
    data.frame(feature = x,
               pval = wilcox.p$p.value,
               coef = coef)
  }))
lsarp.rf.models.importances.df = merge(lsarp.rf.models.importances.wilcox,
                                       lsarp.rf.models.importances.df, by="feature")
lsarp.rf.models.importances.df$sig = ifelse(lsarp.rf.models.importances.df$pval < 0.05, "*", "")
lsarp.rf.models.importances.df = lsarp.rf.models.importances.df %>% arrange(mean.imp)

# clean feature names (remove long COG names)
lsarp.rf.models.importances.df$clean.feature = gsub(" \\(.*", "", lsarp.rf.models.importances.df$feature)

lsarp.rf.loocv.imp.plot = ggplot(lsarp.rf.models.importances.df,
                                 aes(x=mean.imp, y=reorder(clean.feature, mean.imp)))+
  geom_segment(aes(x=imp.low, xend=imp.high, y=clean.feature, yend=clean.feature), color="black")+
  geom_point(shape=21, aes(fill=(coef)), size=3.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("white", "black"))+
  theme_classic()+theme(legend.position="right",
                        axis.title.y = element_blank(),
                        #axis.text.y = element_blank(),
                        # axis.ticks.y=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  guides(color=FALSE)+
  labs(x="Mean Decrease in Accuracy", fill="Log2FC")+
  facet_wrap(~"Feature Importance")
lsarp.rf.loocv.imp.plot


# :: RF ttest plots -------------------------------------------------------

# Just do all features

# loop through top features and plot; extract microbiome values
lsarp.rf.models.importances.ttest = #do.call(rbind, lapply(slice_max(lsarp.rf.models.importances.df, mean.imp, n=10)$feature, function(x){
  do.call(rbind, lapply(colnames(lsarp.omics.baseline[,!colnames(lsarp.omics.baseline) %in% c("RS_Name", "flare.group", "standard.name")]), function(x){
    print(x)
    # add meta
    new.data = lsarp.omics.baseline
    new.data = merge(new.data,
                     metadata.lsarp.stool.asv[,c("Group", "standard.name","ave.fiber", "richness", "shannon", "but.i", "but.ii", "fd","load.asv","load.mgx","starch","mucin","starch.mucin")],
                     by="standard.name")
    new.data = subset(new.data, Group == "RS")
    new.data$HM = substr(new.data$standard.name, 1, 6)
    new.data$standard.name = NULL
    new.data = merge(lsarp.metadata.responders[,c("HM", "flare.group")] %>% distinct(),
                     new.data,
                     by="HM")
    
    # extract values 
    data.subset = new.data[,c(x, "flare.group")]
    data.subset$feature = colnames(data.subset)[1]
    data.subset$sample = rownames(data.subset)
    colnames(data.subset)[1] = "value"
    data.subset
  }))
# add stats
lsarp.rf.models.importances.ttest = merge(lsarp.rf.models.importances.ttest,
                                          lsarp.rf.models.importances.df, by="feature")
lsarp.rf.models.importances.ttest.p = lsarp.rf.models.importances.ttest[,c("feature", "pval")] %>% distinct()

# now clip names to 30 char
lsarp.rf.models.importances.ttest$short.feature = ifelse(nchar(lsarp.rf.models.importances.ttest$feature)>30, paste(substr(lsarp.rf.models.importances.ttest$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.ttest$feature)
lsarp.rf.models.importances.ttest.p$short.feature = ifelse(nchar(lsarp.rf.models.importances.ttest.p$feature)>30, paste(substr(lsarp.rf.models.importances.ttest.p$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.ttest.p$feature)
lsarp.rf.models.importances.df$short.feature = ifelse(nchar(lsarp.rf.models.importances.df$feature)>30, paste(substr(lsarp.rf.models.importances.df$feature, start = 1, stop = 30), "...", sep=""), lsarp.rf.models.importances.df$feature)

lsarp.rf.models.importances.ttest$short.feature = gsub("shannon", "Shannon Diversity", lsarp.rf.models.importances.ttest$short.feature)
lsarp.rf.models.importances.ttest$short.feature = gsub("but.ii", "Kircher Butyrogens", lsarp.rf.models.importances.ttest$short.feature)
lsarp.rf.models.importances.ttest$short.feature = gsub("but.i", "Butyrogens", lsarp.rf.models.importances.ttest$short.feature)

lsarp.rf.models.importances.ttest.p$short.feature = gsub("shannon", "Shannon Diversity", lsarp.rf.models.importances.ttest.p$short.feature)
lsarp.rf.models.importances.ttest.p$short.feature = gsub("but.ii", "Kircher Butyrogens", lsarp.rf.models.importances.ttest.p$short.feature)
lsarp.rf.models.importances.ttest.p$short.feature = gsub("but.i", "Butyrogens", lsarp.rf.models.importances.ttest.p$short.feature)



lsarp.rf.models.importances.ttest.plot = ggplot(lsarp.rf.models.importances.ttest %>%
                                                  # clean names
                                                  mutate(Group = flare.group) %>%
                                                  # add indicator of not present
                                                  group_by(feature) %>%
                                                  mutate(pseudo = ifelse(value == min(value), "pseudo", "real")) %>%
                                                  # reorder taxa based on imp
                                                  mutate(short.feature = factor(short.feature, levels=arrange(distinct(lsarp.rf.models.importances.df[,c("short.feature", "mean.imp")]),-mean.imp)$short.feature)),
                                                aes(x=Group, y=value))+
  geom_boxplot(width=0.3, outlier.shape=NA)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  #scale_y_log10()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=Group, alpha=pseudo), size=2)+
  scale_alpha_manual(values=c(0.2, 1))+
  theme_classic()+theme(legend.position="none",
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=8),
                        strip.background = element_rect(color="black"))+
  geom_text(data=lsarp.rf.models.importances.ttest.p%>%
              # clean names
              mutate(short.feature = factor(short.feature, levels=arrange(distinct(lsarp.rf.models.importances.df[,c("short.feature", "mean.imp")]),-mean.imp)$short.feature)),
            x=1.5, y=Inf, vjust=1.2, 
            aes(label = paste("p =", round(pval, digits=3))),
            size=2.5)+
  labs(x="", y="Feature Abundance")+
  facet_wrap(~short.feature, ncol=2, scales="free")
lsarp.rf.models.importances.ttest.plot

# Sulfur Metabolism!

# :: RF Interactions ------------------------------------------------------

# check interactions of most important features
# (otherwise, intractable to compute, even on compute canada)

lsarp.loocv.models.interactions = rf.function(data.types = c("Multi-Omic"),
                                              iters = 15,
                                              output = "interactions")

# analyze
lsarp.loocv.models.interactions.df = 
  # watch this
  rbind(lsarp.loocv.models.interactions %>% mutate(var3 = var1,
                                                   var1 = var2,
                                                   var2 = var3) %>% dplyr::select(-var3),
        lsarp.loocv.models.interactions) %>% as.data.frame() %>%
  # take mean
  group_by(features) %>%
  mutate(value = mean(Difference))%>%
  dplyr::select(`var1`, `var2`, value)

lsarp.loocv.models.interactions.df = reshape2::acast(lsarp.loocv.models.interactions.df, var1~var2, value.var="value", fun.aggregate=mean)%>%
  as.matrix() %>% reshape2::melt() %>%
  mutate(feature.1 = as.character(Var1),
         feature.2 = as.character(Var2))

lsarp.loocv.models.interactions.df$value %>% range()

# clip
lsarp.loocv.models.interactions.df$short.feature.1 = ifelse(nchar(lsarp.loocv.models.interactions.df$feature.1)>30, paste(substr(lsarp.loocv.models.interactions.df$feature.1, start = 1, stop = 30), "...", sep=""), lsarp.loocv.models.interactions.df$feature.1)
lsarp.loocv.models.interactions.df$short.feature.2 = ifelse(nchar(lsarp.loocv.models.interactions.df$feature.2)>30, paste(substr(lsarp.loocv.models.interactions.df$feature.2, start = 1, stop = 30), "...", sep=""), lsarp.loocv.models.interactions.df$feature.2)


# plot
lsarp.rf.loocv.interactions.plot = lsarp.loocv.models.interactions.df %>%
  #mutate(ratio = log10(Paired - Additive)) %>%
  group_by(Var2) %>%
  mutate(sum.1 = max(na.omit(value)))%>%
  ungroup()%>%
  subset(Var1 %in% lsarp.rf.models.importances.df$feature) %>% distinct()%>%
  #slice_max(order_by=na.omit(value), n=10) %>%
  subset(value > 0) %>%
  dplyr::select(feature.1, short.feature.2, sum.1, value) %>% 
  rbind(data.frame(feature.1 = lsarp.rf.models.importances.df$feature,
                   short.feature.2 = "Marginal Importance",
                   sum.1 = 11,
                   value = lsarp.rf.models.importances.df$mean.imp)) %>%
  mutate(clean.feature.1 = factor(feature.1, levels=(lsarp.rf.models.importances.df$feature))) %>%
  mutate(short.feature.2 = factor(short.feature.2, levels=rev(c(lsarp.rf.models.importances.df$short.feature,"Marginal Importance")))) %>%
  ggplot(
    aes(y=clean.feature.1, x=short.feature.2))+
  geom_tile(aes(fill=scale(value)), color="white", size=2)+
  geom_vline(xintercept=1.5, color="black")+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interactions")+
  theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1),
                        strip.text=element_text(size=10),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
  labs(fill="Scaled\nImportance")
lsarp.rf.loocv.interactions.plot
# fixed!
# Difference = Conditional Importance - sum(Marginal Importance)
# Looking at abs(Difference) would be anti-conditional interactions


