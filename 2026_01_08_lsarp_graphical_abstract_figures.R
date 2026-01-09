### 2025_11_20  Graphical abstract figures

# butyrogen response


# plot
lsarp.butyrogen.ga.plot = ggplot(lsarp.metadata.responders.asv %>%
                                          subset(!is.na(but.i))%>%
                                          subset(flare !="flare"& !HM %in% excluded.by.dave),
                                        aes(x=lsarp.days, 
                                            y=but.i*100))+
  #scale_y_log10()+
  #annotate("rect", xmin=0, xmax=max(subset(metadata.lsarp.stool.asv, phase=="treatment")$lsarp.days),
  #         ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  #geom_line(aes(group=HM, color=flare.group), linetype=2, alpha=0.5, linewidth=0.3)+
  #geom_point(aes(fill=flare.group, shape=baseline), size=3)+
  scale_shape_manual(values=c(23,21))+
  geom_smooth(method="lm", se=T, aes(fill=flare.group), color="white")+
  geom_smooth(method="lm", se=F, aes(fill=flare.group, color=flare.group))+
  #scale_fill_manual(values=c("grey",labelcolors$cols[c(1,1,4,5,5,8,9)]))+
  labs(x="Days since starting Product", y="Fecal Butyrogens %")+
  facet_wrap(~group)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        #strip.background = element_rect(color="black"),
                        strip.background = element_blank())
lsarp.butyrogen.ga.plot



lsarp.roc.ga.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  #scale_fill_manual(values = omics.colors)+
  # Mean ROC curve
  geom_path(data = lsarp.rf.selected.output.roc_summary, 
            aes(x = fpr, y = mean_sens), 
            color="black", linewidth=1.5)+
  
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("AUC: ", round(lsarp.rf.selected.output.stats$mean.auc, digits=2)),
           size=5)+
  theme_classic()+theme(strip.text=element_text(size=12),
                        strip.background = element_blank(),
                        legend.position="none")+
  facet_wrap(~"Multi-Omic Random Forest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "False Positive Rate", y = "True Positive Rate") 
lsarp.roc.ga.plot


(lsarp.butyrogen.ga.plot+
  lsarp.roc.ga.plot+
  patchwork::plot_layout(widths=c(2,1)))%>%
  ggsave(filename="./lsarp_plots/2026_01_08_lsarp_graphical_abstract.pdf",
         width=9, height=3,device = cairo_pdf)

## RS selection

(data.frame(RS_Name = factor(rs.names,levels=rs.names),
           RS_type = factor(c("Potato", "Potato", "Potato", "Green banana", "Corn", "Corn", "Cross-linked", "Cross-linked", "Cross-linked potato"),
                            levels=c("Potato", "Green banana", "Corn", "Cross-linked", "Cross-linked potato")),
           but = c(0.55, 0.45, 0.5, 0.8, 0.1, 0.13, 0.22, 0.25, 0.30)) %>%
  ggplot(aes(x=RS_Name, y=but*100))+
  geom_point(shape=21, aes(fill=RS_type), size=5)+
  scale_fill_manual(values=labelcolors$cols[c(1,4,5,7,9)])+
    scale_y_continuous(limits=c(0,100))+
  theme_classic(base_size=14)+
  theme(legend.position="none",
        axis.text.x = element_blank())+
    guides(fill=guide_legend(nrow=5,byrow=TRUE))+
  labs(y = "Butyrogens %", x=NULL, fill=NULL))%>%
  ggsave(filename="./lsarp_plots/lsarp_graphical_abstract_rapidaim.pdf",
         width=4, height=2,device = cairo_pdf)

  
