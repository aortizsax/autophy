##########
#ploting Random Forest from python
#NARCH 
gini_BSH <- read.csv('../data/RF/2021-03-01_RF_BSH.csv')
g_narch_b <- ggplot(data = gini_BSH, aes(x= reorder(X, Gini.importance),y=Gini.importance))+
  geom_bar(stat = "identity",fill="red",width = 0.5)+ coord_flip()+
  ylim(0,0.22)+ylab('Random Forest Gini Importance')+xlab('Alignment position')+theme_classic()+
  ggtitle('Alignment position by subfamily Gini importance ',)+ 
  theme(plot.title = element_text(size = 35),
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 35))
g_narch_b


gini_NARCH_bact <- read.csv('../data/RF/RF_plot_bact_NARCH_hiLoAvgPD.csv')
gini_NARCH_bact_15 <- read.csv('../data/RF/2020-08-27_15_RF_plot_bact_NARCH_hiLoAvgPD.csv')
gini_NARCH_bact_cyto <- read.csv('../data/RF/RF_plot_bact_cyto_NARCH_hiLoAvgPD.csv')
gini_NARCH_bact_cyto_15 <- read.csv('../data/RF/2020-08-27_15_RF_plot_bact_cyto_NARCH_hiLoAvgPD.csv')

#CLOROX
gini_clorox_bact_15 <- read.csv('../mSysdata/RF/2020-08-28_15_RF_plot_bact_Clorox_AC.csv')
gini_clorox_bact_metab_15 <- read.csv('../mSysdata/RF/2020-08-28_15_RF_plot_bact_metab_Clorox_AC.csv')

g_narch_b <- ggplot(data = gini_NARCH_bact_15, aes(x= reorder(Feature, Gini.importance),y=Gini.importance))+
  geom_bar(stat = "identity",fill="steelblue",width = 0.5)+ coord_flip()+
  ylim(0,0.2)+xlab('')+ylab('')+theme_classic()+
  ggtitle('16S Low/High Pocket Depth NARCH',)+ 
  theme(plot.title = element_text(size = 12))

g_narch_b_c <- ggplot(data = gini_NARCH_bact_cyto_15, aes(x= reorder(Feature, Gini.importance),y=Gini.importance))+
  geom_bar(stat = "identity",fill="steelblue",width = 0.5)+
  coord_flip()+
  ylim(0,0.2)+xlab('')+
  ylab('')+theme_classic()+
  ggtitle('16S+Cytokine Low/High Pocket Depth NARCH')+ 
  theme(plot.title = element_text(size = 12))

g_clorox_b <- ggplot(data = gini_clorox_bact_15, aes(x= reorder(Feature, Gini.importance),y=Gini.importance))+
  geom_bar(stat = "identity",fill="steelblue",width = 0.5)+ 
  coord_flip()+
  ylim(0,0.2)+xlab('')+
  ylab('')+theme_classic()+
  ggtitle('16S A/C Disease Class Clorox')+ 
  theme(plot.title = element_text(size = 12))

g_clorox_b_m <- ggplot(data = gini_clorox_bact_metab_15, aes(x= reorder(Feature, Gini.importance),y=Gini.importance))+
  geom_bar(stat = "identity",fill="steelblue",width = 0.5)+ coord_flip()+
  ylim(0,0.2)+xlab('')+ylab('')+theme_classic()+
  ggtitle('16S+Metabolites A/C Disease Class Clorox')+ 
  theme(plot.title = element_text(size = 12))

RF_figure<- plot_grid(g_narch_b, g_narch_b_c,g_clorox_b, g_clorox_b_m,
                      labels = paste(LETTERS[1:4],')',sep = ''))

RF_figure <- ggdraw(add_sub(add_sub(RF_figure, "Gini Importance",size=12),'Figure #',x=0,hjust = -0.1,size = 12))
# ggraw(add_sub(,'Figure #'))

#save figures
ggsave('../data/RF/2020-10-09_RF_gini_improtance_FIGURE_#14.pdf',width=10.5, height=7, plot = RF_figure)
ggsave('../data/RF/2020-10-09_bact_narch gini_hilo.pdf',width = 5,height = 5, plot = g_narch_b)
ggsave('../data/RF/2020-10-09_bact_cyto_narch gini_hilo.pdf',width = 5,height = 5, plot = g_narch_b_c)
ggsave('../mSysdata/RF/2020-10-09_bact_clorox_gini_AC.pdf',width = 5,height = 5, plot = g_clorox_b)
ggsave('../mSysdata/RF/2020-10-09_bact_metab_clorox_gini_AC.pdf',width = 5,height = 5, plot = g_clorox_b_m)



