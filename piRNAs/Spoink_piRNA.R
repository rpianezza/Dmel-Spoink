library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

#gdl supplement
df <- read.table("/Users/ascarpa/Spoink/piRNA_Spoink/gdl/sam/bam/Spoink-graph-on-TEforR")
names(df) <- c("Run", "te","pos","pirna")
df_Spoink <- subset(df, te=="Spoink")
df_Spoink <- subset(df_Spoink, Run!="SRR3123288_Spoink")

ylow=min(df_Spoink$pirna)
yhigh=max(df_Spoink$pirna)

g_1 <- ggplot(df_Spoink,aes(x=pos,y=pirna))+
  geom_segment(aes(xend=pos),yend=0)+
  ylim(ylow,yhigh)+
  ylab("piRNA abundance [pmp]")+
  xlab("position of piRNA (5' end)")+
  facet_wrap("Run", labeller = labeller(Run = c("SRR5638228_Spoink" = "I06",
                                                "SRR3123284_Spoink" = "B10",
                                                "SRR5638231_Spoink" = "B12",
                                                "SRR5638224_Spoink" = "T05",
                                                "SRR5638225_Spoink" = "T07",
                                                "SRR5638226_Spoink" = "N10",
                                                "SRR5638227_Spoink" = "N16",
                                                "SRR3123310_Spoink" = "ZW155",
                                                "SRR3123314_Spoink" = "ZW184")))+
  theme(text = element_text(size=24))

plot(g_1)


#Ping-pong signature
df2 <- read.table("/Users/ascarpa/Spoink/piRNA_Spoink/gdl/sam/bam/pingpongforR")
names(df2) <- c("Run", "te", "sense", "pos", "frequency")
df2 <- subset(df2, te=="Spoink")
df2 <- subset(df2, Run!="SRR3123288_Spoink")

df2s <- subset(df2, sense == "s")
df2as <- subset(df2, sense == "as")


g_2 <- ggplot(df2s,aes(x=pos,y=frequency))+
  geom_col()+
  ggtitle("Sense piRNAs")+
  ylab("ping-pong signature")+
  xlab("overlap")+
  facet_wrap("Run", labeller = labeller(Run = c("SRR5638228_Spoink" = "I06",
                                                "SRR3123284_Spoink" = "B10",
                                                "SRR5638231_Spoink" = "B12",
                                                "SRR5638224_Spoink" = "T05",
                                                "SRR5638225_Spoink" = "T07",
                                                "SRR5638226_Spoink" = "N10",
                                                "SRR5638227_Spoink" = "N16",
                                                "SRR3123310_Spoink" = "ZW155",
                                                "SRR3123314_Spoink" = "ZW184")))+
  theme(text = element_text(size=24))

plot(g_2)


g_2as <- ggplot(df2as,aes(x=pos,y=frequency))+
  geom_col()+
  ggtitle("Antisense piRNAs")+
  ylab("ping-pong signature")+
  xlab("overlap")+
  facet_wrap("Run", labeller = labeller(Run = c("SRR5638228_Spoink" = "I06",
                                                "SRR3123284_Spoink" = "B10",
                                                "SRR5638231_Spoink" = "B12",
                                                "SRR5638224_Spoink" = "T05",
                                                "SRR5638225_Spoink" = "T07",
                                                "SRR5638226_Spoink" = "N10",
                                                "SRR5638227_Spoink" = "N16",
                                                "SRR3123310_Spoink" = "ZW155",
                                                "SRR3123314_Spoink" = "ZW184")))+
  theme(text = element_text(size=24))

plot(g_2as)


#Iso1
df_Iso1 <- read.table("/Users/ascarpa/Desktop/RNAovaries/noadapt/trim/bam/SRR3123324-graph-on-TE.forR")
names(df_Iso1) <- c("Run", "te","pos","pirna")
df_Spoink_Iso1 <- subset(df_Iso1, te=="Spoink" | te=="1360" | te == "Beagle" )

ylow=min(df_Spoink_Iso1$pirna)
yhigh=max(df_Spoink_Iso1$pirna)

g_1_Iso1 <- ggplot(df_Spoink_Iso1,aes(x=pos,y=pirna))+
  geom_segment(aes(xend=pos),yend=0)+
  ylim(ylow,yhigh)+
  ylab("piRNA abundance [pmp]")+
  xlab("position of piRNA (5' end)")+
  facet_wrap("te")

plot(g_1_Iso1)


#Ping-pong signature
df2_Iso1 <- read.table("/Users/ascarpa/Spoink/piRNA_Spoink/Iso-1/RNAovaries/noadapt/trim/bam/SRR3123324.pps")
names(df2_Iso1) <- c("Run", "te", "sense", "pos", "frequency")
df2_Iso1 <- subset(df2_Iso1, te=="Spoink" | te=="1360" | te == "Beagle" )

df2s_Iso1 <- subset(df2_Iso1, sense == "s")
df2as_Iso1 <- subset(df2_Iso1, sense == "as")


g_2_Iso1 <- ggplot(df2s_Iso1,aes(x=pos,y=frequency))+
  geom_col()+
  ggtitle("Sense piRNAs")+
  ylab("ping-pong signature")+
  xlab("overlap")+
  facet_wrap("te")

plot(g_2_Iso1)



#SRR3123288_Spoink = I17
#SRR5638228_Spoink = I06
#SRR12831808 = Lausanne-S
#SRR3123324 = Iso-1,

#Combined
dfc <- read.table("/Users/ascarpa/Spoink/piRNA_Spoink/combined/SpoinkTEforR")
names(dfc) <- c("Run", "te","pos","pirna")
df_Spoinkc <- subset(dfc, te=="Spoink" | te=="Beagle")
df_Spoinkc <- subset(df_Spoinkc, Run=="SRR12831808" | Run=="SRR3123288_Spoink")
df_Spoinkc <- df_Spoinkc %>%
  mutate(te = ifelse(te == "Beagle", "HMS Beagle", te))


unique(df_Spoinkc$te)
unique(df_Spoinkc$Run)

ylow=min(df_Spoinkc$pirna)
yhigh=max(df_Spoinkc$pirna)

g_1c <- ggplot(df_Spoinkc,aes(x=pos,y=pirna))+
  geom_segment(aes(xend=pos),yend=0)+
  ylim(ylow,yhigh)+
  ylab("piRNA abundance [pmp]")+
  xlab("position")+
  facet_grid( Run  ~ te, scales="free_x", labeller = labeller(Run = 
                                                 c("SRR12831808" = "Lausanne-S",
                                                   "SRR3123288_Spoink" = "I17")))+
  theme(text = element_text(size=24))

plot(g_1c)

#Ping-pong signature combined
df2_c <- read.table("/Users/ascarpa/Spoink/piRNA_Spoink/combined/PingPongforR")
names(df2_c) <- c("Run", "te", "sense", "pos", "frequency")
df2_combined <- subset(df2_c, te=="Spoink" | te == "Beagle" )
df2_combined <- subset(df2_combined, Run=="SRR3123288_Spoink")
df2_combined <- df2_combined %>%
  mutate(te = ifelse(te == "Beagle", "HMS Beagle", te))


df2s_c <- subset(df2_combined, sense == "s")
df2as_c <- subset(df2_combined, sense == "as")


g_2_c <- ggplot(df2s_c,aes(x=pos,y=frequency))+
  geom_col()+
  ylab("ping-pong signature")+
  xlab("overlap")+
  facet_grid( Run  ~ te, scales="free", labeller = labeller(Run = 
                                                              c("SRR3123288_Spoink" = "I17")))+ 
  theme(text = element_text(size=24))

plot(g_2_c)


ggarrange(g_1c+ theme(strip.text = element_text(size = 36)), g_2_c+ theme(strip.text = element_text(size = 36)),
          ncol = 1, nrow = 2, align = ("v"),
          labels = c("A", "B"), heights = c(2,2), widths = c(2,2)
)
