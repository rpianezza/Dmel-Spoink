library(tidyverse)
theme_set(theme_bw())

df_dNdS <- read.csv("/Users/ascarpa/Desktop/dn_ds.txt", header = TRUE, sep = "\t")

quantile_2.5 <- quantile(df_dNdS$dS, probs = 0.025)

ggplot(df_dNdS, aes(x = dS)) +
  geom_histogram(binwidth = 0.1, fill = "orange", color = "black") +
  geom_vline(xintercept = quantile_2.5, color = "red", linetype = "dashed") +
  labs(title="dS on 140 genes", x = "dS", y = "Density")

df_dNdS$Type_genetic<-"gene"
df_dNdS[nrow(df_dNdS) + 1,] <- c(0.0680, 0.0009, 0.0135, "Spoink")
df_dNdS$dS <- as.numeric(as.character(df_dNdS$dS))



ggplot(df_dNdS, aes(x = dS, fill = Type_genetic)) +
  geom_histogram(binwidth = 0.1, color = "black") +
  geom_vline(xintercept = quantile_2.5, color = "red", linetype = "dashed") +
  labs(title = "dS on 140 genes", x = "dS", y = "Density")
