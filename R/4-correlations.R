library(tidyverse)
library(car)

gaming_df <- read.csv("data/gaming.csv") 
colnames(gaming_df)[1] <- "id"
gaming_df$id <- factor(gaming_df$id)
gaming_df <- gaming_df %>% filter(!id %in% "71")

go_ssrt_df <- read.csv("img/FDSSRTs.csv")
colnames(go_ssrt_df) <- c("id", "goSSRT")
go_ssrt_df$id <- factor(go_ssrt_df$id)

ARI_ssrt_df <- read.csv("img/ari_SSRTs_integrate.csv")
colnames(ARI_ssrt_df) <- c("ARISSRT", "id")
ARI_ssrt_df$id <- factor(ARI_ssrt_df$id)

bad <- c("8", "23", "46", "54", "4", "48", "21")

  
UPPSP <- read.csv("data/UPPSP.csv") %>% filter(!id %in% c("4_old",
                                                         "15_old",
                                                         "42a",
                                                         "71"))
gaming_go <- full_join(gaming_df, go_ssrt_df, by="id")
UPPS_ARI <- full_join(UPPSP, ARI_ssrt_df, by="id")
  
full_cor_df <- full_join(gaming_go, UPPS_ARI, by ="id") %>% filter(!id %in% bad)

full_cor_df <- full_cor_df %>% filter(unweighted<500)

cor.test(full_cor_df$unweighted, full_cor_df$goSSRT)
cor.test(full_cor_df$inhibitory, full_cor_df$goSSRT)

cor.test(full_cor_df$unweighted, full_cor_df$ARISSRT)
cor.test(full_cor_df$inhibitory, full_cor_df$ARISSRT)


fit <- lm(scale(ARISSRT) ~ scale(n_urgency) + 
            scale(lack_premed) + 
            scale(lack_perserverance) + 
            scale(sensation_seeking)
          +scale(p_urgency), data=full_cor_df)
summary(fit)

fit <- lm(scale(goSSRT) ~ scale(n_urgency) + 
            scale(lack_premed) + 
            scale(lack_perserverance) + 
            scale(sensation_seeking)
          +scale(p_urgency), data=full_cor_df)
summary(fit)

cor.test(full_cor_df$goSSRT, full_cor_df$ARISSRT)

plot(full_cor_df$goSSRT, full_cor_df$ARISSRT)


