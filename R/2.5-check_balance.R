#Checks how many we have for 
#each counterbalance 
#accounting for all the mess
#discussed in the readme
load("data/final_dfs.RData")

idnums <- as.numeric(unique(ARdats$id))
pre50 <- idnums[idnums<52]

pre50 <- pre50[!pre50 %in% c(4,48,46,23)]

post50 <- idnums[idnums>51]
post50 <- c(post50, 54, 57)

#Almost a full counterbalance
table((post50-2) %% 4) +  table(pre50%% 4)