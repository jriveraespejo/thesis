# preliminar ####
rm(list=ls())

librerias <- c('stringr','dplyr','tidyverse','reshape2','haven','eeptools')
sapply(librerias, require, character.only=T)
# sapply(librerias, install.packages, character.only=T)

setwd('~/Desktop/application/data')
# getwd()


# 1. load ####
file_save = file.path(getwd(), '3_transformed', '3_N2017_long.csv')
n2017s_long = read_csv(file_save)
# str(n2017s_long)



# 3. add ####

# identifying texts
n2017s_long$text = NA

idx = n2017s_long$item %in% paste0('N17_0',1:5)
n2017s_long$text[idx] = 1

idx = n2017s_long$item %in% c(paste0('N17_0',6:9), 'N17_10')
n2017s_long$text[idx] = 2

idx = n2017s_long$item %in% paste0('N17_',11:15)
n2017s_long$text[idx] = 3

idx = n2017s_long$item %in% paste0('N17_',16:20)
n2017s_long$text[idx] = 4

idx = n2017s_long$item %in% paste0('N17_',21:25)
n2017s_long$text[idx] = 5

# table(n2017s_long$text, useNA='ifany')



# identifying dimensions
n2017s_long$dimension = NA

idx = n2017s_long$item %in% c(paste0('N17_0',c(1,6)), paste0('N17_',c(11,16,18,22,23)))
n2017s_long$dimension[idx] = 1

idx = n2017s_long$item %in% c(paste0('N17_0',c(2:5,7:9)), paste0('N17_',c(12:13,17,19:20)))
n2017s_long$dimension[idx] = 2

idx = n2017s_long$item %in% paste0('N17_',c(10,14:15,21,24:25) )
n2017s_long$dimension[idx] = 3

# table(n2017s_long$dimension, useNA='ifany')




# 3. sample ####
set.seed(56985)
idx_ID = sample( unique(n2017s_long$ID_ind), size=2000)
n2017s_final = n2017s_long[n2017s_long$ID_ind %in% idx_ID, ]

file_save = file.path(getwd(), '4_final', '4_N2017_final.csv')
write_csv(n2017s_final, file=file_save)







# 4. data prep ####

# individuals indexes
IDind_o = unique(n2017s_final$ID_ind)
IDind_o = IDind_o[order(IDind_o)]
extra1 = data.frame(IDind_o=IDind_o, IDind_m=1:length(IDind_o))
n2017s_final = merge(n2017s_final, extra1, by.x='ID_ind', by.y='IDind_o', all.x=T)

# items indexes
IDitem_o = unique(n2017s_final$item)
IDitem_o = IDitem_o[order(IDitem_o)]
extra2 = data.frame(IDitem_o=IDitem_o, IDitem_m=1:length(IDitem_o))
n2017s_final = merge(n2017s_final, extra2, by.x='item', by.y='IDitem_o', all.x=T)

# public experience
extra3 = data.frame(Xpu_o=unique(n2017s_final$public_exp), Xpu_m=4:1)
n2017s_final = merge(n2017s_final, extra3, by.x='public_exp', by.y='Xpu_o', all.x=T)

# private experience
extra4 = data.frame(Xpr_o=unique(n2017s_final$private_exp), Xpr_m=c(2,1,3,4))
n2017s_final = merge(n2017s_final, extra4, by.x='private_exp', by.y='Xpr_o', all.x=T)

# disability
extra5 = data.frame(D_o=unique(n2017s_final$disability), D_m=c(1,2,3,4))
n2017s_final = merge(n2017s_final, extra5, by.x='disability', by.y='D_o', all.x=T)
# table(n2017s_final[, c('disability','D_m')])


# sorting data
idx = with(n2017s_final, order(IDind_m, IDitem_m, text, dimension))
n2017s_final = n2017s_final[idx,]
# str(n2017s_final)


# list
data_app = list(
  
  # evaluation data
  N = nrow(n2017s_final),
  J = length(unique(n2017s_final$ID_ind)),
  K = length(unique(n2017s_final$item)),
  L = length(unique(n2017s_final$text)),
  D = length(unique(n2017s_final$dimension)),
  IDj = n2017s_final$IDind_m,
  IDk = n2017s_final$IDitem_m,
  IDl = n2017s_final$text,
  IDd = n2017s_final$dimension,
  G = ifelse(n2017s_final$sex=='female',2,1),
  A = n2017s_final$age,
  E = with(n2017s_final, ifelse(education=='institute only', 1, 
                                ifelse(education=='university only',2,3))),
  S = with(n2017s_final, ifelse(edu_specialty=='EBR-INI', 1, 
                                ifelse(edu_specialty=='EBR-PRI',2,3))),
  Xpu = n2017s_final$Xpu_m,
  Xpr = n2017s_final$Xpr_m,
  Di = n2017s_final$D_m,
  y = n2017s_final$response
)


individual = unique(data.frame(data_app[c('IDj','G','A','E','S','Xpu','Xpr','Di')]))

data_app$IDind = individual$IDj
data_app$GE = individual$G
data_app$AG = individual$A
data_app$ED = individual$E
data_app$SP = individual$S
data_app$XPpu = individual$Xpu
data_app$XPpr = individual$Xpr
data_app$DI = individual$Di

items = unique(data.frame(data_app[c('IDk','IDl','IDd')]))

data_app$IDitem = items$IDk
data_app$IDtext = items$IDl
data_app$IDdim = items$IDd


save(data_app, file=file.path(getwd(), '4_final', 'data_application.RData') )

