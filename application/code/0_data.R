# preliminar ####
rm(list=ls())

librerias <- c('stringr','dplyr','tidyverse','reshape2','haven','eeptools')
sapply(librerias, require, character.only=T)
# sapply(librerias, install.packages, character.only=T)

setwd('~/Desktop/#Classes/4(24)_G0A16a Master Thesis/data')



# 1. load ####

# standardized assessment
file_dir1 = str_replace(file_dir1, '_N2015_', '_N2017_')
n2017s = read_delim(file_dir1, delim=';', locale=locale(encoding="ISO-8859-1"))
# str(n2017s)

# decentralized stage
file_dir2 = str_replace(file_dir2, '_N2015', '_N2017')
n2017d = read_delim(file_dir2, delim=';', locale=locale(encoding="ISO-8859-1"))
# str(n2017d)


# 2. transformation ####

n2017s_mod = n2017s


## 2.1 sex ####
# nominal dichotomous variable 

n2017s_mod$sex = n2017s_mod$sexo
n2017s_mod$sex = ifelse(n2017s_mod$sex=='FEMENINO', 'female', 'male')
# table(n2017s_mod$sex, useNA='ifany')



## 2.2 age ####
# (at the time of evaluation)
# ordered continuous variable

eval_date = as.Date('28/05/2017', format="%d/%m/%Y")
# https://evaluaciondocente.perueduca.pe/nombramiento-contratacion-2017/cronograma/

# n2017s_mod$fecha_nac = as.Date(n2017s_mod$fecha_nac, format="%d/%m/%Y")
n2017s_mod$age = age_calc(n2017s_mod$fecha_nac, enddate=eval_date, units="years")
n2017s_mod$age = floor(n2017s_mod$age)
# table(n2017s_mod$age, useNA='ifany')




## 2.3 INEI_location ####
# nominal polytomous variable 

n2017s_mod$INEI_location = n2017s_mod$UBIGEO_INEI
# table(n2017s_mod$INEI_location, useNA='ifany')



## 2.4 education ####
# ordered (or nominal) polytomous variable 

n2017s_mod$education = n2017s_mod$procedencia
n2017s_mod$education = str_replace(n2017s_mod$education, 'Solo INS', 'institute only')
n2017s_mod$education = str_replace(n2017s_mod$education, 'Ambos', 'both')
n2017s_mod$education = str_replace(n2017s_mod$education, 'Solo UNIV', 'university only')
# table(n2017s_mod$education, useNA='ifany')



## 2.5 edu_specialty ####
# nominal polytomous variable 

n2017s_mod$edu_specialty = n2017s_mod$mod
# table(n2017s_mod$edu_specialty, useNA = 'ifany')



## 2.6 public_exp ####
# ordered polytomous variable

n2017s_mod$public_exp = n2017s_mod$exp_publica
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'Sin experiencia', 'no experience')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'Menos de 1 año', 'less than 5y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 1 a 2 años', 'less than 5y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 3 a 5 años', 'less than 5y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 6 a 10 años', 'between 6y and 10y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 11 a 15 años', 'more than 10y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 16 a 20 años', 'more than 10y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'De 21 a 25 años', 'more than 10y')
n2017s_mod$public_exp = str_replace(n2017s_mod$public_exp, 'Más de 25 años', 'more than 10y')
# table(n2017s_mod$public_exp, useNA = 'ifany')




## 2.7 private_exp ####
# ordered polytomous variable

# 2017
n2017s_mod$private_exp = n2017s_mod$exp_privada
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'Sin experiencia', 'no experience')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'Menos de 1 año', 'less than 5y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 1 a 2 años', 'less than 5y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 3 a 5 años', 'less than 5y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 6 a 10 años', 'between 6y and 10y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 11 a 15 años', 'more than 10y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 16 a 20 años', 'more than 10y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'De 21 a 25 años', 'more than 10y')
n2017s_mod$private_exp = str_replace(n2017s_mod$private_exp, 'Más de 25 años', 'more than 10y')
# table(n2017s_mod$private_exp, useNA = 'ifany')





## 2.9 disability ####
# nominal polytomous variable

n2017s_mod$disability = n2017s_mod$nombre_discapacidad
n2017s_mod$disability = str_replace(n2017s_mod$disability, 'NINGUNA', 'none')
n2017s_mod$disability = str_replace(n2017s_mod$disability, 'MOTORA', 'motor')
n2017s_mod$disability = str_replace(n2017s_mod$disability, 'CEGUERA', 'blindness')
n2017s_mod$disability = str_replace(n2017s_mod$disability, 'BAJA VISIÓN', 'low vision')
n2017s_mod$disability = str_replace(n2017s_mod$disability, 'AUDITIVA', 'auditory')
# table(n2017s_mod$disability, useNA = 'ifany')



## 2.10 final_condition ####
# nominal (or ordinal) polytomous variable

n2017s_mod$final_condition = n2017s_mod$COND_FINAL_EB.17
idx = str_detect(n2017s_mod$final_condition, 'Ganadores de plaza')
n2017s_mod$final_condition[idx] = 'nominated'
idx = str_detect(n2017s_mod$final_condition, 'Clasificados no ganadores')
n2017s_mod$final_condition[idx] = 'classified'
idx = str_detect(n2017s_mod$final_condition, 'No clasificados')
n2017s_mod$final_condition[idx] = 'non-classified'

idx1 = n2017s_mod$final_condition=='classified'
idx2 = n2017s_mod$contratado2018=='Si'
idx2 = ifelse(is.na(idx2), FALSE, idx2)
idx = idx1 & idx2
n2017s_mod$final_condition[idx] = 'classified-contract'

idx1 = n2017s_mod$final_condition=='non-classified'
idx2 = n2017s_mod$contratado2018=='Si'
idx2 = ifelse(is.na(idx2), FALSE, idx2)
idx = idx1 & idx2
n2017s_mod$final_condition[idx] = 'non-classified-contract'

idx1 = str_detect(n2017s_mod$final_condition, 'No aplica')
idx2 = n2017s_mod$contratado2018=='Si'
idx2 = ifelse(is.na(idx2), FALSE, idx2)
idx = idx1 & idx2
n2017s_mod$final_condition[idx] = 'non-classified-contract'

idx = str_detect(n2017s_mod$final_condition, 'No aplica')
n2017s_mod$final_condition[idx] = 'non-classified'

# table(n2017s_mod$final_condition, useNA = 'ifany')
# table(n2017s_mod$contratado2018, useNA = 'ifany')
# with(n2017s_mod, table(final_condition, contratado2018, useNA = 'ifany'))




## 2.11 evaluation_region ####
# nominal polytomous variable

n2017s_mod$evaluation_region = n2017s_mod$region_eval
# table(n2017s_mod$evaluation_region, useNA = 'ifany')



## 2.12 contract_region ####
# nominal polytomous variable

n2017s_mod$contract_region = n2017s_mod$reg_contratado
# table(n2017s_mod$contract_region, useNA = 'ifany')



# 3. join assessment with decentralized ####
# it includes:
#   - nomination_region: nominal polytomous variable
#   - classObs_score: continuous
#   - classObs_pass: ordinal dichotomous variable

idx = with(n2017d, ganador_plaza=='SI' & !is.na(ganador_plaza) )
data_mom = n2017d[idx, c('id_nomb','region','pje_observacion','aprueba_obsaula')]

n2017s_join = merge(n2017s_mod, data_mom, by='id_nomb', all.x=T)
names(n2017s_join)[90:92] = c('nomination_region','classObs_score','classObs_pass')
n2017s_join$classObs_score = as.numeric(n2017s_join$classObs_score)
# str(n2017s_join)
# str(n2017d)
#
# idx = n2017s_mod$final_condition=='nominated'
# sum(n2017s_mod$id_nomb[idx] %in% data_mom$id_nomb)



# 4. variables of interest ####

## 4.1 renaming items ####

idx_var = str_detect(names(n2017s_join), '^s[:digit:]')
names(n2017s_join)[idx_var] = str_replace(names(n2017s_join)[idx_var], '^s', 'N17_')



## 4.2 selecting variables ####

idx_var = names(n2017s_join)[c(1,79:92,17:41)]
n2017s_clean = n2017s_join[, idx_var]
names(n2017s_clean)[1] = 'ID_ind'
# str(n2017s_clean)
# names(n2017s_clean)
# str(n2017s_join)



## 4.3 export ####

file_save = file.path(getwd(), '3_transformed', '1_N2017_clean.csv')
write_csv(n2017s_clean, file=file_save)




# 5. long format ####


## 5.1 with missings ####

## load
file_save = file.path(getwd(), '3_transformed', '1_N2017_clean.csv')
n2017s_clean = read_csv(file_save)
# str(n2017s_clean)


n2017s_long = melt(n2017s_clean, id.vars=names(n2017s_clean)[1:15],
                   factorsAsStrings=FALSE)
names(n2017s_long)[16:17] = c('item','response')
n2017s_long$item = as.character(as_factor(n2017s_long$item))

file_save = file.path(getwd(), '3_transformed', '2_N2017_long_miss.csv')
write_csv(n2017s_long, file=file_save)

# str(n2017s_clean)
# str(n2017s_long)
# with(n2017s_long, table(item, response, useNA = 'ifany'))






## 5.2 without missings ####

file_save = file.path(getwd(), '3_transformed', '2_N2017_long_miss.csv')
n2017s_long = read_csv(file_save)
# str(n2017s_clean)

n2017s_long = n2017s_long[!is.na(n2017s_long$response),]
file_save = file.path(getwd(), '3_transformed', '3_N2017_long.csv')
write_csv(n2017s_long, file=file_save)
# with(n2017s_long, table(item, response, useNA = 'ifany'))