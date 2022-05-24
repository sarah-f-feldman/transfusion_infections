library(dplyr)
library(ggplot2)
library(tidyr)
library(xlsx)
library(lubridate)
library(stringr)

d <- read.csv2("data/MNH Blood donor spreadsheet.csv")

# d1 <- read.xlsx("data/MNH Blood donor spreadsheet_ys_MAJ20171020.xlsx", sheetIndex = 1)
# d2 <- read.xlsx("data/MNH Blood donor spreadsheet_ys_MAJ20171020.xlsx", sheetIndex = 2)
# d3 <- bind_rows(d1, d2)
# saveRDS(d3, "data/MNH_20171024.rds")

d <- readRDS("data/MNH_20171024.rds")


#===========================
#Data management
#===========================

#-----------------
#doublons
d %>% group_by(Reg) %>% mutate(n = n()) %>% filter(n>=2) %>% arrange(Reg) %>% 
  write.table(., file="clipboard", sep= "\t", row.names = F)

#exclure les doublons (les 2 exemplaires du doublons)
d <- d %>% group_by(Reg) %>% mutate(n = n()) %>% filter(n==1) %>% ungroup

#-----------------
#cleaning Donor.status
d <- d %>% mutate(Donor.status = tolower(Donor.status),
                  Donor.status = gsub(" ", "", Donor.status),
                  Donor.status = ifelse(Donor.status == "", NA, Donor.status),
                  Donor.status = recode(Donor.status, "always" = "repeat", "always(regular)" = "repeat", "donorclub" = "repeat"))
table(d$Donor.status)

#-----------------
#cleaning Universal.donor
d <- d %>% mutate(Universal.donor = as.character(Universal.donor), 
                  Universal.donor = ifelse(Universal.donor == "", NA, Universal.donor))

#exclusion Universal.donor == "Conversion" : les rajouter en analyse de sensibilite : 
d <- d %>% filter(!Universal.donor %in% "Conversion")

#-----------------
#infection et cofinection
#creation column infection oui/non
d <- d %>% mutate(infection = apply(d %>% select(HIV, HBV, HCV, Syphilis), 1, max)) 

#creation column coinfection oui/non
d <- d %>% mutate(coinfection = apply(d %>% select(HIV, HBV, HCV, Syphilis), 1, sum),
                  coinfection = ifelse(coinfection > 1, 1, 0)) #ok car pas de NA pour le status viral
table(d$coinfection)

#pattern de coinf
d <- d %>% mutate(quel_coinf = paste(HIV, HBV, HCV, Syphilis, sep = "-")) 
#----------------------------------
#Dates

d$Date.collected <- as.Date(d$Date.collected,  "%d/%m/%Y")
d$Date.donated <- as.Date(d$Date.donated,  "%d/%m/%Y")

#je supprime les dates donated anterieures à 2016 (erreur de saisie) 
d <- d %>% filter(Date.donated > as_date("2016-01-01")) %>% #exclu en même temps les dates manquantes
  #j'exclus le mois de septembre 2016
  filter(!months(Date.donated) %in% "septembre")

#----------------------------------
#Rajout de Age_cat
d <- d %>% mutate(Age_cat = cut(Age, c(min(Age, na.rm = T), 25, 35, 45, max(Age, na.rm = T)), 
                                include.lowest = T, #pour que 76 soit inclus dans la borne 45,76 (sinon NA)
                                right = F)) #pour que l'intervalle soit ferme à gauche


#saveRDS(d, "data/MNH_20171024_clean.rds")

#version long
dl <- d %>% gather(key = virus, value = status, infection, HIV, HBV, HCV, Syphilis)



#===========================
# Analyse
#===========================
d <- readRDS("data/MNH_20171024_clean.rds")
#20181206 
d <- data.frame(d)
#fin20181206 
virus_vec <- c("infection", "HCV", "HBV", "HIV", "Syphilis")
#-----------------
#Nombre de patients analyses
nrow(d)

#taux d'infection :
round(prop.table(table(d$infection))*100,2) #8%

#coinfection
#population totale
round(prop.table(table(d$coinfection))*100,2) #0.47%
#parmi ceux qui sont infectes
round(prop.table(table(d$coinfection[d$infection==1]))*100,2) #5.55%
#pattern de coinf
table(d$quel_coinf[d$coinfection == 1])

#-----------------
#Description Age et Sexe : table 1 (on n'utilise que la ligne var NA pour age et sexe tot et les lignes universal donor)

my_des <- bind_rows(
  
  #general
  d %>% summarise(n = n(), sex_M = sum(Sex), sex_p = sex_M/n, age_m = mean(Age, na.rm = T), age_sd = sd(Age, na.rm = T)) %>% 
    mutate(var = NA, group = NA, pval_sexe = NA, pval_age = NA),
  #donor staus
  d %>% group_by(Donor.status) %>% 
    summarise(n = n(), sex_M = sum(Sex), sex_p = sex_M/n, age_m = mean(Age, na.rm = T), age_sd = sd(Age, na.rm = T)) %>% 
    mutate(var = "Donor.status", group = Donor.status, 
           pval_sexe = chisq.test(d$Sex, d$Donor.status)$p.value, 
           pval_age = t.test(d$Age[d$Donor.status=="firsttime"],d$Age[d$Donor.status=="repeat"])$p.value),
  #universal donor
  d %>% group_by(Universal.donor) %>% 
    summarise(n = n(), sex_M = sum(Sex), sex_p = sex_M/n, age_m = mean(Age, na.rm = T), age_sd = sd(Age, na.rm = T)) %>% 
    mutate(var = "Universal.donor", group = Universal.donor,
           pval_sexe = chisq.test(table(d$Sex, d$Universal.donor))$p.value, #na non pris en compte 
           pval_age = t.test(d$Age[d$Universal.donor=="0"],d$Age[d$Universal.donor=="1"])$p.value) 
  
) %>% 
  #mise en forme
  select(var, group, n:pval_age) %>% 
  mutate_at(vars(contains("_p")), funs(round(.*100,2))) %>% 
  mutate_at(vars(contains("age")), funs(round(.,2))) %>% 
  mutate_at(vars(contains("pval")), funs(round(.,3)))

write.table(data.frame(my_des), file = "clipboard", sep = "\t", row.names = F)



#---------------------
#first time among family and voluntary
table(d$Donor.status)
prop.table(table(d$Donor.status))
n <- table(d$Donor.status, d$Universal.donor)
pr <- round(prop.table(table(d$Donor.status, d$Universal.donor), 2)*100,1)
chisq.test(n)
data.frame(matrix(paste0(n, " (", pr, "%)"), nrow = 2))
#-----------------
#Prevalence dans la population

prop.table(table(d$HCV))
.s <- sapply(virus_vec, function(virus){
  n <- nrow(d)
  M <- sum(d[ ,virus])
  #20181206 M <- sum(pull(d[ ,virus]))
  btest <- binom.test(x = M, n = n)
  p0 <- as.numeric(btest$estimate)
  IC <- as.numeric(btest$conf.int)
  IC <- round(c(p0, IC)*100, 2) # en %
  return(c(virus, n, M, IC))
  # p0 <- M/n
  # varp0 <-  p0*(1-p0)/n
  # IC_low <-  p0 - 1.96 * sqrt(varp0)
  # IC_up <-  p0 + 1.96 * sqrt(varp0)
  # cond <- all((IC_low * n) >=5 , (IC_up * n) >=5,((1-IC_low) * n) >=5, ((1-IC_up) * n) >=5)
  #IC <- round(c(p0, IC_low, IC_up)*100, 2) # en %
  #return(c(virus, n, M, IC, cond))
})
prev <- data.frame(t(.s))
#colnames(prev) <- c("virus", "Number of Patients", "infected patients (n)", "prevalence(%)", "95IC_l", "95IC_up", "IC_valid") #en %
colnames(prev) <- c("virus", "Number of Patients", "infected patients (n)", "prevalence(%)", "95IC_l", "95IC_up") #en %

write.table(data.frame(prev), file = "clipboard", sep = "\t", row.names = F)
#-----------------
#Prevalence selon universal donor ou non 
#family (universal donor = 0)

.l <- lapply(c(0,1), function(univ){
  .s <- lapply(virus_vec, function(virus){
    tmp <- d[d$Universal.donor==univ & !is.na(d$Universal.donor), ] 
    n <- nrow(tmp)
    M <- sum(tmp[ ,virus], na.rm = T)
    btest <- binom.test(x = M, n = n)
    p0 <- as.numeric(btest$estimate)
    IC <- as.numeric(btest$conf.int)
    IC <- round(c(p0, IC)*100, 2) # en %
    return(c(univ, virus, n, M, IC))
    # p0 <- M/n
    # varp0 <-  p0*(1-p0)/n
    # IC_low <-  p0 - 1.96 * sqrt(varp0)
    # IC_up <-  p0 + 1.96 * sqrt(varp0)
    # cond <- all((IC_low * n) >=5 , (IC_up * n) >=5,((1-IC_low) * n) >=5, ((1-IC_up) * n) >=5)
    # IC <- round(c(p0, IC_low, IC_up)*100, 2) # en %
    # return(c(univ, virus, n, M, IC, cond))
  })
  .s <- data.frame(do.call(rbind,.s))
  return(.s)
})
prev_universal <- data.frame(do.call(cbind, .l))
#colnames(prev_universal) <- rep(c("Universal donor (0: family; 1: voluntary)", "virus", "number of patients (N)", "patient infected (N)", "prevalence(%)", "95IC_l", "95IC_up", "IC_valid"),2)
colnames(prev_universal) <- rep(c("Universal donor (0: family; 1: voluntary)", "virus", "number of patients (N)", "patient infected (N)", "prevalence(%)", "95IC_l", "95IC_up"),2)

write.table(data.frame(prev_universal), file = "clipboard", sep = "\t", row.names = F)

#-----------------------
#prevalence by sex and age
#20181206 prevalence by sex, age, donor.status, universal.donor : colonnes prevalence de la table 2

.m <- lapply(virus_vec, function(virus){
  .l <- lapply(c("Sex", "Age_cat", "Donor.status", "Universal.donor"), function (var){
    my_levels <- levels(as.factor(d[ , var]))
    .s <- lapply(my_levels, function(lev){
      print(c(var, lev))
      tmp <- d[d[ ,var]==lev & !is.na(d[ ,var]), ] 
      n <- nrow(tmp)
      M <- sum(tmp[ ,virus], na.rm = T)
      btest <- binom.test(x = M, n = n)
      p0 <- as.numeric(btest$estimate)
      IC <- as.numeric(btest$conf.int)
      IC <- round(c(p0, IC)*100, 2) # en %
      return(c(virus, var, lev, n, M, IC))
    })
    do.call(rbind, .s)
  })
  do.call(rbind, .l)
})
prev_by <- data.frame(do.call(rbind, .m))
colnames(prev_by) <- c("virus", "variable", "level", "N", "M", "prevalence(%)", "IC95_low", "IC95_up")
write.table(data.frame(prev_by), file = "clipboard", sep = "\t", row.names = F)


#-----------------
# Test : 

#Principal hypothesis : Universal Donors are less infected than family
chisq.test(table(d$infection, d$Universal.donor), correct = F) #significant

#Secondary hypothesis : seperate tests for HBV, HCV, HIV, Syphilis
chisq.test(table(d$HBV, d$Universal.donor), correct = F) #significant
chisq.test(table(d$HCV, d$Universal.donor), correct = F) #not significant
chisq.test(table(d$HIV, d$Universal.donor), correct = F) #significant
chisq.test(table(d$Syphilis, d$Universal.donor), correct = F) #not significant

p_virus <- do.call(rbind, lapply(virus_vec, function(virus){
  tab <- table(d[ ,virus], d$Universal.donor)
  #verif cond
  mg <- addmargins(tab)
  v1 <- mg[row.names(mg)== "Sum",][1:2]
  v2 <- mg[,colnames(mg)== "Sum"][1:2]
  cond <- all(all(v1[1]*v2/max(mg)>5), all(v1[2]*v2/max(mg)>5))
  #chi2 sans correction
  p <- round(chisq.test(tab, correct = F)$p.value,3)
  return(data.frame(virus, cond, p))
}))
p_virus
write.table(data.frame(p_virus), file = "clipboard", sep = "\t", row.names = F)

# Exploratory : Repeated donors are less infected
chisq.test(table(d$infection, d$Donor.status), correct = F) #not significant

do.call(rbind, lapply(virus_vec, function(virus){
  tab <- table(d[ ,virus], d$Donor.status)
  #verif cond
  mg <- addmargins(tab)
  v1 <- mg[row.names(mg)== "Sum",][1:2]
  v2 <- mg[,colnames(mg)== "Sum"][1:2]
  cond <- all(all(v1[1]*v2/max(mg)>5), all(v1[2]*v2/max(mg)>5))
  #chi2 sans correction
  p <- round(chisq.test(tab, correct = F)$p.value,3)
  return(data.frame(virus, cond, p))
}))


#-----------------
#Je merge le test selon universal donors et prevalence
cbind(cbind(prev, prev_universal), p_virus[ ,c("virus", "p")])



#===========================
# Figures


#Virus  
tmp <- dl %>% group_by(virus, status) %>% count() %>% 
  group_by(virus) %>% mutate(perc = nn/sum(nn)) %>% filter(status == 1)

ggplot(tmp, aes(virus, perc, fill = virus)) +
  geom_bar(stat = "identity")


#Donor status
ggplot(d, aes(Donor.status, fill = Donor.status)) +
  geom_bar()

#Donor status en fonction de universal
ggplot(d, aes(Donor.status, fill = Universal.donor)) +
  geom_bar()

ggplot(d, aes(Donor.status, fill = Universal.donor)) +
  geom_bar() +
  facet_wrap(~infection)

#infection function of Donor status 
ggplot(d, aes(Donor.status, fill = as.factor(infection))) +
  geom_bar()

#infection function of Universal or family donor 
ggplot(d, aes(Universal.donor, fill = as.factor(infection))) +
  geom_bar()


#Age_cat en fonction du statut viral pour chaque virus(en pourcentage)
tmp <- dl %>% group_by(virus, status, Age_cat) %>% count() %>% 
  group_by(virus, status) %>% mutate(perc = nn/sum(nn)) 

ggplot(tmp, aes(x=Age_cat, y = perc, fill = Age_cat))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.5))+
  facet_grid(virus~status) +
  geom_text(aes(label = sprintf("%.2f%%", perc * 100)), 
            vjust = -.5)


#Sexe en fonction du statut viral pour chaque virus(en pourcentage)
tmp <- dl %>% group_by(virus, status, Sex) %>% count() %>% 
  group_by(virus, status) %>% mutate(perc = nn/sum(nn)) 

ggplot(tmp, aes(x=as.factor(Sex), y = perc, fill = as.factor(Sex)))+
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1))+
  facet_grid(virus~status) +
  geom_text(aes(label = sprintf("%.2f%%", perc * 100)), 
            vjust = -.5)

#===========================
#Modèle logistique

#----------------------
#univariate logistique

get_univar2 <- function(my_data, my_outcome, my_var){
  #my_data = d2
  #my_var <- "Age_cat"
  #my_var <- factors_notif[2]
  #my_outcome <- "HIV"
  my_data <- as.data.frame(my_data) #sinon probleme pour levels(factor(my_data[ ,my_var]))
  my_data <- na.omit(my_data[ ,c(my_outcome, my_var)])
  print(my_var)
  
  #modele
  my_formula <- as.formula (paste0(my_outcome," ~ ", my_var))
  mod1 <- glm(my_formula, data = my_data, family = "binomial")
  g <- summary(mod1)
  
  #recuperer la pvalue (selon variable categorielle ou non)
  if(length(levels(factor(my_data[ ,my_var])))>2){
    my_formula <- as.formula (paste0(my_outcome," ~ 1"))
    mod0 <- glm(my_formula, data = my_data, family = "binomial")
    pval <- anova(mod1, mod0, test="LRT")
    pval <- pval$`Pr(>Chi)`[2]
  } else {
    tab <- g$coefficients  
    pval <- tab[grepl(my_var, rownames(tab)), "Pr(>|z|)"]
  }
  
  #recuperer OR et IC
  pval <- round(pval, 3)
  OR <- exp(coef(mod1))
  CI <- exp(confint(mod1))
  ORCI <- round(cbind(OR, CI), 2)
  
  #organiser le tableau
  res <- data.frame(outcome = my_outcome, variable = my_var, levels = NA, ORCI, pval, stringsAsFactors = FALSE)
  res <- res[grepl(my_var, rownames(res)), ]
  levels <- gsub(my_var, "", rownames(res))
  levels <- if(any(levels != "")) levels else "1"
  res$levels <- levels
  names(res) <- c("outcome", "variable", "levels", "OR", "IC_l", "IC_u", "pvalue")
  return(res)
}

all_univar <- lapply(virus_vec, function(my_virus){
  lapply(c("Age_cat", "Sex", "Donor.status", "Universal.donor"), function(my_var){
    get_univar2(my_data = d, my_outcome = my_virus, my_var = my_var)  
  }) %>% bind_rows()
}) %>% bind_rows()

write.table(data.frame(all_univar), file = "clipboard", sep = "\t", row.names = F)


#----------------------
#multivariate logistique

#sans forcer de variables : c'est la colonne adjusted OR de la table 2

get_multivar <- function(my_formula, virus){
  my_formula <- formula (my_formula)
  mod1 <- glm(my_formula, data = d, family = "binomial")
  #mod1 <- glm(HBV ~ Age_cat + Sex + Donor.status + Universal.donor, data = d, family = "binomial")
  g <- summary(mod1)
  tab <- data.frame(g$coefficients)
  tab$lev <- NA
  tab$lev[-1] <- sapply(rownames(tab)[-1], function(lev){
    #lev <- rownames(tab)[2]
    my_string <- str_extract(lev, colnames(d))[!is.na(str_extract(lev, colnames(d)))]
    my_string <-names(sapply(my_string, nchar)[sapply(my_string, nchar)==max(as.numeric(sapply(my_string, nchar)))])
    gsub(my_string, "", lev)
  })
  tab$var<-NA
  tab$var[1] <- rownames(tab)[1]
  tab$var[-1] <- sapply(rownames(tab)[-1], function(lev){
    #lev <- rownames(tab)[2]
    my_string <- str_extract(lev, colnames(d))[!is.na(str_extract(lev, colnames(d)))]
    my_string <-names(sapply(my_string, nchar)[sapply(my_string, nchar)==max(as.numeric(sapply(my_string, nchar)))])
    return(my_string)
  })
  rownames(tab) <- NULL
  colnames(tab) <- c("coef", "se","Zvalue", "Pvalue", "level", "variable")
  tab <- tab[-1,]
  tab[, c("coef", "se", "Zvalue", "Pvalue")] <- round(tab[, c("coef", "se", "Zvalue", "Pvalue")], 3)
  cib <- round(cbind(exp(coef(mod1)), exp(confint(mod1))), 2)[-1, ]
  cib <- if(nrow(tab)==1)data.frame(t(cib)) else data.frame(cib) ; colnames(cib) <- c("OR", "95IC_l", "95IC_u")
  tab <- cbind(tab, cib)
  tab <- tab[,c("variable", "level", "OR", "95IC_l", "95IC_u", "coef", "se", "Zvalue", "Pvalue")]
  nevent <- sum(d[!rownames(d)  %in% g$na.action ,virus], na.rm = T)
  cond <- paste0("conditions : number of event = ", nevent, " ; number of variable = ",nrow(tab), " ; ", nrow(tab), " x ", "10 = ", nrow(tab)*10, " have to be inferior to ", nevent)
  print(list(virus, cond))
  tab$virus <- virus
  tab <- tab[,c("virus", "variable", "level", "OR", "95IC_l", "95IC_u", "coef", "se", "Zvalue", "Pvalue")]
  return(tab)
}

.l <- lapply(virus_vec, function(my_virus){
  print(my_virus)
  
  #20181206 my_var <- all_univar %>% filter(virus == my_virus & pval < 0.05) %>% select(variable) %>% pull
  my_var <- all_univar %>% filter(outcome == my_virus & pvalue < 0.05) %>% select(variable) %>% pull
  #fin 20181206
  my_var <- paste(unique(my_var), collapse = " + ")
  print(paste0(my_virus, "~",  my_var))
  get_multivar(paste0(my_virus, "~",  my_var), my_virus)
})
.l <- do.call(rbind, .l)
write.table(.l, file = "clipboard", sep = "\t", row.names = F)



#en forçant Universal.donor
lapply(virus_vec, function(my_virus){
  print(my_virus)
  my_var <- all_univar %>% filter(virus == my_virus & pval < 0.05) %>% select(variable) %>% pull
  my_var <- c("Universal.donor", my_var)
  my_var <- paste(unique(my_var), collapse = " + ")
  get_multivar(paste0(my_virus, "~",  my_var), my_virus)
})



# #avec backward selection ?
# #NON Yusuke me dit qu'il n'utilise la selection que pour les modèles predictifs
# # et pas pour les modèles explicatifs. De plus là comme beaucoup de patients, on prend 0.05 comme seuil.
# my_outcome <- "infection"
# my_outcome <- "HBV"
# my_outcome <- "HCV"
# my_outcome <- "HIV"
# my_outcome <- "Syphilis"
# vars_mm <- all_univar %>% filter(outcome == my_outcome & pvalue <0.2) %>% pull(variable) %>% unique
# vars_mm_paste <- paste(vars_mm, collapse = " + ")
# my_formula <- as.formula(paste0(my_outcome, " ~ ", vars_mm_paste))
# mod1 <- glm(my_formula, data = na.omit(d[,c(my_outcome, vars_mm)]), family = "binomial")
# my_step <- step(object = mod1, direction = "backward")
# g <- summary(my_step)
# select <- lapply(vars_mm, function(myvar) str_extract(names(my_step$coefficients), myvar)) %>% unlist %>% unique %>% na.omit %>% as.vector

#----------------------
#multivariate logistique
#Significance of Age_cat : complète la table 2 : pvalue de age_cat qd age cat dans multivarie

#infection
mod1 <- glm(as.formula("infection ~ Universal.donor + Age_cat "), data = d[!is.na(d$Age_cat), ], family = "binomial")
mod0 <- glm(as.formula("infection ~ Universal.donor"), data = d[!is.na(d$Age_cat), ], family = "binomial")
anova(mod1, mod0, test = "LRT")
L1 <- logLik(mod1)
L0 <- logLik(mod0)
k <- nlevels(d$Age_cat)
pval <- round(1-pchisq(2*(L1-L0),df=k-1), 3)

#HCV : no Age_cat
[1] "HCV~Sex"

#HBV
[1] "HBV~Age_cat + Sex + Donor.status + Universal.donor"
mod1 <- glm(as.formula("HBV ~ Age_cat + Sex + Donor.status + Universal.donor"), data = d[!is.na(d$Age_cat), ], family = "binomial")
mod0 <- glm(as.formula("HBV ~ Sex + Donor.status + Universal.donor"), data = d[!is.na(d$Age_cat), ], family = "binomial")
anova(mod1, mod0, test = "LRT")
L1 <- logLik(mod1)
L0 <- logLik(mod0)
k <- nlevels(d$Age_cat)
pval <- round(1-pchisq(2*(L1-L0),df=k-1), 3)

[1] "HIV~Age_cat + Universal.donor"
mod1 <- glm(as.formula("HIV ~ Age_cat + Universal.donor"), data = d[!is.na(d$Age_cat), ], family = "binomial")
mod0 <- glm(as.formula("HIV ~ Universal.donor"), data = d[!is.na(d$Age_cat), ], family = "binomial")
anova(mod1, mod0, test = "LRT")
L1 <- logLik(mod1)
L0 <- logLik(mod0)
k <- nlevels(d$Age_cat)
pval <- round(1-pchisq(2*(L1-L0),df=k-1), 3)

[1] "Syphilis~Age_cat"
mod1 <- glm(as.formula("Syphilis~Age_cat"), data = d[!is.na(d$Age_cat), ], family = "binomial")
mod0 <- glm(as.formula("Syphilis~1"), data = d[!is.na(d$Age_cat), ], family = "binomial")
anova(mod1, mod0, test = "LRT")
L1 <- logLik(mod1)
L0 <- logLik(mod0)
k <- nlevels(d$Age_cat)
pval <- round(1-pchisq(2*(L1-L0),df=k-1), 3)





