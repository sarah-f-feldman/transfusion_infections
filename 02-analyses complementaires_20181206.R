#===========================
# DM complementaires
#===========================

library(dplyr)
library(ggplot2)
library(tidyr)
library(xlsx)
library(lubridate)
library(stringr)


#20180207

d <- readRDS("data/MNH_20171024_clean.rds")
#d <- readRDS("data/MNH_20171024.rds")

#----------------------------------
#creer variable ds.subgroup(donor.Status sub group)

#first time family donors : ft_fd
#repeat family donors : r_fd
#first time voluntary donors : ft_vd
#repeat voluntary donors : r_vd
#Universal donor (0: family; 1: voluntary)
d$ds.subgroup <- NA
d$ds.subgroup[d$Donor.status == "firsttime" & d$Universal.donor == 0 & !is.na(d$Donor.status) & !is.na(d$Universal.donor)] <- "ft_fd"
d$ds.subgroup[d$Donor.status == "repeat" & d$Universal.donor == 0 & !is.na(d$Donor.status) & !is.na(d$Universal.donor)] <- "r_fd"
d$ds.subgroup[d$Donor.status == "firsttime" & d$Universal.donor == 1 & !is.na(d$Donor.status) & !is.na(d$Universal.donor)] <- "ft_vd"
d$ds.subgroup[d$Donor.status == "repeat" & d$Universal.donor == 1 & !is.na(d$Donor.status) & !is.na(d$Universal.donor)] <- "r_vd"

#verif
table(d$ds.subgroup, useNA = "a")
table(d$Donor.status, d$Universal.donor, useNA = "a")

#-------------------------
#analyses parmi d$Returned.for.results == 1
d$Returned.for.results <- ifelse(is.na(d$Returned.for.results), 0, d$Returned.for.results)

saveRDS(d, "data/MNH_20180207.rds")


#===========================
# Analyse complementaires
#===========================
d <- readRDS("data/MNH_20180207.rds")

#-----------------
#Prevalence selon subgroup : 4 premieres colonnes de la table 3
virus_vec <- c("infection", "HBV", "HCV", "HIV", "Syphilis")

my_prevalence <- function(data){
  #faire un vecteur virus_vec avec le nom des variables contenant l'information infection (0/1)
  prop.table(table(data$HCV))
  .s <- sapply(virus_vec, function(virus){
    n <- nrow(data)
    M <- sum(pull(data[ ,virus]))
    btest <- binom.test(x = M, n = n)
    p0 <- as.numeric(btest$estimate)
    IC <- as.numeric(btest$conf.int)
    IC <- round(c(p0, IC)*100, 2) # en %
    return(c(virus, n, M, IC))
  })
  prev <- data.frame(t(.s))
  #colnames(prev) <- c("virus", "Number of Patients", "infected patients (n)", "prevalence(%)", "95IC_l", "95IC_up", "IC_valid") #en %
  colnames(prev) <- c("virus", "Number of Patients", "infected patients (n)", "prevalence(%)", "95IC_l", "95IC_up") #en %
  return(prev)
}

prev <- my_prevalence(d[d$ds.subgroup == "ft_fd" & !is.na(d$ds.subgroup), ])
prev <- my_prevalence(d[d$ds.subgroup == "r_fd" & !is.na(d$ds.subgroup), ])
prev <- my_prevalence(d[d$ds.subgroup == "ft_vd" & !is.na(d$ds.subgroup), ])
prev <- my_prevalence(d[d$ds.subgroup == "r_vd" & !is.na(d$ds.subgroup), ])

prev$CI <- paste0(prev$`95IC_l`, " ; ", prev$`95IC_up`)
prev$`95IC_l` <- NULL 
prev$`95IC_up` <- NULL
write.table(data.frame(prev), file = "clipboard", sep = "\t", row.names = F)


#verif des resultats precedents
prev <- my_prevalence(d[d$Universal.donor == 1 & !is.na(d$Universal.donor), ])

#-----------------
#Difference de prevalence selon subgroup

# chi2
p_virus <- do.call(rbind, lapply(virus_vec, function(virus){
  tab <- table(pull(d[ ,virus]), d$ds.subgroup)
  #verif cond
  mg <- addmargins(tab)
  v1 <- mg[row.names(mg)== "Sum",][-ncol(mg)]
  v2 <- mg[,colnames(mg)== "Sum"][-nrow(mg)]
  cond <- all(sapply(1:length(v2), function(num)v1*as.numeric(v2[num])) / max(mg)>5)
  cond <- sapply(1:length(v2), function(num)v1*as.numeric(v2[num])) / max(mg)
  all_cond_sup5 <- all(cond>5)
  all_cond_sup3 <- all(cond>3) 
  cond <- ifelse(all_cond_sup5 == TRUE, "sup 5", ifelse(all_cond_sup3== TRUE, "sup3", "inf 3"))
  #cond <- all(all(v1[1]*v2/max(mg)>5), all(v1[2]*v2/max(mg)>5))
  #chi2 sans correction
  my_correct <- ifelse(cond != "sup 5", TRUE, FALSE)
  p <- round(chisq.test(tab, correct = my_correct)$p.value,3)
  return(data.frame(virus, cond, p))
}))
p_virus
write.table(data.frame(p_virus), file = "clipboard", sep = "\t", row.names = F)

# fisher : ce sont les pvalues utilisees pour table 3
p_virus_fisher <- do.call(rbind, lapply(virus_vec, function(virus){
  tab <- table(pull(d[ ,virus]), d$ds.subgroup)
  p <- round(fisher.test(tab)$p.value,3)
  return(data.frame(virus, p))
}))
p_virus_fisher
write.table(data.frame(p_virus_fisher), file = "clipboard", sep = "\t", row.names = F)


#-----------------
#%notif :table 4

my_perc_notify <- function(mydata){
  #faire un vecteur virus_vec avec le nom des variables contenant l'information infection (0/1)
  .s <- sapply(virus_vec, function(virus){
    N <- nrow(mydata[mydata[ ,virus] == 1, ])
    notif <- nrow(mydata[mydata[ ,virus] == 1 & d$Returned.for.results == 1, ])
    btest <- binom.test(x = notif, n = N)
    p0 <- as.numeric(btest$estimate)
    IC <- as.numeric(btest$conf.int)
    IC <- round(c(p0, IC)*100, 2) # en %
    return(c(virus, N, notif, IC))
  })
  prev <- data.frame(t(.s))
  #colnames(prev) <- c("virus", "Number of Patients", "infected patients (n)", "prevalence(%)", "95IC_l", "95IC_up", "IC_valid") #en %
  colnames(prev) <- c("virus", "Number of infected(N)", "notified patients (n)", "prevalence(%)", "95IC_l", "95IC_up") #en %
  return(prev)
}

prev <- my_perc_notify(d)
prev$CI <- paste0(prev$`95IC_l`, " ; ", prev$`95IC_up`)
prev$`95IC_l` <- NULL 
prev$`95IC_up` <- NULL
write.table(data.frame(prev), file = "clipboard", sep = "\t", row.names = F)


#---------------------
#factor associated with notification : # 20181206 pas utilise apparemment (et ne tourne pas...)
factors_notif <- c("Age_cat", "Sex", "Donor.status", "Universal.donor", "HBV","HCV", "HIV", "Syphilis")

table(d$Returned.for.results, useNA = "a")
table(d$infection)
d2 <- d[d$infection == 1, ]

levels <- c("family" = 0, "voluntary" = 1)
d2$Universal.donor <- factor(d2$Universal.donor, levels = levels, labels = names(levels))
levels <- c("Female" = 0, "Male" = 1)
d2$Sex <- factor(d2$Sex, levels = levels, labels = names(levels))


#prevalence by sex and age
my_perc_notify_byvar <- function(mydata, myvar){
  # mydata <- d2
  # myvar <- "Age_cat"
  myvec <- factor(d2[ ,myvar])
  .l <- lapply(levels(myvec), function(lev){
    N <- nrow(d2[myvec == lev, ])
    notif <- nrow(d2[myvec == lev & d2$Returned.for.results == 1, ])
    btest <- binom.test(x = notif, n = N)
    p0 <- as.numeric(btest$estimate)
    IC <- as.numeric(btest$conf.int)
    IC <- round(c(p0, IC)*100, 2) # en %
    return(c(myvar, lev, N, notif, IC))
  })
  prev <- data.frame(do.call(rbind, .l))
  colnames(prev) <- c("variable", "levels", "Number of infected(N)", "notified patients (n)", "prevalence(%)", "95IC_l", "95IC_up") #en %
  return(prev)
}

.l<-lapply(factors_notif, function(var){
  my_perc_notify_byvar(d2, var)  
})
prev_by<-do.call(rbind, .l)
write.table(data.frame(prev_by), file = "clipboard", sep = "\t", row.names = F)

# modele univarie

get_univar <- function(my_data, my_var){
  #my_data = d2
  #my_var <- "Age_cat"
  #my_var <- factors_notif[2]
  my_data <- as.data.frame(my_data) #sinon probleme pour levels(factor(my_data[ ,my_var]))
  print(my_var)
  
  #modele
  my_formula <- as.formula (paste0("Returned.for.results ~ ", my_var))
  mod1 <- glm(my_formula, data = my_data, family = "binomial")
  g <- summary(mod1)
  
  #recuperer la pvalue (selon variable categorielle ou non)
  if(length(levels(factor(my_data[ ,my_var])))>2){
    my_formula <- as.formula (paste0("Returned.for.results ~ 1"))
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
  res <- data.frame(variable = my_var, levels = NA, ORCI, pval, stringsAsFactors = FALSE)
  res <- res[grepl(my_var, rownames(res)), ]
  levels <- gsub(my_var, "", rownames(res))
  levels <- if(any(levels != "")) levels else "1"
  res$levels <- levels
  names(res) <- c("variable", "levels", "OR", "IC_l", "IC_u", "pvalue")
  return(res)
}


.l <- lapply(factors_notif, function(x) get_univar(my_data = d2, my_var = x))
all_univar <- bind_rows(.l)
write.table(all_univar, file = "clipboard", sep = "\t", row.names= FALSE)

# modele multivarie avec les variables significatives
all_univar %>% filter(pvalue<0.2)

selection <- step(object = glm(Returned.for.results ~ Sex + HBV + HIV, data = d2, family = "binomial"), direction = "backward")
names(selection$coefficients)
mod1 <- glm(Returned.for.results ~ Sex + HBV , data = d2, family = "binomial")
g <- summary(mod1)
tab <- g$coefficients  
pval <- tab[, "Pr(>|z|)"]
#recuperer OR et IC
pval <- round(pval, 3)
OR <- exp(coef(mod1))
CI <- exp(confint(mod1))
ORCI <- round(cbind(OR, CI), 2)

#organiser le tableau
res <- data.frame(ORCI, pval, stringsAsFactors = FALSE)
write.table(res, file = "clipboard", sep = "\t", row.names= TRUE)

#-------------
#HIV vs autre infection : # 20181206 pas utilise apparemment 
# d2$HIV comme c'est d2 qui ne contient que les sujets infectes alors HIV = 1 
# signifie sujet infecte par VIH (confinf et monoinf) et HIV = 0 signifie autres infections

mod1 <- step(object = glm(Returned.for.results ~ Sex + HIV, data = d2, family = "binomial"), direction = "backward")
g <- summary(mod1)
tab <- g$coefficients  
pval <- tab[, "Pr(>|z|)"]
#recuperer OR et IC
pval <- round(pval, 3)
OR <- exp(coef(mod1))
CI <- exp(confint(mod1))
ORCI <- round(cbind(OR, CI), 2)
#organiser le tableau
res <- data.frame(ORCI, pval, stringsAsFactors = FALSE)
write.table(res, file = "clipboard", sep = "\t", row.names= TRUE)

#-------------------------
# 20140214
# demande de test for linear trend pour age cat et notification
x = as.numeric(table(d2$Age_cat, d2$Returned.for.results)[,2])
n = as.numeric(table(d2$Age_cat))
prop.trend.test(x, n)

