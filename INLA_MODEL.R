
#################################################################
# LOAD PACKAGE 
#################################################################
library(tidyverse)## Versatile package for data analysis
library(survey)
library(srvyr)
library(purrr)
library(haven)
set.seed(1) # Set the seed is important for getting reproducible reports 
## Clean the envorment 
# rm(list=ls())
options(scipen=4)

#################################################################
# LOAD DATA
#################################################################
library(foreign)
data <- read_dta("ZA7562_v1-0-0.dta")
#data <- read.dta("ZA7562_v1-0-0.dta")
data <- zap_labels(data)


data.f <- data %>%
  dplyr::select(
    caseid, 
    country, 
    isocntry, 
    d11, 
    q1_1:q1_30,
    qc1_1:qc12_8,
    p7be:p7hu_r,
    d1:d77,
    nuts, nutslvl,
    w1:w86
  )

data.f <- as.data.frame(data.f)  

#################################################################
# Create the vaccine hesitancy variable
#################################################################

# Question qc7 
# c7_1 Vaccines overload and weaken the immune system 1==TRUE
# c7_1 Vaccines can cause the disease against which they protect 1==TRUE
# c7_1 Vaccines can often produce serious side-effects 1==TRUE
# c7_1 Vaccines are rigorously tested before being authorised for use 2==FALSE
data.f <- data.f %>%
  mutate (
    # At least one positive hesitancy response (high or low vaccine hesitancy)
    hesitancy_vac1 = (
      ifelse( qc7_1==1 | qc7_2==1 | qc7_3==1 | qc7_4==2, 1,0)
    ),
    # All positive hesitancy response (higher vaccine hesitancy)
    hesitancy_vac2 = (
      ifelse( qc7_1==1 & qc7_2==1 & qc7_3==1 & qc7_4==2, 1,0)
    )
  )

#################################################################
# SET SURVEY DATA 
#################################################################

library(survey)
dstrat_srvyr <- svydesign(ids=~1,
                          weights=~w23 ,
                          strata=~nuts, 
                          data=data.f)


#################################################################
# Small area estimation 
#################################################################


data.f <- subset(data, !is.na(data.f$hesitancy_vac2))
data.f <- subset(data, !is.na(data.f$nuts))
names(data.f)[names(data.f)=='wex'] <- 'strata'
n.area <- length(unique(data.f$nuts))

########################
## Weighted estimates
########################
p.i <- svyby(~hesitancy_vac2,
             ~nuts,dstrat_srvyr,
             svymean)$hesitancy_vac2

dv.i <- svyby(~hesitancy_vac2,
              ~nuts,
              dstrat_srvyr,
              svymean)$se^2

logit.pi <- log(p.i/(1-p.i))
v.i <- dv.i/(p.i^2*(1-p.i)^2)

########################
## Construct data frame for INLA
########################

data <- matrix(NA,nrow=n.area,ncol=1)
data <- as.data.frame(data)
colnames(data)[1] <- "unstruct"
data$hracode <- as.character(unique(data.f$nuts))
data$p.i <- p.i
data$dv.i <- dv.i
data$v.i <- v.i
data$logit.pi <- logit.pi
data$logit.prec <- 1/v.i
data <- data[order(data$hracode),]
data$unstruct <- 1:(n.area)
data$struct <- 1:(n.area)

########################
## LOAD MAP
########################
EU_NUTS <- readOGR(dsn = "NUTS_2.shp", layer = "NUTS_RG_01M_2016_4326_LEVL_2")

########################
## SET INLA
########################

library(spdep)
EU_NUTS.neigh <- poly2nb(EU_NUTS)
library(INLA)
nb2INLA('EU_NUTS.neigh.graph',EU_NUTS.neigh)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.10")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("Rgraphviz")

library(Rgraphviz)


g = inla.read.graph('EU_NUTS.neigh.graph')
summary(g)
plot(g) 

########################
## Fit global/local spatial smoothing model
########################

formula = logit.pi ~ 1 + 
  f(struct, model='besag', 
    adjust.for.con.comp=TRUE, 
    constr=TRUE, 
    graph='EU_NUTS.neigh.graph') + 
  f(unstruct,model='iid', 
    param=c(0.5,0.008))

mod.smooth <- inla(formula, 
                   family = "gaussian", 
                   data = data, 
                   control.predictor = list(compute = TRUE),
                   control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                   scale=logit.prec)

### DOES NOT CONVERGE 

# Question there are some NUTS that have ZERO prevalence, does that affect the model stability ?


########################
## THE END
########################
