#################################################################
#################################################################
########## SHOCKS x FARM SIZE ##########
#################################################################
#################################################################

##############
### SET UP ###
##############
getwd()

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(nloptr)
library(lme4)
library(ggeffects)
library(matrixStats) # For weightedMedian
library(DHARMa)
library(mgcViz)
library(brms)

#################
### FUNCTIONS ###
#################

afn.findssu <- function(data_allcountries) {
  
  SSU <- data_allcountries$LVL3 # start at the smallest possible unit
  index <- SSU=="" # for all hh that do not have LVL3, move up to LVL2
  SSU <- replace(SSU, index, data_allcountries$LVL2[which(index)]) # replace by LVL2
  index <- SSU=="" # repeat
  SSU <- replace(SSU, index, data_allcountries$LVL1[which(index)]) # replace by LVL1
  index <- SSU==""
  SSU <- replace(SSU, index, data_allcountries$LVL0.5[which(index)])
  index <- SSU==""
  SSU <- replace(SSU, index, data_allcountries$LVL0[which(index)])
  
  data_allcountries <- cbind(data_allcountries, SSU)
  
  return(data_allcountries)
  
}

nlopt <- function(par, fn, lower, upper, control) {
  .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                            opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                        maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
  list(par = res$solution,
       fval = res$objective,
       conv = if (res$status > 0) 0 else res$status,
       message = res$message
  )
}

# Function to calculate weighted quantiles
weighted_quantile <- function(x, w, probs) {
  if (length(x) != length(w)) {
    stop("'x' and 'w' must have the same length")
  }
  o <- order(x)
  x <- x[o]
  w <- w[o]
  Sw <- cumsum(w)
  totalW <- sum(w)
  cum_weights <- Sw / totalW
  sapply(probs, function(p) {
    if (any(cum_weights >= p)) {
      min(x[cum_weights >= p])
    } else {
      NA  # Handle cases where p is greater than all cum_weights
    }
  })
}

weighted_mean <- function(x, w) {
  return(sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE))
}


#################
### LOAD DATA ###
#################

##need to grab data from repo and set up public repo##

# microdata <- readRDS("microshock/2.SPEI/2.Output/microdata_spei_merge.RDS")
microdata <- readRDS("microdata_spei_merge.RDS")
num_samples <- 20 # there are 1122 unique lvl2's for Colombia..
set.seed(1234)

#################
### CHECK DATA ###
#################

getdistros<-function(dat){
  
  # Function to calculate weighted mean
  
  # Assuming your data frame is named 'microdata'
  microdata_prepared <- dat %>%
    mutate(HHWGT = ifelse(is.na(HHWGT) | HHWGT == 0, 1, HHWGT)) # Set NA or 0 weights to 1
  
  summary_data <- microdata_prepared %>%
    group_by(LVL0) %>%
    summarize(
      weighted_median_fsize = weightedMedian(FSIZE, w = HHWGT, na.rm = TRUE),
      weighted_mean_fsize = weighted_mean(FSIZE, w = HHWGT),
      min_fsize = min(FSIZE, na.rm = TRUE), # min/max are not typically weighted
      max_fsize = max(FSIZE, na.rm = TRUE),
      p10_fsize = weighted_quantile(FSIZE, w = HHWGT, probs = 0.10),
      p90_fsize = weighted_quantile(FSIZE, w = HHWGT, probs = 0.90)
    ) %>%
    mutate(
      range_text = paste0("Range: ", format(min_fsize, digits = 2), " - ", format(max_fsize, digits = 2), " ha"),
      weighted_median_text = paste0("Median: ", format(weighted_median_fsize, digits = 2), " ha"),
      weighted_mean_text = paste0("Mean: ", format(weighted_mean_fsize, digits = 2), " ha"),
      pctl_text = paste0("P10-P90: ", format(p10_fsize, digits = 2), " - ", format(p90_fsize, digits = 2), " ha")
    )
  
  d<-ggplot(microdata_prepared, aes(x = FSIZE)) +
    geom_histogram(aes(y = ..density..), bins = 30) +
    facet_wrap(~LVL0, scales = "free_x") +
    scale_x_log10() +
    theme_bw() +
    labs(x = "Farm Size (ha)", y = "Density") +
    theme(legend.position = "none") +
    geom_text(
      data = summary_data,
      aes(x = Inf, y = Inf, label = weighted_mean_text),
      hjust = 1, vjust = 2, size = 3
    ) +
    geom_text(
      data = summary_data,
      aes(x = Inf, y = Inf, label = weighted_median_text),
      hjust = 1, vjust = 3.5, size = 3
    ) +
    geom_text(
      data = summary_data,
      aes(x = Inf, y = Inf, label = range_text),
      hjust = 1, vjust = 5, size = 3
    ) +
    geom_text(
      data = summary_data,
      aes(x = Inf, y = Inf, label = pctl_text),
      hjust = 1, vjust = 6.5, size = 3
    )
  print(d)
  return(summary_data)
  
}

getdistros(microdata)

#some countries seem very off. number of reasons why this could be. Sampling variation...
#Also, not all these surveys may capture all very small farms (!) or commercial large ones (!)
#Here we correct 4 countries that are clearly wrong: Guatemala, Bangledesh, Ethiopia", "Burkina Faso

#first remove rows from microdatathat are Guatemala and also FSIZE >=1000
#this seems to fix this country 
microdata.gfix <- microdata %>% filter(!(LVL0 %in% c("Guatemala") & FSIZE >= 1000))
getdistros(microdata.gfix) #yes Guatemala seems more on point

#WCA mean FS https://www.fao.org/faostat/en/#data/WCAD in ha
LVL0<-c("Bangladesh", "Ethiopia", "Burkina Faso")
meanfs<-c(0.52, 1.03, 4.18	)
year<-c(2019, 2000, 2010)
wcamean<-data.frame(LVL0, meanfs, year)

#then divide the farm size by the mean farm size
sm<-getdistros(microdata.gfix)
ratio_off<-sm %>% filter(LVL0 %in% c("Bangladesh", "Ethiopia", "Burkina Faso"))  %>%
  left_join(wcamean, by = "LVL0") %>%
  mutate(ratio=weighted_mean_fsize/meanfs)  %>%
  select (LVL0,  ratio) 

#now run through all of microdata.gfix and divide the FSIZE by the ratio for
#these three countries in ratio_off 
microdata.allfix<-microdata.gfix %>% left_join(ratio_off, by = "LVL0") %>%
  mutate(FSIZE=ifelse(!is.na(ratio), FSIZE/ratio, FSIZE)) %>%
  select(-ratio)

#ok seems to fix.
#some countries prob still don't match census, 
#but at least in sensible range
getdistros(microdata.allfix)

######################
### BINARY COMBINED###
######################

#filter data
binary_loss <- microdata.allfix %>%
  select(HHID, HHWGT, SRAU, LVL0, LVL0.5, LVL1, LVL2, LVL3, FSIZE, FGROUP, PRODF, PRODD) %>%
  pivot_longer( cols = starts_with("PROD"),
    names_to = "Event",
    values_to = "Loss" ) %>%
  mutate(Event = ifelse(Event == "PRODF", "Flood", "Drought"),
  Loss=Loss*1) %>%
  filter(!is.na(Loss)) %>%
  mutate(
    FSIZE=ifelse(FSIZE==0,0.0001,FSIZE), #this is to avoid log(0)
    LogFSIZE=log(FSIZE),
    HHWGT=round(HHWGT,0), #need integer for modelling
    HHWGT=ifelse(is.na(HHWGT) |HHWGT==0,1,HHWGT) #here the HH, but some of these seem super off tbh..
    )


#Colombia option...prob exclude b/c of large size of data
#could return to this later, with prop random sampling across outcome
binary_loss<-afn.findssu(binary_loss)
binary_loss_nc<-binary_loss %>% filter(LVL0 != "Colombia")
binary_loss_sc<-afn.samplecolombia(binary_loss)

#check out continuous functions
model_binary_con <- glmer(Loss ~LogFSIZE*Event+ (1|LVL0/SRAU),
                      family=binomial,
                      data=binary_loss_nc,
                      control = glmerControl(optimizer = "nlopt", calc.derivs = TRUE))

#save model
saveRDS(model_binary_con, "Out/model_binary_con.RDS")

#plot
eff_bin <- ggpredict(model_binary_con, terms = c("LogFSIZE", "Event"))
plot(eff_bin)
#save this plot
ggsave("Out/model_binary_con.png", width = 8, height = 6)

#some checks.. dist seems off..
simulationOutput <- simulateResiduals(fittedModel = model_binary, plot = T)
testOutliers(model_binary_con, type = 'bootstrap')

######################
###CONT COMBINED###
######################

#filter data
cont_loss <- microdata.allfix %>%
  select(HHID,HHWGT,SRAU,LVL0,LVL0.5,LVL1,LVL2,LVL3,FSIZE,FGROUP,HARVF,LOSTF,REVPF,HARVD,LOSTD,REVPD) %>%
  pivot_longer(
    cols = starts_with("REV"),
    names_to = "Event",
    values_to = "Loss" ) %>%
  mutate(
    Event = ifelse(Event == "REVPF", "Flood", "Drought"),
    Loss=Loss*1) %>%
  filter(!is.na(Loss)) %>%
  filter(!(FSIZE %in% c(0, 451973.9))) %>%
  mutate(
    LogFSIZE=log(FSIZE),
    HHWGT=round(HHWGT,0),
    HHWGT=ifelse(is.na(HHWGT) |HHWGT==0,1,HHWGT) #here the HH, but some of these seem super off tbh..for example ETH has 1's but clearly not full microdata
  )

cont_loss$Lost.Revenue<-as.numeric(recode(as.character(cont_loss$Loss), 
                                                        "0" = "0.00001", "1"= "0.9999")) # revalue upper and lower
#cut up with or without colombia...
cont_loss_nc<-cont_loss %>% filter(LVL0 != "Colombia")
cont_loss_sc<-afn.samplecolombia(cont_loss)

#check out continuous functions
model_con <- glmmTMB(Lost.Revenue ~LogFSIZE*Event+ (1|LVL0/SRAU/SSU),
                          family=beta_family(link="logit"),
                       # weights=HHWGT,
                   cont_loss_nc)
summary(model_con )
saveRDS(model_con, "Out/model_con.RDS")

#plot
eff_con <- ggpredict(model_con, terms = c("LogFSIZE", "Event"))
plot(eff_con)
#save this plot
ggsave("Out/model_con.png", width = 8, height = 6)

#Model check https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
#distributional assumptions incorrect.
simulationOutput <- simulateResiduals(fittedModel = model_con, plot = T)
testOutliers(model_con, type = 'bootstrap') #not massively problematic

#Get country lists
cont_loss_nc %>% group_by(LVL0, Event) %>% summarise()

binary_loss_nc %>% group_by(LVL0, Event) %>% summarise()%>%print(n=30)

######################
###QGAM example###
######################

#check to see if overall effects are massively different...
#convert Event, LVL0 and SRAU to factors all together
cont_loss_nc$LVL0<-as.factor(cont_loss_nc$LVL0)
cont_loss_nc$SRAU<-as.factor(cont_loss_nc$SRAU)
cont_loss_nc$Event<-as.factor(cont_loss_nc$Event)

cont_qgam<-qgamV(Lost.Revenue ~ LogFSIZE+
               LogFSIZE* Event + 
               s(LVL0, bs = "re")+
               s(SRAU, bs = "re"),
               qu = 0.5,
              data = cont_loss_nc)


#plot effects/ run diag
print(plot(cont_qgam, allTerms = TRUE),
      pages = 1)

eff_con_qgam <- ggpredict(cont_qgam, terms = c("LogFSIZE", "Event"))
plot(eff_con_qgam)


#save this plot
ggsave("Out/model_con_gam.png", width = 8, height = 6)
saveRDS(cont_qgam, "Out/cont_qgam.RDS")
o<-getViz(cont)
check(o) + l_gridQCheck2D(qu = 0.5)

check1D(o, c("LogFSIZE")) +
  l_gridQCheck1D(qu = 0.5)

######################
######brms########
######################

library(brms)
brms_model_binary_con <- brm(
  formula = Loss ~ LogFSIZE * Event + (1 | LVL0 / SRAU),
  family = bernoulli(), # Bernoulli family for binary outcome (0/1)
  data = binary_loss_nc,
  control = list(adapt_delta = 0.8)
)

# You can then examine the model summary
summary(brms_model_binary_con)

# And get credible intervals for the parameters
posterior_summary(brms_model_binary_con)

saveRDS(brms_model_binary_con, "Out/brms_bin.RDS")

# You can also use ggpredict with the brms object
library(ggeffects)
predicted_values <- ggpredict(brms_model_binary_con, terms = c("LogFSIZE", "Event"))
plot(predicted_values)
ggsave("Out/model_bin_brms.png", width = 8, height = 6)
