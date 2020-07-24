# this does the simulation
#rm(list=ls(all=TRUE))
library(rstan)
library(tibble)
library(tidyverse)
library(readxl)
library(parallel)
library(foreach)
library(doParallel)
library(abind)
cl = makeCluster(8)
registerDoParallel(cl)
##Covid-data Nigeria
#source("Covid_data_processing.R")
source("Covid_data_processing_May7.R")

#covert date to integer for all cases
int_date <- as.numeric(covid_nig_data_all$my_date)
int_date_count <- c(int_date[1],diff(int_date))
int_date_count[1] = 0
covid_nig_data_all$day = cumsum(int_date_count)

##get local cases
local_dat <- covid_nig_data %>% filter(source == "Local")


##get import cases
import_dat <- covid_nig_data %>% filter(source == "Imported")
## find day integer for corressponding date
day_count<- function(data_withonlydate, data_withinteger){
  len_var <- length(data_withinteger$my_date)
  result = rep(0,len_var)
  for (i in seq(1,len_var)){
    if(sum(data_withinteger$my_date[i] == data_withonlydate$my_date)>=1){
      result[i] = data_withonlydate$Incidence[data_withinteger$my_date[i] == 
                                                data_withonlydate$my_date]
    }
  }
  return(result)
}

local_inc <- day_count(local_dat,covid_nig_data_all)
import_inc <- day_count(import_dat,covid_nig_data_all)
##daily report is missing and we adjust for it by
max_day <- max(covid_nig_data_all$day)
Tmain <-  covid_nig_data_all$day
T_seq = seq(0,max_day)
daily_report<-function(Tmain,Tseq,local_inc,import_inc){
  daily_local = rep(0,length(Tseq))
  daily_import = rep(0,length(Tseq))
  mycount = rep(0,length(Tmain))
  counter =1
  for (i in seq(1,length(Tseq))){
    if(sum(Tseq[i]==Tmain)>=1){
      daily_local[i] = local_inc[Tseq[i]==Tmain]
      daily_import[i] = import_inc[Tseq[i]==Tmain]
      mycount[counter] = Tseq[i]+1
      counter = counter + 1
    }
  }
  return(list(daily_local = daily_local,daily_import = daily_import,mycounter = mycount))
}

mydaily_report <-daily_report(Tmain,T_seq,local_inc,import_inc)

loc_case = mydaily_report$daily_local
import_case =  mydaily_report$daily_import


###death data
###readin death datat
###
#dat_death<- read_xlsx("Under_reporting_Nig.xlsx")
dat_death<- read_xlsx("Under_reporting_Nig_May7.xlsx")
dat_death<-arrange(dat_death, dateRep)

start_date<- covid_nig_data_all$my_date[1]
date_starter_sim<-as.Date(c("15/03/2020"),"%d/%m/%Y")

date_starter<-as.Date(c("10/03/2020"),"%d/%m/%Y")
my_death_dat <- dat_death%>% filter(dateRep >= date_starter) %>%
  mutate(cumdeath = cumsum(deaths),daycount= seq(0,(length(deaths)-1)))
my_death_rep <- my_death_dat %>% filter(dateRep >= date_starter_sim)
#make it end on the 26 of April
#my_death_rep = my_death_rep[-length(my_death_rep$dateRep),]  Not needed now
my_death_rep$dateRep<- as.Date(my_death_rep$dateRep)
# my_death_dat %>% ggplot(aes(x=d
my_date_loc <-start_date+ 0:(length(loc_case)-1)
dyn_data <- data.frame(my_date = my_date_loc,local_case = loc_case,
                       cumloccase = cumsum(loc_case),cumimport = cumsum(import_case))

my_dyn_data <- dyn_data %>% filter(my_date >= date_starter_sim)
my_dyn_data<- mutate(my_dyn_data, cumdeath = cumsum(my_death_rep$deaths))
N = length(my_dyn_data$my_date) - 1
ts = 1:N
y_me_init = c(1.98e8-2,0,0,0,3,0,0,0,2,1)

y_me = as.matrix(my_dyn_data[2:(N + 1), 3:5])

y_me = cbind(y_me[ , 3], y_me[ , 2],y_me[ , 1]); # here,data
covid_fitdata = list(N=N, ts=ts, y_me_init=y_me_init, y_me=y_me)
#define initial values for each chain 
initf <- function() {
  list(theta = c(0.344,0.938,0.977,0.13),sigma=c(0.2,64,65))
}
fit = stan(file='modelstan.stan', data=covid_fitdata, iter=5000,
           chains = 4,control=list(adapt_delta=0.99),init=initf)
# 
print(fit, pars=c("theta", "sigma"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
# print(fit, pars=c("theta", "sigma"),
#       probs=c(0.025,0.5,0.975), digits = 3)
# 
# theta= extract(fit,'theta')[[1]]
# sigma = extract(fit,'sigma')[[1]]
# library(ggplot2)
# qplot(theta[,1], fill=I("blue"))


###########solve ode 
model_covid<-function(t,z,theta) {
   s = z[1]
   e1 = z[2]
   e2 = z[3]
   i1 = z[4]
   i2 = z[5]
   i3 = z[6]
   r =  z[7]
   d =  z[8]
   c =  z[9]
   v =  z[10]
    
   beta = theta[1]
   sigma1 = 0.3125
   sigma2 = 0.5
   rho1 = 0.5
   rho2 = theta[2]
   rho3 = theta[3]
   gamma1 = 0.5
   gamma2 = 0.176
   gamma3 = 0.5
   k1 = theta[4]
   n1  = 198000000
    
   ds_dt = - beta*(i1+2*i2+3*i3)*s/n1
   de1_dt = beta*(i1+2*i2+3*i3)*s/n1 - sigma1*e1
   de2_dt = sigma1*e1 - sigma2*e2
   di1_dt = rho1*sigma2*e2 - gamma1*i1
   di2_dt = k1*c + gamma1*i1-gamma2*i2
   di3_dt = (1-rho1)*sigma2*e2 - gamma3*i3
   dr_dt = rho2*gamma2*i2 + gamma3*i3;
   dd_dt = (1-rho2)*gamma2*i2;
   dc_dt = k1*c;
   dv_dt = rho3*(k1*c + gamma1*i1);
   dxdt = c(ds_dt,de1_dt,de2_dt,di1_dt,di2_dt,di3_dt,dr_dt,dd_dt,dc_dt,dv_dt)
    return (list(dxdt))
  }
library(deSolve)
params = c(0.09, 0.946, 0.966, 0.848, 0.965)
times <- 0:N # returns a sequence

out <- as.data.frame(
  ode(
    func=model_covid,
    y=y_me_init,
    times=times,
    parms=params
  )
)

names(out) = c("time","S","E1","E2","I1","I2","I3","R","D","C","V")
load('RstanMay7.RData')# I run got this result from running the code from clusters

print(fit, pars=c("theta", "sigma"),
      probs=c(0.025,0.5,0.975), digits = 3)
theta= rstan::extract(fit,'theta')[[1]]
#pairs(fit,pars=c("theta", "sigma"))
plot(fit,pars=c("theta", "sigma"))
n_size = dim(theta)[1]
my_time = 0:(N+21)
k_R0<- foreach (j = 1:n_size)%dopar%{
   library(deSolve)
  
   out <- as.data.frame(
      ode(
         func=model_covid,
         y=y_me_init,
         times=my_time,
         parms=theta[j,]
      )
   )
}

###get quantiles
###get confidence interval
#range
ci_cal_new<-function(timet,data_me){
   ci_out <-matrix(0, nrow =timet,ncol =4)
   for (j in seq(from=1, to=timet,by=1)){
      yboot =quantile(data_me[j,],c(0.025,0.25,0.75,0.975),na.rm = TRUE)
      ci_out[j,1] = yboot[[1]]
      ci_out[j,2] = yboot[[2]]
      ci_out[j,3] =yboot[[3]]
      ci_out[j,4] = yboot[[4]]
   }
   return(ci_out)
}

###cumulative local trans getter
CI_new<-function(k_travel, sim_num, time_len,det){
   #det determines if I want death or cumulative local cases
   time_dim = length(time_len)
   cal_data_cum = matrix(0,nrow = time_dim, ncol = sim_num)
   lower_est=rep(0,time_dim)
   median_est= rep(0,time_dim)
   upper_est= rep(0,time_dim)
   for (sim in seq(1,sim_num)){ #k is simulation output
      if (det == 1){
        cal_data_cum[,sim] = k_travel[[sim]]$'8'}
      else if(det == 2){
         cal_data_cum[,sim] = k_travel[[sim]]$'10'}
   }
   i_cal<- ci_cal_new(time_dim,cal_data_cum)
   lower5p =i_cal[,1]
   lower25p = i_cal[,2]
   upper75p= i_cal[,3]
   upper95p= i_cal[,4]
   return( data.frame(days = time_len,lower5p = lower5p,lower25p = lower25p,
                      upper75p=upper75p,upper95p=upper95p))
}
######
sim_result_D <-CI_new(k_R0,n_size,my_time,1)
sim_result_L <-CI_new(k_R0,n_size,my_time,2)
#get time in date for plot
my_date = date_starter_sim + my_time
sim_result_D = mutate(sim_result_D, my_date = my_date)
sim_result_L = mutate(sim_result_L, my_date = my_date)
###Plotting####
####Scenario analysis
require(scales)

p1<- ggplot(dyn_data, aes(x=my_date, y=cumloccase))+
   geom_point(size = 3, colour = "black")+
   geom_ribbon(data = sim_result_L,aes(ymin=lower5p,ymax=upper95p,
                                     x=my_date,fill= "r1"),inherit.aes = FALSE,alpha=0.3)+
   geom_ribbon(data = sim_result_L,aes(ymin=lower25p,ymax=upper75p,
                                     x=my_date,fill= "r2"),inherit.aes = FALSE,alpha=0.3)+
   geom_vline(xintercept=as.Date(c("7/05/2020"),"%d/%m/%Y"),
              linetype="dashed", color = "red",size=1)+
   scale_x_date(date_breaks = "2 week", date_labels = "%b %d")+
   #scale_colour_manual(values=c("dodgerblue2"), name="",breaks=c("r1"))+
   scale_fill_manual(values=c("dodgerblue1","dodgerblue2"), name="",breaks=c("r1","r2"))+
   labs(x = "Days", y=paste("Cumulative reported cases"))+
   coord_cartesian(ylim = c(1,NA),expand = FALSE)+
   theme_classic()+
   theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
         axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
         legend.position = c(0.2,0.8),legend.key = element_rect(colour = "transparent", fill = "white"),
         legend.background = element_rect(fill = "transparent",colour =NA),
         axis.line = element_line(color = "black", 
                                 size = 1, linetype = "solid"),
         axis.ticks = element_line(size = 1))+
   guides(colour= FALSE,fill=FALSE)

p2<- ggplot( my_death_rep, aes(x=dateRep, y=cumsum(deaths)))+
   geom_point(size = 3, colour = "black")+
   geom_ribbon(data = sim_result_D,aes(ymin=lower5p,ymax=upper95p,
                                       x=my_date,fill= "r1"),inherit.aes = FALSE,alpha=0.3)+
   geom_ribbon(data = sim_result_D,aes(ymin=lower25p,ymax=upper75p,
                                       x=my_date,fill= "r2"),inherit.aes = FALSE,alpha=0.3)+
   geom_vline(xintercept=as.Date(c("07/05/2020"),"%d/%m/%Y"),
              linetype="dashed", color = "red",size=1)+
   scale_x_date(date_breaks = "2 week", date_labels = "%b %d")+
   #scale_colour_manual(values=c("dodgerblue2"), name="",breaks=c("r1"))+
   scale_fill_manual(values=c("dodgerblue1","dodgerblue2"), name="",breaks=c("r1","r2"))+
   labs(x = "Days", y=paste("Cumulative death cases"))+
   coord_cartesian(ylim = c(1,NA),expand = FALSE)+
   theme_classic()+
   theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
         axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
         legend.position = c(0.2,0.8),legend.key = element_rect(colour = "transparent", fill = "white"),
         legend.background = element_rect(fill = "transparent",colour =NA),
         axis.line =element_line(color = "black", 
                                 size = 1, linetype = "solid"),
         axis.ticks = element_line(size = 1))+
   guides(colour= FALSE,fill=FALSE)
library(Rmisc)
library(lattice)
multiplot(p1,p2, cols =2)

