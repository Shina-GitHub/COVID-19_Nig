# editted on July 15 2020
# #libraries needed
library(tidyverse)
library(dplyr)

#import data
covid_nig_dat <- read_csv("Covid_Data_Nig_May7.csv")
str(covid_nig_dat)
covid_nig_dat$Date =as.Date(covid_nig_dat$Date,"%d/%m/%Y")
#covert date to integer for all cases
distinct_date = distinct(covid_nig_dat, Date)
#my_date = as.Date(distinct_date$Date,"%d/%m/%Y")
int_date = as.numeric(distinct_date$Date)
int_date = sort(int_date)
int_date_count <- c(int_date[1],diff(int_date))
int_date_count[1] = 0
#Distinct date with correspinding serial number
data_uniq_date = data.frame(Tday =cumsum(int_date_count), my_date =sort(distinct_date$Date))

##get local cases
local_dat <- covid_nig_dat %>% filter(Sources == "Local")
##get import cases
import_dat <- covid_nig_dat %>% filter(Sources == "Imported")


## find day integer for corressponding date
day_count_upd<- function(data_withoutinteger, data_withinteger){
  len_var <- length(data_withoutinteger$Date)
  result = rep(0,len_var)
  for (i in seq(1,len_var)){
    result[i] = data_withinteger$Tday[data_withoutinteger$Date[i]== data_withinteger$my_date]
  }
  return(result)
}

local_Tday <- day_count_upd(local_dat,data_uniq_date)
import_Tday<- day_count_upd(import_dat,data_uniq_date)
local_dat$Tday <-local_Tday 
import_dat$Tday<-import_Tday
##daily report is missing and we adjust for it by
max_day <- max(local_Tday)
T_seq = seq(0,max_day)
###Generate daily  date
daily_date = data_uniq_date$my_date[1]
daily_date = daily_date + T_seq 
#################################

daily_report<-function(Tseq,mysource){
  #mysource is daily report by Sources (local vs imported)
  daily_type = rep(0,length(Tseq))
  for (i in seq(1,length(Tseq))){
    
    if(sum(Tseq[i] == mysource$Tday)>=1){ #check if the serial number is present
      daily_type[i] = mysource$Count[Tseq[i] == mysource$Tday]
    }
  }
  return(list(tday = Tseq, cases = daily_type))
}

mydaily_local <-daily_report(T_seq,local_dat)
mydaily_import <-daily_report(T_seq,import_dat)

###################################

###My functions
#predict_local
predict_local <- function(Rt_est,T_seq, gen_time,loc_case,all_count){
  #ij is counter the determines number of likelihood needs
  loc_pred = rep(0,length(T_seq)-1)
  rate_loc = rep(0,length(T_seq)-1)
  for (ik in seq(1,length(T_seq)-1)){
    imp_cond = gen_time[1:ik]*all_count[ik:1]
    rate_loc[ik] = Rt_est[ik]*(sum(imp_cond))
    #print(length(T_seq))
    #print(rate[ik])
    loc_pred[ik] = round(rate_loc[ik])
  }
  return(loc_pred)
}
#likelihood
lik <- function(params,T_seq, gen_time,loc_case,all_count, log= FALSE){
  #ij is counter the determines number of likelihood needs
  like = rep(0,length(T_seq))
  rate = rep(0,length(T_seq))
  for (ik in seq(1,(length(T_seq)-1))){
    imp_cond = gen_time[1:ik]*all_count[ik:1]
    rate[ik] = params[ik]*(sum(imp_cond))
    #print(sum(imp_cond))
    #print(length(T_seq))
    #print(rate[ik])
    # if (rate[ik]== 0){like[ik] = -2^(31)}
    # else{like[ik] = dpois(loc_case[ik], rate[ik], log = TRUE)}
    like[ik] = dpois(loc_case[ik+1], rate[ik], log = TRUE)
  }
  #print(gen_time[1:length(T_seq)])
  #print(like)
  loglik <- sum(like)
  
  
  if  (log){ return (loglik)}
  else {return (exp(loglik))} 
}
#MCMC
g.all <- function (params) {
  # #######################
  sig = 0.005
  which <- sample(1:length(params),1)
  proposed <- params
  proposed[which] <- -1
  while (proposed[which]<0){  #No negative values
    proposed[which] <- rnorm(1, mean = params[which], sd = sig*(1+loc_case[which]) )
  }
  return (list(proposed = proposed, ratio = 1))
  # #############
  # #newValue <- -1
  # #while (newValue < 0){
  #   
  # #Edit to ensure single iteration and corrected for negative valaues
  #      proposed<- mvrnorm(1, mu = params, Sigma = sig)
  #      proposed = ifelse(proposed < 0, 0, proposed)
  #      #test <- proposed >= 0  # Ensuring positive values is selected
  #      #if (sum(test) > length(params)-1){ newValue = 1}
  # #}
  # 
  # return (list(proposed = proposed, ratio = 1))
}
#Prior
dprior_exp <-function (params, log = FALSE){
  logprior <-  sum(dunif(params, min = 0, max = 5, log = TRUE))
  # dgamma(params["Rt"], shape=7.5, scale = 1, log = TRUE)
  if  (log){ return (logprior)}
  else {return (exp(logprior))}
}
####Reff
#Initialise R
Rt_est = matrix(0,nrow = (length(T_seq)-1), ncol = 1)
Rt_est_lower = matrix(0,nrow=(length(T_seq)-1), ncol = 1)
Rt_est_upper =matrix(0,nrow=(length(T_seq)-1), ncol = 1)
loc_case = mydaily_local$cases
all_count = mydaily_local$cases + round(40*mydaily_import$cases/100)
##Examining the data visually to allow the possibility of Infection process
T_seq_new = T_seq
T_seq = seq(0, (max(T_seq_new)-20))
loc_case = loc_case[21:71]
all_count = all_count[21:71]
daily_date_new = daily_date
daily_date = daily_date[21:71]
#Re-Initialise R
Rt_est = matrix(0,nrow = (length(T_seq)-1), ncol = 1)
Rt_est_lower = matrix(0,nrow=(length(T_seq)-1), ncol = 1)
Rt_est_upper =matrix(0,nrow=(length(T_seq)-1), ncol = 1)

  #define serial interval
  gen_time = dlnorm(T_seq[-1], meanlog = 1.39, sdlog = 0.57, log = FALSE)
  # plot(T_seq,gen_time, type = "l", lwd = 3, xlab = "Days",
  #      ylab = "serial interval")
  #
  #initial values
  #
  #Update initial Reff
  #Update initial Reff
  for(iv in seq(1,(length(T_seq)-1))){
    if (loc_case[iv+1]>0){Rt_est[iv]=1}
  }
  
  #define params
  params= Rt_est 
  
  #number of iterations
  mc <-1000000
  
  # storage for output
  mcmc.out <- matrix (NA, mc, length(params)+1)
  #the main loop
  for (i in 1:mc){
    prop <- g.all(params)
    new <-prop$proposed
    #print(new)
    #print(lik(new, data_baye, log = TRUE))
    
    log.accept.ratio <- ( dprior_exp(new, log =TRUE) - dprior_exp(params, log =TRUE)+
                            lik(new,T_seq, gen_time,loc_case,all_count, log= TRUE)-
                            lik(params,T_seq, gen_time,loc_case,all_count, log= TRUE)-
                            log(prop$ratio)
    )
    accept.ratio <- min(1,exp(log.accept.ratio))
    #print(accept.ratio)
    #print(accept.ratio)
    if (runif(1)< accept.ratio){
      params <- new
      
    }
    mcmc.out[i, ] <- c(params, lik(params,T_seq, gen_time,loc_case,all_count, log= TRUE)+
                         dprior_exp(params, log =TRUE))
  }
  ##burn out 10%
  burnin <- seq(100000, 1000000, by = 1)
  accept <- 1-mean(duplicated(mcmc.out[burnin,]))
  
  #first para
  data_median = mcmc.out[burnin,-(length(T_seq))] #No plus one because length of Rt is Length (T_seq)-1
  para1<- apply(data_median,2,median)
  bootCI1<-apply(data_median,2, quantile,c(0.025,0.25,0.75,0.975))
  boot_adj = matrix(bootCI1,nrow=4)
  Rt_est = para1
  Rt_est_lower5p=boot_adj[1,]
  Rt_est_lower25p=boot_adj[2,]
  Rt_est_upper75p=boot_adj[3,]
  Rt_est_upper95p=boot_adj[4,]
  
  ###Calculate DIC 
  likBIC <- lik(Rt_est,T_seq, gen_time,loc_case,all_count, log= TRUE)
  BIC_me <- length(T_seq)*length(Rt_est)-2*likBIC  
  
  library(pracma)
  
  Rt_est = movavg(Rt_est, 4, "s")
  Rt_est_lower5p= movavg(Rt_est_lower5p, 4, "s")
  Rt_est_lower25p= movavg(Rt_est_lower25p, 4, "s")
  Rt_est_upper75p = movavg(Rt_est_upper75p, 4, "s")
  Rt_est_upper95p = movavg(Rt_est_upper95p, 4, "s")
  
  #Save Reff
  # write.csv(Rt_est,"Reff.csv_April26.csv")
  # 
  # plot(T_seq,Rt_est,type="l",ylim=c(0,10))
  # lines(T_seq,Rt_est_lower,type="l",col = "red")
  # lines(T_seq,Rt_est_upper,type="l",col = "red")
  # abline(h = 1, col="red", lwd=3, lty=2)
  # predict cases
  pred<-predict_local(Rt_est,T_seq, gen_time,loc_case,all_count)
  pred_lower<-predict_local(Rt_est_lower5p,T_seq, gen_time,loc_case,all_count)
  pred_upper<-predict_local(Rt_est_upper95p,T_seq, gen_time,loc_case,all_count)
  # Nice Ploting
  #
  data_plot1 = data.frame(days = daily_date, import = (all_count-loc_case), local_cases = loc_case)
  data_plot2 = data.frame(days =  daily_date[-1], R_eff = Rt_est,
                          R_lower5p =Rt_est_lower5p, R_lower25p =Rt_est_lower25p,
                          R_upper75p =Rt_est_upper75p,R_upper95p =Rt_est_upper95p,
                          pred_local = pred,pred_local_lower = pred_lower,
                          pred_local_upper = pred_upper)
  
  

  ggplot(data_plot2)+
    geom_line(data=data_plot2, aes(x=days,y=R_eff),color = "blue",
              linetype= 1,size = 2)+
    geom_ribbon(data = data_plot2,aes(ymin=R_lower5p,ymax=R_upper95p,
                                      x=days,fill= "r1"),inherit.aes = FALSE,alpha=0.3)+
    geom_ribbon(data = data_plot2,aes(ymin=R_lower25p,ymax=R_upper75p,
                                      x=days,fill= "r2"),inherit.aes = FALSE,alpha=0.3)+
    geom_hline(yintercept=1, linetype="dashed", color = "red",size=2)+
    scale_x_date(date_breaks = "2 week", date_labels = "%b %d")+
    #scale_colour_manual(values=c("dodgerblue2"), name="",breaks=c("r1"))+
    scale_fill_manual(values=c("dodgerblue1","dodgerblue2"), name="",breaks=c("r1","r2"))+
    labs(x = "Days", y=paste("R(t)"))+
    coord_cartesian( ylim = c(0,5), expand = TRUE)+
    theme_classic()+
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
          legend.position = c(0.2,0.8),legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent",colour =NA),
          axis.line =element_line(color = "black", 
                                  size = 1, linetype = "solid"),
          axis.ticks = element_line(size = 1))+
    guides(colour= FALSE,fill=FALSE)

# p2<-ggplot(data= data_plot1,aes(x=days,y=import, fill = 'r1'))+
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_bar(data_plot1,mapping=aes(x=days,y=local_cases, fill = 'r2'),
#            stat="identity", color="black", position=position_dodge())+
#   geom_line(data=data_plot2, aes(x=days,y=pred_local, color="r3"),
#             linetype= 2,size = 2)+
#   geom_ribbon(data = data_plot2,aes(ymin=pred_local_lower,ymax=pred_local_upper,
#                                    x=days,fill= "r4"),show.legend = FALSE,inherit.aes = FALSE,alpha=0.3)+
#   scale_x_date(date_breaks = "1 week", date_labels = "%b %d")+
#   scale_colour_manual(values=c("black"), name="",breaks=c("r3"),
#                       label = c(""))+
#   scale_fill_manual(values=c("deeppink","deepskyblue","black"),
#                     label = c("importation","local cases",  "Model fit"),
#                     name="COVID-19",breaks= c("r1","r2","r4"))+
#   labs(x = "Days", y=paste("SARS-CoV-2 incidence", paste0("(","NSW",")")))+
#   coord_cartesian( expand = TRUE)+
#   theme_bw()+
#   theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
#         legend.position = c(0.2,0.8),legend.key = element_rect(colour = "transparent", fill = "white"),
#         legend.background = element_rect(fill = "transparent",colour =NA))+
#   guides(color = FALSE)
# library(Rmisc)
# library(lattice)
# multiplot(p1,p2, cols =2)