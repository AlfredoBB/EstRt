library(TTR)
library(dplyr)
library(tidyr)
#library(R6)
library(EpiEstim)
library(incidence)
library(MASS) #fitdistr
library(openxlsx)
library(ISOcodes)
library(ggplot2)
library(truncnorm)
library(parallel)
library(readr) #read_csv
offline<-TRUE

#debug<-TRUE 

mcmc_length<-200

dbfs<-"Distribute backwards before smoothing"

jhu_str <-"CSSE John Hopkins University" # All OK 09/05/2021
jhu_source <-"JHU CSSE COVID-19 Dataset, https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data"
jhu_path<-"https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"


world_pop<-read.csv("WPP2019_TotalPopulationBySex.csv") %>% dplyr::filter(Variant=="Medium",Time==2020)
world_pop<-rbind(world_pop,world_pop[world_pop$Location=="United States of America",])
world_pop[nrow(world_pop),"Location"]<-"US"
world_pop<-rbind(world_pop,world_pop[world_pop$Location=="Russian Federation",])
world_pop[nrow(world_pop),"Location"]<-"Russia"
world_pop<-rbind(world_pop,world_pop[world_pop$Location=="Iran (Islamic Republic of)",])
world_pop[nrow(world_pop),"Location"]<-"Iran"
world_pop$PopTotal <-world_pop$PopTotal*1000


 
fdf_to_inci<-function(dataframe,dateformat){ # Filtered dataframe to incidence
  #accepts a dataframe passed like %>% rename (Dates=,Cum=) 
  return (dataframe %>% tidyr::drop_na()  %>% dplyr::mutate(dates=as.Date(dates,dateformat)) %>%
            dplyr::arrange(dates) %>% dplyr::mutate(I=Cum - lag (Cum,order_by=dates)) %>%
            dplyr::select (dates,I) %>% dplyr::slice (-1))# %>% tibble::column_to_rownames(var="dates"))
  #incidf<-data.frame(dates=incidf[2:nrow(incidf),"FECHA"],I=diff(incidf$Cum))
}


#################################################### SI DISTRIBUTIONS ####################################
DuXuWu<-read.xlsx("https://raw.githubusercontent.com/MeyersLabUTexas/COVID-19/master/Table%20S5.xlsx",startRow=2);
China_si_data<-c(DuXuWu["Seconday.-.symptom.onset.date"]-DuXuWu["Index.-.symptom.onset.date"])[[1]]
Du_et_al_distr<-fitdistr(China_si_data[China_si_data>0],"lognormal")


Zhao_mean<-4.9
Zhao_sd<-4.4

Zhao_meanlog<-log(Zhao_mean^2/sqrt(Zhao_sd^2+Zhao_mean^2))
Zhao_sdlog<-sqrt(log((Zhao_sd^2/Zhao_mean^2+1)))
#sd(rlnorm(10000,Zhao_meanlog,Zhao_sdlog))
#Zhao_mean<-(4.75)^2/3.96
#Zhao_shape<-3.96/Zhao_scale


#Zhao_scale<-(4.75)^2/3.96
#Zhao_shape<-3.96/Zhao_scale

#zhao_means<-rnorm(1000,4.75,0.43/2)
#zhao_sds_<-rnorm(1000,4.75)
#Du_et_al_distr_mom<-c(mean=exp(Du_et_al_distr$estimate["meanlog"]+(Du_et_al_distr$estimate["meanlog"]^2)/2),
#                      sd= 

get_last_R<-function(dataframe,mcmc=NA,meanlog=SI_fit_clever_X@ests["meanlog","est"],sdlog=SI_fit_clever_X@ests["sdlog","est"],config){
  if (!(is.na(mcmc)))
  return(estimate_R(dataframe,method="si_from_sample",si_sample=mcmc,config=config)$R[nrow(dataframe)-9,`Mean.R.`])
  else {
    return(estimate_R(dataframe,method="non_parametric_si",config=config)$R[nrow(dataframe)-9,`Mean.R.`])
    }
}

Nishiura_si<-function(){ #license non-specified
  Nishiura<-read.csv("https://raw.githubusercontent.com/aakhmetz/nCoVSerialInterval2020/master/data/supplemetary%20table.csv",stringsAsFactors = TRUE) %>% dplyr::mutate(type=0) #%>% dplyr::rename(ER=EL,EL=ER)
  SI_fit_clever_Nishiura<<-coarseDataTools::dic.fit.mcmc(dat=Nishiura, dist= "L", init.pars = init_mcmc_params(Nishiura, "L"), burnin = 2000, n.samples=10000)
  check_cdt_samples_convergence(SI_fit_clever_Nishiura@samples)
  si_sample_mcmc_N<-coarse2estim(SI_fit_clever_Nishiura,thin=5,dist="L")
  save(si_sample_mcmc_N,SI_fit_clever_Nishiura,file="C:\\Users\\HP\\OneDrive\\Documentos\\myR\\si_sample_mcmc_N_T2k.R")
}


Du_si<-function(){
  dat<-data.frame(EL=0,ER=0,SL=DuXuWu["Seconday.-.symptom.onset.date"]-DuXuWu["Index.-.symptom.onset.date"],SR=DuXuWu["Seconday.-.symptom.onset.date"]-DuXuWu["Index.-.symptom.onset.date"],type=2)
  colnames(dat)<-c("EL","ER","SL","SR","type")
  SI_fit_clever_X<-coarseDataTools::dic.fit.mcmc(dat=dat[dat$SR>0,], dist= "L", init.pars = init_mcmc_params(dat[dat$SR>0,], "L"), burnin = 2000, n.samples=10000)
  check_cdt_samples_convergence(SI_fit_clever_X@samples)
  #SI_fit_naive<-coarseDataTools::dic.fit.mcmc(dat=dat[dat$SR>0,], dist= "L", burnin = 10000, n.samples=50000)
  #check_cdt_samples_convergence(SI_fit_naive@samples)
  si_sample_mcmc_X<-coarse2estim(SI_fit_clever_X, thin = 5,dist="L")
  save(si_sample_mcmc_X,SI_fit_clever_X,file="C:\\Users\\HP\\OneDrive\\Documentos\\myR\\si_sample_mcmc_X_T2k.R")
}

Zhao_si<-function(){
  Zhao<-read.xlsx("https://assets.researchsquare.com/files/rs-18805/v3/dataset.xlsx") %>% mutate(mintime=min(as.Date(Infector.date.lwr,"%m/%d/%Y"))) %>%
    mutate (EL=as.numeric(as.Date(Infector.date.lwr,"%m/%d/%Y")-mintime), ER=as.numeric(as.Date(Infector.date.upr,"%m/%d/%Y")-mintime),
            SL=as.numeric(as.Date(Infectee.date,"%m/%d/%Y")-mintime), SR=as.numeric(as.Date(Infectee.date,"%m/%d/%Y")-mintime),type=2 - (EL!=ER)) %>%
    dplyr::select(c(EL,ER,SL,SR,type))
  SI_fit_clever_Zhao<<-coarseDataTools::dic.fit.mcmc(dat=Zhao, dist= "L", init.pars = init_mcmc_params(Zhao, "L"), burnin = 2000, n.samples=10000)
  check_cdt_samples_convergence(SI_fit_clever_Zhao@samples)
  si_sample_mcmc_Zhao<-coarse2estim(SI_fit_clever_Zhao, thin = 5,dist="L")
  save(si_sample_mcmc_Zhao,SI_fit_clever_Zhao,file="C:\\Users\\HP\\OneDrive\\Documentos\\myR\\si_sample_mcmc_Zhao_T2k.R")
}

Ali_si<-function(){ #MIT License 15/10/2020
  Ali<-read.csv("https://raw.githubusercontent.com/PDGLin/COVID19_EffSerialInterval_NPI/master/raw_data/TableS1_1407TransPairs.csv",encoding="UTF-8") %>%
    mutate(infectee_onsetDate=as.Date(infectee_onsetDate), infector_onsetDate=as.Date(infector_onsetDate), infector_returnDate_fromOtherCity = as.Date(infector_returnDate_fromOtherCity),infector_isolateDate_beforeSymptom=as.Date(infector_isolateDate_beforeSymptom),infector_isolateDate_afterSymptom=as.Date(infector_isolateDate_afterSymptom)) %>%
    mutate(mintime=min(infector_onsetDate, na.rm=TRUE)) %>%
    mutate(EL = as.numeric(infector_onsetDate - mintime), ER = as.numeric(infector_onsetDate - mintime), SL = as.numeric(infectee_onsetDate - mintime), SR = as.numeric(infectee_onsetDate - mintime),type=2) %>%
    dplyr::filter((!is.na(EL)),(!is.na(SR)),SL>ER) %>% dplyr::select(c(EL,ER,SL,SR,type))

  #as.numeric(max(infector_onsetDate,min(infector_isolateDate_beforeSymptom,infector_isolateDate_afterSymptom,na.rm=TRUE),na.rm=TRUE) - mintime)
  SI_fit_clever_Ali<<-coarseDataTools::dic.fit.mcmc(dat=Ali, dist= "L", init.pars = init_mcmc_params(Ali, "L"), burnin = 2000, n.samples=10000)
  check_cdt_samples_convergence(SI_fit_clever_Ali@samples)
  #SI_fit_naive<-coarseDataTools::dic.fit.mcmc(dat=Ali, dist= "L", burnin = 10000, n.samples=50000)
  #check_cdt_samples_convergence(SI_fit_naive@samples)
  si_sample_mcmc_Ali<-coarse2estim(SI_fit_clever_Ali, thin = 5,dist="L")
  save(si_sample_mcmc_Ali,SI_fit_clever_Ali,file="C:\\Users\\HP\\OneDrive\\Documentos\\myR\\si_sample_mcmc_Ali_T2k.R")
}

#Nishiura_si()
#Du_si()
#Zhao_si()
#Ali_si()



gamma_m_to_par<-function(mean,sd){
  theta<-sd^2/mean
  kappa<-mean/theta
  return(list(theta=theta,kappa=kappa))
}

#mode_sing_mean<-5.20*3-3.78-6.78
#mode_sing_sd<-1.72*3-0.91-3.93
#mode_tian_mean<-3.95*3-3.01-4.91
#mode_tian_sd<-1.51*3-0.74-2.97

sing_m<-gamma_m_to_par(mean=5.2,sd=1.72)
tian_m<-gamma_m_to_par(mean=3.95,sd=1.51)
liu_m<-gamma_m_to_par(mean=5.81,sd=3.24)
li_m<-gamma_m_to_par(mean=7.5,sd=3.4)
ali_early_isolation<-list(mean=3.3,mean_ci_low=2.7,mean_ci_high=3.8,sd=4.5,sd_ci_low=4.1,sd_ci_high=4.9)
ali_late_isolation<-list(mean=6.8,mean_ci_low=6.2,mean_ci_high=7.3,sd=5.3,sd_ci_low=4.9,sd_ci_high=5.7)


G_sing_distr<-matrix(0,ncol=mcmc_length,nrow=31)
G_tian_distr<-matrix(0,ncol=mcmc_length,nrow=31)
Ali_late_distr<-matrix(0,ncol=mcmc_length,nrow=31)
Ali_early_distr<-matrix(0,ncol=mcmc_length,nrow=31)

sing_distros<-gamma_m_to_par(mean=mc2d::rpert(ncol(G_sing_distr),min=3.78,max=6.78,mode=5.20,shape=4),
                             sd=mc2d::rpert(ncol(G_sing_distr),min=0.91,max=3.93,mode=1.72,shape=4))

tian_distros<-gamma_m_to_par(mean=mc2d::rpert(ncol(G_tian_distr),min=3.01,max=4.91,mode=3.95,shape=4),
                             sd=mc2d::rpert(ncol(G_tian_distr),min=0.74,max=2.97,mode=1.51,shape=4))
ali_early_distros<-list(mean=mc2d::rpert(ncol(Ali_early_distr),min=ali_early_isolation$mean_ci_low,max=ali_early_isolation$mean_ci_high,
                               mode=ali_early_isolation$mean,shape=4),
                        sd=mc2d::rpert(ncol(Ali_early_distr),min=ali_early_isolation$sd_ci_low,max=ali_early_isolation$sd_ci_high,
                                       mode=ali_early_isolation$sd,shape=4))
ali_late_distros<-list(mean=mc2d::rpert(ncol(Ali_late_distr),min=ali_late_isolation$mean_ci_low,max=ali_late_isolation$mean_ci_high,
                                        mode=ali_late_isolation$mean,shape=4),
                       sd=mc2d::rpert(ncol(Ali_late_distr),min=ali_late_isolation$sd_ci_low,max=ali_late_isolation$sd_ci_high,
                                      mode=ali_late_isolation$sd,shape=4))

for (i in 1:ncol(G_sing_distr)){ 
  distro<-dgamma(seq(0,30),shape=sing_distros$kappa[i],#kappa
                 scale=sing_distros$theta[i])#theta
  G_sing_distr[,i]<-distro/sum(distro)
}
for (i in 1:ncol(G_tian_distr)){ 
  distro<-dgamma(seq(0,30),shape=tian_distros$kappa[i],#kappa
                 scale=tian_distros$theta[i])#theta
  G_tian_distr[,i]<-distro/sum(distro)
}
for (i in 1:ncol(Ali_early_distr)){ 
  distro<-dtruncnorm(seq(0,30),mean=ali_early_distros$mean[i],#mean
                 sd=ali_early_distros$sd[i],a=1e-5)#sd
  Ali_early_distr[,i]<-distro/sum(distro)
}
#Ali_early_distr[2,]<-Ali_early_distr[2,] + Ali_early_distr[1,]
#Ali_early_distr[1,]<-0

for (i in 1:ncol(Ali_late_distr)){ 
  distro<-dtruncnorm(seq(0,30),mean=ali_late_distros$mean[i],#mean
                sd=ali_late_distros$sd[i],a=1e-5)#sd
  Ali_late_distr[,i]<-distro/sum(distro)
}
#Ali_late_distr[2,]<-Ali_late_distr[2,] + Ali_late_distr[1,]
#Ali_late_distr[1,]<-0


sing_distros<-gamma_m_to_par(mean=mc2d::dpert(ncol(G_sing_distr),min=3.78,max=6.78,mode=5.20,shape=4),
                             sd=mc2d::dpert(ncol(G_sing_distr),min=0.91,max=3.93,mode=1.72,shape=4))

tian_distros<-gamma_m_to_par(mean=mc2d::dpert(ncol(G_tian_distr),min=3.01,max=4.91,mode=3.95,shape=4),
                             sd=mc2d::dpert(ncol(G_tian_distr),min=0.74,max=2.97,mode=1.51,shape=4))

si_distr_boot<-matrix(0,ncol=600,nrow=31)
for (i in 1:ncol(si_distr_boot)){ 
  boot_fit<-fitdistr(sample(China_si_data[China_si_data>0],replace=TRUE),"Lognormal")
  #if(i==1) plot(density(rlnorm(500,boot_fit$estimate)))
  #lines(density(rlnorm(1000,meanlog=rnorm(1,boot_fit$estimate["meanlog"],boot_fit$sd["meanlog"]),
  #                     sdlog=rnorm(1,boot_fit$estimate["sdlog"],boot_fit$sd["sdlog"]))))
  distro<-dlnorm(seq(0,30),meanlog=rnorm(1,boot_fit$estimate["meanlog"],boot_fit$sd["meanlog"]),
                 sdlog=rnorm(1,boot_fit$estimate["sdlog"],boot_fit$sd["sdlog"]))
  si_distr_boot[,i]<-distro/sum(distro)
}
 
##############################################

narr_str<- "Narrativa Coronavirus Data API" # All OK 09/05/2021
narr_source<-"Narrativa Coronavirus Data API, https://covid19tracking.narrativa.com/index_en.html"

ecdc_str <- "European CDC"

casal_str <- "Spain (R Casal)"

uk_str<-"UK (T White)"
uk_source<-"T White, UK Covid 19 data, https://github.com/tomwhite/covid-19-uk-data, UK government"

#alvaro_str<-("Spain (@alvariteus)")
#alvaro_source<-"Twitter: @alvariteus, Spain data, https://docs.google.com/spreadsheets/d/1aOSfDXbqawSCngJwPHUiRDQEA2QLFUyZ7XJ6P1nl6LI/edit#gid=230974174. See the original for specific sources."

MSC_path<-"https://cnecovid.isciii.es/covid19/resources/agregados.csv"
MSC_str<-"Spain cumulative: (MSC, deprecated)"

CAT_z_path<-("https://analisi.transparenciacatalunya.cat/api/views/xuwf-dxjd/rows.csv?accessType=DOWNLOAD&sorting=true")
CAT_m_path<-("https://analisi.transparenciacatalunya.cat/api/views/jj6z-iyrp/rows.csv?accessType=DOWNLOAD&sorting=true")

CAT_zc_str<-("Catalonia Official - Basic health areas - Conf (PCR)")
CAT_zcr_str<-("Catalonia Official - Basic health areas - Conf (PCR + rapid test)")
CAT_zs_str<-("Catalonia Official - Basic health areas - Suspected")
CAT_zsc_str<-("Catalonia Official - Basic health areas - Susp and/or PCR+")
CAT_zrap_str<-("Catalonia Official - Basic health areas - Rapid test")
CAT_zabs_str<-("Catalonia Official - Basic health areas - Antibodies")
CAT_mc_str<-("Catalonia Official - Town/Municipio - Confirmed (PCR)")
CAT_mcr_str<-("Catalonia Official - Town/Municipio - Confirmed (PCR + rapid test)")
CAT_ms_str<-("Catalonia Official - Town/Municipio - Suspected")
CAT_msc_str<-("Catalonia Official - Town/Municipio - Suspected and/or PCR+")
CAT_mrap_str<-("Catalonia Official - Town/Municipio - Rapid test")
CAT_mabs_str<-("Catalonia Official - Town/Municipio - Antibodies")

CAT_z_source<-("Registre de casos de COVID-19 realitzats a Catalunya. Segregació per sexe i àrea bàsica de salut (ABS). Publisher: dades obertes catalunya. Author: Departament de Salut. Generalitat de Catalunya. https://analisi.transparenciacatalunya.cat/Salut/Registre-de-casos-de-COVID-19-realitzats-a-Catalun/xuwf-dxjd")
CAT_m_source<-("Registre de casos de COVID-19 realitzats a Catalunya. Segregació per sexe i municipi. Publisher: dades obertes catalunya. Author: Departament de Salut. Generalitat de Catalunya. https://analisi.transparenciacatalunya.cat/Salut/Registre-de-casos-de-COVID-19-realitzats-a-Catalun/jj6z-iyrp")

USA_path<-("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
USAs_path<-"https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"
USA_str<-("USA by State and County (The New York Times)")
USA_source<-("Coronavirus (Covid-19) Data in the United States. Data from The New York Times, based on reports from state and local health agencies. https://github.com/nytimes/covid-19-data")

setClass("CD",slots=list(path="character",name="character",CountryName="character",df="data.frame",source="character",provinces="character",states="character",getData="function",
                         getStates="function",getProvinces="function",includeAll="logical",date_obtained="Date",source_type="character"))
setGeneric("getData", function(x) standardGeneric("getData") )
setMethod("getData","CD",function(x){df<-read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE);
  if(x@provinces!="X"){df[,x@provinces]<- NAtoAll(df[,x@provinces])};
  if(length(x@CountryName)!=0){df["Country"]<-x@CountryName};
  if(x@states!="X"){df[,x@states]<- NAtoAll(df[,x@states])};
return(df)})
setGeneric("getStates", function(x) standardGeneric("getStates") )
setMethod("getStates","CD",function(x){vector_of_states<-c(if (x@includeAll) "* All *" else NULL, (if (is.factor(x@df[,x@states]))
  as.character(levels(x@df[,x@states])) else unique(x@df[,x@states])));
if (is.list(vector_of_states)) list_of_states <- list(vector_of_states %>%  purrr::flatten_chr())
else list_of_states <- list(vector_of_states)
names(list_of_states)<-x@states
return(list_of_states)
})
setGeneric("getProvinces", function(x,statefilter=NA) standardGeneric("getProvinces") )
setMethod("getProvinces","CD",function(x,statefilter=NA){
  p_list<-x@df[if (length(statefilter==1) & (is.na(statefilter[1]) | statefilter[1]=="* All *")) TRUE else (as.character(unlist(x@df[,x@states])) %in% statefilter),x@provinces]#  %>% purrr::flatten_chr()
  vector_of_provinces<-c(if (x@includeAll) "* All *" else NULL,if(is.factor(p_list)) as.character(unique(p_list)) else unique(p_list))
  if (is.list(vector_of_provinces)) list_of_provinces <- list(vector_of_provinces %>%  purrr::flatten_chr())
  else list_of_provinces <- list(vector_of_provinces)
  names(list_of_provinces)<-x@provinces;
  return(list_of_provinces)
})
setGeneric("getDate", function(x) standardGeneric("getDate") )
setMethod("getDate","CD",function(x){return(Sys.Date())})
datasources<-list(
  #test with
  #datasources[["Brazil_cities"]]@df<-datasources[["Brazil_cities"]]@getData(datasources[["Brazil_cities"]])
  Italy= new("CD",path="https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-province/dpc-covid19-ita-province.csv", # OK 08/05/2021
           name="Italy: official data",
           source="Dati COVID-19 Italia. Dipartimento della Protezione Civile - Emergenza Coronavirus: la risposta nazionale",
           states="denominazione_regione",
           provinces="denominazione_provincia",
           getData=getMethod("getData","CD"),
           getStates=getMethod("getStates","CD"),
           getProvinces=getMethod("getProvinces","CD"),
           date_obtained=getMethod("getDate","CD")(),
           includeAll=FALSE,
           source_type="Country_official"
  ), Dhub= new("CD",path="https://datahub.io/core/covid-19/r/time-series-19-covid-combined.csv", # OK 08/05/2021
            name="DataHub COVID-19 time series",
            source="DataHub Novel Coronavirus 2019 https://datahub.io/core/covid-19",
            states="Country.Region",
            provinces="Province.State",
            getData=getMethod("getData","CD"),
            getStates=getMethod("getStates","CD"),
            getProvinces=getMethod("getProvinces","CD"),
            date_obtained=getMethod("getDate","CD")(),
            includeAll=FALSE,
            source_type="World"
  ), IND= new("CD",path="https://www.isibang.ac.in/~athreya/incovid19/data/summarymohfw1update.csv", # OK 08/05/2021
           name="India: Athreya compilation",
           source="COVID-19 India-Timeline an understanding across States and Union Territories",
           states="Country.Region",
           provinces="X",
           getData=function(x){
             df<-read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE);
             colnames(df)<-paste(colnames(df),unlist(df[1,]))
             df<-(df[-1,] %>% tidyr::pivot_longer(colnames(df[,-1]),names_to="dates_",values_to="Cum"))
             df$mes<-sapply(substr(df$dates_,5,7),function(x) {return(which(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")==x))})
             df$dates<-as.Date(paste("2020-",df$mes,"-",substr(df$dates_,2,3),sep=""),"%Y-%m-%d")
             df<-df[sapply(df$dates_,function(x){return(strsplit(x," ")[[1]][2])}) %in% c("TCIN","TCFN"),] %>% group_by(`X `,dates) %>% summarize (Cum=sum(as.numeric(Cum))) %>%
               rename(`State`=`X `)
             return(df)
           },
           getStates=function(x){return(c("India"))},
           getProvinces=function(x,country){return(list(states=as.character(levels(x@df$State))))},
           date_obtained=getMethod("getDate","CD")(),
           includeAll=FALSE,
           source_type="Country_semiofficial"
  ),MSC_cum= new("CD",path="agregados.csv",
                name=MSC_str,
                source="Spain cumulative incidence data (previously official, now deprecated). https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos",
                states="Country",
                provinces="CCAA",
                getData=function(x){return(subset(read.csv(x@path,encoding="UTF-8",na.strings="NANA",stringsAsFactors = TRUE),(FECHA!="") & (!(is.na(`PCR.`)))) %>%
                                             mutate(Country="Spain")
                                           )},
                getStates=function(x){return(c("Spain"))},
                getProvinces=function(x,country){return(list(provinces=c("* All *",as.character(unique(x@df[,x@provinces])))))},
                date_obtained=getMethod("getDate","CD")(),
                includeAll=FALSE,
                source_type="Country_official"
  ), Spain_provinces= new("CD",path="https://cnecovid.isciii.es/covid19/resources/casos_tecnica_provincia.csv", #OK 2021-05-08
                name="Spain RENAVE (epi. surv.) provinces",
                source="Spain provinces COVID-19 data from the Red Nacional de Vigilancia Epidemiológica (Surveillance). Instituto de Salud Carlos III (ISCIII). Ministerio de Sanidad, Consumo y Bienestar Social. Gobierno de España. https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos   https://datos.gob.es/es/catalogo/e05070101-evolucion-de-enfermedad-por-el-coronavirus-covid-19",
                states="Country",
                provinces="provincia_names",
                getData=function(x){
                  return(read.csv(x@path,encoding="UTF-8",na.strings="NANA") %>%
                           dplyr::mutate(Country="Spain",provincia_iso=paste("ES-",provincia_iso,sep="")) %>%
                           merge(ISO_3166_2 %>% dplyr::select(Code,Name),by.x="provincia_iso",by.y="Code",all.x=TRUE,all.y=FALSE) %>%
                           dplyr::mutate(provincia_names=ifelse(is.na(Name),provincia_iso,Name)))
                         
                },
                getStates=function(x){return(c("Spain"))},
                getProvinces=function(x,country){return(list(province_iso=c("* All *",as.character(unique(x@df[,x@provinces])))))},
                date_obtained=getMethod("getDate","CD")(),
                includeAll=FALSE,
                source_type="Country_official"
  ), Spain_regions= new("CD",path="https://cnecovid.isciii.es/covid19/resources/casos_tecnica_ccaa.csv", #OK 2021-05-08
                name="Spain RENAVE (epi. surv.) regions",
                source="Spain regions COVID-19 data from the Red Nacional de Vigilancia Epidemiológica (Surveillance). Instituto de Salud Carlos III (ISCIII). Ministerio de Sanidad, Consumo y Bienestar Social. Gobierno de España. https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos   https://datos.gob.es/es/catalogo/e05070101-evolucion-de-enfermedad-por-el-coronavirus-covid-19",
                states="Country",
                provinces="ccaa_names",
                getData=function(x){
                  ISO_3166_2[ISO_3166_2$Code=="ES-NC","Name"]<-"Navarra, Comunidad Foral de"
                  ISO_3166_2[ISO_3166_2$Code=="ES-PV","Name"]<-"País Vasco"
                  ISO_3166_2[ISO_3166_2$Code=="ES-VC","Name"]<-"Valenciana, Comunidad"
                  ISO_3166_2[ISO_3166_2$Code=="ES-CT","Name"]<-"Cataluña"
                  return(read.csv(x@path,encoding="UTF-8",na.strings="NANA",stringsAsFactors = TRUE) %>%
                           mutate(Country="Spain",ccaa_iso=paste("ES-",ccaa_iso,sep="")) %>%
                           merge(ISO_3166_2 %>% dplyr::select(Code,Name),by.x="ccaa_iso",by.y="Code",all.x=TRUE,all.y=FALSE) %>%
                           dplyr::mutate(ccaa_names=ifelse(is.na(Name),ccaa_iso,Name)));
                },
                getStates=function(x){return(c("Spain"))},
                getProvinces=function(x,country){return(list(ccaa_names=c("* All *",as.character(unique(x@df[,x@provinces])))))},
                date_obtained=getMethod("getDate","CD")(),
                includeAll=TRUE,
                source_type="Country_official"
  ), BEL = new("CD",path="https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv", #OK 2021-05-08
            name="Belgium (Sciensano), by province",
            source="Belgium official (Sciensano), by province. https://epistat.wiv-isp.be/covid/",
            states="REGION",
            provinces="PROVINCE",
            getData=function(x){
              return(read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE) %>% dplyr::filter(!is.na(DATE) & !is.na(REGION)));
            },
            getStates=getMethod("getStates","CD"),
            getProvinces=getMethod("getProvinces","CD"),
            date_obtained=getMethod("getDate","CD")(),
            includeAll=TRUE,
            source_type="Country_official"
  ), PT = new("CD",path="https://raw.githubusercontent.com/dssg-pt/covid19pt-data/master/data_concelhos.csv", #OK all 2021-05-08
           name="Portugal (DSSG-pt)",
           source="Portugal (Data Science for Social Good), by concelho. GPL-3.0 license. https://github.com/dssg-pt/covid19pt-data",
           states="Country",
           provinces="Concelho",
           getData=function(x){
             return(read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE) %>% pivot_longer(-1,names_to="Concelho",values_to="Cum")
                    %>% rename (`dates`=data) %>% mutate(Cum=replace_na(Cum,0),Country="Portugal"))
           },
           getStates=function(...){return("Portugal")},
           getProvinces=getMethod("getProvinces","CD"),
           date_obtained=getMethod("getDate","CD")(),
           includeAll=FALSE,
           source_type="Country_semiofficial"
  ), Germany = new("CD",path="https://raw.githubusercontent.com/jgehrcke/covid-19-germany-gae/master/cases-rki-by-state.csv", #OK all 2021-05-08
            name="Germany (RKI, J Gehrcke)",
            source="Germany (Robert Koch Institute, compiled by J Gehrcke), by state. MIT license. https://github.com/jgehrcke/covid-19-germany-gae/",
            states="Country",
            provinces="State",
            getData=function(x){
              ger<-read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE);
              return(ger %>% rename (`Total`=sum_cases,`dates`=time_iso8601) %>% pivot_longer(cols=-1,names_to="State",values_to="Cum") %>% mutate(Country="Ger"));
            },
            getStates=function(...){return("Ger")},
            getProvinces=getMethod("getProvinces","CD"),
            date_obtained=getMethod("getDate","CD")(),
            includeAll=FALSE,
            source_type="Country_semiofficial"
  ), GerAgs = new("CD",path="https://raw.githubusercontent.com/jgehrcke/covid-19-germany-gae/master/cases-rki-by-ags.csv", #OK 2021-05-08
               name="Germany by Ags (RKI, J Gehrcke)",
               source="Germany (Robert Koch Institute, compiled by J Gehrcke), by state AND AGS. MIT license. https://github.com/jgehrcke/covid-19-germany-gae/",
               states="state",
               provinces="ags_name",
               getData=function(x){
                 agsjs<-jsonlite::read_json("https://raw.githubusercontent.com/jgehrcke/covid-19-germany-gae/master/ags.json", simplifyVector=TRUE)
                 ags<-as.data.frame(do.call(rbind,agsjs)) %>% tibble::rownames_to_column("ags") %>% mutate(ags=paste("X",ags,sep=""))
                 ger<-read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE);
                 return(ger %>% dplyr::select(-sum_cases) %>% rename (`dates`=time_iso8601) %>% pivot_longer(cols=-1,names_to="ags",values_to="Cum") %>%
                          merge(ags,by="ags") %>% dplyr::rename(ags_name=name));
               },
               getStates=getMethod("getStates","CD"),
               getProvinces=getMethod("getProvinces","CD"),
               date_obtained=getMethod("getDate","CD")(),
               includeAll=TRUE,
               source_type="Country_semiofficial"
  ), Canada = new("CD",path="https://health-infobase.canada.ca/src/data/covidLive/covid19.csv", #OK 2021-05-08
              name="Canada, by province/territory",
              source="Public Health Infobase - Data on COVID-19 in Canada. https://open.canada.ca/data/en/dataset/261c32ab-4cfd-4f81-9dea-7b64065690dc. License: Open Government Licence - Canada",
              CountryName="Canada",
              states="Country",
              provinces="prname",
              getData=getMethod("getData","CD"),
              #return(read.csv(x@path,encoding="UTF-8") %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
              getStates=getMethod("getStates","CD"),
              getProvinces=getMethod("getProvinces","CD"),
              date_obtained=getMethod("getDate","CD")(),
              includeAll=FALSE,
              source_type="Country_official"
  ), Chile = new("CD",path="https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19_std.csv", #OK 2021-05-08
              name="Chile MinCiencia, by region & comuna",
              source="Chile official COVID-19 epidemiological reports, by comuna: DP1 - Casos totales por comuna incremental. (https://github.com/MinCiencia/Datos-COVID19/tree/master/output/producto1). Licences: http://www.minciencia.gob.cl/sites/default/files/1771596.pdf",
              states="Region",
              provinces="Comuna",
              getData=getMethod("getData","CD"),
              #return(read.csv(x@path,encoding="UTF-8") %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
              getStates=getMethod("getStates","CD"),
              getProvinces=getMethod("getProvinces","CD"),
              date_obtained=getMethod("getDate","CD")(),
              includeAll=FALSE,
              source_type="Country_official"
  ), Chile_Regions = new("CD",path="https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto11/bulk/producto4.csv", #OK 2021-05-08
                 name="Chile MinCiencia, by region",
                 source="Chile official COVID-19 MINSAL epidemiological reports, by region: DP11 - Enriquecimiento del Data Product 4. (https://github.com/MinCiencia/Datos-COVID19/tree/master/output/producto11). License: MIT. http://www.minciencia.gob.cl/sites/default/files/1771596.pdf",
                 states="Region",
                 provinces="X",
                 getData=getMethod("getData","CD"),
                 #return(read.csv(x@path,encoding="UTF-8") %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
                 getStates=getMethod("getStates","CD"),
                 getProvinces=function(...){return(c("* All *"))},
                 date_obtained=getMethod("getDate","CD")(),
                 includeAll=FALSE,
                 source_type="Country_official"
  ), Mex= new("CD",path="https://raw.githubusercontent.com/carranco-sga/Mexico-COVID-19/master/Mexico_COVID19_CTD.csv", #OK 2021-05-08
             name="Mexico (Secretaría de Salud Federal, compiled)",
             source="Mexico (Secretaría de Salud Federal, compiled by carranco-sga). (https://github.com/carranco-sga/Mexico-COVID-19/).",
             states="State",
             provinces="X",
             getData=function(x){
               mex<-read.csv(x@path,encoding="UTF-8",stringsAsFactors = TRUE);
               return((mex[,sapply(strsplit(colnames(mex),"_"),function(x) length(x)<2)]) %>% dplyr::select(-c("Susp","Recovered","Deceased")) %>% rename (`Total`=Pos) %>% #"Pos_L","Pos_rep","Susp_rep","Neg_rep","IRAG_Test","Tested_tot",
                       pivot_longer(cols=-1,names_to="State",values_to="Cum"));
             },
             #return(read.csv(x@path,encoding="UTF-8") %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
             getStates=getMethod("getStates","CD"),
             getProvinces=function(...){return(c("* All *"))},
             date_obtained=getMethod("getDate","CD")(),
             includeAll=FALSE,
             source_type="Country_semiofficial"
  ), Brazil_cities= new("CD",path="https://raw.githubusercontent.com/wcota/covid19br/master/cases-brazil-cities-time.csv.gz", #OK all 2021-05-08
                name="Brazil cities (W Cota)",
                source="Brazil, cases by city. (Wesley Cota), (https://covid19br.wcota.me/).",
                states="state",
                provinces="city",
                getData=function(x){
                  return(read_csv(file=x@path, progress=FALSE, col_types=cols(),guess_max=10^6))},# %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
                getStates=getMethod("getStates","CD"),
                getProvinces=getMethod("getProvinces","CD"),
                date_obtained=getMethod("getDate","CD")(),
                includeAll=FALSE,
                source_type="Country_semiofficial"
  ), Brazil_states= new("CD",path="https://raw.githubusercontent.com/wcota/covid19br/master/cases-brazil-states.csv", #OK all 2021-05-08
                 name="Brazil states (W Cota)",
                 source="Brazil (cases by state). (Wesley Cota), (https://covid19br.wcota.me/).",
                 states="country",
                 provinces="state",
                 getData=getMethod("getData","CD"),
                 getStates=getMethod("getStates","CD"),
                 getProvinces=getMethod("getProvinces","CD"),
                 date_obtained=getMethod("getDate","CD")(),
                 includeAll=FALSE,
                 source_type="Country_semiofficial"
  ),Peru= new("CD",path="https://raw.githubusercontent.com/jmcastagnetto/covid-19-peru-data/main/datos/covid-19-peru-data.csv", #OK 2021-05-08
              name="Peru (Minsal, JM Castagnetto)",
              source="Peru (Minsal, JM Castagnetto. (https://github.com/jmcastagnetto/covid-19-peru-data).",
              states="country",
              provinces="region",
              getData=getMethod("getData","CD"),
              getStates=function(...){return("Peru")},
              getProvinces=getMethod("getProvinces","CD"),
              date_obtained=getMethod("getDate","CD")(),
              includeAll=FALSE,
              source_type="Country_semiofficial"
  ), Col= new("CD",path="https://www.datos.gov.co/api/views/gt2j-8ykr/rows.csv?accessType=DOWNLOAD", #OK 2021-05-08
             name="Colombia Official",
             source="Colombia. Casos positivos de COVID-19 en Colombia. Instituto Nacional de la Salud. https://www.datos.gov.co/Salud-y-Protecci-n-Social/Casos-positivos-de-COVID-19-en-Colombia/gt2j-8ykr",
             states="Departamento.o.Distrito",
             provinces="Ciudad.de.ubicación",
             getData=function(x){
               return(read.csv(x@path,encoding="UTF-8") %>% group_by(Fecha.diagnostico,Ciudad.de.ubicación,Departamento.o.Distrito) %>% summarize(I=length(ID.de.caso)))
             },
             getStates=getMethod("getStates","CD"),
             getProvinces=getMethod("getProvinces","CD"),
             date_obtained=getMethod("getDate","CD")(),
             includeAll=TRUE,
             source_type="Country_official"
  ))


deactivated_datasources=list(Madrid_ZBS= new("CD",path="http://datos.comunidad.madrid/catalogo/dataset/b3d55e40-8263-4c0b-827d-2bb23b5e7bab/resource/43708c23-2b77-48fd-9986-fa97691a2d59/download/covid19_tia_zonas_basicas_salud_s.csv", #OK 2021-05-08
                       name="Madrid (ZBS) Official",
                       source="Covid 19 -TIA Zonas Básicas de Salud. Madrid. http://datos.comunidad.madrid/catalogo/dataset/covid19_tia_zonas_basicas_salud. Obtained via vertical union of data reported daily and weekly.",
                       states="Region",
                       provinces="zona_basica_salud",
                       getData=function(x){
                         return(rbind(read.csv2(x@path,encoding="UTF-8",fileEncoding="latin1")[,c("zona_basica_salud","fecha_informe","casos_confirmados_totales")] %>%
                                        dplyr::mutate (fecha_informe=as.Date(fecha_informe,"%Y/%m/%d %H:%M:%S")) %>% dplyr::filter(fecha_informe>=as.Date("2020-07-01")),
                                      read.csv2("http://datos.comunidad.madrid/catalogo/dataset/b3d55e40-8263-4c0b-827d-2bb23b5e7bab/resource/b7b9edb4-0c70-47d3-9c64-8c4913830a24/download/covid19_tia_zonas_basicas_salud.csv",
                                               encoding="UTF-8",fileEncoding="latin1")[,c("zona_basica_salud","fecha_informe","casos_confirmados_totales")]  %>%
                                        dplyr::mutate (fecha_informe=as.Date(fecha_informe,"%Y/%m/%d %H:%M:%S")) %>% dplyr::filter(fecha_informe<as.Date("2020-07-01"))
                                      
                                      
                                      ) %>% mutate (Region="Madrid"))
                         
                       },
                       #return(read.csv(x@path,encoding="UTF-8") %>% dplyr::pivot_longer(cols=-c("Region","Codigo.region","Comuna","Codigo.comuna","Poblacion")));
                       getStates=getMethod("getStates","CD"),
                       getProvinces=getMethod("getProvinces","CD"),
                       date_obtained=getMethod("getDate","CD")(),
                       includeAll=TRUE,
                       source_type="Region_official"
  )
)

#datasources[["Madrid_ZBS"]]@df<-datasources[["Madrid_ZBS"]]@getData(datasources[["Madrid_ZBS"]])

class_sources<-lapply(datasources,function(x) return(x@name))
Country_official_sources<-lapply(datasources,function(x) {if (x@source_type=="Country_official") return(x@name)})
Regional_official_sources<-lapply(datasources,function(x) {if (x@source_type=="Region_official") return(x@name)})
Country_semiofficial_sources<-lapply(datasources,function(x) {if (x@source_type=="Country_semiofficial") return(x@name)})
World_sources<-lapply(datasources,function(x) {if (x@source_type=="World") return(x@name)})
sources_str<-append(list(jh=jhu_str,ecdc=ecdc_str,narr=narr_str,#casal=casal_str,alv=alvaro_str, 
                         uk=uk_str,catmc=CAT_mc_str,catmsc=CAT_msc_str,catzc=CAT_zc_str,catzsc=CAT_zsc_str,usa=USA_str),class_sources)

sources_list<-list(World_data=c(jhu_str,ecdc_str,narr_str,datasources[["Dhub"]]@name),
                   Official_country_data=unlist(c(Country_official_sources,use.names=FALSE)),
                   Semiofficial_country_data=unlist(c(uk_str,USA_str,Country_semiofficial_sources,use.names=FALSE)),
                   Official_regional_data=unlist(c(Regional_official_sources,CAT_mc_str,CAT_msc_str,CAT_zc_str,CAT_zsc_str,use.names=FALSE))
)


#---------------
Ln_bs_du_str<-"SI: Lognormal + Bootstrap (Du et al) (much slower)"
Ln_mc_du_str<-"SI: Lognormal + MCMC (Du et al) (much slower)"
Ln_du_str<-"SI: Lognormal (Du et al)"
Ln_mc_ni_str<-"SI: Lognormal + MCMC (Nishiura et al)"
Ln_ni_str<-"SI: Lognormal (Nishiura et al)"
Ln_mc_zh_str<-"SI: Lognormal + MCMC (Zhao et al)"
Ln_zh_str<-"SI: Lognormal (Zhao et al)"
Ln_mc_ali_str<-"SI: Lognormal + MCMC (Ali et al)"
Ln_ali_str<-"SI: Lognormal (Ali et al)"
Tn_bs_ali_early_str<-"SI: Trunc. normal + Bootstrap (Ali et al, early isolation) (much slower)"
Tn_ali_early_str<-"SI: Trunc. normal (Ali et al, early isolation)"
Tn_bs_ali_late_str<-"SI: Trunc. normal + Bootstrap (Ali et al, late isolation) (much slower)"
Tn_ali_late_str<-"SI: Trunc. normal (Ali et al, late isolation)"

G_sing_bs_str<-"GT: Gamma (Singapore) + BS (Ganyani et al)"
G_tian_bs_str<-"GT: Gamma (Tianjin) + BS (Ganyani et al)"
G_sing_str<-"GT: Gamma (Singapore) (Ganyani et al)"
G_tian_str<-"GT: Gamma (Tianjin) (Ganyani et al)"


App_si_equivalences<-list(      `Du Point`=Ln_du_str,
                                `Du Bs`=Ln_bs_du_str,
                                Du=Ln_mc_du_str,
                                Nishiura = Ln_mc_ni_str,
                                `Nishiura Pt`= Ln_ni_str,
                                Zhao = Ln_mc_zh_str,
                                `Zhao Pt`= Ln_zh_str,
                                Ali = Ln_mc_ali_str,
                                `Ali Pt` = Ln_ali_str,
                                `Ali early` = Tn_bs_ali_early_str,
                                `Ali early Pt` = Tn_ali_early_str,
                                `Ali late` = Tn_bs_ali_late_str,
                                `Ali late Pt` = Tn_ali_late_str,
                                `Ganyani Singapore` = G_sing_bs_str,
                                `Ganyani Singapore Pt` = G_sing_str,
                                `Ganyani Tianjin` = G_tian_bs_str,
                                `Ganyani Singapore Pt` = G_tian_str)


si_strs<-unlist(App_si_equivalences,use.names=FALSE)
#si_strs<-c(Ln_bs_du_str, Ln_mc_du_str, Ln_du_str, Ln_mc_ni_str, Ln_ni_str,Ln_mc_zh_str,Ln_zh_str,G_sing_bs_str,G_sing_str,G_tian_bs_str,G_tian_str)

NAtoAll<-function(vector){
  if(is.factor(vector)){
    vector<-factor(vector,levels=c("* All *",as.character(levels(vector))))
    vector[vector==""]<-"* All *"
    vector[is.na(vector)]<-"* All *"
  }
  return(vector)
}

#--------------------------------------------------


get_cols<-function(df,source,country="* All *",province = "* All *",warnings=""){
  if (country!=""){
    if (source==jhu_str) {
      if (province!=""){
        incidf<-fdf_to_inci(df[(df$Country.Region==country) & (df$Province.State==province),] %>%
                              dplyr::select(-c("Province.State","Country.Region","Lat","Long")) %>% t() %>% as.data.frame() %>%
                              dplyr::rename (Cum=1) %>% tibble::rownames_to_column("dates"),"X%m.%d.%y")
      }
    }
    else if (source==ecdc_str){
      incidf<-df[df$countriesAndTerritories==country,c("dateRep","cases")] %>%
        dplyr::rename (dates=dateRep,I=cases)  %>% dplyr::mutate(dates=as.Date(dates,"%d/%m/%Y")) %>% dplyr::arrange(dates)
    }
    else if (source==narr_str){ # All OK 09/05/2021
      incidf<-fdf_to_inci(df[(df$Country_EN==country) & (df$Region==province),] %>%
                            dplyr::select(-c("Country_EN","Country_ES","Country_IT","Region")) %>% t() %>% as.data.frame() %>%
                            dplyr::rename (Cum=1) %>% tibble::rownames_to_column("dates"),"X%Y.%m.%d")
    }
    else if (source==casal_str){
      incidf<-df[df$ccaa==province,c("fecha","nuevos")] %>%
        dplyr::rename(dates=fecha,I=nuevos) %>% dplyr::mutate(dates=as.Date(dates))
    }
    else if (source==alvaro_str){
      incidf<-fdf_to_inci(data.frame(t(df[province,])) %>% dplyr::rename(Cum=1) %>%
                            tibble::rownames_to_column("dates"),"X%d.%m.%Y")
    }
    else if (source %in% c("MSC_cum",MSC_str,datasources[["MSC_cum"]]@name)){
      incidf<-fdf_to_inci(df[if (province %in% c("* All *","")) TRUE else (df$CCAA==province),c("PCR.","FECHA")] %>%
                            dplyr::group_by(FECHA) %>% dplyr::summarize(Cum=sum(PCR.)) %>% dplyr::rename(dates="FECHA"), "%d/%m/%Y")
    }
    else if (source==uk_str){
      incidf<-fdf_to_inci(df[(df$Country==country) & (df$Area==province),c("Date","TotalCases")] %>%
                            rename(dates="Date",Cum="TotalCases"),"%Y-%m-%d")}
    else if (source==USA_str){
      incidf<-fdf_to_inci(df[(if (country %in% c("* All *","")) TRUE else df$state==country) &
                               (if (province %in% c("* All *","")) TRUE else(df$county==province)),]  %>%
                            rename(dates="date") %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(cases)),"%Y-%m-%d") 
    } 
    else if (source %in% c(CAT_zc_str,CAT_zs_str,CAT_zrap_str,CAT_zabs_str,CAT_zsc_str,CAT_zcr_str)){ # All OK 09/05/2021
      incidf<-df[(if (country %in% c("* All *","")) TRUE else (df$RegioSanitariaDescripcio==country)) &
                   (if (province %in% c("* All *","")) TRUE else(df$ABSDescripcio==province)) &
                   (if(source==CAT_zc_str) ((df$TipusCasDescripcio=="Positiu PCR") || (df$TipusCasDescripcio=="Positiu TAR")) else
                     (if(source==CAT_zs_str) ((as.character(df$TipusCasDescripcio) == "Sospitós") ||(as.character(df$TipusCasDescripcio) == "Epidemiològic") ||(as.character(df$TipusCasDescripcio) == "PCR probable")) else
                       (if(source==CAT_zrap_str) (df$TipusCasDescripcio =="Positiu per Test Ràpid") else
                         (if(source==CAT_zabs_str) (df$TipusCasDescripcio =="Positiu per ELISA") else
                           (if(source==CAT_zsc_str) (as.character(df$TipusCasDescripcio) %in% c("Sospitós","Positiu PCR")) else
                             (if(source==CAT_zcr_str) (as.character(df$TipusCasDescripcio) %in% c("Positiu PCR","Positiu per Test Ràpid")) else FALSE)))))),
                 c("TipusCasData","NumCasos")] %>%
        dplyr::group_by(TipusCasData) %>% dplyr::summarize(I= sum(NumCasos,na.rm=TRUE)) %>% 
        dplyr::mutate(dates=as.Date(TipusCasData,"%d/%m/%Y")) %>% tidyr::drop_na() %>%
        dplyr::select(c("dates","I"))# %>% tibble::column_to_rownames(var="dates")
    }
    else if (source %in% c(CAT_mc_str,CAT_ms_str,CAT_mrap_str,CAT_mabs_str,CAT_msc_str,CAT_mcr_str)){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else (df$ComarcaDescripcio==country)) &
                   (if (province %in% c("* All *","")) TRUE else(df$MunicipiDescripcio==province)) &
                   (if(source==CAT_mc_str) (df$TipusCasDescripcio=="Positiu PCR") else
                     (if(source==CAT_ms_str) (df$TipusCasDescripcio =="Sospitós") else
                       (if(source==CAT_mrap_str) (df$TipusCasDescripcio ==("Positiu per Test Ràpid")) else
                         (if(source==CAT_mabs_str) (df$TipusCasDescripcio ==("Positiu per ELISA")) else
                           (if(source==CAT_msc_str) (as.character(df$TipusCasDescripcio) %in% c("Sospitós","Positiu PCR")) else
                             (if(source==CAT_mcr_str) (as.character(df$TipusCasDescripcio) %in% c("Positiu PCR","Positiu per Test Ràpid")) else FALSE)))))),
                 c("TipusCasData","NumCasos")] %>%
        dplyr::group_by(TipusCasData) %>% dplyr::summarize(I= sum(NumCasos,na.rm=TRUE)) %>% 
        dplyr::mutate(dates=as.Date(TipusCasData,"%d/%m/%Y")) %>% tidyr::drop_na() %>% dplyr::arrange(dates) %>%
        dplyr::select(c("dates","I"))# %>% tibble::column_to_rownames(var="dates")
    }
    else if (source %in% c("Italy",datasources[["Italy"]]@name)){
      incidf<-fdf_to_inci(df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Italy"]]@states]==country) &
                               (if (province %in% c("* All *","")) TRUE else(df[,datasources[["Italy"]]@provinces]==province)),]  %>%
                            rename(dates="data") %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(totale_casi)),"%Y-%m-%dT%T") 
    }
    else if (source %in% datasources[["Dhub"]]@name) {
      incidf<-fdf_to_inci(df[(df$Country.Region==country) & 
                               (if (province %in% c("* All *","")) TRUE else df[,datasources[["Dhub"]]@provinces]==province),] %>%
                            rename(dates="Date") %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(Confirmed)),"%Y-%m-%d")
    }
    else if (source %in% datasources[["IND"]]@name){
      incidf<-fdf_to_inci(df[(df$State==province),c("dates","Cum")],"%Y-%m-%d")
    }
    else if (source %in% datasources[["Spain_provinces"]]@name){
      incidf<-df[if (province=="* All *") TRUE else (df[,datasources[["Spain_provinces"]]@provinces]==province),c("fecha","num_casos")] %>% dplyr::group_by (fecha) %>%
        dplyr::summarize(I=sum(num_casos)) %>% dplyr::rename(dates=fecha) %>% dplyr::mutate(dates=as.Date(dates,"%Y-%m-%d")) %>% dplyr::arrange(dates)
    }
    else if (source %in% c("Spain_regions",datasources[["Spain_regions"]]@name)){
      incidf<-df[if (province=="* All *") TRUE else (df[,datasources[["Spain_regions"]]@provinces]==province),c("fecha","num_casos")] %>% dplyr::group_by (fecha) %>%
        dplyr::summarize(I=sum(num_casos)) %>% dplyr::rename(dates=fecha) %>% dplyr::mutate(dates=as.Date(dates,"%Y-%m-%d")) %>% dplyr::arrange(dates)
    }
    else if (source %in% datasources[["BEL"]]@name){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["BEL"]]@states]==country) &
                   (if (province %in% c("* All *","")) TRUE else(df[,datasources[["BEL"]]@provinces]==province)),]  %>%
        mutate(dates=as.Date(DATE,"%Y-%m-%d")) %>% dplyr::group_by(dates) %>% dplyr::summarize(I= sum(CASES)) %>%
        tidyr::drop_na()  %>% dplyr::arrange(dates) %>% dplyr::select (dates,I)
    }
    else if (source %in% c("PT",datasources[["PT"]]@name)){
      incidf<-fdf_to_inci(df[df[,datasources[["PT"]]@provinces]==province,],"%d-%m-%Y")
    }
    else if (source %in% c("Germany",datasources[["Germany"]]@name)){
      incidf<-fdf_to_inci(df[df[,datasources[["Germany"]]@provinces]==province,],"%Y-%m-%dT%H:%M:%S%z")
    }
    else if (source %in% c("GerAgs",datasources[["GerAgs"]]@name)){
      incidf<-fdf_to_inci(df[(if (country %in% c("* All *","")) TRUE else (df[,datasources[["GerAgs"]]@states]==country)) &
                               (if (province %in% c("* All *","")) TRUE else(df[,datasources[["GerAgs"]]@provinces]==province)),] %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum=sum(Cum)),"%Y-%m-%dT%H:%M:%S%z")
    }
    else if (source %in% c("Canada",datasources[["Canada"]]@name)){
      incidf<-fdf_to_inci(df[(df[,datasources[["Canada"]]@provinces]==province),]  %>%
                            dplyr::rename(dates=date) %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(numconf)),"%d-%m-%Y")
    }
    else if (source %in% c("Chile",c(datasources[["Chile"]]@name))){
      incidf<-fdf_to_inci(df[(if (country %in% c("* All *","")) (!(province %in% c("Total"))) else df[,datasources[["Chile"]]@states]==country) &
                               (if (province %in% c("* All *","")) (!(province %in% c("Total"))) else(df[,datasources[["Chile"]]@provinces]==province)),]  %>%
                            dplyr::rename(dates=Fecha) %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(Casos.confirmados)),"%Y-%m-%d")
    }
    else if (source %in% c("Chile_Regions",datasources[["Chile_Regions"]]@name)){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Chile_Regions"]]@states]==country),]  %>%
        mutate(dates=as.Date(Fecha,"%Y/%m/%d")) %>% dplyr::group_by(dates) %>% dplyr::summarize(I= sum(Casosnuevostotales)) %>%
        tidyr::drop_na(c(dates,I))  %>% dplyr::arrange(dates) %>% dplyr::select (dates,I)
    }
    else if (source %in% c("Mexico",datasources[["Mex"]]@name)){
      incidf<-fdf_to_inci(df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Mex"]]@states]==country),]  %>%
                            dplyr::mutate(dates=Fecha) %>%
                            dplyr::group_by(dates) %>% dplyr::summarize(Cum= sum(Cum)) %>% 
                            dplyr::select (c(dates,Cum)),"%Y-%m-%d")
    }
    else if (source %in% datasources[["Brazil_states"]]@name){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Brazil_states"]]@states]==country) &
                   (if (province %in% c("* All *","")) TRUE else(df[,datasources[["Brazil_states"]]@provinces]==province)),]  %>%
        mutate(dates=as.Date(date,"%Y-%m-%d")) %>% dplyr::group_by(dates) %>% dplyr::summarize(I= sum(newCases)) %>%
        tidyr::drop_na()  %>% dplyr::arrange(dates) %>% dplyr::select (dates,I)
    }
    else if (source %in% datasources[["Brazil_cities"]]@name){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Brazil_cities"]]@states]==country) &
                   (if (province %in% c("* All *","")) TRUE else(df[,datasources[["Brazil_cities"]]@provinces]==province)),]  %>%
        mutate(dates=as.Date(date,"%Y-%m-%d")) %>% dplyr::group_by(dates) %>% dplyr::summarize(I= sum(newCases)) %>%
        tidyr::drop_na()  %>% dplyr::arrange(dates) %>% dplyr::select (dates,I)
    }
    else if (source %in% c("Peru",datasources[["Peru"]]@name)){
      incidf<-fdf_to_inci(df[df[,datasources[["Peru"]]@provinces]==province,]  %>%
                            rename(dates=date) %>% tidyr::drop_na(confirmed) %>% dplyr::group_by(dates) %>% dplyr::summarize(Cum=sum(confirmed,na.rm=TRUE)),"%Y-%m-%d")
    }
    else if (source %in% c("Col",datasources[["Col"]]@name)){
      incidf<-df[(if (country %in% c("* All *","")) TRUE else df[,datasources[["Col"]]@states]==country) &
                   (if (province %in% c("* All *","")) TRUE else(df[,datasources[["Col"]]@provinces]==province)),]  %>%
        mutate(dates=as.Date(Fecha.diagnostico,"%Y-%m-%d")) %>% dplyr::group_by(dates) %>% dplyr::summarize(I= sum(I)) %>%
        tidyr::drop_na()  %>% dplyr::arrange(dates) %>% dplyr::select (dates,I)
    }
    else if (source %in% c("Madrid_ZBS",datasources[["Madrid_ZBS"]]@name)){
      incidf<-fdf_to_inci(df[(if (province %in% c("* All *","")) TRUE else df[,datasources[["Madrid_ZBS"]]@provinces]==province),]  %>%
                            rename(dates=fecha_informe) %>% dplyr::group_by(dates)  %>% dplyr::summarize(Cum=sum(casos_confirmados_totales,na.rm=TRUE)) %>%
        tidyr::drop_na(Cum),"%Y/%m/%d %H:%M:%S")
    }
  }
  return(incidf)
}



incidf_prepro<-function(incidf,undetected=0,negatives=dbfs,positives="Do nothing",smooth_method,smooth_value,ignore=3,warnings=""){
  # Preprocessing
  if (!("dates" %in% colnames(incidf))){
    if ("date" %in% colnames(incidf))
      incidf <- incidf %>% dplyr::rename(dates=date)
  }
  if(is.factor(incidf$dates)) incidf$dates <- as.Date(incidf$dates)

  if (nrow(incidf)==0) {return (list(df=data.frame(incidf)))}
  incidf<-incidf[order(incidf$dates),] %>% tidyr::drop_na()
  incidf$I<-incidf$I/(1-undetected)
  incidf$I_pre_smooth<-incidf$I
  
  lsmooth<-NA
  
  if (max(incidf$dates)-min(incidf$dates)>=nrow(incidf)){
    #print(paste("non-contiguous dates: ",max(incidf$dates)-min(incidf$dates),"days vs",nrow(incidf),"rows"))
    warnings<-paste(warnings,"Non-contiguous dates: ",max(incidf$dates)-min(incidf$dates)+1,"days in the interval but",nrow(incidf),"days are available in the data; the data has been completed with zeros.")
    zeros<-data.frame(dates=seq(as.Date(min(incidf$dates)),as.Date(max(incidf$dates)),1),I=0,I_pre_smooth=0)
    zeros<-zeros[!(zeros$dates %in% incidf$dates),]
    incidf<-rbind(incidf,zeros)
    incidf<-incidf[order(incidf$dates),]
  }
  
  if (negatives!="Do nothing"){
    if (negatives=="Remove before smoothing") {incidf<-incidf[incidf$I>=0,]}
    else if (negatives=="to 0 before smoothing"){incidf[incidf$I<0,"I"]<-0}
    else if (negatives==dbfs){
      for (i_bis in 1:4){
        if(any(incidf$I<0)){
          last_distribution<- 1
          for (i in (1:nrow(incidf))) {
            if ((incidf[i,"I"]<0)){
              #print(paste("i",i,"=",incidf[i,"I"]))
              j<-1
              while(((i+j)<nrow(incidf)) & (incidf[i+j,"I"]<0)) {
                #print(paste("j",j))
                incidf[i,"I"]<-(incidf[i,"I"]+incidf[i+j,"I"])
                incidf[i+j,"I"]<-0
                j<-j+1
              }
              
              #print(incidf[(last_distribution:(i-1)),"I"])
              #print(incidf[(last_distribution:(i-1)),"I"] *
              #        (sum(incidf[(last_distribution:i),"I"]) / sum(incidf[(last_distribution:(i-1)),"I"])))
              incidf[(last_distribution:(i-1)),"I"]<-incidf[(last_distribution:(i-1)),"I"] *
                (sum(incidf[(last_distribution:i),"I"]) / sum(incidf[(last_distribution:(i-1)),"I"]))
              incidf[i,"I"]<-0
              last_distribution<-i
              #if(sum(incidf$I_pre_smooth)!=sum(incidf$I)){print("warning: total sum was changed after backwards distribution")}
            }
          }
        }
      }
    }
  }
  if (positives=="Distribute backwards before smoothing") # incidf<-get_cols(msc_csv,source=MSC_str,province="CT")
  {
    last_distribution<-1
    for (i in (3:(nrow(incidf)-3))) {
      if ((incidf[i,"I"]>200) & (all((incidf[(i-2):(i+2),"I"])>75)) & all((incidf[i,"I"]/c(incidf[c((i-2),(i-1),(i+1),(i+2)),"I"]))>3))# & (incidf[i,"I"]/incidf[i-1,"I"]>3) & (incidf[i,"I"]/incidf[i+1,"I"]>3) & (incidf[i,"I"]/incidf[i+2,"I"]>3))
      {
        distribute_this<-incidf[i,"I"]-mean(unlist(incidf[(i-2):(i+2),"I"]))
        incidf[i,"I"]<-round(mean(unlist(incidf[(i-2):(i+2),"I"])))
        incidf[(last_distribution:(i-1)),"I"]<-round(incidf[(last_distribution:(i-1)),"I"]*
                                                       (sum(incidf[(last_distribution:(i-1)),"I"],distribute_this) / sum(incidf[(last_distribution:(i-1)),"I"])))
        last_distribution<-i
      }
    }
  }
  if (smooth_value>0){
    if (smooth_value>=nrow(incidf)) smooth_value=nrow(incidf)-1
    #print("smoothing")
    smooth_past<-ifelse(substr(smooth_method,5,11)=="(past)",TRUE,FALSE)
    if (smooth_method=="SMA (centered)"){
      incidf$I<-as.numeric(forecast::ma(x=incidf$I,order=smooth_value,centre=TRUE))                    
    }
    else if (substr(smooth_method,1,3)=="SMA"){
      lsmooth<-TTR::SMA(x=incidf[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1),"I"], n=smooth_value)
      incidf$I<-lsmooth[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1)]
    }
    if (substr(smooth_method,1,3)=="EMA"){
      lsmooth<-TTR::EMA(x=incidf[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1),"I"], n=smooth_value)
      incidf$I<-lsmooth[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1)]
    }
    if (substr(smooth_method,1,3)=="WMA"){
      lsmooth<-TTR::WMA(x=incidf[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1),"I"], n=smooth_value)
      incidf$I<-lsmooth[if (smooth_past) c(1:nrow(incidf)) else c(nrow(incidf):1)]    
    }                 
    else if (smooth_method=="Lowess"){
      lsmooth<-lowess(x=incidf$dates, y=incidf$I, f=smooth_value/nrow(incidf))
      incidf$I<-lsmooth$y
      incidf[incidf$I<0,"I"]<-0
    }
    else if (smooth_method=="Loess"){
      lsmooth<-loess(data=incidf,formula=I~dates, span=smooth_value/nrow(incidf))
      incidf$I<-lsmooth$y
    }
    incidf<-incidf[!is.na(incidf$I),] # MA smooth creates NAs.
  if (negatives!="Do nothing"){
  }
    if (negatives=="Remove after smoothing") {incidf<-incidf[incidf$I>=0,]}
    else if (negatives=="to 0 after smoothing"){incidf[incidf$I<0,"I"]<-0}
    #else if (negatives=="Extrapolate"){ for (i in input) incidf[incidf$I<0,"I"]<-0}  }
  }
  incidf$I<-round(incidf$I)
  
  while(nrow(incidf)>2 & sum(incidf[1:min(nrow(incidf),3),"I"])<ignore){incidf<-incidf[2:nrow(incidf),]}
  return (list(df=data.frame(incidf),smooth=lsmooth,warnings=warnings))
  
}


#names(inci)<-gsub(".","-",names(inci),fixed=TRUE)
#names(inci)<-gsub("X","",names(inci),fixed=TRUE)
#names(inci)<-substr(paste("2020-0",names(inci),sep=""),1,10)
#names(inci)<-sapply(names(inci),FUN =function(x){
#    if(substr(x,str_length(x),str_length(x))=="-")
#    {paste(substr(x,1,str_length(x)-2),"0",substr(x,str_length(x)-1,str_length(x)-1),sep="")
#    }else{x}
#})

# ------------------------------------------------------------------------------------
#clusterExport(cl=cl,ls())
load("si_sample_mcmc_N_Tk.R")
load("si_sample_mcmc_X_Tk.R")
load("si_sample_mcmc_Zhao_Tk.R")
load("si_sample_mcmc_Ali_Tk.R")




calculate_R<-function(incidf=Spain,prepro=TRUE,undetected=0.1,negatives="Distribute backwards before smoothing",smooth=5,smooth_method="SMA (centered)",distribution="Du",window=7,ignore=15,return_object=FALSE,mean_prior=5,std_prior=5){
  if (prepro){incidf<-incidf_prepro(incidf=incidf,undetected=undetected,negatives=negatives,smooth_method=smooth_method,smooth_value=smooth,ignore=ignore)$df[,c("dates","I")]}
  t_start<-seq(2, nrow(incidf)-window-1)
  config<-make_config(list(
    mean_prior=mean_prior,
    std_prior=std_prior,
    t_start =t_start,
    t_end = t_start + window - 1))
  incidf<-incidf[,c("dates","I")]
  if (distribution %in% c("Du","Du Bs","Nishiura","Zhao","Ali","Ali early","Ali late","Ganyani Singapore","Ganyani Tianjin")){
    r_est<-estimate_R(incidf, method="si_from_sample",si_sample=switch(distribution,
                                                                       `Du Bs`=si_distr_boot[,1:mcmc_length],
                                                                       Du=si_sample_mcmc_X$si_sample[,1:mcmc_length],
                                                                       Nishiura = si_sample_mcmc_N$si_sample[,1:mcmc_length],
                                                                       Zhao = si_sample_mcmc_Zhao$si_sample[,1:mcmc_length],
                                                                       Ali = si_sample_mcmc_Ali$si_sample[,1:mcmc_length],
                                                                       `Ali early` = Ali_early_distr[,1:mcmc_length],
                                                                       `Ali late` = Ali_late_distr[,1:mcmc_length],
                                                                       `Ganyani Singapore` = G_sing_distr[,1:mcmc_length],
                                                                       `Ganyani Tianjin` = G_tian_distr[,1:mcmc_length]),config=config)  }
  else {
    if (distribution %in% c("Ganyani Singapore Pt","Ganyani Tianjin Pt","Li","Liu")){
      config$si_distr<-dgamma(seq(0,30),shape=switch(distribution,
                                                      `Ganyani Singapore Pt`=sing_m$kappa, `Ganyani Tianjin Pt`=tian_m$kappa,
                                                      `Li`=li_m$kappa,`Liu`=liu_m$kappa),#kappa
                              scale=switch(distribution, `Ganyani Singapore Pt`=sing_m$theta, `Ganyani Tianjin Pt`=tian_m$theta,
                                           `Li`=li_m$theta,`Liu`=liu_m$theta))#theta
    }
    else if (distribution %in% c("Du Point","Nishiura Pt","Zhao Pt", "Ali Pt")){
      config$si_distr<-dlnorm(seq(0,30),meanlog=switch(distribution,
                                                       `Du Point`=SI_fit_clever_X@ests["meanlog","est"],
                                                       `Nishiura Pt`=SI_fit_clever_Nishiura@ests["meanlog","est"],
                                                       `Zhao Pt`=SI_fit_clever_Zhao@ests["meanlog","est"],
                                                       `Ali Pt`=SI_fit_clever_Ali@ests["meanlog","est"]),
                              sdlog=switch(distribution,
                                           `Du Point`=SI_fit_clever_X@ests["sdlog","est"],
                                           `Nishiura Pt`=SI_fit_clever_Nishiura@ests["sdlog","est"],
                                           `Zhao Pt`=SI_fit_clever_Zhao@ests["sdlog","est"],
                                           `Ali Pt`=SI_fit_clever_Ali@ests["sdlog","est"]))
      }
      else if (distribution %in% c("Ali early Pt","Ali late Pt")){
        config$si_distr<-dtruncnorm(seq(0,30),a=1e-6, mean=switch(distribution,
                                                        `Ali early Pt`=ali_early_isolation$mean,
                                                        `Ali late Pt`=ali_late_isolation$mean),
                                sd=switch(distribution,
                                          `Ali early Pt`=ali_early_isolation$sd,
                                          `Ali late Pt`=ali_late_isolation$sd))
      }
    config$si_distr<-config$si_distr/sum(config$si_distr)
    r_est<-estimate_R(incidf, method="non_parametric_si",config=config)
  }  
  if (return_object) {config<-NA;gc();return(r_est)}
  else{
  resultado<-data.frame(date_start=r_est$dates[r_est$R$t_start],date_end=r_est$dates[r_est$R$t_end],r_est$R[ ,c("Mean(R)","Quantile.0.025(R)", "Quantile.0.975(R)")]);
  incidf<-config<-r_est<-NA;
  gc();
  return(resultado)
    }

}

calculate_Empirical_R<-function(incidf,mode=""){
  for (i in (7:(nrow(incidf)-1))){
    incidf[i,"eR"]<-sum(incidf[(i-1):(i+1),"I"])/sum(incidf[(i-6):(i-4),"I"])
    incidf[i-3,"rho7"]<-sum(incidf[((i-6):(i)),"eR"])/7
  }
  gc();
  return(incidf[,c("dates","rho7","Country")])
}

calculate_IRR<-function(incidf,IRR_c=c(0,6)){ # Ratio de tasas
  result<-data.frame(dates=c(),IRR=c(),IRR_N=c(),Country=c())
  for (IRR_N in IRR_c){
    incidf <- incidf %>% mutate(IRR_N=IRR_N,IRR=NA)
      for (i in ((2*(IRR_N+1)):(nrow(incidf)))){
        incidf[i,"IRR"]<-sum(incidf[(i):(i-IRR_N),"I"])/sum(incidf[(i-(IRR_N+1)):(i-(2*IRR_N+1)),"I"])
      }
    result<-rbind(result,incidf[,c("dates","IRR","IRR_N","Country")])
  }
  gc();
  return(result)
}




calculate_Ind<-function(orig_incidf,n_list=14,population=1,pop_factor=10^5){
  #orig_incidf<-incidf
  results<-data.frame(dates=c(),Ind=c(),Country=c(),n=c())
  for (n in n_list){
    incidf<-orig_incidf %>% dplyr::mutate(n=n)
    if (!(("Country") %in% colnames(orig_incidf))){incidf$Country<-""}
    incidf[(n:nrow(incidf)),"Ind"]<-0
    for(i in (0:(n-1))){
      incidf[(n:nrow(incidf)),"Ind"]<-incidf[(n:nrow(incidf)),"Ind"]+incidf[(n-i):(nrow(incidf)-i),"I"]
    }
    if((length(population)==0) || (is.na(population))) {incidf["Ind"]<-(0-1)}
    else {incidf["Ind"]<-incidf["Ind"]/(population)*(pop_factor)}
    results<-rbind(results,incidf[,c("dates","Ind","Country","n")])
    gc();
  }
  gc();
  return(results %>% dplyr::mutate(population=population))
}

calculate_poor_R <- function(incidf,distribution){
  if (distribution %in% c("Ganyani Singapore","Ganyani Singapore Pt","Ganyani Tianjin","Ganyani Tianjin Pt","Liu","Li")){
  single_si_distr<-dgamma(seq(0,30),shape=switch(distribution, 
                                                 `Ganyani Singapore`=sing_m$kappa,`Ganyani Singapore Pt`=sing_m$kappa,
                                                 `Ganyani Tianjin`=tian_m$kappa,`Ganyani Tianjin Pt`=tian_m$kappa,
                                                 `Liu`=liu_m$kappa,`Li`=li_m$kappa),#kappa
                          scale=switch(distribution,
                                       `Ganyani Singapore`=sing_m$theta,`Ganyani Singapore Pt`=sing_m$theta,
                                       `Ganyani Tianjin`=tian_m$theta,`Ganyani Tianjin Pt`=tian_m$theta,
                                       `Liu`=liu_m$theta,`Li`=li_m$theta))#theta
}else if (distribution %in% c("Du","Nishiura","Zhao","Du Point","Nishiura Pt","Zhao Pt","Ali","Ali Pt","Ali")){
  single_si_distr<-dlnorm(seq(0,30),meanlog=switch(distribution,
                                                   Du=SI_fit_clever_X@ests["meanlog","est"],`Du Point`=SI_fit_clever_X@ests["meanlog","est"], `Du Bs`=SI_fit_clever_X@ests["meanlog","est"],
                                                   Nishiura=SI_fit_clever_Nishiura@ests["meanlog","est"],`Nishiura Pt`=SI_fit_clever_Nishiura@ests["meanlog","est"],
                                                   Zhao=SI_fit_clever_Zhao@ests["meanlog","est"],`Zhao Pt`=SI_fit_clever_Zhao@ests["meanlog","est"],
                                                   Ali=SI_fit_clever_Ali@ests["meanlog","est"], `Ali Pt`=SI_fit_clever_Ali@ests["meanlog","est"]),
                          sdlog=switch(distribution,
                                       Du=SI_fit_clever_X@ests["sdlog","est"],`Du Point`=SI_fit_clever_X@ests["sdlog","est"],`Du Bs`==SI_fit_clever_X@ests["sdlog","est"],
                                       Nishiura=SI_fit_clever_Nishiura@ests["sdlog","est"],`Nishiura Pt`=SI_fit_clever_Nishiura@ests["sdlog","est"],
                                       Zhao=SI_fit_clever_Zhao@ests["sdlog","est"],`Zhao Pt`=SI_fit_clever_Zhao@ests["sdlog","est"],
                                       Ali=SI_fit_clever_Ali@ests["meanlog","est"], `Ali Pt`=SI_fit_clever_Ali@ests["meanlog","est"]))
} else if (distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt")){
  single_si_distr<-dtruncnorm(seq(0,30),a=1e-6, mean=switch(distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt"),
                                                            ali_early_isolation$mean,ali_early_isolation$mean,
                                                            ali_late_isolation$mean,ali_late_isolation$mean),
                              sd=switch(distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt"),
                                        ali_early_isolation$sd,ali_early_isolation$sd,
                                        ali_late_isolation$sd,ali_late_isolation$sd))
}

}

calculate_Rc<-function(incidf,distribution="Du",method="TD"){
  td_distribution<-switch(distribution,`Du Point`="Du",`Du Bs`="Du",`Nishiura Pt`="Nishiura",`Zhao Pt`="Zhao",`Ali Pt`="Ali",`Ali early Pt`="Ali early",`Ali late Pt`="Ali late",`Ganyani Tianjin Pt`="Ganyani Tianjin", `Ganyani Sinagpore Pt`="Ganyani Singapore",distribution)
  #print(incidf)
  est.r<-tryCatch(estimate.R(epid=incidf$I, t=incidf$dates, begin=1,end=nrow(incidf)-1,method = method,GT = switch(td_distribution,
                                                                                                                   Du=generation.time("empirical",val= diff(plnorm(c(0,0.5+c(0:31)),
                                                                                                                                                                   meanlog=SI_fit_clever_X@ests[c("meanlog"),"est"],
                                                                                                                                                                   sdlog=SI_fit_clever_X@ests[c("sdlog"),"est"]))),
                                                                                                                   Nishiura=generation.time("empirical",val= diff(plnorm(c(0,0.5+c(0:31)),
                                                                                                                                                                         meanlog=SI_fit_clever_Nishiura@ests[c("meanlog"),"est"],
                                                                                                                                                                         sdlog=SI_fit_clever_Nishiura@ests[c("sdlog"),"est"]))),
                                                                                                                   Zhao=generation.time("empirical", val = diff(plnorm(c(0,0.5+c(0:31)),
                                                                                                                                                                       meanlog=SI_fit_clever_Zhao@ests[c("meanlog"),"est"],
                                                                                                                                                                       sdlog=SI_fit_clever_Zhao@ests[c("sdlog"),"est"]))),
                                                                                                                   Ali=generation.time("empirical", val = diff(plnorm(c(0,0.5+c(0:31)),
                                                                                                                                                                       meanlog=SI_fit_clever_Ali@ests[c("meanlog"),"est"],
                                                                                                                                                                       sdlog=SI_fit_clever_Ali@ests[c("sdlog"),"est"]))),
                                                                                                                   Liu=generation.time("empirical", val = diff(pgamma(c(0,0.5+c(0:34)),
                                                                                                                                                                      shape=liu_m$kappa,
                                                                                                                                                                      scale=liu_m$theta))),
                                                                                                                   Li=generation.time("empirical", val = diff(pgamma(c(0,0.5+c(0:38)),
                                                                                                                                                                      shape=li_m$kappa,
                                                                                                                                                                     scale=li_m$theta))),
                                                                                                                   `Ali early`=generation.time("empirical", val = diff(ptruncnorm(c(0,0.5+c(0:31)),a=1e-6,
                                                                                                                                                                                  mean=ali_early_isolation$mean,
                                                                                                                                                                                  sd=ali_early_isolation$sd))),
                                                                                                                   `Ali late`=generation.time("empirical", val = diff(ptruncnorm(c(0,0.5+c(0:31)),a=1e-6,
                                                                                                                                                                                 mean=ali_late_isolation$mean,
                                                                                                                                                                                 sd=ali_late_isolation$sd))),
                                                                                                                   `Ganyani Singapore`=generation.time("empirical", val = diff(pgamma(c(0,0.5+c(0:31)),
                                                                                                                                                                                      shape=sing_m$kappa,
                                                                                                                                                                                      scale=sing_m$theta))),
                                                                                                                   `Ganyani Tianjin`=generation.time("empirical", val = diff(pgamma(c(0,0.5+c(0:31)),
                                                                                                                                                                                    shape=tian_m$kappa,
                                                                                                                                                                                    scale=tian_m$theta)))
  )),error=function(e) {print("error"); print(e);NA})
  #gc();
  return(est.r)
  
}



comparemodes<-function(countries,mode,moderesults=list(data_to_plot=data.frame(Country=c(),Method=c(),distr=c(),Mode=c(),date_start=c(),date_end=c(),Mean.R.=c(),`Quantile.0.025.R.`=c(),`Quantile.0.975.R.`=c()),
                                                       OA=data.frame(Country=c(),distr=c(),Mode=c(),dates=c(),`Overall Infectivity`=c(),OA_ratio=c()),
                                                       incidences=data.frame(dates=c(),Incidences=c(),I_pre_smooth=c(),Country=c(),Mode=c()),
                                                       empiricalR=data.frame(dates=c(),rho7=c(),Country=c(),Mode=c()),
                                                       `n_cum`=data.frame(dates=c(),Ind=c(),Country=c(),n=c(),Mode=c(),population=c()),
                                                       irr=data.frame(dates=c(),IRR=c(),IRR_N=c(),Country=c(),Mode=c())),
                       distributions=c("Du", "Nishiura","Zhao","Ganyani Singapore","Ganyani Tianjin","Liu"),TD=TRUE,windows=c(7),
                       with_inci=FALSE,with_rho7=FALSE,with_n_cum=FALSE, with_irr=FALSE, pop_data=world_pop, loc_column="Location",pop_column="PopTotal"){
  if (offline){
  cl<-makeCluster(3,outfile="C:\\TEMP\\CLSUTER.TXT")
  clusterEvalQ(cl,{library(R0);library(EpiEstim);library(dplyr);library(truncnorm)})
} 

clusterExport(cl=cl,list("mcmc_length",#"pop_data","pop_column","loc_column",
                         "tian_distros","li_m","liu_m","Ali_early_distr","ali_early_isolation","Ali_late_distr","ali_late_isolation",
                         #"Ali_early_distr","Ali_late_distr",
                         "G_sing_distr","G_tian_distr","tian_m","sing_m",
                         "calculate_R","calculate_Rc","calculate_Empirical_R","calculate_Ind", "calculate_IRR",
                         "si_sample_mcmc_X","si_sample_mcmc_Ali","si_sample_mcmc_N","si_sample_mcmc_Zhao",
                         "SI_fit_clever_X","SI_fit_clever_Nishiura","SI_fit_clever_Zhao","SI_fit_clever_Ali"))

  #if(with_rho7) moderesults[["empiricalR"]]<-data.frame(dates=c(),rho7=c(),Country=c())
  #if(with_n_cum>0) moderesults[["n_cum"]]<-data.frame(dates=c(),Ind=c(),Country=c())
  for(i in names(countries)) countries[[i]]<-countries[[i]] %>% dplyr::mutate(Country=i)
  
  #mcmc_length<-20
  par_comparison<-parLapply(cl,X=countries,fun=calculate_a_country,with_inci=with_inci,with_rho7=with_rho7,with_n_cum=with_n_cum,
                            distributions=distributions,windows=windows,mode=mode,TD=TD,with_irr=with_irr, pop_data=pop_data, loc_column=loc_column, pop_column=pop_column)
  #par_comparison<-lapply(X=countries,FUN=calculate_a_country,with_inci=with_inci,with_rho7=with_rho7,with_n_cum=with_n_cum,
  #                          distributions=distributions,windows=windows,mode=mode,TD=TD,with_irr=with_irr)
  #print("calc ended")
  for(key in names(moderesults)){
    list_of_dfs<-list()
    for (country in names(par_comparison)){
      if(((length(par_comparison[[country]]))>0) && (length(par_comparison[[country]][[key]])>0) && (nrow(par_comparison[[country]][[key]])>0)) moderesults[[key]]<-rbind(moderesults[[key]],par_comparison[[country]][[key]])
      #list_of_dfs[[country]]<-par_comparison[[country]][[key]]
    }
    #moderesults[[key]]<-dplyr::bind_rows(moderesults[[key]],list_of_dfs)
  }
  stopCluster(cl)
  return(moderesults)
}



calculate_a_country<-function(incidf,with_inci,with_rho7,with_n_cum,distributions,windows,mode,TD,with_irr,pop_data=world_pop,pop_column,loc_column){
  #print("caa")
  country<-unique(incidf[,"Country"])
  #incidf <-  incidf %>% dplyr::select(-Country)
  moderesults<-list(data_to_plot=data.frame(Country=c(),Method=c(),distr=c(),Mode=c(),date_start=c(),date_end=c(),Mean.R.=c(),`Quantile.0.025.R.`=c(),`Quantile.0.975.R.`=c()),
                    OA=data.frame(Country=c(),distr=c(),Mode=c(),dates=c(),`Overall Infectivity`=c(),OA_ratio=c()),
                    incidences=data.frame(dates=c(),Incidences=c(),I_pre_smooth=c(),Country=c(),Mode=c()),
                    #empiricalR=data.frame(dates=c(),rho7=c(),Country=c(),Mode=c()),
                    #n_cum=data.frame(dates=c(),Ind=c(),Country=c(),Mode=c(),population=c()),
                    irr=data.frame(dates=c(),IRR=c(),IRR_N=c(),Country=c(),Mode=c()))
  if (with_rho7 & (nrow(incidf)>8)) moderesults$empiricalR<-calculate_Empirical_R(incidf) %>% dplyr::mutate(Mode=mode)
  if ((with_n_cum) & (nrow(incidf)>with_n_cum)) moderesults$n_cum <- calculate_Ind(incidf,n=with_n_cum,population = if (country %in% pop_data[,loc_column]) {pop_data[pop_data[,loc_column]==country,pop_column]} else 1)  %>% dplyr::mutate(Mode=mode) 
  if (any(with_irr) & (nrow(incidf)>max(with_irr))) moderesults$irr <- calculate_IRR(incidf=incidf,IRR_c=with_irr)  %>% dplyr::mutate(Mode=mode)
  if (with_inci)
    moderesults$incidences <- incidf  %>% dplyr::rename(Incidences=I) %>% dplyr::mutate(Mode=mode)
  if (nrow(incidf)==0) {fix(incidf);stop(paste("Incidence empty",country, mode, dim(incidf)))}
  if (sum(incidf$I<0)>0) {fix(incidf);stop(paste("Incidence negative ",country, mode))}
  gc();
  for (distribution in distributions){
    #print(incidf)
    if (distribution %in% c("Ganyani Singapore","Ganyani Singapore Pt","Ganyani Tianjin","Ganyani Tianjin Pt","Liu","Li")){
      single_si_distr<-dgamma(seq(0,30),shape=switch(distribution, 
                                                     `Ganyani Singapore`=sing_m$kappa,`Ganyani Singapore Pt`=sing_m$kappa,
                                                     `Ganyani Tianjin`=tian_m$kappa,`Ganyani Tianjin Pt`=tian_m$kappa,
                                                     `Liu`=liu_m$kappa,`Li`=li_m$kappa),#kappa
                              scale=switch(distribution,
                                           `Ganyani Singapore`=sing_m$theta,`Ganyani Singapore Pt`=sing_m$theta,
                                           `Ganyani Tianjin`=tian_m$theta,`Ganyani Tianjin Pt`=tian_m$theta,
                                           `Liu`=liu_m$theta,`Li`=li_m$theta))#theta
    }else if (distribution %in% c("Du","Nishiura","Zhao","Du Point","Nishiura Pt","Zhao Pt","Ali","Ali Pt","Ali")){
      single_si_distr<-dlnorm(seq(0,30),meanlog=switch(distribution,
                                                       Du=SI_fit_clever_X@ests["meanlog","est"],`Du Point`=SI_fit_clever_X@ests["meanlog","est"], `Du Bs`=SI_fit_clever_X@ests["meanlog","est"],
                                                       Nishiura=SI_fit_clever_Nishiura@ests["meanlog","est"],`Nishiura Pt`=SI_fit_clever_Nishiura@ests["meanlog","est"],
                                                       Zhao=SI_fit_clever_Zhao@ests["meanlog","est"],`Zhao Pt`=SI_fit_clever_Zhao@ests["meanlog","est"],
                                                       Ali=SI_fit_clever_Ali@ests["meanlog","est"], `Ali Pt`=SI_fit_clever_Ali@ests["meanlog","est"]),
                              sdlog=switch(distribution,
                                           Du=SI_fit_clever_X@ests["sdlog","est"],`Du Point`=SI_fit_clever_X@ests["sdlog","est"],`Du Bs`==SI_fit_clever_X@ests["sdlog","est"],
                                           Nishiura=SI_fit_clever_Nishiura@ests["sdlog","est"],`Nishiura Pt`=SI_fit_clever_Nishiura@ests["sdlog","est"],
                                           Zhao=SI_fit_clever_Zhao@ests["sdlog","est"],`Zhao Pt`=SI_fit_clever_Zhao@ests["sdlog","est"],
                                           Ali=SI_fit_clever_Ali@ests["meanlog","est"], `Ali Pt`=SI_fit_clever_Ali@ests["meanlog","est"]))
    } else if (distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt")){
      single_si_distr<-dtruncnorm(seq(0,30),a=1e-6, mean=switch(distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt"),
                                                                ali_early_isolation$mean,ali_early_isolation$mean,
                                                                ali_late_isolation$mean,ali_late_isolation$mean),
                                  sd=switch(distribution %in% c("Ali early","Ali early Pt","Ali late","Ali late Pt"),
                                            ali_early_isolation$sd,ali_early_isolation$sd,
                                            ali_late_isolation$sd,ali_late_isolation$sd))
    }
    single_si_distr<-single_si_distr/sum(single_si_distr)
    moderesults$OA<-rbind(moderesults$OA,
                          data.frame(Country=country,distr=distribution,Mode=mode,dates=incidf$dates,
                                     `Overall Infectivity`=overall_infectivity(incidf[,c("dates","I")],si_distr=single_si_distr)) %>%
                            mutate (OA_ratio=if (country %in% pop_data[,loc_column]) (Overall.Infectivity/(pop_data[pop_data[loc_column]==country,pop_column])*(10^6))
                          else NA))    
    for (window in windows){
      if (nrow(incidf)>=window){
        #print(distribution)
        #r_est<-estimate_R(incidf, method="si_from_sample",si_sample=si_distribution,config=config)
        
        moderesults$data_to_plot<-rbind(moderesults$data_to_plot,calculate_R(incidf=incidf[,c("dates","I")],prepro=FALSE,distribution=distribution,window=window) %>%
          mutate(Method=paste("Cori, w=",window),distr=distribution, Mode=mode, Country=country))
        

      gc()
      } 
    }
    if (((length(TD)==1) & (TD!=0)) | length(TD)>1){
      #print(TD)
      if ((length(TD)==1) & (TD==TRUE)) TD<-c("TD","SB")
      for (method in TD){
        est.r<-calculate_Rc(incidf=incidf,distribution=distribution,method=method)
        if(!is.na(est.r)){
          moderesults$data_to_plot<-rbind(moderesults$data_to_plot,
                                        data.frame(Country=country,Method=if (method=="TD") "Wallinga Teunis" else if (method=="SB") "Bettencourt Ribeiro",distr=distribution,
                                             Mode=mode,
                                             date_start=est.r$estimates[[method]]$epid$t[est.r$estimates[[method]]$begin.nb :if (method=="SB") est.r$estimates[[method]]$end.nb-1 else est.r$estimates[[method]]$end.nb],
                                             date_end=est.r$estimates[[method]]$epid$t[(if (method=="SB") est.r$estimates[[method]]$begin.nb+1 else  est.r$estimates[[method]]$begin.nb):est.r$estimates[[method]]$end.nb],
                                             Mean.R.=est.r$estimates[[method]]$R,
                                             if (method=="TD") est.r$estimates[[method]]$conf.int %>%
                                               rename(Quantile.0.025.R.=lower, Quantile.0.975.R.=upper) else
                                                 (est.r$estimates[[method]]$conf.int %>%
                                                    rename(Quantile.0.025.R.=CI.lower., Quantile.0.975.R.=CI.upper.))))
        
          remove(est.r)
          gc()
        }
      }                          
    }
  }
  gc()
  #print(paste(country," ok"))
  return(moderesults)
}


#------------------------------------
mondaysdf<-data.frame(monday=seq(as.Date("2020-01-06"),as.Date("2021-01-04"),7),friday=seq(as.Date("2020-01-10"),as.Date("2021-01-08"),7))
add_weekdays<-function(plot,mondays=mondaysdf,mindate=min(mondaysdf$monday),maxdate=max(mondaysdf$monday),fill="pink",alpha=0.03,border=TRUE){
  df_selected<-subset(mondays,(monday<maxdate) & (monday>mindate))
  if (border)
    for (i in (1:nrow(df_selected))){plot<-plot+geom_rect(aes_string(xmin=df_selected[i,"monday"],xmax=df_selected[i,"friday"],ymin=0,ymax=Inf),fill=fill,alpha=alpha)}
  else
    for (i in (1:nrow(df_selected))){plot<-plot+geom_rect(aes_string(xmin=df_selected[i,"monday"],xmax=df_selected[i,"friday"],ymin=0,ymax=Inf),fill=fill,alpha=alpha,color=NA)}
  return(plot)
}

#--------------------------------------
my_theme_12<- ggplot2::theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12,angle=45),
                             axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12,angle=45),
                             strip.text = element_text(size=12))
my_theme<-  ggplot2::theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14,angle=45),
                           axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
                           legend.text=element_text(size=14),legend.title=element_text(size=14),strip.text = element_text(size = 14))

my_theme_13<-  ggplot2::theme(axis.title.x = element_text(size = 13), axis.text.x = element_text(size = 13,angle=45),
                              axis.title.y = element_text(size = 13), axis.text.y = element_text(size = 13),
                              legend.text=element_text(size=13),legend.title=element_text(size=13),strip.text = element_text(size = 13))


my_theme_15<-  ggplot2::theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15,angle=45),
                              axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15),
                              legend.text=element_text(size=15),legend.title=element_text(size=14),strip.text = element_text(size = 15))
theme_bottom_legend<-ggplot2::theme(legend.position="bottom",legend.key.width =  unit(4, "lines")
)



#--------------------------------------

incidencePlot <- function(x,comparemodes,shiny=TRUE){
  #print(comparemodes$incidences %>% filter(Country==x,Incidences>0,complete.cases(Incidences)))
  #print(paste("Rendering",x))
  df<-comparemodes$raw_incidences %>% dplyr::filter((country==x),I>1,I_pre_smooth>0,complete.cases(I,I_pre_smooth)) %>%
    dplyr::mutate(I_pre_smooth=ifelse(I_pre_smooth<1,NA,I_pre_smooth))

    the_plot<-add_weekdays(ggplot(df, aes(x = dates, y = I+0.01, group = factor(country))), mindate=min(comparemodes$raw_incidences$dates),maxdate=max(comparemodes$raw_incidences$dates)) +
      geom_line(size=1.7,alpha=0.8) +
      geom_point(aes(y = I_pre_smooth+0.01),size=2,alpha=0.8) +
      #geom_text(aes(x = as.Date("2020-06-21"),y=1000,label="2020-06-21"),text=element_text(size=14),vjust=-1,hjust=1)+
      ggtitle(paste("Daily incidence in",x))+#s in Catalonia, Spain",strftime(Sys.Date(),"%d-%b-%Y")))+ 
      geom_vline(xintercept=as.Date("2020-06-21"),size=1)+
      labs(x=element_blank(),y=element_blank()) + my_theme_15 +
      scale_y_continuous(trans='log10')
    if(shiny=="FALSE") return(the_plot)
    return(renderCachedPlot({the_plot
          },sizePolicy=sizeGrowthRatio(540,480,1.1),cache="app",cacheKeyExpr = {comparemodes$data_to_plot %>% filter(Country==x)}))
  }

transmissionPlot <-function(x,comparemodes,shiny=TRUE){
  df<-comparemodes$data_to_plot[complete.cases(comparemodes$data_to_plot) &
                                  ((comparemodes$data_to_plot$Country==x)),]# %>% dplyr::mutate(date_end=as.Date(date_end))
  the_plot<- add_weekdays(ggplot(df, aes(x = date_end, y = as.numeric(`Mean.R.`)+0.01 )),mindate=min(df$date_end),maxdate=max(df$date_end)) +
    geom_line(aes(y = Mean.R.),size=1.3) +
    geom_ribbon(aes(ymin = `Quantile.0.025.R.`, ymax = `Quantile.0.975.R.`),alpha=0.05) +
    geom_line(aes(y = `Quantile.0.025.R.`),linetype=2,alpha=0.4) +
    geom_line(aes(y = `Quantile.0.975.R.`),linetype=2,alpha=0.4) +
    geom_hline(yintercept=1, linetype="dashed",  color = "black", size=1.05) +
    #ggtitle("Estimated R")+scale_fill_discrete(name = "Country", labels = EU) +
    lims(y=c(0,min(max(df$Mean.R.),4))) + geom_line(aes(y = Mean.R.),color="black",size=0.1,alpha=0.8) +
    annotate(geom="label",x = max(df$date_end)+2,y=tail(df,1)$Mean.R.+0.5,label=round(tail(df,1)$Mean.R.,3),fontface ="bold", color = "black", size = 5) +
    my_theme_15 +
    #geom_vline(xintercept=as.Date("2020-06-21"),size=1)+
    #geom_vline(xintercept=as.Date("2020-07-04"),size=1)+
    ggtitle(paste("estimated Rt, ",x,", ",Sys.Date(),sep=""))+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="bottom")
  if(shiny=="FALSE") return(the_plot)
  return(renderCachedPlot({the_plot
      #geom_rect(mapping=aes(xmin=min(comparemodes$data_to_plot$dates),xmax=min(comparemodes$data_to_plot$dates),ymin=0,ymax=4,fill="Weekdays"),alpha=0.3)
  },sizePolicy=sizeGrowthRatio(540,480,1.1),cache="app",cacheKeyExpr = {comparemodes$data_to_plot %>% filter(Country==x)}))
}

OAPlot<-function(x,comparemodes){
  df<-comparemodes$OA[complete.cases(comparemodes$OA) &
                                  ((comparemodes$OA$Country == x)),]# %>% dplyr::mutate(date_end=as.Date(date_end))
  return(renderCachedPlot({
  add_weekdays(ggplot(df, aes(x = dates, y = as.numeric(`Overall.Infectivity`))),mindate=min(df$dates),maxdate=max(df$dates)) +
    geom_line(size=1.3) +
    #geom_line(aes(y = Mean.R.,color=Regions),size=1.3,alpha=0.8) +
    #ggtitle("Estimated R")+scale_fill_discrete(name = "Country", labels = EU) +
    ggtitle(paste("Overall Infectivity,",x))+ 
      labs(x=element_blank(),y=element_blank())  +  #ylim(c(0.2,4.8)) + geom_line(aes(y = Mean.R.),color="black",size=0.1,alpha=0.8) +
    my_theme_15

  
  },sizePolicy=sizeGrowthRatio(540,480,1.1),cache="app",cacheKeyExpr = {comparemodes$OA %>% filter(Country==x)}))
}