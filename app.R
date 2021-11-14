#10/07/2021
library(shiny)
library(R0)
library(EpiEstim)# GPL >2
#library(incidence)
#library(stringr)
library(ggplot2)
#library(googlesheets)
library(DT)
source("libraryP.R",encoding="UTF-8")
library(openxlsx)
library(MASS)
library(dplyr)
library(shinythemes)
library(shinycssloaders) #GPL 3
#library(countrycode)
# Packages: shiny, R0, EpiEstim, ggplot2, DT, openxlsx,MASS, dplyr,tidyr,shinythemes,shinycssloaders, TTR, incidence
debug<- FALSE


par(mar=c(3.3,1.3,1.1,1.1))

ccovid19<-read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",stringsAsFactors = TRUE)
#dcovid19<-read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
#rcovid19<-read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")

#casal<-load(url("https://raw.githubusercontent.com/rubenfcasal/COVID-19/master/acumula22.RData"))

dhub_csv<-USAs_csv<-USA_csv<-cat_m<-cat_z<-msc_csv<-incidfl<-narr_csv<-covid19<-alvaro<-ukcsv<-NA
cat_m_date<-cat_z_date<-NA
acc_plot<-NA 
warnings<-""



mcmc_limit<<-1:500

ccovid19$Province.State<-NAtoAll(ccovid19$Province.State)
#dcovid19$Province.State<-NAtoAll(dcovid19$Province.State)
#rcovid19$Pro vince.State<-NAtoAll(rcovid19$Province.State)


load("si_sample_mcmc_Zhao.R")
load("si_sample_mcmc_X.R")
load("si_sample_mcmc_N.R")
load("si_sample_mcmc_Ali.R")
#Nishiura



# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("journal"),
    tags$head(
      tags$meta(charset="UTF-8"),
      tags$meta(name="description", content="COVID-19 time dependent reproductive ratio (Rt-R0) estimation tool"),
      tags$meta(name="google-site-verification", content="E2IN_HT-O2cqnhFRkSM14_SpxUrTkHOQ2Q0bK3kU4Ys"),
      tags$meta(name="keywords", content="R, Shiny, COVID-19, SARS-CoV-2, coronavirus, reproductive ratio, R0, estimation, COVID, transmission"),
      tags$meta(name="author", content="Luis Alfredo Bautista Balbas"),
      tags$meta(name="viewport", content="width=device-width, initial-scale=1.0")
    ),
    # Application title
    titlePanel("COVID-19 Rt estimator"),
    
    # Sidebar with a slider input for number of bins 
    #tabsetPanel(
    #shinythemes::themeSelector(),  
    # tabPanel("Custom estimation",
    sidebarLayout(
      sidebarPanel=sidebarPanel(
        selectInput(inputId="datar", label="Data Repository", choices=sources_list),
        selectInput(inputId='country', label='Country', choices="",selected="US"),
        selectInput(inputId='province', label='State/Province', choices=c("")),
        selectInput("smooth_method", label="Smooth method",choices=c("SMA (centered)","Lowess","SMA (past)","SMA (future)","EMA (past)","EMA (future)","WMA (past)","WMA (future)","Lowess","No smoothing"),selected="SMA (centered)"),
        sliderInput("smooth_value", label="Number of days to apply smooth",
                    min = 0, max = 21, value = 2, step = 1),
        sliderInput("undetected", label="Proportion of undetected cases",
                    min = 0, max = 0.99, value = 0.0),
        sliderInput("ignore", "Start after N values/3 days",
                    min = 0, max = 20, value = 3, step =1),
        dateInput("start_date", "Start date",
                  min = "2020-02-15", max = Sys.Date(), value = "2020-03-05"),
        dateInput("end_date", "End date",
                  min = "2020-02-15", max = Sys.Date(), value = "2021-11-14"),
        selectInput("negatives","Negative values",choices=c("Distribute backwards before smoothing","Do nothing","Remove before smoothing","Remove after smoothing","to 0 before smoothing","to 0 after smoothing"),selected="to 0 before smoothing"),        
        selectInput("positives","Spurious positive peaks",choices=c("Distribute backwards before smoothing","Do nothing"),selected="Do nothing"),
        selectInput("method",label="Estimation method", choices=c("A Cori et al (2013)","Wallinga and Teunis (2004)"),selected="A Cori et al (2013)"),#"Wallinga and Teunis (2004)","Bettencourt and Ribeiro (2008)
        sliderInput("window", label="Estimation window (days) (A Cori et al)",
                    min = 1, max = 30, value = 7, step =1),
        selectInput("si_distr", label="Serial interval distribution",
                    choices=si_strs,selected=Ln_ni_str)),
      mainPanel=mainPanel(
        h2(textOutput("country_province")),
        tabsetPanel(type="tabs",
                    tabPanel("Plot",textOutput("warnings"), textOutput("report"),
                             plotOutput("plot_smooth",height="600px") %>% withSpinner(color="darkred"),
                             plotOutput("miplot",height="600px") %>% withSpinner(color="darkred"),
                             plotOutput("cioi",height="600px") %>% withSpinner(color="darkred"),
                             hr(),
                             textOutput("summary"),hr(),
                             textOutput("source")),
                    tabPanel("Rt Table", p("Estimated Rt values"),dataTableOutput('tabla'),downloadButton("downloadRt", "Download Table")),
                    tabPanel("Data Table", dataTableOutput('data')),
                    tabPanel("Active prevalence table",
                             p("The indicators of prevalence of active cases include: Cumulative incidence rates (CIR) of the previous 14 days, CIR of the previous 7 days, overall infectivity as calculated via EpiEstim package"),
                             dataTableOutput("cioidf"),
                             downloadButton("downloadcioi", "Download Table")),
                    tabPanel("Serial interval distribution",
                             p("The serial interval is the interval between case onsets of a primary case and another subsequent secondary case. The distribution is of this interval is assumed to be the following:"),
                             plotOutput('distro'))
                    #tabPanel("Provinces summary"),
                    #plotOutput("plot_provinces",  actionButton("provinces_btn", "Estimate R0 for all provinces in that country"), plotOutput('provinces_plot'),dataTableOutput('provinces_table'))),
                    
        )))
    #)#,tabPanel("Predefined maps",plotOutput("USA",height="800px"))
    #  ),#Layout panels 
    ,hr(),
    p("This tool uses publicly available data from third party sources related to the Coronavirus pandemia. The preprocessing steps include filtering by dates, estimating a real value (assuming a specified rate of asymptomatic cases), negative values management and time series smoothing. Several functions (EMA, SMA, WMA) can perform smoothing, and moving averages can be applied using past values (the usual financial markets technique) or future values (a deviation from standard technique that might make more sense for an underreported disease). The EpiEstim package by Anne Cori is used to estimate the time-dependent R0 values. The serial interval is the interval between symptoms-onset in an infector and symptoms-onset in an infectee, and is required for R estimation. Available distributions are:"),
    tags$div(
      tags$ul(
        tags$li("Nishiura, Linton and Akhmetzhanov (2020) (lognormal distribution with mean = 4.7 days and SD = 2.9 days), available as single distribution or multiple via MCMC"),
        tags$li("Du, Xu, Wu et al (2020) (lognormal distribution, available as single distribution or multiple via bootstrapping or via MCMC)"),
        tags$li("Zhao et al (2020) (lognormal distribution), available as single distribution or multiple via MCMC"),
        tags$li("Ali et al (2020) (all samples (available as single distribution or multiple via MCMC), or only those with late or early isolation(single distributions))"),
        tags$li("Ganyani et al (2020) (gamma distribution, that estimates the generation time distribution, obtained from data from Tianjin or from Singapore"),
      )
    ),
    p("Multiple distributions (up to 200 in this app) can be obtained via Markov-Chain Monte Carlo methods; these are used to estimate the R_t. The EpiEstim default weekly sliding window is used to estimate the R0."),
    hr(),
    p("Options explanation:"),
    hr(),
    p("Smooth methods include: Lowess regression (R lowess function) and moving averages (EMA, SMA, WMA). When the lowess function is selected, the numbeer of days that influence every day is specified. When a moving average the number of days used to calculate the average (before (=past) or after (=future) a given day is selected. If you apply smoothing, you are assuming the real COVID-19 incidence corresponds to the smoothed time series, for example due to the reduced diagnostic workload on weekends."),
    hr(),
    p("Another parameter is the proportion of undetected cases, set at 40% by default (asymptomatic cases in the Iceland Serology study). Every daily incidence is divided by the given proportion of asymptomatic cases to obtain an estimation of total cases. Finally, three day periods are sequentially tested from the start until a period with more than a minimum number of cases is found"),
    hr(),
    p("Time series smoothing and preprocessing are not recommended for all cases, and should be applied on the basis of specific knowledge of the issues present in the time series; and under the assumption that modified values better reflect the real incidences than the original 'raw' values in the datasets."),
    hr(),
    p("Known limitions: Imported cases are not taken into account. Locations are asumed to be isolated"),
    hr(),
    p("2020. Luis Alfredo Bautista Balbas. You can contact the author at luisalfredo d0t bautistabalbas at gmail d0t com."),
    tags$p(tags$a(href="https://www.medrxiv.org/content/10.1101/2020.07.15.20154039v2","Preprint description of the methodology")),
    #p("  //  Data sources: John Hopkins University CSSU, European CDC, Ministerio de Sanidad y Consumo, ISCIII, UK Government, @alvariteus, Narrativa"),
    p("  //  EpiEstim package. Anne Cori et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013), DOI: 10.1093/aje/kwt133DOI . DOI: 10.1016/j.epidem.2019.100356"),
    p("Packages used: shiny, R0, EpiEstim, ggplot2, DT, openxlsx, MASS, dplyr,tidyr,shinythemes, shinycssloaders, TTR, incidence, forecast, ISOcodes, truncnorm. Their respective authors, their respective licenses.")
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  countries <- reactive({
    mydata <- input$datar
    if (debug) print(mydata)
    if (mydata==jhu_str) {return(unique(ccovid19$Country.Region))}
    else if (mydata==ecdc_str) {
      if (is.na(covid19))
      {covid19<<-read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", stringsAsFactors = TRUE)}
      return(unique(covid19$countriesAndTerritories))}
    else if (mydata==narr_str) { # All OK 09/05/2021
      if (is.na(narr_csv)) {
        narr_csv<<-read.csv("https://covid19tracking.narrativa.com/csv/confirmed.csv",encoding="UTF-8",stringsAsFactors = TRUE);
        narr_csv$Region<<-NAtoAll(narr_csv$Region)}
      return(levels(narr_csv$Country_EN))
    }
    #else if (mydata==casal_str || mydata==alvaro_str) {   return("Spain")}
    else if (mydata==uk_str) {
      if (is.na(ukcsv)) {ukcsv<<-read.csv("https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-cases-uk.csv", stringsAsFactors = TRUE)}
      return(unique(ukcsv$Country));}
    else if (mydata==USA_str){
      {if (is.na(USA_csv)){USA_csv<<-read.csv(USA_path,encoding="UTF-8", stringsAsFactors = TRUE)};return(c("* All *",levels(USA_csv$state)))}
    }
    else if (mydata %in% class_sources){
      source<-names(which(class_sources==mydata))
      if (nrow(datasources[[source]]@df)==0) {datasources[[source]]@df<<-datasources[[source]]@getData(datasources[[source]])}
      return(datasources[[source]]@getStates(datasources[[source]]))
    }
    else if (mydata %in% c(CAT_mc_str,CAT_msc_str)){if (is.na(cat_m)){cat_m<<-read.csv(CAT_m_path,encoding="UTF-8", stringsAsFactors = TRUE);cat_m_date<<-Sys.Date()};return(c("* All *",levels(cat_m$ComarcaDescripcio)))}
    else if (mydata %in% c(CAT_zc_str,CAT_zsc_str)){if(is.na(cat_z)){cat_z<<-read.csv(CAT_z_path,encoding="UTF-8", stringsAsFactors = TRUE);cat_z_date<<-Sys.Date()};return(c("* All *",levels(cat_z$RegioSanitariaDescripcio)))}
    #else if (mydata==cam_)
  }) 
  provinces <- reactive({
    mycountry <- input$country
    if (debug) print(mycountry)
    if (input$datar==jhu_str) {return(unique(ccovid19[ccovid19$Country.Region==mycountry,"Province.State"]))}
    else if (input$datar==ecdc_str) {return("")}
    else if (input$datar==narr_str) {
      return (unique(narr_csv[narr_csv$Country_EN==mycountry,"Region"]))
    }
    #else if (input$datar==casal_str) {
    #  if (!exists("acumula2")){load(url("https://raw.githubusercontent.com/rubenfcasal/COVID-19/master/acumula2.RData"))};
    #  if (!exists("acumulados"))load(url("https://raw.githubusercontent.com/rubenfcasal/COVID-19/master/acumulados.RData"));
    #  return(unique(acumulados$ccaa))}
    #else if (input$datar==alvaro_str) {
    #  if (is.na(alvaro)) {
    #    alvaro<-read.csv(url("https://docs.google.com/spreadsheets/d/1aOSfDXbqawSCngJwPHUiRDQEA2QLFUyZ7XJ6P1nl6LI/export?format=csv"),skip=1,row.names=1,encoding="UTF-8");
    #    alvaro<<-alvaro[,substr(colnames(alvaro),1,1)=="X"];}
    #  return(row.names(alvaro))}
    else if (input$datar==uk_str) { return(unique(ukcsv[ukcsv$Country==mycountry,"Area"]));}
    else if (input$datar==USA_str) { return(c("* All *",levels(USA_csv$county)[unique(USA_csv[if (mycountry!="* All *") USA_csv$state==mycountry else TRUE,"county"])]));}
    #else if (input$datar==MSC_str) {msc_csv<<-subset(read.csv(MSC_path),(FECHA!="") & (!(is.na(`PCR.`))));return(c("* All *",as.character(unique(msc_csv$CCAA))))}
    else if (input$datar %in% c(CAT_mc_str,CAT_msc_str)){return(c("* All *",as.character(unique(cat_m[(if (mycountry %in% c("* All *","")) TRUE else cat_m$ComarcaDescripcio==mycountry),"MunicipiDescripcio"]))))}
    else if (input$datar %in% c(CAT_zc_str,CAT_zsc_str)){return(c("* All *",as.character(unique(cat_z[(if (mycountry %in% c("* All *","")) TRUE else cat_z$RegioSanitariaDescripcio==mycountry),"ABSDescripcio"]))))}
    else if (input$datar %in% class_sources){
      source<-names(which(class_sources==input$datar))
      if (nrow(datasources[[source]]@df)==0) {datasources[[source]]@df<<-datasources[[source]]@getData(datasources[[source]])}
      return(datasources[[source]]@getProvinces(datasources[[source]],input$country))
    }
  })
  observe({
    if (debug) print("Actualizando provincias")
    updateSelectInput(session, "province",choices = provinces())
  })
  observe({
    if (debug) print("Actualizando countries")
    updateSelectInput(session, "country",choices = countries())
  })
  incidf <-reactive({
    if (debug) print("Getting incidence")
    wait_for_values()
    if (input$country!=""){
      incidf<-get_cols(df=switch(names(which(input$datar==sources_str)), jh=ccovid19, ecdc=covid19, narr=narr_csv,#casal=acumula2,alv=alvaro,
                                 MSC=msc_csv, uk=ukcsv, catmc = cat_m,catmsc=cat_m, catzc=cat_z,catzsc=cat_z,usa=USA_csv,datasources[[names(which(class_sources==input$datar))]]@df),
                       source=input$datar, country=input$country,province=input$province)
      #preprocess
      #incidf[is.na(incidf),"I"]<-0
      total_cases_original_data<-sum(incidf$I)
      if (nrow(incidf)>0){
        incidf<-incidf[(incidf$dates>=as.Date(input$start_date)) &
                         (incidf$dates<=as.Date(input$end_date)),]
        incidfl<<-c(incidf_prepro(incidf=incidf,undetected=input$undetected, negatives=input$negatives,positives=input$positives,
                                  smooth_method=input$smooth_method,smooth_value=input$smooth_value,
                                  ignore=input$ignore,warnings=warnings),list(total_cases_original_data=total_cases_original_data,
                                                                              source=input$datar,country=input$country,province=input$province))
        incidfl$date_obtained=if (input$datar %in% class_sources) datasources[[names(which(class_sources==input$datar))]]@date_obtained else
          if (input$datar %in% c(CAT_mc_str,CAT_msc_str)) cat_m_date else if (input$datar %in% c(CAT_zc_str,CAT_zsc_str)) cat_z_date else Sys.Date()
        return(incidfl)
      }
      else {
        if (debug) print("Error getting incidence")
        incidfl$incidfl$df<-incidfl$incidfl$df[c(),]
        return (incidfl)} # Return the last incidfl
    }
  })
  wait_for_values<-reactive({return(c(input$datar,input$province,input$state))})
  
  si_distr<-reactive({
    if (debug)print("si_distr")
    distribution<-names(which(App_si_equivalences==input$si_distr))
    if (distribution %in% c("Du","Du Bs","Nishiura","Zhao","Ali","Ali early","Ali late","Ganyani Singapore","Ganyani Tianjin")) {
      distr<-switch(distribution,
                       `Du Bs`=si_distr_boot[,mcmc_limit],
                       Du=si_sample_mcmc_X$si_sample[,mcmc_limit],
                       Nishiura = si_sample_mcmc_N$si_sample[,mcmc_limit],
                       Zhao = si_sample_mcmc_Zhao$si_sample[,mcmc_limit],
                       Ali = si_sample_mcmc_Ali$si_sample[,mcmc_limit],
                       `Ali early` = Ali_early_distr[,mcmc_limit],
                       `Ali late` = Ali_late_distr[,mcmc_limit],
                       `Ganyani Singapore` = G_sing_distr[,mcmc_limit],
                       `Ganyani Tianjin` = G_tian_distr[,mcmc_limit])
      return(distr)
    }
    else if (distribution %in% c("Ganyani Singapore Pt","Ganyani Tianjin Pt","Li","Liu")){
      distr<-dgamma(seq(0,30),shape=switch(distribution,
                                                     `Ganyani Singapore Pt`=sing_m$kappa, `Ganyani Tianjin Pt`=tian_m$kappa,
                                                     `Li`=li_m$kappa,`Liu`=liu_m$kappa),#kappa
                              scale=switch(distribution, `Ganyani Singapore Pt`=sing_m$theta, `Ganyani Tianjin Pt`=tian_m$theta,
                                           `Li`=li_m$theta,`Liu`=liu_m$theta))#theta
    }
    else if (distribution %in% c("Du Point","Nishiura Pt","Zhao Pt", "Ali Pt")){
      distr<-dlnorm(seq(0,30),meanlog=switch(distribution,
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
      distr<-dtruncnorm(seq(0,30),a=1e-6, mean=switch(distribution,
                                                                `Ali early Pt`=ali_early_isolation$mean,
                                                                `Ali late Pt`=ali_late_isolation$mean),
                                  sd=switch(distribution,
                                            `Ali early Pt`=ali_early_isolation$sd,
                                            `Ali late Pt`=ali_late_isolation$sd))
    }
    return(distr/sum(distr))
  })
  config<-reactive({
    distribution<-
    if (debug)print("configurando ventana")
    t_start <- seq(2, nrow(incidf()$df)-input$window-1)
    return(make_config(list(
      si_distr=if(names(which(App_si_equivalences==input$si_distr)) %in% c("Du","Du Bs","Nishiura","Zhao","Ali","Ali early","Ali late","Ganyani Singapore","Ganyani Tianjin")) si_distr() else NA,
      #si_parametric_distr =if(input$si_distr %in% c("Lognormal + MCMC (Du et al) (much slower)","Lognormal + MCMC (Nisiura et al) (much slower)")) "L" else NA,
      #mcmc_control=if(input$si_distr=="Lognormal + MCMC (Du et al) (much slower)") make_mcmc_control(init_pars=c(4,2)) else NA,
      t_start =t_start,
      t_end = t_start + input$window -1)))
    
  })
  estR0<-reactive({
    if (debug)print("calculando")
    #if (input$datar!=incidf()$country | input$country!= incidf()$country | input$)  
    wait_for_values()
    static_incidf<-incidf()
    if (input$method=="A Cori et al (2013)") {
      if(debug) print(incidf())
      return(calculate_R(incidf=static_incidf$df, prepro = FALSE,window = input$window,distribution =
                    names(which(input$si_distr ==App_si_equivalences)),return_object=TRUE))
    }
    else if (input$method %in% c("Wallinga and Teunis (2004)","Bettencourt and Ribeiro (2008)")) {
      if (input$si_distr %in% c(Ln_du_str,Ln_ni_str,Ln_zh_str,Ln_ali_str,Tn_ali_early_str,Tn_ali_late_str,G_tian_str,G_sing_str)){
        rtc<-calculate_Rc(incidf=static_incidf$df,method = (if (input$method=="Wallinga and Teunis (2004)") "TD" else "SB"),
        distribution=names(which(input$si_distr ==App_si_equivalences)))
        return(if(is.na(rtc)) {c("Error. Estimating too early in the epidemic")} else rtc)
        }
    else if (input$method %in% c("Empirical R")) {
      if(length(static_incidf) & (nrow(static_incidf)>0))
        return(calculate_Empirical_R(incidf=static_incidf))
      }
    }
  })
  #infectividad<-reactive({
  #  if (debug)print("infectividad")
  #  overall_infectivity(incidf()$df,si_distr())
  #})
  output$plot_smooth <- renderCachedPlot({
    if (debug)print("plot_smooth")
    inci_df<-incidf()
    if (!is.null(incidf()) & (length(incidf)>0) & (nrow(inci_df$df)>0)){
      
      plot(y=(inci_df$df)$I_pre_smooth,x=(inci_df$df)$dates, xlab="Dates", ylab="Daily incidence",cex.lab=1.7, cex.axis=1.8,cex=1.7)
      lines(y=(inci_df$df)$I_pre_smooth,x=(inci_df$df)$dates,col="black",type="b",lty=1,cex=1.6)
      if(input$smooth_value>0) lines(inci_df$smooth)
      points(y=(inci_df$df)$I,x=(inci_df$df)$dates,col="darkgreen",pch=4,cex=1.2)
      lines(y=(inci_df$df)$I,x=(inci_df$df)$dates,col="darkgreen",type="b",lty=2,cex=1.2)
      legend('topleft',legend=c("Before smoothing","After smoothing"), col=c("black", "darkgreen"), lty=1:2, pch=c(1,4))
      grid()
    } else {print("No data")}
  }, cache="app", cacheKeyExpr={list(input$datar,input$country,input$province,input$start_date,input$end_date,input$ignore,input$undetected,input$smooth_method,input$smooth_value,input$negatives,input$positives,Sys.Date())}
  #,sizePolicy = sizeGrowthRatio(width =1200, height = 800, growthRate = 1.1)
  )
  cioidf<-reactive({
    if(debug)print("cioidf")
    static_incidf<-incidf()$df
    static_si_distr<-si_distr()
    si_distr_oi<-if (length(dim(static_si_distr))==2) apply(static_si_distr,MARGIN=1,FUN=mean) else static_si_distr
    return(rbind(data.frame(dates=static_incidf$dates,Ind=overall_infectivity(static_incidf[,c("dates","I")],si_distr_oi)*10, Series="Overall Infectivity * 10"),
               calculate_Ind(orig_incidf=static_incidf, n_list=c(7,14),population=1,pop_factor=1)[,c("dates","Ind","n")] %>% dplyr::mutate(Series=paste("CIR_",n,sep="")) %>%
                 dplyr::select(-n)))
  })
  output$cioidf<-reactive({return(cioidf())})
  output$cioi <- renderPlot({
    if(debug)print("cioi")
    static_cioi<-cioidf()
    if((length(cioidf)>0)){
      ggplot(static_cioi, aes(x = dates, y = as.numeric(`Ind`),color=Series)) +
        geom_line(size=1.3) +
        #geom_line(aes(y = Mean.R.,color=Regions),size=1.3,alpha=0.8) +
        #ggtitle("Estimated R")+scale_fill_discrete(name = "Country", labels = EU) +
        ggtitle(paste("Prevalence indicators"))+ 
        ggplot2::theme(plot.title = element_text(size=20), legend.title = element_text(size = 18), legend.text = element_text(size = 18),
          #line = element_line(size=18),
          axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20,angle = 45, vjust = 0.5),
          axis.text.y = element_text(size = 20),
          legend.position="bottom", legend.box = "horizontal")+
        labs(x=element_blank(),y=element_blank())    #ylim(c(0.2,4.8)) + geom_line(aes(y = Mean.R.),color="black",size=0.1,alpha=0.8) +
    }
  })#, cache="app", cacheKeyExpr={list(input$datar,input$country,input$province,input$start_date,input$end_date,input$ignore,input$undetected, input$smooth_method, input$smooth_value, input$negatives, input$positives, input$Sys.Date(), input$si_distr, input$window)}
  #)
  output$miplot <- renderCachedPlot({
    static_estR0<-estR0()
    if (input$method %in% c("Wallinga and Teunis (2004)","Bettencourt and Ribeiro (2008)")){
      if (debug)print("plot wt")
      update_geom_defaults("line", list(size = 1.1))
      if(static_estR0==c("Error. Estimating too early in the epidemic")) {plot(c(0,1),c(0,1));text(0.5,0.5,"Too few cases in initial incidences\n",cex=2)}
      else ggplot(data=dtable() %>% tibble::rownames_to_column("dates") ,aes(x=as.Date(dates),y=R))+geom_point()+geom_line()+
        geom_line(aes(y=conf.int.lower))+geom_line(aes(y=conf.int.upper))+
        geom_ribbon(aes(ymin = `conf.int.lower`, ymax = `conf.int.upper`), fill = "grey70",alpha=0.3)+
        ggplot2::theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20,angle=45),
                       axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 20))
    }
    else {update_geom_defaults("line", list(size = 1.7));
      if (debug)print("plot r0")
      if(length(static_estR0)>0){
        plot(static_estR0,what="R") + scale_x_date(date_breaks = "2 month",date_labels = "%Y-%b-%d")+ ggplot2::theme(
        plot.title = element_text(size=24), legend.text = element_text(size = 24),
        #line = element_line(size=18),
        axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 22,angle = 45, vjust = 0.5),
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 24))}
      #estimate_R_plots(list(estR0(),wate()),options_R = list(col = c("blue", "red")),what="R")
    }
  }, cache="app", cacheKeyExpr=list(input$datar,input$country,input$province,input$start_date,input$end_date,input$ignore,input$undetected,input$smooth_method,input$smooth_value,input$negatives,input$positives,input$method,input$window,input$si_distr,Sys.Date())
  #,sizePolicy = sizeGrowthRatio(width =900, height = 700, growthRate = 1.1)
  )
  output$distro <- renderCachedPlot({
    if (debug)print("plot distro")
    if(input$method %in% c("Wallinga and Teunis (2004)","Bettencourt and Ribeiro (2008)"))
      plot(calculate_R(incidf=incidf()$df, prepro = FALSE,window = input$window,distribution =
                                            (names(which(input$si_distr ==App_si_equivalences))),return_object=TRUE),
           what="SI") + ggplot2::theme(
             axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24),
             axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 24))
    else{  
      plot(estR0(),what="SI") + ggplot2::theme(
        axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 24))
      }
  },cache="app", cacheKeyExpr = list(input$si_distr) )
  output$country_province <- reactive({
    return (paste(input$datar,"-",input$country,"-",input$province))
  })
  output$USA <- renderCachedPlot({
    if (is.na(USAs_csv[1])) {
      USAs_csv<<-read.csv(USAs_path,encoding="UTF-8", stringsAsFactors = TRUE)}
    if (nrow(USAs_csv)){
      USA<-lapply(levels(USAs_csv$state),FUN = function(a_state){
        df<-fdf_to_inci(USAs_csv %>% filter(state==a_state) %>% rename(dates="date",) %>% dplyr::group_by(dates) %>%
                          dplyr::summarize(Cum= sum(cases)) ,"%Y-%m-%d")
        df<-incidf_prepro(incidf=df,undetected=0.4, negatives="Distribute backwards before smoothing",
                          smooth_method="Lowess",smooth_value=6,
                          ignore=20,warnings="")$df %>% dplyr::select (c("dates","I"))
        
        if (nrow(df)<10) return (c(a_state,NA))
        Rs<-(estimate_R(df,method="si_from_sample",si_sample=si_sample_mcmc_X$si_sample[,1])$R)
        if (diff(unlist(Rs[nrow(Rs),c("Quantile.0.025(R)","Quantile.0.975(R)")]))<1){
          return (c(a_state,Rs[nrow(Rs),"Mean(R)"]))
        }
        else return (c(a_state,NA))
      })
    }
    USA_df<-as.data.frame(do.call(rbind,USA),stringsAsFactors = FALSE) %>% rename (state=1,R=2) %>% mutate (R=as.numeric(R))
    plot_usmap(data=USA_df,values="R",labels=TRUE) + scale_fill_continuous(low="white",high="red")
  },cache="app",cacheKeyExpr = {format(Sys.time(), "%Y-%m-%d")})
  output$report <- reactive({
    if(debug)print("report")
    #inci_df<-incidf()
    return(paste(incidf()$total_cases_original_data,"COVID-19 cases were reported in",if ((input$province!="* All *") & (input$datar!=ecdc_str)) paste (input$province,",",sep="") else "", input$country))
  })
  output$warnings<- reactive(incidf()$warnings)
  output$data <- renderDataTable(incidf()$df)
  dtable <- reactive({
    if (debug)print("dtable")
    if(input$method=="Wallinga and Teunis (2004)"){
      if (estR0()==c("Error. Estimating too early in the epidemic")) return(data.frame(c("Error","Estimating too early in the epidemic")))
      else return(data.frame(R=estR0()$estimates$TD$R, conf.int.lower=estR0()$estimates$TD$conf.int$lower, conf.int.upper=estR0()$estimates$TD$conf.int$upper))
    }
    else if (input$method=="A Cori et al (2013)"){
      tabla<-estR0()$R
      tabla[,"date_start"]<-as.Date(estR0()$dates[tabla[,"t_start"]],origin="1970-01-01")
      tabla[,"date_end"]<-as.Date(estR0()$dates[tabla[,"t_end"]],origin="1970-01-01")
      return(tabla[,c("date_start","date_end",colnames(tabla)[3:11])])
    }
    else if (input$method=="Empirical_R"){
      tabla<-estR0()
      return(tabla)
    }
  })
  output$tabla <- renderDataTable({
    if (input$method=="A Cori et al (2013)") dtable() %>% datatable() %>% formatRound(c(3:11),digits=3) else
      dtable() %>% datatable() %>% formatRound(c(1:3),digits=3)
  })
  output$source <- reactive({
    data_source<-NA
    if (input$datar==jhu_str) data_source<-jhu_source
    else if (input$datar==ecdc_str) data_source<-(paste("European Centre for Disease Prevention and Control, https://data.europa.eu/euodp/en/data/dataset/covid-19-coronavirus-data"))
    else if (input$datar==narr_str) data_source<-(narr_source)
    #else if (input$datar==casal_str) data_source<-(paste("Rubén F Casal, https://rubenfcasal.github.io/COVID-19/, Ministerio de Sanidad y Consumo, Instituto de Salud Carlos III. (file: acumula2.RData"))
    else if (input$datar==uk_str) data_source<-(uk_source)
    #else if (input$datar==alvaro_str) data_source<-(alvaro_source)
    else if (input$datar==USA_str) data_source<-(USA_source)
    else if (input$datar %in%  c(CAT_mc_str,CAT_msc_str)){data_source<-(CAT_m_source)}
    else if (input$datar %in%  c(CAT_zc_str,CAT_zsc_str)){data_source<-(CAT_z_source)}
    else if (input$datar %in% class_sources){
      source<-names(which(class_sources==input$datar))
      data_source<-(datasources[[source]]@source)
    }
    return (paste (data_source,". Date obtained: ",strftime(incidf()$date_obtained,"%Y-%m-%d")))
  })
  
    output$downloadcioi <- downloadHandler(
    filename = function() {
      paste(input$country, "_prev.csv", sep = "")
    },
    content = function(file) {
      write.csv(cioidf(), file, row.names = FALSE)
    })
  output$downloadRt <- downloadHandler(
    filename = function() {
      paste(input$country, "_rt.csv", sep = "")
    },
    content = function(file) {
      write.csv(dtable(), file, row.names = FALSE)
    })
  output$summary <- reactive({
    #colnames(estR0()$R)#t_start,t_end,Mean(R),Std(R),Quantile.0.025(R),Quantile.0.05(R),Quantile.0.25(R),Median(R),Quantile.0.75(R),Quantile.0.95(R),Quantile.0.975(R)
    if(input$method=="Wallinga and Teunis (2004)"){
      lastdata<-tail(dtable(),1)
      return(paste("Last effective R",row.names(lastdata[1,]), "=>",round(lastdata[3,"R"],3),"(",round(lastdata[3,"conf.int.lower"],3),
                   "-",round(lastdata[3,"conf.int.upper"],3),")"))
    }
    else{
      lastdata<-tail(dtable(),3)
      return(paste("Last effective R, between",lastdata[3,"date_start"],"-",lastdata[3,"date_end"], "=>",round(lastdata[3,"Median(R)"],3),"(",round(lastdata[3,"Quantile.0.025(R)"],3),
                   "-",round(lastdata[3,"Quantile.0.975(R)"],3),")"))
    }
  })
  provinces_summary<- eventReactive(input$provinces_btn, {
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
