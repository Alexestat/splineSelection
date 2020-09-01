####################################################################################
###### Obtem melhor ajuste para os dados do covid para um modelo de
###### regressão com splines cúbicos naturais, com número e posições
###### variáveis para os knots, utilizando métodos de seleção de
###### modelos e detecção de pontos de mudança
######
###### Autores: Alex Rodrigo, Florencia Leonardi e Magno Severino
###### Instituto de Matemática e Estatística - Universidade de São Paulo
###### Data: 03/06/2020
####################################################################################


library(lubridate)
library(dplyr)
library(ggplot2)
library(splines)
library(stringr)
library(parallel)


##### INÍCIO FUNÇÕES

##### Best fit for a fixed number of knots  ######

bestfitcovidfix = function(covdata,nknots,delta,degree,
                           b.knots=c(min(covdata$index),max(covdata$index)),
                           bspline=TRUE){
  
  if( (b.knots[2] - b.knots[1]) / delta < (nknots+1) ){
    print("nknots is too large with respect to delta: returning NA for this number of nknots\n")
    return(list(knots=NA,min.eqm=NA,model=NULL))
  }
  
  posknots <- covdata$index[!(covdata$index %in% c(min(covdata$index), max(covdata$index)))] %>% 
    as.integer() %>% sort()

  posknots <- posknots[posknots >= b.knots[1] & posknots <= b.knots[2]]
  
  if( nknots > 0){
   
      allposknots = rbind(min(covdata$index),
                         combn(posknots,nknots),
                         max(covdata$index))
      
      allposknots = allposknots[,apply(apply(allposknots,FUN=dist,2),
                                       FUN=min,2)>delta]
      
      if( ncol(allposknots) == 0){
          print("allposknots with 0 columns: decrease nknots or delta, or aumentgs sample size\n")
          return(list(knots=NA,min.eqm=NA,model=NULL))
      }
      
      EQM = NA

      ret = mclapply(1:ncol(allposknots), EQM_cal, bspline, covdata,
                     allposknots, nknots, degree, b.knots,
                     mc.cores=detectCores())
      EQM = simplify2array(ret)      

      ## for(i in 1:ncol(allposknots)){
      ##     if ( isTRUE(bspline) ) 
      ##         fit=lm(cases ~ bs(index,knots=c(allposknots[2:(nknots+1),i]),
      ##                           degree=degree), data = covdata )
      ##     else 
      ##         fit=lm(cases ~ ns(index,knots=c(allposknots[2:(nknots+1),i]),
      ##                           Boundary.knots=b.knots), data = covdata )
      ##     EQM[i] = mean(fit$residuals^2)
      ## }

      chosenknots = allposknots[2:(nknots+1),which.min(EQM)]
      
      if ( isTRUE(bspline) )
          bestfit=lm(cases ~ bs(index,knots = c(chosenknots),degree=degree),
                     data = covdata )
      else
          bestfit=lm(cases ~ ns(index,knots = c(chosenknots),
                                Boundary.knots=b.knots), data = covdata ) 
  }
  else {
      chosenknots = NA
      if ( isTRUE(bspline) )
          bestfit=lm(cases ~ bs(index,knots=NULL,degree=degree),
                     data = covdata )
      else
          bestfit=lm(cases ~ ns(index,knots=NULL,Boundary.knots=b.knots),
                     data = covdata ) 
      EQM = mean(bestfit$residuals^2)
  } 
    
    return(list(knots=chosenknots,min.eqm=min(EQM),model=bestfit))
    
}

#### Best fit for all number of knots (until a limit nlimknot) ####

bestfitcovid = function(covdata,nlimknot,lambda,delta,degree,bspline=TRUE,
                        b.knots=NULL,nome='nome'){
    
    fixknots = matrix(0,nlimknot+1,nlimknot+1)
    if ( is.null(b.knots) ) {
        if ( isTRUE(bspline) ) {
            first.knot = sort(covdata$index)[round(length(covdata$index) *
                                                   0.25)]
            last.knot = max(covdata$index)
            b.knots = c(first.knot, last.knot)
            print(b.knots)
        }
        else {
            first.knot = sort(covdata$index)[round(length(covdata$index) *
                                                   0.10)] 
            last.knot = sort(covdata$index)[round(length(covdata$index) *
                                                  0.95)]
            b.knots = c(first.knot, last.knot)
        }
    }
    
    EQMpen = vector()  
    for(i in 1:(nlimknot+1)){
        a = bestfitcovidfix(covdata,i-1,delta,degree,b.knots,bspline)
        EQMpen[i] = a$min.eqm + lambda*i
        fixknots[i,1:length(a$knots)] = a$knots
    }
    
    min.nknot = which.min(EQMpen)
    if ( min.nknot == 1 ) {bestknots = NULL
    }else {bestknots = fixknots[which.min(EQMpen),1:(min.nknot-1)]}
    
    if ( isTRUE(bspline) ){
        bestfitall=lm(cases ~ bs(index,knots = c(bestknots),degree=degree), data = covdata )
    }else{
        bestfitall=lm(cases ~ ns(index,knots = c(bestknots),Boundary.knots=b.knots), data = covdata ) }
    
    return(list(bestknots, bestfitall))
}

## Calcula EQM da i-ésima combinação 
EQM_cal <- function(i, bspline, covdata, allposknots, nknots, degree,
                    b.knots) {
    if ( isTRUE(bspline) ) 
              fit=lm(cases ~ bs(index,knots=c(allposknots[2:(nknots+1),i]),
                                degree=degree), data = covdata )
          else 
              fit=lm(cases ~ ns(index,knots=c(allposknots[2:(nknots+1),i]),
                                Boundary.knots=b.knots), data = covdata )
    mean(fit$residuals^2)
}

####   Passar filtro suavizando os dados  ####
suavizar_dados <- function(covdata, janela=7){
  
  n <- nrow(covdata)
  new_data_cases <- c()
  new_data_deaths <- c()
  
  for(i in 1:(n-janela+1)){
    mean_cases <- mean(covdata$cases[i:(i+janela-1)], na.rm = T)
    mean_deaths <- mean(covdata$deaths[i:(i+janela-1)], na.rm = T)
    new_data_cases <- c(new_data_cases, ifelse(is.nan(mean_cases), 0, mean_cases))
    new_data_deaths <- c(new_data_deaths, ifelse(is.nan(mean_deaths), 0, mean_deaths))
  }
  
  new_data_cases <- append(new_data_cases, rep(0, janela-1))
  new_data_deaths <- append(new_data_deaths, rep(0, janela-1))
  
  covdata <- covdata %>% mutate(cases_smooth = new_data_cases, deaths_smooth = new_data_deaths)
  
  return(covdata)
}


##### FIM FUNÇÕES


##### INÍCIO ANÁLISE DE DADOS DE PAÍSES

data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                 na.strings = "", fileEncoding = "UTF-8-BOM")

data <- data %>% 
  mutate(date = dmy(dateRep))

min_date_paises = min(data$date)

data <- data %>% 
  mutate(index = date - min_date_paises) %>% 
  select(countriesAndTerritories, date, index, cases, deaths)


# obter pontos de mudanca paises --------------------------------------
dicionario_paises <- data.frame(country = c('Brazil', 'Argentina', 'Chile', 'Mexico', 'Colombia', 
                                            'Peru', 'Ecuador', 'Bolivia', 'Guatemala', 'Uruguay', 
                                            'United_States_of_America', 'Spain', 'Italy', 
                                            'United_Kingdom', 'Germany', 'France', 'Austria', 
                                            'Bangladesh', 'Belarus', 'Belgium', 'Canada', 'China', 
                                            'Denmark', 'Dominican_Republic', 'Egypt', 'India', 
                                            'Indonesia', 'Iran', 'Ireland', 'Israel', 'Japan', 
                                            'Kuwait', 'Netherlands', 'Pakistan', 'Philippines', 
                                            'Poland', 'Portugal', 'Qatar', 'Romania', 'Russia', 
                                            'Saudi_Arabia', 'Serbia', 'Singapore', 'South_Africa', 
                                            'South_Korea', 'Sweden', 'Switzerland', 'Turkey', 
                                            'Ukraine', 'United_Arab_Emirates'),
                                pais = c("Brasil", "Argentina", "Chile", "México",
                                         "Colômbia", "Peru", "Equador", "Bolívia",
                                         "Guatemala", "Uruguai",
                                         "Estados Unidos", "Espanha", "Itália",
                                         "Reino Unido", "Alemanha", "França",
                                         'Áustria', 'Bangladesh', 'Belarus', 'Bélgica', 
                                         'Canadá', 'China', 'Dinamarca', 'República Dominicada', 
                                         'Egito', 'Índia', 'Indonésia', 'Irã', 'Irlanda', 
                                         'Israel', 'Japão', 'Kuwait', 'Holanda', 'Paquistão', 
                                         'Filipinas', 'Polônia', 'Portugal', 'Qatar', 'Romênia', 
                                         'Rússia', 'Arábia Saudita', 'Sérvia', 'Singapura', 
                                         'África do Sul', 'Coreia do Sul', 'Suécia', 'Suiça', 
                                         'Turquia', 'Ucrânia', 'Emirados Árabes Unidos'),
                                stringsAsFactors = FALSE)

df_cases <- data.frame(country = character(), 
                       date = as.Date(character()), 
                       cases = numeric(), 
                       deaths = numeric(), 
                       pred = numeric(),
                       ic_max = numeric(),
                       ic_min = numeric())

df_pred <- data.frame(country = character(), 
                      date = as.Date(character()), 
                      pred = numeric(),
                      ic_max = numeric(),
                      ic_min = numeric())

df_changepoints <- data.frame(country = character(),
                              points = character(),
                              dates = character())


inicio <- Sys.time()
for(i in 1:nrow(dicionario_paises)){ 
  print(paste(dicionario_paises$country[i], ": ", i))
  covdata <- data[data$countriesAndTerritories == dicionario_paises$country[i],]
  first.obs <- max(which(covdata$cases>0))
  covdata <- covdata[1:first.obs,]
  
  period <- seq(min(covdata$date), max(covdata$date), by="1 days")
  df <- data.frame(date = period,
                   inputed = !(period %in% covdata$date)) #marca aqueles que foram imputados
  
  covdata <- covdata %>% full_join(df, by = c("date"="date")) %>%
    ungroup() %>%
    mutate(state = covdata$state[1],
           countriesAndTerritories = covdata$countriesAndTerritories[1],
           index = ifelse(is.na(index), date - min_date_paises, index)) %>%
    arrange(desc(date))
  
  covdata$cases_raw <- covdata$cases
  covdata$deaths_raw <- covdata$deaths

  covdata <- suavizar_dados(covdata, janela=7)
 
  first.obs_new <- min(max(which(covdata$cases_smooth > 0.01*max(covdata$cases_smooth,na.rm=TRUE))),75)
  
  covdata <- covdata[1:first.obs_new,] # %>% 
 
  covdata$cases <- covdata$cases_smooth
 
  if(sum(is.nan(covdata$cases))>0)
    next
  
  resultados_smooth <- bestfitcovid(covdata, nlimknot = 4, lambda = 1, 
                                    delta = 5, degree = 3, bspline=FALSE,
                                    b.knots=c(min(covdata$index)+5,max(covdata$index)-5),
                                    nome = dicionario_paises$country[i])
  
  bestknots_new <- resultados_smooth[[1]]
  bestfitall_new <- resultados_smooth[[2]]
  
  covdata <- covdata %>% full_join(data.frame(index = covdata$index, 
                                                      pred = predict(bestfitall_new)),
                                           by = c("index"="index"))
  
  previsao_new <- seq(max(covdata$index)+1, max(covdata$index)+7)
  
  futuro_new <- data.frame(index = previsao_new,
                           pred = predict(bestfitall_new,data.frame(index=previsao_new)))
  
  newx_new = previsao_new
  
  newy_new = predict(bestfitall_new,data.frame(index=newx_new),interval="prediction")
  preds.data_new = predict(bestfitall_new,data.frame(index=covdata$index),interval="prediction")
  
  se.bands1_new=newy_new[,2:3]
  se.bands2_new=preds.data_new[,2:3]
  
  df_cases <- rbind(df_cases, 
                    cbind(covdata[, c("countriesAndTerritories", "date", 
                                          "cases", "deaths", "pred", "cases_raw", "deaths_raw")],
                          se.bands2_new))
  
  df_changepoints <- rbind(df_changepoints,
                           data.frame(country = dicionario_paises$country[i],
                                      dates = bestknots_new + min(data$date)))
  
  df_pred <- rbind(df_pred, 
                   data.frame(country = dicionario_paises$country[i], 
                              date = futuro_new$index + min(data$date),
                              pred = futuro_new$pred,
                              ic_max = se.bands1_new[,2],
                              ic_min = se.bands1_new[,1]))
}
fim <- Sys.time()
fim-inicio
names(df_cases) <- c("country", "date", "cases", "deaths", "pred", "cases_raw", "deaths_raw", 
                     "ic_min", "ic_max")

dicionario_paises <- dicionario_paises %>% filter(country %in% unique(df_cases$country))

# salvar rdata ------------------------------------------------------------

save(df_changepoints, df_cases, df_pred, dicionario_paises,
     file = paste("paises",max(data$date),"_suavizados_cubicos.RData", sep=""))


##### FIM ANÁLISE DE DADOS DE PAÍSES


##### INÍCIO ANÁLISE DE DADOS DO BRASIL

library(tidyverse)
dfbrasil <- read_csv("https://raw.githubusercontent.com/wcota/covid19br/master/cases-brazil-cities-time_changesOnly.csv")
names(dfbrasil) <- c('date', 'country', 'state', 'city', 'ibge_id', 'new_deaths', 'deaths', 
                     'new_cases', 'total_cases', 'deaths_per_100k_inhabitants', 
                     'total_cases_per_100k_inhabitants', 'deaths_by_total_cases', 'source')

cidades_1000_casos <- dfbrasil %>% filter(city!="TOTAL") %>%  group_by(city) %>%
  summarise(cases = sum(new_cases), n = n()) %>% arrange(desc(cases)) %>%
  filter(cases >= 1000 & n > 40) %>% select(city) %>% pull()
capitais <- c("Rio Branco/AC",
              "Maceió/AL",
              "Macapá/AP",
              "Manaus/AM",
              "Salvador/BA",
              "Fortaleza/CE",
              "Brasília/DF",
              "Vitória/ES",
              "Goiânia/GO",
              "São Luís/MA",
              "Cuiabá/MT",
              "Campo Grande/MS",
              "Belo Horizonte/MG",
              "Belém/PA",
              "João Pessoa/PB",
              "Curitiba/PR",
              "Recife/PE",
              "Teresina/PI",
              "Rio de Janeiro/RJ",
              "Natal/RN",
              "Porto Alegre/RS",
              "Porto Velho/RO",
              "Boa Vista/RR",
              "Florianópolis/SC",
              "São Paulo/SP",
              "Aracaju/SE",
              "Palmas/TO"#,
              # "Blumenau/SC"
)

cidades <- sort(unique(c(capitais, 
                         cidades_1000_casos[!str_detect(cidades_1000_casos, "CASO SEM LOCALIZAÇÃO DEFINIDA")])))

min_date_cidades = min(dfbrasil$date)
brasil <- dfbrasil %>% 
  filter(city %in% cidades,
         new_cases>=0) %>% 
  mutate(index = date -min_date_cidades,
         cases = new_cases) %>% 
  select(city, state, date, total_cases, cases, new_deaths, index) %>% 
  mutate(deaths=new_deaths) %>% 
  select(city, state, date, total_cases, cases, deaths, index) %>% 
  group_by(city) %>% 
  arrange(city, desc(date))

# obter pontos de mudanca cidades brasil --------------------------------------
df_cases_brasil <- data.frame(city = character(),
                              date = as.Date(character()),
                              cases = numeric(),
                              deaths = numeric(),
                              pred = numeric(),
                              ic_max = numeric(),
                              ic_min = numeric())

df_pred_brasil <- data.frame(city = character(),
                             date = as.Date(character()),
                             pred = numeric(),
                             ic_max = numeric(),
                             ic_min = numeric(),
                             se = numeric())

df_changepoints_brasil <- data.frame(city = character(),
                                     points = character(),
                                     dates = character())

brasil <- brasil %>% filter(city !="CASO SEM LOCALIZAÇÃO DEFINIDA/CE") 

inicio <- Sys.time()
for(cty in sort(unique(brasil$city))){
  print(cty)

  covdata <- brasil[brasil$city == cty,] %>% as.data.frame()
  first.obs <- max(which(covdata$cases>0))
  covdata <- covdata[1:first.obs,]

  period <- seq(min(covdata$date), max(covdata$date), by="1 days")
  df <- data.frame(date = period,
                   inputed = !(period %in% covdata$date)) 

  covdata <- covdata %>% full_join(df, by = c("date"="date")) %>%
    ungroup() %>%
    mutate(state = covdata$state[1],
           city = covdata$city[1],
           index = ifelse(is.na(index), date-min_date_cidades, index)) %>%
    arrange(desc(date))
 
  covdata$cases_raw <- covdata$cases
  covdata$deaths_raw <- covdata$deaths
  
  covdata <- suavizar_dados(covdata, janela=7)
  
  first.obs_new <- min(max(which(covdata$cases_smooth > 0.01*max(covdata$cases_smooth,na.rm=TRUE))),75)
  
  covdata <- covdata[1:first.obs_new,] 
 
  covdata$cases <- covdata$cases_smooth

  if(sum(is.nan(covdata$cases))>0)
    next

  new_covdata <- covdata
  
  resultados_smooth <- bestfitcovid(new_covdata, nlimknot = 4, lambda = 1,
                                    delta = 5, degree = 3, bspline=FALSE,
                                    b.knots=c(min(covdata$index)+5,max(covdata$index)-5),
                                    nome = dicionario_paises$country[i])

  bestknots_new <- resultados_smooth[[1]]
  bestfitall_new <- resultados_smooth[[2]]

  new_covdata <- new_covdata %>% full_join(data.frame(index = new_covdata$index,
                                                      pred = predict(bestfitall_new)),
                                           by = c("index"="index"))

  previsao_new <- seq(max(new_covdata$index)+1, max(new_covdata$index)+7)

  futuro_new <- data.frame(index = previsao_new,
                           pred = predict(bestfitall_new,data.frame(index=previsao_new)))

  newx_new = previsao_new

  newy_new = predict(bestfitall_new,data.frame(index=newx_new),interval="prediction")
  preds.data_new = predict(bestfitall_new,data.frame(index=new_covdata$index),interval="prediction")

  se.bands1_new=newy_new[,2:3]
  se.bands2_new=preds.data_new[,2:3]

  df_cases_brasil <- rbind(df_cases_brasil,
                           cbind(new_covdata[, c("city", "date", "cases", "deaths", "pred",
                                                 "cases_raw", "deaths_raw")],
                                 se.bands2_new))


  if(!is.null(bestknots_new)){
    df_changepoints_brasil <- rbind(df_changepoints_brasil,
                                    data.frame(city = cty,
                                               points = bestknots_new,
                                               dates = bestknots_new + min(brasil$date)))
  }

  df_pred_brasil <- rbind(df_pred_brasil,
                          data.frame(city = cty,
                                     date = futuro_new$index + min(brasil$date),
                                     pred = futuro_new$pred,
                                     ic_max = se.bands1_new[,2],
                                     ic_min = se.bands1_new[,1]))
}
fim <- Sys.time()
fim-inicio

# salvar dados ------------------------------------------------------------

names(df_cases_brasil) <- c("city", "date", "cases", "deaths", "pred", "cases_raw", "deaths_raw", 
                            "ic_min", "ic_max")

save(df_changepoints_brasil, df_cases_brasil, df_pred_brasil,
     file = paste("cidades",max(brasil$date),"_suavizados_cubicos.RData", sep=""))



# obter dados por estados -------------------------------------------------
min_date_estados <- min(dfbrasil$date)
estados <- dfbrasil %>%
  filter(state != "TOTAL") %>%
  group_by(state, date) %>%
  summarise(cases = sum(new_cases),
            total_cases = sum(total_cases),
            new_deaths = sum(new_deaths)) %>%
  mutate(index = date - min_date_estados) %>%
  group_by(state) %>%
  arrange(state, desc(date)) %>%
  mutate(deaths = new_deaths) %>%
  select(state, date, cases, deaths, index)

# obter pontos de mudanca estados ------------------------------------------
df_cases_estados <- data.frame(state = character(),
                               date = as.Date(character()),
                               cases = numeric(),
                               deaths = numeric(),
                               pred = numeric(),
                               ic_max = numeric(),
                               ic_min = numeric())

df_pred_estados <- data.frame(state = character(),
                              date = as.Date(character()),
                              pred = numeric(),
                              ic_max = numeric(),
                              ic_min = numeric())

df_changepoints_estados <- data.frame(state = character(),
                                      points = character(),
                                      dates = character())

inicio <- Sys.time()
for(est in unique(estados$state)){
  print(est)

  covdata <- estados[estados$state == est,] %>% as.data.frame()
  first.obs <- max(which(covdata$cases>0))
  covdata <- covdata[1:first.obs,]
  
  period <- seq(min(covdata$date), max(covdata$date), by="1 days")
  df <- data.frame(date = period,
                   inputed = !(period %in% covdata$date)) 

  covdata <- covdata %>% full_join(df, by = c("date"="date")) %>%
    ungroup() %>%
    mutate(state = covdata$state[1],
           city = covdata$city[1],
           index = ifelse(is.na(index), date-min_date_estados, index)) %>%
    arrange(desc(date))
  
  covdata$cases_raw <- covdata$cases
  covdata$deaths_raw <- covdata$deaths
 
  covdata <- suavizar_dados(covdata, janela=7)
  
  first.obs_new <- min(max(which(covdata$cases_smooth > 0.01*max(covdata$cases_smooth,na.rm=TRUE))),75)
  
  covdata <- covdata[1:first.obs_new,] 
  
  covdata$cases <- covdata$cases_smooth
  
  new_covdata <- covdata
  
  resultados_smooth <- bestfitcovid(new_covdata, nlimknot = 4, lambda = 1,
                                    delta = 5, degree = 3, bspline=FALSE,
                                    b.knots=c(min(covdata$index)+5,max(covdata$index)-5),
                                    nome = dicionario_paises$country[i])

  bestknots_new <- resultados_smooth[[1]]
  bestfitall_new <- resultados_smooth[[2]]

  new_covdata <- new_covdata %>% full_join(data.frame(index = new_covdata$index,
                                                      pred = predict(bestfitall_new)),
                                           by = c("index"="index"))

  previsao_new <- seq(max(new_covdata$index)+1, max(new_covdata$index)+7)
  futuro_new <- data.frame(index = previsao_new,
                           pred = predict(bestfitall_new,data.frame(index=previsao_new)))

  newx_new = previsao_new

  newy_new = predict(bestfitall_new,data.frame(index=newx_new),interval="prediction")
  preds.data_new = predict(bestfitall_new,data.frame(index=new_covdata$index),interval="prediction")

  se.bands1_new=newy_new[,2:3]
  se.bands2_new=preds.data_new[,2:3]


  df_cases_estados <- rbind(df_cases_estados,
                            cbind(new_covdata[, c("state", "date", "cases", "deaths", "pred", 
                                                  "cases_raw", "deaths_raw")],
                                  se.bands2_new))

  if(!is.null(bestknots_new)){
    df_changepoints_estados <- rbind(df_changepoints_estados,
                                     data.frame(state = est,
                                                points = bestknots_new,
                                                dates = bestknots_new + min_date_estados))
  }

  df_pred_estados <- rbind(df_pred_estados,
                           data.frame(state = est,
                                      date = futuro_new$index + min_date_estados,
                                      pred = futuro_new$pred,
                                      ic_min = se.bands1_new[,1],
                                      ic_max = se.bands1_new[,2]))
}
fim <- Sys.time()
fim-inicio

names(df_cases_estados) <- c("state", "date", "cases", "deaths", "pred", "cases_raw", "deaths_raw", 
                             "ic_min", "ic_max")

dicionario_estados <- data.frame(sigla = c("AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", 
                                           "MA", "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", 
                                           "RJ", "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO"),
                                 nome = c("Acre", "Alagoas", "Amazonas", "Amapá", "Bahia", 
                                          "Ceará", "Distrito Federal", "Espirito Santo", "Goiás", 
                                          "Maranhão", "Minas Gerais", "Mato Grosso do Sul", 
                                          "Mato Grosso", "Pará", "Paraíba", "Pernambuco", 
                                          "Piauí", "Paraná", "Rio de Janeiro", "Rio Grande do Norte",
                                          "Rondonia", "Roraima", "Rio Grande do Sul", "Santa Catarina", 
                                          "Sergipe", "São Paulo", "Tocantins"),
                                 stringsAsFactors = FALSE)

save(df_changepoints_estados, df_cases_estados, df_pred_estados, dicionario_estados,
     file = paste("estados",max(estados$date),"_suavizados_cubicos.RData", sep=""))


