library(lubridate)
library(dplyr)
library(ggplot2)
library(splines)
library(stringr)

data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
                 na.strings = "", fileEncoding = "UTF-8-BOM")

data <- data %>%
  mutate(date = dmy(dateRep))

min_date_paises = min(data$date)

data <- data %>%
  mutate(index = date - min_date_paises) %>%
  select(countriesAndTerritories, date, index, cases, deaths)

covid19data = data

usethis::use_data(covid19data,compress="xz")





                       
