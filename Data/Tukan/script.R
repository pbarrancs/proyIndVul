###########################


# This is the helper script to summarise the infrastructure data 
# of Mexico's public health system.
# Author: TUKAN MX
#Last update: 03/29/20


##########################


############
# Load libraries
###########

library(tidyverse)

setwd('') #Here goes the directory on which you're working on.


#Let's load the data

health_2018 <- read.csv('raw_data/Recursos_Salud_2018.csv', encoding = 'utf-8', stringsAsFactors = F)
health_2017 <- read.csv('raw_data/Recursos_Salud_2017.csv', encoding = 'utf-8', stringsAsFactors = F)
health_2016 <- read.csv('raw_data/Recursos_Salud_2016.csv', encoding = 'utf-8', stringsAsFactors = F)


#First let's merge the three datasets together 

health <- rbind(health_2018[,(names(health_2017) %in% names(health_2018))], 
                health_2017[,(names(health_2017) %in% names(health_2018))],
                health_2016[,(names(health_2017) %in% names(health_2018))])

#Now let's select the variables of interest

health <- health %>% 
  select(AÑO, SECTOR, CLAVE.INSTITUCIÓN, Clave.Estado, Nombre.Estado, Tipo.de.Establecimiento,
         X.Cuenta.con.Unidad.de.Cuidados.Intensivos., X.Cuenta.con.área.de.Hospitalización., TOTAL.CAMAS.AREA.HOSPITALIZACIÓN, 
         TOTAL.CAMAS.EN.OTRAS.AREAS..NO.CONSIDERA.HOSPITALIZACIÓN., Número.de.camas.de.aislados, Número.de.camas.en.área.de.Urgencias,
         Número.de.camas.en.la.Unidad.de.Cuidados.Intensivos..incluye.pediátricas.y.adulto., Número.de.camas.en.la.Unidad.de.Cuidados.Intermedios..incluye.pediátricas.y.adulto.,
         Total.médicos.generales.y.especialistas, Total.enfermeras.en.contacto.con.el.paciente, Total.médicos.en.formación)

names(health) <- c('YEAR', 'SECTOR', 'INSTITUTION', 'STATE_ID', 'STATE_NAME', 'TYPE', 'INTENSIVE_CARE', 'HOSPITAL_AREA', 'HOSPITAL_AREA_BEDS',
                   'OTHER_AREAS_BEDS', 'ISOLATION_BEDS', 'ER_BEDS', 'INTENSIVE_CARE_BEDS', 'MEDIUM_CARE_BEDS', 'DOCTORS', 'NURSES', 'DOCTORS_IN_STUDY')

#We fix inconsistencies
health$TYPE[health$TYPE %in% c('DE CONSULTA EXTERNA', 'CONSULTA EXTERNA')] <- 'CONSULTA EXTERNA'
health$TYPE[health$TYPE %in% c('DE ASISTENCIA SOCIAL', 'ASISTENCIA SOCIAL')] <- 'ASISTENCIA SOCIAL'
health$TYPE[health$TYPE %in% c('ESTABLECIMIENTO DE APOYO', 'DE APOYO')] <- 'ESTABLECIMIENTO DE APOYO'
health$TYPE[health$TYPE %in% c('HOSPITALIZACIÓN', 'DE HOSPITALIZACIÓN')] <- 'HOSPITALIZACIÓN'

health$INSTITUTION[health$INSTITUTION %in% c('ISSSTE', 'IST')] <- 'ISSSTE'
health$INSTITUTION[health$INSTITUTION %in% c('PEMEX', 'PMX')] <- 'PEMEX'
health$INSTITUTION[health$INSTITUTION %in% c('SEDENA', 'SDN')] <- 'SEDENA'
health$INSTITUTION[health$INSTITUTION %in% c('SEMAR', 'SMA')] <- 'SEMAR'
health$INSTITUTION[health$INSTITUTION %in% c('IMSS PROSPERA', 'IMSS-PROSPERA')] <- 'IMSS PROSPERA'
health$INSTITUTION[!(health$INSTITUTION %in% c('SEMAR', 'ISSSTE', 'PEMEX', 'SEDENA', 'SEMAR', 'IMSS PROSPERA', 'SSA', 'SME', 'SMA'))] <- 'OTROS'


#Now let's summarise the data by state, by institution and by type.

health_summary <- health %>% 
  group_by(YEAR, SECTOR, INSTITUTION, STATE_ID, STATE_NAME, TYPE) %>% 
  summarise(BEDS = sum(HOSPITAL_AREA_BEDS, na.rm = T) + sum(OTHER_AREAS_BEDS, na.rm = T),
            CRITICAL_BEDS = sum(ISOLATION_BEDS, na.rm = T) + sum(ER_BEDS, na.rm = T) + sum(INTENSIVE_CARE_BEDS, na.rm = T) + sum(MEDIUM_CARE_BEDS, na.rm = T),
            PERSONNEL = sum(DOCTORS, na.rm = T) + sum(NURSES, na.rm = T) + sum(DOCTORS_IN_STUDY, na.rm = T),
            DOCTORS = sum(DOCTORS, na.rm = T),
            NURSES = sum(NURSES, na.rm = T),
            INTENSIVE_CARE = sum(INTENSIVE_CARE, na.rm = T),
            HOSPITAL_AREA = sum(HOSPITAL_AREA, na.rm = T))

write.csv(file = 'csv_output/health_summary.csv', x = health_summary, row.names = F)
