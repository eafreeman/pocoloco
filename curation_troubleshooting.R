library(data.table)
library(behavr)
library(damr)
library(sleepr)
library(ggetho)
library(zeitgebr)
library(circatools)
library(signal)
library(pracma)
library(ggstatsplot)
library(cowplot)
library(shiny)
library(tidyverse)
library(plotly)


#### Metadata for manual curation of dead animals 



#computer date and time were messed up. Computer recorded the year as 2009 and the time was 1 hour and 35 minutes behind
##this works out to 1:25 pm dusk

manual_metadata <- data.table(file = rep(c("Data/EF_220223_5681/Monitor82.txt",
                                    "Data/EF_220223_5681/Monitor81.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("AGAP005681_GFP"), each = 64),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_31012023/Monitor92.txt"), each =64),
                       incubator = rep(c("1_Behavioural_room"), each = 64),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 64),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2009-02-13 01:25:00", times = 64)),                     
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2009-02-21 15:00:00", times = 64)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 2),
                       ## the genotypes used in the experiment
                       genotype = rep(c("5681",
                                        "GFP"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=64),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 64),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("alive")))

#add in which mosquitoes died 
manual_metadata[c(2,3,5,6,8,9,13,17,20,28,33,36,37,38,40,42,43,50,59,61,48), status := "dead"]


##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

manual_metadata <- link_dam_metadata(manual_metadata, result_dir = data_dir)

##load raw monitor files

manual_dt <- load_dam(manual_metadata[status == "alive"], FUN = sleepr::sleep_dam_annotation)


ggetho(manual_dt, aes(z=activity)) +
  stat_tile_etho() +
  stat_ld_annotations()

ggetho(manual_dt, aes(x=t, z=moving)) + stat_bar_tile_etho()

summary(manual_dt)



###### Metadata for automatic curation of dead animals

metadata <- data.table(file = rep(c("Data/EF_220223_5681/Monitor82.txt",
                                    "Data/EF_220223_5681/Monitor81.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("AGAP005681_GFP"), each = 64),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_31012023/Monitor92.txt"), each =64),
                       incubator = rep(c("1_Behavioural_room"), each = 64),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 64),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2009-02-13 01:25:00", times = 64)),                     
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2009-02-21 15:00:00", times = 64)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 2),
                       ## the genotypes used in the experiment
                       genotype = rep(c("5681",
                                        "GFP"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=64),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 64),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("alive")))

##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

metadata <- link_dam_metadata(metadata, result_dir = data_dir)

##load raw monitor files

dt <- load_dam(metadata, FUN = sleepr::sleep_dam_annotation)

##add simple unique id (uid) and map back to id
dt[, uid := 1 : .N, meta = TRUE]
dt[, .(id, uid) , meta = TRUE]

## curate dead animals or when animals die
dt_curated <- curate_dead_animals(dt, prop_immobile = 0.0001, resolution = 72)
summary(dt_curated)


##curate animals that did not live the whole time
# we make a summary table of all lifespan for each animals
lifespan_dt <- dt_curated[, .(lifespan = max(t)), by=id]
# we filter this table for lifespan>2 and we keep the id
valid_ids <- lifespan_dt[lifespan > days(6), id]
# we apply this filter
dt_curated <- dt_curated[id %in% valid_ids]
summary(dt_curated)

ggetho(dt_curated, aes(z=activity)) +
  stat_tile_etho() +
  stat_ld_annotations()

## see which flies were removed
setdiff(dt[, id, meta=T],
        dt_curated[, id, meta=T])



