############################################
##                                        ##
##   Trikinetics Analysis Pipeline        ##
##                                        ##
##      Created by Jason Somers           ##
##    modified by Elizabeth Freeman       ##
##                                        ##
############################################

#have a look at rethomics documentation: https://rethomics.github.io/

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
library(plotly)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(stringr)
#can add other packages but DO NOT add "Lubridate"

################# Metadata #############################

## use this section to to create a metadata structure that includes all the important information about your experiment
## and any metavariable you may want to analyse across all your different experimental groups

metadata <- data.table(file = rep(c("Data/EF_100723/Monitor84.txt",
                                    "Data/EF_100723/Monitor83.txt",
                                    "Data/EF_100723/Monitor82.txt",
                                    "Data/EF_100723/Monitor81.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("Normality"), each = 64),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_100723/Monitor94.txt"), each =64),
                       incubator = rep(c("1_Behavioural_room"), each = 64),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 64),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2023-07-02 02:00:00", times = 64)),   #ZT0                  
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2023-07-10 10:00:00", times = 64)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 2),
                       ## the genotypes used in the experiment
                       genotype = rep(c("Injected",
                                        "Not injected"), each = 64),
                       ## sex of the individuals
                       sex = rep(c("M"), each=64),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 64),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("dead")))

#add in which mosquitoes died 
metadata[c(1,5,6,7,8,9,10,11,12,14,15,16,17,18,20,21,23,26,27,29,30,32,
           33,35,36,37,38,39,40,42,
           65,66,67,68,67,70,71,72,73,
           97,98,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,
           116,117,118,119,120,122,123,124,125,126,127,128), status := "alive"]

## can flip the status of "dead" or "alive depending on what is easiest for your data

## add as many other variables as you like you should end up with a data.table that has 1 row for each individual activity
## tube in the experiment with a column for each of the metavariables you include.

################# Data Cleaning ###########################

##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

metadata <- link_dam_metadata(metadata, result_dir = data_dir)

##load raw monitor files

dt <- load_dam(metadata[status == "alive"], FUN = sleepr::sleep_dam_annotation)

##add simple unique id (uid) and map back to id
dt[, uid := 1 : .N, meta = TRUE]
dt[, .(id, uid) , meta = TRUE]

##visualize all monitor data to see abnormalities
ggetho(dt, aes(z=activity)) +
  stat_tile_etho() +
  stat_ld_annotations()

## add experiment phase information to each segment of experiment (calling this dt_curated)
dt_curated <- dt[, phase := ifelse(t %between% c(days(0), days(3)), "LD",
                                           ifelse(t %between% c(days(3), days(6)), "FR",
                                                                "Not-used"))]
dt_curated[str_detect(id, "Monitor83"), genotype := "injected"]
dt_curated[str_detect(id, "Monitor84"), genotype := "injected"]
dt_curated[str_detect(id, "Monitor82"), genotype := "not injected"]
dt_curated[str_detect(id, "Monitor81"), genotype := "not injected"]

##interactively plot data and adjust phase days as necessary
#runactPlottR()

##add boolean rhythm columns to exclude dead/arr flies
rhythmcols()

##low pass filter design to remove high frequency activity components
bpfilt <- butter(n = 2, W = c((1/hours(8))/((1/60)/2)), type = "low", plane = "z")

##filter activity across entire experiment by individual
dt_curated[, bpfiltered := as.vector(filtfilt(bpfilt, x = activity)),
           by = c("id", "phase")]

##will return calculated phase peaks for all flies in all phases of experiments
##set filters hours and minimum distance between peaks
dt_peaks <- peak_returnR(dt_curated, filterHours = 16, minpeakdist = 18)

##use to remove dead flies
#runpeakPlottR()

##set key and link to metadata
setkeyv(dt_peaks, "id")
setmeta(dt_peaks, dt_curated[, meta = TRUE])

#create bins for activity data

dt_curated[, total_activity := sum(activity), by = id]

dt_curated[, day := floor(t/86400)] #create day variable
dt_curated[, hour := floor(t/3600)] #hour
dt_curated[, half := floor(t/1800)] #half housr

dt_curated[, daily_activity := sum(activity), by = .(id,dt_curated$day)]

dt_curated[, hourly_activity := sum(activity), by = .(id,dt_curated$hour)]
dt_curated[, avg_hourly_activity := mean(hourly_activity), by = .(genotype,dt_curated$hour)]

dt_curated[, halfhourly_activity := sum(activity), by = .(id,dt_curated$half)]
dt_curated[, avg_halfhourly_activity := mean(halfhourly_activity), by = .(genotype,dt_curated$half)]


################ Activity Plots ##########################

#ggetho functions as a ggpplot item and can be manipulated in same way

#basic ggetho plots 

ggetho(dt, aes(x=t, y=moving)) + 
  stat_pop_etho() +
  facet_grid(genotype ~ .) #prop moving

ggetho(dt_curated, aes(x=t, y=activity)) + #beam breaks
  stat_pop_etho() +
  facet_grid(genotype ~ .) +
  stat_ld_annotations() +
  labs(title = "Injected Activity", x = "Day", y = "Beam crosses (per minute)")

ggetho(dt, aes(x=t, z=moving)) + stat_bar_tile_etho()

#stylized plot example
c <- c("#326B9A", "#533F8D")
ggetho(dt_curated, aes(x=t, y=avg_halfhourly_activity, color = genotype)) +
  stat_ld_annotations() +
  stat_ld_annotations(height=1, alpha=0.2, outline = NA) +
  stat_pop_etho() +
  ylim(0, 150) +
  labs(y = "Average Half Hourly Activity") +
  scale_color_manual(values = c) +
  facet_wrap(~genotype) +
  theme_few()

################ Activity Statistics #########################

############# Peak and Period Analysis #######################

peak_summary <- rejoin(dt_peaks)
peak_summary[, peak_no := as.factor(peak_no)]

grouped_ggbetweenstats(
  data = peak_summary[FR_rhythmic == TRUE & phase %in% c("FR")],
  x = entrainment,
  y = peak_time,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR"
)


dt_peak_sum_id <- peak_summary[phase %in% c("FR"), .(meantime = mean(peak_time),
                                                     medtime = median(peak_time),
                                                     n = length(peak_time),
                                                     sdphase = sd(peak_time)),
                               by = c("genotype", "sex", "entrainment", "phase", "id")]

dt_peak_sum <- peak_summary[phase %in% c("FR"), .(meantime = mean(peak_time),
                                                  medtime = median(peak_time),
                                                  n = length(peak_time),
                                                  sdphase = sd(peak_time)),
                            by = c("genotype", "sex", "entrainment", "phase")]


grouped_ggbetweenstats(
  data = dt_peak_sum_id[phase %in% c("FR")],
  x = entrainment,
  y = medtime,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR"
)


## autocorrelation

filterHours <- 16
minpeakdist <- 18

#calculate DAM sample rate to pass to bp filter design
sampRate <- diff(dt_curated[1:2, t])

#low pass filter design to remove high frequency activity components
bpfilt <- butter(n = 2, W = c((1/hours(filterHours))/((1/sampRate)/2)), type = "low", plane = "z")

#filter activity across entire experiment by individual
dt_curated[, bpfiltered := as.vector(filtfilt(bpfilt, x = activity)),
           by = c("id", "phase")]

dt_curated <- dt_curated[, autoco := acf(bpfiltered, lag.max = length(bpfiltered), plot = FALSE)$acf,
                         by = c("id", "phase")]


#find daily peak in activity for each fly and remap it to same metadata
dt_auto <- dt_curated[, data.table(findpeaks(autoco,
                                             zero = "0",
                                             minpeakdistance = 18*60,
                                             peakpat = "[+]{1,}[0]*[-]{1,}",
                                             npeaks = floor(length(t)/1440))),
                      by = c("id", "phase")]

setnames(dt_auto, c("V1", "V2", "V3", "V4"), c("height", "peak", "start", "end"))
setorderv(dt_auto, cols = c("id", "phase",  "peak"))
setkeyv(dt_auto, "id")
metadata <- dt_curated[, meta = TRUE]
setmeta(dt_auto, metadata)
dt_auto[, uid := xmv(uid)]

dt_auto <- dt_auto[, peak_no := rank(peak),
                   by = c("id", "phase")]
dt_auto[, peak_phase := peak%%1440]
dt_auto[, peak_time := peak_phase/1440*24]

setkeyv(dt_auto, "id")
setmeta(dt_auto, dt_curated[, meta = TRUE])


auto_summary <- rejoin(dt_auto)


dt_auto_sum_id <- auto_summary[phase %in% c("FR") &
                                 peak_no == 3, .(medheight = median(height)),
                               by = c("genotype", "sex", "entrainment", "phase", "id")]

dt_auto_sum <- auto_summary[phase %in% c("FR") &
                              peak_no == 3, .(medheight = median(height)),
                            by = c("genotype", "sex", "entrainment", "phase")]

grouped_ggbetweenstats(
  data = dt_auto_sum_id[phase %in% c("FR")],
  x = entrainment,
  y = medheight,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR1"
)


## 

dt_summary <- merge.data.table(dt_peak_sum_id, dt_auto_sum_id)

dt_sum <- merge.data.table(dt_peak_sum, dt_auto_sum)


##periodogram - FR

dt_pgram_FR1 <- periodogram(moving, dt_curated[phase == "FR"], FUN = chi_sq_periodogram, resample_rate = 1/mins(1))

dt_pgram_FR1 <-find_peaks(dt_pgram_FR1)

ggperio(dt_pgram_FR1, aes(y = power - signif_threshold, colour=genotype)) + 
  stat_pop_etho() +
  geom_hline(yintercept = 0) +
  facet_wrap(genotype ~ entrainment)

dt_pgram_FR1_sum <- rejoin(dt_pgram_FR1[peak == 1])
dt_pgram_FR1_sum[, period_h := period/hours(1)]

grouped_ggbetweenstats(
  data = dt_pgram_FR1_sum,
  x = entrainment,
  y = period_h,
  grouping.var = genotype,
  title.text = "Chi^2 periodogram FR1"
)


