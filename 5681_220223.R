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

## use this file to to create a metadata structure that includes all the important information about your experiment
## and any metavariable you may want to analyse across all your different experimental groups

#computer date and time were messed up. Computer recorded the year as 2009 and the time was 1 hour and 35 minutes behind
##this works out to 1:25 pm dusk

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

#add in which mosquitoes died 
metadata[c(2,3,5,6,8,9,13,17,20,28,33,36,37,38,40,42,43,50,59,61,48), status := "dead"]

## now change the status of "dead" or "alive

## add as many other variables as you like you should end up with a data.table that has 1 row for each individual activity
## tube in the experiment with a column for each of the metavariables you include.

#add in statements of who died or survived for future automated stats 
gfp_dead <- 
  dsRNA_dead <- 
  
  
  
  ###########################################################

#Analysis 


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

#change dt to dt_curated because we manually curated when loading in data

dt_curated <- dt

ggetho(dt_curated, aes(z=activity)) + #should be the same as above
  stat_tile_etho() +
  stat_ld_annotations()


## add experiment phase information to each segment of experiment
dt_curated <- dt_curated[, phase := ifelse(t %between% c(days(0), days(3)), "LD1",
                                           ifelse(t %between% c(days(3), days(6)), "FR1",
                                                  ifelse(t %between% c(days(6), days(12)), "VS",
                                                         ifelse(t %between% c(days(12.5), days(20)), "FR2",
                                                                "Not-used"))))]


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


peak_summary <- rejoin(dt_peaks)
peak_summary[, peak_no := as.factor(peak_no)]

grouped_ggbetweenstats(
  data = peak_summary[FR1_rhythmic == TRUE & phase %in% c("FR1")],
  x = entrainment,
  y = peak_time,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR1"
)



dt_peak_sum_id <- peak_summary[phase %in% c("FR1", "FR2"), .(meantime = mean(peak_time),
                                                             medtime = median(peak_time),
                                                             n = length(peak_time),
                                                             sdphase = sd(peak_time)),
                               by = c("genotype", "sex", "entrainment", "phase", "id")]

dt_peak_sum <- peak_summary[phase %in% c("FR1", "FR2"), .(meantime = mean(peak_time),
                                                          medtime = median(peak_time),
                                                          n = length(peak_time),
                                                          sdphase = sd(peak_time)),
                            by = c("genotype", "sex", "entrainment", "phase")]


grouped_ggbetweenstats(
  data = dt_peak_sum_id[phase %in% c("FR1")],
  x = entrainment,
  y = medtime,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR1"
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


dt_auto_sum_id <- auto_summary[phase %in% c("FR1", "FR2") &
                                 peak_no == 3, .(medheight = median(height)),
                               by = c("genotype", "sex", "entrainment", "phase", "id")]

dt_auto_sum <- auto_summary[phase %in% c("FR1", "FR2") &
                              peak_no == 3, .(medheight = median(height)),
                            by = c("genotype", "sex", "entrainment", "phase")]

grouped_ggbetweenstats(
  data = dt_auto_sum_id[phase %in% c("FR1")],
  x = entrainment,
  y = medheight,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR1"
)


## 

dt_summary <- merge.data.table(dt_peak_sum_id, dt_auto_sum_id)

dt_sum <- merge.data.table(dt_peak_sum, dt_auto_sum)


##periodogram - FR1

dt_pgram_FR1 <- periodogram(moving, dt_curated[phase == "FR1"], FUN = chi_sq_periodogram, resample_rate = 1/mins(1))

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



#final visualizations 

ggetho(dt, aes(x=t, y=moving)) + stat_pop_etho() +
  facet_grid(genotype ~ .) #prop moving

ggetho(dt, aes(x=t, y=activity)) + stat_pop_etho() +
  facet_grid(genotype ~ .) #beam breaks 

ggetho(dt, aes(x=t, z=moving)) + stat_bar_tile_etho()

#difference in periods 

per_xsq_dt <- periodogram(activity, 
                          dt_curated,
                          FUN = chi_sq_periodogram)
per_xsq_dt

per_xsq_dt <- find_peaks(per_xsq_dt)
per_xsq_dt

summary_dt <- rejoin(per_xsq_dt[peak==1])
summary_dt

ggplot(summary_dt, aes(genotype, period, fill= genotype)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(size=power -  signif_threshold), alpha=.5) +
  scale_y_hours(name = "Period") 

ggperio(per_xsq_dt) + geom_line(aes(group = id, colour=genotype))

pairwise.wilcox.test(summary_dt$period, summary_dt$genotype) #no difference in GFP and 5681 per


##differences in total activity 

#create new variable

dt_curated[str_detect(id, "Monitor81"), genotype := "GFP"]
dt_curated[str_detect(id, "Monitor82"), genotype := "5681"]

dt_curated[, total_activity := sum(activity), by = id]

pairwise.wilcox.test(dt_curated$total_activity, dt_curated$genotype)

ggplot(data = dt_curated, aes(x = id, y = total_activity, color = genotype)) +
  geom_point() +
  geom_hline(yintercept = mean(dt_curated$total_activity))

##differences in activity at ZT13

ggetho(dt_curated, aes(x=t, y=activity)) + 
  stat_pop_etho() +
  stat_ld_annotations() +
  facet_grid(genotype ~ .)


#average activity for each genotype at each time interval (hour?)
#add in as geom_point to this plot. then do stats tests? 

ggperio(dt_pgram_FR1, aes(period, power, colour=genotype)) + 
  geom_line() +
  geom_peak(colour="blue") +
  facet_wrap( ~ uid) 

## differences in activity at swarming time 

dt_curated[, hour := t/360]



