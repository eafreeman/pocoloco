############################################
##                                        ##
##   Trikinetics Analysis Pipeline        ##
##                                        ##
##      Created by Jason Somers           ##
##    modified by Elizabeth Freeman       ##
##                                        ##
############################################

## all 5681 trials, big ole' dataset version

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
library(mdthemes)
library(viridis)
library(dplyr)
#can add other packages but DO NOT add "Lubridate"

################# Metadata 1 #############################

## use this section to to create a metadata structure that includes all the important information about your experiment
## and any metavariable you may want to analyse across all your different experimental groups

metadata1 <- data.table(file = rep(c("Data/EF_210524_5681/Monitor81.txt",
                                    "Data/EF_210524_5681/Monitor82.txt",
                                    "Data/EF_210524_5681/Monitor83.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("5681"), each = 96),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_210524_5681/Monitor94.txt"), each = 96),
                       incubator = rep(c("1_Behavioural_room"), each = 96),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 96),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2024-05-13 02:00:00", times = 96)),   #ZT0                  
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2024-05-21 10:00:00", times = 96)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 3),
                       ## the genotypes used in the experiment
                       genotype = rep(c("not_inj",
                                        "ds5681",
                                        "dsgfp"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=96),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 96),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("alive")))

#add in which mosquitoes died 
metadata1[c(5,7,19,23,25,26,28,32,38,40,41,46,56,58,62,63,64,65,66,68,72,74,81,82,85,86,88,89,93,
           94,95), status := "dead"]

## can flip the status of "dead" or "alive depending on what is easiest for your data

## add as many other variables as you like you should end up with a data.table that has 1 row for each individual activity
## tube in the experiment with a column for each of the metavariables you include.

################# Data Cleaning ###########################

##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

metadata1 <- link_dam_metadata(metadata1, result_dir = data_dir)

##load raw monitor files

dt1 <- load_dam(metadata1[status == "alive"], FUN = sleepr::sleep_dam_annotation)

dt1[str_detect(id, "Monitor81"), genotype := "not injected"]
dt1[str_detect(id, "Monitor82"), genotype := "ds5681"]
dt1[str_detect(id, "Monitor83"), genotype := "dsGFP"]

##add simple unique id (uid) and map back to id
dt1[, uid := 1 : .N, meta = TRUE]
dt1[, .(id, uid) , meta = TRUE]




################# Metadata 2#############################

## use this section to to create a metadata structure that includes all the important information about your experiment
## and any metavariable you may want to analyse across all your different experimental groups

metadata2 <- data.table(file = rep(c("Data/EF_080524_5681/Monitor81.txt",
                                    "Data/EF_080524_5681/Monitor82.txt",
                                    "Data/EF_080524_5681/Monitor83.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("5681"), each = 96),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_080524_5681/Monitor94.txt"), each = 96),
                       incubator = rep(c("1_Behavioural_room"), each = 96),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 96),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2024-04-29 02:00:00", times = 96)),   #ZT0                  
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2024-05-07 10:00:00", times = 96)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 3),
                       ## the genotypes used in the experiment
                       genotype = rep(c("dsgfp",
                                        "ds5681",
                                        "not_inj"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=96),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 96),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("alive")))

#add in which mosquitoes died 
metadata2[c(4,7,8,10,11,12,14,22,23, 24, 26,28,29,31,32,35,39,40,41,44,45,46,47,48,49,52,56,58,59,60,64,
           71,81,82,85,86,87,90,91,92,96), status := "dead"]

## can flip the status of "dead" or "alive depending on what is easiest for your data

## add as many other variables as you like you should end up with a data.table that has 1 row for each individual activity
## tube in the experiment with a column for each of the metavariables you include.

################# Data Cleaning ###########################

##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

metadata2 <- link_dam_metadata(metadata2, result_dir = data_dir)

##load raw monitor files

dt2 <- load_dam(metadata2[status == "alive"], FUN = sleepr::sleep_dam_annotation)

##add simple unique id (uid) and map back to id
dt2[, uid := 1 : .N, meta = TRUE]
dt2[, .(id, uid) , meta = TRUE]


#create bins for activity data

dt2[str_detect(id, "Monitor81"), genotype := "dsGFP"]
dt2[str_detect(id, "Monitor82"), genotype := "ds5681"]
dt2[str_detect(id, "Monitor83"), genotype := "not injected"]


############### rep 3 ############################


metadata3 <- data.table(file = rep(c("Data/EF_110624_5681/Monitor81.txt",
                                    "Data/EF_110624_5681/Monitor82.txt",
                                    "Data/EF_110624_5681/Monitor83.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("5681"), each = 96),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_110624_5681/Monitor94.txt"), each = 96),
                       incubator = rep(c("1_Behavioural_room"), each = 96),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 96),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2024-06-03 02:00:00", times = 96)),   #ZT0                  
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2024-06-11 10:00:00", times = 96)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 3),
                       ## the genotypes used in the experiment
                       genotype = rep(c("not_inj",
                                        "ds5681",
                                        "dsgfp"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=96),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 96),
                       ## whether or not the individual survived to the end: d = dead, a = alive
                       status = rep(c("alive")))

#add in which mosquitoes died 
metadata3[c(2,8,13,23,35,43,46,51,52,56,60,67,70,80,85,90,92,94,96), status := "dead"]

data_dir <- "./data"

##link meta and raw monitor files

metadata3 <- link_dam_metadata(metadata3, result_dir = data_dir)

dt3 <- load_dam(metadata3[status == "alive"], FUN = sleepr::sleep_dam_annotation)

##add simple unique id (uid) and map back to id
dt3[, uid := 1 : .N, meta = TRUE]
dt3[, .(id, uid) , meta = TRUE]

dt3[str_detect(id, "Monitor81"), genotype := "not injected"]
dt3[str_detect(id, "Monitor82"), genotype := "ds5681"]
dt3[str_detect(id, "Monitor83"), genotype := "dsGFP"]

################## ALL TOGETHER NOW ##########################

dt_1 <- full_join(dt1, dt2)
dt <- full_join(dt_1, dt3)

## add experiment phase information to each segment of experiment (calling this dt_curated)
dt_curated <- dt[, phase := ifelse(t %between% c(days(0), days(5)), "LD",
                                   ifelse(t %between% c(days(5), days(8)), "FR",
                                          "Not-used"))]

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

#Creating time and activity variables

dt_curated[, total_activity := sum(activity), by = id]

dt_curated[, day := floor(t/86400)] #create day variable
dt_curated[, hour := floor(t/3600)] #hour
dt_curated[, half := floor(t/1800)] #half hour
dt_curated[, min := floor(t/60)]

dt_curated[, daily_activity := sum(activity), by = .(id,dt_curated$day)]

dt_curated[, avg_min_activity := mean(activity), by = .(genotype,min)]

dt_curated[, hourly_activity := sum(activity), by = .(id,dt_curated$hour)]
dt_curated[, avg_hourly_activity := mean(hourly_activity), by = .(genotype,dt_curated$hour)]

dt_curated[, halfhourly_activity := sum(activity), by = .(id,dt_curated$half)]
dt_curated[, avg_halfhourly_activity := mean(halfhourly_activity), by = .(genotype,dt_curated$half)]

dt_curated[, ZT := floor(hour - (day*24))]
dt_curated[, ZT_min := floor(t - (day*86400))]

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
  labs(title = "Activity", x = "Day", y = "Beam crosses (per minute)")

ggetho(dt, aes(x=t, z=moving)) + stat_bar_tile_etho() #individuals

ggetho(dt, aes(x=t, z=moving, y = genotype)) + #grouped
  stat_bar_tile_etho() + 
  stat_ld_annotations() +
  stat_ld_annotations(x_limits = days(c(5,9)), 
                      ld_colours = c("darkgrey", "black" )) +
  labs(y = "", title = "Actogram") +
  theme_few()

#Beam breaks by hour
c <- c("#f768a1","#c51b8a", "#7a0177")
anno_hour <- data.frame(day = c(2:7), x1 = rep(c(12, 0), c(3,3)), x2 = rep(c(24), 6), y1 = rep(c(0), 6), y2 = rep(c(300), 6))
#anno_hour <- data.frame(day = c(7), x1 = rep(c(0)), x2 = rep(c(24)), y1 = rep(c(0)), y2 = rep(c(50)))
ggplot() +
  geom_rect(data = anno_hour, aes(xmin = x1, xmax= x2 , ymin = y1, ymax = y2), 
            fill = "grey85", color = "grey75") +
  geom_bar(data = dt_curated %>% filter(day %in% 2:7), 
           aes(x = ZT, y = avg_hourly_activity, fill = genotype), 
           stat = "summary", fun = "mean", position= "dodge") +
  facet_wrap(~ day) +
  scale_x_continuous(breaks = seq(0, 24, by= 4), expand = c(0,0)) +
  scale_y_continuous(limits=c(0, 300), expand = c(0,0)) +
  labs(y = "Average Hourly Beam Breaks", x = "Hour", fill = "Treatment", 
       title = "Hourly activity", legend = "Treatment") +
  scale_fill_manual(values = c) +
  theme_few()

#Beam breaks by minute

anno_min <- data.frame(day = c(2:7), x1 = rep(c(43200, 0), c(3,3)), x2 = rep(c(86400), 6), y1 = rep(c(0), each = 6), y2 = rep(c(12), each = 6))

ggplot() +
  geom_rect(data = anno_min, aes(xmin = x1, xmax= x2 , ymin = y1, ymax = y2), fill = "grey85", color = "grey75") +
  geom_line(data = dt_curated %>% filter(day %in% 2:7), aes(x = ZT_min, y = avg_min_activity, color = genotype, group = genotype), linewidth = .35) +
  facet_wrap(~ day) +
  scale_x_continuous(breaks = seq(0, 86400, by = 14400), labels = seq(0,24, by = 4), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2), expand = c(0,0)) +
  labs(y = "Average Beam Breaks (per minute)", x = "ZT", color = "Treatment", title = "Activity per Minute") +
  scale_color_brewer(palette = "Set2") +
  theme_few()

#double actograms per genotype

dt_5681 <- filter(dt_curated, genotype == "ds5681")
dt_gfp <- filter(dt_curated, genotype == "dsGFP")
dt_not <- filter(dt_curated, genotype == "not injected")

a <- ggetho(dt_5681, aes(x=t, z=moving), multiplot = 2) + #5681
  stat_bar_tile_etho() +
  ggtitle("*dsAgDopEcR*") +
  md_theme_minimal()

b <- ggetho(dt_gfp, aes(x=t, z=moving), multiplot = 2) + #GFP
  stat_bar_tile_etho() +
  ggtitle("*dsGFP*") +
  md_theme_minimal()

c <- ggetho(dt_not, aes(x=t, z=moving), multiplot = 2) + #not injected
  stat_bar_tile_etho() +
  ggtitle("Not Injected") +
  md_theme_minimal()

ggarrange(a,b,c, labels = c("A", "B", "C"))

################ Activity Statistics #########################

#these are wrong, df should only be 2! 
kruskal.test(dt_curated$daily_activity, dt_curated$genotype) #all "genotypes"

dt_curated_inj <- filter(dt_curated, genotype != "not injected") #filters for only injected
kruskal.test(dt_curated_inj$daily_activity, dt_curated_inj$genotype) #just injected

dt_curated_inj_LD <- filter(dt_curated_inj, day == 2:5)
kruskal.test(dt_curated_inj_LD$daily_activity, dt_curated_inj_LD$genotype)

dt_curated_inj_DD <- filter(dt_curated_inj, day == 6:8)
kruskal.test(dt_curated_inj_DD$daily_activity, dt_curated_inj_DD$genotype)

median(dt_curated_inj_DD$daily_activity[dt_curated_inj_DD$genotype == "ds5681"])
median(dt_curated_inj_LD$daily_activity[dt_curated_inj_LD$genotype == "ds5681"])

median(dt_curated_inj_DD$daily_activity[dt_curated_inj_DD$genotype == "dsGFP"])
median(dt_curated_inj_LD$daily_activity[dt_curated_inj_LD$genotype == "dsGFP"])

t.test(dt_curated_inj_DD$daily_activity ~ dt_curated_inj_DD$genotype)

#trying something else for stats

#want to make a dataset with means from each replicate, 5 df total
#want to break it down by daily activity LD and DD of days 2-4 (and total FR)

#unique dataset

dt_curated_inj_LD <- filter(dt_curated_inj, day == 2:4)
unique_LDdays <- dt_curated_inj_LD %>% distinct(daily_activity, .keep_all = TRUE)

kruskal.test(unique_LDdays$daily_activity, unique_LDdays$genotype)

dt_curated_inj_DD <- filter(dt_curated_inj, day == 6:8)
unique_DDdays <- dt_curated_inj_DD %>% distinct(daily_activity, .keep_all = TRUE)

kruskal.test(unique_DDdays$daily_activity, unique_DDdays$genotype)

t.test(unique_LDdays$daily_activity ~ unique_LDdays$genotype)
t.test(unique_DDdays$daily_activity ~ unique_DDdays$genotype)

#plotting daily averages 

dt_curated_inj <- dt_curated_inj %>%
  group_by(genotype, day) %>%
  mutate(mean_da = mean(daily_activity)) %>%
  mutate(se_da = sd(daily_activity)/sqrt(30)) %>% #need to fix length for sqrt
  ungroup()

unique_days <- dt_curated_inj %>% distinct(daily_activity, .keep_all = TRUE)

unique_days %>%
ggplot(aes(x = day, y = mean_da, fill = genotype)) +
  geom_bar(stat="identity", position = position_dodge(), color = "darkgray") +
  geom_errorbar(aes(ymin = mean_da-se_da, ymax = mean_da+se_da), position = position_dodge(.9), width = .2, color = "gray34") +
  theme_few()
  

#Scaling down df

hours <- dt_curated %>% distinct(hour, genotype, .keep_all = T)
write.csv(hours, "avg_hourly_activity_270624.csv")

mins <- dt_curated %>% distinct(min, genotype, .keep_all = T)
mins <-  filter(mins, day == 2:7)


c <- c("#f768a1","#c51b8a", "#7a0177")
anno_hour <- data.frame(day = c(2:7), x1 = rep(c(12, 0), c(3,3)), x2 = rep(c(24), 6), y1 = rep(c(0), 6), y2 = rep(c(300), 6))
#anno_hour <- data.frame(day = c(7), x1 = rep(c(0)), x2 = rep(c(24)), y1 = rep(c(0)), y2 = rep(c(50)))
ggplot() +
  geom_rect(data = anno_hour, aes(xmin = x1, xmax= x2 , ymin = y1, ymax = y2), 
            fill = "grey85", color = "grey75") +
  geom_bar(data = hours %>% filter(day %in% 2:7), 
           aes(x = ZT, y = avg_hourly_activity, fill = genotype), 
           stat = "summary", fun = "mean", position= "dodge") +
  facet_wrap(~ day) +
  scale_x_continuous(breaks = seq(0, 24, by= 4), expand = c(0,0)) +
  scale_y_continuous(limits=c(0, 300), expand = c(0,0)) +
  labs(y = "Average Hourly Beam Breaks", x = "Hour", fill = "Treatment", 
       title = "Hourly activity", legend = "Treatment") +
  scale_fill_manual(values = c) +
  theme_few()
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


dt_peak_sum_id[str_detect(id, "04-29 02:00:00|Monitor81"), genotype := "dsgfp"]
dt_peak_sum_id[str_detect(id, "04-29 02:00:00|Monitor83"), genotype := "not_inj"]
dt_peak_sum_id[str_detect(id, "06-03 02:00:00|Monitor81"), genotype := "dsgfp"]
dt_peak_sum_id[str_detect(id, "06-03 02:00:00|Monitor83"), genotype := "not_inj"]
dt_peak_sum_id[str_detect(id, "Monitor82"), genotype := "ds5681"]

dt_peak_sum <- peak_summary[phase %in% c("FR"), .(meantime = mean(peak_time),
                                                  medtime = median(peak_time),
                                                  n = length(peak_time),
                                                  sdphase = sd(peak_time)),
                            by = c("genotype", "sex", "entrainment", "phase")]


grouped_ggbetweenstats(
  data = dt_peak_sum_id[phase %in% c("FR")],
  x = genotype,
  y = medtime,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR"
)


ggbetweenstats(
  data = dt_peak_sum_id[phase %in% c("FR")],
  x = genotype,
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

dt_pgram_FR1_sum[str_detect(id, "04-29 02:00:00|Monitor81"), genotype := "dsgfp"]
dt_pgram_FR1_sum[str_detect(id, "04-29 02:00:00|Monitor83"), genotype := "not_inj"]
dt_pgram_FR1_sum[str_detect(id, "06-03 02:00:00|Monitor81"), genotype := "dsgfp"]
dt_pgram_FR1_sum[str_detect(id, "06-03 02:00:00|Monitor83"), genotype := "not_inj"]
dt_pgram_FR1_sum[str_detect(id, "Monitor82"), genotype := "ds5681"]

grouped_ggbetweenstats(
  data = dt_pgram_FR1_sum,
  x = entrainment,
  y = period_h,
  grouping.var = genotype,
  title.text = "Chi^2 periodogram FR1",
)

#I think this is the proper one 
ggbetweenstats(
  data = dt_pgram_FR1_sum,
  x = genotype,
  y = period_h,
  title = "Chi^2 periodogram",
  ylab = "Period",
  xlab = "Genotype",
  pairwise.display = "significant"
)


