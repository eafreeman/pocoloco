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

## use this file to to create a metadata structure that includes all the important information about your experiment
## and any metavariable you may want to analyse across all your different experimental groups

metadata <- data.table(file = rep(c("Data/EF_31012023/Monitor81.txt",
                                    "Data/EF_31012023/Monitor79.txt"), each = 32),
                       ## a unique experiment ID if you want to analyse across experiments
                       exp_ID = rep(c("AGAP005681_GFP"), each = 64),
                       ## a link to the environmental monitor file use in the experiment
                       env_monitor = rep(c("Data/EF_31012023/Monitor92.txt"), each =64),
                       incubator = rep(c("1_Behavioural_room"), each = 64),
                       ##entrainment = do these for each ZT 
                       entrainment = rep(c("LD_DD"), each = 64),
                       ##startdatetime at lights on day before 
                       start_datetime = 
                         c(rep("2023-01-24 02:00:00", times = 64)),                     
                       ## experimental stop time
                       stop_datetime = 
                         c(rep("2023-01-31 14:05:00", times = 64)),
                       ## this is the tube number position in each monitor
                       region_id = rep(1:32, 2),
                       ## the genotypes used in the experiment
                       genotype = rep(c("GFP",
                                        "5681"), each = 32),
                       ## sex of the individuals
                       sex = rep(c("M"), each=64),
                       ## temperature of experiment
                       temp = rep(c("28"), each = 64))

## this creates another metavariable "treatment" for specific combinations of varibales that might be interesting
## here it creates one for genotype ~ sex
metadata[, treatment := paste(genotype, sex, sep='_')]

## add as many other variables as you like you should end up with a data.table that has 1 row for each individual activity
## tube in the experiment with a column for each of the metavariables you include.


###########################################################

#Analysis 

##set directory of DAM monitor files

data_dir <- "./data"

##link meta and raw monitor files

metadata <- link_dam_metadata(metadata, result_dir = data_dir)

##load raw monitor files

dt <- load_dam(metadata, FUN = sleepr::sleep_dam_annotation)

##add simple unique id (uid) and map back to id
dt[, uid := 1 : .N, meta = TRUE]
dt[, .(id, uid) , meta = TRUE]


##visualize all monitor data to see abnormalities
ggetho(dt, aes(z=activity)) +
  stat_tile_etho() +
  stat_ld_annotations()


## curate dead animals or when animals die
dt_curated <- curate_dead_animals(dt, prop_immobile = 0.005, resolution = 24)

summary(dt_curated)

ggetho(dt_curated, aes(x=t, z=moving)) + stat_bar_tile_etho()

ggetho(dt_curated, aes(z=activity)) +
  stat_tile_etho() +
  stat_ld_annotations()

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


##more viz

ggetho(dt_curated, aes(x=t, y=moving)) + stat_pop_etho() +
  facet_grid(genotype ~ .) #prop moving

ggetho(dt_curated, aes(x=t, y=activity)) + stat_pop_etho() +
  facet_grid(genotype ~ .) #beam breaks 

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

## add experiment phase information to each segment of experiment
dt_curated <- dt[, phase := ifelse(t %between% c(days(1), days(3)), "LD1",
                                   ifelse(t %between% c(days(3), days(8)), "FR1",
                                          "Not-used"))]


##interactively plot data and adjust phase days as necessary
runactPlottR()

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
runpeakPlottR()

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


ggplot(data = dt_peak_sum_id[phase %in% "FR2"], aes(x = phase, y = medtime, colour = entrainment)) +
  geom_point() +
  geom_segment(data = dt_peak_sum[phase %in% "FR2"], aes(x = 0, y = medtime, yend = medtime, xend = phase), linewidth = 1.5) +
  coord_polar(theta = "y", start = 0) +
  scale_y_continuous(limits = c(0,24), breaks = seq(0, 24, 6)) +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  labs(x = "Rhythmicity (A.U.)", y = "Peak phase (H)") +
  facet_wrap(. ~ genotype) +
  theme_minimal_grid()




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

grouped_ggbetweenstats(
  data = dt_auto_sum_id[phase %in% c("FR2")],
  x = entrainment,
  y = medheight,
  grouping.var = genotype,
  mean.plotting = FALSE,
  title.text = "peak phase FR2"
)

## 

dt_summary <- merge.data.table(dt_peak_sum_id, dt_auto_sum_id)

dt_sum <- merge.data.table(dt_peak_sum, dt_auto_sum)

ggplot(data = dt_summary[phase %in% "FR2"], aes(x = medheight, y = medtime, colour = entrainment)) +
  geom_point(alpha = 0.5) +
  geom_segment(data = dt_sum[phase %in% "FR2"], aes(x = 0, y = medtime, yend = medtime, xend = medheight), linewidth = 1.5) +
  coord_polar(theta = "y", start = 0) +
  scale_y_continuous(limits = c(0,24), breaks = seq(0, 24, 6)) +
  scale_x_continuous(limits = c(0,0.4), breaks = seq(0, 0.5, 0.1)) +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  labs(x = "Rhythmicity (A.U.)", y = "Peak phase (H)") +
  #facet_wrap(. ~ genotype) +
  theme_minimal_grid()



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

