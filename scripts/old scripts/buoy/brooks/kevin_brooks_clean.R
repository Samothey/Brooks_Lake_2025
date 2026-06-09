## kevins script for cleaning brooks lake buoy data 

#Load required packages
library(ggplot2)
library(dplyr)
library(scales)
library(lubridate)
library(caTools)
library(zoo)


# Set "dir" to location of raw buoy data files
dir <- "G:/.shortcut-targets-by-id/1JM2QoknF8B3Pc75S9OrH4J_-4ONf7H81/WYACT Team Aquatic Shared/Project Data/Buoys/2025/Brooks"

# FOR DATA FROM CAMPBELL DATALOGGER ON BUOY

# Read in Brooks Lake buoy data, add column names, format datetime, and select desired columns and order them
brooks_wq <- read.csv(paste0(dir,"/Raw/Brooks_2025_WQ_raw_final.csv"))
brooks_wq$X <- NULL

# Add 0:00 timestamp where it is needed: R has issues reading/writing midnight timestamps sometimes
brooks_wq$DATETIME_UTC <- as.POSIXct(
  ifelse(
    grepl(":", brooks_wq$DATETIME_UTC),                     # has time already, leave alone
    brooks_wq$DATETIME_UTC,
    paste(brooks_wq$DATETIME_UTC, "00:00:00")), tz = "UTC") # add midnight

brooks_wq$DATETIME_MST <- with_tz(brooks_wq$DATETIME_UTC,"MST") # add column for local time (MST)


# MINIDOT and MINIPAR DATA

# Load data (9 m depth for PAR, 15 m depth for minidot)
brooks_DO <- read.csv(paste0(dir,"/Raw/Brooks_DO_303968_15m_2025.csv"),skip=8)
brooks_PAR <- read.csv(paste0(dir,"/Raw/Brooks_PAR_666312_9m_2025.csv"),skip=6)

# Add column names
colnames(brooks_DO) <- c("TIMESTAMP_UNIX","DATETIME_UTC","DATETIME_MST","batV","WTemp15m","odomgL_15m","odosat_15m","Q")
colnames(brooks_PAR) <- c("TIMESTAMP_UNIX","DATETIME_UTC","DATETIME_MST","batV","WTemp9m","PAR9m","acc_x","acc_y","acc_z")

# Format datetime column
brooks_DO$DATETIME_UTC <- as.POSIXct(brooks_DO$DATETIME_UTC,format="%Y-%m-%d %H:%M:%S",tz="UTC")
brooks_PAR$DATETIME_UTC <- as.POSIXct(brooks_PAR$DATETIME_UTC,format="%Y-%m-%d %H:%M:%S",tz="UTC")

# Round timestamps to nearest 15 min interval
output <- brooks_DO %>%
  mutate(bin = floor_date(DATETIME_UTC, "15 mins")) %>% # change ceiling to floor if timestamps should be rounded down instead of up
  group_by(bin) %>%
  summarise(DATETIME_UTC = bin,
            odomgL_15m = mean(odomgL_15m),
            odosat_15m = mean(odosat_15m),
            WTemp15m = mean(WTemp15m), .groups = "drop")

output2 <- brooks_PAR %>%
  mutate(bin = floor_date(DATETIME_UTC, "15 mins")) %>%
  group_by(bin) %>%
  summarise(DATETIME_UTC = bin,
            PAR9m = mean(PAR9m),
            WTemp9m = mean(WTemp9m), .groups = "drop")

# Remove unwanted column
output$bin <- NULL
output2$bin <- NULL

# Join output files together
output <- left_join(output,output2,by=c("DATETIME_UTC"))

# Convert to 10 minute interval data before merging with buoy data  **Not ran for 2025 due to timestamp correction

# Create 10-minute sequence from first to last timestamp
#new_times <- data.frame(
#  DATETIME_UTC = seq(
#    from = floor_date(min(output$DATETIME_UTC), "10 minutes"),
#    to   = ceiling_date(max(output$DATETIME_UTC), "10 minutes"),
#    by   = "10 min"))

# Merge and use linear interpolation to fill NAs
#output_10min <- new_times %>%
#  left_join(output, by = "DATETIME_UTC") %>%
#  mutate(odomgL_15m = na.approx(odomgL_15m, x = DATETIME_UTC, na.rm = FALSE),
#         odosat_15m = na.approx(odosat_15m, x = DATETIME_UTC, na.rm = FALSE),
#         WTemp15m   = na.approx(WTemp15m, x = DATETIME_UTC, na.rm = FALSE),
#         PAR9m      = na.approx(PAR9m, x = DATETIME_UTC, na.rm = FALSE),
#         WTemp9m    = na.approx(WTemp9m, x = DATETIME_UTC, na.rm = FALSE))


# Trim to correct start/end times based on buoy deployment and retrieval
output <- output %>% filter(DATETIME_UTC>=as.POSIXct('2025-06-10 16:30:00',tz="UTC"))
output <- output %>% filter(DATETIME_UTC<=as.POSIXct('2025-10-14 22:10:00',tz="UTC"))


# Add minidot data to final dataframe
brooks_wq <- brooks_wq %>%
  distinct(DATETIME_UTC, .keep_all = TRUE) %>%     # fix duplicate keys
  full_join(output %>% distinct(DATETIME_UTC, .keep_all = TRUE),
            by = "DATETIME_UTC",
            suffix = c("", "_out")) %>%
  mutate(WTemp15m = coalesce(WTemp15m, WTemp15m_out)) %>%
  select(-WTemp15m_out) %>%
  arrange(DATETIME_UTC)

# Add in missing timestamps in MST
brooks_wq$DATETIME_MST <- with_tz(brooks_wq$DATETIME_UTC,"MST") # add column for local time (MST)

# Re-order desired columns
brooks_wq <- brooks_wq %>% select(DATETIME_UTC,DATETIME_MST,WTempC,ConduS,SpConduS,CHLrfu,CHLugL,
                                  BGAPCrfu,BGAPCugL,odomgL,odosat,WTemp1m,WTemp4m,WTemp7m,WTemp9m,WTemp10m,
                                  WTemp13m,WTemp15m,odomgL_15m,odosat_15m,PAR9m)



# PLOTTING

# Plot all variables for visual inspection
brooks_wq %>% tidyr::gather("id", "value", 3:21) %>% 
  ggplot(., aes(DATETIME_MST, value)) + geom_point() + facet_wrap(~id,scales="free_y") + theme_bw()


# Re-plot all variables for visual inspection
#brooks_wq %>% tidyr::gather("id", "value", 2:11) %>% 
#  ggplot(., aes(DATETIME_MST, value)) + geom_point() + facet_wrap(~id)

# Start saving plots for quick reference later
ggsave(paste0(dir,"/Figures/allvars.png")) # save most recent plot in viewer to desired location

# Plot EXO temp for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, WTempC)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/temp_EXO.png")) # save plot

# Plot conductivity for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, ConduS)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/cond.png")) # save plot

# Plot spec conductivity for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, SpConduS)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/speccond.png")) # save plot

# Plot chl-a (rfu) for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, CHLrfu)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/CHLrfu.png")) # save plot

# Plot BGA (rfu) for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, BGAPCrfu)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/BGArfu.png")) # save plot

# Plot DO (mg/L) for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, odomgL)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/DOmgL.png")) # save plot

# Plot DO (% sat) for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, odosat)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/DOsat.png")) # save plot

# Plot DO (mg/L) at 15 m for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, odomgL_15m)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/DOmgL_15m.png")) # save plot

# Plot DO (% sat) at 15 m for visual inspection and save
ggplot(brooks_wq, aes(DATETIME_MST, odosat_15m)) + geom_point(na.rm=TRUE) + theme_bw()
ggsave(paste0(dir,"/Figures/DOsat_15m.png")) # save plot

# Plot PAR at 9 m for visual inspection and save
p <- ggplot(brooks_wq, aes(DATETIME_UTC, PAR9m)) + geom_point(na.rm=TRUE) + theme_bw()
ggplotly(p)
ggsave(paste0(dir,"/Figures/PAR_9m.png")) # save plot

# Plot salinity for visual inspection and save
#ggplot(brooks_wq, aes(DATETIME_MST, salpsu)) + geom_point(na.rm=TRUE) + theme_bw()
#ggsave("salinity.png") # save plot

# Plot total dissolved solids for visual inspection and save
#ggplot(brooks_wq, aes(DATETIME_MST, TDS)) + geom_point(na.rm=TRUE) + theme_bw()
#ggsave("TDS.png") # save plot

# Plot pH for visual inspection and save
#ggplot(brooks_wq, aes(DATETIME_MST, pH)) + geom_point(na.rm=TRUE) + theme_bw()
#ggsave("pH.png") # save plot

# Plot ORP for visual inspection and save
#ggplot(brooks_wq, aes(DATETIME_MST, ORP)) + geom_point(na.rm=TRUE) + theme_bw()
#ggsave("ORP.png") # save plot

# Plot temperature string timeseries
legend_colors <- c("1 m" = "black",   # set unique colors for each depth
                   "4 m" = "orange", 
                   "7 m" = "sky blue",
                   "9 m" = "#009E73",
                   "10 m" = "#D55E00",
                   "13 m" = "#CC79A7",
                   "15 m" = "#0072B2")

p1 <- ggplot(brooks_wq) + geom_line(aes(DATETIME_UTC, WTemp1m, col="1 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp4m, col="4 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp7m, col="7 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp9m, col="9 m")) +
  geom_line(aes(DATETIME_UTC, WTemp10m, col="10 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp13m, col="13 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp15m, col="15 m")) + ggtitle("Brooks Lake Water Temp 2025") + 
  xlab("Date") + ylab("Temp (C)") + ylim(4,20) +
  scale_color_manual(name="",values=legend_colors,breaks=c("1 m","4 m","7 m","9 m","10 m","13 m","15 m")) +
  theme_bw()
ggplotly(p1)
ggsave(paste0(dir,"/Figures/tempstring.png")) # save plot


# Plot temp profile heat map for Brooks Lake
heat <- brooks_wq[,c(1,12,13,14,15,16,17,18)]
colnames(heat) <- c("datetime","var_1","var_4","var_7","var_9","var_10","var_13","var_15")
heat$datetime <- format(heat$datetime, "%Y-%m-%d %H:%M:%S")
heat <- heat %>% na.omit(heat)
write.table(heat,paste0(dir,"/Figures/brookstemp2025_rlakeanalyzer_format.txt"),sep="\t",
            row.names = FALSE,
            quote = FALSE)
data2=load.ts(paste0(dir,"/Figures/brookstemp2025_rlakeanalyzer_format.txt"))
wtr.heat.map(data2) # save plot manually from viewer



###################
# Remove outliers #
###################
brooks_wq$BGAPCrfu_corr <- brooks_wq$BGAPCrfu

# Calculate rolling average using surrounding 30 values using runmean in caTools package
brooks_wq$moving_av_30obs <- runmean(brooks_wq$BGAPCrfu_corr,30,alg="C",endrule="mean",align="center")

# Calculate rolling standard deviation using surrounding 30 values using runsd in caTools package
brooks_wq$moving_sd_30obs <- runsd(brooks_wq$BGAPCrfu_corr,30,center=runmean(brooks_wq$BGAPCrfu_corr,30),endrule="sd",align="center") 

# Calculate z score based on rolling mean and sd
brooks_wq$zscore <- abs((brooks_wq$BGAPCrfu_corr-brooks_wq$moving_av_30obs)/brooks_wq$moving_sd_30obs)

# Set values to NA if the z score is >4
for (i in 1:nrow(brooks_wq)) {
  if (!is.na(brooks_wq$zscore[i]) && brooks_wq$zscore[i] > 4) {
    brooks_wq$BGAPCrfu_corr[i] <- NA
  }}

ggplot() + geom_point(data=brooks_wq,aes(DATETIME_MST,BGAPCrfu)) + geom_point(data=brooks_wq,aes(DATETIME_MST,BGAPCrfu_corr),color="red") + xlab("Time") + ylab("BGA (rfu)")

brooks_wq$BGAPCugL <- ifelse(is.na(brooks_wq$BGAPCrfu_corr), NA, brooks_wq$BGAPCugL)
brooks_wq$BGAPCrfu <- brooks_wq$BGAPCrfu_corr

brooks_wq$BGAPCrfu_corr <- NULL
brooks_wq$moving_av_30obs <- NULL
brooks_wq$moving_sd_30obs <- NULL
brooks_wq$zscore <- NULL


## Chl-a
brooks_wq$CHLrfu_corr <- brooks_wq$CHLrfu

# Calculate rolling average using surrounding 30 values using runmean in caTools package
brooks_wq$moving_av_30obs <- runmean(brooks_wq$CHLrfu_corr,30,alg="C",endrule="mean",align="center")

# Calculate rolling standard deviation using surrounding 30 values using runsd in caTools package
brooks_wq$moving_sd_30obs <- runsd(brooks_wq$CHLrfu_corr,30,center=runmean(brooks_wq$CHLrfu_corr,30),endrule="sd",align="center") 

# Calculate z score based on rolling mean and sd
brooks_wq$zscore <- abs((brooks_wq$CHLrfu_corr-brooks_wq$moving_av_30obs)/brooks_wq$moving_sd_30obs)

# Set values to NA if the z score is >4
for (i in 1:nrow(brooks_wq)) {
  if (!is.na(brooks_wq$zscore[i]) && brooks_wq$zscore[i] > 4) {
    brooks_wq$CHLrfu_corr[i] <- NA
  }}

ggplot() + geom_point(data=brooks_wq,aes(DATETIME_MST,CHLrfu)) + geom_point(data=brooks_wq,aes(DATETIME_MST,CHLrfu_corr),color="red") + xlab("Time") + ylab("Chla (rfu)")

brooks_wq$CHLugL <- ifelse(is.na(brooks_wq$CHLrfu_corr), NA, brooks_wq$CHLugL)
brooks_wq$CHLrfu <- brooks_wq$CHLrfu_corr

brooks_wq$CHLrfu_corr <- NULL
brooks_wq$moving_av_30obs <- NULL
brooks_wq$moving_sd_30obs <- NULL
brooks_wq$zscore <- NULL


brooks_wq$DATETIME_UTC <- format(
  brooks_wq$DATETIME_UTC,
  "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

brooks_wq$DATETIME_MST <- format(
  brooks_wq$DATETIME_MST,
  "%Y-%m-%d %H:%M:%S",
  tz = "MST"
)


################################################################################
# Correct temperatures after 7/23/25 17:10 at 7, 9, 10, and 13 m depths due to 
#   the temp string (below 4 m) and PAR sensor shifting depths; use vertical
#   profile data to offset buoy temperatures

# Read in Brooks Lake buoy and profile temp data for comparison
d1 <- read.csv(paste0(dir,"/Raw/Brooks_2025_buoy_profile_temp_comp.csv"))
d1$DATETIME_UTC <- as.POSIXct(d1$DATETIME_UTC,format="%m/%d/%Y %H:%M",tz="UTC")

# Subset comparison data by depth
d7 <- d1 %>% filter(Depth==7)
d9 <- d1 %>% filter(Depth==9)
d10 <- d1 %>% filter(Depth==10)
d13 <- d1 %>% filter(Depth==13)


# Get timestamps to fill with corrected temperature data
raw7 <- brooks_wq %>% select(DATETIME_UTC, WTemp7m) %>% filter(DATETIME_UTC>as.POSIXct('2025-07-23 17:10:00',tz="UTC"))

# Add in offset values at corresponding datetimes for a specific depth
raw7 <- left_join(raw7,d7,by=c("DATETIME_UTC"))

# Linearly interpolate the rest of the time offset values in between vertical profiles
raw7$offset_interp <- approx(
  x    = raw7$DATETIME_UTC,
  y    = raw7$Offset,
  xout = raw7$DATETIME_UTC,   # ← interpolate onto FULL raw time series
  rule = 2
)$y

# Calculate corrected temps based on offset at each timestamp
raw7$WTemp7m_corrected <- raw7$WTemp7m - raw7$offset_interp
raw7 <- raw7 %>% select(DATETIME_UTC, WTemp7m_corrected)

# Add corrected temps back into original dataframe
brooks_wq <- left_join(brooks_wq,raw7,by=c("DATETIME_UTC"))

# Use original temperature values for timestamps before 7/23/25 17:10
brooks_wq <- brooks_wq %>% mutate(WTemp7m_corrected=coalesce(WTemp7m_corrected, WTemp7m))

# Plot original vs corrected data
p7 <- ggplot(brooks_wq, aes(x = DATETIME_UTC)) +
  geom_point(aes(y = WTemp7m, color = "Original"),size=1) +
  geom_point(aes(y = WTemp7m_corrected, color = "Corrected"),size=1) +
  scale_color_manual(
    values = c(
      Original = "grey50",
      Corrected = "blue")) + theme_bw()
ggplotly(p7)



## Repeat above code for other depths 
# Get timestamps to fill with corrected temperature data
raw9 <- brooks_wq %>% select(DATETIME_UTC, WTemp9m) %>% filter(DATETIME_UTC>as.POSIXct('2025-07-23 17:10:00',tz="UTC"))

# Add in offset values at corresponding datetimes for a specific depth
raw9 <- left_join(raw9,d9,by=c("DATETIME_UTC"))

# Linearly interpolate the rest of the time offset values in between vertical profiles
raw9$offset_interp <- approx(
  x    = raw9$DATETIME_UTC,
  y    = raw9$Offset,
  xout = raw9$DATETIME_UTC,   # ← interpolate onto FULL raw time series
  rule = 2
)$y

# Calculate corrected temps based on offset at each timestamp
raw9$WTemp9m_corrected <- raw9$WTemp9m - raw9$offset_interp
raw9 <- raw9 %>% select(DATETIME_UTC, WTemp9m_corrected)

# Add corrected temps back into original dataframe
brooks_wq <- left_join(brooks_wq,raw9,by=c("DATETIME_UTC"))

# Use original temperature values for timestamps before 9/23/25 19:10
brooks_wq <- brooks_wq %>% mutate(WTemp9m_corrected=coalesce(WTemp9m_corrected, WTemp9m))

# Plot original vs corrected data
p9 <- ggplot(brooks_wq, aes(x = DATETIME_UTC)) +
  geom_point(aes(y = WTemp9m, color = "Original"),size=1) +
  geom_point(aes(y = WTemp9m_corrected, color = "Corrected"),size=1) +
  scale_color_manual(
    values = c(
      Original = "grey50",
      Corrected = "blue")) + theme_bw()
ggplotly(p9)


# Get timestamps to fill with corrected temperature data
raw10 <- brooks_wq %>% select(DATETIME_UTC, WTemp10m) %>% filter(DATETIME_UTC>as.POSIXct('2025-07-23 17:10:00',tz="UTC"))

# Add in offset values at corresponding datetimes for a specific depth
raw10 <- left_join(raw10,d10,by=c("DATETIME_UTC"))

# Linearly interpolate the rest of the time offset values in between vertical profiles
raw10$offset_interp <- approx(
  x    = raw10$DATETIME_UTC,
  y    = raw10$Offset,
  xout = raw10$DATETIME_UTC,   # ← interpolate onto FULL raw time series
  rule = 2
)$y

# Calculate corrected temps based on offset at each timestamp
raw10$WTemp10m_corrected <- raw10$WTemp10m - raw10$offset_interp
raw10 <- raw10 %>% select(DATETIME_UTC, WTemp10m_corrected)

# Add corrected temps back into original dataframe
brooks_wq <- left_join(brooks_wq,raw10,by=c("DATETIME_UTC"))

# Use original temperature values for timestamps before 10/23/25 110:10
brooks_wq <- brooks_wq %>% mutate(WTemp10m_corrected=coalesce(WTemp10m_corrected, WTemp10m))

# Plot original vs corrected data
p10 <- ggplot(brooks_wq, aes(x = DATETIME_UTC)) +
  geom_point(aes(y = WTemp10m, color = "Original"),size=1) +
  geom_point(aes(y = WTemp10m_corrected, color = "Corrected"),size=1) +
  scale_color_manual(
    values = c(
      Original = "grey50",
      Corrected = "blue")) + theme_bw()
ggplotly(p10)


# Get timestamps to fill with corrected temperature data
raw13 <- brooks_wq %>% select(DATETIME_UTC, WTemp13m) %>% filter(DATETIME_UTC>as.POSIXct('2025-07-23 17:10:00',tz="UTC"))

# Add in offset values at corresponding datetimes for a specific depth
raw13 <- left_join(raw13,d13,by=c("DATETIME_UTC"))

# Linearly interpolate the rest of the time offset values in between vertical profiles
raw13$offset_interp <- approx(
  x    = raw13$DATETIME_UTC,
  y    = raw13$Offset,
  xout = raw13$DATETIME_UTC,   # ← interpolate onto FULL raw time series
  rule = 2
)$y

# Calculate corrected temps based on offset at each timestamp
raw13$WTemp13m_corrected <- raw13$WTemp13m - raw13$offset_interp
raw13 <- raw13 %>% select(DATETIME_UTC, WTemp13m_corrected)

# Add corrected temps back into original dataframe
brooks_wq <- left_join(brooks_wq,raw13,by=c("DATETIME_UTC"))

# Use original temperature values for timestamps before 13/23/25 113:13
brooks_wq <- brooks_wq %>% mutate(WTemp13m_corrected=coalesce(WTemp13m_corrected, WTemp13m))

# Plot original vs corrected data
p13 <- ggplot(brooks_wq, aes(x = DATETIME_UTC)) +
  geom_point(aes(y = WTemp13m, color = "Original"),size=1) +
  geom_point(aes(y = WTemp13m_corrected, color = "Corrected"),size=1) +
  scale_color_manual(
    values = c(
      Original = "grey50",
      Corrected = "blue")) + theme_bw()
ggplotly(p13)


# Plot full corrected temp string timeseries
p2 <- ggplot(brooks_wq) + geom_line(aes(DATETIME_UTC, WTemp1m, col="1 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp4m, col="4 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp7m_corrected, col="7 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp9m_corrected, col="9 m")) +
  geom_line(aes(DATETIME_UTC, WTemp10m_corrected, col="10 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp13m_corrected, col="13 m")) + 
  geom_line(aes(DATETIME_UTC, WTemp15m, col="15 m")) + ggtitle("Brooks Lake Water Temp 2025") + 
  xlab("Date") + ylab("Temp (C)") + ylim(4,20) +
  scale_color_manual(name="",values=legend_colors,breaks=c("1 m","4 m","7 m","9 m","10 m","13 m","15 m")) +
  theme_bw()
ggplotly(p2)
ggsave(paste0(dir,"/Figures/tempstring.png")) # save plot



# Load DEQ profile data for comparison
brooks_deq <- read.csv("G:/.shortcut-targets-by-id/1JM2QoknF8B3Pc75S9OrH4J_-4ONf7H81/WYACT Team Aquatic Shared/Project Data/DEQ/Brooks_2025_profiles_for_comp.csv")
brooks_deq$DATETIME_MST<-as.POSIXct(brooks_deq$DATETIME_MST, format="%m/%d/%Y %H:%M")


# Plot DEQ profile data in red against buoy data
p3 <- ggplot(brooks_wq) + geom_line(aes(DATETIME_MST, WTemp1m, col="1 m"), na.rm=TRUE) +
  geom_line(aes(DATETIME_MST, WTemp4m, col="4 m"), na.rm=TRUE) + 
  geom_line(aes(DATETIME_MST, WTemp7m_corrected, col="7 m"), na.rm=TRUE) + 
  geom_line(aes(DATETIME_UTC, WTemp9m_corrected, col="9 m")) +
  geom_line(aes(DATETIME_MST, WTemp10m_corrected, col="10 m"), na.rm=TRUE) + 
  geom_line(aes(DATETIME_MST, WTemp13m_corrected, col="13 m"), na.rm=TRUE) + 
  geom_line(aes(DATETIME_MST, WTemp15m, col="15 m"), na.rm=TRUE) + 
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp1m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp4m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp7m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp9m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp10m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp13m),shape=21,size=3,fill="red",color="black") +
  geom_point(data=brooks_deq,aes(DATETIME_MST, WTemp15m),shape=21,size=3,fill="red",color="black") +
  xlab("Date") + ylab("Temp (C)") + ylim(4,18) +
  scale_color_manual(name="",values=legend_colors,breaks=c("1 m","4 m","7 m","10 m","13 m","15 m")) +
  theme_bw()
ggplotly(p3) # Looks goog
ggsave(paste0(dir,"/Raw/Time_Corr_Figures/DEQ_temp_comp.png")) # save plot


# Remove uncorrected temperature data
brooks_wq$WTemp7m <- brooks_wq$WTemp7m_corrected
brooks_wq$WTemp9m <- brooks_wq$WTemp9m_corrected
brooks_wq$WTemp10m <- brooks_wq$WTemp10m_corrected
brooks_wq$WTemp13m <- brooks_wq$WTemp13m_corrected

brooks_wq <- brooks_wq %>% select(DATETIME_UTC,DATETIME_MST,WTempC,ConduS,SpConduS,CHLrfu,CHLugL,
                                  BGAPCrfu,BGAPCugL,odomgL,odosat,WTemp1m,WTemp4m,WTemp7m,WTemp9m,WTemp10m,
                                  WTemp13m,WTemp15m,odomgL_15m,odosat_15m,PAR9m)


# Save cleaned datatable
write.csv(brooks_wq,paste0(dir,'/Cleaned/Brooks_WQ_2025_cleaned.csv'))



# Load cleaned data
#brooks_wq <- read.csv(paste0(dir,'/Cleaned/Brooks_WQ_2025_cleaned.csv'))
#brooks_wq$X <- NULL
#brooks_wq$DATETIME_UTC <- as.POSIXct(brooks_wq$DATETIME_UTC,format="%Y-%m-%d %H:%M:%S",tz="UTC")


