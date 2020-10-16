library(ggpmisc)
library(tidyverse)
options(scipen = 999)


spikes <- c("0uM", "0.5uM", "1.0uM", "2.5uM")
vial_types <- c("Durapore", "GF75", "Omnipore", "Vial")

GBT_StdCurves <- read.csv("data_raw/GBT-Fate_THAA_T0-Tlong.csv", stringsAsFactors = FALSE) %>% 
  select(Replicate.Name, Precursor.Ion.Name, Area) %>%
  mutate(runtype = ifelse(str_detect(Replicate.Name, "Blk"), "Blank", "Sample")) %>%
  filter(!str_detect(Replicate.Name, "-1"))

# Subtract the blanks from the areas
GBT_StdCurves_BlksSub <- GBT_StdCurves %>%
  mutate(Replicate.Name = ifelse(runtype == "Blank", 
                                 substr(Replicate.Name, 1, nchar(Replicate.Name)-2), Replicate.Name)) %>%
  group_by(Precursor.Ion.Name, runtype) %>%
  mutate(Area = ifelse(str_detect(Replicate.Name, "Blk"), 
                       round(mean(Area, na.rm = TRUE)), Area)) %>%
  unique() %>%
  ungroup() %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(blankArea = Area[which(runtype == "Blank")]) %>%
  mutate(Area_noBlank = Area - blankArea)

# Isolate Concentrations
GBT_StdCurves_wConcentration <- GBT_StdCurves_BlksSub %>%
  separate(Replicate.Name, into = c("Date", "run", "TempConcentration", "SampID", "Replicate"), sep = "_") %>%
  filter(TempConcentration %in% spikes) %>%
  mutate(TempConcentration = substr(TempConcentration, 1, nchar(TempConcentration)-2),
         TempConcentration = as.numeric(TempConcentration)) %>%
  mutate(Concentration_uM = TempConcentration) %>%
  unite(Date, run, TempConcentration, SampID, Replicate, col = "Replicate.Name")


# plots
ggplot(GBT_StdCurves_wConcentration, aes(x=Concentration_uM, y=Area_noBlank, group = Precursor.Ion.Name)) +
  facet_wrap(~Precursor.Ion.Name) +
  geom_point() + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  stat_poly_eq(aes(label=..eq.label..),
               geom="label", alpha=0.33, formula=(y ~ x),
               label.y = 0.9 * max(GBT_StdCurves_wConcentration$Area_noBlank), 
               label.x = 0.5 * max(GBT_StdCurves_wConcentration$Concentration_uM),
               size = 2.5, parse=TRUE) +
  theme(text=element_text(size=10)) 

# Extract slope values
Slope.Values <- GBT_StdCurves_wConcentration %>% 
  group_by(Precursor.Ion.Name) %>% 
  do({
    mod = lm(Area_noBlank ~ Concentration_uM, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

All.Standards <- GBT_StdCurves_wConcentration %>%
  left_join(Slope.Values) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area_noBlank, Concentration_uM, Intercept, Slope)

Final.Concentrations <- read.csv("data_raw/GBT-Fate_THAA_T0-Tlong.csv") %>%
  left_join(All.Standards %>% select(Precursor.Ion.Name, Slope, Intercept)) %>%
  unique() %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Calculated.Concentration = (Area - Intercept) / Slope)

## Calculation using only the slope 
Final.Concentrations.SlopeCalc <- Final.Concentrations %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Slope, Intercept, Calculated.Concentration) %>%
  filter(!str_detect(Replicate.Name, "A-1|Blk")) %>%
  mutate(Runtype = ifelse(str_detect(Replicate.Name, "A-2"), "StandardCurve", "Sample"),
         Conc.Type = ifelse(Runtype == "Standard"|str_detect(Replicate.Name, "_0uM_"), 0, NA)) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Zero.Concentration.Area = Area[which(Conc.Type == 0)],
         # (x2 - x1) = (y2 - y1) / slope
         NewConc = (Area - Zero.Concentration.Area)/(Slope),
         Conc.Diff = abs(NewConc - Calculated.Concentration))
  
