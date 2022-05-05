library(readxl)
library(ggplot2)
require(dplyr)
require(tidyr)
library(RColorBrewer)
library(scales)
library(ggfortify)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(patchwork)
library(stringr)
library(zoo)

data_path = "~/networkDrives/smb-share:server=s455fm120,share=share/DIS/LABS/Environmental/Environmental_Microbiology/wastewater_covid/"
save_path = "~/networkDrives/smb-share:server=s455fm120,share=share/DIS/LABS/Environmental/Environmental_Microbiology/wastewater_covid/wastewater_plot/"


# WWTP RNA copies/L -------------------------------------------------------

dfw = read.csv(paste0(data_path, 'dep_wastewater.csv'), header = T, stringsAsFactors = F) %>%
  mutate(Date = as.Date(Date,  format=c("%m/%d/%Y"))) %>%
  gather(key = WWTP, value = conc, -Date) %>%
  mutate(WWTP = replace(WWTP, WWTP == "X26W", "26W")) %>%
  mutate(conc = as.numeric(conc)) %>%
  mutate(WWTP = case_when(
    WWTP == '26W' ~ '26th Ward',
    WWTP == 'CI' ~ 'Coney Island',
    WWTP == 'OH' ~ 'Owls Head',
    WWTP == 'RH' ~ 'Red Hook',
    WWTP == 'BB' ~ 'Bowery Bay',
    WWTP == 'JA' ~ 'Jamaica',
    WWTP == 'RK' ~ 'Rockaway',
    WWTP == 'TI' ~ 'Tallman Island',
    WWTP == 'PR' ~ 'Port Richmond',
    WWTP == 'OB' ~ 'Oakwood Beach',
    WWTP == 'NC' ~ 'Newtown Creek',
    WWTP == 'WI' ~ 'Wards Island',
    WWTP == 'NR' ~ 'North River',
    WWTP == 'HP' ~ 'Hunts Point'
  ))

P1 = ggplot(dfw, aes(x = Date, y = conc)) + 
  theme_bw() +
  geom_point(size=1, alpha = 0.7, color='forestgreen') + 
  labs(title = paste0("SARS-CoV2 RNA - Copies/L (",
                      format(min(dfw$Date), "%m/%d/%Y"),
                      " to ",
                      format(max(dfw$Date), "%m/%d/%Y"),
                      ")"), x = "Date", y = "Copies/L") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = 12,
                                                                       angle = 45,
                                                                       hjust = 1)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
        axis.text.y = element_text(hjust = 1, size = 12)) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 12)) +
  theme(strip.text = element_text(size= 12)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%m/%Y") +
  theme(strip.text = element_text(size = 18)) +
  facet_wrap(WWTP ~., nrow = 4)  
P1

ggsave(filename = paste0(save_path, "Wastwater_covid.png"), 
       width = 40, height = 20, units = "cm", plot = P1)


# clinical ----------------------------------------------------------------
dfc = read.csv(paste0(data_path, 'dep_clinical.csv'), header = T, stringsAsFactors = F) %>%
  mutate(Date = as.Date(Date,  format=c("%d/%m/%Y"))) %>%
  gather(key = WWTP, value = conc, -Date) %>%
  mutate(WWTP = replace(WWTP, WWTP == "X26W", "26W")) %>%
  mutate(conc = as.numeric(conc)) %>%
  group_by(WWTP) %>%
  mutate(seven_avg= rollmean(conc, 7,
                             align="right",
                             fill=0)) %>%
  filter(seven_avg != 0) %>% 
  mutate(type = "clinical") %>%
  select(-conc) %>%
  rename(conc = seven_avg) %>%
  ungroup()


dfw2 = read.csv(paste0(data_path, 'dep_wastewater_2.csv'), header = T, stringsAsFactors = F) %>%
  mutate(Date = as.Date(Date,  format=c("%m/%d/%Y"))) %>%
  gather(key = WWTP, value = conc, -Date) %>%
  mutate(WWTP = replace(WWTP, WWTP == "X26W", "26W")) %>%
  mutate(conc = as.numeric(conc)) %>% 
  mutate(conc = conc / 10^7) %>%
  mutate(type = "wastewater")

dfp = rbind(dfc, dfw2) %>%
  mutate(WWTP = case_when(
    WWTP == '26W' ~ '26th Ward',
    WWTP == 'CI' ~ 'Coney Island',
    WWTP == 'OH' ~ 'Owls Head',
    WWTP == 'RH' ~ 'Red Hook',
    WWTP == 'BB' ~ 'Bowery Bay',
    WWTP == 'JA' ~ 'Jamaica',
    WWTP == 'RK' ~ 'Rockaway',
    WWTP == 'TI' ~ 'Tallman Island',
    WWTP == 'PR' ~ 'Port Richmond',
    WWTP == 'OB' ~ 'Oakwood Beach',
    WWTP == 'NC' ~ 'Newtown Creek',
    WWTP == 'WI' ~ 'Wards Island',
    WWTP == 'NR' ~ 'North River',
    WWTP == 'HP' ~ 'Hunts Point'
  ))


P2 = ggplot(dfp, aes(x = Date, y = conc, color=type)) + 
  theme_bw() +
  geom_line(data = filter(dfp, type == "clinical"), 
            size=1, alpha = 0.7) + 
  geom_point(data = filter(dfp, type == "wastewater"),
             size=1.5, alpha = 0.5) + 
  scale_color_manual(values = c("clinical" = "brown", "wastewater" = "forestgreen"),
                     labels = c("New COVID-19 Cases", "SARS-CoV-2 Wastewater")) +
  scale_y_continuous(sec.axis = sec_axis(~.* 10^7, name="N1 GC/day/population")) + 
  labs(title = "", x = "Date", y = "Cases/100,000 people \n (7-day Average)", color = "") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  theme(axis.title.x=element_text(size= 18), axis.text.x = element_text(size = 10,
                                                                        angle = 45,
                                                                        hjust = 1)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
        axis.text.y = element_text(hjust = 1, size = 12),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15),
                                          colour = "forestgreen"),
        axis.text.y.right = element_text(colour = "forestgreen", size = 10),
        axis.title.y.left = element_text(colour = "brown"),
        axis.text.y.left = element_text(colour = "brown", size = 10)) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 12)) +
  theme(strip.text = element_text(size= 13)) +
  guides(colour = guide_legend(override.aes =  list(linetype = c("solid", "blank")))) +
  scale_x_date(date_breaks = "3 month", date_labels =  "%m/%Y") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black')) +
  facet_wrap(WWTP ~., nrow = 4)  
P2

ggsave(filename = paste0(save_path, "Wastwater vs clinical_7_day_average.png"), 
       width = 40, height = 20, units = "cm", plot = P2)














