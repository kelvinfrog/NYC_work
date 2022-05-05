library(purrr)
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

data_path = "~/networkDrives/smb-share:server=s455fm120,share=share/DIS/LABS/Environmental/Environmental_Microbiology/wastewater_covid/"
save_path = "~/networkDrives/smb-share:server=s455fm120,share=share/DIS/LABS/Environmental/Environmental_Microbiology/wastewater_covid/wastewater_plot/"


# clinical data -----------------------------------------------------------
df = read.csv(paste0(data_path, 'clinical_data_2.csv'), header = F, stringsAsFactors = F)
empty_pos = which(!apply(df == "", 1, all) == FALSE)
ww_name = df[c(1, empty_pos+1), 1]
index = c(0, empty_pos)
len_data = empty_pos[1] -2

table_1 = list()
for (i in c(1:length(index))){
  A = index[i]+1
  B = index[i]+len_data
  db = df[c(A : B),]
  db$name = ww_name[i]
  names(db) <- db[2,]
  db = db[c(-1, -2),-12]
  column_names = names(db)
  column_names[c(1,length(column_names))] = c("Date", "WWTP")
  names(db) = column_names
  db = db %>%
    mutate(WWTP = case_when(
      WWTP == '26th Ward' ~ '26W',
      WWTP == 'Coney Island' ~ 'CI',
      WWTP == "Owl's Head" ~ 'OH',
      WWTP == 'Red Hook' ~ 'RH',
      WWTP == 'Bowery Bay' ~ 'BB',
      WWTP == 'Jamaica' ~ 'JA',
      WWTP == 'Rockaway' ~ 'RK',
      WWTP == 'Tallman Island' ~ 'TI',
      WWTP == 'Port Richmond' ~ 'PR',
      WWTP == 'Oakwood Beach' ~ 'OB',
      WWTP == 'Newtown Creek' ~ 'NC',
      WWTP == "Ward's Island" ~ 'WI',
      WWTP == 'North River' ~ 'NR',
      WWTP == 'Hunts Point' ~ 'HP'
    ))
  
  db =  db %>% 
    mutate_at(vars(Alpha, Beta, Delta, Gamma, Kappa, 
                              Lambda, Mu, `N/A`, Omicron, Other), 
                         as.numeric) %>%
    replace(is.na(.), 0) %>%
    mutate(Other = Other + `N/A`) %>%
    select(-`N/A`) %>%
    rename(`Theta/Mu` = Mu)
  
  table_1[[i]] <- db
  
}

df1 = do.call(rbind, table_1)
df1[c('L452R', 'E484K', 'S477N',
      'WNY1', 'WNY2', 'WNY3',
      'WNY4', 'WNY5', 'Mixed')] <- 0
dfc = df1 %>%
  mutate(Date = do.call("c", lapply(Date, function(x) mean.Date(as.Date(unlist(str_split(x, " - ")), 
                                                                        format=c("%m/%d/%Y"))) ))) %>%
  mutate(Type = "Clinical") %>%
  mutate(row_sum = rowSums(across(where(is.numeric)))) %>%
  modify_if(is.numeric,  `/`, .$row_sum) %>% 
  select(-row_sum)


# WW data -----------------------------------------------------------------
WW = c("NR", "26W", "CI", "OH", "RH", "HP", "BB", "JA", "RK", "TI", "OB", "PR", "NC", "WI")

dfw = read.csv(paste0(data_path, 'wastewater_data.csv'), header = T) %>%
  mutate(Date = as.Date(Date,  format=c("%m/%d/%Y"))) %>%
  mutate(Date = Date - 1) %>%
  select(!c(Flow.rate..MGD., iSeq.MiSeq, SRA.Accession, number_reads, Code, 
            Borough, Virus.Concentration..virus.copies.L.)) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  mutate(WWTP = toupper(WWTP)) %>%
  mutate(WWTP = factor(WWTP, levels = WW)) %>%
  mutate(Type = "Wastewater")

names(dfw) = gsub("Theta", "Theta/Mu", gsub("\\..*", "", names(dfw)))


# Plot data ---------------------------------------------------------------
strain_order = c("Alpha", "Beta", "Gamma", "Delta", "Kappa", "Lambda", "Omicron","Theta/Mu", "L452R", "E484K", "S477N", 
                 "WNY1", "WNY2", "WNY3", "WNY4", "WNY5", "Mixed", "Other")

dfp = dplyr::bind_rows(dfc, dfw) %>%
  gather("strain", "abundance", -c(WWTP, Type, Date)) %>%
  mutate(Borough = case_when(
    WWTP == '26W' ~ 'Brooklyn',
    WWTP == 'CI' ~ 'Brooklyn',
    WWTP == 'OH' ~ 'Brooklyn',
    WWTP == 'RH' ~ 'Brooklyn',
    WWTP == 'BB' ~ 'Queens',
    WWTP == 'JA' ~ 'Queens',
    WWTP == 'RK' ~ 'Queens',
    WWTP == 'TI' ~ 'Queens',
    WWTP == 'PR' ~ 'Staten Is.',
    WWTP == 'OB' ~ 'Staten Is.',
    WWTP == 'NC' ~ 'Man. Qu. & BL',
    WWTP == 'WI' ~ 'Man. & Bronx',
    WWTP == 'NR' ~ 'Manhattan',
    WWTP == 'HP' ~ 'Bronx'
  )) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  mutate(strain = factor(strain, levels = strain_order)) %>%
  mutate(WWTP = factor(WWTP, levels = WW))

ww_min_date = format(min(dfp$Date), "%m/%d/%Y")
ww_max_date = format(max(dfp$Date), "%m/%d/%Y")

P1 = ggplot(dfp, aes(x = Date, y = abundance, color = as.factor(Type), group = as.factor(Type))) + 
  theme_bw() +
  geom_line(size = 1.0, alpha = 0.7) + 
  geom_point(size = 1.0, alpha = 0.7) + 
  scale_color_manual(values=c("#CC6666", "blue")) +
  labs(title = paste0("Clinical and Wastewater SARS-CoV-2 VOC/VOI and Cryptic Variants (", ww_min_date, 
                      " to ", ww_max_date, ")"), x = "Date", y = "Abundance",
       color = "") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = 10,
                                                                       angle = 90,
                                                                       vjust = 0.3)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
        axis.text.y = element_text(hjust = 1, size = 7)) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 10)) +
  theme(strip.text = element_text(size=11)) +
  facet_grid(WWTP~strain)

P1

ggsave(filename = paste0(save_path, 'Clinical vs WW Abundance_all_variants', '.png'),
       width = 40, height = 20, units = "cm", plot = P1)



P2 = ggplot(dfp %>%
              filter(strain %in% c("Alpha", "Beta", "Gamma", "Delta", "Kappa", "Lambda", "Theta/Mu", "Omicron")), 
            aes(x = Date, y = abundance, color = as.factor(Type), group = as.factor(Type))) + 
  theme_bw() +
  geom_line(size = 1.0, alpha = 0.7) + 
  geom_point(size = 1.0, alpha = 0.7) + 
  scale_color_manual(values=c("#CC6666", "blue")) +
  labs(title = paste0("Clinical and Wastewater SARS-CoV-2 Major Variants (", ww_min_date, 
                      " to ", ww_max_date, ")"), x = "Date", y = "Abundance",
       color = "") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = 10,
                                                                       angle = 90,
                                                                       vjust = 0.3)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
        axis.text.y = element_text(hjust = 1, size = 7)) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 10)) +
  theme(strip.text = element_text(size=11)) +
  facet_grid(strain~WWTP)

P2

ggsave(filename = paste0(save_path, 'Clinical vs WW Abundance_major_variants', '.png'),
       width = 40, height = 20, units = "cm", plot = P2)
