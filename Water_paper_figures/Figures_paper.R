# Figures for Rowan et al., 2026 

# Load packages
library(tidyverse)
library(scico)
library(lubridate)
library(ggpubr)
#library(Hmisc)
library(stringr)
library(tidyr)
library(broom)


########### Data Loading ###########

# Load data trace elements 
te <- read.csv("data/Rowan_2026_Trace_elements.csv")
te$month_date <- my(te$Month)  # parses date format automatically
# Clean up element names
clean_names <- str_extract(names(te), "(?<=\\.\\.)[A-Za-z]{1,2}(?=\\.\\.)")
names(te) <- ifelse(is.na(clean_names), names(te), clean_names)

# Load data Mg and Sr ratios for Sinclair test
sin <- read.csv("data/Rowan_2026_Mg_Sr_Ca_Sinclair.csv")
sin$month_date <- dmy(sin$Month)  # parses date format automatically

# Load TOC data
toc <- read.csv("data/Rowan_2026_TOC.csv")
toc$month_date <- my(toc$Month)  # parses date format automatically

# Load DOC data
doc <- read.csv("data/Rowan_2026_DOC_radiocarbon.csv")
doc$month_date <- my(doc$Month)  # parses date format automatically
doc <- doc %>% rename(Sample.Name = sample.label)

# Load DIC data
dic <- read.csv("data/Rowan_2026_DIC_radiocarbon.csv")
dic$month_date <- my(dic$Month)  # parses date format automatically

# Load precipitation data
prec <- read.csv("data/Rowan_2026_precipitation_data.csv")
prec$month_date <- my(prec$Month)  # parses date format automatically

# Load cave river discharge data
river <- read.csv("data/Rowan_2026_Milandre_river_discharge.csv")
river$month_date <- dmy(river$Date)  # parses date format automatically



########### CARBON AND HYDROLOGY TIME SERIES ######

###########  Figure 3 sites and averages over time ########### 
# Calculate monthly mean of each group of samples and plot
toc_avg <- toc %>% 
  filter(Corrected.TOC.mg.L != 'is.na') %>% 
  group_by(Group, month_date) %>% 
  summarise(avg = mean(Corrected.TOC.mg.L))

doc_avg <- doc %>% 
  filter(F14C != 'is.na') %>% 
  group_by(Group, month_date) %>% 
  summarise(avg = mean(F14C))

dic_avg <- dic %>% 
  filter(F14C != 'is.na') %>% 
  group_by(Group, month_date) %>% 
  summarise(avg = mean(F14C))

# Here we focus on elements Sr, S and Zn, calculate and normalise to average across all sites

elements <- c("Sr", "S", "Zn")

element_avg <- te %>%
  group_by(Group) %>% #first calculate total average by sample group
  mutate(across(all_of(elements),
                ~ mean(.x, na.rm = TRUE),
                .names = "tot_{.col}")) %>%
  
  group_by(Group, month_date) %>% #now calculate monthly average by sample group
  summarise(
    across(all_of(elements),
           list(avg = ~ mean(.x, na.rm = TRUE),
                sd  = ~ sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    across(starts_with("tot_"), first),
    .groups = "drop")
  
# Normalise for each element
element_avg <- element_avg %>% 
  mutate(Sr_norm = Sr_avg / tot_Sr,
         S_norm = S_avg / tot_S, 
         Zn_norm = Zn_avg / tot_Zn)


# Plot averages and individual timeseries 
fig3_toc <-ggplot(data = toc %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = Corrected.TOC.mg.L, 
                                colour = Group))+
                geom_point(alpha=0.4)+
                geom_line(data = toc_avg %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = avg, 
                                colour = Group), linewidth = 1)+
                xlab("") + ylab("TOC (mg/L)") +  scale_x_date(
                  limits = as.Date(c("2021-12-01", "2024-01-01")),
                  date_breaks = "2 months", date_labels = "%b %Y")+ 
                scale_colour_scico_d(name = 'Sample group', alpha = 1,begin = 0,
                                 end = 1, direction = 1, palette = "batlow") +
                theme_classic() + theme(axis.text.x = element_blank())

fig3_doc_dic <- ggplot(data = doc, aes(x = as.Date(month_date), y = F14C, colour = Group))+
                geom_point(alpha=0.4)+
                geom_line(data = doc_avg, aes(x = as.Date(month_date), y = avg, colour = Group), linewidth = 1)+
                geom_line(data = dic_avg, aes(x = as.Date(month_date), y = avg, colour = Group), 
                          linewidth = 1, linetype = 2)+
                xlab("") + ylab("F14C") + scale_x_date(
                  limits = as.Date(c("2021-12-01", "2024-01-01")),
                  date_breaks = "2 months", date_labels = "%b %Y")+ 
                scale_colour_scico_d(name = 'Sample group', alpha = 1,begin = 0,
                                      end = 1, direction = 1, palette = "batlow") +
                theme_classic() + theme(axis.text.x = element_blank())

#Pivot dataframe for te plot
element_avg_p <- element_avg %>% 
                     pivot_longer(
                       cols = c(Sr_avg, S_avg, Zn_avg, Sr_sd, S_sd, Zn_sd),
                       names_to = c("element", ".value"),
                       names_pattern = "(.*)_(avg|sd)") %>%
                     pivot_longer(
                       cols = starts_with("tot_"),
                       names_to = "element_tot",
                       values_to = "total") %>%
                     mutate(element_tot = sub("tot_", "", element_tot)) %>%
                     filter(element == element_tot) %>%
                     mutate(
                       norm = avg / total,
                       sd_norm = sd / total) %>%
                     select(-element_tot)

# Define a scaling factor 
scale_factor <- 0.5

fig3_te <- ggplot(data = element_avg_p %>% filter(Group == 'Drip'), aes(x = as.Date(month_date),
                       y = ifelse(element == "Zn", norm * scale_factor, norm),colour = element)) +
           geom_line(linewidth = 0.7) + 
           geom_ribbon(aes(ymin = ifelse(element == "Zn",
                    (norm - sd_norm) * scale_factor,
                    norm - sd_norm),
                          ymax = ifelse(element == "Zn",
                    (norm + sd_norm) * scale_factor,
                    norm + sd_norm),fill = element),alpha = 0.2)+ 
                    scale_y_continuous(name = "Concentration (mmol/mol)",
                    sec.axis = sec_axis(~ . / scale_factor, name = "Zn (mmol/mol)"))+ xlab("") +
            scale_x_date(limits = as.Date(c("2021-12-01", "2024-01-01")),
             date_breaks = "2 months", date_labels = "%b %Y") +
            scale_colour_scico_d(name = 'Sample group', palette = "imola") +
            scale_fill_scico_d(name = 'Sample group', palette = "imola") +
            theme_classic() +
            theme(axis.text.x = element_blank())


# Add monthly precipitation and daily river discharge 
fig3_prec <- ggplot(data = prec, aes(x = as.Date(month_date), y = Precipitation))+
                geom_bar(stat="identity", colour = 'lightblue', fill = 'lightblue')+
                xlab("") + ylab("Monthly precipitation (mm)") +
                scale_x_date(
                  limits = as.Date(c("2021-12-01", "2024-01-01")),
                  date_breaks = "2 months", date_labels = "%b %Y")+
                theme_classic() + theme(axis.text.x = element_blank())

fig3_disc <- ggplot(data = river, aes(x = as.Date(month_date), y = Discharge))+
                geom_line(linewidth = 0.7, colour = 'blue2')+
                xlab("Collection date") + ylab("Discharge (L s-1)") +
                scale_x_date(
                  limits = as.Date(c("2021-12-01", "2024-01-01")),
                  date_breaks = "2 months", date_labels = "%b %Y")+ 
                theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange in one plot 
ggpubr::ggarrange(fig3_toc, fig3_doc_dic, fig3_te, fig3_prec, fig3_disc,  
                  ncol=1, nrow = 5, heights = c(1.5,1.5,1.5, 1,1), labels = c('A', 'B', 'C', 'D'),
                  common.legend = FALSE, legend = "bottom", align = "v") 



# Save plot as png
ggsave(filename = file.path("output","Fig3_timeseries.pdf"), width = 12, height = 30, units = "cm")


########### Supplementary Figure 1 with all the sites individually ########### 
sfig1_toc <- ggplot(data = toc %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = Corrected.TOC.mg.L, 
                                                                colour = Site.ID, shape= Group), size=3)+
  geom_point()+
  geom_line(data = toc%>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = Corrected.TOC.mg.L, 
                                                       colour = Site.ID))+
  xlab("Collection date") + ylab("TOC (mg/L)") +
  theme(legend.position="none") + scale_colour_scico_d(name = 'Site ID', alpha = 1,
                                                       begin = 0, end = 1, direction = 1,palette = "batlow") +
  theme_bw() 

sfig1_doc <- ggplot(data = doc, aes(x = as.Date(month_date), y = F14C, colour = Site.ID, shape= Group), size=3)+
  geom_point()+
  geom_line(data = doc, aes(x = as.Date(month_date), y = F14C, colour = Site.ID))+
  xlab("Collection date") + ylab("DOC F14C") +
  theme(legend.position="none") + scale_colour_scico_d(name = 'Site ID', alpha = 1,
                                                       begin = 0, end = 1, direction = 1,palette = "batlow") +
  theme_bw()  

sfig1_dic <- ggplot(data = dic, aes(x = as.Date(month_date), y = F14C, colour = Site.ID, shape= Group), size=3)+
  geom_point()+
  geom_line(data = dic, aes(x = as.Date(month_date), y = F14C, colour = Site.ID))+
  xlab("Collection date") + ylab("DIC F14C") +
  theme(legend.position="none") + scale_colour_scico_d(name = 'Site ID', alpha = 1,
                                                       begin = 0, end = 1, direction = 1,palette = "batlow") +
  theme_bw()  

# Arrange in one plot 
ggpubr::ggarrange(sfig1_toc, sfig1_doc, sfig1_dic,  
                  ncol=1, nrow = 3, heights = c(1,1), widths = c(1,1), 
                  common.legend = TRUE, legend = "bottom", align = "v") 

# Save plot as png
ggsave(filename = file.path("output","SupplFig1_C-timeseries.png"), width = 15, height = 20, units = "cm")

#Calculate average for each drip
drip_avg_doc14C <- doc %>% 
  filter(F14C != 'is.na') %>% 
  group_by(Site.ID) %>% 
  summarise(avg = mean(F14C))

########### Supplementary figure 2 - Trace element concentrations over time and by site ######
Mg <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Mg) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Mg) , colour = Site.ID))+
  xlab("") + ylab("Mg") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Na <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Na) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Na) , colour = Site.ID))+
  xlab("") + ylab("Na") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Al <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Al) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Al) , colour = Site.ID))+
  xlab("") + ylab("Al") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Si <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Si) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Si) , colour = Site.ID))+
  xlab("") + ylab("Si") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw()

P <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(P) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(P) , colour = Site.ID))+
  xlab("") + ylab("P") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw()

S <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(S) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(S) , colour = Site.ID))+
  xlab("") + ylab("S") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

K <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(K) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(K) , colour = Site.ID))+
  xlab("") + ylab("K") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Ca <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Ca) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Ca) , colour = Site.ID))+
  xlab("") + ylab("Ca") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Zn <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Zn) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Zn) , colour = Site.ID))+
  xlab("") + ylab("Zn") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Sr <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Sr) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Sr) , colour = Site.ID))+
  xlab("Date") + ylab("Sr") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Ba <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Ba) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Ba) , colour = Site.ID))+
  xlab("Date") + ylab("Ba") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

Fe <- ggplot(data = te %>% filter(Group != 'Soil') , aes(x = as.Date(month_date), y = as.double(Fe) , colour = Site.ID)) +
  geom_point() +
  geom_line(data = te %>% filter(Group != 'Soil'), aes(x = as.Date(month_date), y = as.double(Fe) , colour = Site.ID))+
  xlab("Date") + ylab("Fe") +
  theme(legend.position="none") + scale_colour_discrete(name = 'Site ID') +
  theme_bw() 

# Arrange in one plot 
ggpubr::ggarrange(Mg, Na, Al, Si, P, S, K, Ca, Zn, Sr, Ba, Fe,  
                  ncol=3, nrow = 4, heights = c(1,1), widths = c(1,1), 
                  common.legend = TRUE, legend = "bottom", align = "v") 

# Save plot as png
ggsave(filename = file.path("output","SupplFig2_TE-concentrations.png"), width = 25, height = 30, units = "cm")


########### SCATTER PLOTS, SINCLAIR TEST ######

########### Figure 4 - TOC vs F14C relationship ######

# Merge TOC and F14C in one dataframe
temp <- merge(toc, doc, by=c("month_date", "Site.ID"))

# Calculate Spearman correlation for the whole dataset
corr_global <- temp %>%
  filter(grepl('Drip_GF|Drip_GR', Group.x)) %>%
  summarise(
    test = list(cor.test(Corrected.TOC.mg.L, F14C,
                         method = "spearman",
                         exact = FALSE))
  ) %>%
  mutate(
    rho = test[[1]]$estimate,
    p = test[[1]]$p.value
  )

label_text <- paste0( #label for plotting later
  "\u03C1 = ", round(corr_global$rho, 2),
  ", p = ", signif(corr_global$p, 2)
)


# TOC vs F14C by site

fig4_corr_site <- temp %>% filter(grepl('Drip_GF|Drip_GR', Group.x)) %>%  
  ggplot(aes(x = Corrected.TOC.mg.L, y = F14C, colour = Site.ID))+
  geom_point()+ #sample points
  geom_smooth(aes(group = Site.ID), method = "lm", #regressions by site  
              se = FALSE,linewidth = 0.4, alpha = 0.3)+
  geom_smooth(method = "lm", colour = "black", linewidth = 1.2,se = FALSE) +  #global regression         
  annotate("text", x = 1.05, y = 0.25, label = label_text, colour = "black", hjust = 0) + #global Spearman correlation     
  xlab("TOC (mg/L)") + ylab("F14C") + 
  scale_colour_scico_d(name = 'Site', alpha = 1,begin = 0,
                       end = 1, direction = 1, palette = "batlow") +
  theme_classic() 


fig4_corr_time <- temp %>% filter(grepl('Drip_GF|Drip_GR', Group.x)) %>%
  mutate(month_date = as.Date(month_date)) %>% 
  ggplot(aes(x = Corrected.TOC.mg.L, y = F14C, colour = factor(month_date)))+
  geom_point()+ #sample points
  geom_smooth(aes(group = factor(month_date)), method = "lm", #regressions by month
              se = FALSE,linewidth = 0.4, alpha = 0.3)+
  geom_smooth(method = "lm", colour = "black", linewidth = 1.2,se = FALSE) +  #global regression   
  annotate("text", x = 1.05, y = 0.1, label = label_text, colour = "black", hjust = 0) + #global Spearman correlation     
  xlab("TOC (mg/L)") + ylab("F14C") + ylim(0, 1.05)+
  scale_colour_scico_d(name = 'Collection month', alpha = 1,begin = 0,
                       end = 1, direction = 1, palette = "batlow",
                       labels = function(x) format(as.Date(x), "%b %y")) +
  theme_classic() 




# Arrange in one plot 
ggpubr::ggarrange(fig4_corr_site, fig4_corr_time,  
                  ncol=1, nrow = 2, heights = c(1,1), labels = c('A', 'B'),
                  common.legend = FALSE, legend = "bottom", align = "v") 



# Save plot as png
ggsave(filename = file.path("output","Fig4_TOC-F14C.pdf"), width = 15, height = 25, units = "cm")


########### Supplementary Fig. 3: Sinclair Test ###########
# Calculate individual correlations to have as a table

# Correlations per month
corr_by_month_spearman <- temp %>%
  filter(grepl('Drip_GF|Drip_GR', Group.x)) %>%
  mutate(month_date = as.Date(month_date)) %>%
  group_by(month_date) %>%
  summarise(
    test = list(cor.test(Corrected.TOC.mg.L, F14C,
                         method = "spearman",
                         exact = FALSE)),
    n = sum(complete.cases(Corrected.TOC.mg.L, F14C))
  ) %>%
  mutate(
    rho = sapply(test, \(x) x$estimate),
    p_value = sapply(test, \(x) x$p.value)
  ) %>%
  select(-test)

corr_by_month_spearman %>%
  mutate(
    month = format(month_date, "%b %y"),
    rho = round(rho, 2),
    p_value = signif(p_value, 2)
  ) %>%
  select(month, rho, p_value, n)


# Correlations per site
corr_by_site_spearman <- temp %>%
  filter(grepl('Drip_GF|Drip_GR', Group.x)) %>%
  group_by(Site.ID) %>%
  summarise(
    test = list(cor.test(Corrected.TOC.mg.L, F14C,
                         method = "spearman",
                         exact = FALSE))  ,
    n = sum(complete.cases(Corrected.TOC.mg.L, F14C))
  ) %>%
  mutate(
    rho = sapply(test, \(x) x$estimate),
    p_value = sapply(test, \(x) x$p.value)
  ) %>%
  select(-test)

corr_by_site_spearman %>%
  mutate(
    site = Site.ID,
    rho = round(rho, 2),
    p_value = signif(p_value, 2)
  ) %>%
  select(site, rho, p_value, n)


########### Sinclair test ######
# Sinclair test by site
fig6_Sinclair_site <- sin %>% filter(grepl('Drip_GF|Drip_GR|RW', Group)) %>%  
  ggplot(aes(x = ln.Mg.Ca., y = ln.Sr.Ca., colour = Site.ID))+
  geom_point()+ #sample points
  geom_smooth(aes(group = Site.ID), method = "lm", #regressions by site  
              se = FALSE,linewidth = 0.4, alpha = 0.3)+
  xlab("ln(Mg/Ca)") + ylab("ln(Sr/Ca)") + 
  scale_colour_scico_d(name = 'Site', alpha = 1,begin = 0,
                       end = 1, direction = 1, palette = "batlow") +
  theme_classic() 

# Save plot as png
ggsave(filename = file.path("output","SupplFig3_Sinclair.pdf"), width = 15, height = 15, units = "cm")


####### Supplementary figure 4 - Figure Sr, Mg vs DO14C by site ########

# Merge TE and DOC data by collection date
d <- merge(doc, te, by=c("month_date", "Sample.Name"))

# Calculate stats for regression
stats_sr <- d %>%
  group_by(Site.ID.x) %>%
  do({
    model <- lm(F14C ~ as.double(Sr), data = .)
    
    data.frame(
      r2 = summary(model)$r.squared,
      p = tidy(model)$p.value[2]  # p-value for Sr
    )
  }) %>%
  ungroup()

# Annotation labels (for plot)
stats_sr <- stats_sr %>%
  mutate(
    label = paste0(
      "R² = ", round(r2, 2),
      "\n p = ", signif(p, 2)))
stats_sr_sig <- stats_sr %>% filter(p < 0.05)

# Positions for plot
stats_sr_sig <- stats_sr_sig %>%
  mutate(x = -Inf, y = Inf)

# Plot for Sr, only show labels for significant correlations
Sr_doc <- ggplot(d, aes(x = as.double(Sr), y = F14C)) +
                geom_point(aes(color = Site.ID.x)) +
                geom_smooth(method = "lm", fill = NA, linewidth = 0.5) +
                facet_wrap(~ Site.ID.x, scales = "free_x") +
                geom_text(data = stats_sr_sig, aes(x = x, y = y, label = label),
                           hjust = -0.1,
                           vjust = 1.1,
                           size = 3,
                           inherit.aes = FALSE) +
                xlab("Sr (mmol/mol)") +
                ylab("F14C DOC") + ylim(0, 1) +
                theme_bw() + theme(legend.position = "none")
Sr_doc 

# Save plot as png
ggsave(filename = file.path("output","SupplFig4_Sr-DOC.pdf"), width = 20, height = 20, units = "cm")


