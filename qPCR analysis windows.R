
library(readxl)
library(tidyr)
library(dplyr)
library(TTR)
library(zoo)
library(sciplot)
library(ggplot2)
library(xlsx)



#########################ENTER NUMBERS and NAMES!!!#######################

bio_rep <- 5 #ENTER how many biological replicates
treatments <- 2 #ENTER how many different treatments eg. NC and Prog: 2
treatment_names <- 
  c('SC', 'siRNA') #ENTER treatment names eg. SC, siRNA1...
dose_count <- 1 #ENTER how many different doses used
dose_names <- c('10 nM') #ENTER dose names eg. 10nM, 20nM
excel_path <- 'C:/Users/ismaeldc/OneDrive - Nexus365/Documents/Data/qPCR Analysed/HepG2 POR KD n5.xlsx'
append <- TRUE # False if no Excel for this Analysis, if Excel already exists set to TRUE
figure_title <- 'PPARGC1'
#if one is left empty comment out with: #

#########################CHOOSE GENE TO ANALYSE!!!########################
readIn_goi <- read_excel('C:/Users/ismaeldc/OneDrive - Nexus365/Documents/Data/qPCR Analysed/PPARGC1 HepG2 POR KD n5.xlsx', sheet = 4) #gene of interest Excel

#######################CHOOSE GENE TO NORMALIZE!!!########################
readIn_hkg <- read_excel('C:/Users/ismaeldc/OneDrive - Nexus365/Documents/Data/qPCR Analysed/18S HepG2 POR KD n5.xlsx', sheet = 4) #house keeping gene Excel



#########################################################################




total_samples <- sum(!(readIn_goi[25:132, 'Av CT'] > 40.5))
samples_per_treat_dose_rep <-    total_samples/(treatments*bio_rep*dose_count)
#how many samples per treatment and dose in each replicate
paste('Did you use', samples_per_treat_dose_rep, 'samples per treatment, dose and replicate?')
if (samples_per_treat_dose_rep%%1==0){
  print('Yes? Samples are complete!')
  print('No? check number of treatments and biological replicates entered.')
} else {
  warning('Samples are missing! Check in Excel')
  samples_per_treat_dose_rep <- round(samples_per_treat_dose_rep)
}


#Cleaning Data

###Adding house keeping gene
readIn_goi$'house keeping delta CT' <- readIn_hkg$`delta CT (to min each Panel)`

### Formatting gene of interest

treatment_col <- rep(rep(rep(treatment_names, each =
                               samples_per_treat_dose_rep), bio_rep),dose_count)
treatment_dose <- rep(dose_names, each =
                        samples_per_treat_dose_rep*bio_rep*treatments)

goi_samples <- readIn_goi[25:(24+length(treatment_col)),]
goi_samples$`delta CT` <- goi_samples$`delta CT (to min each Panel)`
cleaned_goi <- goi_samples[, c('probe and plate Name', 'sample',
                               'Av CT', 'delta CT',
                               'house keeping delta CT')]



cleaned_goi$Treatment <- treatment_col
cleaned_goi$Dose <- treatment_dose


### Removing outliers (only works if 3 samples per treatment and dose)

for (i in seq(1, nrow(cleaned_goi), 1)) {
  if (cleaned_goi[i, 'Av CT'] > (colMeans(
    cleaned_goi[1:nrow(cleaned_goi), 'Av CT']) + 5)){
    if (i%%3==1){
      x <- i + 1
      y <- i + 2
    } else if (i%%3==2){
      x <- i - 1
      y <- i + 1
    } else if (i%%3==0){
      x <- i - 1
      y <- i - 2
    } 
    cleaned_goi[i, c('delta CT')] <- as.list(colMeans(
      rbind(cleaned_goi[(x), c('delta CT')],
            cleaned_goi[(y), c('delta CT')])))
    cleaned_goi[i, c('house keeping delta CT')] <- as.list(colMeans(
      rbind(cleaned_goi[(x), c('house keeping delta CT')],
            cleaned_goi[(y), c('house keeping delta CT')])))
    warning('Removing and imputing outlier in row ')
    warning(i)
    warning(' column "sample": ')
    warning(c(cleaned_goi[i, 'sample']))
  } 
}


### Normalizing delta CT of gene of interest
complete_goi <- cleaned_goi
complete_goi$norm_delta_CT <- cleaned_goi$'delta CT' / cleaned_goi$'house keeping delta CT'

print(complete_goi[,2:8], n = nrow(complete_goi))


biological_rep <- rep(rep(1:bio_rep, each = samples_per_treat_dose_rep*treatments), dose_count)

complete_goi$bio_rep <- biological_rep

outlier_check <- complete_goi

goi_av_ct <- complete_goi[seq(samples_per_treat_dose_rep,
                              nrow(complete_goi),
                              samples_per_treat_dose_rep),
                          c('probe and plate Name', 'bio_rep',
                            'Treatment', 'Dose')]

goi_av_ct$mean_delta_CT <- runMean(complete_goi$norm_delta_CT,
                                   samples_per_treat_dose_rep)[seq(
                                     samples_per_treat_dose_rep,
                                     nrow(complete_goi),
                                     samples_per_treat_dose_rep)]

print(goi_av_ct, n = nrow(goi_av_ct))

write.xlsx(goi_av_ct, excel_path, sheetName = figure_title, append = append)



summary_stats <- goi_av_ct %>%
  group_by(Dose, Treatment) %>%
  summarise(mean = mean(mean_delta_CT),
            se = se(mean_delta_CT))

summary_stats

#Look at individual Values if there are outliers -> complete_goi or outlier_check



#Figures
figure_boxplot <- ggplot(goi_av_ct, aes(Treatment, mean_delta_CT)) +
  facet_wrap(~Dose) +
  geom_boxplot() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = 'mRNA levels
(relative expression ratio)',
       title = figure_title)

figure_boxplot_by_treatment <- ggplot(
  goi_av_ct, aes(Dose, mean_delta_CT)) +
  facet_wrap(~Treatment) +
  geom_boxplot() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = 'mRNA levels
(relative expression ratio)',
       title = figure_title)
par(mfrow=c(1,2))
figure_boxplot
figure_boxplot_by_treatment

figure_column <- ggplot(goi_av_ct, aes(Treatment, mean_delta_CT)) +
  facet_wrap(~Dose) +
  geom_col() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = 'mRNA levels
(relative expression ratio)',
       title = figure_title)
#figure_column



#ANOVA
if ((treatments + dose_count) > 3) {
  res.aov <- aov(mean_delta_CT ~ Treatment + Dose, data = goi_av_ct)
  summary(res.aov)
} else if (treatments == 2) {
  t.test(mean_delta_CT ~ Treatment, data = goi_av_ct, paired = TRUE)
} else if (dose_count == 2) {
  t.test(mean_delta_CT ~ Dose, data = goi_av_ct, paired = TRUE)
}


#TukeyHSD
if ((treatments + dose_count) > 3) {
  TukeyHSD(res.aov)
}