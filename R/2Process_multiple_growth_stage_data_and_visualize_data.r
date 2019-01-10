##Load required packages
library(tidyverse)
library(broom)
library(plotrix)

##Read into the two datasets, two datasets has the shared Gen id
##dataset1
data1 <- read_tsv("/XXX/data1")
head(data1)
dim(data1)

##some of the 'Gen' has multiple measurements for each of the trait, following code is try to get the mean for each Gen of each trait
data1f <- data1 %>%
  select(Gen,NDVI_06132017, NDVI_06202017, NDVI_07062017, NDVI_07192017, NDVI_08302017) %>%
  gather(Trait, Value, -Gen) %>%
  group_by(Trait, Gen) %>%
  summarise(Avg_Value = mean(Value,na.rm = TRUE)) %>%
  ungroup() %>% 
  spread(key = Trait, value = Avg_Value)

##dataset2 has somemissing values recorded as '.', which can be taken care of when processing the data
data2 <- read_tsv("C:/Users/jinyuw/Box Sync/Project/Highthrougput Image analysis/Image_data_analysis/AmesPanel/Data/Processed/AP_data_long_Traits_1-5_from_Laura.txt") %>%
  mutate(Value = as.numeric(ifelse(Value == ".",NA,Value)))

data2f <- data2 %>%
  filter(Year == 2017 & Trait == 1) %>%
  select(Gen, Trait, Value)  %>%
  group_by(Trait, Gen)  %>%
  summarise(Avg_Value = mean(Value,na.rm = TRUE)) %>%
  ungroup() %>% 
  spread(Trait,Avg_Value,sep = "_")

##Combine the two datasets together and check the pairwise correlation 
combinedtest <- data2f %>%
  arrange(Gen) %>%
  na.omit() %>% 
  left_join(data1f,by = "Gen") %>%
  select(-Gen) %>%
  do(tidy(cor(.)))

datacombined <- data2f %>%
  arrange(Gen) %>%
  na.omit() %>% 
  left_join(data1f,by = "Gen") 

colnames(datacombined) <- c("Gen", "FT", "37", "44", "60", "73", "115")

ngen <- length(datacombined$Gen)

##Add a new variable Genor which orders based on the 'FT' of the data
datacf <- datacombined %>%
  arrange(FT) %>%
  mutate(Genor = c(1:ngen)) %>% 
  gather(Date, value, -Gen, -Genor, -FT) %>% 
  mutate(Date = as.numeric(Date)) %>%
  group_by(Date) %>% 
  mutate(Avg_value = mean(value, na.rm = T)) %>% 
  ungroup() 

##Extract the unique value of each data frame and plot the  
datacf1 <- datacf %>% 
  select(Date, Avg_value) %>% 
  unique()

##Create a gradient color vector

colfunc <- colorRampPalette(c("limegreen","yellow", "orangered" ))
colorsn <- colfunc(ngen)

x_lims <- c(35, 120)
y_lims <- c(95, 145)



tiff("/XXX/Fig1.tiff", width = 8, height = 12, units = 'in', res = 400)
layout(matrix(1:2, ncol=1, byrow = FALSE), widths= 4, heights = c(3,3 ))

##Initiate  a plot
par(mar =c(3, 3.5, 0.5, 0.5))
plot( 0, 0,  xlab = "", ylab = "", ylim = y_lims, xlim = x_lims,  col = "white",  xaxt = "n", yaxt="n")

##Use forloop plot out the growth trend line for each genotpe 
for(i in 1:ngen){
  datap <- datacf[datacf$Genor == i, ]
  points(datap$Date, datap$value, col = colorsn[i], pch = 19)
  lines(datap$Date, datap$value, col = colorsn[i], lwd = 2,  lty = 1, pch = 19)
}

##plot the trend of the population mean
points(datacf1$Date, datacf1$Avg_value, cex = 1.4, pch = 19, col = "blue")
lines(datacf1$Date, datacf1$Avg_value,lty =1, lwd = 2,   col = "blue")

##Add xaxis and yaxis labels
x<-seq(from = round(x_lims[1], digits = 0), to= round(x_lims[2],digits = 0), by= round((x_lims[2]-x_lims[1])/4, digits = 0))
axis(1, mgp=c(3, .1, 0), at=x, labels=x, cex.axis=1.3, tck=.01, las=0)
mtext("Days after planting", side=1, line=1.5, cex=1.5,las=0)

y<- seq(from = round(y_lims[1], digits = 0), to= round(y_lims[2], digits = 0), by= round((y_lims[2]-y_lims[1])/4, digits = 0))
axis(2, mgp=c(3, .1, 0), at=y, labels=y, cex.axis=1.3, tck=.010)
mtext("NDVI (Pixel)", side=2, line=1.5, cex=1.5,las=0)

box(lwd = 2)

########################################################################################################
##plot reversely to review the hidden pattern
##Initiate  a plot
par(mar =c(3, 3.5, 0.5, 0.5))

plot( 0, 0,  xlab = "", ylab = "", ylim = y_lims, xlim = x_lims,  col = "white",  xaxt = "n", yaxt="n")

##Use forloop plot out the growth trend line for each genotpe 
for(i in ngen:1){
  datap <- datacf[datacf$Genor == i, ]
  points(datap$Date, datap$value, col = colorsn[i], pch = 19)
  lines(datap$Date, datap$value, col = colorsn[i], lwd = 2,  lty = 1, pch = 19)
}

##plot the trend of the population mean
points(datacf1$Date, datacf1$Avg_value, cex = 1.4, pch = 19, col = "blue")
lines(datacf1$Date, datacf1$Avg_value,lty =1, lwd = 2,   col = "blue")

##Add xaxis and yaxis labels
x<-seq(from = round(x_lims[1], digits = 0), to= round(x_lims[2],digits = 0), by= round((x_lims[2]-x_lims[1])/4, digits = 0))
axis(1, mgp=c(3, .1, 0), at=x, labels=x, cex.axis=1.3, tck=.01, las=0)
mtext("Days after planting", side=1, line=1.5, cex=1.5,las=0)

y<- seq(from = round(y_lims[1], digits = 0), to= round(y_lims[2], digits = 0), by= round((y_lims[2]-y_lims[1])/4, digits = 0))
axis(2, mgp=c(3, .1, 0), at=y, labels=y, cex.axis=1.3, tck=.010)
mtext("NDVI (Pixel)", side=2, line=1.5, cex=1.5,las=0)

box(lwd = 2)

dev.off()

