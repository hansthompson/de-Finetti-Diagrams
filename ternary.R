#install.packages("V:/WORK/Hans/hierfstat_0.04-14.zip", repos = NULL, type = "binary")
library(hierfstat)
library(reshape2)
library(ggtern)
library(ggplot2)
library(dplyr)

x <- seq(0, 1, by = 0.01)
aa <- x^2
AA <- rev(x^2)
Aa <-  -(2 * (seq(-0.5, 0.5, by = 0.01)^2) ) + .5
HWEline <-data.frame(AA, Aa, aa)

files <- list.files(pattern = ".csv")
dat <- read.csv(files[1], skip = 15)


dat <- dat %>% group_by(Assay, Final) %>% summarize( n  =  length(Final))

dat <- dat %>% 
    filter(! Final %in% c("Invalid", "No Call", "NTC")) 

chi_square <- c()
for(i in seq(levels(dat$Assay))) {
    freqs <- dat[dat$Assay == levels(dat$Assay)[i],]
    
    Q <- ((freqs[freqs$Final == "XX",]$n * 2) + freqs[freqs$Final == "YX",]$n) /
        (2 * sum(freqs$n))
    
    P <- 1 - Q 
    
    expAA <- P^2 * sum(freqs$n)
    expAa <- 2 * P * Q * sum(freqs$n)
    expaa <- Q^2 * sum(freqs$n)
    
    
    chi <-     ((freqs[freqs$Final == "YY",]$n - expAA)^2) / expAA  + 
        ((freqs[freqs$Final == "YX",]$n - expAa)^2) / expAa  + 
        ((freqs[freqs$Final == "XX",]$n - expaa)^2) / expaa   
    
    chi_square <- c(chi_square, chi)
    
}

dat <- dat %>% 
    group_by(Assay) %>% 
    mutate(n_percent = n/sum(n))
dat <- dcast(dat, Assay ~ Final)
dat <- dat[complete.cases(dat),]

colors <- rgb(dat$XX, dat$YX, dat$YY)
dat <- cbind(dat, colors)

dat$XX =
    
    p <- ggtern() + tern_limits(1,1,1) + geom_path(data = HWEline, aes(x = AA, y = Aa, z = aa)) 
p + geom_text(data = dat, aes(x = XX, y = YX, z = YY, label = Assay, colour = colors), size = 6) + 
    scale_colour_identity() +
    theme(axis.title=element_text(colour="red"))
