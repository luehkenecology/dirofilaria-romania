# set working directory--------------------------------------------------------------------------------
RPROJ <- list(PROJHOME = normalizePath(getwd()))
attach(RPROJ)
rm(RPROJ)
setwd(PROJHOME)

# libaries
library(tidyr)
library(plyr)
library(binGroup)
library(ggplot2)

# read data
data <- read.table (file = "data/data_final_final.csv",
                    header=TRUE,sep=";")

#diro_screen <- read.table (file = "data/screening_new.csv",
#                           header=TRUE,sep=";")
data_2 <- data %>% gather(species, specimens, Culex.modestus:Cx..spec, na.rm = T)


#data_merge <- merge(data_2, diro_screen, by = "Sample.ID", all.x = T)

#data_merge_2 <- na.omit(data_merge)



data_merge_3 <- subset(data_2, !(data_2$species == "Ochlerotatus.dorsalis"))

data_merge_3$immitis <- ifelse(data_merge_3$diro == "D. immitis" | data_merge_3$diro == "D. immitis + D. repens", 
                               1, 0)
data_merge_3$repens <- ifelse(data_merge_3$diro == "D. repens" | data_merge_3$diro == "D. immitis + D. repens", 
                              1, 0)
data_merge_3$both <- ifelse(data_merge_3$diro == "D. immitis" | data_merge_3$diro == "D. immitis + D. repens"|data_merge_3$diro == "D. repens", 
                            1, 0)

data_merge_3$repens[is.na(data_merge_3$repens)] <- 0
data_merge_3$immitis[is.na(data_merge_3$immitis)] <- 0
data_merge_3$both[is.na(data_merge_3$both)] <- 0

FULL_repens <- cbind(data_merge_3[,c(4:8,10:11)], diro = "repens", val = data_merge_3$repens)
FULL_immitis <- cbind(data_merge_3[,c(4:8,10:11)], diro = "immitis", val = data_merge_3$immitis)
FULL_full <- rbind(FULL_repens, FULL_immitis)

source("R/infection_rate.R")

#anzahl <- ddply(FULL_full,.(species),
#                summarize,
#                sum_val1=sum(val),
#                sum_val2=sum(val)/length(individuals)*100)

FULL_full_2 <- subset(FULL_full, specimens > 0)

dateaa <- paste(FULL_full_2$day, FULL_full_2$month, FULL_full_2$year)
FULL_full_2$isoweek <- lubridate::isoweek(as.POSIXct(dateaa, format = "%d%m%Y"))

tztztztzt <- ddply(FULL_full_2,.(diro, species, specimens),
                   summarize,
                   sum=sum(val),
                   freq=length(val))

colnames(tztztztzt) <- c("diro","species","full","sum","freq")
colnames(tztztztzt)
AA <- data.frame(ddply(tztztztzt,.(diro, species),
                       function(x) infection_rate(x$sum,x$full,x$freq)))
AA[,3] <- round(as.numeric(AA[,3])*100, 3)
AA[,4] <- round(as.numeric(AA[,4])*100, 3)
AA[,5] <- round(as.numeric(AA[,5])*100, 3)

#

sum(FULL_full_2$specimens)/2
dateaa <- paste(FULL_full_2$day, FULL_full_2$month, FULL_full_2$year)
FULL_full_2$isoweek <- lubridate::isoweek(as.POSIXct(dateaa, format = "%d%m%Y"))
colnames(FULL_full_2)
dfdfdf <-ddply(FULL_full_2, .(isoweek, TrappingSite), summarize, sum_diro = sum(val), sum_spec = sum(specimens))
plot(dfdfdf[,1], dfdfdf[,2], type = "l")
ggplot(dfdfdf)+
  geom_point(aes(x = isoweek, y = sum_spec)) +
  geom_point(aes(x = isoweek, y = sum_diro*1000))+
  geom_line(aes(x = isoweek, y = sum_spec)) +
  geom_line(aes(x = isoweek, y = sum_diro*1000)) +
  facet_wrap(~TrappingSite, scales = "free_y") +
  theme_bw()
cor.test(dfdfdf[,2], dfdfdf[,3], method = "spearman")  

ee <- ddply(FULL_full_2, .(species),
            summarize,
            sum_val1 = sum(val))

FULL_full_3 <- subset(FULL_full_2,
                      FULL_full_2$species %in% ee[ee[,2] > 0, 1])

summarize_diro <- ddply(FULL_full_3, .(diro, species, isoweek),
                        summarize,
                        sum_val1 = sum(val))

summarize_diro_all <- ddply(FULL_full_3,.(diro, species),
                            summarize,
                            sum_val1 = sum(val))

diro_merge <- merge(summarize_diro, summarize_diro_all, by = c("diro", "species"))

diro_merge$perc <- diro_merge[,4]/diro_merge[,5]*100

gtzu <- ddply(FULL_full_3,.(species),
              summarize,
              sum_val1 = sum(specimens)/2)

phno <- ddply(FULL_full_3, .(species, isoweek), 
              summarize, ind = sum(specimens, na.rm = T)/2)

gl <- merge(phno, gtzu, by = "species")
gl$perc <- gl[,3]/(gl[,4])*100
diro_merge[,1]
gl$imm_perc <- diro_merge$perc[1:130] 
gl$rep_perc <- diro_merge$perc[131:nrow(diro_merge)] 


te <- ggplot(gl) +
  geom_point(aes(x=isoweek,y=perc)) +
  geom_line(aes(x=isoweek,y=perc)) +
  geom_point(aes(x=isoweek,y=imm_perc), colour = "red") +
  geom_line(aes(x=isoweek,y=imm_perc), colour = "red") +
  ylab("percentage") +
  geom_point(aes(x=isoweek,y=rep_perc), colour = "blue") +
  geom_line(aes(x=isoweek,y=rep_perc), colour = "blue") +
  theme_bw() +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  ggtitle("red = perecentage of D. immitis positive pools")


print_resutl <- data.frame(species = AA[1:13,2],
                           specimens = ddply(FULL_full,.(species),
                                             summarize,
                                             sum_val1=sum(individuals))[,2],
                           n_pools = n_pools,
                           n_pos_pools = paste(positive_pools, " (", positive_pools/n_pools*100,")",sep=""),
                           n_pos_repens = paste(positive_pools_species[1:13]," (", positive_pools_species[1:13]/n_pools*100,")",sep=""),
                           infection_rate_diro = paste(AA[1:13,3], " (", AA[1:13,4], "-", AA[1:13,5], ")", sep = ""),
                           n_pos_imm = paste(positive_pools_species[14:26]," (", positive_pools_species[14:26]/n_pools*100,")",sep=""),
                           infection_rate_imm = paste(AA[14:26,3], " (", AA[14:26,4], "-", AA[14:26,5], ")", sep = ""))

write.table(print_resutl,"output/AAA.txt",sep="\t")

#
AA$num <- ddply(FULL_full_2,.(species),
                summarize,
                sum_val1=sum(specimens))[,2]

limits1 <- aes(ymax = as.numeric(AA$V3), ymin=as.numeric(AA$V2))

AA$species <- revalue(AA$species, c("Ae...Oc..spec" = "Aedes spec.",
                                    "Aedes.cinereus" = "Aedes cinereus",
                                    "Aedes.vexans" = "Aedes vexans",
                                    "An..hyrcanus" = "Anopheles hyrcanus",
                                    "Anopheles.algeriensis" = "Anopheles algeriensis",
                                    "Anopheles.maculipennis.s.l." = "Anopheles maculipennis s.l.",
                                    "Coquillettidia.richiardii" = "Coquillettidia richiardii",
                                    "Culex.pipiens.s.l...torrentium"="Culex pipiens s.l./torrentium",
                                    "Culex.modestus" = "Culex modestus",
                                    "Ochlerotatus.caspius" = "Aedes caspius",
                                    "Culiseta.annulata" = "Culiseta annulata",
                                    "Ochlerotatus.detritus" = "Aedes detritus",
                                    "Ochlerotatus.flavescens" = "Aedes flavescens",
                                    "Ochlerotatus.hungaricus" = "Aedes hungaricus",
                                    "unbestimmbar" = "unidentified",
                                    "Cx..spec" = "Cx. spec",
                                    "Uranotaenia.unguiculata" = "Uranotaenia unguiculata"))

#
p2 <- ggplot(AA,aes(x=reorder(species,-num),y=(as.numeric(V1))))+
  geom_point()+
  geom_errorbar(limits1, width=0.2)+
  theme_bw()+
  #theme_linedraw()+
  xlab("")+
  ylab("infection rate")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0,hjust=1))+
  ## facet_wrap(~HOST.GROUP,scales="free")+
  scale_fill_manual(guide = guide_legend(title = "Dirofilaria infection"),values=c("gray75","gray50","black","gray25"))+
  #scale_y_continuous(expand=c(0,0),limits=c(0,105))+
  scale_y_log10()+
  theme(axis.text.x =  element_text(face = "italic"))+
  facet_wrap(~diro, ncol = 1)

png(file = "figure/inf_rate.jpg",width = 8, height=5, units = 'in', res = 1000)
plot(p2)
dev.off()
