################
#
#    hora
#    2021.11
#
################

#boxplot

library(ggplot2)
#
stat3 <- as.data.frame(t(dat['Stat3',]))
colnames(stat3) <- 'Stat3'
stat3$Group <- c(rep('a',4),rep('b',4))

ggplot(stat3,aes(x=Group,y = Stat3,fill=Group))+
  geom_boxplot()+
  geom_jitter()
#
tlr4 <- as.data.frame(t(dat['Tlr4',]))
colnames(stat3) <- 'Tlr4'
tlr4$Group <- c(rep('a',4),rep('b',4))

ggplot(tlr4,aes(x=Group,y = Tlr4,fill=Group))+
  geom_boxplot()+
  geom_jitter()
#
Hmox1 <- as.data.frame(t(dat['Hmox1',]))
colnames(stat3) <- 'Hmox1'
Hmox1$Group <- c(rep('a',4),rep('b',4))

ggplot(Hmox1,aes(x=Group,y = Hmox1,fill=Group))+
  geom_boxplot()+
  geom_jitter()
#
Hif1a <- as.data.frame(t(dat['Hif1a',]))

Hif1a$Group <- c(rep('a',4),rep('b',4))

ggplot(Hif1a,aes(x=Group,y = Hif1a,fill=Group))+
  geom_boxplot()+
  geom_jitter()

#
Cybb <- as.data.frame(t(dat['Cybb',]))
Cybb$Group <- c(rep('a',4),rep('b',4))

ggplot(Cybb,aes(x=Group,y = Cybb,fill=Group))+
  geom_boxplot()+
  geom_jitter()

#miRNA
miR_485 <- as.data.frame(t(GSE19890['hsa-miR-485-5p/ mmu-miR-485/ rno-miR-485',]))
colnames(miR_485) <- 'miR485'
miR_485$Group <- c(rep('b',5),rep('a',5))

ggplot(miR_485,aes(x=Group,y = miR485,fill=Group))+
  geom_boxplot()+
  geom_jitter() 

miR383 <- as.data.frame(t(GSE19890['hsa-miR-383',]))
colnames(miR383) <- 'miR383'
miR383$Group <- c(rep('b',5),rep('a',5))

ggplot(miR383,aes(x=Group,y = miR383,fill=Group))+
  geom_boxplot()+
  geom_jitter()


miR127 <- as.data.frame(t(GSE19890['hsa-miR-127-3p/ mmu-miR-127/ rno-miR-127',]))
colnames(miR127) <- 'miR127'
miR127$Group <- c(rep('b',5),rep('a',5))

ggplot(miR127,aes(x=Group,y = miR127,fill=Group))+
  geom_boxplot()+
  geom_jitter()

miR92 <- as.data.frame(t(GSE19890['hsa-miR-92b*',]))
colnames(miR92) <- 'miR92'
miR92$Group <- c(rep('b',5),rep('a',5))

ggplot(miR92,aes(x=Group,y = miR92,fill=Group))+
  geom_boxplot()+
  geom_jitter()