# Pain_Relief_using_R
Magnetic fields have been shown to have an effect on living tissue as early as the 1930's. Plants have been shown to have an improved growth rate when raised in a magnetic field (Mericle et al., 1964). More recently, doctors and physical therapists have used either static or fluctuating magnetic fields to aid in pain management, most comonly for broken bones. In the case study presented here, Carlos Vallbona and his colleagues sought to answer the question "Can the chronic pain experienced by postpolio patients be relieved by magnetic fields applied directly over an identified pain trigger point?"  Mericle, R. P. et al. "Plant Growth Responses". in : Biological effects of magnetic fields. New York: Plenium Press 1964. p183-195.  Vallbona, Carlos et. al., "Response of pain to static Magnetic fields in postpolio patients, a double blind Pilot study" Archives of Physical Medicine and Rehabilitation. Vol 78, American congress of rehabilitaion medicine, p 1200-1203.
Pain_Relief = read.csv("C:/Users/eleni/Desktop/Pain_Relief.txt", sep="")
View(Pain_Relief)
attach(Pain_Relief)
names(Pain_Relief)
stats.d(Pain_Relief$Score_1)


par(mar=c(4,5,4,5))
qqnorm(Pain_Relief$Score_1)
stats.d(Pain_Relief$Score_1)
t.test(Pain_Relief$Score_1,alternative="two.sided",mu=9,conf.level=0.95)
#το 95% διάστημα εμπιστοσύνης προκύπτει 9.356961-9.803039
t.test(Pain_Relief$Score_2,alternative="two.sided",mu=9,conf.level=0.95)
#το 95% διάστημα εμπιστοσύνης για τον μέσο του δείγματος μετά τη θεραπεία  5.131603 7.028397
t.test(Pain_Relief$Change,alternative = "two.sided",mu=9,conf.level = 0.95)
#το 95% διάστημαεμπιστοσύνης για την διαφορά του πριν και μετά 2.545695 4.454305

plot(Pain_Relief,Pain_Relief$Score_1)
boxplot(Pain_Relief$Active,Pain_Relief$Score_1)

x=Pain_Relief$Score_1
tr.x=log(x);tr.x=bc.t(x,0)## με την bc.t() από MY.LECTURE.FUNCTIONS
stats.d(tr.x)

tr.x1=bc.t(x,0.2)
stats.d(tr.x1)

tr.x2=bc.t(x,0.25)
stats.d(tr.x2)

library(MASS)

boxcox(x~1)## library(MASS)
title('Box-Cox μετασχηματισμός για Scores1')

aux=boxcox(x~1,lambda=seq(0,0.4,1/100),plotit=T)
aux$x[which(aux$y==max(aux$y))]


tr.x3=bc.t(x,0.4)
stats.d(tr.x3)
#απορρίπτεται πάλι η κανονικότητα με την επιλογή της τιμής λ=0.4 

#δημιουργία νεας μεταβλητής για ν

Elite =rep ("No",nrow(Pain_Relief ))
Elite[Pain_Relief$Score_1<8]="Yes"
Elite=as.factor(Elite)
Elite<-factor(Elite,labels=c(0,1))#κωδικοποιούμε την μεταβλητή ελιτ με 0 για όσους είαι >8 και 1 για όσους το score1<8
pain1=data.frame(Pain_Relief,Elite)#προσθέτουμε τη νεα μεταβλητή στο data frame



# δεν ξεχνάμε να μετασχηματίσουμε την Η0 !!!
#mu.0=bc.t(25,0.4);mu.0
#t.test(tr.x,alternative="two.sided",mu=mu.0,conf.level=0.95)
#αυτο δεν το καταλαβα

##### Προσδιορισμός μεγέθους δείγματος για δ=s/4, α=0.05, β=0.20,(power=1-β=0.80)
# Έλεγχος μ=bc.t(25,0.2) vs μ!=bc.t(25,0.2)
## n μέγεθος δείγματος
power.t.test(n=NULL,delta=0.15,sd=sd(tr.x3,na.rm=TRUE),sig.level=0.05,
             power=0.80,type=c("one.sample"),alternative=c("two.sided"))
#Το μέγεθος του δείγματος που απαιτείται είναι 18 άτομα

## δ: μέγεθος επίδρασης
power.t.test(n=length(na.omit(Pain_Relief$Score_1)),delta=NULL,sd=sd(tr.x,na.rm=TRUE),sig.level=0.05,
             power=0.80,type=c("one.sample"),alternative=c("two.sided"))

### ΚΑΜΠΥΛΗ ΙΣΧΥΟΣ
my.delta=seq(0,2.5,length.out=30)
aux=power.t.test(n=length(na.omit(x)),delta=my.delta,sd=sd(na.omit(x)),sig.level=0.05,
                 power=NULL,type=c("one.sample"),alternative=c("two.sided"))
plot(my.delta,aux$power,type='l',main="Καμπύλη Ισχύος",xlab="δέλτα",ylab="ισχύς")

#_________________________________________________________________________________________________________________

## Μη-παραμετρικός, sign test,έλεγχος για διάμεσο
no.sucs=sum(ifelse(na.omit(x)>25,1,0))
no.trials=length(na.omit(x))
binom.test(no.sucs,no.trials)
#απορρίπτεται η μηδενική υπόθεση άρα η κατανονή όχι συμμετρική


#_________________________________________________________________________________________________________________
## Bootstrap δε,  Efron and Tibshirani (1993)
library(bootstrap)
my.mean=function(x){mean(x,na.rm=TRUE)}
set.seed(100)
boot.res=bootstrap(Pain_Relief$Score_1,2000,my.mean)

show.bootstrap.results(Pain_Relief$Score_1,boot.res,my.mean)## από MY.LECTURE.FUNCTIONS
hist(boot.res$thetastar)

stats.d(boot.res$thetastar)
     

set.seed(100)
res.b=bcanon(Pain_Relief$Score_1,nboot=2000,my.mean,alpha=c(0.025,0.05,0.95,0.975))
res.b$confpoints
#το 95% διάστημα εμπιστοσύνης boostrap διαφέρει από το προηγούμενο και είναι 9.34-9.72

my.mean=function(x){mean(x,na.rm=TRUE)}
set.seed(100)
boot.res=bootstrap(Pain_Relief$Score_2,2000,my.mean)

show.bootstrap.results(Pain_Relief$Score_2,boot.res,my.mean)## από MY.LECTURE.FUNCTIONS
hist(boot.res$thetastar)

stats.d(boot.res$thetastar)


set.seed(100)
res.b=bcanon(Pain_Relief$Score_2,nboot=2000,my.mean,alpha=c(0.025,0.05,0.95,0.975))
res.b$confpoints
#_________________________________________________________________________________________________________________
library(gmodels)
CrossTable(Pain_Relief$Active,Pain_Relief$Score_1)

Pain.active<-split(Pain_Relief$Score_1,Pain_Relief$Active,drop=False)#1os tropos να χωρίσουμε το Score1 σε Αctive και Inactive

Origin.F<-Pain_Relief$Active
Origin.F<-factor(Origin.F,labels = c("Active","Inactive"))#2os tropos
Scores1.Group<-split(x,Origin.F)
Scores2.Group<-split(Pain_Relief$Score_2,Origin.F)
stats.d(Scores1.Group[[1]])#από τη συνάρτηση προκύπτουν όλα τα περιγραφικά μέτρα για τους ασθενείς πρίν την ακτινοθεραπεία
stats.d(Scores1.Group[[2]])#και για τους ασθενεις πριν τη θεραπεία της ομάδας ελέγχου
stats.d(Scores2.Group[[1]])
stats.d(Scores2.Group[[2]])




summary(Scores1.Group)
summary(Scores2.Group)

Active=rep("Active",nrow(Pain_Relief))#φτίαχνουμε μια καινουργια στήλη με την μεταβλητη active
Active[Pain_Relief$Active==2]="ακτινοθεραπεία"
Active

library('psych')
describeBy(Pain_Relief$Score_1,Active)
describeBy(Pain_Relief$Score_2,Active)

Change.group<-split(Pain_Relief$Change,Origin.F)
stats.d(Change.group[[1]])
stats.d(Change.group[[2]])

library("ggplot2")
par(mfrow=c(1,2))
qplot(Scores1.Group[[1]])
qplot(Scores1.Group[[2]])

stats.d(Change.group[[1]])
stats.d(Change.group[[2]])

hist(Scores2.Group[[1]])
hist(Scores2.Group[[2]])
par(mfrow=c(1,1))

