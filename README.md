# Pain_Relief_using_R
Magnetic fields have been shown to have an effect on living tissue as early as the 1930's. Plants have been shown to have an improved growth rate when raised in a magnetic field (Mericle et al., 1964). More recently, doctors and physical therapists have used either static or fluctuating magnetic fields to aid in pain management, most comonly for broken bones. In the case study presented here, Carlos Vallbona and his colleagues sought to answer the question "Can the chronic pain experienced by postpolio patients be relieved by magnetic fields applied directly over an identified pain trigger point?"  Mericle, R. P. et al. "Plant Growth Responses". in : Biological effects of magnetic fields. New York: Plenium Press 1964. p183-195.  Vallbona, Carlos et. al., "Response of pain to static Magnetic fields in postpolio patients, a double blind Pilot study" Archives of Physical Medicine and Rehabilitation. Vol 78, American congress of rehabilitaion medicine, p 1200-1203.
Pain_Relief = read.csv("C:/Users/eleni/Desktop/Pain_Relief.txt", sep="")
View(Pain_Relief)
attach(Pain_Relief)
names(Pain_Relief)
stats.d(Pain_Relief$Score_1)

########################################################
#
# R-πακέτα και συναρήσεις για τα μαθήματα στατιστικής
# Κάντε run το αρχείο, και σώστε στο workspace
#
########################################################


### ΠΑΚΕΤΑ ΓΙΑ ΜΑΘΗΜΑΤΑ ΣΤΑΤΙΣΤΙΚΗΣ

my.packages=c("TeachingDemos","MASS","pwr","car","lmtest","TSA","UsingR","SMIR","bootstrap",
"VGAM","class","e1071","lattice")
install.packages(my.packages,repos='http://cloud.r-project.org')


### ΣΥΝΑΡΤΗΣΕΙΣ ΓΙΑ ΜΑΘΗΜΑΤΑ ΣΤΑΤΙΣΤΙΚΗΣ

stats.d=function(x,s.trim=0.05,my.adjust=2,my.factor=2,...)
{# descriptive statistics
## x        : is a numeric vector
## my.factor: is a jitter() argument
## s.trim   : is a mean() argument
## my.adjust: sets the adjust argument in density() to a different default

  options(warn=-1)## no warnings 
	xname=paste(deparse(substitute(x),500),collapse="\n")
#check if data set name has been missed
	if(missing(x)) stop("ΠΡΟΣΟΧΗ! Δεν ορίσατε διάνυσμα δεδομένων")	
	#check the data set to be numeric
	if(!is.numeric(x)) stop("ΠΡΟΣΟΧΗ! Τα δεδομένα δεν είναι αριθμητικό διάνυσμα. Δοκιμάστε as.numeric()")	
	#check the sample size  
	n0=length(x)
	#count the number of NAs
	n=length(na.omit(x));NMs=n0-n
	if(n < 2) stop("ΠΡΟΣΟΧΗ! Λιγότερες από δύο έγκυρες παρατηρήσεις")

aux=paste("      ΠΕΡΙΓΡΑΦΙΚΑ ΣΤΑΤΙΣΤΙΚΑ για ",xname);names(aux)=c("");print(aux,quote=F)
##υπολογίζει μέσο και διορθωμένες για μέσο παρατηρήσεις
	x=na.omit(x)
	mean.x=mean(x)
	c.x=x - mean.x
	## υπολογίζει ροπές γύρω από μέσο 
	x2=c.x^2; x2=sum(x2)
	x3=c.x^3; x3=sum(x3)
	x4=c.x^4; x4=sum(x4)
	##υπολογίζει λειασμένο μέσο 
 	tmean.x=mean(x,trim=s.trim)	
	##υπολογίζει δειγματική διακύμανση και μέτρα μορφής 
	s.var.x=var(x)
	sd.x=sqrt(s.var.x)
	se.mean.x=sd.x/sqrt(n)
	skewness=(x3/n)/(x2/n)^1.5;  se.skewness <- sqrt(6/n)
	kurtosis=((x4/n)/(x2/n)^2)-3;  se.kurtosis <- sqrt(24/n)
	se.s.var.x=sqrt(((kurtosis+2)*s.var.x^2/n))
	if(abs(mean.x/sd.x)==0){coef.var=NaN} else {coef.var=(1+1/(4*n))*sd.x/mean.x}
	min.x=min(x);	max.x=max(x)
	median.x=median(x)
	a=quantile(x,probs=c(0.25,0.75),na.rm=TRUE,names=FALSE)
	Lqua.x=a[1]; Uqua.x=a[2]
	range.x=max.x-min.x
	SW.prob=shapiro.test(x)[[2]]## έλεγχος Shapiro-Wilk 

### plots
	par(oma=c(0,0,2,0))
	layout(matrix(c(1,2,3,4,5,5),nr=3,byrow=TRUE)) ## δημιουργεί 2X3 graphs
	aa=hist(x,plot=FALSE,...)
	x.axis.low=min(mean.x-3*sd.x,(min(x)-0.5*sd.x))
	x.axis.upper=max(mean.x+3*sd.x,min(x)+0.5*sd.x)
	my.xlim=c(x.axis.low,x.axis.upper)	
	my.height=aa$density
	my.ylim=c(0,1.2*max(my.height,dnorm(mean.x,mean.x,sd.x)))

	plot(aa,freq=FALSE,xlim=my.xlim,ylim=my.ylim,main="Ιστόγραμμα και Κανονική Καμπύλη",
	xlab="Τιμές Δείγματος",ylab="Πυκνότητα",col="grey95");box(which="plot",lty ="solid")
	aux.1=seq(x.axis.low,x.axis.upper,length.out=300)
	lines(aux.1,dnorm(aux.1,mean.x,sd.x),lwd=1,lty=1)
	my.jitter=jitter(x,my.factor);	rug(my.jitter)
	abline(h=0,lty=1)
	
	plot(aa,freq=FALSE,xlim=my.xlim,ylim=my.ylim,main="Ιστόγραμμα και Εξομαλυσμένη Καμπύλη",
	xlab="Τιμές Δείγματος",ylab="Πυκνότητα",col="grey95");box(which="plot",lty ="solid")
	par(new=TRUE)
	plot(density(x,adjust=my.adjust,...),type="l",xlim=my.xlim,ylim=my.ylim,ann=FALSE,bty="n")
	rug(my.jitter)
	abline(h=0,lty=1)

	boxplot(x,pars=list(pch=8),main="Θηκόγραμμα",xlab="Τιμές Δείγματος",ylim=my.xlim,horizontal=TRUE,col="grey95")
	s.x <- c.x/sd.x
	qqnorm(s.x,main="qq-plot για Τυπική Κανονική",xlab="Θεωρητικά Ποσοστημόρια",ylab="Δειγματικά Ποσοστημόρια",pch=16,col="grey50")
	qqline(s.x)
	plot(x,main="Δείκτης - Τιμή Δείγματος",xlab="Δείκτης",ylab="Τιμές Δείγματος",pch=16,col="grey50")
	lines(c(-1,n+2),c(mean.x,mean.x),type="l",lty=2)
	title(main=paste("ΓΡΑΦΙΚΑ για",xname),line=-0.6,out=TRUE,cex.main=2.0,lwd=0.8)
	par(mfrow=c(1,1),oma=c(0,0,0,0))## επαναφέρει σε μονό διάγραμμα 

##αποτελέσματα
aux.name=paste(s.trim," λειασμένος μέσος")
Όνομα=c("αριθ. παρατηρήσεων","αριθ. τιμών που λείπουν", "μέσος",aux.name,"διακύμανση","συν.μεταβλητότητας",
"ελάχιστο","μέγιστο","διάμεσος","Q1","Q3","εύρος","ασυμμετρία","κύρτωση(διορθωμένη)","ελ-χος Shapiro-Wilk(p-τιμή)")
τιμή=c(n0,NMs,mean.x,tmean.x,s.var.x,coef.var,min.x,max.x,median.x,Lqua.x,Uqua.x,range.x,skewness,kurtosis,SW.prob)
τιμή=round(τιμή,digits=6)
τυπ.σφάλμα=c(NA,NA,se.mean.x,NA,se.s.var.x,NA,NA,NA,NA,NA,NA,NA,se.skewness,se.kurtosis,NA)
τυπ.σφάλμα=round(τυπ.σφάλμα,digits=6)
d.names.1=rep("",length(τιμή))
d.names.2=c("     ΟΝΟΜΑ","  ΤΙΜΗ","ΤΥΠ.ΣΦΑΛΜΑ")
d.names=list(d.names.1,d.names.2)
aux=cbind(Όνομα,τιμή,τυπ.σφάλμα)
dimnames(aux)=d.names

if(abs(mean.x/se.mean.x)>2){print(aux,na.print="",quote=FALSE)}
else {print(aux,na.print="",quote=FALSE)
aux1='***ΠΡΟΣΟΧΗ: ο μέσος είναι σχεδόν μηδέν, η τιμή  του συν.μεταβλητότητας δεν είναι έγκυρη'
names(aux1)=c("")
print(aux1,quote=FALSE)} }


bc.t=function(x,l,inv=FALSE,...)
{# bc.t performs a Box-Cox trans on the vector x with parameter l
# with inv set to TRUE, performs the inverse transformation
my.epsilon=.Machine$double.neg.eps*1e2
ifelse((x<=0 & inv==FALSE),return('STOP, non-positive values'),{
bc=function(x,l){if(abs(l)<my.epsilon) x=log(x) else x=(x^l-1)/l}
bc.i=function(x,l){if(abs(l)<my.epsilon) x=exp(x) else x=(l*x+1)^(1/l)}
if(inv==FALSE) bc(x,l) else bc.i(x,l)})}



my.distr.plot.2=function(x,factor=NULL,my.amount=NULL,...)
{
## Χρησιμοποιείστε την my.distr.plot.2 για γραφική παράσταση των αποτελεσμάτων της anova
## x is a numeric vector 
## factor, if defined, is a factor object
## my.amount is a jitter() argument, see ?jitter

xname=paste(deparse(substitute(x),500),collapse="\n")
if(is.null(factor)){aux.0=1;aux.2=list(x=x)} else 
{if(is.factor(factor)) {aux.0=length(levels(factor));aux.2=split(x,factor)}
else stop("factor should be either NULL or a factor")}

G.mean.x=mean(x,na.rm=TRUE)

## define auxiliary functions
my.fu.1=function(x)
{mean.x=mean(x,na.rm=TRUE)
sd.x=sd(x,na.rm=TRUE)
max.h=dnorm(mean.x,mean.x,sd.x)
c(mean.x,sd.x,max.h)}

my.fu.2=function(v)
{mean.x=v[1];sd.x=v[2]
x.axis.low=mean.x-3*sd.x
x.axis.upper=mean.x+3*sd.x
x.axis=seq(x.axis.low,x.axis.upper,length=100)
y.axis=dnorm(x.axis,mean.x,sd.x)
lines(x.axis,y.axis,lwd=2,lty=2,col="black")
x1=c(mean.x,mean.x);y1=c(0,dnorm(mean.x,mean.x,sd.x))
points(x1,y1,type='h',lwd=3,lty=1,col="black")}

res.s=matrix(nrow=aux.0,ncol=3)
for(i in 1:aux.0){y=aux.2[[i]];res.s[i,1:3]=my.fu.1(y)}
aa=hist(x,plot=FALSE,...)
my.height=aa$density
max.y=max(c(my.height,res.s[,3]))
min.x=min(res.s[,1]-3*res.s[,2])
max.x=max(res.s[,1]+3*res.s[,2])

## plot  the histogram
my.ylim=c(0,1.2*max.y)
my.xlim=c(min.x-0.15*min.x,max.x+0.15*max.x)
plot(aa,freq=FALSE,ylim=my.ylim,xlim=my.xlim,main="",
col="bisque",xlab=xname,ylab="πυκνότητα",...)
rug(jitter(x,amount=my.amount))
for(i in 1:aux.0){v=res.s[i,1:3];my.fu.2(v)}

points(G.mean.x,0,lwd=2,pch=25,bg="red")
text(x=G.mean.x,y=0.05*max.y,labels="Grand Mean",col="red")

if(aux.0==1) title(main=paste("Κατανομή ",xname),line=-1,out=TRUE,cex.main=1.5,lwd=1.2) else title(main=paste("Κατανομή και μέσοι για ",xname," και ομάδες"),line=-1,out=TRUE,cex.main=1.5,lwd=1.2)}



my.distr.plot.3=function(x,factor=NULL,my.amount=NULL,my.adjust=2,...)
{
## Χρησιμοποιείστε την my.distr.plot.3 για γραφική παράσταση των αποτελεσμάτων της anova
## x is a numeric vector
## factor, if defined, is a factor object
## my.amount is a jitter() argument, see ?jitter
## my.adjust sets the adjust argument in density(),see ?density

xname=paste(deparse(substitute(x),500),collapse="\n")
if(is.null(factor)){k=1;aux.1=list(x=x)} else 
{if(is.factor(factor)) {k=length(levels(factor));aux.1=split(x,factor)}
else stop("factor should be either NULL or a factor object")}

options(warn=-1)

x.matrix=matrix(nrow=k,ncol=512)
y.matrix=matrix(nrow=k,ncol=512)
min.x=numeric(k)
max.x=numeric(k)
max.y=numeric(k)

for(i in 1:k){aux.temp=density(na.omit(aux.1[[i]]),adjust=my.adjust,...)
x.matrix[i,]=aux.temp$x
y.matrix[i,]=aux.temp$y
min.x[i]=min(aux.temp$x)
max.x[i]=max(aux.temp$x)
max.y[i]=max(aux.temp$y)}

G.mean.x=mean(x,na.rm=TRUE)
mean.x=sapply(aux.1,mean,na.rm=TRUE)

## Ιστόγραμμα
hist.plot=hist(x,plot=FALSE,...)
my.height=max(hist.plot$density)

## set plot limits
my.ylim=c(0,1.1*max(max.y,my.height))
min.x=min(na.omit(x),min.x)
min.x=min.x-0.1*abs(min.x)
max.x=max(na.omit(x),max.x)
max.x=max.x-0.1*abs(max.x)
my.xlim=c(min.x,max.x)

## plot histogram
plot(hist.plot,freq=FALSE,xlim=my.xlim,ylim=my.ylim,main="",xlab=xname,ylab="πυκνότητα",col="bisque")
rug(jitter(x,amount=my.amount))

## plot smooth curves and group averages
for (i in (1:k)){
lines(x.matrix[i,],y.matrix[i,],col=grey(i*0.55/k),lwd=4,...)
x1=c(mean.x[i],mean.x[i])
aux.x=min(abs(x.matrix[i,]-mean.x[i]))
index=which(abs(x.matrix[i,]-mean.x[i])==aux.x)
y1=c(0,y.matrix[i,index[1]])
points(x1,y1,type="h",col=grey(i*0.55/k),lwd=3,lty=2)}

## plot Grand mean, etc
points(G.mean.x,0,lwd=2,pch=25,bg="red")
text(x=G.mean.x,y=0.05*max(max.y,my.height),labels="Grand Mean",col="red")

if(k==1) title(main=paste("Κατανομή ",xname),line=-1,out=TRUE,cex.main=1.5,lwd=1.2)
else title(main=paste("Κατανομή και μέσοι για ",xname," και ομάδες"),
line=-1,out=TRUE,cex.main=1.5,lwd=1.2) }



my.plot.pre=function(x,y,res.reg,res.pred,new.data,tr.function=function(x){x},...){
## produces a plot for point and conf.limit predictions for an atomic value
## res.reg: is the result of lm() on a simple linear regression
## of the type y=b.0+b.1*x
## res.pred: is the result of predict.lm() on new.data & res.reg,
## of the form > predict.lm(res.reg,new.data,type="response",interval="prediction")
## tr.function: is a function of the type f(y,...),
## it is the inverse tranformation of that applied on y

xname=paste(deparse(substitute(x),500),collapse="\n")
yname=paste(deparse(substitute(y),500),collapse="\n")

aux1=as.numeric(new.data[,1])
aux2=as.numeric(tr.function(res.pred,...)[,1])
aux3=as.numeric(tr.function(res.pred,...)[,2:3])
plot(aux1,aux2,xlab=xname,ylab=yname,
pch=4,col="red",lwd=2,xlim=c(min(x,aux1,na.rm=TRUE),max(x,aux1,na.rm=TRUE)),
ylim=c(min(y,aux2,aux3,na.rm=TRUE),max(y,aux2,aux3,na.rm=TRUE)))
text(aux1,aux2,labels="P",pos=4,offset=0.6,cex =1.2,col="red")
title(main=paste("Πρόβλεψη της ",yname," με την ",xname),line=-1,out=TRUE,cex.main=1.5,lwd=1.2)
title(sub="Προβλέψεις και 95% Διαστήματα για μια Ειδική Τιμή")
points(x,y)

b.0=as.numeric(res.reg[[1]])[1]
b.1=as.numeric(res.reg[[1]])[2]
data.aug=c(x,aux1)
data.aug.min=min(data.aug,na.rm=TRUE);data.aug.max=max(data.aug,na.rm=TRUE)
x.aug=seq(data.aug.min,data.aug.max,length.out=300)
y=b.0+b.1*x.aug
lines(x.aug,tr.function(y,...),lty=1)

for(i in 1:length(aux1))
{x.axis=rep(aux1[i],2);y=res.pred;y.axis=as.numeric(tr.function(y,...)[i,2:3])
lines(x.axis,y.axis,lwd=2,col="red")}}



plot.cum.1=function(x,y=NULL,ylab1="πιθανότητα",my.main=xname.1,d=0.1,...)
{# plot.cum.1 plots a cumulative density plot
# x is a numeric vector that contains the x-axis values
# y is a numeric vector that contains weights related to x, most usually a d.f
# ylab1 is the default y-axis label
# my.main is the default main title
# d controls the line length before and after min and max x

xname=paste(deparse(substitute(x),500),collapse="\n")
xname.1=paste("Αθροιστική Συνάρτηση Κατανομής της ",xname)

x=na.omit(x)

if(is.null(y)){x=as.numeric(names(table(x)));y=as.numeric(table(x))} else {O.I=order(x)
x=x[O.I];y=y[O.I]}

if(is.null(y)==FALSE & min(y)<=0) stop("Αρνητικές ή Μηδενικές Συχνότητες στο x")

freq.x=cumsum(y)/sum(y)
max.x=max(x);min.x=min(x)
range.x=max.x-min.x
aux.1=d*range.x
low.x=min.x-aux.1;high.x=max.x+aux.1

x.axis=c(low.x,x,high.x)
y.axis=c(0,freq.x,1)
plot(x.axis,y.axis,type="s",xlab=xname,ylab=ylab1,ylim=c(0,1),...)

points(x,freq.x,pch=19,...)
title(main=my.main,line=-2,out=TRUE,...) }


plot.cum.2=function(x,y,delta=700,ylab1="c.d.f.",...)
{#x,y are numerical vectors
# x contains the x-axis values, some ordered figures
# y some sort of non-decreasing weights related to x, most usually a c.d.f
# delta controls the resolution of the plot
# ylab1 is the default y-axis label
max.x=max(x)
min.x=min(x)
n=length(x)
d=numeric(n-1)
points=numeric(n-1)

step=(max.x-min.x)/delta
for(i in 2:n){d[i-1]=x[i]-x[i-1]
		points[i-1]=floor(d[i-1]/step)}

xx=matrix(ncol=max(points),nrow=n-1)
yy=matrix(ncol=max(points),nrow=n-1)
i=1
for(i in 1:(n-1)){j=1
			while(j<=points[i])
			{xx[i,j]=x[i]+(j-1)*step;yy[i,j]=y[i];j=j+1}}

k.max=sum(points)
x.axis=numeric();y.axis=numeric()

# Create some values before x(1)
back=50
for (k in 1:back){x.axis[k]=x[1]-step*(back-k);y.axis[k]=0}
# Note thet k practically starts from back
for(i in 1:(n-1)){for (j in 1:points[i]){k=k+1;x.axis[k]=xx[i,j];y.axis[k]=yy[i,j]}}

# Create some values after x(n)
forth=50
for (kk in 1:forth){k=k+1;x.axis[k]=x[n]+kk*step;y.axis[k]=y[n]}

#plot(x.axis,y.axis,type="l",xlab="x-values",ylab="cum probs")
plot(x.axis,y.axis,pch=".",xlab="x-values",ylab=ylab1,...)
points(x,y,pch=19)

lines(c(x[1],x[1]),c(0,y[1]),...)
for (i in 2:(n)){lines(c(x[i],x[i]),c(y[i-1],y[i]),...)} }




colorie=function (x, y1, y2, ...)
 {###EXAMPLE of use
# x.values=seq(-5,5,length=1000)
# y.values=dnorm(x.values)
# aa=list(x=x.values,y=y.values)
# i=-z.p<aa$x&aa$x<z.p
# colorie(aa$x[i],aa$y[i],rep(0,sum(i)),col=grey(0.9))
polygon( c(x, x[length(x):1]), c(y1, y2[length(y2):1]), ... )}




ts.dec.plot=function(res.dec,my.main=NULL,...)
{## plots the results of decompose()
## makes use of  season() from package TSA
## res.dec: is an object produced by decompose
## my.main: is a character string to be used as title

s.detd=res.dec$x-res.dec$trend
my.frequency=frequency(res.dec$x)
my.length=length(res.dec$x)
season.names=season(res.dec$x)[1:my.frequency]

# Διαγράμματα εποχικότητας
par(mfrow=c(2,2))
#Ραβδόγραμμα για εποχικούς δείκτες
my.width=c(rep(1,my.frequency))
Heights=res.dec$figure
barplot(Heights,my.width,space=0.2,col=grey(0.9),xlab="season",
ylab="seasonal index", main="seasonal index",
axis.lty=1,names.arg=season.names,cex.main=1,ylim=c(1.1*min(Heights,0),1.1*max(Heights,0)))
abline(h=c(0,0),lty=1,col="black",lwd=1)

#boxplots for seasonal factor
season.x=factor(rep(1:my.frequency,len=my.length),labels=season.names)
boxplot(s.detd~season.x,col=grey(0.9),xlab="season",ylab="seasonal component",
main="seasonal component variation",cex.main=1)

#boxplots for random
boxplot(res.dec$random~season.x,xlab="season",ylab="random component",
main="random component variation",col=grey(0.9),cex.main=1)

# plot res.dec$x and trend
ts.plot(res.dec$x,res.dec$trend,lty=c(1,2),main="series and trend",xlab="time",
ylab="series",gpars=list(cex.main=1))

par(mfrow=c(1,1))

title(my.main,outer=TRUE,line=-1,cex.main=1.2)}




show.bootstrap.results=function (obs.x, res.b, my.statistic, ...) 
{   theta.o = my.statistic(obs.x)
    theta.b = mean(res.b$thetastar)
    bias.b = theta.b - theta.o
    se.b = sd(res.b$thetastar)
    min.b = min(res.b$thetastar)
    max.b = max(res.b$thetastar)
    names = c("obs.statistic", "boot.statistic", "bias.b", "s.error.b", 
        "min.b", "max.b")
    values = c(theta.o, theta.b, bias.b, se.b, min.b, max.b)
    show.b = data.frame(matrix(values, ncol = 6))
    colnames(show.b) = names
    rownames(show.b) = c("")
    aa = hist(res.b$thetastar, plot = FALSE)
    i = 2
    while (theta.o > aa$breaks[i]) i = i + 1
    h1 = aa$counts[i - 1]
    while (theta.b > aa$breaks[i]) i = i + 1
    h2 = aa$counts[i - 1]
    plot(aa, main = "bootstrap replications", xlab = "bootstrap values, blue:observed stat, red:boot stat", 
        ylab = "freq")
    lines(c(theta.o, theta.o), c(0, h1), col = "blue", lwd = 2, 
        lty = 3)
    lines(c(theta.b, theta.b), c(0, h2), col = "red", lwd = 2, 
        lty = 6)
    show.b
}



psRsqr=function(y,my.fit,k=2,...)
{### ψευδο-R^2, για μη-γραμμικά μοντέλα
### y: απόκριση, αρχική κλίμακα
### my.fit: πρόβλεψη για απόκριση y, αρχική κλίμακα
### k=2, αριθ ερμηνευτικών μεταβλητών (προκ τιμή=2)
e=y-my.fit
s2.a=sum(e^2)/(length(e)-k)
SSRes=sum(e^2)
SST=sum((y-mean(my.y))^2)
Q_R.2=1-(SSRes/SST)
res=c(sqrt(s2.a),Q_R.2)
names(res)=c('Residual standard error:','pseudo R-squared:')
res}



my.panel.cor=function(x,y,digits=3,prefix="",method="spearman",cex.cor,...){
## put (absolute) correlations on the upper panesl,
## with size proportional to the correlations
##pairs HELP, MODIFIED
usr=par("usr");on.exit(par(usr))
par(usr=c(0,1,0,1))
r=cor(x,y,method=method,use ="pairwise.complete.obs")
txt=format(c(r,0.123456789),digits=digits)[1]
txt=paste(prefix,"ρ=",txt,sep="")
if(missing(cex.cor)) {if(abs(r)<=0.5)cex.cor=3.0 else cex.cor=(2.5*(abs(r)-0.5)+3)}
if(abs(r)<=0.5) aux.1=2.0 else aux.1=10*(abs(r)-0.5)+2.0
text(0.5,0.5,txt,cex=cex.cor,col=grey(sqrt(1/aux.1)))}



my.panel.line=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex=1,col.smooth ="black", span = 2/3, iter = 3,lty.smooth=2,col.line="black",...) 
{    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	my.lm=lm(y~x,...)
	abline(my.lm,col=col.line,...)
    ok <- is.finite(x) & is.finite(y)
    lines(lowess(x[ok],y[ok],f=span,iter=iter),col=col.smooth,lty=lty.smooth,...) }


##pairs HELP  
## put histograms on the diagonal
panel.hist <- function(x,...)
{   usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
	rug(jitter(x))}





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

