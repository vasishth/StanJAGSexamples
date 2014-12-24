d <- as.data.frame(read.delim(file="Results/data.txt",sep=" ",col.names=c('subj','expt', 'item','cond','Wpos','word','region','rt')))

 ### look at dataset and understand its structure
#######################################################################
str(d)
summary(d)
nrow(d)

# exclude practice trials (because reading times on them reflect
# probably mostly the training itself
#######################################################################
d <- d[-grep("^practice", d$expt),]
d$logRT <- log(d$rt)
d$Wlen <- nchar(as.character(d$word), type="chars")

abs(scale(unlist(lapply(split(d$logRT, as.factor(as.character(d$subj))), mean)))) < 2.5

## remove non-natives

d <- subset(d,subj!=2 & subj!=19)

first.cut <- function(dat, rt.name='rt', lower=100, upper=2000) {
    rt <- dat[[rt.name]]
    num.too.low <- length(rt[rt < lower])
    num.too.high <- length(rt[rt > upper])
    cat('There are', num.too.low, 'RTs under', lower, 'and',
        num.too.high, 'over', upper, '.\n')
    return(rt < lower | rt > upper)
}


ans <- subset(d,region=="1"|region=="0")
ans$Wpos <- NULL
ans$word <- NULL
ans$ct <- ans$rt
ans$rt <- NULL
ans$logRT <- NULL
ans$Wlen <- NULL
ans$region <- ans$region[drop=TRUE]
ans$correct <- as.numeric(ans$region) - 1
ans$region <- NULL


pos <- function(dat,  subj.name='subj',item.name='item',expt.name='expt') {
    dat <- data.frame(subj = dat[[subj.name]], item = dat[[item.name]],expt=dat[[expt.name]])
     k <- 0
     m <- 0
    ipos <- matrix(0,nrow= nrow(dat), ncol=1, dimnames=list(1:nrow(dat), paste("Lpos")))
    ipos[1,] <- 1
    for(i in 1:nrow(dat)) {
    	   if (dat$subj[i] != k)  {  # if the subject number is different than the last start counter over
		k <- dat$subj[i]
		counter <- 1
		ipos[i,] <- counter }
		  
          else if ((dat$item[i] != dat$item[i-1]) || (dat$expt[i] != dat$expt[i-1]) ) { ## if item number or expt name is different, increment counter
          	 m <- dat$item[i]
	    	counter <- counter + 1
	       ipos[i,] <- counter }
	         
	  else { 
		ipos[i,] <- counter } 
            
}
    return(ipos)
}


d$Lpos <- NA
d$Lpos <- pos(d)

d <- merge(d,ans,by=c('subj','expt','item','cond'))

uans <- subset(ans,expt=="univ")



## remove low accuracy subjects, i.e. < 75% on all trials - not done in final analysis

### d <- subset(d,subj!=8 & subj!=17 & subj!=21 & subj!=26 & subj!=28 & subj!=43 & subj!=45)

## those who scored below 60% on the critical trials - not done in final analysis

### d <- subset(d,subj!=8 & subj!=17 & subj!=21 & subj!=55)


spr <- subset(d,region=="-")
spr <- subset(spr,word!="+")

spr$Wpos <- spr$Wpos[drop=TRUE]
spr$Wpos <- as.character(spr$Wpos)
spr$Wpos <- as.numeric(spr$Wpos)

to.remove <- first.cut(spr,rt.name="rt",upper=5000)
spr$rt[to.remove=='TRUE'] <- NA
spr$logRT[to.remove=='TRUE'] <- NA

cutout <- spr[which(is.na(spr$logRT)), ]

spr <- spr[which(!is.na(spr$logRT)), ]

library(lme4)
l <- lmer(logRT ~ Wlen + log(Lpos) +  (1 | subj), spr)
spr$logRTresidual <- residuals(l)
cutout$logRTresidual <- NA
 
e <- rbind(spr, cutout)
e <- e[order(as.numeric(rownames(e))),]

univ <- subset(e,expt=="univ")
univ$expt <- univ$expt[drop=TRUE]
univ$cond <- univ$cond[drop=TRUE]
univ$Wpos <- univ$Wpos[drop=TRUE]

univ$region <- NULL

univ$region[univ$Wpos==7 & univ$cond=="complex_simple"] <- 0
univ$region[univ$Wpos==8 & univ$cond=="complex_simple"] <- 1
univ$region[univ$Wpos==9 & univ$cond=="complex_simple"] <- 2
univ$region[univ$Wpos==10 & univ$cond=="complex_simple"] <- 3
univ$region[univ$Wpos==11 & univ$cond=="complex_simple"] <- 4
univ$region[univ$Wpos==12 & univ$cond=="complex_simple"] <- 5
univ$region[univ$Wpos==13 & univ$cond=="complex_simple"] <- 6
univ$region[univ$Wpos==14 & univ$cond=="complex_simple"] <- 7
univ$region[univ$Wpos==15 & univ$cond=="complex_simple"] <- 8
univ$region[univ$Wpos==16 & univ$cond=="complex_simple"] <- 9
univ$region[univ$Wpos==17 & univ$cond=="complex_simple"] <- 10

univ$region[univ$Wpos==5 & univ$cond=="simple_simple"] <- 0
univ$region[univ$Wpos==6 & univ$cond=="simple_simple"] <- 1
univ$region[univ$Wpos==7 & univ$cond=="simple_simple"] <- 2
univ$region[univ$Wpos==8 & univ$cond=="simple_simple"] <- 3
univ$region[univ$Wpos==9 & univ$cond=="simple_simple"] <- 4
univ$region[univ$Wpos==10 & univ$cond=="simple_simple"] <- 5
univ$region[univ$Wpos==11 & univ$cond=="simple_simple"] <- 6
univ$region[univ$Wpos==12 & univ$cond=="simple_simple"] <- 7
univ$region[univ$Wpos==13 & univ$cond=="simple_simple"] <- 8
univ$region[univ$Wpos==14 & univ$cond=="simple_simple"] <- 9
univ$region[univ$Wpos==15 & univ$cond=="simple_simple"] <- 10

univ$region[univ$Wpos==9 & univ$cond=="complex_complex"] <- 0
univ$region[univ$Wpos==10 & univ$cond=="complex_complex"] <- 1
univ$region[univ$Wpos==11 & univ$cond=="complex_complex"] <- 2
univ$region[univ$Wpos==12 & univ$cond=="complex_complex"] <- 3
univ$region[univ$Wpos==13 & univ$cond=="complex_complex"] <- 4
univ$region[univ$Wpos==14 & univ$cond=="complex_complex"] <- 5
univ$region[univ$Wpos==15 & univ$cond=="complex_complex"] <- 6
univ$region[univ$Wpos==16 & univ$cond=="complex_complex"] <- 7
univ$region[univ$Wpos==17 & univ$cond=="complex_complex"] <- 8
univ$region[univ$Wpos==18 & univ$cond=="complex_complex"] <- 9
univ$region[univ$Wpos==19 & univ$cond=="complex_complex"] <- 10


univ$region[univ$Wpos==7 & univ$cond=="simple_complex"] <- 0
univ$region[univ$Wpos==8 & univ$cond=="simple_complex"] <- 1
univ$region[univ$Wpos==9 & univ$cond=="simple_complex"] <- 2
univ$region[univ$Wpos==10 & univ$cond=="simple_complex"] <- 3
univ$region[univ$Wpos==11 & univ$cond=="simple_complex"] <- 4
univ$region[univ$Wpos==12 & univ$cond=="simple_complex"] <- 5
univ$region[univ$Wpos==13 & univ$cond=="simple_complex"] <- 6
univ$region[univ$Wpos==14 & univ$cond=="simple_complex"] <- 7
univ$region[univ$Wpos==15 & univ$cond=="simple_complex"] <- 8
univ$region[univ$Wpos==16 & univ$cond=="simple_complex"] <- 9
univ$region[univ$Wpos==17 & univ$cond=="simple_complex"] <- 10

uans <- subset(uans,!is.na(uans$region))

uans$f1[uans$cond=="complex_complex"] <- 1
uans$f1[uans$cond=="complex_simple"] <- 1
uans$f1[uans$cond=="simple_complex"] <- 0
uans$f1[uans$cond=="simple_simple"] <- 0

uans$f2[uans$cond=="complex_complex"] <- 1
uans$f2[uans$cond=="complex_simple"] <- 0
uans$f2[uans$cond=="simple_complex"] <- 1
uans$f2[uans$cond=="simple_simple"] <- 0


univ$region <- as.factor(univ$region)
puniv <- univ

## cutout <- subset(univ,is.na(logRTresidual))
## univ <- subset(univ,!is.na(logRTresidual))
## univ <- rbind(univ,cutout)
## univ <- univ[order(univ$subj,univ$item,univ$f1,univ$f2),]

univ <- transform(univ, subject.index=as.integer(as.factor(univ$subj)))
r7 <- subset(univ,region==7)
r7 <- r7[order(r7$subj,r7$item,r7$f1,r7$f2),]

r8 <- subset(univ, region==8)
r8 <- r8[order(r8$subj,r8$item,r8$f1,r8$f2),]

r9 <- subset(univ, region==9)
r9 <- r9[order(r9$subj,r9$item,r9$f1,r9$f2),]

puniv <- subset(univ,logRTresidual!="NA")
punivc <- subset(puniv,correct==1)


params <- c("beta_f1","beta_f2","beta_f1_X_f2")

buildtable <- function(dat, params) {

	for (i in params) {

		b <- sum(dat[[i]] < 0)/ length(dat[[i]]) ;
		cat(sprintf("%s : %f\n", i, b));
	}

}