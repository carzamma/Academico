######################################################################
###Datos con sobredispersión###
rpoisdisp<-function(lambda=10, 
                    N=1000000, 
                    rep=2, 
                    phi)
{
    if(rep==0)
        stop("El número de réplicas debe 
          ser mayor a cero")
    if(lambda==0)
        stop("El valor de lambda debe 
          ser mayor a cero")
    if(phi==1)
    {
        y<-rpois(N,lambda)
    }
    else
    {
        ni<-N/rep
        k<-lambda/(phi-1)
        y<-rnbinom(n=N, size=k, mu=lambda)
        x<-rep(1:rep, ni)
        mod<-glm(y~factor(x),
                 family=quasipoisson)
    }
    list(Dispersion=summary(mod)$disp, 
         Y=y, X=x)
}
######################################################################

######################################################################
###Bootstrap para el error estándar###

##Bootstrap de Poisson

library(boot)
set.seed(27)

bootpoiss<-function(d,i)
{
    modb<-glm(Y~X, family=poisson, data=d[i,])
    return(summary(modb)$coef[,1])
}

i<-1
eeb0<-NULL
eeb1<-NULL
eepoisb0<-NULL
eepoisb1<-NULL
eequasib0<-NULL
eequasib1<-NULL
phiesp<-NULL
phihat<-NULL
for(n in c(20,40,60,80,100))
{
 for(p in seq(1.1,4, by=0.1))
 {
  #Se generan datos aleatorios    
  dat<-rpoisdisp(N=n, lambda=15, phi=p) 
  #Crea el DF para los modelos
  dat1<-data.frame(Y=dat$Y, X=dat$X)
  #Modelo Poisson clásico
  mod1<-glm(Y~X, family=poisson, data=dat1) 
  #Error estándar beta 0 del mod1
  eepoisb0[i]<-summary(mod1)$coef[1,2]
  #Error estándar beta 1 del mod1
  eepoisb1[i]<-summary(mod1)$coef[2,2]
  #Modelo clásico Quasipoisson
  mod2<-glm(Y~X, family=quasipoisson, data=dat1)
  #Error estándar beta 0 del mod2
  eequasib0[i]<-summary(mod2)$coef[1,2]
  #Error estándar beta 1 del mod2
  eequasib1[i]<-summary(mod2)$coef[2,2]
  #Valor esperado de phi
  phiesp[i]<-p 
  #Valor estimado de Phi
  phihat[i]<-summary(mod2)$disp 
  #Modelo Poisson por Bootstrap
  modboot<-boot(dat1,bootpoiss,1000)
  #Error estándar beta 0 Poisson Bootstrap
  eeb0[i]<-sd(modboot$t[,1]) 
  #Error estándar beta 1 Poisson Bootstrap
  eeb1[i]<-sd(modboot$t[,2]) 
  #Indexador
  i<-i+1       
 }
}

difpoiss<-round(cbind(n=rep(c(20,40,60,80,100), each=30),
                      phiesp,
                      phihat,
                      eeb0,
                      eepoisb0,
                      eequasib0,
                      eeb1,
                      eepoisb1,
                      eequasib1),4)
library(dplyr)
datpoiss<-tbl_df(difpoiss)
datpoiss<-mutate(datpoiss, 
                 deeb0p=abs(eeb0-eepoisb0),
                 deeb0q=abs(eeb0-eequasib0),
                 deeb1p=abs(eeb1-eepoisb1),
                 deeb1q=abs(eeb1-eequasib1))
datpoiss<-datpoiss[,-c(4:9)]

datpoiss<-mutate(datpoiss,
                 difb0=abs(deeb0p-deeb0q),
                 difb1=abs(deeb1p-deeb1q))

##Bootstrap Quasipoisson

library(boot)
set.seed(72)

bootquasi<-function(d,i)
{
    modb<-glm(Y~X, family=quasipoisson, data=d[i,])
    return(summary(modb)$coef[,1])
}

i<-1
eeb0<-NULL
eeb1<-NULL
eepoisb0<-NULL
eepoisb1<-NULL
eequasib0<-NULL
eequasib1<-NULL
phiesp<-NULL
phihat<-NULL
for(n in c(20,40,60,80,100))
{
    for(p in seq(1.1,4, by=0.1))
    {
        #Se generan datos aleatorios    
        dat<-rpoisdisp(N=n, lambda=15, phi=p) 
        #Crea el DF para los modelos
        dat1<-data.frame(Y=dat$Y, X=dat$X)
        #Modelo Poisson clásico
        mod1<-glm(Y~X, family=poisson, data=dat1) 
        #Error estándar beta 0 del mod1
        eepoisb0[i]<-summary(mod1)$coef[1,2]
        #Error estándar beta 1 del mod1
        eepoisb1[i]<-summary(mod1)$coef[2,2]
        #Modelo clásico Quasipoisson
        mod2<-glm(Y~X, family=quasipoisson, data=dat1)
        #Error estándar beta 0 del mod2
        eequasib0[i]<-summary(mod2)$coef[1,2]
        #Error estándar beta 1 del mod2
        eequasib1[i]<-summary(mod2)$coef[2,2]
        #Valor esperado de phi
        phiesp[i]<-p 
        #Valor estimado de Phi
        phihat[i]<-summary(mod2)$disp 
        #Modelo Poisson por Bootstrap
        modboot<-boot(dat1,bootquasi,1000)
        #Error estándar beta 0 Poisson Bootstrap
        eeb0[i]<-sd(modboot$t[,1]) 
        #Error estándar beta 1 Poisson Bootstrap
        eeb1[i]<-sd(modboot$t[,2]) 
        #Indexador
        i<-i+1       
    }
}

difquasi<-round(cbind(n=rep(c(20,40,60,80,100), each=30),
                      phiesp,
                      phihat,
                      eeb0,
                      eepoisb0,
                      eequasib0,
                      eeb1,
                      eepoisb1,
                      eequasib1),4)
library(dplyr)
datquasi<-tbl_df(difquasi)
datquasi<-mutate(datquasi, 
                 deeb0p=abs(eeb0-eepoisb0),
                 deeb0q=abs(eeb0-eequasib0),
                 deeb1p=abs(eeb1-eepoisb1),
                 deeb1q=abs(eeb1-eequasib1))
datquasi<-datquasi[,-c(4:9)]

datquasi<-mutate(datquasi,
                 difb0=abs(deeb0p-deeb0q),
                 difb1=abs(deeb1p-deeb1q))
######################################################################
###Potencias###
potencia<-function(I, muestra, Phi)
{  
 p1<-NULL
 p2<-NULL
 for(i in 1:I)
 {
  g<-rpoisdisp(N=muestra, phi=Phi)
  mod1<-glm(Y~X, 
            family=poisson, 
            data=g)
  mod2<-glm(Y~X, 
            family=quasipoisson, 
            data=g)
  p1[i]<-summary(mod1)$coef[2,4]
  p2[i]<-summary(mod2)$coef[2,4]
 }
 pot1<-sum(1*(p1<0.05))/I
 pot2<-sum(1*(p2<0.05))/I
 potencia<-c(pot1,pot2, pot1-pot2)
 return(c(round(muestra,0),round(g$Dispersion, 2),potencia))
}

set.seed(2772)

potencias<-NULL
i<-1
for(j in c(20,40,60,80,100))
{
 for(k in seq(1.1,4, by=0.1))
 {
  potencias[[i]]<-potencia(100, j, k)
  i<-i+1
 }
}

pot<-data.frame(t(sapply(potencias, round, 2)))
colnames(pot)<-c("Muestra","Phi", "PPoisson", "PQuasi", "Dif")
pot<-tbl_df(pot)
pot
######################################################################


######################################################################

###Tablas###

library(dplyr)

tabla<-function(d,data, cd, c1, c2)
{
    tab<-NULL
    j<-1
    for(i in d)
    {
        dat<-filter(data, data[, cd]<=i)
        dat<-as.matrix(dat)
        tab[[j]]<-t(sapply(tapply(dat[, c1], 
                                  dat[, c2], 
                                  quantile,
                                  probs=c(97.5/100)),round, 2))
        j<-j+1
    }
    names(tab)<-paste(d)
    tab
}

tablas<-tabla(c(.025, .05, .1), pot, 5, 2, 1)
library(xtable)

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{2}{c}{Cuantiles} \\\\\n",
                      "Muestra & 2.5% & 97.5% \\\\\n")
print(xtable(tablas[[3]], 
             caption="Cuantiles del parámetro de dispersión con diferencia de proporciones menor igual a 0.1 para cada tamaño de muestra"), 
             caption.placement = 'top',
             add.to.row = addtorow, include.colnames = FALSE)
######################################################################


######################################################################

###Multiplot###

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
######################################################################
