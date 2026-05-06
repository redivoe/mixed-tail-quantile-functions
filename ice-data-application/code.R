library(mixqev)

load("preprocessed-data.RData")

y=dati.tutti$ICE_EXT/10^13
nina=ifelse(substr(dati.tutti$NINO_Category,1,2)=="La",1,0)
nino=ifelse(substr(dati.tutti$NINO_Category,1,2)=="El",1,0)
n=length(y)
ti2=((1:n)/60)^2
ti=(1:n)/60
#X=as.matrix(data.frame(,ti2=((1:n)/60)^2,nino=ifelse(substr(dati.tutti$NINO_Category,1,2)=="El",1,0),nina=nina,inter=nina*(1:n)/60))

X=as.matrix(data.frame(ti=ti,ti2=ti2))

thresh=2
y=y[X[,1]>=thresh]
X=X[X[,1]>=thresh,]
X[,1]=X[,1]-thresh
n=length(y)

#fit <- mqgam(y~s(ti, k=20, bs="ad") + nino + nina + ultimo.nino, data = x, qu = c(0.1,0.9))
#summary(fit$fit[[1]])
#summary(fit$fit[[2]])


probs=c(0.9,0.95,0.99,0.999)

out_mixq <- fit_mixq_ls_X(y,X,it_max=15000, which_type=c(1,2,3))
#altri=rep(NA,6)
#altri[1]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(1,2,3))$obj)
#altri[2]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(1,2))$obj)
#altri[3]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(2,3))$obj)
#altri[4]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(1))$obj)
#altri[5]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(2))$obj)
#altri[6]=min(fit_mixq_ls_X(y,X,it_max=15000, which_type=c(3))$obj)

pred=mixq_pred(0.999,X,out_mixq)
pred2=mixq_pred(0.99,X,out_mixq)

# gf
perc=1:99/100
for(i in 1:99) {perc[i]=mean(mixq_pred(perc[i],X,out_mixq)>y)}
pdf("gf.pdf")
plot(1:99/100,perc,xlab="u",ylab="uhat")
abline(a=0,b=1)
dev.off()

but=matrix(NA,500,ncol(X)+6)
butPrev=butPrev2=matrix(NA,500,length(pred))
but[1,]=c(out_mixq$beta,out_mixq$mu,out_mixq$sigma,out_mixq$alpha2,out_mixq$xi2,out_mixq$w[1:2])
butPrev[1,]=pred
butPrev[2,]=pred2
for(j in 2:nrow(but)) {
w=sample(n,n,replace=T)
out=fit_mixq_ls_X(y[w],X[w,],it_max=15000, which_type=c(1,2,3))
but[j,]=c(out$beta,out$mu,out$sigma,out$alpha2,out$xi2,out$w[1:2])
butPrev[j,]=mixq_pred(0.999,X,out)
butPrev2[j,]=mixq_pred(0.99,X,out)
}

ps=rep(NA,ncol(X))

for(j in 1:ncol(X)) {
ps[j]=mixq_bootstrap(out_mixq, y=y, X=X, which_beta=j, B=1000)$pvalue
cat(j, "\t")
}

pdf("pred.pdf")
plot(X[,1],y,type="l",ylim=c(1.5,3.75),xlab="Time",ylab="Ice Extent")
lines(X[,1],pred,col="red")
lines(X[,1],pred2,col="blue")
lines(X[,1],apply(butPrev,2,quantile,0.975),col="red",lty=2)
lines(X[,1],apply(butPrev,2,quantile,0.025),col="red",lty=2)
lines(X[,1],apply(butPrev2,2,quantile,0.975,na.rm=T),col="blue",lty=2)
lines(X[,1],apply(butPrev2,2,quantile,0.025,na.rm=T),col="blue",lty=2)
dev.off()
