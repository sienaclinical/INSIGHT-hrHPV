
# ##########################################################
# PLOT FUNCTION
# ##########################################################
plotZV <- function(V,Z,a1,c1,thetar,title,char){ 
  Y = thetar
  c1_3 <- c1*3  
  V1_0 <- a1 / c1_3
  V1_a <- a1 / c1; Z1_a <- a1*2
  
  maxV <- V1_a * 1.1  
  maxZ <- Z1_a * 1.1
  
  mtitle <- paste(title," with Y = ",round(Y,2)," / ",round(1/Y,2))
  plot(V,Z,axes="FALSE",xlab="",ylab="",xlim=c(0,60),ylim=c(-40,40),pch=char,cex=1)
  axis(1,pos=0)
  axis(2,pos=-maxV*0.05,las=1)
  text(59,40*0.1,"V")
  text(0,40*0.9,"Z")
  lines(c(0,maxV),c(0,0))
  lines(c(0,V1_a),c(a1,Z1_a),lty=1) 
  lines(c(0,V1_0),c(a1,0),lty=2)
  lines(c(V1_0,V1_a),c(0,-Z1_a),lty=1) 
  lines(c(0,V1_0),c(-a1,0),lty=2)
  lines(c(V1_0,V1_a),c(0,Z1_a),lty=1)
  lines(c(0,V1_a),c(-a1,-Z1_a),lty=1)
  return
} 
# ##########################################################


# ##########################################################
# TAU2 ESTIMATIONS
# ##########################################################

estimate.tau2 = function(y, si2, k, method)
{
  #replace si2=0 with small value
  id=which(si2==0)
  if((length(id))>0){si2[id]=0.01}
  
  if(method == "dl") {
    wi0 = 1/si2
    yw0 = sum(wi0 * y) / sum(wi0)
    Q = sum(wi0 * (y - yw0)^2)
    tau2 = max(0, (Q - (k-1))/(sum(wi0) - sum(wi0^2) / sum(wi0)))    
  } 
  else if (method == "dl2"){
    DL = estimate.tau2(y,si2,k,"dl")
    wd = 1/(DL +si2);
    mw = sum(wd*y)/sum(wd)
    numer = (sum(wd*((y-mw)^2))) - (sum(wd*si2) - sum((wd^2)*si2)/sum(wd)) 
    denom = sum(wd) - (sum(wd^2)/sum(wd))
    tau2 = max(0,numer/denom) 
  }
  else if (method == "vc") {
    tau2 = max(0, (var(y) - mean(si2)))
  } 
  else if (method == "mv") {
    if (sum(((y-mean(y))^2)/k) == 0) {mv = 0.01} else {mv = sum(((y-mean(y))^2)/k)}
    ri = si2/mv
    vi = ri + 1
    Y_bar_v = sum(y/vi)/sum((1/vi))
    tau2 = sum((1/vi)*(y - Y_bar_v)^2)/(k-1)    
  } 
  else if (method == "mvvc") {
    tau2ca = max(0, (var(y) - mean(si2)))
    if (tau2ca==0) {tau2ca_MVvc = 0.01} else {tau2ca_MVvc = tau2ca}
    ri_vc = si2/tau2ca_MVvc
    vi_vc = ri_vc + 1
    Y_bar_vc = sum(y/vi_vc)/sum((1/vi_vc))
    tau2 = sum((1/vi_vc)*(y - Y_bar_vc)^2)/(k-1)
  } 
  else if (method == "pm") {
    wi0 = 1/si2
    yw0 = sum(wi0 * y) / sum(wi0)
    Q = sum(wi0 * (y - yw0)^2)
    tau2pr = 0
    if ((Q - (k-1)) <= 0) {tau2 = 0} else
      repeat{      
        wi = 1 / (tau2pr + si2)
        yw = sum(wi * y) / sum(wi)
        Ftau2 = sum(wi * (y - yw)^2) - (k-1)
        if (Ftau2 < 0) {tau2 = 0; break}
        if (abs(Ftau2) < 0.0001) {tau2 = tau2pr; break}
        if (Ftau2 > 0) {
          deltatau2 = Ftau2 / (sum(wi^2 * (y-yw)^2))
          tau2pr = tau2pr + deltatau2
        }
      }
  }
  else {
    tau2.re = 0
    repeat {
      wi <- 1 / (tau2.re + si2)
      teta_w <- sum(wi*y)/sum(wi)
      tau2 <- (sum(wi^2 * ((y-teta_w)^2 - si2)) / sum(wi^2)) + (1/sum(wi))
      if (tau2<0) {tau2=0}
      if (tau2 - tau2.re == 0) {break}
      tau2.conv = (abs(tau2-tau2.re))/ (1+tau2.re)
      if (tau2.conv <= 0.00001) {break}
      if (tau2.conv > 0){tau2.re = tau2}
    }
  }
  
  return(tau2)
  
}
# ################################################################


# ################################################################
# DECISION FUNCTION
# ################################################################
decision = function(a,c,v,z){
  k=dim(v)[1]; p=dim(v)[2]
  dec.mat = matrix(NA,k,p)
  
  for(i in 1:k){
    for(g in 1:p){
      dec=NA 
      vi =sum(v[1:i,g]);  zi =sum(z[1:i,g]); 
      
      # acceptance v limits
      vb0= a/(c*3); vb1=a/c;
      
      if(zi>0) {
        zb.r = a + c*vi
        zb.a = -a+(3*c*vi)
        if(zi>zb.r){dec=1}
        else if(zi<zb.r && vi<vb0){dec=0}
        else if(zi<zb.r && zi>zb.a){dec=0}
        else if(zi<zb.r && zi<zb.a){dec=-1}
        else if(zi<zb.r && vi>vb1){dec=-1}     
      }
      else if(zi<0){
        zb.r = -a - c*vi
        zb.a = a+(3*c*vi)
        if(zi<zb.r){dec=1}
        else if(zi>zb.r && vi<vb0){dec=0}
        else if(zi>zb.r && zi<zb.a){dec=0}
        else if(zi>zb.r && zi>zb.a){dec=-1}
        else if(zi>zb.r && vi>vb1){dec=-1}
      }
      dec.mat[i,g]=dec
    }
  }
  return(dec.mat)
}
# ################################################################


# ####################################################################### #
#                 MATHEW'S CORRELATION COEFICIENT                         #
# ####################################################################### #
mcc = function(model){
  k = length(model);   M = matrix()
  for(i in 1:k){
    a = model[[i]]
    Y = a@y
    Yhat = a@yhat
    c0_true = which(Y==0);       c1_true = which(Y==1)
    c0_hat = which(Yhat==0);     c1_hat = which(Yhat==1)
    tp = length(intersect(c0_true,c0_hat))
    tn = length(intersect(c1_true,c1_hat))
    fp = length(intersect(c1_true,c0_hat))
    fn = length(intersect(c0_true,c1_hat))
    # CORRECTION
    if(tp+fp==0 | tp+fn==0 | tn+fp==0 | tn+fn==0){
      TP=tp+0.5; FP=fp+0.5; TN=tn+0.5; FN=fn+0.5
    }else{
      TP=tp; TN=tn; FP=fp; FN=fn}
    
    M[i] = (TP*TN - FP*FN) / (sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  }
  return(M)
}
# ########################################################################


# ####################################################################### #
#                                 ACCURACY                                #
# ####################################################################### #
acc = function(model,ytrue){
  ncor = 0
  Yhat = as.character(model@yhat);
  n = length(Yhat)
  for(i in 1:n){
    if(Yhat[i]==ytrue[i]){ncor=ncor+1}else{ncor=ncor}}
  ACC = ncor/n
  return(ACC)
}
# ########################################################################



# ########################################################################
#                          FUNCTIONS FOR DA MODEL                        #
# ########################################################################
build.da = function(dataX, dataY, train, dam.k, model.eval){
  da.method = get(paste(dam.k,"CMA",sep=""))
  model.da = classification(X=dataX, y=dataY, learningsets=train, classifier=da.method);
  if(model.eval=="MCC"){
    result = summary(mcc(model.da))[4]  
  }else{ #accuracy
    result = 1-summary(evaluation(model.da, measure="misclassification"))[4]
  } 
  return(result)
}

test.discr= function(evalv.k,dam.k,topk,data,label,train, model.eval){
  if(model.eval=="MCC"){
    ngene = topk[which((evalv.k)==max(evalv.k))[1]];
    mcc.da = build.da(data[,1:ngene], label, train, dam.k, model.eval) 
    result = cbind(mcc.da,ngene)
  }else{
    ngene = topk[which((1-evalv.k)==min(1-evalv.k))[1]];
    acc = build.da(data[,1:ngene], label, train, dam.k, model.eval)
    result = cbind(acc,ngene)
  }
  return(result)
}
# ########################################################################

