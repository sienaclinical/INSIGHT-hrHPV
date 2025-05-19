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
