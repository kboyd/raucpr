### Copyright (c) 2013 Kendrick Boyd. This is free software. See
### LICENSE for details.

### Precision-Recall Analysis Stuff for R
### 


### Precision-Recall curves
###
### positive (cases) scores (outputs, probabilities, etc.) are random
### variable Y negative (controls) scores are random variable X
###
### pi (prevalence) =  # positives / (# positives + # negatives)


######################################################################
### Binormal
### Assume scores are normally distributed.

### parameters:
### X ~ Normal(0,1)
### Y ~ Normal(pos.mean,pos.sd)


## Calculate true precision for a threshold under binormal distribution.
binormal.precision <- function(c,pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  r = (pi*pnorm(c,pos.mean,pos.sd,lower.tail=FALSE))/(pi*pnorm(c,pos.mean,pos.sd,lower.tail=FALSE) + (1-pi)*pnorm(c,mean=0,sd=1,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from a denominator of 0 should be precision 1
  return(r)
}


## Calculate true recall for a threshold under binormal distribution.
binormal.recall <- function(c,pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  pnorm(c,pos.mean,pos.sd,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under binormal distribution.
binormal.precision.at.recall <- function(recall, pi=0.5,pos.mean=1.0,pos.sd=1.0) {
  return(binormal.precision(qnorm(recall,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd))
}

# Plot the true binormal PR curve.
binormal.plot <- function(pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  plot(function(x) {binormal.precision(qnorm(x,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd)},xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",sub=paste("pi=",pi,", pos.mean=",pos.mean,", pos.sd=",pos.sd))
}


# Calculate true area under binormal PR curve using numeric integration.
binormal.area <- function(pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  f <- function(q) { binormal.precision(qnorm(q,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd) }
  ## monte carlo integration  
  #mean(f(runif(samples)))

  ## integrate method of R, much faster (although might have problems with x=0)
  r = integrate(f,0,1)
  return(r$value)
}


## Generate sample from binormal distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
binormal.sample <- function(n=1, pi=0.5, pos.mean=1.0, pos.sd=1.0, pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)    
  }
  num.neg = n-num.pos

  pos.values = rnorm(num.pos,mean=pos.mean,sd=pos.sd)
  neg.values = rnorm(num.neg,mean=0,sd=1)
  return(list(pos.values=pos.values,neg.values=neg.values))
}




######################################################################
### Bibeta
### Assumes scores are from two Beta distributions.
### Y ~ Beta(pos.a,pos.b)
### X ~ Beta(neg.a,neg.b)


## Calculates true precision for a threshold under bibeta distribution.
bibeta.precision <- function(c,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  ## X ~ Beta(neg.a,neg.b), Y ~ Beta(pos.a,pos.b)
  r = (pi*pbeta(c,pos.a,pos.b,lower.tail=FALSE))/(pi*pbeta(c,pos.a,pos.b,lower.tail=FALSE) + (1-pi)*pbeta(c,neg.a,neg.b,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from a denominator of 0 should have precision of 1
  return(r)
}

## Calculates true recall for a threshold under bibeta distribution.
bibeta.recall <- function(c,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  pbeta(c,pos.a,pos.b,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under bibeta distribution.
bibeta.precision.at.recall <- function(recall,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  return(bibeta.precision(qbeta(recall,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b))
}

## Plot the true bibeta PR curve.
bibeta.plot <- function(pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  plot(function(x) {bibeta.precision(qbeta(x,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b)},
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="Recall",
       ylab="Precision",
       sub=paste("pi=",pi,",X (neg) ~ Beta(",neg.a,",",neg.b,"), Y (pos) ~ Beta(",pos.a,",",pos.b,")")
       )
  
}

## Calculate the true area under the bibeta PR curve using numeric integration.
bibeta.area <- function(pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  f <- function(q) { bibeta.precision(qbeta(q,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b) }

  ## use R's integrate method
  res = integrate(f,0,1)
  return(res$value)
}

## Generate sample from bibeta distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
bibeta.sample <- function(n=1,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1,pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)
  }
  
  num.neg = n-num.pos

  ## Generates a sample of observed values for neg.values (X) ~
  ## Beta(neg.a,neg.b) and pos.values (Y) ~ Beta(pos.a,pos.b).

  pos.values = rbeta(num.pos,pos.a,pos.b)
  neg.values = rbeta(num.neg,neg.a,neg.b)

  return(list(pos.values=pos.values,neg.values=neg.values))
}


######################################################################
### Offset Uniform
### Assumes the scores are drawn from two uniform distributions
### spanning different ranges

### X ~ Uniform(0,1)
### Y ~ Uniform(a,a+b)


## Calculate true precision for a threshold under offset uniform distributions.
offsetuniform.precision <- function(c,pi=0.5,a=0.5,b=1) {
  r = (pi*punif(c,min=a,max=a+b,lower.tail=FALSE))/(pi*punif(c,min=a,max=a+b,lower.tail=FALSE) + (1-pi)*punif(c,min=0,max=1,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from denominator of 0 should have precision of 1
  return(r)
}

## Calculate true recall for a threshold under offset uniform distributions.
offsetuniform.recall <- function(c,pi=0.5,a=0.5,b=1) {
  punif(c,min=a,max=a+b,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under offset uniform distribution.
offsetuniform.precision.at.recall <- function(recall, pi=0.5,a=0.5,b=1) {
  return(offsetuniform.precision(qunif(recall,a,a+b,lower.tail=FALSE),pi=pi,a=a,b=b))
}

## Plot the true offset uniform PR curve.
offsetuniform.plot <- function(pi=0.5,a=0.5,b=1) {
  plot(function(x) {offsetuniform.precision(qunif(x,min=a,max=a+b,lower.tail=FALSE),
                                            pi=pi,
                                            a=a,
                                            b=b)},
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="Recall",
       ylab="Precision",
       sub=paste("pi=",pi," a=",a,"b=",b))  
}

## Calculate true area under offset uniform PR curve using numeric integration.
offsetuniform.area <- function(pi=0.5,a=0.5,b=1) {
  f <- function(q) { offsetuniform.precision(qunif(q,min=a,max=a+b,lower.tail=FALSE),
                                             pi=pi,
                                             a=a,
                                             b=b)}
  r = integrate(f,0,1)
  return (r$value)
}

## Generate sample from offset uniform distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
offsetuniform.sample <- function(n=1,pi=0.5,a=0.5,b=1,pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)
  }
  num.neg = n-num.pos

  ## X (neg) ~ Uniform(0,1)
  ## Y (pos) ~ Uniform(a,a+b)

  pos.values = runif(n=num.pos,min=a,max=a+b)
  neg.values = runif(n=num.neg,min=0,max=1)

  return(list(pos.values=pos.values,neg.values=neg.values))  
}




aucpr.conf.int <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,method="binomial",bootstrap.replicates=1000) {
  ## Calculates confidence interval for an AUCPR estimate.
  ## method=c("binomial","expit")
  m = match.arg(method,c("binomial","expit","bootstrap"))
  
  if (m=="binomial") {
    return(aucpr.conf.int.binomial(estimate,num.pos=length(pos.values),num.neg=length(neg.values),conf.level=conf.level))
  }
  else if (m=="expit") {
    return(aucpr.conf.int.expit(estimate,num.pos=length(pos.values),num.neg=length(neg.values),conf.level=conf.level))
  }
  else if (m=="bootstrap") {
    return(aucpr.conf.int.bootstrap(estimate,pos.values,neg.values,prcurve.estimator,conf.level,bootstrap.replicates))
  }
}

aucpr.conf.int.binomial <- function(estimate,num.pos,num.neg,conf.level=0.95) {
  ## Calculates confidence interval for an AUCPR estimate under
  ## binomial assumptions. Uses num.pos as the sample size (might not
  ## be best?).
  
  ci = estimate+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*sqrt(estimate*(1-estimate)/num.pos)
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "binomial"
  return(ci)
}

aucpr.conf.int.expit <- function(estimate,num.pos,num.neg,conf.level=0.95) {
  ## Calculates confidence interval for an AUCPR estimate using expit.
  
  ## convert to logit scale
  est.logit = log(estimate/(1-estimate))
  ## standard error (from Kevin Eng)
  se.logit = sqrt(estimate*(1-estimate)/num.pos)*(1/estimate + 1/(1-estimate))
  ## confidence interval in logit
  ci.logit = est.logit+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*se.logit

  ## back to original scale
  ci = exp(ci.logit)/(1+exp(ci.logit))
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "expit"
  return(ci)
}

aucpr.conf.int.bootstrap <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,replicates=1000) {
  areas = rep(0,replicates)

  for (i in 1:replicates) {
    p.v = sample(pos.values,replace=TRUE)
    n.v = sample(neg.values,replace=TRUE)

    res = prcurve.estimator(p.v,n.v)
    areas[i] = res$area
  }

  q = quantile(areas,c((1-conf.level)/2,(1+conf.level)/2))
  
  ci = c(q[1],q[2])
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "bootstrap"
  attr(ci,"median") = median(areas)
  attr(ci,"mean") = mean(areas)
  
  return (ci)
}

aucpr.conf.int.crossvalidation <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,folds=10) {
  areas = rep(0,folds)

  pos = sample(pos.values)
  neg = sample(neg.values)

  pos.counts = rep(c(as.integer(length(pos.values)/folds)+1,as.integer(length(pos.values)/folds)),c(length(pos.values)%%folds,10-length(pos.values)%%folds))
  neg.counts = rep(c(as.integer(length(neg.values)/folds)+1,as.integer(length(neg.values)/folds)),c(length(neg.values)%%folds,10-length(neg.values)%%folds))
  
  pos.index = 1
  neg.index = 1
  for (k in 1:folds) {
    p = pos[pos.index:(pos.index+pos.counts[k]-1)]
    n = neg[neg.index:(neg.index+neg.counts[k]-1)]

    areas[k] = prcurve.estimator(p,n)$area
    
    pos.index = pos.index + pos.counts[k]
    neg.index = neg.index + neg.counts[k]
  }
    

  ## normal approximation
  ##ci = mean(areas) + qnorm(c((1-conf.level)/2,(1+conf.level)/2))*sd(areas)/sqrt(length(areas))
  ## use t-distribution
  ci = mean(areas) + qt(c((1-conf.level)/2,(1+conf.level)/2),df=length(areas)-1)*sd(areas)/sqrt(length(areas))
  
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "crossvalidation, t-dist"
  attr(ci,"values") = areas
  attr(ci,"mean") = mean(areas)
  
  return(ci)
}





# return confusion matrix info for >= threshold
# assume labels are 0 and 1
confusion.matrix <- function(values,labels,threshold)
{
  counts=table(factor(values>=threshold,c(FALSE,TRUE)),labels)
  ## print(counts)
  
  TP = counts[2,2]
  TN = counts[1,1]
  FP = counts[2,1]
  FN = counts[1,2]
  
  if (TP+FP == 0)
    return(list(tp=TP,tn=TN,fp=FP,fn=FN,precision=NaN,sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),accuracy=(TP+TN)/(TP+TN+FP+TN)))
  else	
    return(c(tp=TP,tn=TN,fp=FP,fn=FP,precision=TP/(TP+FP),sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),accuracy=(TP+TN)/(TP+TN+FP+TN)))
}


## Create all confusion matrices from set of pos.values and neg.values
make.confusion.matrices <- function(pos.values,neg.values) {
  ## thresholds, sorting to have increasing recall
  thresholds = sort(c(-Inf,unique(c(pos.values,neg.values)),Inf),decreasing=TRUE)
  
  ## 0 - negative (control)
  ## 1 - positive (case)
  labels = c(rep(1,length(pos.values)),rep(0,length(neg.values)))
  values = c(pos.values,neg.values)
  
  # create set of confusion matrices
  d=sapply(thresholds,function(t) { confusion.matrix(values,labels,t) })

  return(d)
}

precisions.recalls.fast <- function(pos.values,neg.values) {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]

  tp = 0
  fp = 0

  z = 1
  precisions = rep(0,length(indices)+1)
  recalls = rep(0,length(indices)+1)

  tp = 0
  fp = 0

  for (i in 1:length(indices)) {
    if (labels[i]==1) {
      tp = tp + 1      
    }
    else {
      fp = fp + 1
    }

    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can put a threshold

      ## if first, insert (r=0,p=1) point
      if (z==1) {
        recalls[z] = 0
        precisions[z] = 1
        z = z + 1
      }
      recalls[z] = tp/(length(pos.values))
      precisions[z] = tp/(tp+fp)
      z = z + 1
    }
  }

  ## truncate recals and precisions
  recalls = recalls[1:(z-1)]
  precisions = precisions[1:(z-1)]
  
  return (list(recalls=recalls,precisions=precisions))
}


# return c(v1[1],v2[1],v1[2],v2[2],...)  works with uneven vectors,
# appending the extra elements from the longer vector so return always
# has length of length(v1)+length(v2) from
# http://tolstoy.newcastle.edu.au/R/help/06/03/22717.html
interleave <- function(v1,v2) {
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}


#################################################################
## Estimators
#################################################################
# Under, connects lowest precision on left side to highest precision
# on the right
prcurve.lowertrap <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions
  

  # get max and min precision for each unique recall
  r.unique = unique(recalls)
  p.min = sapply(r.unique,function(r) { min(precisions[recalls==r])})
  p.max = sapply(r.unique,function(r) { max(precisions[recalls==r])})

  # use min on left side and max on right side of each area between known recalls
  ids=2:length(r.unique)
  area = as.double((r.unique[ids]-r.unique[ids-1]) %*% (p.max[ids] + p.min[ids-1]) / 2)

  rs = rep(r.unique,each=2)
  ps = interleave(p.max,p.min)

  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.lowertrap,
    conf.level=conf.level,
    method=conf.int.method)

  return(list(area=area,
              x=rs,
              y=ps,
              conf.int=ci))  
}

# connects highest precisions at each recall, definitely an
# over-estimate
prcurve.uppertrap <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions

  # get max and min precision for each unique recall
  r.unique = unique(recalls)
  p.min = sapply(r.unique,function(r) { min(precisions[recalls==r])})
  p.max = sapply(r.unique,function(r) { max(precisions[recalls==r])})

  # use max precision for all recalls (NOT LEGITIMATE!)
  
  ids=2:length(r.unique)
  area = as.double((r.unique[ids]-r.unique[ids-1]) %*% (p.max[ids] + p.max[ids-1]) / 2)


  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.uppertrap,
    conf.level=conf.level,
    method=conf.int.method)
  

  return(list(area=area,
              x=r.unique,
              y=p.max,
              conf.int=ci))
}


## estimate PR curve and area under PR curve using average precision
## method, mean precision at each positive example
prcurve.ap.slow <- function(pos.values, neg.values, conf.level=0.95, conf.int.method="binomial") {
  ## thresholds, sorting to have increasing recall
  ## only use thresholds for positives
  thresholds = sort(pos.values,decreasing=TRUE)
  
  labels = c(rep(1,length(pos.values)),rep(0,length(neg.values)))
  values = c(pos.values,neg.values)
  
  ## create set of confusion matrices
  d=sapply(thresholds,function(t) { confusion.matrix(values,labels,t) })
  recalls = unlist(d[6,])
  precisions = unlist(d[5,])
  
  area = mean(precisions)
  
  
  y = rep(precisions,each=2)
  x = c(0,rep(recalls[1:(length(recalls)-1)],each=2),recalls[length(recalls)])

  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.ap.slow,
    conf.level=conf.level,
    method=conf.int.method)
  
  return (list(area=area,
               x=x,
               y=y,
               conf.int=ci))    
}

prcurve.ap <- function(pos.values, neg.values, conf.level=0.95, conf.int.method="binomial") {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]

  
  tp = 0
  fp = 0

  rs = rep(0,length(pos.values))
  ps = rep(0,length(pos.values))
  z = 1
  cur.tp = 0
  cur.fp = 0
  for (i in 1:length(indices)) {
    if (labels[i]==1) {
      cur.tp = cur.tp + 1      
    }
    else {
      cur.fp = cur.fp + 1
    }
    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can put a threshold here
      tp = tp + cur.tp
      fp = fp + cur.fp

      if (cur.tp > 0) {
        for (j in 1:cur.tp) {
          rs[z] = tp/length(pos.values)
          ps[z] = tp/(tp+fp)
          z = z + 1              
        }
      }
      cur.tp = 0
      cur.fp = 0
    }
    
  }

  
  if (z != length(pos.values)+1) {
    cat("WARNING: did not fill recall and precision arrays correctly (z: received=",z,", expected=",length(pos.values)+1,")\n",sep="")
  }
  area = mean(ps)

  y = rep(ps,each=2)
  x = c(0,rep((rs[1:(length(rs)-1)]),each=2),rs[length(rs)])

  
  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.ap,
    conf.level=conf.level,
    method=conf.int.method)
  
  return (list(area=area,
               x=x,
               y=y,
               conf.int=ci))  
}

## estimate PR curve and area under PR curve with MLE parameters for
## positive and negative normal distributions
prcurve.binormal <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  ## Assume X ~ N(0,1), i.e. negatives have standard normal distribution.
  ## So only have 2 df to calculate.

  if (length(pos.values)<=1) {
    cat("ERROR: cannot use binormal estimate with fewer than 2 positive samples\n")
  }
  if (length(neg.values)<=1) {
    cat("ERROR: cannot use binormal estimate with fewer than 2 negative samples\n")
  }
  mean.hat = (mean(pos.values)-mean(neg.values))/sd(neg.values)
  sd.hat = sd(pos.values)/sd(neg.values)

  pi = length(pos.values)/(length(pos.values)+length(neg.values))
  
  area = binormal.area(pi=pi,pos.mean=mean.hat,pos.sd=sd.hat)

  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.binormal,
    conf.level=conf.level,
    method=conf.int.method)
  
  return(list(area=area,pi.hat=pi,pos.mean.hat=mean.hat,pos.sd.hat=sd.hat, conf.int=ci))
}


## Use linear interpolation in ROC space, translate to PR space and
## find the precision for the specified recall (r) for interpolating
## between PR points (r1,p1) and (r2,p2)
interpolate.point <- function(r1,p1,r2,p2,r) {
  res = r / (r*(1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))) + (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1)))
  res[r==0] = p2
  return(res)
## if (r==0) {
##   ## hmm, use limit which is just p2
##   return(p2)
## }
## else {
##   return(r / (r*(1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))) + (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1))))
## }
}

interpolate.curve <- function(recalls,precisions,num.samples=1000) {

  ## sort by increasing recall
  indices = order(recalls,decreasing=FALSE)
  rs = recalls[indices]
  ps = precisions[indices]

  xs = rep(0,num.samples)
  ys = rep(0,num.samples)
  ## point lies between rs[index-1] and rs[index]
  index = 1
  for (i in 1:num.samples) {
    ## evenly spaced in [0,1], including end-points
    recall = (i-1)/(num.samples-1)
    xs[i] = recall
    
    while (rs[index]<recall) {
      index = index + 1
    }

    if (index==1) {
      ## use (0,1) as beginning interp point (creates horizontal line
      ## from smallest (r,p) at y=p
      ys[i] = interpolate.point(0,1,rs[index],ps[index],recall)
    }
    else {
      ys[i] = interpolate.point(rs[index-1],ps[index-1],rs[index],ps[index],recall)
    }
  }
  return(list(x=xs,y=ys))
}

## Use ROC interpolation to calculate areas between PR points
## Assumes recalls has no duplicates and sorted in ascending order
## only calculates area from min(recalls) to max(recalls)
interpolate.area <- function(recalls, precisions) {
  area = 0.0

  for (i in 2:length(recalls)) {
    ## calculate area for recalls[i-1] to recalls[i]

    r1 = recalls[i-1]
    p1 = precisions[i-1]
    r2 = recalls[i]
    p2 = precisions[i]

    if (r1 == r2) {
      cat("WARNING: recalls passed to interpolate.area are not unique\n")
    }
    if (r1 > r2) {
      cat("WARNING: recalls passed to interpolate.area are not sorted ascending\n")
      cat("recalls: ",recalls,"\n")
      cat("precisions: ",precisions,"\n")
    }
    
    if (r1==0) {
      ## definite integral is undefined
      ## use rectangle with height p2 (as interpolate.curve creates)
      area = area + p2*(r2-r1)
      
    }
    else {
      ## formula between these is p' = r' / (a + b*r')
      a = (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1))
      b = 1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))
      
      ## indefinite integral is (bx - a log (a + bx))/b^2
      ## area contribution is definite integral from r1 to r2
      area = area + (b*r2 - a*log(a+b*r2))/(b*b) - (b*r1 - a*log(a+b*r1))/(b*b)
    }
    
  }
  return(area)
}

## estimate PR curve and AUCPR by interpolating between the max
## precisions for each recall
prcurve.interpolate <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial",aggregator = max) {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions

  r.unique = unique(recalls)
  p.max = sapply(r.unique,function(r) { aggregator(precisions[recalls==r])})

  ## sort
  indices = order(r.unique,decreasing=FALSE)
  rs = r.unique[indices]
  ps = p.max[indices]

  area = interpolate.area(rs,ps)

  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator= function(p.v,n.v,c.l,c.i.m) { prcurve.interpolate(pos.values=p.v,neg.values=n.v,conf.level=c.l,conf.int.method=c.i.m,aggregator) },
    conf.level=conf.level,
    method=conf.int.method)

  curve = interpolate.curve(rs,ps)
  
  return(list(area=area,x=curve$x,y=curve$y,conf.int=ci))
  
}


roccurve.points <- function(pos.values,neg.values) {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]

  tp = 0
  fp = 0

  z = 1
  tprs = rep(0,length(indices)+1)
  fprs = rep(0,length(indices)+1)

  for (i in 1:length(indices)) {
    if (labels[i]==1) { tp = tp + 1 }
    else { fp = fp + 1 }

    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can place threshold

      ## if first, insert the all negative point (tpr=0,fpr=0)
      if (z==1) {
        tprs[z] = 0
        fprs[z] = 0
        z = z + 1
      }
      tprs[z] = tp/length(pos.values)
      fprs[z] = fp/length(neg.values)
      z = z + 1
    }      
  }

  ## truncate to used indices
  tprs = tprs[1:(z-1)]
  fprs = fprs[1:(z-1)]

  return (list(tprs=tprs,fprs=fprs))
}

## Estimate PR curve by taking convex hull in ROC space, convert only
## points on convex hull in ROC to PR, then use interpolation to go
## between.
prcurve.interpolate.convex <- function(pos.values, neg.values, conf.level=0.95,conf.int.method="binomial") {
  pi = length(pos.values)/(length(pos.values)+length(neg.values))
  
  t = roccurve.points(pos.values,neg.values)

  ## add the point (1,0) to make convex hull calculation clean (won't
  ## pick up any worse than random points)
  points = matrix(c(t$fprs,1,t$tprs,0),ncol=2)

  hull.indices = chull(points)

  ## remove the dummy point (0,1) we added, it's always the last element
  ## with index length(points)/2
  indices = hull.indices[!hull.indices==length(points)/2]
  points = points[indices,]

  ## due to bugs in chull and just to be safe sort points ascending by
  ## recall (second column)
  points = points[order(points[,2],decreasing=FALSE),]
    
  ## convert to PR space
  precisions = rep(0,length(indices))
  recalls = rep(0,length(indices))
  for (i in 1:length(indices)) {
    if (points[i,1]==0 & points[i,2]==0) {
      ## all negative
      precisions[i] = 1.0
      recalls[i] = 0.0
    }
    else {
      recalls[i] = points[i,2]
      precisions[i] = pi*points[i,2]/(pi*points[i,2] + (1-pi)*points[i,1])
    }
  }

  ## remove duplicate recalls, most likely having multiple recall=1.0
  ## points coming from the y=1 horizontal line in ROC space
  recalls.unique = unique(recalls)
  precisions.unique = sapply(recalls.unique,function(r) { max(precisions[recalls==r])})

  ## check that recalls are sorted ascending, having problems with non-sorted recalls
  if (is.unsorted(recalls.unique)) {
    cat("pos.values = ",pos.values,"\n")
    cat("neg.values = ",neg.values,"\n")
    cat("recalls = ",recalls,"\n")
    cat("precisions = ",precisions,"\n")
    cat("recalls.unique = ",recalls.unique,"\n")
    cat("precisions.unique = ",precisions.unique,"\n")
    
  }
  
  area = interpolate.area(recalls.unique,precisions.unique)
  
  ci = aucpr.conf.int(area,
    pos.values=pos.values,
    neg.values=neg.values,
    prcurve.estimator=prcurve.interpolate.convex,
    conf.level=conf.level,
    method=conf.int.method)

  curve = interpolate.curve(recalls.unique,precisions.unique)

  return(list(area=area,x=curve$x,y=curve$y,conf.int=ci))
  
}

## stratified bootstrap
## resample positive and negative separately
bootstrap <- function(pos.values,neg.values,prcurve.function,replicates=1000) {
  areas = rep(0,replicates)
  
  for (i in 1:replicates) {
    p.v = sample(pos.values,replace=TRUE)
    n.v = sample(neg.values,replace=TRUE)

    res = prcurve.function(p.v,n.v)
    areas[i] = res$area
  }

  return(areas)
}




