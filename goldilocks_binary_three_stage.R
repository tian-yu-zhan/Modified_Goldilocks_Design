
### three stage design
#scen.ind = 7

# setwd("/Users/zhantx/Documents/code/goldilocks/")
setwd("~/goldilocks")
library(rmutil)
library(doParallel)
n.cluster = 79

###### sample from beta-binomial
beta.bin.sample <- function(prior.a, prior.b, n, x, n.out){
  post.a = x + prior.a
  post.b = n - x + prior.b
  out.result = 
    rbeta(n.out, shape1 = post.a, shape2 = post.b, ncp = 0)
  return(out.result)
}

###### indicator of prop.test test
sign.prop.test = function(x.in, n.in){
  if (sum(x.in)==0) x.in = c(1,1)
  temp = prop.test(x = x.in, n = n.in, alternative = "less",correct = FALSE)$p.value
  # temp = pvalue.chi(c(x.in[1], n.in[1]-x.in[1]), c(x.in[2], n.in[2]-x.in[2]))
  return(as.numeric(temp<=alpha))
}

###### p-value of prop.test test
p.prop.test = function(x.in, n.in){
  if (sum(x.in)==0) x.in = c(1,1)
  temp = prop.test(x = x.in, n = n.in, alternative = "less", correct = FALSE)$p.value
  # temp = pvalue.chi(c(x.in[1], n.in[1]-x.in[1]), c(x.in[2], n.in[2]-x.in[2]))
  return(temp)
}

###### calculate p-value from fisher
p.cal.fisher = function(p1.in, p2.in){
  return(p1.in*p2.in - p1.in*p2.in*log(p1.in*p2.in))
}

##### indicator of significant result in chi-square test
sign.chi = function(n.pbo, n.trt, alpha){
  temp = pvalue.chi(n.pbo, n.trt)
  #temp = chisq.test(temp.in, correct = FALSE)$p.value
  dec = FALSE
  if (is.na(temp)) temp = 1
  
  ## one-sided test
  if (temp< alpha) dec = TRUE
  #if ((n.pbo[1])/sum(n.pbo)>(n.trt[1])/sum(n.trt)) dec = FALSE
  return(dec)
}

##### pvalue chi-square test
pvalue.chi = function(n.pbo, n.trt){
  temp.in = rbind(n.pbo, n.trt)
  chi.fit = chisq.test(temp.in, correct = TRUE)
  
  if ((n.pbo[1])/sum(n.pbo)>(n.trt[1])/sum(n.trt)){
    stat = - sqrt(chi.fit$statistic)
  } else {
    stat = sqrt(chi.fit$statistic)
  }
  
  p.value = 1-pnorm(stat)
  if (is.na(p.value)) p.value = 1
  return(p.value)
}

##### pvalue chi-square test
pvalue.fisher = function(n.pbo, n.trt){
  temp.in = rbind(n.pbo, n.trt)
  temp = fisher.test(temp.in, alternative = "two.sided")$p.value
  return(temp)
}

#### beta-binomial density
beta.bin.den = function(prior.a, prior.b, n, k){
  temp = choose(n, k)*beta(k+prior.a, n-k+prior.b)/beta(prior.a, prior.b)
  return(temp)
}

###### calculate P_n
pred.prob <- function(nc.1.un, nt.1.un, prior.a, prior.b, xc.1, xt.1, nc.1, nt.1, alpha){
  
  pn.mat = matrix(0, nrow = nc.1.un+1, ncol = nt.1.un+1)
  for (i in 0:nc.1.un){
    for (j in 0:nt.1.un){
      ind = sign.chi(n.pbo = c(xc.1+i, nc.1+nc.1.un-xc.1-i),
                     n.trt = c(xt.1+j, nt.1+nt.1.un-xt.1-j), alpha)
      #ind = TRUE
      if (!ind) next
      pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
        beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)

    }
  }

  
  # time1 = Sys.time()
  # temp.fun = function(x, y){
  #   ind = sign.chi(n.pbo = c(xc.1+x, nc.1+nc.1.un-xc.1-x),
  #                                 n.trt = c(xt.1+y, nt.1+nt.1.un-xt.1-y), alpha)
  # 
  #   if (ind){
  #     temp.out = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1,
  #                  n = nc.1.un, k=x)*
  #       beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=y)
  #   } else{
  #     temp.out = 0
  #   }
  #   return(temp.out)
  # }
  # 
  # temp = outer(0:nc.1.un, 0:nt.1.un, Vectorize(temp.fun))
  # print(Sys.time()- time1)
  #  print(sum(pn.mat))
  #  print(sum(temp))
  
  return(sum(pn.mat))
  
  #print(sum(pn.mat)==sum(pn.mat.1))
  
  # pn.mat = matrix(0, nrow = nc.1.un+1, ncol = nt.1.un+1)
  # for (i in 0:nc.1.un){
  #   ind = FALSE
  #   for (j in 0:nt.1.un){
  #     if (!ind)
  #       ind = sign.chi(n.pbo = c(xc.1+i, nc.1+nc.1.un-xc.1-i),
  #                      n.trt = c(xt.1+j, nt.1+nt.1.un-xt.1-j), alpha)
  # 
  #     if (!ind) next
  #     pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
  #       beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)
  # 
  #   }
  # }
  
}

#### solve nonlinear equation of critical value c
solve.c = function(fun, end){
  cand.vec = seq(0.000001, end-0.0000001, length.out = 10^7)
  out.vec = fun(cand.vec)
  out = cand.vec[which.min(abs(out.vec))]
  return(out)
}

#### solve c for inverse normal combination function: tol = 10E-6 
solve.c.inverse = function(w1, w2, alpha, a0, tol){
  inverse.num = function(c, w1, w2, a0, precision){
    n = precision
    p1 = seq(0, a0, length.out = n)
    out = sum(1-pnorm((qnorm(1-c)-w1*qnorm(1-p1))/w2))*a0/n
    return(out)
  }
  
  # x.cand = seq(0.00000001, 0.9999999, length.out = 10^4)
  # y.cand = rep(0, 10^4)
  # 
  # for (i in 1:10^4){
  #   print(i)
  #   x = x.cand[i]
  #   y.cand[i] = inverse.num(c= x, w1=w1, w2=w2, a0=a0, precision = 10^5)
  # }
  
  
  start = 0.000000001
  end = a0 - 0.000000001
  ind = TRUE
  
  while(ind){
    mid = (start+end)/2
    out.start = inverse.num(start, w1, w2, a0, 1/tol) - alpha
    out.end = inverse.num(end, w1, w2, a0, 1/tol) - alpha
    out.mid = inverse.num(mid, w1, w2, a0, 1/tol) - alpha
    if (out.mid>0){
      end = mid
    } else {
      start = mid
    }
    if (abs(out.start - out.end)<tol) ind = FALSE
  }
  
  return(mid)
}

#### function to calculate futility p value
solve.a0 = function(input){
  alpha = input$alpha
  prior.a = input$prior.a
  prior.b = input$prior.b
  
  pn.cut = input$pn.cut
  p.max.cut = input$p.max.cut
  
  nc.1 = input$nc.1 ## number of pbo in check 1
  nc.1.un = input$nc.1.un ## unobserved observations
  nc.max = input$nc.max
  nc.2 = nc.max - nc.1 - nc.1.un
  
  nt.1 = input$nt.1 ## number of trt in check 1
  nt.1.un = input$nt.1.un
  nt.max = input$nt.max
  nt.2 = nt.max - nt.1 - nt.1.un
  
  rate.c = input$rate.c ## rate in control
  rate.t = input$rate.t ## rate in trt
  
  a0.pmax.mat = a0.out.mat = a0.mat = matrix(NA, nrow = nc.1+1, ncol = nt.1+1)
  for (i in 0:nt.1){
    print(paste("solve a0", i))
    j.left = i
    ind.left = TRUE
    while ((ind.left)&(j.left>=0)){
      p.max.temp = pred.prob(nc.1.un = nc.1.un + nc.2,
                        nt.1.un = nt.1.un + nt.2,
                        prior.a,
                        prior.b,
                        i,
                        j.left,
                        nc.1,
                        nt.1,
                        alpha)
      if (p.max.temp < p.max.cut){
        a0.mat[(i+1), (j.left+1)] = pvalue.chi(c(i, nc.1-i), c(j.left, nt.1-j.left))
        a0.pmax.mat[(i+1), (j.left+1)] = p.max.temp
        j.left = j.left-1
      } else {
        ind.left=FALSE
      }
      
    }
    if ((j.left)>=0) a0.out.mat[(i+1), (j.left+1)] = pvalue.chi(c(i, nc.1-i), c(j.left, nt.1-j.left))
      
  
    
    j.right = i+1
    ind.right = TRUE
    while ((ind.right)&(j.right<=nt.1)){
      p.max.temp = pred.prob(nc.1.un = nc.1.un + nc.2,
                             nt.1.un = nt.1.un + nt.2,
                             prior.a,
                             prior.b,
                             i,
                             j.right,
                             nc.1,
                             nt.1,
                             alpha)
      if (p.max.temp < p.max.cut){
        a0.mat[(i+1), (j.right+1)] = pvalue.chi(c(i, nc.1-i), c(j.right, nt.1-j.right))
        a0.pmax.mat[(i+1), (j.right+1)] = p.max.temp
        j.right = j.right+1
      } else {
        ind.right=FALSE
      }
      
    }
    if ((j.right)<=nt.1) a0.out.mat[(i+1), (j.right+1)] = pvalue.chi(c(i, nc.1-i), 
                                                                       c(j.right, nt.1-j.right))
    
  }
  
  newlist = 
    list(a0 = max(a0.out.mat, na.rm=TRUE), a0.out.mat = a0.out.mat, a0.mat = a0.mat, a0.pmax.mat = a0.pmax.mat)
  
  return(newlist)

}

#### function to get pn table
solve.pn.table = function(nc.1, nt.1, nc.1.un, nt.1.un, nc.max, nt.max, prior.a, prior.b,
                          alpha){
  
  pn.mat = pmax.mat = matrix(NA, nrow = nc.1+1, ncol = nt.1+1)
  
  for (i in 0:nc.1){
    for (j in 0:nt.1){
      print(paste("solve pn table i", i, "j", j))
      pn = pred.prob(nc.1.un = nc.1.un,
                     nt.1.un = nt.1.un,
                     prior.a,
                     prior.b,
                     i,
                     j,
                     nc.1,
                     nt.1,
                     alpha)
      
      pn.mat[i+1, j+1] = pn
      
      p.max = pred.prob(nc.1.un = nc.max - nc.1,
                        nt.1.un = nt.max - nt.1,
                        prior.a,
                        prior.b,
                        i,
                        j,
                        nc.1,
                        nt.1,
                        alpha)
      
      pmax.mat[i+1, j+1] = p.max
      
    }
  }
  
  write.csv(pn.mat, 
            paste0("results/cutoff/same/pn_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                   "_nc.max_", nc.max, ".csv"), row.names = FALSE)
  write.csv(pmax.mat, 
            paste0("results/cutoff/same/pmax_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                   "_nc.max_", nc.max, ".csv"), row.names = FALSE)
}

#### function to calculate power
power.cal = function(input){
  alpha = input$alpha
  prior.a = input$prior.a
  prior.b = input$prior.b
  
  pn.cut = input$pn.cut
  p.max.cut = input$p.max.cut
  
  nc.1 = input$nc.1 ## number of pbo in check 1
  nc.1.un = input$nc.1.un ## unobserved observations
  nc.2 = input$nc.2
  nc.max = input$nc.max

  nt.1 = input$nt.1 ## number of trt in check 1
  nt.1.un = input$nt.1.un
  nt.2 = input$nt.2
  nt.max = input$nt.max
  
  rate.c = input$rate.c ## rate in control
  rate.t = input$rate.t ## rate in trt
  
  #######################################################################
  ## solve combination function
  c.fisher = solve.c(function(x){1/2*x*(log(x))^2+x-x*log(x)-alpha}, end = alpha)
  
  c.com.vec = c(0.015, 0.016, 0.017, 0.018, 0.019, 0.02)
  n.sim.inner = 1000

  ############################################ simulations
  n.itt = input$n.itt

  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  para.fit = foreach(itt=1:n.itt) %dopar% {
    
    ###### indicator of prop.test test
    sign.prop.test = function(x.in, n.in){
      if (sum(x.in)==0) x.in = c(1,1)
      temp = prop.test(x = x.in, n = n.in, alternative = "less",correct = FALSE)$p.value
      # temp = pvalue.chi(c(x.in[1], n.in[1]-x.in[1]), c(x.in[2], n.in[2]-x.in[2]))
      return(as.numeric(temp<=alpha))
    }
    
    ###### p-value of prop.test test
    p.prop.test = function(x.in, n.in){
      if (sum(x.in)==0) x.in = c(1,1)
      temp = prop.test(x = x.in, n = n.in, alternative = "less", correct = FALSE)$p.value
      # temp = pvalue.chi(c(x.in[1], n.in[1]-x.in[1]), c(x.in[2], n.in[2]-x.in[2]))
      return(temp)
    }
    
    ###### calculate p-value from fisher
    p.cal.fisher = function(p1.in, p2.in){
      return(p1.in*p2.in - p1.in*p2.in*log(p1.in*p2.in))
    }
    
    # ##### indicator of significant result in chi-square test
    # sign.chi = function(n.pbo, n.trt, alpha){
    #   temp = pvalue.chi(n.pbo, n.trt)
    #   #temp = chisq.test(temp.in, correct = FALSE)$p.value
    #   dec = FALSE
    #   if (is.na(temp)) temp = 1
    #   
    #   ## one-sided test
    #   if (temp< alpha) dec = TRUE
    #   #if ((n.pbo[1])/sum(n.pbo)>(n.trt[1])/sum(n.trt)) dec = FALSE
    #   return(dec)
    # }
    # 
    # ##### pvalue chi-square test
    # pvalue.chi = function(n.pbo, n.trt){
    #   temp.in = rbind(n.pbo, n.trt)
    #   chi.fit = chisq.test(temp.in, correct = TRUE)
    #   
    #   if ((n.pbo[1])/sum(n.pbo)>(n.trt[1])/sum(n.trt)){
    #     stat = - sqrt(chi.fit$statistic)
    #   } else {
    #     stat = sqrt(chi.fit$statistic)
    #   }
    #   
    #   p.value = 1-pnorm(stat)
    #   if (is.na(p.value)) p.value = 1
    #   return(p.value)
    # }
    
    
    library(rmutil)
    p1 = p2 = p3 = NULL
    
    # for (itt in 1:n.itt){
    # if ((itt %% 10)==0) print(paste("scen", scen.ind, 'itt', itt))
    ## observed events in stage 1. 
    xc.1 = rbinom(1, nc.1, rate.c)
    xt.1 = rbinom(1, nt.1, rate.t)
    
    p1 = p.prop.test(c(xc.1, xt.1), c(nc.1, nt.1))
    
    ## first stage pn and pmax
    pn.pmax.func = function(pbo.vec.in, trt.vec.in, 
                            n.unobs.in,
                            prior.a.in, prior.b.in,
                            n.sim.inner){
      pn.vec = sapply(1:n.sim.inner, function(x){
        
        pbo.new.sample =rbetabinom(n = 1, size = n.unobs.in, 
                                   m = ((prior.a.in+pbo.vec.in[1]))/(prior.a.in+prior.b.in+pbo.vec.in[2]), 
                                   s = (prior.a.in+prior.b.in+pbo.vec.in[2]))
        trt.new.sample = rbetabinom(n = 1, size = n.unobs.in, 
                                    m = ((prior.a.in+trt.vec.in[1]))/(prior.a.in+prior.b.in+trt.vec.in[2]), 
                                    s = (prior.a.in+prior.b.in+trt.vec.in[2]))
        
        
        x.new.in = c(pbo.vec.in[1], trt.vec.in[1])+c(pbo.new.sample, trt.new.sample)
        n.new.in = c(pbo.vec.in[2], trt.vec.in[2])+n.unobs.in
        return(sign.prop.test(x.new.in, n.new.in))
      })
      return(mean(pn.vec))
    }
    
    pn.1 = pn.pmax.func(pbo.vec.in = c(xc.1, nc.1),
                        trt.vec.in = c(xt.1, nt.1),
                        n.unobs.in = nc.1.un,
                        prior.a.in = prior.a,
                        prior.b.in = prior.b,
                        n.sim.inner = n.sim.inner)
    
    pmax.1 = pn.pmax.func(pbo.vec.in = c(xc.1, nc.1),
                          trt.vec.in = c(xt.1, nt.1),
                          n.unobs.in = nc.max - nc.1,
                          prior.a.in = prior.a,
                          prior.b.in = prior.b,
                          n.sim.inner = n.sim.inner)
    
    ind.futility = 1
    status.out = "stage_1_stop_futility"
    p.naive.out = 1
    if (pmax.1>=p.max.cut){
      # stage 1 does not stop for futility
      
      ind.futility = 0
      xc.1.un = rbinom(1, nc.1.un, rate.c)
      xt.1.un = rbinom(1, nt.1.un, rate.t)
      
      if (pn.1>=pn.cut){
        # stop recruiting
        
        status.out = "stage_1_stop_recruiting"
        p2 = p.prop.test(c(xc.1.un, xt.1.un), 
                                     c(nc.1.un, nt.1.un))
        
        p.naive.out = p.prop.test(c(xc.1+xc.1.un, xt.1+xt.1.un), 
                                  c(nc.1+nc.1.un, nt.1+nt.1.un))
        
      } else {
        # continue recruiting
        xc.2.left = rbinom(1, nc.2 - nc.1 - nc.1.un, rate.c)
        xt.2.left = rbinom(1, nt.2 - nt.1 - nt.1.un, rate.t)
        
        p2 = p.prop.test(c(xc.1.un+xc.2.left, xt.1.un+xt.2.left), 
                                       c(nc.2 - nc.1, nt.2 - nt.1))
        
        pn.2 = pn.pmax.func(pbo.vec.in = c(xc.1 + xc.1.un + xc.2.left, nc.2),
                            trt.vec.in = c(xt.1 + xt.1.un + xt.2.left, nt.2),
                            n.unobs.in = nc.1.un,
                            prior.a.in = prior.a,
                            prior.b.in = prior.b,
                            n.sim.inner = n.sim.inner)
        
        pmax.2 = pn.pmax.func(pbo.vec.in = c(xc.1 + xc.1.un + xc.2.left, nc.2),
                              trt.vec.in = c(xt.1 + xt.1.un + xt.2.left, nt.2),
                              n.unobs.in = nc.max - nc.2,
                              prior.a.in = prior.a,
                              prior.b.in = prior.b,
                              n.sim.inner = n.sim.inner)
        
        ## second stage check
        if (pmax.2<p.max.cut){
          status.out = "stage_2_stop_futility"
          ind.futility = 1
          p.naive.out = 1
          # p3=1
        } else {
          xc.2.un = rbinom(1, nc.1.un, rate.c)
          xt.2.un = rbinom(1, nt.1.un, rate.t)
          
          if (pn.2>=pn.cut){
            # stop recruiting
            status.out = "stage_2_stop_recruiting"
            p3 = p.prop.test(c(xc.2.un, xt.2.un), 
                                         c(nc.1.un, nt.1.un))
            p.naive.out = p.prop.test(c(xc.1+xc.1.un+xc.2.left+xc.2.un, 
                                        xt.1+xt.1.un+xt.2.left+xt.2.un), 
                                      c(nc.2+nc.1.un, 
                                        nt.2+nt.1.un))
            
          } else {
            
            status.out = "stage_2_max"
            xc.3.left = rbinom(1, nc.max - nc.2 - nc.1.un, rate.c)
            xt.3.left = rbinom(1, nt.max - nt.2 - nt.1.un, rate.t)
            
            p3 = p.prop.test(c(xc.2.un+xc.3.left, xt.2.un+xt.3.left), 
             
                                           c(nc.max - nc.2, nt.max - nt.2))
            
            p.naive.out = p.prop.test(c(xc.1+xc.1.un+xc.2.left+xc.2.un+xc.3.left, 
                                        xt.1+xt.1.un+xt.2.left+xt.2.un+xt.3.left), 
                                      c(nc.max, 
                                        nt.max))
            
          }
        }
      }
    }
    
    # print(p.vec)
    ind.futility.out = ind.futility
    # p.vec.list[[itt]] = p.vec
    
    return(list("p.vec" = c(p1, p2, p3),
                "ind.futility" = ind.futility,
                "p.naive" = p.naive.out))
    ## finish trial  
  }
  
  stopCluster(cl)
  

  sign.mat = matrix(NA, nrow= n.itt, ncol = 2+length(c.com.vec))
  colnames(sign.mat) = c("fisher","inverse",c.com.vec)
  
  for (itt in 1:n.itt){
    para.fit.temp = para.fit[[itt]]
    
    p.vec.temp = para.fit.temp$p.vec
    ind.futility.temp = para.fit.temp$ind.futility
    p.naive.temp = para.fit.temp$p.naive
    
    if (length(p.vec.temp)==1){
      p.fish = p.inv = as.numeric(p.vec.temp)
    } else if (length(p.vec.temp)==2){
      p.fish = p.cal.fisher(p.vec.temp[1], p.vec.temp[2])
      p.inv = pnorm(sum(qnorm(p.vec.temp,lower.tail = FALSE)*sqrt(c(1/2, 1/2))),lower.tail = FALSE)
    } else {
      p.fish = p.cal.fisher(p.vec.temp[1],
                            p.cal.fisher(p.vec.temp[2], p.vec.temp[3]))
      p.inv = pnorm(sum(qnorm(p.vec.temp,lower.tail = FALSE)*sqrt(c(1/2, 1/4, 1/4))),lower.tail = FALSE)
    }
    
    sign.mat[itt, 1] = (p.fish<=alpha)&(ind.futility.temp==0)
    sign.mat[itt, 2] = (p.inv<=alpha)&(ind.futility.temp==0)
    # sign.mat[itt, 2] = p.vec.temp[1]
    for (sign.mat.in in 1:length(c.com.vec)){
      sign.mat[itt, 2+sign.mat.in] = (p.naive.temp<=c.com.vec[sign.mat.in])&(ind.futility.temp==0)
    }
    

  }
  print(apply(sign.mat, 2, function(x){mean(x,na.rm = TRUE)}))
  
  return(apply(sign.mat, 2, function(x){mean(x,na.rm = TRUE)}))
}

#######################################################################################################

## type 1 error simulations
time.temp = Sys.time()
output.table = NULL

for (scen.ind in c(1:18)){

  #random.seed = NA
  random.seed = 1
  set.seed(random.seed)

  #print(scen.ind)
  n.itt.in = 10^5
  #n.itt.in = 1000

  p.max.cut.in = 0
  
  if (scen.ind<=6){
    pn.cut.in = 0.8

      if (scen.ind==1) {rate.c.in = 0.2; delta.in = 0}
      if (scen.ind==2) {rate.c.in = 0.4; delta.in = 0}
      if (scen.ind==3) {rate.c.in = 0.6; delta.in = 0}
      if (scen.ind==4) {rate.c.in = 0.8; delta.in = 0}

      if (scen.ind==5) {rate.c.in = 0.6; delta.in = 0.1}
      if (scen.ind==6) {rate.c.in = 0.6; delta.in = 0.15}

    
  } else if (scen.ind<=12){
    pn.cut.in = 0.85

      if (scen.ind==7) {rate.c.in = 0.2; delta.in = 0}
      if (scen.ind==8) {rate.c.in = 0.4; delta.in = 0}
      if (scen.ind==9) {rate.c.in = 0.6; delta.in = 0}
      if (scen.ind==10) {rate.c.in = 0.8; delta.in = 0}

      if (scen.ind==11) {rate.c.in = 0.6; delta.in = 0.1}
      if (scen.ind==12) {rate.c.in = 0.6; delta.in = 0.15}
    
  } else {
    pn.cut.in = 0.9

      if (scen.ind==13) {rate.c.in = 0.2; delta.in = 0}
      if (scen.ind==14) {rate.c.in = 0.4; delta.in = 0}
      if (scen.ind==15) {rate.c.in = 0.6; delta.in = 0}
      if (scen.ind==16) {rate.c.in = 0.8; delta.in = 0}

      if (scen.ind==17) {rate.c.in = 0.6; delta.in = 0.1}
      if (scen.ind==18) {rate.c.in = 0.6; delta.in = 0.15}
    
  }
  
  

  
  rate.t.in = rate.c.in + delta.in
  
  nc.1.in = 100 # first stage observed
  nc.1.un.in = 40 # first stage unobserved
  nc.2.in = 200 # second stage observed (> (nc.1.in + nc.1.un.in))
  
  input = list(alpha = 0.025, # chi square test is for two-sided test
               prior.a = 1,
               prior.b = 1,
               pn.cut = pn.cut.in,
               p.max.cut = p.max.cut.in,
               nc.1 = nc.1.in,
               nc.1.un = nc.1.un.in,
               nc.2 = nc.2.in, 
               nc.max = 300,
               nt.1 = nc.1.in,
               nt.1.un = nc.1.un.in,
               nt.2 = nc.2.in, 
               nt.max = 300,
               rate.c = rate.c.in,
               rate.t = rate.t.in,
               n.itt = n.itt.in)

  power.fit = power.cal(input)

  output.table = rbind(output.table,
                       power.fit
                       )
  rownames(output.table)[scen.ind] =
    paste("rate_c", rate.c.in, "rate_t", rate.t.in, "Sn", pn.cut.in, "Fn", p.max.cut.in,
                                           "n1", nc.1.in, "n1.un", nc.1.un.in)
  
  write.csv(output.table, paste0("results/subset/three_test", ".csv"))
}

# colnames(output.table) = c("fisher", "inverse", "goldilocks", "fisher unif", "inverse unif", "stop for futility",
#                            "stop recruiting", "continue recruiting", "a0","min p1 coh 1",
#                            "max p1 coh 2", "random seed", "n.itt")

print(Sys.time() - time.temp)
print(output.table)

