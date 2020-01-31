
#add ppn^prime
#setwd("/Users/zhantx/Documents/code/goldilocks/")
# supp table 1 and 2
setwd("~/goldilocks")

###### sample from beta-binomial
beta.bin.sample <- function(prior.a, prior.b, n, x, n.out){
  post.a = x + prior.a
  post.b = n - x + prior.b
  out.result = 
    rbeta(n.out, shape1 = post.a, shape2 = post.b, ncp = 0)
  return(out.result)
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
  chi.fit = chisq.test(temp.in, correct = FALSE)
  
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
pred.prob <- function(nc.1.un, nt.1.un, prior.a, prior.b, xc.1, xt.1, nc.1, nt.1, alpha, c){
  
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
  
  return(sum(pn.mat))
  
}

###### calculate P_n^\prime for fisher
pred.prob.fisher <- function(nc.1.un, nt.1.un, prior.a, prior.b, xc.1, xt.1, nc.1, nt.1, alpha,
                             c.fisher){
  
  pn.mat = matrix(0, nrow = nc.1.un+1, ncol = nt.1.un+1)
  for (i in 0:nc.1.un){
    for (j in 0:nt.1.un){
      p1 = pvalue.chi(n.pbo = c(xc.1, nc.1-xc.1), n.trt = c(xt.1, nt.1 - xt.1))
      p2 = pvalue.chi(n.pbo = c(i, nc.1.un-i), n.trt = c(j, nt.1.un - j))
      
      ind = FALSE
      if (p1*p2<c.fisher) ind = TRUE
      
      #ind = TRUE
      if (!ind) next
      pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
        beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)
      
    }
  }
  
  return(sum(pn.mat))
  
}

###### calculate P_n^\prime for inverse
pred.prob.inverse <- function(nc.1.un, nt.1.un, prior.a, prior.b, xc.1, xt.1, nc.1, nt.1, alpha, c.inverse){
  
  w1.inv = w2.inv = sqrt(0.5)
  
  pn.mat = matrix(0, nrow = nc.1.un+1, ncol = nt.1.un+1)
  for (i in 0:nc.1.un){
    for (j in 0:nt.1.un){
      p1 = pvalue.chi(n.pbo = c(xc.1, nc.1-xc.1), n.trt = c(xt.1, nt.1 - xt.1))
      p2 = pvalue.chi(n.pbo = c(i, nc.1.un-i), n.trt = c(j, nt.1.un - j))
      
      p1 = max(0.0000001, min(p1, 0.9999999))
      p2 = max(0.0000001, min(p2, 0.9999999))
      
      #print(c.inverse)
      ind = FALSE
      #print(p1)
      #print(p2)
      if ((1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2))) <c.inverse) ind = TRUE
      
      #ind = TRUE
      if (!ind) next
      pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
        beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)
      
    }
  }
  
  return(sum(pn.mat))
  
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
solve.a0 = function(input, pred.prob.fun, c.value){
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
      p.max.temp = pred.prob.fun(nc.1.un = nc.1.un + nc.2,
                        nt.1.un = nt.1.un + nt.2,
                        prior.a,
                        prior.b,
                        i,
                        j.left,
                        nc.1,
                        nt.1,
                        alpha, c.value)
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
      p.max.temp = pred.prob.fun(nc.1.un = nc.1.un + nc.2,
                             nt.1.un = nt.1.un + nt.2,
                             prior.a,
                             prior.b,
                             i,
                             j.right,
                             nc.1,
                             nt.1,
                             alpha, c.value)
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
                          alpha, pred.prob.fun, pred.name, c.value){
  
  pn.mat = pmax.mat = matrix(NA, nrow = nc.1+1, ncol = nt.1+1)
  
  for (i in 0:nc.1){
    for (j in 0:nt.1){
      print(paste("solve pn table i", i, "j", j))
      pn = pred.prob.fun(nc.1.un = nc.1.un,
                     nt.1.un = nt.1.un,
                     prior.a,
                     prior.b,
                     i,
                     j,
                     nc.1,
                     nt.1,
                     alpha, c.value)
      
      pn.mat[i+1, j+1] = pn
      
      p.max = pred.prob.fun(nc.1.un = nc.max - nc.1,
                        nt.1.un = nt.max - nt.1,
                        prior.a,
                        prior.b,
                        i,
                        j,
                        nc.1,
                        nt.1,
                        alpha, c.value)
      
      pmax.mat[i+1, j+1] = p.max
      
    }
  }
  
  write.csv(pn.mat, 
            paste0("results/cutoff/diff/pn_", pred.name, "_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                   "_nc.max_", nc.max, ".csv"), row.names = FALSE)
  write.csv(pmax.mat, 
            paste0("results/cutoff/diff/pmax_", pred.name, "_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
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
  nc.max = input$nc.max
  nc.2 = nc.max - nc.1 - nc.1.un
  
  nt.1 = input$nt.1 ## number of trt in check 1
  nt.1.un = input$nt.1.un
  nt.max = input$nt.max
  nt.2 = nt.max - nt.1 - nt.1.un
  
  rate.c = input$rate.c ## rate in control
  rate.t = input$rate.t ## rate in trt
  
  c.fisher.a0 = solve.c(function(x) x+x*log(1)-x*log(x)-alpha, end = 1)
  w1.inv = w2.inv = sqrt(0.5)
  
  if (!file.exists(paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                          "_w2_", w2.inv, "_a0_", 1, ".csv"))){
    c.inverse.a0 = solve.c.inverse(w1.inv, w2.inv, alpha, a0 = 1, 10^(-7))
    write.csv(c.inverse.a0, paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                                "_w2_", w2.inv, "_a0_", 1, ".csv"), 
              row.names = FALSE)
  }
  
  c.inverse.a0 = as.numeric(read.csv(paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                              "_w2_", w2.inv, "_a0_", 1, ".csv")))
  print(c.inverse.a0)
  
  
  
  if (p.max.cut==0){
    a0.gold = a0.fisher = a0.inverse = 1
  } else {
    if (!file.exists(paste0("results/cutoff/diff/a0_gold_nc.1_", nc.1, 
                            "_nc.max_", nc.max, "_pmax_cut_", p.max.cut, ".csv"))){
      
      a0.gold = solve.a0(input, pred.prob.fun = pred.prob, 1)$a0
      a0.fisher = solve.a0(input, pred.prob.fun = pred.prob.fisher, c.fisher.a0)$a0
      a0.inverse = solve.a0(input, pred.prob.fun = pred.prob.inverse, c.inverse.a0)$a0
      write.csv(a0.gold, paste0("results/cutoff/diff/a0_gold_nc.1_", nc.1, 
                           "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv"), 
                row.names = FALSE)
      write.csv(a0.fisher, paste0("results/cutoff/diff/a0_fisher_nc.1_", nc.1, 
                                "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv"), 
                row.names = FALSE)
      write.csv(a0.inverse, paste0("results/cutoff/diff/a0_inverse_nc.1_", nc.1, 
                                "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv"), 
                row.names = FALSE)
    }
    
    a0.gold = as.numeric(read.csv(paste0("results/cutoff/diff/a0_gold_nc.1_", nc.1, 
                                 "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv")))
    a0.fisher = as.numeric(read.csv(paste0("results/cutoff/diff/a0_fisher_nc.1_", nc.1, 
                                         "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv")))
    a0.inverse = as.numeric(read.csv(paste0("results/cutoff/diff/a0_inverse_nc.1_", nc.1, 
                                         "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv")))

    #a0=0.6273726
  }
  
  print(a0.gold)
  print(a0.fisher)
  print(a0.inverse)
  
  w1.inv = sqrt(0.5)
  w2.inv = sqrt(0.5)
  c.fisher = solve.c(function(x) x+x*log(a0.fisher)-x*log(x)-alpha, end = a0.fisher)
  print(c.fisher)
  c.com = 0.022
  
  if (!file.exists(paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                          "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"))){
    c.inverse = solve.c.inverse(w1.inv, w2.inv, alpha, a0 = a0.inverse, 10^(-7))
    write.csv(c.inverse, paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                                "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"), 
              row.names = FALSE)
  }
  
  c.inverse = read.csv(paste0("results/cutoff/diff/cinverse_w1_", w1.inv, 
                       "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"))
  c.inverse = as.numeric(c.inverse)
  print(c.inverse)
  
  ############################################ simulations
  n.itt = input$n.itt
  p1.fisher.vec = p2.fisher.vec = p1.inverse.vec = p2.inverse.vec = p.gold.vec = rep(NA, n.itt)
  dec.fisher.mat = dec.inverse.mat = dec.gold.mat = rep(NA, n.itt)
  
  sign.mat = matrix(0, nrow= n.itt, ncol = 3)
  colnames(sign.mat) = c("fisher","inverse","goldilocks")
  
  ## get pn pmax table
  if (!file.exists(paste0("results/cutoff/diff/pn_gold_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                         "_nc.max_", nc.max, ".csv"))){
    solve.pn.table(input$nc.1, input$nt.1, input$nc.1.un, input$nt.1.un,
                   input$nc.max, input$nt.max, input$prior.a, input$prior.b, input$alpha,
                   pred.prob.fun = pred.prob, pred.name = "gold", 1)
    solve.pn.table(input$nc.1, input$nt.1, input$nc.1.un, input$nt.1.un,
                   input$nc.max, input$nt.max, input$prior.a, input$prior.b, input$alpha,
                   pred.prob.fun = pred.prob.fisher, pred.name = "fisher", c.fisher.a0)
    solve.pn.table(input$nc.1, input$nt.1, input$nc.1.un, input$nt.1.un,
                   input$nc.max, input$nt.max, input$prior.a, input$prior.b, input$alpha,
                   pred.prob.fun = pred.prob.inverse, pred.name = "inverse", c.inverse.a0)
  }
  
  pn.c.gold.table = read.csv(paste0("results/cutoff/diff/pn_gold_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                               "_nc.max_", nc.max, ".csv"))
  pmax.c.gold.table = read.csv(paste0("results/cutoff/diff/pmax_gold_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                 "_nc.max_", nc.max, ".csv"))
  
  pn.c.fisher.table = read.csv(paste0("results/cutoff/diff/pn_fisher_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                    "_nc.max_", nc.max, ".csv"))
  pmax.c.fisher.table = read.csv(paste0("results/cutoff/diff/pmax_fisher_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                      "_nc.max_", nc.max, ".csv"))
  
  pn.c.inverse.table = read.csv(paste0("results/cutoff/diff/pn_inverse_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                      "_nc.max_", nc.max, ".csv"))
  pmax.c.inverse.table = read.csv(paste0("results/cutoff/diff/pmax_inverse_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                        "_nc.max_", nc.max, ".csv"))
  
  for (itt in 1:n.itt){
    if ((itt %% 1000)==0) print(paste("scen", scen.ind, 'itt', itt))
    ## observed events in stage 1. 
    xc.1 = rbinom(1, nc.1, rate.c)
    xt.1 = rbinom(1, nt.1, rate.t)
    
    pn.gold = max(min(1, pn.c.gold.table[xc.1+1, xt.1+1]), 0)
    pmax.gold = max(min(1, pmax.c.gold.table[xc.1+1, xt.1+1]), 0)
    
    pn.fisher = max(min(1, pn.c.fisher.table[xc.1+1, xt.1+1]), 0)
    pmax.fisher = max(min(1, pmax.c.fisher.table[xc.1+1, xt.1+1]), 0)
    
    pn.inverse = max(min(1, pn.c.inverse.table[xc.1+1, xt.1+1]), 0)
    pmax.inverse = max(min(1, pmax.c.inverse.table[xc.1+1, xt.1+1]), 0)
    
    p1 = pvalue.chi(c(xc.1, nc.1-xc.1), c(xt.1, nt.1-xt.1))
    p1.fisher.vec[itt] = p1.inverse.vec[itt] = p1
    
    ### naive adaptation
    # pn = 0.1
    # ptemp = pvalue.chi(c(xc.1, nc.1-xc.1), c(xt.1, nt.1-xt.1))
    # if (ptemp<0.05) pn = 0.9
    # p.max = 0.5

    ## Make interim look decisions for gold
    if (pmax.gold<p.max.cut){
      # stop for futility
      dec.gold.mat[itt] = 1
    } else {
      if (pn.gold>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un, rate.c)
        xt.2 = rbinom(1, nt.1.un, rate.t)
        
        dec.gold.mat[itt] = 2
        
        p.com = pvalue.chi(c(xc.1 + xc.2, nc.1+nc.1.un -xc.1-xc.2), 
                           c(xt.1 + xt.2, nt.1+nt.1.un -xt.1-xt.2))
        p.gold.vec[itt] = p.com
        
      } else {
        # continue recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un + nc.2, rate.c)
        xt.2 = rbinom(1, nt.1.un + nt.2, rate.t)
        
        dec.gold.mat[itt] = 3
        
        p.com = pvalue.chi(c(xc.1 + xc.2, nc.1+nc.1.un+nc.2 -xc.1-xc.2), 
                           c(xt.1 + xt.2, nt.1+nt.1.un+nt.2 -xt.1-xt.2))
        p.gold.vec[itt] = p.com
      }
      
      if (p.com < c.com) sign.mat[itt, 3] = 1
    }
    
    ## Make interim look decisions for fisher
    if (pmax.fisher<p.max.cut){
      # stop for futility
      dec.fisher.mat[itt] = 1
    } else {
      if (pn.fisher>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un, rate.c)
        xt.2 = rbinom(1, nt.1.un, rate.t)

        p2 = pvalue.chi(c(xc.2, nc.1.un - xc.2), c(xt.2, nt.1.un-xt.2))
        p2.fisher.vec[itt] = p2
        dec.fisher.mat[itt] = 2

      } else {
        # continue recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un + nc.2, rate.c)
        xt.2 = rbinom(1, nt.1.un + nt.2, rate.t)

        p2 = pvalue.chi(c(xc.2, nc.1.un + nc.2-xc.2), c(xt.2, nt.1.un + nt.2-xt.2))
        p2.fisher.vec[itt] = p2
        dec.fisher.mat[itt] = 3
        
      }

      if (p1*p2< (c.fisher*0.98)) sign.mat[itt, 1] = 1
      #if (1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2)) <c.inverse) sign.mat[itt, 2] = 1
    }
    
    ## Make interim look decisions for inverse
    if (pmax.inverse<p.max.cut){
      # stop for futility
      dec.inverse.mat[itt] = 1
    } else {
      if (pn.inverse>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un, rate.c)
        xt.2 = rbinom(1, nt.1.un, rate.t)
        
        p2 = pvalue.chi(c(xc.2, nc.1.un - xc.2), c(xt.2, nt.1.un-xt.2))
        dec.inverse.mat[itt] = 2
        
      } else {
        # continue recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un + nc.2, rate.c)
        xt.2 = rbinom(1, nt.1.un + nt.2, rate.t)
        
        p2 = pvalue.chi(c(xc.2, nc.1.un + nc.2-xc.2), c(xt.2, nt.1.un + nt.2-xt.2))
        dec.inverse.mat[itt] = 3
        
      }
      
      p1 = max(0.0000001, min(p1, 0.9999999))
      p2 = max(0.0000001, min(p2, 0.9999999))
      p2.inverse.vec[itt] = p2
      
      #if (p1*p2< c.fisher) sign.mat[itt, 1] = 1
      if ((1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2))) <(c.inverse*0.99)) sign.mat[itt, 2] = 1
    }
    
  }
  
  ## uniform (0, 1) p1 and p2
  p1.unif = runif(n.itt, 0, 1)
  p2.unif = runif(n.itt, 0, 1)
  fisher.unif = mean((p1.unif<a0.fisher)&(p1.unif*p2.unif<c.fisher))
  inverse.unif = mean((p1.unif<a0.inverse)&((1-pnorm(w1.inv*qnorm(1-p1.unif)+w2.inv*qnorm(1-p2.unif))) <c.inverse))
  dec.unif = c(fisher.unif, inverse.unif)
  
  newlist = list(dec.gold.mat = dec.gold.mat, dec.fisher.mat = dec.fisher.mat, 
                 dec.inverse.mat = dec.inverse.mat, sign.mat = sign.mat, 
                 p1.fisher.vec = p1.fisher.vec, p1.inverse.vec = p1.inverse.vec,
                 p2.fisher.vec = p2.fisher.vec, p2.inverse.vec = p2.inverse.vec,
                 p.gold.vec = p.gold.vec,
                 a0.gold = a0.gold, a0.fisher = a0.fisher, a0.inverse = a0.inverse,
                 c.fisher = c.fisher, c.inverse = c.inverse, c.com = c.com, dec.unif = dec.unif)
  return(newlist)
}

#######################################################################################################

## type 1 error simulations
time.temp = Sys.time()
output.table = NULL

for (scen.ind in 1:24){

  #random.seed = NA
   random.seed = 1
   set.seed(random.seed)

  #print(scen.ind)
  n.itt.in = 10^6
  #n.itt.in = 1000

  if (scen.ind<=6){
    pn.cut.in = 0.8
    p.max.cut.in = 0.1
    nc.1.in = 80
    rate.in = 0.6
    if (scen.ind==1) nc.1.un.in = 15
    if (scen.ind==2) nc.1.un.in = 20
    if (scen.ind==3) nc.1.un.in = 40
    if (scen.ind==4) nc.1.un.in = 60
    if (scen.ind==5) nc.1.un.in = 80
    if (scen.ind==6) nc.1.un.in = 100
  } else if (scen.ind<=12){
    pn.cut.in = 0.9
    p.max.cut.in = 0
    nc.1.in = 80
    rate.in = 0.4
    if (scen.ind==7) nc.1.un.in = 15
    if (scen.ind==8) nc.1.un.in = 20
    if (scen.ind==9) nc.1.un.in = 40
    if (scen.ind==10) nc.1.un.in = 60
    if (scen.ind==11) nc.1.un.in = 80
    if (scen.ind==12) nc.1.un.in = 100
  } else if (scen.ind<=18){
    pn.cut.in = 0.8
    p.max.cut.in = 0.1
    nc.1.in = 80
    nc.1.un.in = 20
    if (scen.ind==13) rate.in = 0.1
    if (scen.ind==14) rate.in = 0.3
    if (scen.ind==15) rate.in = 0.4
    if (scen.ind==16) rate.in = 0.6
    if (scen.ind==17) rate.in = 0.8
    if (scen.ind==18) rate.in = 0.9
  } else if (scen.ind<=24){
    pn.cut.in = 0.8
    p.max.cut.in = 0
    rate.in = 0.6
    nc.1.un.in = 20
    if (scen.ind==19) nc.1.in = 20
    if (scen.ind==20) nc.1.in = 40
    if (scen.ind==21) nc.1.in = 60
    if (scen.ind==22) nc.1.in = 80
    if (scen.ind==23) nc.1.in = 100
    if (scen.ind==24) nc.1.in = 120

  }

  input = list(alpha = 0.025, # chi square test is for two-sided test
               prior.a = 1,
               prior.b = 1,
               pn.cut = pn.cut.in,
               p.max.cut = p.max.cut.in,
               nc.1 = nc.1.in,
               nc.1.un = nc.1.un.in,
               nc.max = 200,
               nt.1 = nc.1.in,
               nt.1.un = nc.1.un.in,
               nt.max = 200,
               rate.c = rate.in,
               rate.t = rate.in,
               n.itt = n.itt.in)

  power.fit = power.cal(input)

  output.table = rbind(output.table,
                       c(apply(power.fit$sign.mat, 2, mean), power.fit$dec.unif,
                         mean(power.fit$dec.gold.mat==1),
                         mean(power.fit$dec.gold.mat==2), mean(power.fit$dec.gold.mat==3),
                         mean(power.fit$dec.fisher.mat==1),
                         mean(power.fit$dec.fisher.mat==2), mean(power.fit$dec.fisher.mat==3),
                         mean(power.fit$dec.inverse.mat==1),
                         mean(power.fit$dec.inverse.mat==2), mean(power.fit$dec.inverse.mat==3),
                         power.fit$a0.gold, power.fit$a0.fisher, power.fit$a0.inverse,
                         min(power.fit$p1.fisher.vec[power.fit$dec.fisher.mat==1], na.rm=TRUE),
                         max(power.fit$p1.fisher.vec[!power.fit$dec.fisher.mat==1], na.rm=TRUE),
                         min(power.fit$p1.inverse.vec[power.fit$dec.inverse.mat==1], na.rm=TRUE),
                         max(power.fit$p1.inverse.vec[!power.fit$dec.inverse.mat==1], na.rm=TRUE),
                         random.seed, input$n.itt
                         )
                       )
  rownames(output.table)[scen.ind] =
    paste("rate", rate.in, "Sn", pn.cut.in, "Fn", p.max.cut.in,
                                           "n1", nc.1.in, "n1.un", nc.1.un.in)
}

colnames(output.table) = c("fisher", "inverse", "goldilocks", "fisher unif", "inverse unif",
                           "gold stop futility", "gold stop recruiting", "gold continue",
                           "fisher stop futility", "fisher stop recruiting", "fisher continue",
                           "inverse stop futility", "inverse stop recruiting", "inverse continue",
                           "gold a0", "fisher a0", "inverse a0",
                           "fisher min p1 coh1", "fisher max p1 coh2",
                           "inverse min p1 coh1", "inverse max p1 coh2",
                           "random seed", "n.itt")

print(Sys.time() - time.temp)
#print(output.table)
write.csv(output.table, paste0("results/subset/prime_type_1_error_com", ".csv"))
#write.csv(output.table, paste0("results/subset/type_1_error_", scen.ind, ".csv"))

## type 1 error simulation ends
###############################################################################################

## power simulations
time.temp = Sys.time()
output.table = NULL

for (scen.ind in 1:33){

  random.seed = 1
  set.seed(random.seed)

  print(scen.ind)
  n.itt.in = 10^6

  if (scen.ind<=11){
    pn.cut.in = 0.8
    p.max.cut.in = 0.1
    nc.1.in = 80
    nc.1.un.in = 15
    
    if (scen.ind <=5 ){ # check type 1 error
      rate.delta = 0
      if (scen.ind==1) rate.c.in = 0.1
      if (scen.ind==2) rate.c.in = 0.3
      if (scen.ind==3) rate.c.in = 0.5
      if (scen.ind==4) rate.c.in = 0.7
      if (scen.ind==5) rate.c.in = 0.9
    } else {
      if (scen.ind==6){rate.c.in = 0.3; rate.delta = 0.143 }
      if (scen.ind==7){rate.c.in = 0.3; rate.delta = 0.17 }
      if (scen.ind==8){rate.c.in = 0.3; rate.delta = 0.188 }
      if (scen.ind==9){rate.c.in = 0.7; rate.delta = 0.127 }
      if (scen.ind==10){rate.c.in = 0.7; rate.delta = 0.145 }
      if (scen.ind==11){rate.c.in = 0.7; rate.delta = 0.16 }
    }
    
  } else if (scen.ind<=22){
    
    pn.cut.in = 0.8
    p.max.cut.in = 0.1
    nc.1.in = 80
    nc.1.un.in = 40
    
    if (scen.ind <=16 ){ # check type 1 error
      rate.delta = 0
      if (scen.ind==12) rate.c.in = 0.1
      if (scen.ind==13) rate.c.in = 0.3
      if (scen.ind==14) rate.c.in = 0.5
      if (scen.ind==15) rate.c.in = 0.7
      if (scen.ind==16) rate.c.in = 0.9
    } else {
      if (scen.ind==17){rate.c.in = 0.4; rate.delta = 0.15 }
      if (scen.ind==18){rate.c.in = 0.4; rate.delta = 0.175 }
      if (scen.ind==19){rate.c.in = 0.4; rate.delta = 0.195 }
      if (scen.ind==20){rate.c.in = 0.6; rate.delta = 0.14 }
      if (scen.ind==21){rate.c.in = 0.6; rate.delta = 0.16 }
      if (scen.ind==22){rate.c.in = 0.6; rate.delta = 0.18 }
    }
  } else if (scen.ind<=33){
    
    pn.cut.in = 0.8
    p.max.cut.in = 0.1
    nc.1.in = 80
    nc.1.un.in = 80
    
    if (scen.ind <=27 ){ # check type 1 error
      rate.delta = 0
      if (scen.ind==23) rate.c.in = 0.1
      if (scen.ind==24) rate.c.in = 0.3
      if (scen.ind==25) rate.c.in = 0.5
      if (scen.ind==26) rate.c.in = 0.7
      if (scen.ind==27) rate.c.in = 0.9
    } else {
      if (scen.ind==28){rate.c.in = 0.4; rate.delta = 0.149 }
      if (scen.ind==29){rate.c.in = 0.4; rate.delta = 0.172 }
      if (scen.ind==30){rate.c.in = 0.4; rate.delta = 0.192 }
      if (scen.ind==31){rate.c.in = 0.6; rate.delta = 0.137 }
      if (scen.ind==32){rate.c.in = 0.6; rate.delta = 0.159 }
      if (scen.ind==33){rate.c.in = 0.6; rate.delta = 0.178 }
    }
  }

  rate.t.in = rate.c.in + rate.delta

  input = list(alpha = 0.025, # chi square test is for two-sided test
               prior.a = 1,
               prior.b = 1,
               pn.cut = pn.cut.in,
               p.max.cut = p.max.cut.in,
               nc.1 = nc.1.in,
               nc.1.un = nc.1.un.in,
               nc.max = 200,
               nt.1 = nc.1.in,
               nt.1.un = nc.1.un.in,
               nt.max = 200,
               rate.c = rate.c.in,
               rate.t = rate.t.in,
               n.itt = n.itt.in)

  power.fit = power.cal(input)

  
  output.table = rbind(output.table,
                       c(apply(power.fit$sign.mat, 2, mean), power.fit$dec.unif, 
                         mean(power.fit$dec.gold.mat==1),
                         mean(power.fit$dec.gold.mat==2), mean(power.fit$dec.gold.mat==3),
                         mean(power.fit$dec.fisher.mat==1),
                         mean(power.fit$dec.fisher.mat==2), mean(power.fit$dec.fisher.mat==3),
                         mean(power.fit$dec.inverse.mat==1),
                         mean(power.fit$dec.inverse.mat==2), mean(power.fit$dec.inverse.mat==3),
                         power.fit$a0.gold, power.fit$a0.fisher, power.fit$a0.inverse, 
                         min(power.fit$p1.fisher.vec[power.fit$dec.fisher.mat==1], na.rm=TRUE),
                         max(power.fit$p1.fisher.vec[!power.fit$dec.fisher.mat==1], na.rm=TRUE),
                         min(power.fit$p1.inverse.vec[power.fit$dec.inverse.mat==1], na.rm=TRUE),
                         max(power.fit$p1.inverse.vec[!power.fit$dec.inverse.mat==1], na.rm=TRUE),
                         random.seed, input$n.itt
                       )
  )
  rownames(output.table)[scen.ind] =
    paste("rate c", rate.c.in, "delta", rate.delta, "Sn", pn.cut.in, "Fn", p.max.cut.in,
          "n1", nc.1.in, "n1.un", nc.1.un.in)
}

colnames(output.table) = c("fisher", "inverse", "goldilocks", "fisher unif", "inverse unif", 
                           "gold stop futility", "gold stop recruiting", "gold continue",
                           "fisher stop futility", "fisher stop recruiting", "fisher continue",
                           "inverse stop futility", "inverse stop recruiting", "inverse continue",
                           "gold a0", "fisher a0", "inverse a0",
                           "fisher min p1 coh1", "fisher max p1 coh2",
                           "inverse min p1 coh1", "inverse max p1 coh2",
                           "random seed", "n.itt")

print(Sys.time() - time.temp)
write.csv(output.table, "results/subset/prime_output_power_com.csv")
# 
# ## power simulation ends


























