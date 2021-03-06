
#scen.ind = 7
setwd("~/goldilocks/")
library(gridExtra)
library(ggplot2)
library(grid)
library(Rmisc)
library(fields)
library(scatterplot3d)

#setwd("~/goldilocks")

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
pred.prob <- function(nc.1.un, nt.1.un, prior.a, prior.b, xc.1, xt.1, nc.1, nt.1, alpha){
  
  # nc.1.un = nt.1.un = 50
  # prior.a = 1
  # prior.b = 1
  # xc.1 = 90
  # xt.1 = 60
  # nc.1 = nt.1 = 100
  # alpha = 0.05
  # 
  # start = sign.chi(n.pbo = c(xc.1, nc.1+nc.1.un-xc.1),
  #                  n.trt = c(xt.1, nt.1+nt.1.un-xt.1), alpha)
  # c1 = c2 = NA
  # if (!start){
  #   c1 = c2 = 0
  # } else if (xc.1/nc.1>xt.1/nt.1){
  #   for (j in 0:nt.1.un){
  #     ind= sign.chi(n.pbo = c(xc.1, nc.1+nc.1.un-xc.1),
  #                   n.trt = c(xt.1+j, nt.1+nt.1.un-xt.1-j), alpha)
  #     if (!ind) {
  #       c1 = 0
  #       c2 = j
  #       break
  #     }
  #     c1 = 0
  #     c2 = j+1
  #   }
  #   
  # } else if (xc.1/nc.1<xt.1/nt.1){
  #   
  #   for (i in 0:nc.1.un){
  #     ind = sign.chi(n.pbo = c(xc.1+i, nc.1+nc.1.un-xc.1-i),
  #                   n.trt = c(xt.1, nt.1+nt.1.un-xt.1), alpha)
  #     if (!ind) {
  #       c1 = i
  #       c2 = 0
  #       break
  #     }
  #     c1 = i+1 
  #     c2 = 0
  #   }
  #   
  # }
  # 
  # 
  # #chisq.c = qchisq(1-alpha, df=1)
  # pn.mat = pn.ind.mat = matrix(0, nrow = nc.1.un+1, ncol = nt.1.un+1)
  # 
  # for (i in 0:nc.1.un){
  #     if (c1 > 0){
  #       if (i<=c1){
  #         j.right = 0
  #         j.left = -1
  #       } else {
  #         j.right = j.left = i-c1
  #       }
  #     } else if (c2 > 0) {
  #       if ((i+c2)<=nt.1.un){
  #         j.right = j.left = i + c2
  #       } else {
  #         j.left = nt.1.un
  #         j.right = nt.1.un+1
  #       }
  #       
  #       
  #     } else if ((c1==0)&(c2==0)){
  #       j.right = j.left = i
  #     }
  #   
  #     ind.right = ind.left = FALSE
  # 
  #     while ((j.right<=nt.1.un)&(!ind.right)){
  #       ind.right = sign.chi(n.pbo = c(xc.1+i, nc.1+nc.1.un-xc.1-i),
  #                      n.trt = c(xt.1+j.right, nt.1+nt.1.un-xt.1-j.right), alpha)
  #       if (!ind.right) j.right = j.right+1
  #     }
  # 
  #     while ((j.left>=0)&(!ind.left)){
  #       ind.left = sign.chi(n.pbo = c(xc.1+i, nc.1+nc.1.un-xc.1-i),
  #                            n.trt = c(xt.1+j.left, nt.1+nt.1.un-xt.1-j.left), alpha)
  #       if (!ind.left) j.left = j.left-1
  #     }
  # 
  #     if (j.left >=0){
  #       for (j in 0:j.left){
  #         pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
  #           beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)
  #       }
  #     }
  # 
  #     if (j.right<=nt.1.un){
  #       for (j in j.right:nt.1.un){
  #         pn.mat[i+1, j+1] = beta.bin.den(prior.a+xc.1, prior.b+nc.1-xc.1, n = nc.1.un, k=i)*
  #           beta.bin.den(prior.a+xt.1, prior.b+nt.1-xt.1, n = nt.1.un, k=j)
  #       }
  #     }
  # 
  # }
  # 
  # pn.mat.1 = pn.mat

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
  nc.max = input$nc.max
  nc.2 = nc.max - nc.1 - nc.1.un
  
  nt.1 = input$nt.1 ## number of trt in check 1
  nt.1.un = input$nt.1.un
  nt.max = input$nt.max
  nt.2 = nt.max - nt.1 - nt.1.un
  
  rate.c = input$rate.c ## rate in control
  rate.t = input$rate.t ## rate in trt
  
  if (p.max.cut==0){
    a0=1
  } else {
    if (!file.exists(paste0("results/cutoff/same/a0_nc.1_", nc.1, 
                            "_nc.max_", nc.max, "_pmax_cut_", p.max.cut, ".csv"))){
      a0 = solve.a0(input)$a0
      write.csv(a0, paste0("results/cutoff/same/a0_nc.1_", nc.1, 
                           "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv"), 
                row.names = FALSE)
    }
    
    a0 = read.csv(paste0("results/cutoff/same/a0_nc.1_", nc.1, 
                                 "_nc.max_", nc.max, "_pmax_cut_", p.max.cut,".csv"))
    a0 = as.numeric(a0)
    #a0=0.6273726
  }
  
  print(a0)
  
  w1.inv = sqrt(0.5)
  w2.inv = sqrt(0.5)
  c.fisher = solve.c(function(x) x+x*log(a0)-x*log(x)-alpha, end = a0)
  print(c.fisher)
  c.com = 0.022
  
  if (!file.exists(paste0("results/cutoff/same/cinverse_w1_", w1.inv, 
                          "_w2_", w2.inv, "_a0_", a0, ".csv"))){
    c.inverse = solve.c.inverse(w1.inv, w2.inv, alpha, a0 = a0, 10^(-7))
    write.csv(c.inverse, paste0("results/cutoff/same/cinverse_w1_", w1.inv, 
                                "_w2_", w2.inv, "_a0_", a0, ".csv"), 
              row.names = FALSE)
  }
  
  c.inverse = read.csv(paste0("results/cutoff/same/cinverse_w1_", w1.inv, 
                       "_w2_", w2.inv, "_a0_", a0, ".csv"))
  c.inverse = as.numeric(c.inverse)
  print(c.inverse)
  
  ############################################ simulations
  n.itt = input$n.itt
  p1.vec = p2.vec = p.com.vec = dec.mat = pn.vec = pmax.vec = rep(NA, n.itt)
  sign.mat = matrix(0, nrow= n.itt, ncol = 3)
  colnames(sign.mat) = c("fisher","inverse","goldilocks")
  
  ## get pn pmax table
  if (!file.exists(paste0("results/cutoff/same/pn_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                         "_nc.max_", nc.max, ".csv"))){
    solve.pn.table(input$nc.1, input$nt.1, input$nc.1.un, input$nt.1.un,
                   input$nc.max, input$nt.max, input$prior.a, input$prior.b, input$alpha)
  }
  
  pn.c.table = read.csv(paste0("results/cutoff/same/pn_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                               "_nc.max_", nc.max, ".csv"))
  pmax.c.table = read.csv(paste0("results/cutoff/same/pmax_nc.1_", nc.1, "_nc.1.un_", nc.1.un, 
                                 "_nc.max_", nc.max, ".csv"))
  
  for (itt in 1:n.itt){
    if ((itt %% 1000)==0) print(paste("scen", scen.ind, 'itt', itt))
    ## observed events in stage 1. 
    xc.1 = rbinom(1, nc.1, rate.c)
    xt.1 = rbinom(1, nt.1, rate.t)
    
    pn = max(min(1, pn.c.table[xc.1+1, xt.1+1]), 0)
    p.max = max(min(1, pmax.c.table[xc.1+1, xt.1+1]), 0)
    
    p1 = p2 = NA
    p1 = pvalue.chi(c(xc.1, nc.1-xc.1), c(xt.1, nt.1-xt.1))
    p1.vec[itt] = p1
    
    ### naive adaptation
    # pn = 0.1
    # ptemp = pvalue.chi(c(xc.1, nc.1-xc.1), c(xt.1, nt.1-xt.1))
    # if (ptemp<0.05) pn = 0.9
    # p.max = 0.5

    ## Make interim look decisions
    if (p.max<p.max.cut){
      # stop for futility
      dec.mat[itt] = 1
    } else {

      if (pn>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un, rate.c)
        xt.2 = rbinom(1, nt.1.un, rate.t)

        p2 = pvalue.chi(c(xc.2, nc.1.un - xc.2), c(xt.2, nt.1.un-xt.2))
        p2.vec[itt] = p2
        dec.mat[itt] = 2
        
        p.com = pvalue.chi(c(xc.1 + xc.2, nc.1+nc.1.un -xc.1-xc.2), 
                           c(xt.1 + xt.2, nt.1+nt.1.un -xt.1-xt.2))
        p.com.vec[itt] = p.com

      } else {
        # continue recruiting
        ## observed events in stage 2.
        xc.2 = rbinom(1, nc.1.un + nc.2, rate.c)
        xt.2 = rbinom(1, nt.1.un + nt.2, rate.t)

        p2 = pvalue.chi(c(xc.2, nc.1.un + nc.2-xc.2), c(xt.2, nt.1.un + nt.2-xt.2))
        p2.vec[itt] = p2
        dec.mat[itt] = 3
        
        p.com = pvalue.chi(c(xc.1 + xc.2, nc.1+nc.1.un+nc.2 -xc.1-xc.2), 
                           c(xt.1 + xt.2, nt.1+nt.1.un+nt.2 -xt.1-xt.2))
        p.com.vec[itt] = p.com
      }

      ## check if significant
      #print(p1)
      #print(p2)
      
      if (p1*p2< c.fisher) sign.mat[itt, 1] = 1
      if (1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2)) <c.inverse) sign.mat[itt, 2] = 1
      if (p.com < c.com) sign.mat[itt, 3] = 1
    }
    
  }
  
  ## uniform (0, 1) p1 and p2
  p1.unif = runif(n.itt, 0, 1)
  p2.unif = runif(n.itt, 0, 1)
  fisher.unif = mean((p1.unif<a0)&(p1.unif*p2.unif<c.fisher))
  inverse.unif = mean((p1.unif<a0)&((1-pnorm(w1.inv*qnorm(1-p1.unif)+w2.inv*qnorm(1-p2.unif))) <c.inverse))
  dec.unif = c(fisher.unif, inverse.unif)
  
  newlist = list(dec.mat = dec.mat, sign.mat = sign.mat, p1 = p1.vec, 
                 p2 = p2.vec, a0 = a0, pn.vec = pn.vec, pmax.vec = pmax.vec,
                 c.fisher = c.fisher, c.inverse = c.inverse, c.com = c.com, dec.unif = dec.unif)
  return(newlist)
}

#######################################################################################################

## power simulations
time.temp = Sys.time()
output.table = NULL

for (scen.ind in c(6, 28)){

  
  print(scen.ind)
  n.itt.in = 10^4

  if (scen.ind<=11){
    pn.cut.in = 0.8
    p.max.cut.in = 0
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
      if (scen.ind==6){rate.c.in = 0.3; rate.delta = 0.15 }
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
    p.max.cut.in = 0
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
      if (scen.ind==28){rate.c.in = 0.3; rate.delta = 0.15 }
      if (scen.ind==29){rate.c.in = 0.4; rate.delta = 0.172 }
      if (scen.ind==30){rate.c.in = 0.4; rate.delta = 0.192 }
      if (scen.ind==31){rate.c.in = 0.6; rate.delta = 0.137 }
      if (scen.ind==32){rate.c.in = 0.6; rate.delta = 0.159 }
      if (scen.ind==33){rate.c.in = 0.6; rate.delta = 0.178 }
    }
  }

  rate.t.in = rate.c.in + rate.delta

  random.seed = 1
  set.seed(random.seed)
  
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
                       c(apply(power.fit$sign.mat, 2, mean), mean(power.fit$dec.mat==1),
                         mean(power.fit$dec.mat==2), mean(power.fit$dec.mat==3),
                         power.fit$a0, min(power.fit$p1[power.fit$dec.mat==1], na.rm=TRUE),
                         max(power.fit$p1[!power.fit$dec.mat==1], na.rm=TRUE),
                         random.seed, input$n.itt
                       )
  )
  print(output.table)

  ## make power plot
  p1.vec = power.fit$p1
  p2.vec = power.fit$p2
  fisher.dec = power.fit$sign.mat[,1]
  inverse.dec = power.fit$sign.mat[,2]
  c.fisher = power.fit$c.fisher
  c.inverse = power.fit$c.inverse
  
  p1.plot = p2.plot = c(seq(0, 0.1, length=4000), seq(0.10000001, 1, length=1000))
  
  library(fields)
  library(scatterplot3d)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  
  fisher.f = function(x, y) { x*y }
  fisher.f.reverse = function(x, com){com/x}
  p.fisher = outer(p1.plot, p2.plot, fisher.f)
  cut.fisher = cbind(p1.plot, fisher.f.reverse(p1.plot, c.fisher))
  cut.fisher = data.frame(cut.fisher)
  
  inverse.f = function(x, y) {
    w1 = w2 = sqrt(1/2)
    z = w1*qnorm(1-x) + w2*qnorm(1-y)
    pnorm(z, lower.tail = FALSE)
  }
  
  inverse.f.reverse = function(x, com){
    w1 = w2 = sqrt(1/2)
    z = (qnorm(1-com) - w1*qnorm(1-x))/w2
    pnorm(z, lower.tail = FALSE)
  }
  
  p.inverse = outer(p1.plot, p2.plot, inverse.f)
  cut.inverse = cbind(p1.plot, inverse.f.reverse(p1.plot, c.inverse))
  cut.inverse = data.frame(cut.inverse)
  
  # 1: both non-sig, 2: fisher sig, 3: inverse sig, 4: both sig
  plot.DF = data.frame(p1 = p1.vec, p2 = p2.vec, legend = 3)
  plot.DF[(fisher.dec==1)&(inverse.dec==0),3] = 6
  plot.DF[(fisher.dec==0)&(inverse.dec==1),3] = 5
  plot.DF[(fisher.dec==1)&(inverse.dec==1),3] = 4
  
  cut.DF = data.frame(p1 = c(cut.fisher$p1, cut.inverse$p1), p2 = c(cut.fisher$V2, cut.inverse$V2),
                      legend = c(rep(1, dim(cut.fisher)[1]), rep(2, dim(cut.inverse)[1])))
  all.DF = rbind(plot.DF, cut.DF)
  all.DF$legend = factor(all.DF$legend)
  
  order.vec = c(3, 4, 5, 1, 2, 6)
  all.DF = all.DF[order(match(all.DF$legend, order.vec)),]
  

  gg.title = paste0("Fisher only: ", mean(plot.DF[,3]==6)*100,
                    "%, Inverse only: ", mean(plot.DF[,3]==5)*100, "%")
  
  if (scen.ind==6) show.legend.ind = TRUE
  if (scen.ind==28) show.legend.ind = TRUE
  
  gg = ggplot(all.DF, aes(x=p1, y=p2, color=legend, size=legend))+
    geom_point(show.legend = show.legend.ind) +
    scale_color_manual(labels = c("Fisher c", "Inverse c", "Both not significant", 
                                  "Both significant", "Inverse significant", "Fisher significant"
    ),
    values=c("springgreen2", "red", 'snow3',  "tan4", "yellow1", "blue"))+
    scale_size_manual(values=c(rep(0.8, 2), rep(3, 4)), guide=FALSE)+
    xlim(0, 1) + ylim(0, 1) +
    theme_bw() + 
    #ggtitle(gg.title)+
    theme(text = element_text(size=43),
          axis.text.x = element_text(colour="black",size=50,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=50,angle=0,hjust=1,vjust=0,face="plain"),
          axis.title.x = element_text(colour="black",size=50,angle=0,hjust=.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="black",size=50,angle=0,hjust=.5,vjust=.5,face="plain"),
          legend.text = element_text(colour="black", size = 35, face = "plain"),
          legend.title = element_text(colour="black", size = 35, face = "plain"),
          legend.key.size = unit(2,"line"),
          legend.position="bottom", plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(1,1,1.5,1), "cm")
          )+
    guides(colour = guide_legend(override.aes = list(size=6)))
  
  if (scen.ind==6) {gg1 = gg; plot.DF.1 = plot.DF}
  if (scen.ind==28) {gg2 = gg; plot.DF.2 = plot.DF}
  #print(gg)
  
}

gglist = list(gg1, gg2)

png("results/combination_plot.png", width = 2400, height = 1200)
multiplot(plotlist = gglist, cols=2)
dev.off()

print(Sys.time() - time.temp)


## power simulation ends




































