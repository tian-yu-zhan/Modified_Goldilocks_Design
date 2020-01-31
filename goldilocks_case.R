
# case study: TTE with OS and log rank test
#add ppn^prime
# setwd("/Users/zhantx/Documents/code/goldilocks/")
setwd("~/goldilocks")
library(survival)
library(doParallel)


n.cluster = 119

###### sample survival time from haz
surv.sample = function(n, haz){
  temp = runif(n, 0, 1)
  t.temp = -log(temp)/haz
}

###### log rank p value
log.rank.p.func = function(surv.fit){
  log.rank = sqrt(surv.fit$chisq)
  z1.diff= surv.fit$obs - surv.fit$exp
  if (z1.diff[1]<=z1.diff[2]) log.rank = -log.rank
  
  return(pnorm(log.rank, lower.tail = FALSE))
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

#### function to calculate power
power.cal = function(input){
  alpha = input$alpha
  
  # prior.c.a = input$prior.c.a
  # prior.c.b = input$prior.c.b
  
  pn.cut = input$pn.cut
  p.max.cut = input$p.max.cut
  
  nc.1 = input$nc.1 ## number of pbo in check 1
  enrol.rate.c = input$enrol.rate.c ## unobserved observations
  nc.max = input$nc.max
  
  nt.1 = input$nt.1 ## number of trt in check 1
  enrol.rate.t = input$enrol.rate.t
  nt.max = input$nt.max
  
  haz.c = input$haz.c ## rate in control
  haz.t = input$haz.t ## rate in trt
  
  c.fisher.a0 = solve.c(function(x) x+x*log(1)-x*log(x)-alpha, end = 1)
  w1.inv = w2.inv = sqrt(0.5)
  
  if (!file.exists(paste0("results/case/cinverse_w1_", w1.inv, 
                          "_w2_", w2.inv, "_a0_", 1, ".csv"))){
    c.inverse.a0 = solve.c.inverse(w1.inv, w2.inv, alpha, a0 = 1, 10^(-7))
    write.csv(c.inverse.a0, paste0("results/case/cinverse_w1_", w1.inv, 
                                "_w2_", w2.inv, "_a0_", 1, ".csv"), 
              row.names = FALSE)
  }
  
  c.inverse.a0 = as.numeric(read.csv(paste0("results/case/cinverse_w1_", w1.inv, 
                              "_w2_", w2.inv, "_a0_", 1, ".csv")))
  print(c.inverse.a0)
  
  
  # a0.gold = a0.fisher = a0.inverse = 1
  a0.gold = 1
  
  # haz.c.pull = (haz.c*nc.1 + haz.t*nt.1)/(nc.1+nt.1)
  
  if (p.max.cut==0){
    a0.fisher = a0.inverse = 1
  } else {
    
    if (!file.exists(paste0("results/case/a0_fisher_nc.1_", nc.1, 
                            "_nc.max_", nc.max, "_enrol.rate.c_", enrol.rate.c, 
                            "_pmax_cut_", p.max.cut, 
                            "_haz.c_", haz.c, "_haz.t_", haz.t, ".csv"))){
      
    pmax.fisher.a0.vec = pmax.inverse.a0.vec =
      p1.a0.vec = rep(0, 10^4)
  
    cl <- makeCluster(n.cluster)
    registerDoParallel(cl)
    pred = foreach(i = 1:length(pmax.fisher.a0.vec)) %dopar% {
      
      ###### sample survival time from haz
      surv.sample = function(n, haz){
        temp = runif(n, 0, 1)
        t.temp = -log(temp)/haz
      }
      
      ###### log rank p value
      log.rank.p.func = function(surv.fit){
        log.rank = sqrt(surv.fit$chisq)
        z1.diff= surv.fit$obs - surv.fit$exp
        if (z1.diff[1]<=z1.diff[2]) log.rank = -log.rank
        
        return(pnorm(log.rank, lower.tail = FALSE))
      }
      
    #for (i in 1:length(pmax.fisher.a0.vec)){
      # print(i)
      library(survival)
      haz.c.1 = haz.c
      haz.t.1 = haz.c
      
      # temp = rep(0, 10^4)
      # nc.max = nt.max = 300
      # nc.1 = nt.1 = 200
      # 
      # for (i in 1:length(temp)){
      # print(i)
      
      surv.sample.c = surv.sample(nc.max, haz.c.1)
      surv.sample.t = surv.sample(nt.max, haz.t.1)
      
      data.stage = data.frame("start" = c((1:nc.max)/nc.max*(nc.max/enrol.rate.c),
                                          (1:nt.max)/nt.max*(nt.max/enrol.rate.t)),
                              "TTE" = c(surv.sample.c, surv.sample.t),
                              "group" = c(rep("C", nc.max), rep("T", nt.max)),
                              "stage.1.ind" = c(rep(1, nc.1), rep(0, nc.max-nc.1),
                                                rep(1, nt.1), rep(0, nt.max-nt.1))
      )
      
      data.stage$TTE_calen = data.stage$start + data.stage$TTE
      data.stage$cen = data.stage$start + 12
      data.stage$t_calend = pmin( data.stage$TTE_calen, data.stage$cen)
      
      data.stage$TTE_1 = pmin(data.stage$TTE, 12)
      data.stage$event_1 = as.numeric(data.stage$TTE<=12)
      
      data.stage$cut = nc.1/enrol.rate.c
      data.stage$exp = pmin(data.stage$t_calend, data.stage$cut)-data.stage$start
      data.stage$obs_1 = data.stage$t_calend<=data.stage$cut
      data.stage$event_cut = (data.stage$obs_1==TRUE)&(data.stage$event_1==1)
      
      # data.stage$cen_1 = pmin(data.stage$cen, nc.1/enrol.rate.c)
      
      
      ####################################################################
      ## stage 1 pvalue
      data.stage.1 = data.stage[(data.stage$stage.1.ind==1)&(data.stage$obs_1==TRUE),]
      data.stage.cut.1 = data.stage[data.stage$stage.1.ind==1,]
      
      log.rank.fit.1 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.stage.1, rho=0)
      #log.cox.fit.1 = coxph(Surv(TTE_1, event_1) ~ group, data = data.stage.1)
      #p1 = pnorm(summary(log.cox.fit.1)$coefficients[4], lower.tail = FALSE)
      
      log.rank.p.1 = log.rank.p.func(log.rank.fit.1)
      log.rank.p.1 = max(0.0000001, min(log.rank.p.1, 0.9999999))
      # p1.fisher.vec[itt] = p1.inverse.vec[itt] = log.rank.p.1
      # p1.a0.vec[i] = log.rank.p.1
      
      pmax.vec = pmax.fisher.vec = pmax.inv.vec = rep(0, 10^3)
      for (pmax.vec.ind in 1:length(pmax.vec)){
        # print(pmax.vec.ind)
        haz.c.pred = rgamma(1, shape=haz.c*10+sum(data.stage.1$event_cut[(data.stage.1$group=="C")]),
                            rate = 10+sum(data.stage.1$exp[data.stage.1$group=="C"]))
        haz.t.pred = rgamma(1, shape=haz.t*10+sum(data.stage.1$event_cut[data.stage.1$group=="T"]),
                            rate = 10+sum(data.stage.1$exp[data.stage.1$group=="T"]))
        
        surv.sample.c.pred = surv.sample(nc.max, haz.c.pred)
        surv.sample.t.pred = surv.sample(nt.max, haz.t.pred)
        
        data.stage$TTE_pred_init = c(surv.sample.c.pred, surv.sample.t.pred)
        # data.stage$TTE_pred = data.stage$TTE_pred_init + data.stage$TTE
        data.stage$TTE_pred = data.stage$TTE_pred_init
        
        data.stage$TTE_pred_f = pmin(data.stage$TTE_pred, 12)
        data.stage$event_pred = data.stage$TTE_pred <= 12
        
        data.stage$TTE_pred_f[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
          data.stage$TTE_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]
        data.stage$event_pred[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
          data.stage$event_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]
        
        ##### pmax pmax for GD
        data.pmax = data.stage
        log.rank.pmax = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.pmax)
        log.rank.p.pmax = log.rank.p.func(log.rank.pmax)
        if (log.rank.p.pmax<=0.025) pmax.vec[pmax.vec.ind]=1
        
        ## pmax for fisher and inverse
        data.com.pmax = data.stage[(data.stage$stage.1.ind==0)|(data.stage$obs_1=="FALSE"),]
        log.rank.com.pmax = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.com.pmax)
        log.rank.p2.pmax = log.rank.p.func(log.rank.com.pmax)
        
        log.rank.p2.pmax = max(0.0000001, min(log.rank.p2.pmax, 0.9999999))
        
        if (log.rank.p.1*log.rank.p2.pmax< (c.fisher.a0)) pmax.fisher.vec[pmax.vec.ind] = 1
        #if (log.rank.p.1 > a0.fisher) pmax.fisher.vec[pmax.vec.ind] = 0
        
        if ((1-pnorm(w1.inv*qnorm(1-log.rank.p.1)+w2.inv*qnorm(1-log.rank.p2.pmax))) <
            (c.inverse.a0)) pmax.inv.vec[pmax.vec.ind] = 1
        # if (log.rank.p.1 > a0.inverse) pmax.inv.vec[pmax.vec.ind] = 0
      }
      # pmax.gold = mean(pmax.vec)
      # pmax.fisher.a0.vec[i] = mean(pmax.fisher.vec)
      # pmax.inverse.a0.vec[i] = mean(pmax.inv.vec)
      
      newlist = list("p1.a0" = log.rank.p.1,
                     "pmax.fisher.a0" = mean(pmax.fisher.vec),
                     "pmax.inverse.a0" = mean(pmax.inv.vec))
      return(newlist)
    }
    
    stopCluster(cl)
    
    for (i in 1:length(pmax.fisher.a0.vec)){
      
      pred.temp = pred[[i]]
      p1.a0.vec[i] = pred.temp$p1.a0
      pmax.fisher.a0.vec[i] = pred.temp$pmax.fisher.a0
      pmax.inverse.a0.vec[i] = pred.temp$pmax.inverse.a0
    }
    
    a0.fisher.1 = max(p1.a0.vec[pmax.fisher.a0.vec>p.max.cut])
    a0.fisher.2 = min(p1.a0.vec[pmax.fisher.a0.vec<=p.max.cut])
    
    a0.fisher = a0.fisher.1
    #a0.fisher = mean(c(a0.fisher.1, a0.fisher.2))
    # a0.fisher = max(a0.fisher.2 +  
    #           mean((p1.a0.vec>a0.fisher.2)&(pmax.fisher.a0.vec>p.max.cut)),
    #                 a0.fisher.1 - 
    #   mean((p1.a0.vec<=a0.fisher.1)&(pmax.fisher.a0.vec>p.max.cut)) )
                
    
    a0.inverse.1 = max(p1.a0.vec[pmax.inverse.a0.vec>p.max.cut])
    a0.inverse.2 = min(p1.a0.vec[pmax.inverse.a0.vec<=p.max.cut])
    
    a0.inverse =  a0.inverse.1
    # a0.inverse = mean(c(a0.inverse.1, a0.inverse.2))
    # a0.inverse = max(a0.inverse.2 + 
    #               mean((p1.a0.vec>a0.inverse.2)&(pmax.inverse.a0.vec>p.max.cut)),
    #                  a0.inverse.1 - 
    #   mean((p1.a0.vec<=a0.inverse.1)&(pmax.inverse.a0.vec>p.max.cut)) )
      
    
    #a0.fisher = mean(pmax.fisher.a0.vec>p.max.cut)
    #a0.inverse = mean(pmax.inverse.a0.vec>p.max.cut)
    
    # a0.fisher.init = max(p1.a0.vec[pmax.fisher.a0.vec>p.max.cut])
    # a0.inverse.init = max(p1.a0.vec[pmax.inverse.a0.vec>p.max.cut])

    write.csv(a0.fisher, paste0("results/case/a0_fisher_nc.1_", nc.1, 
                                "_nc.max_", nc.max, "_enrol.rate.c_", enrol.rate.c, 
                                "_pmax_cut_", p.max.cut, 
                                "_haz.c_", haz.c, "_haz.t_", haz.t, ".csv"), 
              row.names = FALSE)
    write.csv(a0.inverse, paste0("results/case/a0_inverse_nc.1_", nc.1, 
                                 "_nc.max_", nc.max, "_enrol.rate.c_", enrol.rate.c, 
                                 "_pmax_cut_", p.max.cut, 
                                 "_haz.c_", haz.c, "_haz.t_", haz.t, ".csv"), 
              row.names = FALSE)
    }
    
    a0.fisher = as.numeric(read.csv(paste0("results/case/a0_fisher_nc.1_", nc.1, 
                                           "_nc.max_", nc.max, "_enrol.rate.c_", enrol.rate.c, 
                                           "_pmax_cut_", p.max.cut, 
                                           "_haz.c_", haz.c, "_haz.t_", haz.t, ".csv")))
    a0.inverse = as.numeric(read.csv(paste0("results/case/a0_inverse_nc.1_", nc.1, 
                                            "_nc.max_", nc.max, "_enrol.rate.c_", enrol.rate.c, 
                                            "_pmax_cut_", p.max.cut, 
                                            "_haz.c_", haz.c, "_haz.t_", haz.t, ".csv")))

  }
  
  print("Done: a0.fisher a0.inverse")
    
  # if (input$a0.fisher.force.ind){
  #   a0.fisher = input$a0.fisher.force
  # }
  # if (input$a0.inverse.force.ind){
  #   a0.inverse = input$a0.inverse.force
  # }
  

  print(a0.gold)
  print(a0.fisher)
  print(a0.inverse)
  
  w1.inv = sqrt(0.5)
  w2.inv = sqrt(0.5)
  c.fisher = solve.c(function(x) x+x*log(a0.fisher)-x*log(x)-alpha, end = a0.fisher)
  print(c.fisher)
  c.com = 0.021
  # if (enrol.rate.c==10) c.com = 0.02
  # if (enrol.rate.c==5) c.com = 0.02
  
  if (!file.exists(paste0("results/case/cinverse_w1_", w1.inv, 
                          "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"))){
    c.inverse = solve.c.inverse(w1.inv, w2.inv, alpha, a0 = a0.inverse, 10^(-7))
    write.csv(c.inverse, paste0("results/case/cinverse_w1_", w1.inv, 
                                "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"), 
              row.names = FALSE)
  }
  
  c.inverse = read.csv(paste0("results/case/cinverse_w1_", w1.inv, 
                       "_w2_", w2.inv, "_a0_", a0.inverse, ".csv"))
  c.inverse = as.numeric(c.inverse)
  print(c.inverse)
  
  print("Done: c.inverse")
  
  ############################################ simulations
  n.itt = input$n.itt
  p1.fisher.vec = p2.fisher.vec = p1.inverse.vec = p2.inverse.vec = p.gold.vec = rep(NA, n.itt)
  dec.fisher.mat = dec.inverse.mat = dec.gold.mat = rep(NA, n.itt)
  
  sign.mat = matrix(0, nrow= n.itt, ncol = 3)
  colnames(sign.mat) = c("fisher","inverse","goldilocks")
  
  # n.itt = 10
  # for (itt in 1:n.itt){
  cl <- makeCluster(n.cluster)
  registerDoParallel(cl)
  pred = foreach(itt = 1:n.itt) %dopar% {
    library(survival)
    
    ###### sample survival time from haz
    surv.sample = function(n, haz){
      temp = runif(n, 0, 1)
      t.temp = -log(temp)/haz
    }
    
    ###### log rank p value
    log.rank.p.func = function(surv.fit){
      log.rank = sqrt(surv.fit$chisq)
      z1.diff= surv.fit$obs - surv.fit$exp
      if (z1.diff[1]<=z1.diff[2]) log.rank = -log.rank
      
      return(pnorm(log.rank, lower.tail = FALSE))
    }
    
    sign.mat.cluster = rep(0, 3)
    # print(itt)
    # if ((itt %% 1000)==0) print(paste("scen", scen.ind, 'itt', itt))
    ## observed events in stage 1. 

    # haz.c = 0.1
    # haz.t = 0.07
    
    # haz.c.1 = rgamma(1, shape=haz.c*10, rate = 10)
    # haz.t.1 = rgamma(1, shape=haz.t*10, rate = 10)
    
    haz.c.1 = haz.c
    haz.t.1 = haz.t
    
    # temp = rep(0, 10^4)
    # nc.max = nt.max = 300
    # nc.1 = nt.1 = 200
    # 
    # for (i in 1:length(temp)){
    # print(i)
    
    surv.sample.c = surv.sample(nc.max, haz.c.1)
    surv.sample.t = surv.sample(nt.max, haz.t.1)
    
    data.stage = data.frame("start" = c((1:nc.max)/nc.max*(nc.max/enrol.rate.c),
                                          (1:nt.max)/nt.max*(nt.max/enrol.rate.t)),
                              "TTE" = c(surv.sample.c, surv.sample.t),
                              "group" = c(rep("C", nc.max), rep("T", nt.max)),
                              "stage.1.ind" = c(rep(1, nc.1), rep(0, nc.max-nc.1),
                                                rep(1, nt.1), rep(0, nt.max-nt.1))
                                )
    
    data.stage$TTE_calen = data.stage$start + data.stage$TTE
    data.stage$cen = data.stage$start + 12
    data.stage$t_calend = pmin( data.stage$TTE_calen, data.stage$cen)
    
    data.stage$TTE_1 = pmin(data.stage$TTE, 12)
    data.stage$event_1 = as.numeric(data.stage$TTE<=12)
    
    data.stage$cut = nc.1/enrol.rate.c
    data.stage$exp = pmin(data.stage$t_calend, data.stage$cut)-data.stage$start
    data.stage$obs_1 = data.stage$t_calend<=data.stage$cut
    data.stage$event_cut = (data.stage$obs_1==TRUE)&(data.stage$event_1==1)
    
    # data.stage$cen_1 = pmin(data.stage$cen, nc.1/enrol.rate.c)

    
    ####################################################################
    ## stage 1 pvalue
    data.stage.1 = data.stage[(data.stage$stage.1.ind==1)&(data.stage$obs_1==TRUE),]
    data.stage.cut.1 = data.stage[data.stage$stage.1.ind==1,]
    
    print(dim(data.stage.1))
    
    log.rank.fit.1 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.stage.1, rho=0)
    #log.cox.fit.1 = coxph(Surv(TTE_1, event_1) ~ group, data = data.stage.1)
    #p1 = pnorm(summary(log.cox.fit.1)$coefficients[4], lower.tail = FALSE)
    
    log.rank.p.1 = log.rank.p.func(log.rank.fit.1)
    log.rank.p.1 = max(0.0000001, min(log.rank.p.1, 0.9999999))
    p1.fisher.cluster = p1.inverse.cluster = log.rank.p.1
    p1 = log.rank.p.1
    
    # p1 = pchisq(log.rank.fit.1$chisq, df = 1, lower.tail=FALSE)
    
    # temp[i] = p1
    # }
    
    
    
    
    
    
    ###################################################################
    ## pn 
    # pn.gold = pn.fisher = pn.inverse = 1
    # if (p1 <=0.025){
    #   pn.gold = pn.fisher = pn.inverse = 0.8
    # }
    # 
    # pmax.gold = pmax.fisher = pmax.inverse = 1
    # if (p1 >0.4){
    #   pmax.gold = pmax.fisher = pmax.inverse = 0.05
    # }
    
    pn.vec = pn.fisher.vec = pn.inv.vec = rep(0, 10^3)

    for (pn.vec.ind in 1:length(pn.vec)){
      # print(pn.vec.ind)
      haz.c.pred = rgamma(1, shape=haz.c*10+sum(data.stage.1$event_cut[(data.stage.1$group=="C")]),
                          rate = 10+sum(data.stage.1$exp[data.stage.1$group=="C"]))
      haz.t.pred = rgamma(1, shape=haz.t*10+sum(data.stage.1$event_cut[data.stage.1$group=="T"]),
                          rate = 10+sum(data.stage.1$exp[data.stage.1$group=="T"]))

      surv.sample.c.pred = surv.sample(nc.max, haz.c.pred)
      surv.sample.t.pred = surv.sample(nt.max, haz.t.pred)
      
      data.stage$TTE_pred_init = c(surv.sample.c.pred, surv.sample.t.pred)
      # data.stage$TTE_pred = data.stage$TTE_pred_init + data.stage$TTE
      data.stage$TTE_pred = data.stage$TTE_pred_init
      
      data.stage$TTE_pred_f = pmin(data.stage$TTE_pred, 12)
      data.stage$event_pred = data.stage$TTE_pred <= 12

      data.stage$TTE_pred_f[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
        data.stage$TTE_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]
      data.stage$event_pred[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
        data.stage$event_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]

      ##### pn for GD
      data.pn = data.stage[data.stage$stage.1.ind==1,]
      log.rank.pn = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.pn)
      log.rank.p.pn = log.rank.p.func(log.rank.pn)
      if (log.rank.p.pn<=0.025) pn.vec[pn.vec.ind]=1

      ## pn for fisher and inverse
      data.com.pn = data.stage[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="FALSE"),]
      log.rank.com.pn = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.com.pn)
      log.rank.p2.pn = log.rank.p.func(log.rank.com.pn)

      log.rank.p2.pn = max(0.0000001, min(log.rank.p2.pn, 0.9999999))

      if (log.rank.p.1*log.rank.p2.pn< (c.fisher.a0)) pn.fisher.vec[pn.vec.ind] = 1
      # if (log.rank.p.1 > a0.fisher) pn.fisher.vec[pn.vec.ind] = 0

      if ((1-pnorm(w1.inv*qnorm(1-log.rank.p.1)+w2.inv*qnorm(1-log.rank.p2.pn))) <
          (c.inverse.a0)) pn.inv.vec[pn.vec.ind] = 1
      # if (log.rank.p.1 > a0.inverse) pn.inv.vec[pn.vec.ind] = 0

    }
    pn.gold = mean(pn.vec)
    pn.fisher = mean(pn.fisher.vec)
    pn.inverse = mean(pn.inv.vec)

    ## pmax
    pmax.vec = pmax.fisher.vec = pmax.inv.vec = rep(0, 10^3)
    for (pmax.vec.ind in 1:length(pmax.vec)){
      # print(pmax.vec.ind)
      haz.c.pred = rgamma(1, shape=haz.c*10+sum(data.stage.1$event_cut[(data.stage.1$group=="C")]),
                          rate = 10+sum(data.stage.1$exp[data.stage.1$group=="C"]))
      haz.t.pred = rgamma(1, shape=haz.t*10+sum(data.stage.1$event_cut[data.stage.1$group=="T"]),
                          rate = 10+sum(data.stage.1$exp[data.stage.1$group=="T"]))
      
      surv.sample.c.pred = surv.sample(nc.max, haz.c.pred)
      surv.sample.t.pred = surv.sample(nt.max, haz.t.pred)

      data.stage$TTE_pred = c(surv.sample.c.pred, surv.sample.t.pred)
      data.stage$TTE_pred_f = pmin(data.stage$TTE_pred, 12)
      data.stage$event_pred = data.stage$TTE_pred <= 12

      data.stage$TTE_pred_f[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
        data.stage$TTE_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]
      data.stage$event_pred[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")] =
        data.stage$event_1[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="TRUE")]

      ##### pmax pmax for GD
      data.pmax = data.stage
      log.rank.pmax = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.pmax)
      log.rank.p.pmax = log.rank.p.func(log.rank.pmax)
      if (log.rank.p.pmax<=0.025) pmax.vec[pmax.vec.ind]=1

      ## pmax for fisher and inverse
      data.com.pmax = data.stage[(data.stage$stage.1.ind==0)|(data.stage$obs_1=="FALSE"),]
      log.rank.com.pmax = survdiff(Surv(TTE_pred_f, event_pred) ~ group, data = data.com.pmax)
      log.rank.p2.pmax = log.rank.p.func(log.rank.com.pmax)

      log.rank.p2.pmax = max(0.0000001, min(log.rank.p2.pmax, 0.9999999))

      if (log.rank.p.1*log.rank.p2.pmax< (c.fisher.a0)) pmax.fisher.vec[pmax.vec.ind] = 1
      # if (log.rank.p.1 > a0.fisher) pmax.fisher.vec[pmax.vec.ind] = 0

      if ((1-pnorm(w1.inv*qnorm(1-log.rank.p.1)+w2.inv*qnorm(1-log.rank.p2.pmax))) <
          (c.inverse.a0)) pmax.inv.vec[pmax.vec.ind] = 1
      # if (log.rank.p.1 > a0.inverse) pmax.inv.vec[pmax.vec.ind] = 0
    }
    pmax.gold = mean(pmax.vec)
    pmax.fisher = mean(pmax.fisher.vec)
    pmax.inverse = mean(pmax.inv.vec)
    
    
    ## Make interim look decisions for gold
    if (pmax.gold<p.max.cut){
      # stop for futility
      dec.gold.cluster = 1
      p.gold.cluster = NA
    } else {
      if (pn.gold>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        dec.gold.cluster = 2
        
        data.pn.real = data.stage[data.stage$stage.1.ind==1,]
        
        log.rank.pn.final = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pn.real)
        p.com = log.rank.p.func(log.rank.pn.final)
        
        p.gold.cluster = p.com
        
      } else {
        # continue recruiting
        ## observed events in stage 2.
        dec.gold.cluster = 3
        
        data.pmax.real = data.stage
        
        log.rank.pmax.final = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pmax.real)
        p.com = log.rank.p.func(log.rank.pmax.final)
        
        p.gold.cluster = p.com
      }
      
      if (p.com < c.com) sign.mat.cluster[3] = 1
    }
    
    
    
    ## Make interim look decisions for fisher
    if (pmax.fisher<p.max.cut){
      # stop for futility
      dec.fisher.cluster = 1
      p2.fisher.cluster = NA
    } else {
      if (pn.fisher>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        data.pn.2 = data.stage[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="FALSE"),]
        
        log.rank.pn.2 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pn.2)
        p2 = log.rank.p.func(log.rank.pn.2)
        
        p2.fisher.cluster = p2
        dec.fisher.cluster = 2

      } else {
        # continue recruiting
        ## observed events in stage 2.
        data.pmax.2 = data.stage[(data.stage$stage.1.ind==0)|(data.stage$obs_1=="FALSE"),]
        
        log.rank.pmax.2 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pmax.2)
        p2 = log.rank.p.func(log.rank.pmax.2)
        
        p2.fisher.cluster = p2
        dec.fisher.cluster = 3
        
      }

      if (p1*p2< (c.fisher)) sign.mat.cluster[1] = 1
      #if (1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2)) <c.inverse) sign.mat[itt, 2] = 1
    }
    
    ## Make interim look decisions for inverse
    if (pmax.inverse<p.max.cut){
      # stop for futility
      dec.inverse.cluster = 1
      p2.inverse.cluster = NA
    } else {
      if (pn.inverse>pn.cut){
        # stop recruiting
        ## observed events in stage 2.
        data.pn.2 = data.stage[(data.stage$stage.1.ind==1)&(data.stage$obs_1=="FALSE"),]
        
        log.rank.pn.2 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pn.2)
        p2 = log.rank.p.func(log.rank.pn.2)
        
        dec.inverse.cluster = 2
        p2 = max(0.0000001, min(p2, 0.9999999))
        p2.inverse.cluster = p2
        
      } else {
        # continue recruiting
        ## observed events in stage 2.
        data.pmax.2 = data.stage[(data.stage$stage.1.ind==0)|(data.stage$obs_1=="FALSE"),]
        
        log.rank.pmax.2 = survdiff(Surv(TTE_1, event_1) ~ group, data = data.pmax.2)
        p2 = log.rank.p.func(log.rank.pmax.2)
        
        dec.inverse.cluster = 3
        p2 = max(0.0000001, min(p2, 0.9999999))
        p2.inverse.cluster = p2
        
      }
      
      p1 = max(0.0000001, min(p1, 0.9999999))

      
      #if (p1*p2< c.fisher) sign.mat[itt, 1] = 1
      if ((1-pnorm(w1.inv*qnorm(1-p1)+w2.inv*qnorm(1-p2))) <(c.inverse)) sign.mat.cluster[2] = 1
    }
    
    clulist = list(dec.gold.cluster = dec.gold.cluster, 
                   dec.fisher.cluster = dec.fisher.cluster, 
                   dec.inverse.cluster = dec.inverse.cluster, 
                   sign.mat.cluster = sign.mat.cluster, 
                   p1.fisher.cluster = p1.fisher.cluster, 
                   p1.inverse.cluster = p1.inverse.cluster,
                   p2.fisher.cluster = p2.fisher.cluster, 
                   p2.inverse.cluster = p2.inverse.cluster,
                   p.gold.cluster = p.gold.cluster)
    return(clulist)
    
  }
  stopCluster(cl)
  
  
  for (itt in 1:n.itt){
    
    pred.temp = pred[[itt]]
    dec.gold.mat[itt] = pred.temp$dec.gold.cluster
    dec.fisher.mat[itt] = pred.temp$dec.fisher.cluster 
    dec.inverse.mat[itt] = pred.temp$dec.inverse.cluster
    sign.mat[itt,] = pred.temp$sign.mat.cluster
    p1.fisher.vec[itt] = pred.temp$p1.fisher.cluster
    p1.inverse.vec[itt] = pred.temp$p1.inverse.cluster
    p2.fisher.vec[itt] = pred.temp$p2.fisher.cluster
    p2.inverse.vec[itt] = pred.temp$p2.inverse.cluster
    p.gold.vec[itt] = pred.temp$p.gold.cluster

  }
  
  print("Done: simulation iterations")
  
  ## uniform (0, 1) p1 and p2
  p1.unif = runif(10^6, 0, 1)
  p2.unif = runif(10^6, 0, 1)
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

for (scen.ind in 1:8){

  #random.seed = NA
   random.seed = 1
   set.seed(random.seed)

  #print(scen.ind)
  n.itt.in = 10^5
  # n.itt.in = 10

  if (scen.ind<=4){
    pn.cut.in = 0.9
    p.max.cut.in = 0.1
    nc.1.in = 100

    enrol.rate.in = 5 # per month

    # if (scen.ind==0) haz.in = -log(0.1)/12
    if (scen.ind==1) haz.in = -log(0.2)/12
    if (scen.ind==2) haz.in = -log(0.3)/12
    if (scen.ind==3) haz.in = -log(0.5)/12
    if (scen.ind==4) haz.in = -log(0.6)/12
    # if (scen.ind==5) haz.in = -log(0.9)/12
  } else if (scen.ind<=8){
    pn.cut.in = 0.9
    p.max.cut.in = 0.1
    nc.1.in = 100

    enrol.rate.in = 2.5 # per month

    # if (scen.ind==0) haz.in = -log(0.1)/12
    if (scen.ind==5) haz.in = -log(0.2)/12
    if (scen.ind==6) haz.in = -log(0.3)/12
    if (scen.ind==7) haz.in = -log(0.5)/12
    if (scen.ind==8) haz.in = -log(0.6)/12
    # if (scen.ind==5) haz.in = -log(0.9)/12
  }

  input = list(alpha = 0.025, # chi square test is for two-sided test
               prior.c.a = 0.3,
               prior.c.b = 1,
               prior.t.a = 0.3,
               prior.t.b = 1,
               pn.cut = pn.cut.in,
               p.max.cut = p.max.cut.in,
               nc.1 = nc.1.in,
               enrol.rate.c = enrol.rate.in,
               nc.max = 200,
               nt.1 = nc.1.in,
               enrol.rate.t = enrol.rate.in,
               nt.max = 200,
               haz.c = haz.in,
               haz.t = haz.in,
               n.itt = n.itt.in)

  power.fit = power.cal(input)

  output.table.temp =  c(apply(power.fit$sign.mat, 2, mean), power.fit$dec.unif,
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

  output.table.temp = t(as.matrix(output.table.temp))

  output.table = rbind(output.table, output.table.temp)

  colnames(output.table.temp) = c("fisher", "inverse", "goldilocks", "fisher unif", "inverse unif",
                             "gold stop futility", "gold stop recruiting", "gold continue",
                             "fisher stop futility", "fisher stop recruiting", "fisher continue",
                             "inverse stop futility", "inverse stop recruiting", "inverse continue",
                             "gold a0", "fisher a0", "inverse a0",
                             "fisher min p1 coh1", "fisher max p1 coh2",
                             "inverse min p1 coh1", "inverse max p1 coh2",
                             "random seed", "n.itt")

  rownames(output.table)[scen.ind] =
    paste("haz", haz.in, "Sn", pn.cut.in, "Fn", p.max.cut.in,
                                           "n1", nc.1.in, "enrol.rate",enrol.rate.in)

  write.csv(output.table.temp, paste0("results/case/prime_type_1_error_",
                                 "haz_", haz.in, "_Sn_", pn.cut.in, "_Fn_",
                                 p.max.cut.in,
                                 "_n1_", nc.1.in, "_enrol.rate_",enrol.rate.in,
                                 "_n.itt_", n.itt.in,
                                 ".csv"))

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
write.csv(output.table, paste0("results/case/prime_type_1_error_com_n.itt", n.itt.in, ".csv"))
#write.csv(output.table, paste0("results/subset/type_1_error_", scen.ind, ".csv"))
# 
# ## type 1 error simulation ends
###############################################################################################

## power simulations
time.temp = Sys.time()
output.table = NULL

for (scen.ind in 1:6){

  random.seed = 1
  set.seed(random.seed)

  print(scen.ind)
  n.itt.in = 10^5

  if (scen.ind<=3){

     pn.cut.in = 0.9
     p.max.cut.in = 0.1
     nc.1.in = 100

     enrol.rate.in = 5 # per month

     haz.c.in = -log(0.3)/12

    if (scen.ind==1){prob.delta = 0.12 }
    if (scen.ind==2){prob.delta = 0.15 }
    if (scen.ind==3){prob.delta = 0.18 }

  } else if (scen.ind<=6){

    pn.cut.in = 0.9
    p.max.cut.in = 0.1
    nc.1.in = 100

    enrol.rate.in = 2.5 # per month

    haz.c.in = -log(0.3)/12

    if (scen.ind==4){prob.delta = 0.11 }
    if (scen.ind==5){prob.delta = 0.14 }
    if (scen.ind==6){prob.delta = 0.17 }

  } else if (scen.ind<=9){

    pn.cut.in = 0.9
    p.max.cut.in = 0.1
    nc.1.in = 100

    enrol.rate.in = 2.5 # per month

    haz.c.in = -log(0.3)/12

    if (scen.ind==7){prob.delta = 0.15 }
    if (scen.ind==8){prob.delta = 0.2 }
    if (scen.ind==9){prob.delta = 0.25 }

  }

  haz.t.in = -log(0.3+prob.delta)/12

     input = list(alpha = 0.025, # chi square test is for two-sided test
                  prior.c.a = 0.3,
                  prior.c.b = 1,
                  prior.t.a = 0.3,
                  prior.t.b = 1,
                  pn.cut = pn.cut.in,
                  p.max.cut = p.max.cut.in,
                  nc.1 = nc.1.in,
                  enrol.rate.c = enrol.rate.in,
                  nc.max = 200,
                  nt.1 = nc.1.in,
                  enrol.rate.t = enrol.rate.in,
                  nt.max = 200,
                  haz.c = haz.c.in,
                  haz.t = haz.t.in,
                  n.itt = n.itt.in)

  power.fit = power.cal(input)

    output.table.temp =  c(apply(power.fit$sign.mat, 2, mean), power.fit$dec.unif,
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

    output.table.temp = t(as.matrix(output.table.temp))

    output.table = rbind(output.table, output.table.temp)

    colnames(output.table.temp) = c("fisher", "inverse", "goldilocks", "fisher unif", "inverse unif",
                               "gold stop futility", "gold stop recruiting", "gold continue",
                               "fisher stop futility", "fisher stop recruiting", "fisher continue",
                               "inverse stop futility", "inverse stop recruiting", "inverse continue",
                               "gold a0", "fisher a0", "inverse a0",
                               "fisher min p1 coh1", "fisher max p1 coh2",
                               "inverse min p1 coh1", "inverse max p1 coh2",
                               "random seed", "n.itt")

    rownames(output.table)[scen.ind] =
      paste("haz", haz.c.in, "_prob.delta_", prob.delta, "Sn", pn.cut.in,
            "Fn", p.max.cut.in,
                                             "n1", nc.1.in, "enrol.rate",enrol.rate.in)

    write.csv(output.table.temp, paste0("results/case/prime_power_",
                                   "haz_", haz.c.in, "_prob.delta_", prob.delta,
                                   "_Sn_", pn.cut.in, "_Fn_",
                                   p.max.cut.in,
                                   "_n1_", nc.1.in, "_enrol.rate_",enrol.rate.in,
                                   ".csv"))


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
write.csv(output.table, "results/case/prime_output_power_com.csv")

## power simulation ends
################################################################################







































