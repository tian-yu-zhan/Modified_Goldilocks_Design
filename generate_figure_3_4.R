
## plot simulation scenario plot
library(gridExtra)
library(ggplot2)
library(grid)
library(reshape)
library(xtable)
library("ggplot2")
library("reshape")
library("scales")
# library("wesanderson")

setwd("~/goldilocks/")

output.table = read.csv(paste0("results/subset/same_type_1_error_multiple_plot_prime", ".csv"))
output.table = data.frame(output.table)


n.cutoff = 5
n.scen = dim(output.table)[1]/2/n.cutoff
data.in = data.frame(power = c(output.table$fisher[1:n.scen], output.table$inverse[1:n.scen], 
                               output.table$goldilocks[1:(n.scen*n.cutoff)],
                               output.table$fisher[(1:n.scen)+n.scen*n.cutoff],
                               output.table$inverse[(1:n.scen)+n.scen*n.cutoff],
                               output.table$goldilocks[(1:(n.scen*n.cutoff))+n.scen*n.cutoff]),
                      group = c(rep("MGDF Type I error", n.scen), rep("MGDI Type I error", n.scen),
                                rep(paste("GD Type I error with cutoff", 0.025), n.scen), 
                                rep(paste("GD Type I error with cutoff", 0.022), n.scen),
                                rep(paste("GD Type I error with cutoff", 0.02), n.scen),
                                rep(paste("GD Type I error with cutoff", 0.015), n.scen),
                                rep(paste("GD Type I error with cutoff", 0.01), n.scen),
                                rep("MGDF power", n.scen), rep("MGDI power", n.scen),
                                rep(paste("GD power with cutoff", 0.025), n.scen), 
                                rep(paste("GD power with cutoff", 0.022), n.scen),
                                rep(paste("GD power with cutoff", 0.02), n.scen),
                                rep(paste("GD power with cutoff", 0.015), n.scen),
                                rep(paste("GD power with cutoff", 0.01), n.scen)
                                ),
                      n1un = rep((1:n.scen)*5, 2*7),
                      power_type = c(rep("Type I Error", n.scen*7), rep("Power", n.scen*7)),
                      method = c(rep(rep(c("MGDF", "MGDI", paste("GD with cutoff", 0.025),
                                         paste("GD with cutoff", 0.022),
                                         paste("GD with cutoff", 0.02),
                                         paste("GD with cutoff", 0.015),
                                         paste("GD with cutoff", 0.01)), each = n.scen), 2))
                      )
data.in$method = factor(as.character(data.in$method), 
                        levels = c("MGDF", "MGDI", "GD with cutoff 0.025", "GD with cutoff 0.022", 
                                   "GD with cutoff 0.02", "GD with cutoff 0.015", "GD with cutoff 0.01"))
  
  
  
  


## plot
data.in.error.wide = data.in[which(data.in$power_type=="Type I Error"), ]

png("results/simulation_scenario_plot_error.png", width = 1800, height = 1200)
ggplot.fit = ggplot(data.in.error.wide, aes(x = n1un, y = power, color = method)) +
  geom_line(aes(alpha = method), size =2) +  #plot flow
  geom_point(size = 8, aes(shape = method, alpha = method))+
  #scale_linetype_manual(values=c("dotted", "solid")) + 
  scale_shape_manual(values=c(17, 17, rep(16, 5)))+
  scale_alpha_manual(values=c(1, 1, rep(0.7, 5)))+
  scale_color_manual(values=
                       c("black", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#E69F00", "#CC79A7"))+
  # scale_y_continuous(breaks = c(0.75, 0.8, 0.85, 0.9), limits = c(0.65, 0.9), sec.axis = sec_axis(~.*0.5-0.65/2, name = "Type I error", breaks = c(0, 0.01, 0.025, 0.04))) + 
  scale_y_continuous(breaks = c(0, 0.013, 0.025, 0.035), limits = c(0, 0.035)) +
  # scale_y_continuous(sec.axis = sec_axis(~., name = "Type I error")) + 
  scale_x_continuous(breaks = c(5, 20, 35, 50, 65, 80, 95))+
  labs(title = "") +
  ylab ("Type I error rate") + xlab("Number of Subjects with Unobserved Responses") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,1,1,1),units="lines"),
        text = element_text(size=45),
        axis.text.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=45,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=40,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text = element_text(colour="black", size = 37, face = "plain"),
        legend.title = element_text(colour="black", size = 42, face = "plain"),
        legend.key.size = unit(2,"line"),
        legend.position="bottom", plot.title = element_text(hjust = 0.5))

print(ggplot.fit)

dev.off()

data.in.power.wide = data.in[which(data.in$power_type=="Power"), ]

png("results/simulation_scenario_plot_power.png", width = 1800, height = 1200)
ggplot.fit = ggplot(data.in.power.wide, aes(x = n1un, y = power, color = method)) +
  geom_line(aes(alpha = method), size =2) +  #plot flow
  geom_point(size = 8, aes(shape = method, alpha = method))+
  #scale_linetype_manual(values=c("dotted", "solid")) + 
  scale_shape_manual(values=c(17, 17, rep(16, 5)))+
  scale_alpha_manual(values=c(1, 1, rep(0.7, 5)))+
  scale_color_manual(values=
                       c("black", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#E69F00", "#CC79A7"))+
  # scale_y_continuous(breaks = c(0.75, 0.8, 0.85, 0.9), limits = c(0.65, 0.9), sec.axis = sec_axis(~.*0.5-0.65/2, name = "Type I error", breaks = c(0, 0.01, 0.025, 0.04))) + 
  scale_y_continuous(breaks = c(0.75, 0.8, 0.85, 0.9), limits = c(0.75, 0.9)) +
  # scale_y_continuous(sec.axis = sec_axis(~., name = "Type I error")) + 
  scale_x_continuous(breaks = c(5, 20, 35, 50, 65, 80, 95))+
  labs(title = "") +
  ylab ("Power") + xlab("Number of Subjects with Unobserved Responses") +
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,1,1,1),units="lines"),
        text = element_text(size=45),
        axis.text.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=45,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=40,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text = element_text(colour="black", size = 37, face = "plain"),
        legend.title = element_text(colour="black", size = 42, face = "plain"),
        legend.key.size = unit(2,"line"),
        legend.position="bottom", plot.title = element_text(hjust = 0.5))

print(ggplot.fit)

dev.off()

# png("manuscript/v6/simulation_scenario_plot.png", width = 2400, height = 1600)
# ggplot.fit = ggplot(data.in.wide, aes(x = n1un, y = power, color = method)) +
#   geom_line(aes(alpha = method, linetype=power_type), size =2) +  #plot flow
#   geom_point(size = 8, aes(shape = method, alpha = method))+
#   geom_line(aes(y = error,  color = method, alpha = method), 
#             size = 2, linetype="solid")+ 
#   geom_point(size = 8, aes(y = error, shape = method, alpha = method))+
#   scale_linetype_manual(values=c("dotted", "solid")) + 
#   scale_shape_manual(values=c(17, 17, rep(16, 5)))+
#   scale_alpha_manual(values=c(1, 1, rep(0.7, 5)))+
#   scale_color_manual(values=
#                        c("black", "#0072B2", "#56B4E9", "#009E73", "#D55E00", "#E69F00", "#CC79A7"))+
#   scale_y_continuous(breaks = c(0.75, 0.8, 0.85, 0.9), limits = c(0.65, 0.9), sec.axis = sec_axis(~.*0.5-0.65/2, name = "Type I error", breaks = c(0, 0.01, 0.025, 0.04))) + 
#   # scale_y_continuous(sec.axis = sec_axis(~., name = "Type I error")) + 
#   scale_x_continuous(breaks = c(5, 20, 35, 50, 65, 80, 95))+
#   labs(title = "") +
#   ylab ("Power") + xlab("Number of Subjects with Unobserved Responses") +
#   theme_bw()+
#   theme(plot.background = element_rect(fill = "transparent"),
#         plot.margin = unit(c(2,0,1,1),units="lines"),
#     text = element_text(size=45),
#         axis.text.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=.5,face="plain"),
#         axis.text.y = element_text(colour="black",size=45,angle=0,hjust=1,vjust=0,face="plain"),
#         axis.title.x = element_text(colour="black",size=45,angle=0,hjust=.5,vjust=0,face="plain"),
#         axis.title.y = element_text(colour="black",size=45,angle=90,hjust=.5,vjust=.5,face="plain"),
#         legend.text = element_text(colour="black", size = 37, face = "plain"),
#         legend.title = element_text(colour="black", size = 42, face = "plain"),
#         legend.key.size = unit(2,"line"),
#         legend.position="right", plot.title = element_text(hjust = 0.5))
#   
# print(ggplot.fit)
# 
# dev.off()

###########################################################
## GD type I error table
# GD.table = data.in[(data.in$power_type=="Type I Error")&((!(data.in$method=="MGDF"))&
#                                                                  (!data.in$method=="MGDI")),]
# GD.table.wide = matrix(GD.table$power, nrow = 5, ncol = 19, byrow = TRUE)
# rownames(GD.table.wide) = unique(GD.table$method)
# colnames(GD.table.wide) = unique(GD.table$n1un)
# 
# GD.table.final = GD.table.wide[1:4, c(4, 7, 10, 13, 16, 19)]
# print(xtable(GD.table.final, digits = 4, include.rownames=FALSE))


