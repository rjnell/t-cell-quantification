data = list()

data[["CM-1"]] = list()
data[["CM-1"]][["Classic TCF"]] = c(-21.48,-26.89,-16.33)
data[["CM-1"]][["Adjusted TCF"]] = c(8.82,3.66,13.99)
data[["CM-1"]][["TCF positions"]] = c(0.25,1.35,2.325)
data[["CM-1"]][["IHC TCFs"]] = c(8.0,6.1,0.0,0.0,0.8)
data[["CM-1"]][["IHC positions"]] = c(0,0, -0.075, 0.075, 0)
data[["CM-1"]][["IHC combined TCF"]] = 3.4    

data[["CM-2"]] = list()
data[["CM-2"]][["Classic TCF"]] = c(-21.83,-27.96,-16.02)
data[["CM-2"]][["Adjusted TCF"]] = c(11.89,6.06,17.75)
data[["CM-2"]][["TCF positions"]] = c(0.25,1.3,2.275)
data[["CM-2"]][["IHC TCFs"]] = c(9.0,7.7,22.2,6.9,5.9)
data[["CM-2"]][["IHC positions"]] = c(0,0.075,0,-0.075, 0)
data[["CM-2"]][["IHC combined TCF"]] = 10.5  

data[["CM-3"]] = list()
data[["CM-3"]][["Classic TCF"]] = c(5.12,1.01,9.06)
data[["CM-3"]][["Adjusted TCF"]] = c(16.38,12.09,20.69)
data[["CM-3"]][["TCF positions"]] = c(0.35,1.3,2.275)
data[["CM-3"]][["IHC TCFs"]] = c(17.3,40.0,24.8,16.7,15.8)
data[["CM-3"]][["IHC positions"]] = c(0.075,0,0,-0.075, 0)
data[["CM-3"]][["IHC combined TCF"]] = 22.9 

for (id in c("CM-1", "CM-2", "CM-3")) {
  
  file = paste0("experimental-results/Figure-5-",id,".png")
  png(file, res=600, 3000, 2500)  
  par(mar=c(5,5,5,5))
  xlim = c(0,3)
  ylim = c(-30,45)
  
  plot(xlim, 
       ylim, 
       type = "n", 
       axes = F, 
       main= "", 
       xlab = "",
       ylab = "",
       xaxs = "i", 
       yaxs = "i")
  
  xat = seq(xlim[1], xlim[2], by=25)
  yat = seq(ylim[1], ylim[2], by=10)
  axis(side = 1, at = xat, labels = rep("",length(xat)), col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  axis(side = 2, at = yat, labels = paste0(yat,"%"), las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  title(ylab="T-cell fraction", line = 4, col="#333333")
  title(main=id, line = 1, col="#333333")
  
  segments(xat, ylim[1], xat, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
  
  segments(1,ylim[1], 1, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
  
  segments(xlim[1], 0, xlim[2], 0, col="#b1b1b1", lwd=1.4, xpd=T)
  segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
  
  xpos = data[[id]][["TCF positions"]][1]
  tcf = data[[id]][["Classic TCF"]]
  arrows(xpos, tcf[2], xpos, tcf[3], code=3, angle=90, length=0.05, col="#b1b1b1", lwd=1.4)
  points(xpos, tcf[1],pch=15)
  text(xpos+0.05, tcf[1], labels = paste0(round(tcf[1]),"%"), pos=4, col="#333333")
  text(0.5, -32.5, labels = "Classic", pos=1, font=3, xpd=T, col="#333333")
  text(0.5, -37.5, labels = "model", pos=1, font=3, xpd=T, col="#333333")
  
  xpos = data[[id]][["TCF positions"]][2]
  tcf = data[[id]][["Adjusted TCF"]]
  arrows(xpos, tcf[2], xpos, tcf[3], code=3, angle=90, length=0.05, col="#b1b1b1", lwd=1.4)
  points(xpos, tcf[1],pch=15)
  text(xpos+0.05, tcf[1], labels = paste0(round(tcf[1]),"%"), pos=4, col="#333333")
  text(1.5, -32.5, labels = "Adjusted", pos=1, font=3, xpd=T, col="#333333")
  text(1.5, -37.5, labels = "model", pos=1, font=3, xpd=T, col="#333333")
  
  segments(2,ylim[1], 2, ylim[2], col="#b1b1b1", lwd=1.4, xpd=T, lty=3)
  
  xpos = data[[id]][["TCF positions"]][3]
  tcfs = data[[id]][["IHC TCFs"]]
  xposs = data[[id]][["IHC positions"]]
  for (i in 1:length(tcfs)) {
    points(xpos+xposs[i], tcfs[i],pch=16,col="#B1B1B1",xpd=T)  
  }
  segments(xpos-0.125, data[[id]][["IHC combined TCF"]], xpos+0.125, col="#333333", lwd=1.4)
  points(xpos, data[[id]][["IHC combined TCF"]],pch=15)
  text(xpos+0.15, data[[id]][["IHC combined TCF"]], labels = paste0(round(data[[id]][["IHC combined TCF"]]),"%"), pos=4, col="#333333")
  text(2.5, -32.5, labels = "IHC", pos=1, font=1, xpd=T, col="#333333")
  text(2.5, -37.5, labels = "(CD3+)", pos=1, font=1, xpd=T, col="#333333")
  
  dev.off()
  system(paste("open",file))
}