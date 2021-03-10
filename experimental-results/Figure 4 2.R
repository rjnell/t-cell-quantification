

plot_data = function(x,y,yl,yh,file,main) {
  file = paste0(file)
  png(file, res=600, 2500, 3250)  
  par(mar=c(5,5,5,5))
  ylim = c(-50,100)
  xlim = c(0,100)
  plot(x = xlim,
       y = ylim,
       pch = 16,
       xlab = "",
       ylab = "", 
       ylim = ylim,
       xlim = xlim, 
       bty = "l", 
       type="n",
       xaxs = "i", 
       yaxs = "i",
       axes = F)
  
  yat = seq(from=ylim[1],to=ylim[2],by=25)
  xat = seq(from=xlim[1],to=xlim[2],by=25)
  segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
  segments(xat, ylim[1], xat, ylim[2],col="#eeeeee",lwd=1.4,xpd=T)
  
  
  
  segments(xlim[1],ylim[1],xlim[1],ylim[2],xpd=T,col="#B1B1B1",lwd=1.4)
  axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.035,lwd=1.4)
  axis(side = 2,at=yat,las=2,labels=paste0(yat,"%"), lwd=0, col.axis="#333333",line=-0.23)
  
  axis(side = 1,at=xat,las=1,labels=rep("",length(xat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.035,lwd=1.4)
  axis(side = 1,at=xat,las=1,labels=paste0(xat,"%"), lwd=0, col.axis="#333333",line=-0.23, )
  
  mtext(text = "T-cell fraction (mixture)", side = 1, line=2.7,col="#333333")  
  mtext(text = "T-cell fraction (digital PCR)", side = 2, line=3.7,col="#333333")  
  
  regression = lm(y ~ x)
  newx=sort(x)
  pred_interval <- predict(regression, newdata=data.frame(newx), interval="prediction",
                           level = 0.95)
  pred_interval = pred_interval[order(pred_interval[,1]),]
  polygon(c(newx,rev(newx)), c(pred_interval[,2],rev(pred_interval[,3])), col="#eeeeee", border=NA)
  segments(0,0,100,100,col="#b1b1b1",lwd=1.4,xpd=T,lty=3) 
  segments(xlim[1],ylim[1],xlim[2],ylim[1],xpd=T,col="#B1B1B1",lwd=1.4)
  
  for (i in 1:length(x)) {
    arrows(x[i], 
           yl[i], 
           x[i], 
           yh[i], length=0.05, angle=90, code=3, col="#b1b1b1", lwd=1.4, xpd=T)
  }
  points(x = x,
         y = y,
         pch = 15, cex=0.9, col="#333333",xpd=T)
  
 
  
  #segments(min(pbmc_data$`TCF flow cytometry`)-10,
  #         regression$coefficients[1]+regression$coefficients[2]*(min(pbmc_data$`TCF flow cytometry`)-10),
  #         max(pbmc_data$`TCF flow cytometry`)+10,
  #         regression$coefficients[1]+regression$coefficients[2]*(max(pbmc_data$`TCF flow cytometry`)+10))
  
  r = cor.test(x, y)
  
  slope = format(round(regression$coefficients[2]*100)/100, nsmall = 2)
  t = paste0("y = ", slope, "x + ", format(round(regression$coefficients[1]*0.01, digits = 2), nsmall = 2))
  if (regression$coefficients[1] < 0) {
    t = paste0("y = ", slope, "x - ", format(round(regression$coefficients[1]*-0.01, digits = 2), nsmall = 2)) 
  }
  
  
  #lines(newx, pred_interval[,3], col="orange", lty=2)
  
  #t = bquote(italic(y) == .(format(r$estimate^2, digits = 3)) - 1)
  rect(5,97.5,60,87.5,border=NA,col="#FFFFFF")
  text(0, 92.5, labels = t,
       pos=4,col="#333333")
  
  #t = substitute("R"^2," = ", ) 
  v = 
    #t = expression(paste("R"^"2","="),)
    #t = paste0("*km^2~size"~round(r$estimate, digits = 2))
    #t = paste0(expression(R^2),(round(r$estimate, digits = 2)))
    
    t = bquote("R"^2~"="~.(format(r$estimate^2, digits = 3)))
  
  #rect(5,87.5,40,77.5,border=NA,col="#FFFFFF")
  #text(0, 82.5, labels = t,
  #     pos=4,col="#333333")
  
  mtext(text = main, side = 3, line=1.2,col="#333333", font=4)
  
  dev.off()
  system(paste("open",file))
}




mix_data = read_xlsx("D:/SURFDRIVE/ORGANIZED/Projecten/Project 'Multiplexing'/manuscript/data/Book1.xlsx", sheet=5)
plot_data(x=mix_data$`TCF mixed`,
          y=mix_data$`Classic TCF (estimate)`,
          yl=mix_data$`Classic TCF (95%-CI min)`,
          yh=mix_data$`Classic TCF (95%-CI max)`,
          "experimental-results/Figure-C.png", main="Classic model")
plot_data(x=mix_data$`TCF mixed`,
          y=mix_data$`Adjusted TCF (estimate)`,
          yl=mix_data$`Adjusted TCF (95%-CI min)`,
          yh=mix_data$`Adjusted TCF (95%-CI max)`,
          "experimental-results/Figure-D.png", main="Adjusted model")


mean(mix_data$`Adjusted TCF (estimate)` - mix_data$`TCF mixed` )
