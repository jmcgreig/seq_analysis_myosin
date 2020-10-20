library(MASS)
library(foreign)
library(quantreg)  
library(robustbase)


setwd('/home/jake/Documents/myosin_size_var/dnds/plots/')
# pdf("dndsplots_sup.pdf", paper="a4", width = 8.3, height = 11.7)
jpeg(filename = "dnds_sup.jpeg",
     width = 18, height = 25, units = "cm",
     quality = 75, res = 200)

# pdf("dndsplots.pdf", paper="a4", width = 8.3, height = 11.7)

files <- c("2b_mot.csv", "2b_tail.csv",
           "2a_mot.csv", "2a_tail.csv",
           "2x_mot.csv", "2x_tail.csv",
           "alpha_motor.csv", "alpha_tail.csv")

# 
# files <- c("emb_mot.csv", "emb_tail.csv",
#            "nma_mot.csv", "nma_tail.csv")


opar <- par(mfrow = c (4,2), mar = c(0,1,1.5,1), oma = c(3, 3, 0, 0))

for (input in files){
  input
  ni <- gsub("_.*?csv", "", input)
  rtab = read.csv(input, header = T)
  rtab
  mass <- rtab$Mass..kg.
  mass
  cla <- rtab$Clade
  
  rtab$Clade <- as.character(rtab$Clade)
  rtab$Clade[rtab$Clade == "Euarchontoglires"] <- 1
  rtab$Clade[rtab$Clade == "Laurasiatheria"] <- 2
  rtab$Clade[rtab$Clade == "zAfrotheria"] <- 3
  rtab$Clade[rtab$Clade == "Metatheria"] <- 3
  as.factor(rtab$Clade)
  
  rtab$Clade
  ########Robust regression################################
  model1 = lmrob(rtab$dn.ds ~ log10(mass))
  plot(rtab$dn.ds ~ mass, data=rtab,
       type = "p",
       log = "x",
       xlab = "",
       ylab = "",
       col = c("grey36", "black", "black", "grey36")[as.factor(rtab$Clade)],
       pch = c(15, 17, 1)[as.factor(rtab$Clade)],
       bty = "o",
       cex = 1.5,
       ylim = c(0.0,0.06),
       xlim=c(0.001,100000),
       
       xaxt = "n",
       yaxt = "n")
  g <- summary(model1)
  robslop <- g$coefficients[2]
  roberr <- g$coefficients[4]
  g$coefficients[4]
  robR <- g$r.squared
  robR
  
  

    
  axis(1,  at=c(0.01,0.1,1,10,100,1000,10000,100000), labels= FALSE)
  if (grepl('emb', input) | grepl("nma", input)){
    legend("topleft", inset=c(0, 0.02), 
           legend = paste(c("R² = ", round(robR, digits = 3),
                            "\nSlope =", round(robslop, digits = 3),
                            "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  else{
  legend("bottomright", inset=c(0, 0.02), 
  legend = paste(c("R² = ", round(robR, digits = 3),
                   "\nSlope =", round(robslop, digits = 3),
                   "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  
  
    if (grepl('mot', input)){
      if (!(grepl('alpha', input))){
        
      # title (main="Motor")
      legend("topright", legend=toupper(ni), text.font = 2, bty = "n", cex = 2)
      }
      # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  0,0.2,0.4,0.6,0.8))
      axis(2, at= seq(0,0.08, by = 0.02), labels = seq(0,0.08, by = 0.02), cex.axis = 1.2)
      axis(2, at= seq(0,0.8, by = 0.01), labels = F)
      
    }
    if (grepl('tail', input)){
      # title (main="Tail")
      # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  "","","","",""))
      # axis(2,  at=c( 0,0.02,0.04,0.06,0.08), labels=F)
      axis(2, at= seq(0,0.08, by = 0.01), labels = F)
      
    }
    if (grepl('alpha', input)){
      if (grepl('mot', input)){
        
      legend("topright", legend=expression(alpha), text.font = 2, bty = "n", cex = 3)
      }
      axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
    }
  if (grepl('nma', input)){
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }

    
    all = lmrob(rtab$dn.ds ~ log10(mass))
    abline(all, lty="solid", col = 'black')
  }
  
  
  
  
mtext('Motor', side = 1, outer = TRUE, line = -71, adj = 0.25)
mtext('Tail', side = 1, outer = TRUE, line = -71, adj = 0.75)
  
  mtext('Mass (kg)', side = 1, outer = TRUE, line = 2 )
  mtext('dN/dS', side = 2, outer = TRUE, line = 1)
  

dev.off()

# pdf("dnplots_sup.pdf", paper="a4", width = 8.3, height = 11.7)
# pdf("dnplots.pdf", paper="a4", width = 8.3, height = 11.7)
jpeg(filename = "dn_sup2.jpeg",
     width = 18, height = 25, units = "cm",
     quality = 75, res = 200)



opar <- par(mfrow = c (4,2), mar = c(0,1,1.5,1), oma = c(3, 3, 0, 0))

for (input in files){
  input
  ni <- gsub("_.*?csv", "", input)
  rtab = read.csv(input, header = T)
  rtab
  mass <- rtab$Mass..kg.
  mass
  cla <- rtab$Clade
  rtab$Clade <- as.character(rtab$Clade)
  rtab$Clade[rtab$Clade == "Euarchontoglires"] <- 1
  rtab$Clade[rtab$Clade == "Laurasiatheria"] <- 2
  rtab$Clade[rtab$Clade == "zAfrotheria"] <- 3
  rtab$Clade[rtab$Clade == "Metatheria"] <- 3
  as.factor(rtab$Clade)
  ########Robust regression################################
  model1 = lmrob(rtab$dn ~ log10(mass))
  plot(rtab$dn ~ mass, data=rtab,
       type = "p",
       log = "x",
       xlab = "",
       ylab = "",
       col = c("grey36", "black", "black", "grey36")[as.factor(rtab$Clade)],
       pch = c(15, 17, 1)[as.factor(rtab$Clade)],
       bty = "o",
       cex = 1.5,
       ylim = c(0.0,0.04),
       xlim=c(0.001,100000),
       
       xaxt = "n",
       yaxt = "n")
  g <- summary(model1)
  robslop <- g$coefficients[2]
  roberr <- g$coefficients[4]
  g$coefficients[4]
  robR <- g$r.squared
  robR
  
  
  
  
  axis(1,  at=c(0.01,0.1,1,10,100,1000,10000,100000), labels= FALSE)
  
  if (grepl('emb', input) | grepl("nma", input)){
    legend("topleft", inset=c(0, 0.02), 
           legend = paste(c("R² = ", round(robR, digits = 3),
                            "\nSlope =", round(robslop, digits = 3),
                            "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  else{
    legend("bottomright", inset=c(0, 0.02), 
           legend = paste(c("R² = ", round(robR, digits = 3),
                            "\nSlope =", round(robslop, digits = 3),
                            "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  
  if (grepl('mot', input)){
    if (!(grepl('alpha', input))){
    # title (main="Motor")
    legend("topright", legend=toupper(ni), text.font = 2, bty = "n", cex = 2)}
    # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  0,0.2,0.4,0.6,0.8))
    axis(2, at= seq(0,0.08, by = 0.02), labels = seq(0,0.08, by = 0.02), cex.axis = 1.2)
    axis(2, at= seq(0,0.8, by = 0.01), labels = F)
    
  }
  if (grepl('tail', input)){
    # title (main="Tail")
    # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  "","","","",""))
    # axis(2,  at=c( 0,0.02,0.04,0.06,0.08), labels=F)
    axis(2, at= seq(0,0.08, by = 0.01), labels = F)
  }
  if (grepl('alpha', input)){
    if (grepl('mot', input)){
      
    
    legend("topright", legend=expression(alpha), text.font = 2, bty = "n", cex = 3)
    }
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
  
  if (grepl('nma', input)){
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
  
  all = lmrob(rtab$dn ~ log10(mass))
  abline(all, lty="solid", col = 'black')
}





mtext('Motor', side = 1, outer = TRUE, line = -71, adj = 0.25)
mtext('Tail', side = 1, outer = TRUE, line = -71, adj = 0.75)

mtext('Mass (kg)', side = 1, outer = TRUE, line = 2 )
mtext('dN', side = 2, outer = TRUE, line = 1)

dev.off()

# pdf("dsplots_sup.pdf", paper="a4", width = 8.3, height = 11.7)
# pdf("dsplots.pdf", paper="a4", width = 8.3, height = 11.7)
jpeg(filename = "ds_sup.jpeg",
     width = 18, height = 25, units = "cm",
     quality = 75, res = 200)

opar <- par(mfrow = c (4,2), mar = c(0,1,1.5,1), oma = c(3, 3, 0, 0))

for (input in files){
  input
  ni <- gsub("_.*?csv", "", input)
  rtab = read.csv(input, header = T)
  rtab
  mass <- rtab$Mass..kg.
  mass
  cla <- rtab$Clade
  rtab$Clade <- as.character(rtab$Clade)
  rtab$Clade[rtab$Clade == "Euarchontoglires"] <- 1
  rtab$Clade[rtab$Clade == "Laurasiatheria"] <- 2
  rtab$Clade[rtab$Clade == "zAfrotheria"] <- 3
  rtab$Clade[rtab$Clade == "Metatheria"] <- 3
  as.factor(rtab$Clade)
  ########Robust regression################################
  model1 = lmrob(rtab$ds ~ log10(mass))
  plot(rtab$ds ~ mass, data=rtab,
       type = "p",
       log = "x",
       xlab = "",
       ylab = "",
       col = c("grey36", "black", "black", "grey36")[as.factor(rtab$Clade)],
       pch = c(15, 17, 1)[as.factor(rtab$Clade)],
       bty = "o",
       cex = 1.5,
       ylim = c(0.3,1.15),
       xlim=c(0.001,100000),
       
       xaxt = "n",
       yaxt = "n")
  g <- summary(model1)
  robslop <- g$coefficients[2]
  roberr <- g$coefficients[4]
  g$coefficients[4]
  robR <- g$r.squared
  robR
  
  
  
  
  axis(1,  at=c(0.01,0.1,1,10,100,1000,10000,100000), labels= FALSE)
  
  if (grepl('emb', input) | grepl("nma", input)){
    legend("topleft", inset=c(0, 0.02), 
           legend = paste(c("R² = ", round(robR, digits = 3),
                            "\nSlope =", round(robslop, digits = 3),
                            "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  else{
    legend("bottomright", inset=c(0, 0.02), 
           legend = paste(c("R² = ", round(robR, digits = 3),
                            "\nSlope =", round(robslop, digits = 3),
                            "± " , round(roberr, digits = 3)), collapse = " ") , 
           text.font = 1, bty = "n", cex = 1.2)
  }
  
  if (grepl('mot', input)){
    if (!(grepl('alpha', input))){
    # title (main="Motor")
    legend("topright", legend=toupper(ni), text.font = 2, bty = "n", cex = 2)
    # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  0,0.2,0.4,0.6,0.8))
    }
    axis(2, at= seq(0,2, by = 0.2), labels = seq(0,2, by = 0.2), cex.axis = 1.2)
    axis(2, at= seq(0,2, by = 0.1), labels = F)
    
  }
  if (grepl('tail', input)){
    # title (main="Tail")
    # axis(2,  at=c( 0,0.2,0.4,0.6,0.8), labels= c(  "","","","",""))
    # axis(2,  at=c( 0,0.02,0.04,0.06,0.08), labels=F)
    axis(2, at= seq(0,2, by = 0.1), labels = F)
  }
  if (grepl('alpha', input)){
    if (grepl('mot', input)){
    legend("topright", legend=expression(alpha), text.font = 2, bty = "n", cex = 3, inset = c(0, -0.1))
    }
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
  if (grepl('nma', input)){
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
  
  
  all = lmrob(rtab$ds ~ log10(mass))
  abline(all, lty="solid", col = 'black')
}



mtext('Motor', side = 1, outer = TRUE, line = -71, adj = 0.25)
mtext('Tail', side = 1, outer = TRUE, line = -71, adj = 0.75)


mtext('Mass (kg)', side = 1, outer = TRUE, line = 2 )
mtext('dS', side = 2, outer = TRUE, line = 1)


dev.off()







