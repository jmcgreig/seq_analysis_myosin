library(MASS)
library(foreign)
library(quantreg)  
library(MASS)
library(foreign)
library(quantreg)  
library(robustbase)

setwd('/home/jake/Documents/projects/myosin/chimera/scripts/')

# pdf("id_vs_mass.pdf", paper="a4", width = 8.3, height = 11.7)
# jpeg(filename = "id_vs_mass_paper_sup.jpeg",
     # width = 18, height = 25, units = "cm",
     # quality = 75, res = 200)
files <- c("beta_mot.csv")
kount = 0


opar <- par(mfrow = c (3,2), mar = c(0,1,1.5,1), oma = c(4, 4, 0, 0))

for (input in files){
  if (kount == 0){
    
    jpeg(filename = "../figures/id_vs_mass_2c_reduced.jpeg",
         width = 18, height = 25, units = "cm",
         quality = 75, res = 200)
    opar <- par(mfrow = c (3,2), mar = c(0,1,1.5,1), oma = c(4, 4, 0, 0))
    
  }
  print(kount)
  kount = kount +1
  if (kount == 7){
    mtext('Motor', side = 1, outer = TRUE, line = -70, adj = 0.25)
    mtext('Tail', side = 1, outer = TRUE, line = -70, adj = 0.75)
    mtext('Mass (kg)', side = 1, outer = TRUE, line = 2 )
    mtext('Sequence Identity (%)', side = 2, outer = TRUE, line = 2)
    dev.off()
    jpeg(filename = "../figures/id_vs_mass_2c_reduced.jpeg",
         width = 18, height = 25, units = "cm",
         quality = 75, res = 200)
    opar <- par(mfrow = c (3,2), mar = c(0,1,1.5,1), oma = c(4, 4, 0, 0))
    
  }
  ni <- gsub("_.*?csv", "", input)
  rtab = read.csv(input, header = T)
  mass <- rtab$Mass..kg.
  id <- rtab$id

  
  ########Robust regression################################
  slopey_lin <- lm(id ~ log10(mass))$coeff[2]
  slopey_rr <- rlm(id ~ log10(mass), data=rtab)$coeff[2]
  slopey_rr
  fm.orig <- lm(id ~ log10(mass), data=rtab)
  # fm.rq <- rq(id ~ log10(mass), data=rtab)
  fm.rlm <- rlm(id ~ log10(mass), data=rtab)
  
  summ <- summary(fm.orig)
  Rvalue <- summ$r.squared
  summRR <- summary(fm.rlm)
  sloptRR <- summRR$coefficients[[2]]
  errorRR <- summRR$coefficients[[4]]

  
  summ <- summary(fm.orig)
  slopt <- summ$coefficients[[2]]
  
  error <- summ$coefficients[[4]]

  model = lmrob(id ~ log10(mass))
  g <- summary(model)
  robslop <- g$coefficients[2]
  roberr <- g$coefficients[4]
  g$coefficients[4]
  robR <- g$r.squared
  
  rtab$Clade
  rtab$Clade <- as.character(rtab$Clade)
  rtab$Clade[rtab$Clade == "Euarchontoglires"] <- 1
  rtab$Clade[rtab$Clade == "Laurasiatheria"] <- 2
  rtab$Clade[rtab$Clade == "zAfrotheria"] <- 3
  rtab$Clade[rtab$Clade == "Metatheria"] <- 3
  as.factor(rtab$Clade)
  
  

  
  
  plot(id ~ mass, data=rtab,
       type = "p",
       log = "x",
       xlab = "",
       ylab = "",
       col = c("grey36", "black", "black", "grey36")[as.factor(rtab$Clade)],
       pch = c(15, 17, 1)[as.factor(rtab$Clade)],
       bty = "o",
       cex = 1.5,
       ylim = c(90, 100),
       xlim=c(0.001,100000),
       
       xaxt = "n",
       yaxt = "n")
  
  if (grepl('peri?emb?2a?',  input)){
    if (grepl('mot',  input)){
      title (main="Motor")}
    if (grepl('tail',  input)){
      title (main="Tail")}}
  
  # abline(fm.orig, lty="dashed")    # use a dashed line
  #abline(fm.rq, lty = "dotted")
  abline(model, lty="solid")
  
  
  # legend("bottomleft", inset=c(0.025, 0.18), bty="n", 
  #        legend = c("Linear Model Fit ",  "Robust Linear Model Fit" ),
  #        lty = c(2,   1),      # 1 = "solid" ; 2 = "dashed"
  #        col = c("black", "black"),
  #        text.font = 2, cex = 0.8)
  # 
  # 
  # legend("bottomleft", inset=c(-0.03, 0.1), 
  #        legend = paste(c("  R² =", round(Rvalue, digits = 3),
  #                         "        Robust-R² = ", round(robR, digits = 3),
  #                         "\n - - Linear Slope =", round(slopt, digits = 3), 
  #                         "+/- " , round(error, digits = 3),
  #                         "\n — Robust Linear Slope =", round(robslop, digits = 3),
  #                         "+/- " , round(roberr, digits = 3)), collapse = " ") , 
  #        text.font = 2, bty = "n", cex = 0.8)
  # 
  legend("bottomleft", inset=c(-0.03, 0.03), 
         legend = paste(c("Robust-R² = ", round(robR, digits = 3), 
                          "\nRobust Linear Slope =", round(robslop, digits = 3),
                          "+/- " , round(roberr, digits = 3)), collapse = " ") , 
         text.font = 2, bty = "n", cex = 1.1)
  
  print(ni)
  
  
  

  axis(1, at=c(0.01,0.1,1,10,100,1000,10000,100000), labels= FALSE,   tck = 0.02) 
  
  if (grepl('mot',  input)){
    axis(2,  at=c(90, 92, 94, 96, 98, 100), labels= c(90, 92, 94, 96, 98, 100))
    # axis(2,  at=c(80, 84, 88, 92, 96, 100), labels= c(80, 84, 88, 92, 96, 100))
    
  }
  if (grepl('tail',  input)){
    axis(2,  at=c(90, 92, 94, 96, 98, 100), labels= FALSE)
    # axis(2,  at=c(80, 84, 88, 92, 96, 100), labels= FALSE)
    
  }
  
  if (grepl('mot',  input)){
    if (!(grepl('alpha',  input))){
    legend("bottomright", legend=toupper(ni), text.font = 2, bty = "n", cex = 2, inset = c(0, 0.04)  )
    }
    if ((grepl('alpha',  input))){
      legend("bottomright", legend=expression(alpha), text.font = 2, bty = "n", cex = 3)
    }
  }
  
    if (grepl('sm',  input)){
      axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
    }
  if (grepl('nmb',  input)){
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
  if (grepl('alpha',  input)){
    axis(1,  at=c(0.01,0.1,1,10,100,1000,10000, 100000), labels= c("0.01", "0.1", "1", "10", "100",  expression(paste("10"^"3")), expression(paste("10"^"4")), expression(paste("10"^"5"))), cex.axis = 1.2)
  }
}


mtext('Motor', side = 1, outer = TRUE, line = -70, adj = 0.25)
mtext('Tail', side = 1, outer = TRUE, line = -70, adj = 0.75)
mtext('Mass (kg)', side = 1, outer = TRUE, line = 2 )
mtext('Sequence Identity (%)', side = 2, outer = TRUE, line = 2)
dev.off()