    library(Hmisc) 
    library(ncar)
    library(plotrix)
            #-----------------------------------------------
      # Store graphics in a pdf file. Comment out this
      # line incar v0.3.4f you want plots to the screen
      #-----------------------------------------------
    pdf("../figures/p_value_plots_log10.pdf", paper="a4", width = 8.3, height = 11.7)
    # zz <- file("/home/jake/Documents/file.txt","w")
    
    

      #---------------------------------------------------
      # Threshold for running analysis. Analyses are only
      # run if there are at least 2 different amino acids.
      # (Only the first 2 are kept, others are replaced by
      # missing values).
      # Also, the analysis is only run if 
      # freq_2 / freq_1 > threshold
      # where freq_1 and freq_2 are the frequencies of the 
      # commonest and second commonest alleles
      #---------------------------------------------------
    threshold = 0.1
    
      #-------------------------------------------
      # Initial processing of data - read the data 
      # and extract the mass and clade columns
      #-------------------------------------------
    dd <- read.csv("../data/myh7_plot_H_.csv",header=TRUE, 
                    colClasses = c("factor", "factor", "factor", "numeric",
                                 rep("factor", 808)))
    mass = dd$Mass
    clade = dd$Clade
    
      #---------------------------------------------------------------
      # Only consider clades Euarchontoglires and Laurasiatheria
      # Replace clades Afrotheria (4 values) and Metatheria (2 values)
      # by missing values
      #---------------------------------------------------------------
    c2 = clade
    c2
    c2[clade %in% levels(clade)[c(1,4)]] = NA
    clade2 = factor(as.numeric(c2))
    
    levels(clade2) = c("E", "L")
    
      #----------------------------------------------
      # Store the column names (for labelling figure) 
      #----------------------------------------------
    pos = colnames(data)
    # pos
    
      #--------------------------------------------------------------
      # Create data frame (data) that ONLY stores the amino acid data
      # Also keep a copy of this, called data.old
      #--------------------------------------------------------------
    data = dd[,-(1:4)]
    # data
    n.columns = ncol(data)
    # n.columns
    data.old = data
    
      #--------------------------------------------------------
      # For each column, identify the two most common sequences
      # and replace the others by missing values
      #--------------------------------------------------------
    
    logr.pmass = wilc.pmass = logr.prank = wilc.prank = 
                 p.clade = numeric(n.columns)
    
      #--------------------------------
      # Run through each column in turn
      #--------------------------------
    
    for (j in 1:(n.columns)) {
    
        #--------------------------------------------------------
        # For each column, identify the two most common sequences
        # and replace the others by missing values
        #--------------------------------------------------------
      cat(j, pos[j], sort(tabulate(data[,j]), decreasing=TRUE), "\n")
      aas = levels(data[,j])
      freqs = tabulate(data[,j])
      of = order(freqs, decreasing=TRUE)
      tmp = numeric(length(data[,j]))
      tmp[!( data[,j] == aas[of][1])] = 1
      tmp[!( data[,j] == aas[of][2])] = 0
      tmp[!( data[,j] %in% aas[of][1:2])] = NA
      tmp = factor(tmp)
      
      levels(tmp) = aas[of][1:2]
      data[,j] = tmp
    
        # The following section of code only works if there are
        # at least two distinct amino acids at the site
        # (What you have is fine, but can also replace with alternative)
    #  if (length(levels(data[,j])) != 1){
      run.j = FALSE
      if (nlevels(data[,j]) > 1) {
        if (freqs[of][2] / freqs[of][1] > threshold) run.j = TRUE
      }
      if (run.j) {
       cat("Running", j, "\n")
    
          #-------------------------------------------------------
          # x = rank position, v = y-variable, coded as 0 or 1 for
          # logistic regression
          #-------------------------------------------------------
        x = 1:nrow(data)
        v = as.numeric(data[,j]) - 1
    
          #------------------------------------------------------
          # Do a Wilcoxon test to see if rank position differs 
          # significantly between the two most common amino acids
          # Save the p-value for this. Also save the slope 
          # parameter from logistic regression
          #------------------------------------------------------
        
        wilc = wilcox.test(x~data[,j])
        wilc.prank[j] = wilc$p.value
     
        g=glm(v ~ x, family="binomial")
        logr.prank[j] = summary(g)$coef[2,4]

          #------------------------------------------------------
          # Do a Wilcoxon test to see if log(mass) differs 
          # significantly between the two most common amino acids
          # Save the p-value for this. Also save the slope 
          # parameter from logistic regression
          # (May get warning messages about ties)
          #------------------------------------------------------
        x = log10(as.numeric(mass))
        wilc = wilcox.test(x~data[,j])
        wilc.pmass[j] = wilc$p.value
      
      
        g=glm(v ~ x, family="binomial")
        logr.pmass[j] = summary(g)$coef[2,4]
      
          #------------------------------------------------------
          # Create a 2x2 table comparing the two clades (E and L)
          # and the two (commonest) types of amino acid
          # Test for a significant association using Fisher's
          # exact test
          #------------------------------------------------------
        tt = tapply(1:nrow(data), list(clade2, data[,j]), length)
        tt[is.na(tt)] = 0
    
       # print(tt)
        p.clade[j] = fisher.test(tt)$p.value
      }
    }
    
    
    #---------------------------
      # Use a 5 x 3 array of plots
      #---------------------------
    # par(mfrow=c(5,3), par(mar=c(4,2.5,2.5,1) + 0.1))
    k=1
    lst = (k-1)*15 + 1:15
    lst = 1:n.columns
    
      #----------------------------------------------------------------
      # For each column, plot the two commonest amino acids against the
      # log(mass). Include the logistic regression line. The p-value in
      # The p-value in the header is from the Wilcoxon test based on 
      # the log-mass. p-value is rounded to 4 decimal places, so '0'
      # just means 'P<0.00005'
      #----------------------------------------------------------------
    for (i in 1:length(lst)) {
      run.i = FALSE
      freqs = tabulate(data[,lst[i]])
      
      if (nlevels(data[,lst[i]]) == 2) {
        if (freqs[2]/freqs[1] > threshold) run.i = TRUE
      }
      if (run.i) {
        writeLines(paste(i,",", -log10(p.clade[i]) ),con=zz,sep="\n")
        # writeLines(paste(i,",", -log10(wilc.pmass[i]) ),con=zz,sep="\n")
        # writeLines(paste(-log10(p.clade[i]) ,",", -log10(wilc.pmass[i]) ),con=zz,sep="\n")
        
        cat("Plotting", i, lst[i], "\n")
        
        headr = paste(pos[lst[i]], ": P = ", 
                      as.character(round(wilc.pmass[lst[i]],4)), sep="")
        v = as.numeric(data[,lst[i]]) - 1
        
        
        most_freq = levels(data[,lst[i]])[1]
        pen_freq = levels(data[,lst[i]])[2]
        switchers <- c(8,15,65,111,211,331,339,348,354,426,429,435,516,576,580,587,592)
        
        
        
        
       
        x = log10(as.numeric(mass))
        max(x)
        min(x)
        
        clade2
        
        # plot(x, v, pch=c(15,2)[as.numeric(clade2)], yaxt="n", main=headr, xlab="", ylab="" )
        # axis(2, at=c(0,1), labels=levels(data[,lst[i]]), las=2)
        g=glm(v ~ x, family="binomial")
        
        
        
    #    lines(x[!is.na(v)], fitted(g), col="red")
        xx = x[!is.na(v)]
        # lines(xx[order(xx)], fitted(g)[order(xx)], col="red")

        # 
  
        }
      } 
    
    
      #--------------------------------------------------------------------
      # These are the plots of p-values. The top plot shows the p-value
      # based on the Fisher test for association with clade. The bottom
      # plot shows the p-value for the Wilcoxon test for difference
      # in log(mass). P-values are shown on a log-scale (log to the base e)
      # The horizontal lines show various P-values
      #--------------------------------------------------------------------
    # jpeg(filename = "p_value_plots_log10_A.jpg",
    #      width = 2000, height = 2000, units = "px", 
    #      quality = 75, res = 200)
    # 
    
    
    
    xpos = 0.9 * ncol(data)
    par(mfrow=c(2,1), par(mar=c(4,4,2,1) + 0.1))
    plot(1:ncol(data), -log10(p.clade), pch=20,
         ylab = expression("-log"["10"]*"(P-Value Clade)"), xlab="Residue")
    title("A", adj = 0.02)
    # points(1:5, c(2,3,4,5,6), col = c(1,2,3))
    axis(1,  at=c(seq(0, 800, by=100),
         labels= c(seq(0, 800, by=100))))
    minor.tick(nx = 10)
    # abline(h=-log(0.05), col="gray")
    # text(xpos, -log(0.05), "P = 0.05", col="gray")
    # abline(h=-log(0.01), col="blue")
    # text(xpos, -log(0.01), "P = 0.01", col="blue")
    # abline(h=-log(0.001), col="cyan")
    # text(xpos, -log(0.001), "P = 0.001", col="cyan")
    # draw.circle(358, 1, 30,border="orange",lty=1,lwd=2)
    # draw.circle(333, 4.65, 10,border="orange",lty=1,lwd=2)
    # draw.ellipse(432, 1.5, 34,2,border="darkorchid1",lty=1,lwd=2)
    # draw.circle(559, 2.55, 10,border="red",lty=1,lwd=2)
    # draw.ellipse(582, 0.65, 14,1,border="red",lty=1,lwd=2)
    
    # abline(h=-log(0.0009434), col="gray")
    # text(xpos, -log(0.0009434)+ 0.44, "P = 9.434 x 10^-4", col="gray")
    # 
    # abline(h=-log(0.00018868), col="blue")
    # text(xpos, -log(0.00018868)+ 0.44, "P = 1.887 x 10^-4", col="blue")
    # 
    # abline(h=-log(0.00001887), col="cyan")
    # text(xpos, -log(0.00001887)+ 0.44, "P = 1.887 x 10^-5", col="cyan")
    abline(h=-log10(0.0009434), col="gray")
    text(xpos+50, -log10(0.0009434)+ 0.24, "P = Adj 0.05", col="gray")
    # 
    abline(h=-log10(0.00018868), col="blue")
    text(xpos+50, -log10(0.00018868)+ 0.24, "P = Adj 0.01", col="blue")
    
    
    # abline(h=-log10(0.00001887), col="cyan")
    # text(xpos, -log10(0.00001887)+ 0.44, "P = Adj 0.001", col="cyan")
    # abline(h=-log(0.001/ncol(data)), v=log(0.05/ncol(data)), col="green")
    # text(1:ncol(data), -log10(p.clade), labels=1:ncol(data), cex= 0.6, pos=2)
    # 
    points(331, -log10(p.clade)[331], pch = 19, col = "orange")
    points(354, -log10(p.clade)[354], pch = 19, col = "orange")
    points(348, -log10(p.clade)[348], pch = 19, col = "orange")
    points(371, -log10(p.clade)[371],  pch = 19, col = "orange")
    
    points(439, -log10(p.clade)[439], pch = 19, col = "darkorchid1")
    points(429,-log10(p.clade)[429], pch = 19, col = "darkorchid1")
    points(435,-log10(p.clade)[435],  pch = 19, col = "darkorchid1")
    points(426,-log10(p.clade)[426],  pch = 19, col = "darkorchid1")
    
    points(580,-log10(p.clade)[580],  pch = 19, col = "red")
    points(587,-log10(p.clade)[587],  pch = 19, col = "red")
    points(576,-log10(p.clade)[576],  pch = 19, col = "red")
    points(560,-log10(p.clade)[560], pch = 19, col = "red")
    points(580,-log10(p.clade)[580], pch = 19, col = "red")


    
    #add labels to chimera residues
    text(331, -log10(p.clade)[331], labels=331, cex= 0.6, pos=2)
    text(354, -log10(p.clade)[354], labels=354, cex= 0.6, pos=1)
    text(348, -log10(p.clade)[348], labels=348, cex= 0.6, pos=2)
    text(371, -log10(p.clade)[371], labels=371, cex= 0.6, pos=2)
    
    text(439, -log10(p.clade)[439], labels=439, cex= 0.6, pos=2)
    text(429, -log10(p.clade)[429], labels=429, cex= 0.6, pos=4)
    text(435, -log10(p.clade)[435], labels=435, cex= 0.6, pos=2)
    text(426, -log10(p.clade)[426], labels=426, cex= 0.6, pos=2)
    
    text(580, -log10(p.clade)[580], labels=580, cex= 0.6, pos=2)
    text(587, -log10(p.clade)[587], labels=587, cex= 0.6, pos=2)
    text(576, -log10(p.clade)[576], labels=576, cex= 0.6, pos=2)
    text(560, -log10(p.clade)[560], labels=560, cex= 0.6, pos=2)

    # dev.off()
    #########################################################p-value mass - residue###############
    # jpeg(filename = "p_value_plots_log10_B.jpg",
    #      width = 2000, height = 2000, units = "px", 
    #      quality = 75, res = 200)
    
    
    
    plot(1:ncol(data), -log10(wilc.pmass), pch=20, 
         ylab = expression("-log"["10"]*"(P-Value Mass)"), xlab="Residue")
    title("B", adj = 0.02)
    
    # draw.circle(355, 17, 65,border="orange",lty=1,lwd=2)
    # 
    # draw.ellipse(432, 13, 38,3,border="darkorchid1",lty=1,lwd=2)
    # draw.circle(582, 16.75, 20,border="yellow",lty=1,lwd=2)
    # draw.ellipse(570, 12.2, 38,5,border="red",lty=1,lwd=2)
    axis(1,  at=c(seq(0, 800, by=100),
                  labels= c(seq(0, 800, by=100))))
    
    minor.tick(nx = 10)
    # abline(h=-log(0.05), col="gray")
    # text(xpos, -log(0.05), "P = 0.05", col="gray")
    # abline(h=-log(0.01), col="blue")
    # text(xpos, -log(0.01), "P = 0.01", col="blue")
    # abline(h=-log(0.001), col="cyan")
    # text(xpos, -log(0.001), "P = 0.001", col="cyan")
    # 
    
    # abline(h=-log(0.0009434), col="gray")
    # text(xpos, -log(0.0009434) + 0.44, "P = 9.434 x 10^-4", col="gray")
    # 
    # abline(h=-log(0.00018868), col="blue")
    # text(xpos, -log(0.00018868)+ 0.44, "P = 1.887 x 10^-4", col="blue")
    # 
    # abline(h=-log(0.00001887), col="cyan")
    # text(xpos, -log(0.00001887) + 0.44, "P = 1.887 x 10^-5", col="cyan")
    abline(h=-log10(0.0009434), col="gray")
    text(xpos+50, -log10(0.0009434) + 0.24, "P = Adj 0.05", col="gray")

    abline(h=-log10(0.00018868), col="blue")
    text(xpos+50, -log10(0.00018868)+ 0.24, "P = Adj 0.01", col="blue")
    
    # abline(h=-log10(0.00001887), col="cyan")
    # text(xpos, -log10(0.00001887) + 0.44, "P = Adj 0.001", col="cyan")
    
    # abline(h=-log(0.001/ncol(data)), v=log(0.05/ncol(data)), col="green")
    # text(1:ncol(data), -log10(wilc.pmass), labels=1:ncol(data), cex= 0.6, pos=2)
    
    points(331, -log10(wilc.pmass)[331], pch = 19, col = "orange")
    points(354, -log10(wilc.pmass)[354], pch = 19, col = "orange")
    points(348, -log10(wilc.pmass)[348], pch = 19, col = "orange")
    points(371, -log10(wilc.pmass)[371], pch = 19, col = "orange")
    
    points(439, -log10(wilc.pmass)[439], pch = 19, col = "darkorchid1")
    points(429, -log10(wilc.pmass)[429], pch = 19, col = "darkorchid1")
    points(435, -log10(wilc.pmass)[435], pch = 19, col = "darkorchid1")
    points(426, -log10(wilc.pmass)[426], pch = 19, col = "darkorchid1")
    
    points(580, -log10(wilc.pmass)[580], pch = 19, col = "red")
    points(587, -log10(wilc.pmass)[587], pch = 19, col = "red")
    points(576, -log10(wilc.pmass)[576], pch = 19, col = "red")
    points(560, -log10(wilc.pmass)[560], pch = 19, col = "red")
    
    
    #add labels to chimera residues
    # text(331, -log10(wilc.pmass)[331], labels=331, cex= 0.6, pos=2)
    # text(354, -log10(wilc.pmass)[354], labels=354, cex= 0.6, pos=2)
    # text(348, -log10(wilc.pmass)[348], labels=348, cex= 0.6, pos=2)
    # text(371, -log10(wilc.pmass)[371], labels=371, cex= 0.6, pos=1)
    # 
    # text(439, -log10(wilc.pmass)[439], labels=439, cex= 0.6, pos=2)
    # text(429, -log10(wilc.pmass)[429], labels=429, cex= 0.6, pos=2)
    # text(435, -log10(wilc.pmass)[435], labels=435, cex= 0.6, pos=2)
    # text(426, -log10(wilc.pmass)[426], labels=426, cex= 0.6, pos=2)
    # 
    # text(580, -log10(wilc.pmass)[580], labels=580, cex= 0.6, pos=2)
    # text(587, -log10(wilc.pmass)[587], labels=587, cex= 0.6, pos=2)
    # text(576, -log10(wilc.pmass)[576], labels=576, cex= 0.6, pos=2)
    # text(560, -log10(wilc.pmass)[560], labels=560, cex= 0.6, pos=2)
    # 
    # 
    
    # dev.off()
    ########################    PLOT P values against each other ####################################
    # jpeg(filename = "p_value_plots_log10_C.jpg",
    #      width = 2000, height = 2000, units = "px", 
    #      quality = 75, res = 200)
    # 
    
    plot(-log10(p.clade), -log10(wilc.pmass), pch=20, 
         ylab = expression("-log"["10"]*"(P-Value Mass)"), xlab=expression("-log"["10"]*"(P-Value Clade)"))
    title("C", adj = 0.02)
    # abline(h=-log(0.05), col="gray")
    # text(-log(0.05), 17, "P = 0.05", col="gray",  srt=90)
    # # abline(h=-log(0.01), col="blue")
    # # text(-log(0.01), 17, "P = 0.01", col="blue",  srt=90)
    # abline(h=-log(0.001), col="cyan")
    # text(-log(0.001), 17, "P = 0.001", col="cyan",  srt=90)
    # # 
    # abline(v=-log(0.05), col="gray")
    # text(20, -log(0.05), "P = 0.05", col="gray")
    # # abline(v=-log(0.01), col="blue")
    # # text(20, -log(0.01), "P = 0.01", col="blue")
    # abline(v=-log(0.001), col="cyan")
    # text(20, -log(0.001), "P = 0.001", col="cyan")
    
    abline(h=-log10(0.00018868), col="blue")
    text(-log10(0.00018868)-0.14, 7.5, "P = Adj 0.01", col="blue",  srt=90)
    
    abline(h=-log10(0.0009434), col="gray")
    text(-log10(0.0009434) - 0.14, 7.5, "P = Adj 0.05", col="gray", srt=90)
    # abline(h=-log(0.00018868), col="blue")
    # text(-log(0.00018868)-0.44,16, "P = 1.887 x 10^-4", col="blue",srt=90)
    
    # abline(h=-log10(0.00001887), col="cyan")
    # text(-log10(0.00001887)-0.44,16, "P = Adj 0.001", col="cyan",srt=90)
    
    abline(v=-log10(0.00018868), col="blue")
    text(8.7, -log10(0.00018868)+0.24, "P = Adj 0.01", col="blue")
    
    abline(v=-log10(0.0009434), col="gray")
    text(8.7, -log10(0.0009434) + 0.24, "P = Adj 0.05", col="gray")
    # abline(v=-log(0.00018868), col="blue")
    # text(18, -log(0.00018868)+0.44, "P = 1.887 x 10^-4", col="blue")
    
    # abline(v=-log10(0.00001887), col="cyan")
    # text(18, -log10(0.00001887)+0.44, "P = Adj 0.001", col="cyan")
    
    # text(-log10(p.clade), -log10(wilc.pmass), labels=1:ncol(data), cex= 0.6, pos=2)
    points(-log10(p.clade)[331], -log10(wilc.pmass)[331], pch = 19, col = "orange")
    points(-log10(p.clade)[354], -log10(wilc.pmass)[354], pch = 19, col = "orange")
    points(-log10(p.clade)[348], -log10(wilc.pmass)[348], pch = 19, col = "orange")
    points(-log10(p.clade)[371], -log10(wilc.pmass)[371], pch = 19, col = "orange")
  
    points(-log10(p.clade)[439], -log10(wilc.pmass)[439], pch = 19, col = "darkorchid1")
    points(-log10(p.clade)[429], -log10(wilc.pmass)[429], pch = 19, col = "darkorchid1")
    points(-log10(p.clade)[435], -log10(wilc.pmass)[435], pch = 19, col = "darkorchid1")
    points(-log10(p.clade)[426], -log10(wilc.pmass)[426], pch = 19, col = "darkorchid1")
    
    points(-log10(p.clade)[580], -log10(wilc.pmass)[580], pch = 19, col = "red")
    points(-log10(p.clade)[587], -log10(wilc.pmass)[587], pch = 19, col = "red")
    points(-log10(p.clade)[576], -log10(wilc.pmass)[576], pch = 19, col = "red")
    points(-log10(p.clade)[560], -log10(wilc.pmass)[560], pch = 19, col = "red")
    
    text(-log10(p.clade)[331], -log10(wilc.pmass)[331], labels=331,cex= 0.6, pos=4)
    text(-log10(p.clade)[354], -log10(wilc.pmass)[354], labels=354,cex= 0.6, pos=4)
    text(-log10(p.clade)[348], -log10(wilc.pmass)[348], labels=348,cex= 0.6, pos=3)
    text(-log10(p.clade)[371], -log10(wilc.pmass)[371],labels=371, cex= 0.6, pos=3)
    
    text(-log10(p.clade)[439], -log10(wilc.pmass)[439], labels=439,cex= 0.6, pos=3)
    text(-log10(p.clade)[429], -log10(wilc.pmass)[429], labels=429,cex= 0.6, pos=4)
    text(-log10(p.clade)[435], -log10(wilc.pmass)[435], labels=435,cex= 0.6, pos=3)
    text(-log10(p.clade)[426], -log10(wilc.pmass)[426], labels=426,cex= 0.6, pos=4)
    
    text(-log10(p.clade)[580], -log10(wilc.pmass)[580], labels=580,cex= 0.6, pos=4)
    text(-log10(p.clade)[587], -log10(wilc.pmass)[587], labels=587,cex= 0.6, pos=4)
    text(-log10(p.clade)[576], -log10(wilc.pmass)[576], labels=576,cex= 0.6, pos=4)
    text(-log10(p.clade)[560], -log10(wilc.pmass)[560], labels=560,cex= 0.6, pos=3)
    
    text(-log10(p.clade)[15], -log10(wilc.pmass)[15], labels=15,cex= 0.6, pos=4)
    text(-log10(p.clade)[638], -log10(wilc.pmass)[638], labels=638,cex= 0.6, pos=4)
    text(-log10(p.clade)[11], -log10(wilc.pmass)[11], labels=11,cex= 0.6, pos=4)
    text(-log10(p.clade)[52], -log10(wilc.pmass)[52], labels=52,cex= 0.6, pos=3)
    text(-log10(p.clade)[614], -log10(wilc.pmass)[614], labels=614,cex= 0.6, pos=3)
    text(-log10(p.clade)[65], -log10(wilc.pmass)[65], labels=65,cex= 0.6, pos=3)
    
    dev.off()

    
  