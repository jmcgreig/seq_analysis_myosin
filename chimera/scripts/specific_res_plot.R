    library(Hmisc) 
    library(ncar)
    #-----------------------------------------------
    # Store graphics in a pdf file. Comment out this
    # line incar v0.3.4f you want plots to the screen
    #-----------------------------------------------
    # pdf("Specific_Resi.pdf", paper="a4", width = 8.3, height = 11.7)
    # pdf("Supplementary_Resi.pdf", paper="a4", width = 8.3, height = 11.7)
    
    test_count = 0
    jpeg(filename = "../figures/specific_resi.jpg",
         width = 1500, height = 2000, units = "px",
         quality = 75, res = 200)
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
    # by missing valuesswitchers <- c(331,339,348,354)
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
    
    for (j in 1:n.columns) {
      
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
    par(mfrow=c(3,2), par(mar=c(4,2.5,2.5,1) + 0.1))
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
        
        most_freq = levels(data[,lst[i]])[1]
        pen_freq = levels(data[,lst[i]])[2]
        important <- c(77,125,331,339,348,354)
        # important <- c(426, 808)
        
        switchers <- c(8,15,65,111,211,331, 339,348, 354, 426, 429,435,516,576,580,587,592)
        # switchers <- c(331,339,348,354)
        # switchers <- c(808)
        
        if ( (i %in%  important)) {
          
          test_count = test_count + 1
          # if (test_count == 13){
          #   dev.off()
          #   jpeg(filename = "supplementary_resi_1.jpg",
          #        width = 2000, height = 2000, units = "px",
          #        quality = 75, res = 200)
          #   par(mfrow=c(4,3), par(mar=c(4,2.5,2.5,1) + 0.1))
          # }
          # if (test_count == 25){
          #   dev.off()
          #   jpeg(filename = "supplementary_resi_2.jpg",
          #        width = 2000, height = 2000, units = "px",
          #        quality = 75, res = 200)
          #   par(mfrow=c(4,3), par(mar=c(4,2.5,2.5,1) + 0.1))
          # }
          # if (test_count == 37){
          #   dev.off()
          #   jpeg(filename = "supplementary_resi_3.jpg",
          #        width = 2000, height = 2000, units = "px",
          #        quality = 75, res = 200)
          #   par(mfrow=c(4,3), par(mar=c(4,2.5,2.5,1) + 0.1))
          # }
          # if (test_count == 49){
          #   dev.off()
          #   jpeg(filename = "supplementary_resi_4.jpg",
          #        width = 2000, height = 2000, units = "px",
          #        quality = 75, res = 200)
          #   par(mfrow=c(4,3), par(mar=c(4,2.5,2.5,1) + 0.1))
          # }
          most_freq = levels(data[,lst[i]])[1]
          pen_freq = levels(data[,lst[i]])[2]
          if (i %in% switchers) {
            
          x1  = factor(data[,lst[i]], levels=c(pen_freq, most_freq))
          
          cat("Plotting", i, lst[i], "\n")
          
          # headr = paste("426b: P = ",
          #               as.character(round(wilc.pmass[lst[i]],4)), sep="")
          headr = paste(pos[lst[i]], ": P = ",
                        as.character(round(wilc.pmass[lst[i]],4)), sep="")
          v = as.numeric(x1) - 1
          # v = as.numeric(data[,lst[i]]) - 1
          # x1 = factor(data[,lst[i]], levels=c(most_freq, pen_freq))
          
          } 
          else {
            cat("Plotting", i, lst[i], "\n")
            
            headr = paste(pos[lst[i]], ": P = ", 
                          as.character(round(wilc.pmass[lst[i]],4)), sep="")
            v = as.numeric(data[,lst[i]]) - 1
            x1 = factor(data[,lst[i]], levels=c(most_freq, pen_freq))
        }
    
          
          
          x = log10(as.numeric(mass))
          levels(data[,lst[i]])
          levels(x1)
          plot(x, v, pch=c(15,2)[as.numeric(clade2)], yaxt="n", main=headr, xlab="", ylab="" )
          axis(2, at=c(0,1), labels=levels(x1), las=2)
          g=glm(v ~ x, family="binomial")
          
          
          
          #    lines(x[!is.na(v)], fitted(g), col="red")
          xx = x[!is.na(v)]
          lines(xx[order(xx)], fitted(g)[order(xx)], col="red")
          
        }
      
        
      }
    } 
    
    dev.off()
    
