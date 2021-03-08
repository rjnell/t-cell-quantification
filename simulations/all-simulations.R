# Load library 'digitalPCRsimulations'
library(digitalPCRsimulations)

# Set random seed
set.seed(12345)

# Function to calculate the difference between two concentrations,
# obtained by taking a given sample from two given universes.
sample_from_2_universes_diff = function(universe1, universe2, sample, alpha=0.95) {
  
  # What is the total length of the sample?
  n = length(sample)
  
  # As described in the manuscript:
  # concentration difference target1 - target2 =
  # [1] - [2] = -log(p_1_min / p_2_min) / volume

  p_1_min = (n - sum(universe1[sample])) / n
  p_2_min = (n - sum(universe2[sample])) / n
  
  # Translate this percentage to concentration (difference)...
  conc = -log(p_1_min / p_2_min)
  
  # ... and its CI
  # As described in manuscript:
  # Var = (p_1_min^-1 + p_2_min^-1 - 2) / n
  var = (p_1_min^(-1) + p_2_min^(-1) - 2)/n
  
  # ... and s = sqrt(var), following Dube et al. 2008
  s = sqrt(var)
  z = qnorm(1-(1-alpha)/2)
  conc_low = conc - z * s
  conc_high = conc + z * s
  
  # Volume per droplet
  volume = 0.00085
  
  # Return the calculated concentration (difference) with CI
  conc_ci = c(conc/volume, conc_low/volume, conc_high/volume)  
  names(conc_ci) = c("concentration", "concentration_low", "concentration_high")
  return(conc_ci)
}

# Function to simulate multiplex adjusted model
simulate_multiplex_adjusted_model = function(input_ng, 
                                             tcf, 
                                             cnv_ref, 
                                             cnv_region, 
                                             n_droplets, 
                                             n_simulations, 
                                             alpha) {
  
  # The amount of input alleles of the reference depends on the chosen amount of input DNA, 
  # the T-cell fraction and the CNV of the reference.
  # As in T cells (= healthy cells) a default CNV of 2 is found, 
  # the CNV of the reference can only be altered due to CNA in the non-T cells.
  input_ng_ref = input_ng * (tcf * 2 + (1-tcf) * cnv_ref) / 2
  universe_ref = universe(input_ng_ref)
  
  # The amount of input alleles of the regional marker depends on the chosen amount of input DNA, 
  # the T-cell fraction and the regional CNV.
  # In T cells (= healthy cells) a default CNV of 2 is found, 
  # while in non-T cells we observe the regional CNV.
  input_ng_region = input_ng * (tcf * 2 + (1-tcf) * cnv_region) / 2
  universe_region = universe(input_ng_region)
  
  # The amount of input alleles of the T-cell marker depends on the chosen amount of input DNA, 
  # the T-cell fraction and the regional CNV.
  # In T-cells the marker is completely lost, 
  # while in non-T cells the regional CNV determines the absolute presence of the T-cell marker.
  input_ng_tcm = input_ng * (1-tcf) * cnv_region / 2
  universe_tcm = universe(input_ng_tcm)
  
  # Initialize the list of results of the simulations
  results = list()
  results$tcf_true_value = tcf
  results$tcf_classic = matrix(NA, nrow = 0, ncol = 3)
  results$tcf_adjusted = matrix(NA, nrow = 0, ncol = 3)
  
  for (simulation in 1:n_simulations) {
    
    # Simulate duplex and multiplex
    sample = simulate_sample(n_droplets)
    
    # Determine the concentration of the reference 
    conc_ref = sample_from_universe(universe_ref, 
                                    sample, 
                                    alpha)
    
    # Determine the concentration of the regional marker
    conc_region = sample_from_universe(universe_region, 
                                       sample, 
                                       alpha)
    
    # Determine the concentration of the T-cell marker
    conc_tcm = sample_from_universe(universe_tcm, 
                                    sample, 
                                    alpha)
    
    # Determine the concentration of [regional marker - T-cell marker]
    conc_diff_region_tcm = sample_from_2_universes_diff(universe_region, 
                                                        universe_tcm, 
                                                        sample, 
                                                        alpha)
    
    # Calculate the T-cell fraction according the classic model, 
    # using the T-cell marker and the reference
    tcf_classic = 1-calc_ratio(conc_tcm, conc_ref)[c(1,3,2)]
    
    # Calculate the T-cell fraction according the adjusted model, 
    # using the T-cell marker, the regional marker as corrector,
    # and reference 1.
    tcf_adjusted = calc_ratio(conc_diff_region_tcm, conc_ref) 
    
    # Bind obtained result to the simulation results list
    results$tcf_classic = rbind(results$tcf_classic, tcf_classic) 
    results$tcf_adjusted = rbind(results$tcf_adjusted, tcf_adjusted)
    
  }
  
  # Return obtained simulation_results
  return (results)
}

# Shortcut function to plot summary of simulations
plot_simulations_summary = function(main, simulation_results, true_tcf) {
  
  par(mar=c(5.1, 5, 4.1, 2.1))
  ylim = c(0,60)
  xlim = c(-1,22)
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
  
  yat = seq(from=ylim[1],to=ylim[2],by=10)
  xat = c(xlim[1],xlim[2])
  segments(xlim[1],yat,xlim[2],col="#eeeeee",lwd=1.4,xpd=T)
  segments(xlim[1],50,xlim[2],col="#b1b1b1",lwd=1.4,xpd=T,lty=3)
  segments(xlim[1],ylim[1],xlim[1],ylim[2]+3.75,xpd=T,col="#B1B1B1",lwd=1.4)
    axis(side = 2,at=yat,las=2,labels=rep("",length(yat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = -0.035,lwd=1.4)
  axis(side = 2,at=yat,las=2,labels=paste0(yat,"%"), lwd=0, col.axis="#333333",line=-0.23)
  mtext(text = "T-cell fraction", side = 2, line=3.7,col="#333333")  
  
  r = simulation_results
  s = stats(r, true_value=true_tcf)
  
  m = format(round(s$point_estimate_mean*100, 1), nsmall = 1)
  sd = format(round(s$point_estimate_sd*100, 1), nsmall = 1)
  cc= format(round(s$coverage*100, 1), nsmall = 1)
  
  axis(side = 1, at=xat, labels=rep("",length(xat)),col.ticks = "#b1b1b1",col = "#b1b1b1", tck = 0,lwd=1.4)
  mtext(text = paste0("Mean (sd) = ", m, "% (", sd,")"), side = 1, line=1.2,col="#333333")  
  mtext(text = paste0("95%-CI coverage = ", cc, "%"), side = 1, line=2.2,col="#333333")  
  mtext(text = main, side = 3, line=1.2,col="#333333", font=4)  
  
  segments(xlim[1],ylim[1],xlim[2],ylim[1],xpd=T,col="#B1B1B1",lwd=1.4)
  
  results = r[1:20,]*100
  results = results[order(results[,1]),]
  
  lty=rep(1,20)
  col=rep("#b1b1b1",20)
  w=which(results[,2]>50 | results[,3] <50)
  lty[w] = 1
  col[w] = "#FBB4AF"
  
  arrows(1:nrow(results), results[,2], 1:nrow(results), results[,3], length=0.05, angle=90, code=3, col=col,lty=lty,lwd=1.4)
  points(x = 1:nrow(results),
         y = results[,1],
         pch = 15, cex=0.9, col="#333333")
  
}


### SIMULATION 1: 20 ng, copy number stable, 50% TCF

simulation_1 = simulate_multiplex_adjusted_model(
  
  # Input in ng genomic DNA
  input_ng = 20, 
  
  # T-cell fraction (0-1)
  tcf = 0.50, 
  
  # Copy number value of the reference in non-tcf (def = 2)
  cnv_ref = 2, 
  
  # Copy number value of the T-cell marker region in non-tcf (def = 2)
  cnv_region = 2, 
  
  # Number of droplets per experiment
  n_droplets = 20000, 
  
  # Number of simulations
  n_simulations = 10000, 
  
  # Level of significance
  alpha = 0.95)

# Plot summary of simulation results of classic model
png("simulations/fig-2a-1.png", res=600, width=2000, height=2500)  
plot_simulations_summary("Classic model", simulation_1$tcf_classic, 0.5)
dev.off()

# Plot summary of simulation results of adjusted model
png("simulations/fig-2a-2.png", res=600, width=2000, height=2500)  
plot_simulations_summary("Adjusted model", simulation_1$tcf_adjusted, 0.5)
dev.off()


### SIMULATION 2: 20 ng, copy number unstable, 50% TCF

simulation_2 = simulate_multiplex_adjusted_model(
  
  # Input in ng genomic DNA
  input_ng = 20, 
  
  # T-cell fraction (0-1)
  tcf = 0.5, 
  
  # Copy number value of the reference in non-tcf (def = 2)
  cnv_ref = 2, 
  
  # Copy number value of the T-cell marker region in non-tcf (def = 2)
  cnv_region = 3, 
  
  # Number of droplets per experiment
  n_droplets = 20000, 
  
  # Number of simulations
  n_simulations = 10000, 
  
  # Level of significance
  alpha = 0.95)

# Plot summary of simulation results of classic model
png("simulations/fig-2b-1.png", res=600, width=2000, height=2500)  
plot_simulations_summary("Classic model", simulation_2$tcf_classic, 0.5)
dev.off()

# Plot summary of simulation results of adjusted model
png("simulations/fig-2b-2.png", res=600, width=2000, height=2500)  
plot_simulations_summary("Adjusted model", simulation_2$tcf_adjusted, 0.5)
dev.off()