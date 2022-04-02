#
# try to use a convex combination of beta distributions
# to explain the data at hand.
#
# please load the data you want to analyze as a vector of 
# vaccination frequencies (all enties >0 and <1.0),
# stored in the variable "vaccData"
# thereto, please adapt the section "read data"






###############
# 
# read data
#
#################

# perhaps you want to set the work directiory
# setwd("...");


# Please adapt this line:
# vaccDat is a vector with faccination frequencies you want to fit
# vaccDat = XXXXXXX;
infect.x    = "Measles"; # name of the infection

# a greedy-like optimization algorithm used to ft the data:
source("estiReinforceV0.txt");



######################################
#
#  do it - fit the data
#
#  We use three different setting:
#  first run; Full model (no restrictions)
#  second run: reinforcement parameters are identical for both groups
#  third run: reinforcement parameter is taken to zero (beta distribution).
#
#  We use kolmogoroff-smirov test to check for appropriateness
#  and likelihod-ratio test to compare the three settings.
#
#####################################

#
# Reparametrization of the distribution
#
# \theta_1 = skal * theta.hat * psi
# \theta_2 = skal * theta.hat * (1-psi)
# N_1+1    = skal * (1-theta) * nu 
# N_2+1    = skal * (1-theta) * (1-nu)
# where skal>0, and theta.hat, psi, n1\in[0,1].
# the vector para holds the parameters in the sequence
# para = (n1, theta, theta1, skal)



         # fit the beta model as well as the reinforcement model
         # neglecting space
         #

         res.tab = c();
     
         all.para = c();

         myDataEsti  = vaccDat;
   

         mmean = mean(myDataEsti);

         max.x = 1.0;
         hist(myDataEsti, freq = FALSE, nclass=30, 
              main=as.character(e$year),
              xlim=c(0,max.x));

         ###################################################################
         # first run: estimate the full reinforcement model
         ###################################################################
         para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
         lll.last = lll(para.ref);
         cat(lll.last, "\n");
         curve(g(x), add=TRUE, col="blue", lwd=2);

         unrestricted.model = TRUE;      # we aim at the full model
         unrestrict.theta   = TRUE;     
         opti.cylcic(para.ref);
         cat(lll.last, "\n");
         para.unrest = para.last;        # store the result
         lll.unrest  = lll.last;
         theta.res   = theta*skal;

         # Kolmogorov-Smirnov test
         res.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 
 
  

         ###################################################################
         # second run: reinforcement model, force equal reiforcement parameters 
         ###################################################################
         para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
         lll.last = lll(para.ref);
         cat(lll.last, "\n");
         hist(myDataEsti, freq = FALSE, nclass=30, 
              main=as.character(e$year),
              xlim=c(0,max.x));
          curve(g(x), add=TRUE, col="blue", lwd=2);

         unrestricted.model = TRUE;      # we aim at the full model
         unrestrict.theta   = FALSE;     # we want to keep equal parameters for reinforcement
         opti.cylcic(para.ref);
         cat(lll.last, "\n");
         para.halfRestrict = para.last;        # store the result
         lll.halfRestrict  = lll.last;
         theta.res   = theta*skal;

         # Kolmogorov-Smirnov test
         res.halfRestrict.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 

         
         ###################################################################
         # third run: estimate the zealot model (beta-distrib)
         ###################################################################

         unrestricted.model = FALSE;           # we fix all reinfrcement-paras
         unrestrict.theta   = FALSE;
         para.ref = c(mean(myDataEsti), 0.0, 0.5,100);
         hist(myDataEsti, freq = FALSE, nclass=30, 
              main=as.character(e$year),
              xlim=c(0,max.x));
         curve(g(x), add=TRUE, col="blue", lwd=2);

         opti.cylcic(para.ref);
         cat(lll.last, "\n");
         para.restrict = para.last;            # store result
         lll.restrict  = lll.last;

         # komogorov-smirnov-test
         res.restric.ks = ks.test(myDataEsti, function(x){pReinforce(x)});
         
         line = c(e$year, infect.x, 
                  theta.res,
                  para.unrest,
                  lll.unrest,
                  para.halfRestrict,
                  lll.halfRestrict,
                  para.restrict,
                  lll.restrict, 
                  pchisq(2*(lll.unrest-lll.restrict),df=2, 
                                 lower.tail=FALSE), 
                  pchisq(2*(lll.unrest-lll.halfRestrict),df=2, 
                                 lower.tail=FALSE), 
                  res.ks$p.value,
                  res.halfRestrict.ks$p.value,
                  res.restric.ks$p.value, mean(myDataEsti)
                  );

         res.tab = rbind(res.tab, line); 
   

   # names are the rescaled parameters (see above how to obtain the original paramers
   # from the rescaled parameters)
   col.names = c(
           "year", "infect", 
           "Theta1PlusTheta2.unr",
           "nu.unr", "theta.hat.unr", "psi.unr", "skal.unr",
           "lll.unr",
           "nu.halfr", "theta.hat.halfr", "psi.halfr", "skal.halfr",
           "lll.halfr",
           "nu.restr", "theta.hat.restr", "psi.restr", "skal.restr",
           "lll.restr",
           "lll.unr.rest", "lll.unr.halfr",
           "ks.unres", "ks.halfRestr", "ks.restr", "mean");
   dimnames(res.tab)[[2]] =col.names;

   # save the result   
   save(file="estimatedParameters.rSave", res.tab);




