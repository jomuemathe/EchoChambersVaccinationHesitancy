#
# try to use a convex combination of beta distributions
# to explain the data at hand.
#

### Input data ######
# vaccDat: 
#     Vector of data \in (0,1); relative abundance of opinion 1 in district i
# pairings:
#    array, two columns. one row (i,j) indicates that district i is neighbot of
#    district j (where i,j, refer to the i'th resp. j'th entry of the 
#    vector myDataEsti)
#

###############
# 
# read data
#
#################


#
# adapt!
#
# set work directory
# setwd("...");
#
# load vaccination frequencies (vector, entries between 0 and 1)
# vaccDat = ....
#
# load spatial data (pairings, array, two columns, integer values
# (i,j) does mean vaccData[i] and vaccData[j] come form neighboring districts
# pairings = ....
#
infect.x    = "Measles";  # measels or whatever




source("estiReinforceSpatialV0.txt");





######################################
#
#  fit the parameters (Quasi-Likelihood-based)
#
#####################################





   # fit the beta model as well as the reinforcement model;
   # do that based on the spatial model
   #
   res.tab = c();
   

         myDataEsti = vaccData;
        
         max.x = 1.0;
         hist(myDataEsti, freq = FALSE, nclass=30, 
              main=as.character(e$year),
              xlim=c(0,max.x));

         ###################################################################
         # first run: estimate the full reinforcement model
         ###################################################################
         # start with the homogeneous model
         source("estiReinforceV0.txt");
         prepare.esti();
         para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
         lll.last = lll(para.ref);
         cat(lll.last, "\n");
         curve(g(x), add=TRUE, col="blue", lwd=2);

         unrestricted.model = TRUE;      # we aim at the full model
         unrestrict.theta   = TRUE;  
         cat("homogeneous model\n");   
         opti.cylcic(para.ref);
         cat(lll.last, "\n");

         # now turn over to the spatial model
         source("estiReinforceSpatialV0.txt");
         maxIter = 50;
         para.ref = c(para.last, 8);   # define init para
         lll.last = lll(para.ref);

         cat("spatial model\n");
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
         if (1==0){ # not used
            para.ref = c(mean(myDataEsti), 0.5, 0.5,100, 1);   # define init para
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
         }
         
         ###################################################################
         # third run: estimate the zealot model (beta-distrib)
         ###################################################################
         if (1==0) { # not used
            # start with the homogeneous model
            source("estiReinforceV0.txt");
            prepare.esti();
            para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
            lll.last = lll(para.ref);
            cat(lll.last, "\n");
            curve(g(x), add=TRUE, col="blue", lwd=2);

            unrestricted.model = FALSE;           # we fix all reinfrcement-paras
            unrestrict.theta   = FALSE;
            cat("homogeneous model\n");   
            opti.cylcic(para.ref);
            cat(lll.last, "\n");

            # now turn over to the spatial model
            source("estiReinforceSpatialV0.txt");
            para.ref = c(para.last, 8);   # define init para
            lll.last = lll(para.ref);
   
            cat("spatial model\n");
            unrestricted.model = FALSE;           # we fix all reinfrcement-paras
            unrestrict.theta   = FALSE;
            para.ref = c(mean(myDataEsti), 0.0, 0.5,100, 1);
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
         }         

         line = c(e$year, infect.x, 
                  theta.res,
                  para.unrest,
                  lll.unrest,
                 # para.halfRestrict,            # not used
                 # lll.halfRestrict,
                 #  para.restrict,
                 # lll.restrict, 
                 # pchisq(2*(lll.unrest-lll.restrict),df=2, 
                 #                lower.tail=FALSE), 
                 # pchisq(2*(lll.unrest-lll.halfRestrict),df=2, 
                 #                lower.tail=FALSE), 
                  res.ks$p.value
                 # res.halfRestrict.ks$p.value,
                 # res.restric.ks$p.value, mean(myDataEsti)
                  );

         res.tab = rbind(res.tab, line); 
   }

   # names orient themseves ar the supplement II of the paper
   col.names = c(
           "year", "infect", 
           "Theta1PlusTheta2.unr",
           "nu.unr", "theta.hat.unr", "psi.unr", "s.unr","gamma.unr",
           "lll.unr",
           # "nu.halfr", "theta.hat.halfr", "psi.halfr", "s.halfr","gamma.halfr",
           # "lll.halfr",
           # "nu.restr", "theta.hat.restr", "psi.restr", "s.restr","gamma.restr",
           # "lll.restr",
           # "lll.unr.rest", "lll.unr.halfr",
           "ks.unres"
           #"ks.halfRestr", "ks.restr", "mean"
     );
   dimnames(res.tab)[[2]] =col.names;

   # save resulting parameters
   save(file="spatFitResultV1.rSave", res.tab);









