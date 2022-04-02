#
# Estimate the paras for the reinforcement model
# Spatial setting; we use the quasi-likelihood
#
### Input data ######
# myDataEsti: 
#     Vector of data \in (0,1); relative abundance of opinion 1 in district i
# pairings:
#    array, two columns. one row (i,j) indicates that district i is neighbot of
#    district j (where i,j, refer to the i'th resp. j'th entry of the 
#    vector myDataEsti)
#
###### prepare ##############
# call prepare.esti() before the estimation starts
# to set up partivularly the spatial structure in the way required 
# by the present script
#
#
###### parametr of the model ####
#   n1     <<- para[1];
#   theta  <<- para[2];
#   theta1 <<- para[3];
#   skal   <<- min(skal.max,abs(para[4]));
#   ggamma <<- para[5];
#
#
#
#####################################################################
#
#    prepare the data dor the estimation in the spatial reinforcement model
#
#####################################################################
#
#
# Aim: 
#
# my.neighbors: 
#    list of length of the number of districts; 
#    my.neighbors[[i]] is a vector with the neighbors of district i.
# my.neighbors.vals:
#    list of length of the number of districts; 
#    my.neighbors.vals[[i]] is a vector with the vlalues of the neighbors 
#    of district i.
# my.effecive.influence:
#    vector, \sum_{k\im k'} (x^{(k)}-x^{(k')})^2/d_k
# my.degree:
#    vector, includes the degree for each node (order of vactor myDataEsti).


prepare.esti <- function(){
   my.neighbors      = list(c());
   my.neighbors.vals = list(c());
   my.degree         = rep(0, length(myDataEsti));

   for (i in 2:length(myDataEsti)){
      my.neighbors     = c(my.neighbors, list(c()));
      my.neighbors.vals = c(my.neighbors.vals, list(c()));
   }
   for (i in 2:length(pairings[,1])){
      a = pairings[i,1]; b = pairings[i,2];
      my.neighbors[[a]]      = c(my.neighbors[[a]], b);
      my.neighbors[[b]]      = c(my.neighbors[[b]], a);
      my.neighbors.vals[[a]] = c(my.neighbors.vals[[a]], myDataEsti[b]);
      my.neighbors.vals[[b]] = c(my.neighbors.vals[[b]], myDataEsti[a]);
      my.degree[a] = my.degree[a]+1;
      my.degree[b] = my.degree[b]+1;
   }
   
   my.effecive.influence = rep(NA, length(myDataEsti));
   for (i in 1:length(myDataEsti)){
      my.effecive.influence[i] = sum((my.neighbors.vals[[i]]-myDataEsti[i])**2);
      my.effecive.influence[i] = my.effecive.influence[i]/(4*my.degree[i]);
   }
   # store in environment
   my.neighbors          <<- my.neighbors;
   my.neighbors.vals     <<- my.neighbors.vals;
   my.degree             <<- my.degree;
   my.effecive.influence <<- my.effecive.influence;
   mmean                 <<- mean(myDataEsti);

}





#####################################################################
#
#    estimate the paras of the reinforcement model - general functions
#
#####################################################################


##################
#
# define distrib
#
##################

skal.max = 2000;

f.norm <- function(){
   # an additive constant to avoid large numbers
   # in the exponent; this number is chosen in dependence
   # on the mean value of the data at hand.
   x = mmean;
   return(-1*(skal*theta*(0.5*x**2-theta1*x)
          +log(x)*(1-theta)*n1*skal
          +log(1-x)*(1-theta)*(1-n1)*skal)
          );
}

f.loc <- function(x){
  # N1 = (1-theta)*n1+1, N2 = (1-theta)*(1-n1)*scal +1
  # theta \in (0,1), n1 \in (0,1), scal >0
  return(
      exp(  skal*theta*(0.5*x**2-theta1*x)
          +log(x)*(1-theta)*n1*skal
          +log(1-x)*(1-theta)*(1-n1)*skal +f.norm()
          -sum((neighbVals-x)**2)*ggamma/(4*degreeVal)
         )
   );
}
f <- function(x){return(sapply(x, f.loc))}

get.cc <- function(){
   # get the normalizazion constant
   CC=(integrate(f, lower=10**-6,upper=1-10**-6)$value)**(-1);
   return(CC);
}


g.loc<-function(x){
   # density (note: CC is computed using the constant f.norm().
   # That constant cancels at the end of the day.
   return(
      CC*
       exp(skal*theta*(0.5*x**2-theta1*x)
          +log(x)*(1-theta)*n1*skal
          +log(1-x)*(1-theta)*(1-n1)*skal+f.norm()
          -sum((neighbVals-x)**2)*ggamma/(4*degreeVal)
         )
   );
}
g <- function(x){return(sapply(x, g.loc))}

lll.dat <- function(x,n1,theta,theta1,skal, ggamma, CC, districtNo){
   # Quasi-log likeli for one single data point x,
   # given the data parameter, and the normalization constant CC
   # N1 = (1-theta)*n1+1, N2 = (1-theta)*(1-n1)*scal +1
   # theta \in (0,1), n1 \in (0,1), scal >0
   # district with number distrcictNo
  
   return(
           skal*theta*(0.5*x**2-theta1*x)
          +log(x)*(1-theta)*n1*skal
          +log(1-x)*(1-theta)*(1-n1)*skal
          +log(CC)+f.norm()
          -my.effecive.influence[districtNo]*ggamma
   );
}


pReinforce.loc <- function(x) {
   # parameters given by global parameters;
   # particularly, CC is defined.
   res = integrate(g,lower = 0, upper = x);
   return(res$value);
}

pReinforce <- function(x){
   return(sapply(x, pReinforce.loc));
}

#######################
# estimate paras:
# N1, N2, theta1, theta2
########################
theta.max = 1800;
theta.max = 1900;
tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
				     warning = w.handler),
	 warning = W)
}

# compute marginal log likeli for district myDatassmywith number districtNo
lll.district <- function(para, districtNo){
   # the parameter already is defined in the environment
   
   # indicate the district globally
   targetDist <<- districtNo;
   neighbVals <<- my.neighbors.vals[[districtNo]];
   degreeVal  <<- my.degree[districtNo]

   # get normalization constant
   OK = TRUE;
   aa = tryCatch.W.E(get.cc());
   aa <<- aa;
   if (!(is.double(aa$value))>0) return(-10000);
   CC<<- aa$value;

   # return likelihood 
   n1     <- para[1];
   theta  <- para[2];
   theta1 <- para[3];
   skal   <- min(skal.max,abs(para[4]));
   ggamma <- para[5];
   return(lll.dat(myDataEsti[districtNo], n1,theta,theta1,skal,ggamma, CC, districtNo));   
}


# compute log likeli
lll <- function(para){
   # N1 = (1-theta)*n1+1, N2 = (1-theta)*(1-n1)*scal +1
   # theta \in (0,1), n1 \in (0,1), scal >0
   ppara  <<- para;
   n1     <<- para[1];
   theta  <<- para[2];
   theta1 <<- para[3];
   skal   <<- min(skal.max,abs(para[4]));
   ggamma <<- para[5];

   res = 0;
   for (i in 1:length(myDataEsti)){
      x =   lll.district(para, i);
      if (x< -9999) return(x);  # integration / C failed
      res = res + x;
   }
   return(res);    
}


#
# for the optimization: vary only one of the parameters
#
search.p1 <- function(px){
   # para.last gives the framework; we modify parameter 1 only
   pxx = para.last; pxx[1] = px;
   return( lll(pxx) );
}
search.p2 <- function(px){
   # para.last gives the framework; we modify parameter 2 only
   pxx = para.last; pxx[2] = px;
   return( lll(pxx) );
}
search.p3 <- function(px){
   # para.last gives the framework; we modify parameter 3 only
   pxx = para.last; pxx[3] = px;
   return( lll(pxx) );
}
search.p4 <- function(px){
   # para.last gives the framework; we modify parameter 4 only
   pxx = para.last; pxx[4] = px;
   return( lll(pxx) );
}
search.p5 <- function(px){
   # para.last gives the framework; we modify parameter 5 only
   pxx = para.last; pxx[5] = px;
   return( lll(pxx) );
}

############################################
#
# optimization
#
############################################

opti.cylcic <- function(para.init){
   # optimize cylcically the parameters.
   #
   # we have different modes 
   # unrestricted.model == FALSE: fix all paras expect of skal, and n1 => we can take
   #                              the reinfrocement to zero and fit a beta dstribution.
   # unrestrict.theta   == TRUE:  allow theta to vary.
   #                              (if FALSE: we can fix theta=0.5, and in this
   #                               try to find ut where the reinforceent takes place)

   para.last <- para.init; 
   lll.last  <- lll(para.ref);
   last.skal = -1; no.skal.const = 0;last.lll=-1e10;
   fertig = FALSE;

   i = 0;
   if(!exists("maxIter")) {maxIter=1000;}
   while ((i<maxIter)&(fertig==FALSE)){
      i = i+1;
      cat(para.last,"\n");
      lll.x = lll.last;
      para.last <<- para.last;
      res1 = optimize(search.p1, interval=c(0,1), maximum=TRUE);
      para.loc1 = para.last; para.loc1[1]=res1$maximum;
      lll.lok  = lll(para.loc1);
      if (lll.lok>lll.last){
         para.last <- para.loc1;
         lll.last  <- lll.lok;
      }

      if (unrestricted.model){
         para.last <<- para.last;
         res2 = optimize(search.p2, interval=c(0,1), maximum=TRUE);
         para.loc2 = para.last; para.loc2[2]=res2$maximum;
         lll.lok  = lll(para.loc2);
         if (lll.lok>lll.last){
            para.last <- para.loc2;
            lll.last  <- lll.lok;
         }

         if (unrestrict.theta) {
            para.last <<- para.last;
            res3 = optimize(search.p3, interval=c(0,1), maximum=TRUE);
            para.loc3 = para.last; para.loc3[3]=res3$maximum;
            lll.lok  = lll(para.loc3);
            if (lll.lok>lll.last){
               para.last <- para.loc3;
               lll.last  <- lll.lok;
            }
         }
      }

      para.last <<- para.last;
      res5 = optimize(search.p5, interval=c(0,50000), maximum=TRUE);
      para.loc5 = para.last; para.loc5[5]=res5$maximum;
      lll.lok  = lll(para.loc5);
      if (lll.lok>lll.last){
         para.last <- para.loc5;
         lll.last  <- lll.lok;
      } 

      lll.1     = lll(para.last);
      Delta = 0.01;
      if (para.last[4]>10)  Delta=0.1;
      if (para.last[4]>50)  Delta=0.5;
      if (para.last[4]>100) Delta=1;

      lll.p2    = lll(para.last+c(0,0, 0, Delta,0));
      lll.m2    = lll(para.last+c(0,0,0, -Delta,0));
      if (lll.p2>lll.1){
         paral.loc3 = para.last+c(0,0,0, Delta,0);
         para.last  = para.last+c(0,0,0, Delta,0);
         last.lll = lll.p2;
      } else {
         if (lll.m2>lll.1){
            paral.loc3 = para.last+c(0,0,0,-Delta,0);
            para.last  <- para.last+c(0,0,0,-Delta,0);
            last.lll   <- lll.m2
         } else {
            paral.loc3 = para.last+c(0,0,0,0,0);
            para.last  <- para.last+c(0,0,0,0,0);
            last.lll   <- lll.1
         }
      }
      if (last.skal!=para.last[3]){
         last.skal = para.last[3];
         no.skal.const = 0;
         lll.x =last.lll;
      } else {
         no.skal.const = no.skal.const+1;
         loc.lll = lll(para.last);
         if (last.lll>lll.x+1e-6) {
            no.skal.const = 0;
         }
         if(no.skal.const>100) fertig=TRUE;
      }

     para.last <<- para.last;
     lll.last  <<- last.lll;
     cat(i," ", lll(para.last), "\n");     
     curve(g(x), add=TRUE, col="blue");
   }
}
