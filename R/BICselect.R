#' Convolutive regression best model selection BIC
#'
#' Convolutive regression best model selection BIC
#' @param df dataset
#' @param formula.resp response formula
#' @param idx.pred indices of predictors
#' @param dist1 chain of character to idnetify the distribution of variable1
#' @param dist2 chain of character to idnetify the distribution of variable2
#' @param quiet logical. Displaying the progress.
#' @param formulas "main" or "all". All includes the variances.
#' @param parallel logical. Should the function run on parallel threads.
#' @param interactions logical. Should the formulas include interactions.
#' @param debug logical. Working option to debug. Dvp.
#' @param ... Additional arguments fed into convreg().
#' @keywords BIC
#' @export
#' @examples
#' set.seed(123)
#' e=rnorm(n=50,mean=0,sd=0.05)
#' x = rnorm(50,0,0.5)
#' k= rnbinom(50,exp(x),0.5)
#' y= data.frame(obs= k + e +x, f1 = x, f2 = k, f3 = 2*rnorm(n=50,mean=0,sd=0.05))
#' form.r = ~obs
#' df = y
#' idx.pred = 2
#' dist1 = "Nbinom"
#' dist2 = "Gauss"
#' res= BICselect(df = df, formula.resp = form.r, idx.pred = idx.pred, dist1 = dist1, dist2 = dist2)
#' head(res,10)
BICselect = function(df, formula.resp, idx.pred, dist1, dist2, quiet = TRUE, formulas = "main", parallel = FALSE, interactions = FALSE, debug = FALSE, ...){

   idx.resp = match(all.vars(formula.resp), names(df))

   df = df[,c(idx.pred, idx.resp)]

   n = nrow(df)
   p = ncol(df) - 1

   list.comb.nb = lapply(0:(p), function(x) combn(1:p,x, simplify = F))
   list.comb.ga = lapply(0:(p), function(x) combn(1:p,x, simplify = F))



   if(parallel == FALSE){

   if(formulas == "main"){

   summaryStat =NULL

   mat.comb.nb = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)
   mat.comb.ga = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)

   ncol.nb = ncol(mat.comb.nb)
   ncol.ga = ncol(mat.comb.ga)

   #full.comb  = rbind(mat.comb.nb[,rep(1:ncol.ga, each = ncol.ga)], mat.comb.ga[,rep(1:ncol.nb,ncol.nb)])
   mc.nb.rep = matrix(mat.comb.nb[,rep(1:ncol.ga, each = ncol.ga)], nrow = p)
   mc.ga.rep = matrix(mat.comb.ga[,rep(1:ncol.nb,ncol.nb)], nrow = p)

   idx.keep= which(colSums(mc.nb.rep*mc.ga.rep)==0)

   mc.nb.rep = matrix(mc.nb.rep[,idx.keep], nrow = p)
   mc.ga.rep = matrix(mc.ga.rep[,idx.keep], nrow = p)

   message(sprintf('Number of combinations: %s', ncol(mc.nb.rep)))

   for(i in 1:ncol(mc.nb.rep)){
      if(!quiet){
         if((i / round(.05 * ncol(mc.nb.rep))) == round(i/round(.05 * ncol(mc.nb.rep)))){
            message(sprintf('%s / %s achieved.', i, ncol(mc.nb.rep) ))
         }
      }
      idx.nb = which(mc.nb.rep[,i]==1)
      idx.ga = which(mc.ga.rep[,i]==1)
      formula.mu1 = paste0("~ ",paste(names(df)[idx.nb], sep = "", collapse = "+"))
      formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "+"))

      if(interactions){
         formula.mu1 = paste0("~ ",paste(names(df)[idx.nb], sep = "", collapse = "*"))
         formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "*"))
      }

      if(length(idx.nb)==0){formula.mu1 = "~1"}
      if(length(idx.ga)==0){formula.mu2 = "~1"}
      formula.mu1 = as.formula(formula.mu1)
      formula.mu2 = as.formula(formula.mu2)
      res.t = convreg(formula.resp,
                      formula.mu1 = formula.mu1 ,
                      formula.mu2 = formula.mu2 ,
                      data = df, family = "iid",
                      dist1 = dist1, dist2 = dist2, ...)
      BIC = -2 * res.t$loglik + (length(idx.nb) + length(idx.ga) + 2) * log(n)

      summaryStat = rbind(summaryStat, cbind(paste(as.character(formula.mu1), collapse = ""), paste(as.character(formula.mu2), collapse = ""), BIC))

   }


   }
   if(formulas == "all"){

      summaryStat =NULL

      mat.comb.nb = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)
      mat.comb.ga = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)

      ncol.nb = ncol(mat.comb.nb)
      ncol.ga = ncol(mat.comb.ga)

      #full.comb  = rbind(mat.comb.nb[,rep(1:ncol.ga, each = ncol.ga)], mat.comb.ga[,rep(1:ncol.nb,ncol.nb)])

      exp = expand.grid(1:ncol.nb, 1:ncol.nb, 1:ncol.ga)

      mc.nb.rep1 = matrix(mat.comb.nb[,exp[,1]], nrow = p)
      mc.ga.rep =  matrix(mat.comb.ga[,exp[,3]], nrow = p)
      mc.nb.rep2 = matrix(mat.comb.nb[,exp[,2]], nrow = p)

      idx.keep= which(colSums(mc.nb.rep1*mc.ga.rep*mc.nb.rep2)==0)

      mc.nb.rep1 = matrix(mc.nb.rep1[,idx.keep], nrow = p)
      mc.nb.rep2 = matrix(mc.nb.rep2[,idx.keep], nrow = p)
      mc.ga.rep = matrix(mc.ga.rep[,idx.keep], nrow = p)

      message(sprintf('Number of combinations: %s', ncol(mc.nb.rep1)))

      for(i in 1:ncol(mc.nb.rep1)){
         if(!quiet){
            if((i / round(.05 * ncol(mc.nb.rep1))) == round(i/round(.05 * ncol(mc.nb.rep1)))){
               message(sprintf('%s / %s achieved.', i, ncol(mc.nb.rep1) ))
            }
         }
         idx.nb1 = which(mc.nb.rep1[,i]==1)
         idx.nb2 = which(mc.nb.rep2[,i]==1)
         idx.ga = which(mc.ga.rep[,i]==1)
         formula.mu1 = paste0("~ ",paste(names(df)[idx.nb1], sep = "", collapse = "+"))
         formula.sigma1 = paste0("~ ",paste(names(df)[idx.nb2], sep = "", collapse = "+"))
         formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "+"))

         if(interactions){
            formula.mu1 = paste0("~ ",paste(names(df)[idx.nb1], sep = "", collapse = "*"))
            formula.sigma1 = paste0("~ ",paste(names(df)[idx.nb2], sep = "", collapse = "*"))
            formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "*"))
         }


         if(length(idx.nb1)==0){formula.mu1 = "~1"}
         if(length(idx.nb2)==0){formula.sigma1 = "~1"}
         if(length(idx.ga)==0){formula.mu2 = "~1"}

         formula.mu1 = as.formula(formula.mu1)
         formula.sigma1 = as.formula(formula.sigma1)
         formula.mu2 = as.formula(formula.mu2)

         if(debug){browser()}

         res.t = convreg(formula.resp,
                         formula.mu1 = formula.mu1 ,
                         formula.sigma1 = formula.sigma1 ,
                         formula.mu2 = formula.mu2 ,
                         data = df, family = "iid",
                         dist1 = dist1, dist2 = dist2, ...)
         BIC = -2 * res.t$loglik + (length(idx.nb1) + length(idx.nb2) + length(idx.ga) + 1) * log(n)

         summaryStat = rbind(summaryStat, cbind(paste(as.character(formula.mu1), collapse = ""),
                                                paste(as.character(formula.sigma1), collapse = ""),
                                                paste(as.character(formula.mu2), collapse = ""), BIC))

      }

   }

   }

   if(parallel == TRUE){

      # Calculate the number of cores
      no_cores <- detectCores() - 1

      registerDoMC(no_cores)

      mat.comb.nb = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)
      mat.comb.ga = matrix(unlist(sapply(0:p, function(x) combn(1:p, x, tabulate, nbins = p))), nrow = p)

      ncol.nb = ncol(mat.comb.nb)
      ncol.ga = ncol(mat.comb.ga)

      #full.comb  = rbind(mat.comb.nb[,rep(1:ncol.ga, each = ncol.ga)], mat.comb.ga[,rep(1:ncol.nb,ncol.nb)])

      if(formulas == "main"){

         exp = expand.grid(1:ncol.nb, 1:ncol.ga)

         mc.nb.rep = mat.comb.nb[,exp[,1]]
         mc.ga.rep = mat.comb.ga[,exp[,2]]

         idx.keep= which(colSums(mc.nb.rep*mc.ga.rep)==0)

         mc.nb.rep = mc.nb.rep[,idx.keep]
         mc.ga.rep = mc.ga.rep[,idx.keep]


         summaryStat <-foreach(i=1:ncol(mc.nb.rep),.combine=rbind) %do% {

            idx.nb = which(mc.nb.rep[,i]==1)
            idx.ga = which(mc.ga.rep[,i]==1)

            formula.mu1 = paste0("~ ",paste(names(df)[idx.nb], sep = "", collapse = "+"))
            formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "+"))

            if(length(idx.nb)==0){formula.mu1 = "~1"}
            if(length(idx.ga)==0){formula.mu2 = "~1"}

            formula.mu1 = as.formula(formula.mu1)
            formula.mu2 = as.formula(formula.mu2)

            res.t = convreg(formula.resp,
                            formula.mu1 = formula.mu1 ,
                            formula.mu2 = formula.mu2 ,
                            data = df, family = "iid",
                            dist1 = dist1, dist2 = dist2, ...)

            BIC = -2 * res.t$loglik + (length(idx.nb) + length(idx.ga) + 2) * log(n)

            summaryStat = cbind(paste(as.character(formula.mu1), collapse = ""),
                                paste(as.character(formula.mu2), collapse = ""), BIC)
         }

         stopImplicitCluster()
      }

      if(formulas == "all"){

      exp = expand.grid(1:ncol.nb, 1:ncol.nb, 1:ncol.ga)

      mc.nb.rep1 = mat.comb.nb[,exp[,1]]
      mc.ga.rep = mat.comb.ga[,exp[,3]]
      mc.nb.rep2 = mat.comb.nb[,exp[,2]]

      idx.keep= which(colSums(mc.nb.rep1*mc.ga.rep*mc.nb.rep2)==0)

      mc.nb.rep1 = mc.nb.rep1[,idx.keep]
      mc.nb.rep2 = mc.nb.rep2[,idx.keep]
      mc.ga.rep = mc.ga.rep[,idx.keep]


      summaryStat <-foreach(i=1:ncol(mc.nb.rep1),.combine=rbind) %do% {

         idx.nb1 = which(mc.nb.rep1[,i]==1)
         idx.nb2 = which(mc.nb.rep2[,i]==1)
         idx.ga = which(mc.ga.rep[,i]==1)
         formula.mu1 = paste0("~ ",paste(names(df)[idx.nb1], sep = "", collapse = "+"))
         formula.sigma1 = paste0("~ ",paste(names(df)[idx.nb2], sep = "", collapse = "+"))
         formula.mu2 = paste0("~ ",paste(names(df)[idx.ga], sep = "", collapse = "+"))
         if(length(idx.nb1)==0){formula.mu1 = "~1"}
         if(length(idx.nb2)==0){formula.sigma1 = "~1"}
         if(length(idx.ga)==0){formula.mu2 = "~1"}

         formula.mu1 = as.formula(formula.mu1)
         formula.sigma1 = as.formula(formula.sigma1)
         formula.mu2 = as.formula(formula.mu2)

         res.t = convreg(formula.resp,
                         formula.mu1 = formula.mu1 ,
                         formula.sigma1 = formula.sigma1 ,
                         formula.mu2 = formula.mu2 ,
                         data = df, family = "iid",
                         dist1 = dist1, dist2 = dist2, ...)
         BIC = -2 * res.t$loglik + (length(idx.nb1) + length(idx.nb2) + length(idx.ga) + 1) * log(n)

         summaryStat = cbind(paste(as.character(formula.mu1), collapse = ""),
                                                paste(as.character(formula.sigma1), collapse = ""),
                                                paste(as.character(formula.mu2), collapse = ""), BIC)
      }

      stopImplicitCluster()
      }


   }

   idx.ord = sort(as.numeric(summaryStat[,ncol(summaryStat)]), index.return = T)$ix
   if(formulas == "main"){RET = data.frame(formula1 = summaryStat[idx.ord,1], formula2 = summaryStat[idx.ord,2] ,BIC = as.numeric(summaryStat[idx.ord,3]))}
   if(formulas == "all"){RET = data.frame(formula1 = summaryStat[idx.ord,1], formula11 = summaryStat[idx.ord,2] , formula2 = summaryStat[idx.ord,3] ,BIC = as.numeric(summaryStat[idx.ord,4]))}
   return(RET)
}
