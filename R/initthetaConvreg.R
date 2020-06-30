initthetaConvreg = function(
   y,
   dist1,
   dist2,
   formula.mu1, formula.sigma1, formula.mu2, formula.sigma2,
   fixed,
   data,
   theta.name,
   scale,
   scaleInit,
   debug
){

  if(scale){y = y / scaleInit}

   y.r = round((y - min(y)))

if(debug){browser()}

   if(dist1 == "Multinom"){

      n.mult = fixed$value[which(fixed$name == "sigma 1: (Intercept)")]

      y.r = apply(y.r,1,function(x) min(x,n.mult))
   }

   if(dist1=="Nbinom"){model.1 = glm.nb(eval(parse(text=paste("y.r~",formula.mu1[2],sep=""))),data=data)}
   if(dist1=="Pois"){model.1 = glm(eval(parse(text=paste("y.r~",formula.mu1[2],sep=""))),data=data, family = "poisson")}
   if(dist1=="CoMPoisson"){model.1 = glm.cmp(eval(parse(text=paste("y.r~",formula.mu1[2],sep=""))),data=data)}
   if(dist1=="Binom"){model.1 = glm(eval(parse(text=paste("y.r/",fixed$value[which(fixed$name == "sigma 1: (Intercept)")],"~",formula.mu1[2],sep=""))),data=data, family = "quasibinomial")}
   if(dist1=="Multinom"){model.1 = multinom(eval(parse(text=paste("y.r~",formula.mu1[2],sep=""))),data=data)}
   if(dist1=="ZIP"){model.1 = zeroinfl(eval(parse(text=paste("y.r~",formula.mu1[2],"|",formula.sigma1[2],sep=""))),data=data)}
   if(dist1=="HP"){model.1 = hurdle(eval(parse(text=paste("y.r~",formula.mu1[2],"|",formula.sigma1[2],sep=""))),data=data)}


   if(dist1=="Nbinom"){theta.init.mu1 = model.1$coefficients;theta.init.sigma1 = model.1$theta}
   if(dist1=="Pois"){theta.init.mu1 = model.1$coefficients;theta.init.sigma1 = NA}
   if(dist1=="CoMPoisson"){theta.init.mu1 = coef(model.1)$beta;theta.init.sigma1 = coef(model.1)$zeta;}
   if(dist1=="Binom"){theta.init.mu1 = coef(model.1);theta.init.sigma1 = fixed$value[which(fixed$name == "sigma 1: (Intercept)")];}
   if(dist1=="Multinom"){theta.init.mu1 = coef(model.1);theta.init.sigma1 = NA;}
   if(dist1=="ZIP"){theta.init.mu1 = model.1$coefficients$count;theta.init.sigma1 = model.1$coefficients$zero}
   if(dist1=="HP"){theta.init.mu1 = model.1$coefficients$count;theta.init.sigma1 = model.1$coefficients$zero}

   if(dist1 == "Binom" & length(all.vars(formula.mu1)) == 0){names(theta.init.mu1) = "(Intercept)"}

   if(dist1 != "Multinom"){
      #names(theta.init.1)[which(names(theta.init.1) == "(Intercept)")] = "Intercept.mu1"
      #names(theta.init.1b)[which(names(theta.init.1b) == "(Intercept)")] = "Intercept.sigma1"
      names(theta.init.mu1) = paste("mu 1: ",names(theta.init.mu1),sep="")
      if(is.null(names(theta.init.sigma1))){names(theta.init.sigma1) = "(Intercept)"}
      names(theta.init.sigma1) = paste("sigma 1: ",names(theta.init.sigma1),sep="")

      theta.init.mu1[which(abs(theta.init.mu1) > 20)] = 0
      theta.init.sigma1[which(abs(theta.init.sigma1) > 20)] = 0

      idx.th.in.mu1 = na.omit(match(names(theta.init.mu1), theta.name))

      idx.th.in.sigma1 = na.omit(match(names(theta.init.sigma1), theta.name))
   }
   if(dist1 == "Multinom"){
      #theta.init.1 = rbind(rep(0,ncol(theta.init.1)),theta.init.1)

      miss = setdiff(1:n.mult,rownames(theta.init.mu1))

      theta.init.1.f = matrix(0,nrow = n.mult, ncol = ncol(theta.init.mu1))
      colnames(theta.init.1.f)  = colnames(theta.init.mu1)
      rownames(theta.init.1.f)  = 1:n.mult

      theta.init.1.f[rownames(theta.init.mu1),] = theta.init.mu1

      theta.init.1.f[miss,1] = -100

      theta.init.mu1 = theta.init.1.f

      idx.0 = which(colSums(theta.init.mu1^2)==0)

      if(length(idx.0)>0){theta.init.mu1 = theta.init.mu1[,-idx.0]}

      #colnames(theta.init.1)[which(colnames(theta.init.1) == "(Intercept)")] = "Intercept.mu1"

      colnames(theta.init.mu1) = paste("mu 1: ",colnames(theta.init.mu1),sep="")

      idx.th.in.mu1 = match(colnames(theta.init.mu1), theta.name)
   }

   y.2 = y
   if(length(all.vars(formula.mu2))==0 & dist2 == "Gauss"){y.2 = y- round(y - min(y))}

   if(dist2=="Gauss"){model.2 = lm(eval(parse(text=paste("y.2~",formula.mu2[2],"",sep=""))),data)}
   if(dist2=="Lnorm"){model.2 = lm(eval(parse(text=paste("log(y.2)~",formula.mu2[2],"",sep=""))),data)}
   if(dist2=="Gamma"){model.2 = glm(eval(parse(text=paste("y.2~",formula.mu2[2],"",sep=""))),data=data, family=gamma)}

   theta.init.mu2 = model.2$coefficients
   theta.init.sigma2 = rep(0.1, length(grep("sigma 2", theta.name)))

   names(theta.init.mu2) = paste("mu 2: ",names(theta.init.mu2),sep="")
   names(theta.init.sigma2) = theta.name[grep("sigma 2", theta.name)]

   idx.th.in.mu2 = match(names(theta.init.mu2), theta.name)
   idx.th.in.sigma2 = match(names(theta.init.sigma2), theta.name)
   RET = list(theta.mu1 = theta.init.mu1, theta.sigma1 = theta.init.sigma1,
              theta.mu2 = theta.init.mu2, theta.sigma2 = theta.init.sigma2,
              idx1 = idx.th.in.mu1, idx1b = idx.th.in.sigma1,
              idx2 = idx.th.in.mu2, idx2b = idx.th.in.sigma2)


   return(RET)
}
