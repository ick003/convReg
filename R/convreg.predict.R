predict.convreg <- function(obj, newdata=NULL,pred.type = "point", delta= 1, xlim=c(0,25), 
                            with.confidence.interval=FALSE){
   
   if(!(pred.type %in% c("point","hist", "contr2"))){stop("Not the right pred.type call")}
   
   if(is.null(newdata) & pred.type == "point"){
      pred = obj$xdata[,c("y","y.fitted","fitted.up","fitted.low")]
      sample = NULL
   }
   if(pred.type == "hist"){
      
      if(is.null(newdata)){newdata = obj$data}
      
      dist1= unlist(strsplit(obj$distname,"/"))[1]
      dist2= unlist(strsplit(obj$distname,"/"))[2]
      
      X.mu1 = model.matrix(obj$formulas$mu1,model.frame(obj$formulas$mu1,data = newdata,na.action=NULL))
      colnames(X.mu1)[1] <- "Intercept.mu1"
      
      X.sigma1 = model.matrix(obj$formulas$sigma1,model.frame(obj$formulas$sigma1,data = newdata,na.action=NULL))
      colnames(X.sigma1)[1] <- "Intercept.sigma1"
      
      X.mu2 = model.matrix(obj$formulas$mu2,model.frame(obj$formulas$mu2,data = newdata,na.action=NULL))
      colnames(X.mu2)[1] <- "Intercept.mu2"
      
      X.sigma2 = model.matrix(obj$formulas$sigma2,model.frame(obj$formulas$sigma2,data = newdata,na.action=NULL))
      colnames(X.sigma2)[1] <- "Intercept.sigma2"
      
      idx.par2 = match(paste("sigma 1: ",colnames(X.sigma1),sep=""), rownames(obj$estimation))
      sigma1 = obj$transforms[[2]](X.sigma1 %*% obj$estimation$Estimate[idx.par2])
      n.mult = max(sigma1)
      
      if(dist1 != "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1), sep = "")))]}}
      #if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1)," (p=",1:n.mult,")", sep = "")))]}}
      
      if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,which(unlist(lapply(colnames(X.mu1), function(x) !is.null(intersect(rownames(obj$estimation),paste("mu 1: ",x," (p=",1:n.mult,")", sep = ""))))))]}}      
      if(ncol(X.sigma1)>1){X.sigma1 = X.sigma1[,na.omit(match( rownames(obj$estimation),paste("sigma 1: ",colnames(X.sigma1), sep = "")))]}
      
      if(ncol(X.mu2)>1){X.mu2 = X.mu2[,na.omit(match( rownames(obj$estimation),paste("mu 2: ",colnames(X.mu2), sep = "")))]}
      
      if(ncol(X.sigma2)>1){X.sigma2 = X.sigma2[,na.omit(match( rownames(obj$estimation),paste("sigma 2: ",colnames(X.sigma2), sep = "")))]}
      
      if(dist1 == "Multinom"){idx.par1 = match(paste("mu 1: ",rep(colnames(X.mu1),n.mult)," (p=",rep(1:n.mult,each=ncol(X.mu1)),")", sep = ""), rownames(obj$estimation))}
      if(dist1 != "Multinom"){idx.par1 = match(paste("mu 1: ",colnames(X.mu1),sep=""), rownames(obj$estimation))}
      
      idx.par3 = match(paste("mu 2: ",colnames(X.mu2),sep=""), rownames(obj$estimation))
      idx.par4 = match(paste("sigma 2: ",colnames(X.sigma2),sep=""), rownames(obj$estimation))
      
      if(dist1 != "Multinom"){mu1 = matrix(obj$transforms[[1]](X.mu1 %*% obj$estimation$Estimate[idx.par1]),ncol=1)}
      if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.par1],nrow = n.mult)))) / rowSums(exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.par1],nrow = n.mult)))))}
      
      
      mu2 = obj$transforms[[3]](X.mu2 %*% obj$estimation$Estimate[idx.par3])
      sigma2 = as.matrix(obj$transforms[[4]](X.sigma2 %*% obj$estimation$Estimate[idx.par4]))
        
      R.Cmd <- sprintf("sample = matrix(replicate(5000/nrow(mu1), r%s%s(N=nrow(mu1), mu1,sigma1, mu2,sigma2)),ncol=5000/nrow(mu1),byrow=F)", dist1, dist2)
      
      eval(parse(text=R.Cmd))
      
       seq.x2 = seq(xlim[1],xlim[2],delta)
#       
#       if(dist1 != "Multinom"){R.Cmd <- sprintf("hist = apply(cbind(mu1,sigma1,mu2,sigma2),1, function(x) d%s%s(seq.x2, x[1],x[2], x[3],x[4], log=F))", dist1, dist2)}
#       if(dist1 == "Multinom"){R.Cmd <- sprintf("hist = apply(cbind(mu1,sigma1,mu2,sigma2),1, function(x) d%s%s(seq.x2, x[1:ncol(mu1)],x[(ncol(mu1)+ 1):(ncol(mu1)+ncol(sigma1))], x[(ncol(mu1)+ncol(sigma1)+ 1):(ncol(mu1)+ncol(sigma1)+ncol(mu2))],x[(ncol(mu1)+ncol(sigma1)+ncol(mu2)+ 1):(ncol(mu1)+ncol(sigma1)+ncol(mu2)+ncol(sigma2))], log=F))", dist1, dist2)}
#       
#       eval(parse(text=R.Cmd))

      
      
      seq.x3 = seq(xlim[1]-delta/2,xlim[2]+delta/2,delta)
      q = t(apply(apply(sample, 2 , function(x){hist(x,breaks = seq.x3,plot=F)$density}),1,function(y){quantile(y,probs = c(0.025,0.975))}))
      
      hist.f = apply(apply(sample, 2 , function(x){hist(x,breaks = seq.x3,plot=F)$density}),1,function(y){mean(y)})

      
#      for(i in 1:length(seq.x2)){
#         
#         q = rbind(q,quantile(hist[i,], probs = c(0.025,0.975)))
#         
#      }  
#       
#       
      pred = data.frame(x = seq.x2,density = hist.f)
      sample = q
   }
   if(pred.type == "point" & !is.null(newdata)){
      
      dist1= unlist(strsplit(obj$distname,"/"))[1]
      dist2= unlist(strsplit(obj$distname,"/"))[2]
      
      X.mu1 = model.matrix(obj$formulas$mu1,model.frame(obj$formulas$mu1,data = newdata,na.action=NULL))
      colnames(X.mu1)[1] <- "Intercept.mu1"
      
      X.sigma1 = model.matrix(obj$formulas$sigma1,model.frame(obj$formulas$sigma1,data = newdata,na.action=NULL))
      colnames(X.sigma1)[1] <- "Intercept.sigma1"
      
      X.mu2 = model.matrix(obj$formulas$mu2,model.frame(obj$formulas$mu2,data = newdata,na.action=NULL))
      colnames(X.mu2)[1] <- "Intercept.mu2"
      
      X.sigma2 = model.matrix(obj$formulas$sigma2,model.frame(obj$formulas$sigma2,data = newdata,na.action=NULL))
      colnames(X.sigma2)[1] <- "Intercept.sigma2"
      
      idx.par2 = match(paste("sigma 1: ",colnames(X.sigma1),sep=""), rownames(obj$estimation))
      sigma1 = obj$transforms[[2]](X.sigma1 %*% obj$estimation$Estimate[idx.par2])
      n.mult = max(sigma1)
      
      if(dist1 != "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1), sep = "")))]}}
      #if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1)," (p=",1:n.mult,")", sep = "")))]}}
      
      if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,which(unlist(lapply(colnames(X.mu1), function(x) !is.null(intersect(rownames(obj$estimation),paste("mu 1: ",x," (p=",1:n.mult,")", sep = ""))))))]}}      
      if(ncol(X.sigma1)>1){X.sigma1 = X.sigma1[,na.omit(match( rownames(obj$estimation),paste("sigma 1: ",colnames(X.sigma1), sep = "")))]}
      
      if(ncol(X.mu2)>1){X.mu2 = X.mu2[,na.omit(match( rownames(obj$estimation),paste("mu 2: ",colnames(X.mu2), sep = "")))]}
      
      if(ncol(X.sigma2)>1){X.sigma2 = X.sigma2[,na.omit(match( rownames(obj$estimation),paste("sigma 2: ",colnames(X.sigma2), sep = "")))]}
      
      if(dist1 == "Multinom"){idx.par1 = match(paste("mu 1: ",rep(colnames(X.mu1),n.mult)," (p=",rep(1:n.mult,each=ncol(X.mu1)),")", sep = ""), rownames(obj$estimation))}
      if(dist1 != "Multinom"){idx.par1 = match(paste("mu 1: ",colnames(X.mu1),sep=""), rownames(obj$estimation))}
      
      idx.par3 = match(paste("mu 2: ",colnames(X.mu2),sep=""), rownames(obj$estimation))
      idx.par4 = match(paste("sigma 2: ",colnames(X.sigma2),sep=""), rownames(obj$estimation))
      
      if(dist1 != "Multinom"){mu1 = matrix(obj$transforms[[1]](X.mu1 %*% obj$estimation$Estimate[idx.par1]),ncol=1)}
      if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.par1],nrow = n.mult)))) / rowSums(exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.par1],nrow = n.mult)))))}
      
      mu2 = obj$transforms[[3]](X.mu2 %*% obj$estimation$Estimate[idx.par3])
      sigma2 = as.matrix(obj$transforms[[4]](X.sigma2 %*% obj$estimation$Estimate[idx.par4]))
      
      
      
      # FITTED VALUE composed of 2 model
      
      if(dist1 == "Nbinom" | dist1 == "Pois"){exp1 = mu1}
      if(dist1 == "Binom"){exp1 = mu1*sigma1}
      if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
      
      if(dist2 == "Gauss"){exp2 = mu2}
      if(dist2 == "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}
      
      y.pred = exp1 + exp2
      
#       y.pred = rowMeans(sample, na.rm=T)
      
      
      R.Cmd <- sprintf("y.pred.ci = q%s%s(c(0.025, 0.975), mu1,sigma1, mu2,sigma2)", dist1, dist2)
      
      eval(parse(text=R.Cmd))
      
      pred = data.frame(y.fitted = y.pred, fitted.low = y.pred.ci[,1],fitted.up = y.pred.ci[,2])
      sample = NULL
   }
   
   if(pred.type == "contr" & !is.null(newdata)){
      
      dist1= unlist(strsplit(obj$distname,"/"))[1]
      dist2= unlist(strsplit(obj$distname,"/"))[2]
      
      name.contr = names(which(apply(as.matrix(newdata),2,function(x) (x[1] != x[2]))))
      
      idx.contr = (1:nrow(obj$estimation))[-c(grep(name.contr, rownames(obj$estimation)), grep("Intercept",rownames(obj$estimation)))]
      
      estim = obj$estimation$Estimate
      
      #estim[idx.contr] = 0
         
      X.mu1 = model.matrix(obj$formulas$mu1,model.frame(obj$formulas$mu1,data = newdata,na.action=NULL))
      colnames(X.mu1)[1] <- "Intercept.mu1"
      
      X.sigma1 = model.matrix(obj$formulas$sigma1,model.frame(obj$formulas$sigma1,data = newdata,na.action=NULL))
      colnames(X.sigma1)[1] <- "Intercept.sigma1"
      
      X.mu2 = model.matrix(obj$formulas$mu2,model.frame(obj$formulas$mu2,data = newdata,na.action=NULL))
      colnames(X.mu2)[1] <- "Intercept.mu2"
      
      X.sigma2 = model.matrix(obj$formulas$sigma2,model.frame(obj$formulas$sigma2,data = newdata,na.action=NULL))
      colnames(X.sigma2)[1] <- "Intercept.sigma2"
      
      idx.par2 = match(paste("sigma 1: ",colnames(X.sigma1),sep=""), rownames(obj$estimation))
      sigma1 = obj$transforms[[2]](X.sigma1 %*% estim[idx.par2])
      n.mult = max(sigma1)
      
      if(dist1 != "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1), sep = "")))]}}
      #if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1)," (p=",1:n.mult,")", sep = "")))]}}
      
      if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,which(unlist(lapply(colnames(X.mu1), function(x) !is.null(intersect(rownames(obj$estimation),paste("mu 1: ",x," (p=",1:n.mult,")", sep = ""))))))]}}      
      if(ncol(X.sigma1)>1){X.sigma1 = X.sigma1[,na.omit(match( rownames(obj$estimation),paste("sigma 1: ",colnames(X.sigma1), sep = "")))]}
      
      if(ncol(X.mu2)>1){X.mu2 = X.mu2[,na.omit(match( rownames(obj$estimation),paste("mu 2: ",colnames(X.mu2), sep = "")))]}
      
      if(ncol(X.sigma2)>1){X.sigma2 = X.sigma2[,na.omit(match( rownames(obj$estimation),paste("sigma 2: ",colnames(X.sigma2), sep = "")))]}
      
      if(dist1 == "Multinom"){idx.par1 = match(paste("mu 1: ",rep(colnames(X.mu1),n.mult)," (p=",rep(1:n.mult,each=ncol(X.mu1)),")", sep = ""), rownames(obj$estimation))}
      if(dist1 != "Multinom"){idx.par1 = match(paste("mu 1: ",colnames(X.mu1),sep=""), rownames(obj$estimation))}
      
      idx.par3 = match(paste("mu 2: ",colnames(X.mu2),sep=""), rownames(obj$estimation))
      idx.par4 = match(paste("sigma 2: ",colnames(X.sigma2),sep=""), rownames(obj$estimation))
      
      if(dist1 != "Multinom"){mu1 = matrix(obj$transforms[[1]](X.mu1 %*% estim[idx.par1]),ncol=1)}
      if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(estim[idx.par1],nrow = n.mult)))) / rowSums(exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(estim[idx.par1],nrow = n.mult)))))}
      
      mu2 = obj$transforms[[3]](X.mu2 %*% estim[idx.par3])
      sigma2 = as.matrix(obj$transforms[[4]](X.sigma2 %*% estim[idx.par4]))
      
      if(dist1 == "Nbinom" | dist1 == "Pois"){exp1 = mu1}
      if(dist1 == "Binom"){exp1 = mu1*sigma1}
      if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
      if(dist1 == "Multinom"){exp1 = as.numeric(apply(mu1,1,which.max))-1}
      
      if(dist2 == "Gauss"){exp2 = mu2}
      if(dist2 == "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}
      
      y.pred = exp1 + exp2

      R.Cmd <- sprintf("y.pred.ci = q%s%s(c(0.025, 0.975), mu1,sigma1, mu2,sigma2)", dist1, dist2)
      
      eval(parse(text=R.Cmd))
      
      pred = data.frame(y.fitted = y.pred, fitted.low = y.pred.ci[,1],fitted.up = y.pred.ci[,2])
      sample = NULL
   }
   
   if(pred.type == "contr2" & !is.null(newdata)){
   
   dist1= unlist(strsplit(obj$distname,"/"))[1]
   dist2= unlist(strsplit(obj$distname,"/"))[2]
   
   estim = obj$estimation$Estimate
   
   X.mu1 = model.matrix(obj$formulas$mu1,model.frame(obj$formulas$mu1,data = newdata,na.action=NULL))
   colnames(X.mu1)[1] <- "(Intercept)"
   
   X.sigma1 = model.matrix(obj$formulas$sigma1,model.frame(obj$formulas$sigma1,data = newdata,na.action=NULL))
   colnames(X.sigma1)[1] <- "(Intercept)"
   
   X.mu2 = model.matrix(obj$formulas$mu2,model.frame(obj$formulas$mu2,data = newdata,na.action=NULL))
   colnames(X.mu2)[1] <- "(Intercept)"
   
   X.sigma2 = model.matrix(obj$formulas$sigma2,model.frame(obj$formulas$sigma2,data = newdata,na.action=NULL))
   colnames(X.sigma2)[1] <- "(Intercept)"
   
   idx.par2 = match(paste("sigma 1: ",colnames(X.sigma1),sep=""), rownames(obj$estimation))
   sigma1 = obj$transforms[[2]](X.sigma1 %*% estim[idx.par2])
   n.mult = max(sigma1)
   
   if(dist1 != "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1), sep = "")))]}}
   #if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,na.omit(match( rownames(obj$estimation),paste("mu 1: ",colnames(X.mu1)," (p=",1:n.mult,")", sep = "")))]}}
   
   if(dist1 == "Multinom"){if(ncol(X.mu1)>1){X.mu1 = X.mu1[,which(unlist(lapply(colnames(X.mu1), function(x) !is.null(intersect(rownames(obj$estimation),paste("mu 1: ",x," (p=",1:n.mult,")", sep = ""))))))]}}      
   if(ncol(X.sigma1)>1){X.sigma1 = X.sigma1[,na.omit(match( rownames(obj$estimation),paste("sigma 1: ",colnames(X.sigma1), sep = "")))]}
   
   if(ncol(X.mu2)>1){X.mu2 = X.mu2[,na.omit(match( rownames(obj$estimation),paste("mu 2: ",colnames(X.mu2), sep = "")))]}
   
   if(ncol(X.sigma2)>1){X.sigma2 = X.sigma2[,na.omit(match( rownames(obj$estimation),paste("sigma 2: ",colnames(X.sigma2), sep = "")))]}
   
   if(dist1 == "Multinom"){idx.par1 = match(paste("mu 1: ",rep(colnames(X.mu1),n.mult)," (p=",rep(1:n.mult,each=ncol(X.mu1)),")", sep = ""), rownames(obj$estimation))}
   if(dist1 != "Multinom"){idx.par1 = match(paste("mu 1: ",colnames(X.mu1),sep=""), rownames(obj$estimation))}
   
   idx.par3 = match(paste("mu 2: ",colnames(X.mu2),sep=""), rownames(obj$estimation))
   idx.par4 = match(paste("sigma 2: ",colnames(X.sigma2),sep=""), rownames(obj$estimation))
   
   if(dist1 != "Multinom"){mu1 = matrix(obj$transforms[[1]](X.mu1 %*% estim[idx.par1]),ncol=1)}
   if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(estim[idx.par1],nrow = n.mult)))) / rowSums(exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(estim[idx.par1],nrow = n.mult)))))}
   
   mu2 = obj$transforms[[3]](X.mu2 %*% estim[idx.par3])
   sigma2 = as.matrix(obj$transforms[[4]](X.sigma2 %*% estim[idx.par4]))
   
   if(dist1 == "Nbinom" | dist1 == "Pois"){exp1 = mu1}
   if(dist1 == "Binom"){exp1 = mu1*sigma1}
   if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
   if(dist1 == "Multinom"){exp1 = as.numeric(apply(mu1,1,which.max))-1}
   
   if(dist2 == "Gauss"){exp2 = mu2}
   if(dist2 == "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}
   
   R.Cmd <- sprintf("y.pred = r%s%s(nrow(mu1), mu1,sigma1, mu2,sigma2)", dist1, dist2)
   
   eval(parse(text=R.Cmd))
   
   #y.pred = exp1 + exp2
   
   #R.Cmd <- sprintf("y.pred.ci = q%s%s(c(0.025, 0.975), mu1,sigma1, mu2,sigma2)", dist1, dist2)
   
   #eval(parse(text=R.Cmd))
   
   pred = data.frame(y.fitted = y.pred)
   sample = NULL
}

   if(pred.type == "latent" & !is.null(newdata)){
      
      dist1= unlist(strsplit(obj$distname,"/"))[1]
      dist2= unlist(strsplit(obj$distname,"/"))[2]
      
      R.Cmd <- sprintf("p.y = d%s%s(obj$xdata$y, exp(obj$estimation$Estimate[1]), exp(obj$estimation$Estimate[2]),
                   obj$estimation$Estimate[3],exp(obj$estimation$Estimate[4]))", dist1, dist2)
      
      eval(parse(text=R.Cmd))
      
      #p.y = dNbinomGauss(obj$xdata$y, exp(obj$estimation$Estimate[1]), exp(obj$estimation$Estimate[2]),
      #             obj$estimation$Estimate[3],exp(obj$estimation$Estimate[4]))
      
      R.Cmd <- sprintf("p.y.k = sapply(0:50,function(x) dk%s%s(obj$xdata$y, exp(obj$estimation$Estimate[1]), exp(obj$estimation$Estimate[2]),
                            obj$estimation$Estimate[3],exp(obj$estimation$Estimate[4]),x))", dist1, dist2)
      #p.y.k = sapply(0:50,function(x) dkNbinomGauss(obj$xdata$y, exp(obj$estimation$Estimate[1]), exp(obj$estimation$Estimate[2]),
      #                      obj$estimation$Estimate[3],exp(obj$estimation$Estimate[4]),x))
      
      eval(parse(text=R.Cmd))
      
      post.dist = p.y.k / p.y
      
      k.est = apply(post.dist,1,which.max) -1 
      
      plot(k.est,obj$xdata$y )
      
      hist(k.est,breaks = (0:20)/2-0.25, xlim = c(-1,7))
      
      par(mfrow = c(1,1),mar = c(5.1,4.1,4.1,2.1))
      plot(0:20,(p.y.k[3,]/p.y[3]))
      
   }
   
   
   
   
   
   if(with.confidence.interval)
      return(list(pred=pred, dist=sample)) else 
         return(pred$y.fitted)
}
