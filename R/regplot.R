#' Convolutive regression function regression plot
#'
#' Convolutive regression function regression plot
#' @param obj convreg object
#' @param delta discretization parameter
#' @param disp choice of predictors to plot
#' @param cex.pt point size for observations
#' @param ... additional plotting parameters
#' @keywords convreg
#' @export
#' @examples
#' set.seed(123)
#' x1 = rnorm(500,0.75,0.25)
#' x2 = rnorm(500,-2,1)
#' k= rnbinom(500,mu=exp(x2),size=0.5)
#' e= x1 + rnorm(n=500,mean=0,sd=0.25)
#' y= data.frame(obs=k + e , f1 = x1, f2 = x2)
#' res=convreg( ~obs,
#'           formula.mu1 =~ f2,
#'           formula.mu2 =~ f1,
#'           data=y)
#' regplot(res, cex.pt = 0.25)
regplot <- function(obj,delta = 0.5,disp = "all", cex.pt = 0.5, ...)
{


   dist1= unlist(strsplit(obj$distname,"/"))[1]
   dist2= unlist(strsplit(obj$distname,"/"))[2]

   p = predict(obj)

   idx.boot = sample(1:length(p$y), length(p$y), replace=T)

   variables = unlist(lapply(obj$formulas, all.vars))
   variables = variables[variables %in% names(obj$data)]
   names(variables)[grep("mu1",names(variables))] <- paste0(dist1,"")
   names(variables)[grep("mu2",names(variables))] <- paste0(dist2,"")
   names(variables)[grep("sigma1",names(variables))] <- paste0(dist1,"")
   names(variables)[grep("sigma2",names(variables))] <- paste0(dist2,"")

   names = names(variables)

   variables = unique(variables)

   par(mfrow = c(1,1))
   nb.plot.final = length(variables)

   if(disp=="all"){
      par(mfrow=c(max(round(sqrt(nb.plot.final)),1),max(round(sqrt(nb.plot.final+2)),1)))
      if(nb.plot.final==1){par(mfrow = c(1,1))}
   }

   if(length(variables)==0){
      message("No variable to display")
   }else{

      for(ii in 1:length(variables)){

         var.x = obj$data[,which(names(obj$data) %in% variables[ii])]

         bp = is.factor(var.x)

         if(bp){
            init.levels = levels(var.x)
            tb = table(var.x)

            idx.ch = which(var.x %in% names(which(tb <=2)))
            if(length(idx.ch)>0){
               var.x[idx.ch] = sample(setdiff(levels(var.x), names(which(tb <=2))), length(idx.ch))
               var.x = as.factor(as.character(var.x))
            }
            corr.levels = levels(var.x)
         }

         reg = marginalReg(obj, name.reg = variables[ii], bp)

         pred = reg$value

         calc.m.s = meanVarianceConvreg(mu1 = pred$mu1, sigma1 = pred$s1, mu2 = pred$mu2, sigma2 = pred$s2,
                                          dist1 = obj$profile.lik$other[[16]], dist2 = obj$profile.lik$other[[17]])

         if(!bp){

            y.p = calc.m.s$mean

            y.q = eval(parse(text = sprintf("t(matrix(q%s%s(c(0.1,0.25,0.75,0.9), pred$mu1,pred$s1, pred$mu2,pred$s2), ncol=4))", dist1, dist2)))


            plot(pred$x, y.p,type="l", ylim = range(p$y),col = rgb(0,0,0,0.3), xlab = variables[ii], ylab = "")
                 #, ...)
            polygon(c(pred$x,rev(pred$x)),
                    c(y.q[1,],rev(y.q[4,])),col = rgb(0,0,0,0.1),border=NA)
            polygon(c(pred$x,rev(pred$x)),
                    c(y.q[2,],rev(y.q[3,])),col = rgb(0,0,0,0.2),border=NA)
            points(var.x, p$y, cex = cex.pt , pch = 3)
         }

         if(bp){

            keep.levels = which(init.levels %in% corr.levels)
            y.q = eval(parse(text = sprintf("t(matrix(q%s%s(c(0.05,0.25,0.75,0.95), pred$mu1,pred$s1, pred$mu2,pred$s2), ncol=4))", dist1, dist2)))

            mean.y = calc.m.s[[1]][keep.levels]
            min.y = y.q[2,keep.levels]
            max.y = y.q[3,keep.levels]
            min2.y = y.q[1,keep.levels]
            max2.y = y.q[4,keep.levels]

            plot(1,type=  "n",xlab = variables[ii], ylab = "",xaxt = "n",
                 xlim = c(0.5,nlevels(var.x)+0.5),ylim = range(min2.y,max2.y), main = names(variables)[ii])


            axis(1,at = 1:nlevels(var.x), labels = levels(var.x))



            for(jj in 1:nlevels(var.x)){
               points(c(jj-0.2,jj+0.2),rep(mean.y[jj],2), col=rgb(0.5,0.5,.5,0.8), type= "l",lwd=5)

               polygon(c(c(jj-0.2,jj+0.2),rev(c(jj-0.2,jj+0.2))),
                       c(rep(min2.y[jj],2),rev(rep(max2.y[jj],2))),
                       col = rgb(0,0,0,0.1),border="black")
               polygon(c(c(jj-0.2,jj+0.2),rev(c(jj-0.2,jj+0.2))),
                       c(rep(min.y[jj],2),rev(rep(max.y[jj],2))),
                       col = rgb(0,0,0,0.3),border=NA)

            }

         }

         if(disp=="seq"){par(ask = T)}


      }
   }
}
