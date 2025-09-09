#' Plotting  Bayesian B-spline model
#' 
#' @description function for plotting  Bayesian B-spline model for estimating reference series and/or assessing aging bias
#' 
#' @param df data list object from `data_prep` and used to estimate the STAN model
#' @param ext posterior draws extraction list from `extract_draws` function
#' @param probs vector of three probabilities for lower CI, centrality, upper CI
#' @param post.den logical flag to plot posterior densities of certain parameters
#' @param legend logical flag to plot legend
#' 
#' @return figure only, invisible
#' 
#' @export
#' 
#' @examples
#' #REFERENCE-ONLY
#' df <- data_prep(sim_ref)
#' 
#' #estimate model
#' fit <- est_model(df, show_messages = FALSE, show_exceptions = FALSE)
#' 
#' #extract draws
#' draws <- extract_draws(fit)
#' 
#' #plot fit
#' plot_fit(df, draws)
#' 
#' #INTEGRATED
#' \dontrun{
#' df_int <- data_prep(sim_ref, sim_unk)
#' fit_int <- est_model(df_int)
#' draws_int <- extract_draws(fit_int)
#' plot_fit(df_int, draws_int)
#' }
#' 
plot_fit <- function(df, ext, probs=c(0.05,0.5,0.95),post.den=TRUE, legend=TRUE, min.BY, max.BY){
	bias_flag <- df$flag == 'integrated'
	pred.q <- apply(ext$C14_pred, 2, quantile, probs=probs)
	if(bias_flag) obs.q <- apply(ext$BY_adj, 2, quantile, probs=probs)+df$data$BY_mu
	df$data$BY_ref_bck <- df$data$BY_ref+df$data$BY_mu
	if(bias_flag) df$data$BY_obs_bck <- df$data$BY_obs+df$data$BY_mu
	df$data$BY_pred_bck <- df$data$BY_pred + df$data$BY_mu
	if(bias_flag){
		BY.r <- range(c(df$data$BY_ref_bck,
		                df$data$BY_obs_bck,
		                df$data$BY_pred_bck))
	}else{
		BY.r <- range(c(df$data$BY_ref_bck,
		                df$data$BY_pred_bck))
	}
	if(!missing(min.BY)){
		BY.r[1] <- min.BY
	}
	if(!missing(max.BY)){
		BY.r[2] <- max.BY
	}

	

	
	C14.r <- range(c(df$data$C14_ref, df$data$C14_obs))
	
	BY.dr <- abs(diff(BY.r))
	C14.dr <- abs(diff(C14.r))

	par(mar=c(4,4.5,1,1), mgp=c(2.5,1,0))
	
	plot(df$data$BY_ref_bck, df$data$C14_ref,las=1, 
	     xlim=BY.r+c(-0.05,0)*BY.dr,
	     ylim=C14.r+c(-0.02,0.02)*C14.dr,
	     type='n', xlab="Reference/Birth Year",
	     ylab=expression(Delta^14*C~"(\u2030)"))
	points(df$data$BY_ref_bck, df$data$C14_ref,
	       pch=21,bg='gray90',col='gray30')
	lines(df$data$BY_pred_bck,pred.q[2,],lwd=3)
	lines(df$data$BY_pred_bck,pred.q[1,], 
	      lty=3,lwd=2)
	lines(df$data$BY_pred_bck,pred.q[3,], 
	      lty=3,lwd=2)

	peak <- df$data$BY_pred_bck[which.max(pred.q[2,])]

	peak_frac <- (par('usr')[2]-peak)/diff(par('usr')[1:2])

	#bump post.den
	if(peak_frac < 0.8){
		xlow <- peak_frac-0.2
		ylow <- 0.4
		ylow.adj <- 0.1
		legend.loc <- 'topright'
	}else{
		xlow <- 0.1
		ylow <- 0.75
		legend.loc <- 'bottomright'
	}
	xhigh <- xlow + 0.003 * diff(par('usr')[1:2])
	yhigh <- ylow + 0.125
	yhigh.adj <- ylow.adj + 0.125
	sig.y <- ylow + 0.18
	adj.y <- ylow.adj + 0.18
	
	if(bias_flag){
		points(df$data$BY_obs_bck, df$data$C14_obs, 
       	pch=21,bg='firebrick4',col='firebrick1')
		points(obs.q[2,], df$data$C14_obs, pch="|",
		       col='firebrick1')
		segments(y0 = df$data$C14_obs,
		         y1 = df$data$C14_obs,
		         x0 = obs.q[1,],
		         x1 = obs.q[3,],
		         col='firebrick2')
	}
	
	#sigma
		if(bias_flag){
			sigE.r <- range(c(ext$sigma_ref,
		                  ext$sigma_obs))
		}else{
			sigE.r <- range(c(ext$sigma_ref))
		}
		
		sigE.r <- c(-0.1,0.1)*diff(sigE.r)+sigE.r
		d.sigE <- den.sc(ext$sigma_ref,
		                  to.x=par.it(c(xlow,xhigh),1),
		                  to.y=par.it(c(ylow,yhigh),2),
		                  adj=2, from=sigE.r,
		                  probs=probs)
		polygon(x=c(d.sigE$d$x,rev(d.sigE$d$x)),
		        y=c(rep(min(d.sigE$d$y),length(d.sigE$d$y)),
		            rev(d.sigE$d$y)),
		        col=col2rgbA('orange',0.5),
		        border='black',lwd=2)
		if(post.den){
			segments(x0=d.sigE$d$q[,1],
			         y0 = rep(min(d.sigE$d$y),length(probs)),
			         x1=d.sigE$d$q[,1],
			         y1=d.sigE$d$q[,2],
			         col='black',
			         lty=3)
		}
		max_sig_ref <- d.sigE$d$x[which.max(d.sigE$d$y)]
		if(!bias_flag){
			text(max_sig_ref,
			     par.it(sig.y,2),
			     expression(sigma[REF]),
			     cex=2,font=2)
		}
		axis(1, at = d.sigE$ax$p,
		     labels = d.sigE$ax$v,
		     pos = min(d.sigE$d$y),
		     tck = -0.0075,mgp = c(3,0.05,0),
		     cex.axis = 0.6)
		if(bias_flag){
			d.sigE <- den.sc(ext$sigma_obs,
		                  to.x=par.it(c(xlow,xhigh),1),
		                  to.y=par.it(c(ylow,yhigh),2),
		                  adj=2, from = sigE.r,
		                  probs=probs)
			polygon(x=c(d.sigE$d$x,rev(d.sigE$d$x)),
			        y=c(rep(min(d.sigE$d$y),
			                length(d.sigE$d$y)),
			            rev(d.sigE$d$y)),
			        col=col2rgbA('orange',0.5),
			        border='darkorange4',lwd=2)
			if(post.den){
				segments(x0=d.sigE$d$q[,1],
				         y0 = rep(min(d.sigE$d$y),length(probs)),
				         x1=d.sigE$d$q[,1],
				         y1=d.sigE$d$q[,2],
				         col='darkorange4',
				         lty=3)
			}
			max_sig_obs <- d.sigE$d$x[which.max(d.sigE$d$y)]
			if(max_sig_ref < max_sig_obs){
				text(min(d.sigE$d$x),
				     par.it(sig.y,2),
				     expression(sigma[REF]),
				     col = 'black',
				     cex=2,font=2, adj = c(0,0.5))
				text(max(d.sigE$d$x),
				     par.it(sig.y,2),
				     expression(sigma[OBS]),
				     col = 'darkorange4',
				     cex=2,font=2, adj = c(1,0.5))
			}else{
				text(min(d.sigE$d$x),
				     par.it(sig.y,2),
				     expression(sigma[OBS]),
				     col = 'darkorange4',
				     cex=2,font=2, adj = c(0,0.5))
				text(max(d.sigE$d$x),
				     par.it(sig.y,2),
				     expression(sigma[REF]),
				     col = 'black',
				     cex=2,font=2, adj = c(1,0.5))
			}
			
		}
		
	#Age bias
		if(bias_flag){
			d.dmu <- den.sc(ext$adj,
			                  to.x=par.it(c(xlow,xhigh),1),
			                  to.y=par.it(c(ylow.adj,yhigh.adj),2),
			                  adj=2, ret.x=TRUE,
			                  probs=probs)
			polygon(x=c(d.dmu$d$x,rev(d.dmu$d$x)),
			        y=c(rep(min(d.dmu$d$y),length(d.dmu$d$y)),
			            rev(d.dmu$d$y)),
			        col=col2rgbA('gray',0.5),
			        border='black',lwd=2)
			if(post.den){
				segments(x0=d.dmu$d$q[,1],
				         y0 = rep(min(d.dmu$d$y),length(probs)),
				         x1=d.dmu$d$q[,1],
				         y1=d.dmu$d$q[,2],
				         col='black',
				         lty=3)
			}
			text(mean(d.dmu$d$x),
			     par.it(adj.y,2),
			     expression(BY[adj]),
			     cex=1.5,font=2)
			axis(1,at=d.dmu$ax$p,d.dmu$ax$v,
			     pos=min(d.dmu$d$y),tck=-0.0075,mgp=c(3,0.05,0),
			     cex.axis=0.6)
			x.pos <- scales::rescale(c(0),
			                from=range(d.dmu$x),
			                to=par.it(c(xlow,xhigh),1))
			if(x.pos > xlow & x.pos < xhigh){
				segments(x.pos,par.it(c(ylow.adj),2),
			         x.pos,max(d.dmu$d$y),
			         lty=2,col=c('red'))
			}
			
		}
	# Age bias multi
		if(bias_flag & legend){
			legend(legend.loc,
		       legend=c(expression(tilde(y)[REF]),
		                expression(paste('90% CI',phantom()[REF])),
		                expression(tilde(y)[OBS]),
		                expression(paste('90% CI',phantom()[OBS])),
		                'Ref.',
		                'Obs.'), 
		       bty='n',lty=c(1,3,NA,1,NA,NA), lwd=3, 
		       pch=c(NA,NA,3,NA,21,21), y.intersp=1,
		       col=c('black','gray50','firebrick4',
		             'firebrick1','gray50','firebrick1'),
		       pt.bg=c(NA,NA,'firebrick1',NA,'gray90','firebrick4'), 
		       pt.lwd=1, pt.cex=2)
		}else if(legend){
			legend(legend.loc,
		       legend=c(expression(tilde(y)[REF]),
		                expression(paste('90% CI',phantom()[REF])),
		                'Ref.'), 
		       bty='n',lty=c(1,3,NA), lwd=3, 
		       pch=c(NA,NA,21), y.intersp=1,
		       col=c('black','gray50','gray50'),
		       pt.bg=c(NA,NA,'gray90'), 
		       pt.lwd=1, pt.cex=2)
		}
	invisible()
}

den.rel <- function(x,top=NULL,probs,...){
  x <- na.omit(unlist(x))

  d <- density(x, ...)
  d$y <- d$y/max(d$y)
  if(!is.null(top)) d$y <- d$y*top
  if(!missing(probs)){
    q <- quantile(x,probs)
    q.id <- sapply(q, function(x)which.min(abs(d$x-x)))
    d$q <- cbind(d$x[q.id],d$y[q.id])
  }
  return(d)
}
den.sc <- function(x,to.x,to.y,from,ret.x=FALSE,probs,...){
  x <- na.omit(unlist(x))

  if(!missing(from)){
    d <- density(x, from=from[1],to=from[2],...)
  }else{
    d <- density(x, ...)
  }
  dx <- d$x
  d$y <- d$y/max(d$y)
  px <- pretty(d$x)
  if(!missing(to.y)){
    d$y <- scales::rescale(d$y,to=to.y)
  }
  if(!missing(to.x)){
    x.ax <- scales::rescale(px,from=range(dx),to=to.x)
    d$x <- scales::rescale(d$x,to=to.x)
    if(!missing(probs)){
      q <- quantile(x,probs)
      q.id <- sapply(q, function(x)which.min(abs(dx-x)))
      d$q <- cbind(d$x[q.id],d$y[q.id])
    }
    if(ret.x){
      return(list(d=d,ax=data.frame(v=px,
                                p=x.ax),x=dx))
    }else{
      return(list(d=d,ax=data.frame(v=px,
                                p=x.ax)))
    }
  }else{
    if(!missing(probs)){
      q <- quantile(x,probs)
      q.id <- sapply(q, function(x)which.min(abs(dx-x)))
      d$q <- cbind(d$x[q.id],d$y[q.id])
    }
    if(ret.x){
      return(list(d=d,x=dx))
    }else{
      return(list(d=d))
    }
    
  }
}
par.it <- function(p,side){
  if(side==2){
     x <- par('usr')[3]+abs(diff(par('usr')[3:4]))*p
  }else if(side==1){
     x <- par('usr')[1]+abs(diff(par('usr')[1:2]))*p
  }else if(side==3){
     x <- par('usr')[4]+abs(diff(par('usr')[3:4]))*p
  }else if(side==4){
     x <- par('usr')[2]+abs(diff(par('usr')[1:2]))*p
  }
  return(x)
}
col2rgbA<-function(color,transparency){
  rgb(t(col2rgb(color))/255,alpha=transparency)
}