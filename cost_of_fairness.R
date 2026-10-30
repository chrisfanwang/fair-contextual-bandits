# Population cost of Lipschitz fairness and its reduction from L = 0.
# Run with source(file.choose()) or Rscript cost_of_fairness.R.
# --plot-only redraws saved data.
.cost_entry_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  path <- if(length(paths))tail(paths,1L)[[1L]] else args[[1L]]
  if(!file.exists(path))path<-gsub("~+~"," ",path,fixed=TRUE)
  normalizePath(path,mustWork=TRUE)
})
.cost_release_dir <- dirname(.cost_entry_file)
source(file.path(.cost_release_dir, "algorithms.R"))

# These are the coefficients used by all ordinary simulation experiments.
get_beta_setting <- function(regime, d) {
  regime <- match.arg(
    regime,
    c(
      "homogeneous",
      "heterogeneous",
      "heterogeneous_d2_arm1_small",
      "heterogeneous_d2_arm2_small"
    )
  )
  d <- as.integer(d)
  
  if (regime == "homogeneous") {
    beta1_slope <- rep(1 / (2 * sqrt(d)), d)
    beta2_slope <- rep(-1 / (2 * sqrt(d)), d)
    theta_slope <- beta1_slope - beta2_slope
    theta0 <- -0.5 * sum(theta_slope)
    
  } else {
    if (regime %in% c("heterogeneous_d2_arm1_small", "heterogeneous_d2_arm2_small") && d != 2L) {
      stop(
        paste0(
          regime,
          " is a two-dimensional stress-test setting and is only defined for d = 2."
        ),
        call. = FALSE
      )
    }
    
    if (d == 1L) {
      beta1_slope <- c(0.8)
      beta2_slope <- c(-0.4)
      theta_slope <- beta1_slope - beta2_slope
      theta0 <- -0.5 * sum(theta_slope)
      
    } else if (d == 2L) {
      
      if (regime == "heterogeneous") {
        # Baseline heterogeneous d = 2 setting:
        beta1_slope <- c(0.7, 0.4)
        beta2_slope <- c(-0.3, -0.6)
        theta_slope <- beta1_slope - beta2_slope
        theta0 <- -0.5 * sum(theta_slope)
        
      } else if (regime == "heterogeneous_d2_arm1_small") {
        # Arm-1-small stress-test boundary:
        theta0 <- -0.82
        theta_slope <- c(0.27, 1.00)
        beta1_slope <- theta_slope / 2
        beta2_slope <- -theta_slope / 2
        
      } else if (regime == "heterogeneous_d2_arm2_small") {
        # Arm-2-small stress-test boundary:
        beta1_slope <- c(0.7, 0.4)
        beta2_slope <- c(-0.3, -0.6)
        theta_slope <- beta1_slope - beta2_slope
        theta0 <- -0.2
      }
      
    } else if (d == 3L) {
      beta1_slope <- c(0.6, -0.7, 0.15)
      beta2_slope <- c(-0.4, 0.3, 0.15)
      theta_slope <- beta1_slope - beta2_slope
      theta0 <- -0.5 * sum(theta_slope)
      
    } else if (d == 4L) {
      beta1_slope <- c(0.6, 0.35, -0.3, -0.6)
      beta2_slope <- c(-0.4, -0.65, 0.7, 0.4)
      theta_slope <- beta1_slope - beta2_slope
      theta0 <- -0.5 * sum(theta_slope)
      
    } else {
      stop("Heterogeneous setting is only defined for d in {1,2,3,4}.", call. = FALSE)
    }
  }
  
  beta1 <- c(theta0 / 2, beta1_slope)
  beta2 <- c(-theta0 / 2, beta2_slope)
  
  list(beta1 = beta1, beta2 = beta2)
}

# Fixed, unscrambled Halton sequence: skip the first 128 points.
halton_sequence <- function(n, d, skip = 128L) {
  stopifnot(n >= 1L, d %in% 1:4, skip >= 0L)
  bases <- c(2L, 3L, 5L, 7L)
  out <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    index <- seq_len(n) + skip
    fraction <- 1 / bases[j]
    while (any(index > 0)) {
      out[, j] <- out[, j] + fraction * (index %% bases[j])
      index <- index %/% bases[j]
      fraction <- fraction / bases[j]
    }
  }
  out
}

unit_to_context <- function(U, covariate_dist, sigma_X = .18) {
  U <- as.matrix(U)
  if (covariate_dist == "uniform") return(U)
  stopifnot(covariate_dist == "truncated_gaussian", sigma_X > 0)
  lo <- pnorm(0, .5, sigma_X); hi <- pnorm(1, .5, sigma_X)
  matrix(qnorm(lo + U * (hi - lo), .5, sigma_X), nrow(U), ncol(U))
}

generate_contexts <- function(T, d, covariate_dist, sigma_X = .18) {
  unit_to_context(matrix(runif(T * d), nrow = T, ncol = d), covariate_dist, sigma_X)
}

# Shared simulation and population-evaluation settings.
mab_simulation_settings <- function() list(
  seed = 20260924L, horizons = as.integer(2^(10:15)), repetitions = 50L,
  L = .5, nu = 1, lambda = 1, sigma = .3, eta = .05, C_beta = 2,
  ordinary = list(sigma_X = .18),
  evaluation = list(qmc_n_d3 = 16384L, qmc_n_d4 = 32768L, halton_skip = 128L))

# Analytic integration in d = 1 uses moments in x - 1/2.
mab_d1_density <- function(covariate, sigma_X = .18) {
  stopifnot(covariate %in% c("uniform", "truncated_gaussian"))
  list(kind = covariate, sigma_X = sigma_X)
}
mab_d1_power_integral <- function(lo,hi,k) {
  # Midpoint expansion avoids subtraction of nearly equal endpoint powers.
  center<-(lo+hi)/2;half<-(hi-lo)/2
  sum(vapply(seq.int(0L,k,by=2L),function(j)
    choose(k,j)*center^(k-j)*2*half^(j+1)/(j+1),numeric(1)))
}

mab_d1_moments <- function(lo,hi,density) {
  lo<-max(-.5,lo);hi<-min(.5,hi)
  if(hi<=lo)return(c(0,0,0))
  if(density$kind=="uniform")return(vapply(0:2,function(k)mab_d1_power_integral(lo,hi,k),numeric(1)))
  if(density$kind=="truncated_gaussian") {
    sd<-density$sigma_X;a<-lo/sd;b<-hi/sd
    mass<-if(a>=0)pnorm(a,lower.tail=FALSE)-pnorm(b,lower.tail=FALSE) else pnorm(b)-pnorm(a)
    normalizer<-pnorm(.5/sd)-pnorm(-.5/sd)
    return(c(mass,sd*(dnorm(a)-dnorm(b)),sd^2*(mass+a*dnorm(a)-b*dnorm(b)))/normalizer)
  }
  stop("Unsupported density")
}

mab_d1_affine_segments <- function(a,b) {
  knots<-c(-.5,.5)
  if(b!=0)knots<-sort(unique(c(knots,pmax(-.5,pmin(.5,c(-a/b,(1-a)/b))))))
  mid<-(head(knots,-1)+tail(knots,-1))/2;value<-a+b*mid
  data.frame(lo=head(knots,-1),hi=tail(knots,-1),
    a=ifelse(value<=0,0,ifelse(value>=1,1,a)),b=ifelse(value<=0|value>=1,0,b))
}

mab_d1_policy_segments <- function(policy) {
  stopifnot(policy$d==1L)
  p<-policy$params
  if(policy$type=="constant")return(mab_d1_affine_segments(p$prob,0))
  if(policy$type=="population_analytic_oracle") {
    stopifnot(p$kind=="centered_signed_sum")
    return(mab_d1_affine_segments(.5,p$L*p$signs[1]))
  }
  if(policy$type=="analytic_d1_clipped_affine")return(mab_d1_affine_segments(p$a,p$b))
  if(policy$type=="single_sign_exact_oracle")return(mab_d1_affine_segments(p$q,p$L))
  if(policy$type=="boundary_canonical") {
    stopifnot(nrow(p$H_ref)==1L,length(p$q_ref)==1L,p$nu==1)
    b<-sign(p$theta[2])*p$L
    return(mab_d1_affine_segments(p$q_ref+b*(.5-p$H_ref[1,1]),b))
  }
  stop("Unsupported policy for analytic d1 evaluation: ",policy$type)
}

mab_d1_gap_value <- function(segments,theta,density) {
  A<-theta[1]+.5*theta[2];B<-theta[2]
  sum(vapply(seq_len(nrow(segments)),function(i) {
    z<-segments[i,];sum(c(A*z$a,A*z$b+B*z$a,B*z$b)*mab_d1_moments(z$lo,z$hi,density))
  },numeric(1)))
}

mab_d1_policy_value <- function(policy,beta,density) {
  moments<-mab_d1_moments(-.5,.5,density)
  baseline<-(beta$beta2[1]+.5*beta$beta2[2])*moments[1]+beta$beta2[2]*moments[2]
  baseline+mab_d1_gap_value(mab_d1_policy_segments(policy),beta$beta1-beta$beta2,density)
}

mab_d1_regret_rate <- function(oracle,policy,beta,density) {
  a<-mab_d1_policy_segments(oracle);b<-mab_d1_policy_segments(policy)
  knots<-sort(unique(c(a$lo,a$hi,b$lo,b$hi)));mid<-(head(knots,-1)+tail(knots,-1))/2
  ia<-findInterval(mid,c(a$lo,tail(a$hi,1)));ib<-findInterval(mid,c(b$lo,tail(b$hi,1)))
  difference<-data.frame(lo=head(knots,-1),hi=tail(knots,-1),a=a$a[ia]-b$a[ib],b=a$b[ia]-b$b[ib])
  # Integrate the difference directly to avoid cancellation of common rewards.
  mab_d1_gap_value(difference,beta$beta1-beta$beta2,density)
}


# Integrate affine policy regions analytically in d = 2.
mab_d2_square <- function() rbind(c(-.5,-.5),c(.5,-.5),c(.5,.5),c(-.5,.5))

mab_d2_affine <- function(P,b) b[1]+P[,1]*b[2]+P[,2]*b[3]

mab_d2_center <- function(theta) c(theta[1]+.5*sum(theta[-1]),theta[-1])

mab_d2_area <- function(P) {
  if(nrow(P)<3L)return(0)
  Q<-sweep(P,2,P[1,]);j<-c(2:nrow(Q),1L)
  sum(Q[,1]*Q[j,2]-Q[,2]*Q[j,1])/2
}

mab_d2_clip <- function(P,h) {
  if(nrow(P)<3L)return(matrix(numeric(),ncol=2))
  scale<-max(abs(h));if(scale==0)return(P)
  f<-mab_d2_affine(P,h/scale);inside<-f>=0
  if(all(inside))return(P)
  if(!any(inside))return(matrix(numeric(),ncol=2))
  out<-list();k<-0L;n<-nrow(P)
  for(i in seq_len(n)) {
    j<-if(i==n)1L else i+1L
    if(inside[i]){k<-k+1L;out[[k]]<-P[i,]}
    if(inside[i]!=inside[j]){k<-k+1L;out[[k]]<-P[i,]+f[i]/(f[i]-f[j])*(P[j,]-P[i,])}
  }
  Q<-unique(do.call(rbind,out))
  if(nrow(Q)<3L||mab_d2_area(Q)<=0)return(matrix(numeric(),ncol=2))
  Q
}

mab_d2_clip_many <- function(P,H) {
  for(i in seq_len(nrow(H))) {P<-mab_d2_clip(P,H[i,]);if(nrow(P)<3L)break}
  P
}

mab_d2_cell <- function(P,b) list(P=P,b=as.numeric(b))

mab_d2_envelope_update <- function(cells,b,domain,minimum=TRUE) {
  out<-list()
  for(cell in cells) {
    improvement<-if(minimum)cell$b-b else b-cell$b
    # Suppress roundoff-sized duplicate planes, not positive-area policy features.
    if(max(mab_d2_affine(cell$P,improvement))<=32*.Machine$double.eps) {
      out[[length(out)+1L]]<-cell;next
    }
    H<-rbind(domain,improvement)
    if(!nrow(mab_d2_clip_many(cell$P,H))) {out[[length(out)+1L]]<-cell;next}
    remaining<-cell$P
    for(j in seq_len(nrow(H))) {
      f<-mab_d2_affine(remaining,H[j,])
      if(all(f>=0))next
      outside<-mab_d2_clip(remaining,-H[j,])
      if(nrow(outside))out[[length(out)+1L]]<-mab_d2_cell(outside,cell$b)
      remaining<-mab_d2_clip(remaining,H[j,])
      if(!nrow(remaining))break
    }
    if(nrow(remaining))out[[length(out)+1L]]<-mab_d2_cell(remaining,b)
  }
  if(length(out)>getOption("mab.max_d2_regions",100000L))
    stop("Exact d2 policy partition exceeded its region budget; no grid fallback was used.")
  out
}

mab_d2_clipped_regions <- function(b) {
  P<-mab_d2_square();out<-list()
  for(z in list(list(H=matrix(-b,1),v=c(0,0,0)),
    list(H=rbind(b,c(1,0,0)-b),v=b),list(H=matrix(b-c(1,0,0),1),v=c(1,0,0)))) {
    # A constant at a clipping threshold belongs to one region, not two.
    if(all(b[-1]==0))return(list(mab_d2_cell(P,c(clip01(b[1]),0,0))))
    q<-mab_d2_clip_many(P,z$H)
    if(nrow(q))out[[length(out)+1L]]<-mab_d2_cell(q,z$v)
  }
  out
}

mab_d2_trace_regions <- function(p) {
  # Simplify collinear knots for integration and check the reconstructed values.
  old_knots<-p$trace_knots;old_values<-p$q_ref
  retained<-logical(length(old_knots));retained[c(1L,length(old_knots))]<-TRUE
  intervals<-list(c(1L,length(old_knots)))
  while(length(intervals)) {
    ends<-intervals[[length(intervals)]];intervals[[length(intervals)]]<-NULL
    lo<-ends[1];hi<-ends[2];if(hi<=lo+1L)next
    ii<-seq.int(lo+1L,hi-1L)
    fitted<-old_values[lo]+(old_knots[ii]-old_knots[lo])/
      (old_knots[hi]-old_knots[lo])*(old_values[hi]-old_values[lo])
    discrepancy<-abs(old_values[ii]-fitted)
    if(max(discrepancy)<=8*.Machine$double.eps)next
    j<-ii[which.max(discrepancy)];retained[j]<-TRUE
    intervals<-c(intervals,list(c(lo,j),c(j,hi)))
  }
  keep<-which(retained);p$trace_knots<-old_knots[keep];p$q_ref<-old_values[keep]
  stopifnot(max(abs(approx(p$trace_knots,p$q_ref,xout=old_knots)$y-old_values))<=16*.Machine$double.eps)
  p$H_ref<-p$H_ref[keep,,drop=FALSE]
  theta<-mab_d2_center(p$theta);H<-p$H_ref-.5;q<-p$q_ref;L<-p$L
  a<-p$trace_origin-.5;v<-p$trace_direction;t<-p$trace_knots
  result<-list()
  for(side in c(-1,1)) {
    P<-mab_d2_clip(mab_d2_square(),side*theta);if(!nrow(P))next
    cells<-list(mab_d2_cell(P,c(as.numeric(side>0),0,0)))
    for(i in seq_len(nrow(H)))for(s1 in c(-1,1))for(s2 in c(-1,1)) {
      mab_check_deadline();signs<-c(s1,s2)
      domain<-rbind(c(-s1*H[i,1],s1,0),c(-s2*H[i,2],0,s2))
      b<-c(q[i]-side*L*sum(signs*H[i,]),side*L*signs)
      cells<-mab_d2_envelope_update(cells,b,domain,minimum=side>0)
    }
    # Query-dependent knots h_j(t)=x_j are needed in addition to fixed knots.
    for(i in seq_len(length(t)-1L))for(j in which(v!=0)) {
      k<-3L-j;dt<-t[i+1L]-t[i];dq<-(q[i+1L]-q[i])/dt
      b<-c(q[i]+dq*(-a[j]/v[j]-t[i]),0,0);b[j+1L]<-dq/v[j]
      # h_k = a_k + v_k*(x_j-a_j)/v_j.
      diff<-c(-a[k]+v[k]*a[j]/v[j],0,0)
      diff[k+1L]<-1;diff[j+1L]<- -v[k]/v[j]
      lo<-min(H[c(i,i+1L),j]);hi<-max(H[c(i,i+1L),j])
      lower<-upper<-numeric(3);lower[1]<- -lo;upper[1]<-hi
      lower[j+1L]<-1;upper[j+1L]<- -1
      for(sign in c(-1,1)) {
        candidate<-b+side*L*sign*diff
        cells<-mab_d2_envelope_update(cells,candidate,rbind(lower,upper,sign*diff),minimum=side>0)
      }
    }
    result<-c(result,cells)
  }
  result
}

mab_d2_policy_regions <- function(policy) {
  stopifnot(policy$d==2L);p<-policy$params
  if(policy$type=="constant")return(list(mab_d2_cell(mab_d2_square(),c(p$prob,0,0))))
  if(policy$type=="population_analytic_oracle")return(switch(p$kind,
    centered_signed_sum=mab_d2_clipped_regions(c(.5,p$L*p$signs)),
    uniform_arm1=mab_cost_uniform_arm1_regions(p$L),
    arm2_sum=mab_d2_clipped_regions(c(p$q+.8*p$L,p$L,p$L)),
    gaussian_arm1=if(p$L==0)list(mab_d2_cell(mab_d2_square(),c(0,0,0))) else NULL,
    stop("Unsupported analytic d2 oracle: ",p$kind)))
  if(policy$type=="analytic_d2_clipped_affine")return(mab_d2_clipped_regions(p$b))
  if(policy$type=="boundary_canonical_continuous_d2_l1")return(mab_d2_trace_regions(p))
  if(policy$type=="standard_ucb" && (p$alpha==0||identical(p$Vinv[[1]],p$Vinv[[2]]))) {
    gap<-mab_d2_center(p$beta[[1]]-p$beta[[2]])
    if(all(gap[-1]==0))return(list(mab_d2_cell(mab_d2_square(),c(as.numeric(gap[1]>=0),0,0))))
    out<-list()
    for(side in c(-1,1)) {
      P<-mab_d2_clip(mab_d2_square(),side*gap)
      if(nrow(P))out[[length(out)+1L]]<-mab_d2_cell(P,c(as.numeric(side>0),0,0))
    }
    return(out)
  }
  NULL # Nonlinear standard-UCB boundaries retain explicitly labelled numerical evaluation.
}

# Integrate a product of up to three affine factors using simplex moments.

mab_d2_polynomial_integral <- function(P,factors) {
  if(nrow(P)<3L)return(0)
  sum(vapply(2:(nrow(P)-1L),function(i) {
    triangle<-P[c(1L,i,i+1L),,drop=FALSE];area<-mab_d2_area(triangle)
    value<-vapply(seq_len(nrow(factors)),function(j)mab_d2_affine(triangle,factors[j,]),numeric(3))
    sums<-colSums(value)
    average<-switch(nrow(factors),sums[1]/3,
      (prod(sums)+sum(value[,1]*value[,2]))/12,
      (prod(sums)+sum(value[,1]*value[,2])*sums[3]+sum(value[,1]*value[,3])*sums[2]+
        sum(value[,2]*value[,3])*sums[1]+2*sum(value[,1]*value[,2]*value[,3]))/60)
    area*average
  },numeric(1)))
}

mab_d2_bvn <- function(x,y,rho) {
  if(!requireNamespace("pbivnorm",quietly=TRUE))stop("Install pbivnorm for Gaussian d2 polygon moments.")
  if(rho>=1)return(pnorm(pmin(x,y)))
  if(rho<= -1)return(pmax(0,pnorm(x)-pnorm(-y)))
  pbivnorm::pbivnorm(x,y,rho)
}

mab_d2_gaussian_moments <- function(P,sd) {
  Q<-P/sd;M<-numeric(6) # mass, s, t, s^2, s*t, t^2
  for(i in seq_len(nrow(Q))) {
    j<-if(i==nrow(Q))1L else i+1L
    u<-Q[i,1];v<-Q[i,2];du<-Q[j,1]-u;dv<-Q[j,2]-v
    A<-du^2+dv^2;if(A==0)next
    if(du!=0) {
      k<-(v*abs(du)-sign(du)*dv*u)/sqrt(A);rho<- -sign(du)*dv/sqrt(A)
      M[1]<-M[1]-diff(mab_d2_bvn(c(u,u+du),k,rho))
    }
    B<-u*du+v*dv;C<-u^2+v^2;mu<- -B/A;root<-sqrt(A)
    lower<- -mu*root;upper<-(1-mu)*root
    mass<-if(lower>=0)pnorm(lower,lower.tail=FALSE)-pnorm(upper,lower.tail=FALSE) else pnorm(upper)-pnorm(lower)
    J0<-exp(-.5*max(0,C-B^2/A))/sqrt(2*pi*A)*mass
    f0<-dnorm(u)*dnorm(v);f1<-dnorm(u+du)*dnorm(v+dv)
    J1<-mu*J0-(f1-f0)/A
    U<-u*J0+du*J1;V<-v*J0+dv*J1
    M[-1]<-M[-1]+c(-dv*J0,du*J0,-dv*U,(-dv*V+du*U)/2,du*V)
  }
  M[c(4,6)]<-M[c(4,6)]+M[1]
  M*c(1,sd,sd,sd^2,sd^2,sd^2)/(pnorm(.5/sd)-pnorm(-.5/sd))^2
}

mab_d2_region_integral <- function(P,factors,density) {
  if(!nrow(P))return(0)
  if(density$kind=="uniform")return(mab_d2_polynomial_integral(P,factors))
  if(density$kind=="truncated_gaussian") {
    stopifnot(nrow(factors)<=2L)
    a<-factors[1,];b<-if(nrow(factors)==2L)factors[2,] else c(1,0,0)
    coefficients<-c(a[1]*b[1],a[1]*b[2]+a[2]*b[1],a[1]*b[3]+a[3]*b[1],
      a[2]*b[2],a[2]*b[3]+a[3]*b[2],a[3]*b[3])
    return(sum(coefficients*mab_d2_gaussian_moments(P,density$sigma_X)))
  }
  stop("Unsupported density")
}

mab_d2_gap_value <- function(regions,theta,density) {
  gap<-mab_d2_center(theta)
  sum(vapply(regions,function(cell) {
    mab_check_deadline();mab_d2_region_integral(cell$P,rbind(gap,cell$b),density)
  },numeric(1)))
}

mab_d2_policy_value <- function(policy,beta,density,regions=mab_d2_policy_regions(policy)) {
  if(is.null(regions))stop("No analytic d2 integration formula for this nonlinear policy.")
  baseline<-mab_d2_region_integral(mab_d2_square(),matrix(mab_d2_center(beta$beta2),1),density)
  baseline+mab_d2_gap_value(regions,beta$beta1-beta$beta2,density)
}

mab_d2_density <- function(covariate, sigma_X = .18) {
  stopifnot(covariate %in% c("uniform", "truncated_gaussian"))
  list(kind = covariate, sigma_X = sigma_X)
}

mab_cost_uniform_sum <- function(d,L) {
  stopifnot(d%in%1:4,length(L)==1L,is.finite(L),L>=0)
  center <- d/2
  upper <- if(L==0)center else min(center,1/(2*L))
  total <- 0
  for(k in 0:d) {
    lower <- max(0,k-center)
    if(lower>=upper)next
    for(j in 0:(d-1)) {
      coefficient <- (-1)^k*choose(d,k)*choose(d-1,j)*(center-k)^(d-1-j)
      primitive <- function(t).5*t^(j+2)/(j+2)-L*t^(j+3)/(j+3)
      total <- total+coefficient*(primitive(upper)-primitive(lower))
    }
  }
  2*total/factorial(d-1)
}

mab_cost_region_loss <- function(cells,gap,density) {
  sum(vapply(cells,function(cell) {
    positive <- mab_d2_clip(cell$P,gap)
    negative <- mab_d2_clip(cell$P,-gap)
    mab_d2_region_integral(positive,rbind(gap,c(1,0,0)-cell$b),density)+
      mab_d2_region_integral(negative,rbind(-gap,cell$b),density)
  },numeric(1)))
}

mab_cost_uniform_arm1_planes <- function(L) {
  list(L*c(.13,.54,1),c(.5-.185*L,.27*L,L))
}

mab_cost_uniform_arm1_regions <- function(L) {
  if(L==0)return(list(mab_d2_cell(mab_d2_square(),c(0,0,0))))
  # clip(min(L*(y+0.54*x-0.64), 0.5+L*(y+0.27*x-0.82))).
  planes <- mab_cost_uniform_arm1_planes(L)
  cells <- list()
  for(i in 1:2)for(cell in mab_d2_clipped_regions(planes[[i]])) {
    polygon <- mab_d2_clip(cell$P,planes[[3-i]]-planes[[i]])
    if(nrow(polygon))cells[[length(cells)+1L]] <- mab_d2_cell(polygon,cell$b)
  }
  cells
}

mab_cost_symmetric_shape <- function(beta) {
  theta <- beta$beta1-beta$beta2
  slopes <- theta[-1]; active <- which(abs(slopes)>1e-13)
  scale <- max(abs(slopes)); signs <- sign(slopes)
  signs[-active] <- 0
  stopifnot(length(active)>0L,abs(theta[1]+.5*sum(slopes))<1e-13,
    max(abs(abs(slopes[active])-scale))<1e-13)
  list(scale=scale,signs=signs,active=active)
}

mab_cost_symmetric <- function(task,beta,cfg) {
  shape <- mab_cost_symmetric_shape(beta)
  theta <- beta$beta1-beta$beta2
  scale <- shape$scale; signs <- shape$signs; active <- shape$active
  answer <- list(eval_n=0L,analytic=TRUE,error_estimate=NA_real_,root_residual=NA_real_,
    boundary_q=NA_real_,halton_skip=NA_integer_,oracle_method="analytic_centered_signed_sum")
  if(task$covariate=="uniform") {
    answer$cost <- scale*mab_cost_uniform_sum(length(active),task$L)
    answer$cost0 <- scale*mab_cost_uniform_sum(length(active),0)
    answer$reward_method <- "analytic_uniform_sum_moments"
  } else if(task$d==1L) {
    density <- mab_d1_density(task$covariate,cfg$ordinary$sigma_X)
    upper <- if(task$L==0).5 else min(.5,1/(2*task$L))
    moments <- mab_d1_moments(0,upper,density)
    answer$cost <- scale*(moments[2]-2*task$L*moments[3])
    answer$cost0 <- scale*mab_d1_moments(0,.5,density)[2]
    answer$reward_method <- "analytic_normal_moments"
  } else if(task$d==2L) {
    density <- mab_d2_density(task$covariate,cfg$ordinary$sigma_X)
    gap <- mab_d2_center(theta)
    answer$cost <- mab_cost_region_loss(mab_d2_clipped_regions(c(.5,task$L*signs)),gap,density)
    answer$cost0 <- mab_cost_region_loss(mab_d2_clipped_regions(c(.5,0,0)),gap,density)
    answer$reward_method <- "analytic_gaussian_polygon_moments"
  } else {
    n <- cfg$evaluation[[paste0("qmc_n_d",task$d)]]
    skip <- cfg$evaluation$halton_skip
    stopifnot(length(n)==1L,is.finite(n),n>=1,n==as.integer(n),
      length(skip)==1L,is.finite(skip),skip>=0,skip==as.integer(skip))
    X <- unit_to_context(halton_sequence(n,task$d,skip),task$covariate,cfg$ordinary$sigma_X)
    score <- drop((X-.5)%*%signs)
    # Evaluate pointwise loss directly to avoid subtracting close reward means.
    answer$cost <- scale*mean(abs(score)*pmax(0,.5-task$L*abs(score)))
    answer$cost0 <- scale*mean(abs(score))/2
    answer$eval_n <- n
    answer$halton_skip <- skip
    answer$analytic <- FALSE
    answer$reward_method <- "halton_inverse_truncated_normal"
  }
  answer
}

# Serializable pointwise versions of the same analytic policies used above.
# No empirical LP or optimization grid is involved in their construction.

# The two asymmetric oracle formulas are valid for uniform covariates in d = 2.
# Arm 1 small: Delta(x) = x2 - (0.82 - 0.27*x1).
# Arm 2 small: Delta(x) = x1 + x2 - 0.2.
mab_population_oracle_policy <- function(d, parameters, meta = list()) {
  L <- parameters$L; kind <- parameters$kind
  stopifnot(d %in% 1:4, is.finite(L), L >= 0)
  evaluate <- function(X) {
    X <- as.matrix(X); stopifnot(ncol(X) == d)
    if (kind == "centered_signed_sum")
      return(clip01(.5 + L * drop((X - .5) %*% parameters$signs)))
    stopifnot(d == 2L)
    if (kind == "uniform_arm1") {
      h <- .82 - .27 * X[, 1]
      return(clip01(pmin(L * (X[, 2] + 1 - 2*h), .5 + L*(X[, 2] - h))))
    }
    if (kind == "arm2_sum") return(clip01(parameters$q + L*(rowSums(X) - .2)))
    stop("Unknown oracle kind")
  }
  structure(list(type = "population_analytic_oracle", d = d, params = parameters,
    meta = meta, eval_batch = evaluate,
    eval_one = function(x) evaluate(matrix(x, nrow = 1L))[[1L]]),
    class = "estimated_policy")
}

mab_ordinary_oracle <- function(d, regime, L = .5, covariate = "uniform") {
  beta <- get_beta_setting(regime, d)
  parameters <- list(L = L)
  if (regime == "heterogeneous_d2_arm1_small") {
    stopifnot(d == 2L, covariate == "uniform")
    parameters$kind <- "uniform_arm1"
  } else if (regime == "heterogeneous_d2_arm2_small") {
    stopifnot(d == 2L, covariate == "uniform")
    parameters$kind <- "arm2_sum"
    parameters$q <- if (L <= 10/3) 1 - .1*L else
      .5 + L*(.2 - (.2 + sqrt(.04 - 1/(3*L^2)))/2)
  } else {
    parameters$kind <- "centered_signed_sum"
    parameters$signs <- mab_cost_symmetric_shape(beta)$signs
  }
  mab_population_oracle_policy(d, parameters,
    list(fairness_guaranteed = TRUE, oracle_is_numerical_approximation = FALSE))
}

mab_ordinary_cost <- function(task, cfg = mab_simulation_settings()) {
  beta <- get_beta_setting(task$regime, task$d)
  if (grepl("_d2_", task$regime)) {
    stopifnot(task$d == 2L, task$covariate == "uniform")
    density <- mab_d2_density("uniform")
    gap <- mab_d2_center(beta$beta1 - beta$beta2)
    oracle <- mab_ordinary_oracle(2L, task$regime, task$L)
    p0 <- as.numeric(task$regime == "heterogeneous_d2_arm2_small")
    value <- list(cost = mab_cost_region_loss(mab_d2_policy_regions(oracle), gap, density),
      cost0 = mab_cost_region_loss(list(mab_d2_cell(mab_d2_square(), c(p0, 0, 0))), gap, density),
      reward_method = "analytic_polygon_moments", eval_n = 0L)
  } else value <- mab_cost_symmetric(task, beta, cfg)
  tolerance <- MAB_NUMERICAL_TOLERANCE
  stopifnot(is.finite(value$cost), value$cost >= -tolerance,
    value$cost <= value$cost0 + tolerance)
  constant <- max(vapply(beta, function(b) b[1] + .5*sum(b[-1]), numeric(1)))
  data.frame(d = task$d, covariate = task$covariate, regime = task$regime, L = task$L,
    cost_L = value$cost, cost0 = value$cost0,
    cost0_minus_cost_L = value$cost0 - value$cost,
    oracle_value = constant + value$cost0 - value$cost,
    unconstrained_value = constant + value$cost0,
    best_constant_value = constant, evaluation_method = value$reward_method,
    eval_n = value$eval_n)
}

mab_release_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") +
    ggplot2::theme(text = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "#DCDCDC", linewidth = .35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .5),
      panel.spacing = grid::unit(12.7, "pt"),
      strip.background = ggplot2::element_rect(fill = "#EBEBEB", colour = "black", linewidth = .5),
      strip.text = ggplot2::element_text(size = 12, colour = "black", face = "plain",
        margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)),
      axis.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.title = ggplot2::element_text(size = 13, colour = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = .45),
      axis.ticks.length = grid::unit(3, "pt"), legend.position = "bottom",
      legend.direction = "horizontal", legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(14, "pt"),
      legend.spacing.x = grid::unit(5, "pt"), legend.box.spacing = grid::unit(22, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(t = 57, r = 56, b = 10.5, l = 38))
}


mab_cost_plot <- function(values, metric) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  values$dimension <- factor(values$d, levels = 1:4)
  values$covariate <- factor(values$covariate, levels = c("uniform", "truncated_gaussian"))
  values$regime <- factor(values$regime, levels = c("homogeneous", "heterogeneous"))
  values$value <- values[[metric]]
  markers <- do.call(rbind, lapply(split(values, interaction(values$d, values$covariate, values$regime, drop = TRUE)), function(z) {
    z <- z[order(z$L), , drop = FALSE]
    start <- if (z$d[1] == 3L && nrow(z) > 1L) 2L else 1L
    z[seq.int(start, nrow(z), by = 2L), , drop = FALSE]
  }))
  markers$outline <- ifelse(markers$d %in% c(2L, 3L), .35, .5)
  groups <- interaction(values$d, values$covariate, values$regime, drop = TRUE)
  lines <- values[as.integer(table(groups)[groups]) > 1L, , drop = FALSE]
  colors <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#009E73", "4" = "#CC79A7")
  present <- as.character(sort(unique(values$d)))
  sizes <- c("1" = 1.7, "2" = 1.15, "3" = 1.15, "4" = 1.7)
  strokes <- c("1" = .5, "2" = .35, "3" = .35, "4" = .5)
  ylabel <- if (metric == "cost_L") expression(log[10](plain(Cost)(L))) else expression(log[10](plain(Cost)(0)-plain(Cost)(L)))
  ggplot2::ggplot(values, ggplot2::aes(L, value, colour = dimension, shape = dimension, fill = dimension)) +
    ggplot2::geom_line(data = lines, ggplot2::aes(group = dimension), linewidth = .75) +
    ggplot2::geom_point(data = markers, ggplot2::aes(size = dimension, stroke = I(outline))) +
    ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(values = c("1" = 16, "2" = 24, "3" = 22, "4" = 18)) +
    ggplot2::scale_size_manual(values = c("1" = 1.7, "2" = 1.15, "3" = 1.15, "4" = 1.7), guide = "none") +
    ggplot2::scale_x_log10(breaks = c(.1, 1, 10), labels = c("-1", "0", "1"),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = c(.02, .04))) +
    ggplot2::scale_y_log10(breaks = c(.001, .01, .1), labels = c("-3", "-2", "-1"),
      minor_breaks = NULL, expand = ggplot2::expansion(mult = c(.045, .08))) +
    ggplot2::facet_grid(covariate ~ regime, scales = "free_y", labeller = ggplot2::labeller(
      covariate = c(uniform = "Uniform", truncated_gaussian = "Truncated Gaussian"),
      regime = c(homogeneous = "Homogeneous", heterogeneous = "Heterogeneous"))) +
    ggplot2::labs(x = expression(log[10](L)), y = ylabel, colour = "d", shape = "d") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(
      size = unname(sizes[present]), stroke = unname(strokes[present]), fill = unname(colors[present])))) +
    mab_release_plot_theme()
}

run_cost_of_fairness <- function(options = mab_release_options("cost_of_fairness")) {
  output_dir <- options$output_dir
  if (options$worker_count != 1L) stop("Cost calculations use one coordinator; omit --worker-count")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("filelock", quietly = TRUE)) stop("Install filelock")
  lock <- filelock::lock(file.path(output_dir, ".cost.lock"), timeout = 0)
  if (is.null(lock)) return(invisible(NULL))
  on.exit(filelock::unlock(lock), add = TRUE)
  cfg <- mab_simulation_settings()
  L_grid <- exp(seq(log(.02), log(10), length.out = 81L))
  target <- file.path(output_dir, "cost_of_fairness.csv")
  if (!options$plot_only) {
    settings <- list(version = 1L, cfg = cfg, L_grid = L_grid,
      numerical_tolerance = MAB_NUMERICAL_TOLERANCE)
    settings_file <- file.path(output_dir, "settings.rds")
    previous <- mab_read_rds(settings_file)
    if (is.null(previous)) mab_save_rds(settings, settings_file) else if (!identical(previous, settings))
      stop("This output directory contains different cost settings; select a new --output directory")
    tasks <- expand.grid(d = 1:4, covariate = c("uniform", "truncated_gaussian"),
      regime = c("homogeneous", "heterogeneous"), L = L_grid, stringsAsFactors = FALSE)
    saved_path <- file.path(output_dir, "cost_values.rds")
    completed <- mab_read_rds(saved_path,
      default = if (file.exists(target)) read.csv(target) else data.frame())
    key <- function(z) paste(z$d, z$covariate, z$regime, sprintf("%.14g", z$L), sep = "|")
    done <- if (nrow(completed)) key(completed) else character()
    newly_completed <- 0L
    for (i in seq_len(nrow(tasks))) {
      task <- tasks[i, , drop = FALSE]
      if (key(task) %in% done) next
      if (newly_completed >= options$max_tasks) break
      completed <- rbind(completed, mab_ordinary_cost(as.list(task), cfg))
      mab_save_rds(completed, saved_path)
      write.csv(completed, target, row.names = FALSE)
      newly_completed <- newly_completed + 1L
      if (newly_completed %% 16L == 0L || i == nrow(tasks))
        message(sprintf("Cost of fairness: %d/%d values saved", nrow(completed), nrow(tasks)))
    }
  }
  if (!file.exists(target)) stop("No saved cost values; first run without --plot-only")
  values <- read.csv(target)
  for (metric in c("cost_L", "cost0_minus_cost_L")) {
    stem <- if (metric == "cost_L") "cost_of_fairness" else "cost_reduction"
    p <- mab_cost_plot(values, metric)
    ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), p, width = 12, height = 12*1230/1788)
    ggplot2::ggsave(file.path(output_dir, paste0(stem, ".png")), p, width = 12, height = 12*1230/1788, dpi = 180)
  }
  invisible(values)
}

if (!identical(getOption("mab.release.autorun", TRUE), FALSE)) run_cost_of_fairness()
