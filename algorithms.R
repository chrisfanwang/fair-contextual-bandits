# Fair bandit algorithms, policy optimisation and reward evaluation.

MAB_NUMERICAL_TOLERANCE <- 1e-6
MAB_CACHE_KEY_MAX_BYTES <- 9000L
.mab_geometry_loaded <- FALSE
MAB_SOURCE_FILE <- local({
  source_paths <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  file_argument <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  source_path <- if (length(source_paths)) tail(source_paths, 1L)[[1L]] else
    if (length(file_argument)) file_argument[[1L]] else "algorithms.R"
  if(!file.exists(source_path))source_path<-gsub("~+~"," ",source_path,fixed=TRUE)
  normalizePath(source_path, mustWork = TRUE)
})

check_required_packages <- function() {
  pkgs <- c("Rglpk", "slam")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }
}

clip01 <- function(z) {
  out <- pmin(1, pmax(0, z))
  
  z_dim <- dim(z)
  if (!is.null(z_dim)) {
    dim(out) <- z_dim
    dimnames(out) <- dimnames(z)
  }
  
  out
}

add_intercept <- function(X) {
  X <- as.matrix(X)
  cbind(Intercept = 1, X)
}

normalize_nu <- function(nu) {
  if (is.character(nu)) {
    if (tolower(nu) %in% c("inf", "infty", "infinity")) {
      return(Inf)
    }
    nu <- as.numeric(nu)
  }
  nu <- as.numeric(nu)
  if (length(nu) != 1L || is.na(nu) || nu < 1){
    stop("nu must be a scalar in [1, Inf].", call. = FALSE)
  }
  nu
}

norm_nu_rows_to_point <- function(X, h, nu) {
  X <- as.matrix(X)
  h <- as.numeric(h)
  nu <- normalize_nu(nu)
  
  if (ncol(X) != length(h))
    stop("X and h have incompatible dimensions.", call. = FALSE)
  if (nrow(X) == 0L) return(numeric(0))
  
  out <- numeric(nrow(X))
  
  if (is.infinite(nu)) {
    for (k in seq_len(ncol(X)))
      out <- pmax(out, abs(X[, k] - h[k]))
    return(out)
  }
  
  for (k in seq_len(ncol(X)))
    out <- out + abs(X[, k] - h[k])^nu
  
  out^(1 / nu)
}

pairwise_norm_nu <- function(X, Y, nu) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  nu <- normalize_nu(nu)
  
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.", call. = FALSE)
  }
  if (nrow(X) == 0L || nrow(Y) == 0L) {
    return(matrix(0, nrow(X), nrow(Y)))
  }
  
  out <- matrix(0, nrow(X), nrow(Y))
  
  if (is.infinite(nu)) {
    for (k in seq_len(ncol(X))) {
      out <- pmax(out, abs(outer(X[, k], Y[, k], "-")))
    }
    return(out)
  }
  
  for (k in seq_len(ncol(X))) {
    out <- out + abs(outer(X[, k], Y[, k], "-"))^nu
  }
  out^(1 / nu)
}

unique_rows_tol <- function(X, digits = 10L) {
  X <- as.matrix(X)
  if (nrow(X) == 0L) {
    return(X)
  }
  key <- apply(round(X, digits = digits), 1L, paste, collapse = "|")
  X[!duplicated(key), , drop = FALSE]
}

make_constant_policy <- function(prob, d, label = "constant") {
  prob <- as.numeric(prob)
  d <- as.integer(d)
  if (length(prob) != 1L || !is.finite(prob)) {
    stop("prob must be a finite scalar.", call. = FALSE)
  }
  if (d <= 0L) {
    stop("d must be a positive integer.", call. = FALSE)
  }
  p0 <- clip01(prob)
  
  eval_one <- local({
    p_loc <- p0
    d_loc <- d
    function(x) {
      x <- as.numeric(x)
      if (length(x) != d_loc) {
        stop(sprintf("x must have length %d.", d_loc), call. = FALSE)
      }
      p_loc
    }
  })
  
  eval_batch <- local({
    p_loc <- p0
    d_loc <- d
    function(Xnew) {
      Xnew <- as.matrix(Xnew)
      if (ncol(Xnew) != d_loc) {
        stop(sprintf("Xnew must have %d columns.", d_loc), call. = FALSE)
      }
      rep(p_loc, nrow(Xnew))
    }
  })
  
  structure(
    list(
      type = "constant",
      d = d,
      params = list(prob = p0),
      meta = list(label = label, fairness_guaranteed = TRUE, policy_version = "constant_v2",
        canonical_eval = "constant", n_boundary_samples = 0L, all_pairs_certified = TRUE,
        max_pair_violation = 0, initial_max_pair_violation = 0),
      eval_one = eval_one,
      eval_batch = eval_batch
    ),
    class = "estimated_policy"
  )
}

make_d1_boundary_policy <- function(theta, H_ref, q_ref, L, nu = 1,
                                           label = "boundary_canonical",
                                           meta = list()) {
  theta <- as.numeric(theta)
  H_ref <- as.matrix(H_ref)
  q_ref <- clip01(as.numeric(q_ref))
  L <- as.numeric(L)
  nu <- normalize_nu(nu)
  d <- ncol(H_ref)
  if (d != 1L || nrow(H_ref) != 1L)
    stop("The one-dimensional policy requires a single boundary point.", call. = FALSE)
  
  if (length(theta) != d + 1L)
    stop("theta must have length d + 1, where d = ncol(H_ref).", call. = FALSE)
  if (nrow(H_ref) == 0L)
    stop("H_ref must contain at least one boundary reference point.", call. = FALSE)
  if (length(q_ref) != nrow(H_ref))
    stop("q_ref must have length nrow(H_ref).", call. = FALSE)
  if (length(L) != 1L || !is.finite(L) || L < 0)
    stop("L must be a nonnegative finite scalar.", call. = FALSE)
  
  eval_one <- local({
    theta_loc <- theta
    H_loc <- H_ref
    q_loc <- q_ref
    L_loc <- L
    nu_loc <- nu
    d_loc <- d
    
    function(x) {
      x <- as.numeric(x)
      if (length(x) != d_loc)
        stop(sprintf("x must have length %d.", d_loc), call. = FALSE)
      
      Xnew <- matrix(x, nrow = 1L)
      adv <- drop(add_intercept(Xnew) %*% theta_loc)
      D <- norm_nu_rows_to_point(H_loc, x, nu = nu_loc)
      upper_env <- min(q_loc + L_loc * D)
      lower_env <- max(q_loc - L_loc * D)
      
      if (adv >= 0) min(1, upper_env) else max(0, lower_env)
    }
  })
  
  eval_batch <- local({
    theta_loc <- theta
    H_loc <- H_ref
    q_loc <- q_ref
    L_loc <- L
    nu_loc <- nu
    d_loc <- d
    batch_size <- 256L
    
    function(Xnew) {
      Xnew <- as.matrix(Xnew)
      if (ncol(Xnew) != d_loc)
        stop(sprintf("Xnew must have %d columns.", d_loc), call. = FALSE)
      
      n <- nrow(Xnew)
      if (n == 0L) return(numeric(0))
      result <- numeric(n)
      
      for (start in seq.int(1L, n, by = batch_size)) {
        end <- min(start + batch_size - 1L, n)
        idx <- start:end
        X_block <- Xnew[idx, , drop = FALSE]
        adv <- drop(add_intercept(X_block) %*% theta_loc)
        D <- pairwise_norm_nu(X_block, H_loc, nu = nu_loc)
        
        upper_env <- apply(
          sweep(L_loc * D, 2L, q_loc, "+"), 1L, min
        )
        lower_env <- apply(
          sweep(-L_loc * D, 2L, q_loc, "+"), 1L, max
        )
        
        result[idx] <- ifelse(
          adv >= 0,
          pmin(1, upper_env),
          pmax(0, lower_env)
        )
      }
      
      result
    }
  })
  
  structure(
    list(
      type = "boundary_canonical",
      d = d,
      params = list(theta = theta, H_ref = H_ref, q_ref = q_ref,
                    L = L, nu = nu),
      meta = utils::modifyList(meta, list(
        label = label, n_boundary = nrow(H_ref),
        canonical_eval = if (d == 1L) "exact_single_boundary_point" else "legacy_sampled_boundary",
        fairness_guaranteed = d == 1L,
        policy_version = "single_boundary_v1"
      )),
      eval_one = eval_one,
      eval_batch = eval_batch
    ),
    class = "estimated_policy"
  )
}

make_d2_l1_trace_policy <- function(theta, a, v, knots, values, L,
                                    label = "continuous_boundary", meta = list(),
                                    trace_definition = "piecewise_linear_boundary_trace") {
  knots <- as.numeric(knots)
  values <- as.numeric(values)
  if (length(knots) < 2L || length(knots) != length(values) ||
      any(!is.finite(c(knots, values))) || any(diff(knots) <= 0) ||
      knots[1] != 0 || tail(knots, 1L) != 1) {
    stop("Trace knots must strictly increase from 0 to 1.", call. = FALSE)
  }
  if (any(values < -MAB_NUMERICAL_TOLERANCE | values > 1 + MAB_NUMERICAL_TOLERANCE) ||
      any(abs(diff(values)) > L * sum(abs(v)) * diff(knots) + MAB_NUMERICAL_TOLERANCE)) {
    stop("The continuous boundary trace must be [0,1]-valued and L-Lipschitz.", call. = FALSE)
  }
  values <- clip01(values)
  H <- sweep(outer(knots, v), 2L, a, "+")
  trace <- function(t) stats::approx(knots, values, xout = t, rule = 2L)$y
  eval_batch <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    if (ncol(Xnew) != 2L || any(!is.finite(Xnew)))
      stop("Xnew must be a finite matrix with two columns.", call. = FALSE)
    n <- nrow(Xnew)
    out <- numeric(n)
    if (!n) return(out)
    for (first in seq.int(1L, n, by = 256L)) {
      rows <- first:min(n, first + 255L)
      X <- Xnew[rows, , drop = FALSE]
      D <- pairwise_norm_nu(X, H, 1)
      upper <- apply(sweep(L * D, 2L, values, "+"), 1L, min)
      lower <- apply(sweep(-L * D, 2L, values, "+"), 1L, max)
      # On each interval the two objectives are affine. Their only additional
      # knots are where a coordinate of h(t) crosses the query coordinate.
      for (j in which(v != 0)) {
        t <- clip01((X[, j] - a[j]) / v[j])
        h <- sweep(outer(t, v), 2L, a, "+")
        q <- trace(t)
        distance <- rowSums(abs(X - h))
        upper <- pmin(upper, q + L * distance)
        lower <- pmax(lower, q - L * distance)
      }
      out[rows] <- ifelse(hyperplane_value(X, theta) >= 0,
                          pmin(1, upper), pmax(0, lower))
    }
    out
  }
  eval_one <- function(x) eval_batch(matrix(as.numeric(x), nrow = 1L))[1L]
  structure(list(
    type = "boundary_canonical_continuous_d2_l1", d = 2L,
    params = list(theta = theta, H_ref = H, q_ref = values, L = L, nu = 1,
                  trace_origin = a, trace_direction = v, trace_knots = knots),
    meta = utils::modifyList(meta, list(
      label = label, canonical_eval = "exact_continuous_d2_l1",
      policy_version = "continuous_d2_l1_v1", fairness_guaranteed = TRUE,
      n_boundary = length(knots), n_boundary_samples = 0L,
      n_boundary_breakpoints = length(knots),
      boundary_reference_type = "complete_continuous_trace_knots",
      q_definition = trace_definition, uses_boundary_samples = FALSE,
      uses_boundary_trace = TRUE, uses_boundary_canonicalization = TRUE,
      oracle_uses_boundary_samples = FALSE,
      oracle_boundary_sample_source = "none_continuous_boundary",
      policy_form = "continuous_boundary_canonical",
      oracle_policy_form = "continuous_boundary_canonical"
    )), eval_one = eval_one, eval_batch = eval_batch, eval_trace = trace
  ), class = "estimated_policy")
}

make_exact_d2_l1_boundary_policy <- function(theta, X_ref, z_ref, L,
                                             label = "exact_boundary", meta = list()) {
  vertices <- box_hyperplane_points(theta)
  if (nrow(vertices) != 2L)
    stop("The exact d=2 construction requires a nondegenerate boundary segment.", call. = FALSE)
  a <- vertices[1L, ]
  v <- vertices[2L, ] - a
  knots <- c(0, 1)
  values <- c(1, 1)
  for (i in seq_len(nrow(X_ref))) {
    cone <- function(t) {
      h <- sweep(outer(t, v), 2L, a, "+")
      z_ref[i] + L * rowSums(abs(sweep(h, 2L, X_ref[i, ], "-")))
    }
    active <- which(v != 0)
    breaks <- (X_ref[i, active] - a[active]) / v[active]
    t <- sort(unique(c(knots, breaks[breaks > 0 & breaks < 1])))
    delta <- stats::approx(knots, values, xout = t)$y - cone(t)
    left <- head(delta, -1L)
    right <- tail(delta, -1L)
    crossing <- which((left < 0 & right > 0) | (left > 0 & right < 0))
    if (length(crossing)) {
      extra <- t[crossing] + (t[crossing + 1L] - t[crossing]) *
        delta[crossing] / (delta[crossing] - delta[crossing + 1L])
      t <- sort(unique(c(t, extra)))
    }
    values <- pmin(stats::approx(knots, values, xout = t)$y, cone(t))
    knots <- t
  }
  make_d2_l1_trace_policy(theta, a, v, knots, values, L, label, meta,
                         "exact_upper_mcshane_trace_of_whole_space_lp")
}

boundary_l1_two_hop_distances <- function(x, X_ref, theta) {
  d <- ncol(X_ref)
  gamma <- theta[-1L]
  scale <- max(abs(gamma))
  if (!is.finite(scale) || scale == 0)
    stop("A nonzero boundary normal is required.", call. = FALSE)
  gamma <- gamma / scale
  intercept <- theta[1L] / scale
  low <- high <- X_ref
  lower_score <- upper_score <- rep(intercept, nrow(X_ref))
  direct <- numeric(nrow(X_ref))
  for (j in seq_len(d)) {
    low[, j] <- pmin(x[j], X_ref[, j])
    high[, j] <- pmax(x[j], X_ref[, j])
    direct <- direct + high[, j] - low[, j]
    lower_score <- lower_score + gamma[j] * if (gamma[j] >= 0) low[, j] else high[, j]
    upper_score <- upper_score + gamma[j] * if (gamma[j] >= 0) high[, j] else low[, j]
  }
  decrease <- lower_score > 0
  remaining <- pmax(lower_score, -upper_score, 0)
  extra <- numeric(nrow(X_ref))
  for (j in order(abs(gamma), decreasing = TRUE)) {
    if (gamma[j] == 0) next
    capacity_down <- if (gamma[j] > 0) low[, j] else 1 - high[, j]
    capacity_up <- if (gamma[j] > 0) 1 - high[, j] else low[, j]
    capacity <- ifelse(decrease, capacity_down, capacity_up)
    moved <- pmin(remaining / abs(gamma[j]), capacity)
    extra <- extra + moved
    remaining <- pmax(0, remaining - abs(gamma[j]) * moved)
  }
  if (any(remaining > MAB_NUMERICAL_TOLERANCE))
    stop("The boundary does not intersect the unit cube.", call. = FALSE)
  direct + 2 * extra
}

rank1_update <- function(Vinv, b, ell, x, y) {
  x <- as.numeric(x)
  z <- drop(Vinv %*% x)
  den <- as.numeric(1 + crossprod(x, z))
  
  list(
    Vinv = Vinv - (z %o% z) / den,
    b = b + x * y,
    ell = ell + log(den)
  )
}

ucb_value_matrix <- function(Xtilde, beta_hat, Vinv, alpha) {
  Xtilde <- as.matrix(Xtilde)
  mu_hat <- drop(Xtilde %*% beta_hat)
  quad <- rowSums((Xtilde %*% Vinv) * Xtilde)
  mu_hat + alpha * sqrt(pmax(quad, 0))
}

standard_ucb_policy_matrix <- function(Xnew,
                                       beta1_hat,
                                       beta2_hat,
                                       Vinv1,
                                       Vinv2,
                                       alpha) {
  Xtilde <- add_intercept(Xnew)
  U1 <- ucb_value_matrix(Xtilde, beta1_hat, Vinv1, alpha)
  U2 <- ucb_value_matrix(Xtilde, beta2_hat, Vinv2, alpha)
  as.numeric(U1 >= U2)
}

hyperplane_value <- function(X, theta) {
  X <- as.matrix(X)
  theta <- as.numeric(theta)
  drop(add_intercept(X) %*% theta)
}

# Geometric degeneracy thresholds are distinct from accepted solver residuals.
hyperplane_intersects_unit_box <- function(theta, tol = 1e-10) {
  theta <- as.numeric(theta)
  slopes <- theta[-1]
  min_val <- theta[1] + sum(ifelse(slopes < 0, slopes, 0))
  max_val <- theta[1] + sum(ifelse(slopes > 0, slopes, 0))
  min_val <= tol && max_val >= -tol
}

binary_grid <- function(n) {
  n <- as.integer(n)
  if (n == 0L) {
    return(matrix(numeric(0), nrow = 1L, ncol = 0L))
  }
  as.matrix(do.call(expand.grid, rep(list(c(0, 1)), n)))
}

box_hyperplane_points <- function(theta, tol = 1e-10) {
  theta <- as.numeric(theta)
  slopes <- theta[-1]
  d <- length(slopes)
  out <- list()
  ptr <- 1L
  
  if (d == 1L) {
    if (abs(slopes[1]) > tol) {
      x0 <- -theta[1] / slopes[1]
      if (is.finite(x0) && x0 >= -tol && x0 <= 1 + tol) {
        return(matrix(clip01(x0), nrow = 1L, ncol = 1L))
      }
    }
    return(matrix(numeric(0), nrow = 0L, ncol = 1L))
  }
  
  for (free_j in seq_len(d)) {
    fixed_idx <- setdiff(seq_len(d), free_j)
    G <- binary_grid(length(fixed_idx))
    for (r in seq_len(nrow(G))) {
      h <- numeric(d)
      if (length(fixed_idx) > 0L) {
        h[fixed_idx] <- as.numeric(G[r, ])
      }
      if (abs(slopes[free_j]) <= tol) {
        next
      }
      h[free_j] <- -(theta[1] + sum(slopes[fixed_idx] * h[fixed_idx])) / slopes[free_j]
      if (is.finite(h[free_j]) && h[free_j] >= -tol && h[free_j] <= 1 + tol) {
        out[[ptr]] <- clip01(h)
        ptr <- ptr + 1L
      }
    }
  }
  
  if (length(out) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  unique_rows_tol(do.call(rbind, out))
}

mcshane_upper_extension_batch <- function(Xnew,
                                          X_ref,
                                          z_ref,
                                          L,
                                          nu = 1,
                                          batch_size = 5000L) {
  Xnew <- as.matrix(Xnew)
  X_ref <- as.matrix(X_ref)
  z_ref <- clip01(as.numeric(z_ref))
  L <- as.numeric(L)
  nu <- normalize_nu(nu)
  if (nrow(Xnew) == 0L) return(numeric(0))
  if (nrow(X_ref) == 0L) stop("X_ref must be non-empty.", call. = FALSE)
  if (ncol(Xnew) != ncol(X_ref)) stop("Xnew and X_ref have incompatible dimensions.", call. = FALSE)
  if (length(z_ref) != nrow(X_ref)) stop("z_ref must have length nrow(X_ref).", call. = FALSE)
  
  batch_size <- as.integer(batch_size)
  if (!is.finite(batch_size) || batch_size <= 0L) batch_size <- nrow(Xnew)
  out <- numeric(nrow(Xnew))
  for (start in seq.int(1L, nrow(Xnew), by = batch_size)) {
    end <- min(start + batch_size - 1L, nrow(Xnew))
    idx <- start:end
    vals <- rep(Inf, length(idx))
    X_batch <- Xnew[idx, , drop = FALSE]
    for (j in seq_len(nrow(X_ref))) {
      vals <- pmin(vals, z_ref[j] + L * norm_nu_rows_to_point(X_batch, X_ref[j, ], nu = nu))
    }
    out[idx] <- pmin(1, vals)
  }
  clip01(out)
}

mab_validate_neighbor_k <- function(k) {
  if(!is.numeric(k)||length(k)!=1L||!is.finite(k)||k<1||k>.Machine$integer.max||k!=floor(k))
    stop("lp_neighbor_k must be a positive integer.")
  as.integer(k)
}

.mab_geometry_source <- r"---(#include <Rcpp.h>
#include <algorithm>
#include <queue>
#include <numeric>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
struct Node {
  int begin,end,left=-1,right=-1,imax;
  std::vector<double> lo,hi;
  double zlo=0,zhi=0;
};
class Tree {
public:
  NumericMatrix x; NumericVector z; int n,d;
  std::vector<int> ids; std::vector<Node> nodes;
  Tree(NumericMatrix xx, NumericVector zz):x(xx),z(zz),n(xx.nrow()),d(xx.ncol()) {
    ids.resize(n);std::iota(ids.begin(),ids.end(),0);nodes.reserve(2*n);build(0,n);
  }
  int build(int a,int b) {
    Node node;node.begin=a;node.end=b;node.imax=-1;
    node.lo.assign(d,R_PosInf);node.hi.assign(d,R_NegInf);
    node.zlo=R_PosInf;node.zhi=R_NegInf;
    for(int k=a;k<b;k++) {
      int i=ids[k];node.imax=std::max(node.imax,i);
      for(int j=0;j<d;j++) {node.lo[j]=std::min(node.lo[j],x(i,j));node.hi[j]=std::max(node.hi[j],x(i,j));}
      if(z.size()) {node.zlo=std::min(node.zlo,z[i]);node.zhi=std::max(node.zhi,z[i]);}
    }
    int id=nodes.size();nodes.push_back(node);
    if(b-a>24) {
      int axis=0;for(int j=1;j<d;j++)if(node.hi[j]-node.lo[j]>node.hi[axis]-node.lo[axis])axis=j;
      int mid=(a+b)/2;
      std::nth_element(ids.begin()+a,ids.begin()+mid,ids.begin()+b,[&](int i,int k){return x(i,axis)<x(k,axis)||(x(i,axis)==x(k,axis)&&i<k);});
      int l=build(a,mid),r=build(mid,b);nodes[id].left=l;nodes[id].right=r;
    }
    return id;
  }
  double dist(int i,int j)const {double out=0;for(int k=0;k<d;k++)out+=std::abs(x(i,k)-x(j,k));return out;}
  double lower(int i,int id)const {
    double out=0;const Node &v=nodes[id];
    for(int j=0;j<d;j++)out+=std::max(std::max(v.lo[j]-x(i,j),x(i,j)-v.hi[j]),0.0);
    return out;
  }
  typedef std::priority_queue<std::pair<double,int>> Heap;
  void knn(int i,int id,int k,Heap &h)const {
    const Node &v=nodes[id];
    if((int)h.size()==k && lower(i,id)>h.top().first+1e-12)return;
    if(v.left<0) {
      for(int p=v.begin;p<v.end;p++) {
        int j=ids[p];if(j==i)continue;
        std::pair<double,int>a(dist(i,j),j);
        if((int)h.size()<k)h.push(a);else if(a<h.top()){h.pop();h.push(a);}
      }
    } else {
      int l=v.left,r=v.right;if(lower(i,l)>lower(i,r))std::swap(l,r);
      knn(i,l,k,h);knn(i,r,k,h);
    }
  }
  bool dominated(int i,int id,double L)const {
    const Node &v=nodes[id];
    if(v.zlo+L*lower(i,id)>z[i]+1e-12)return false;
    if(v.left<0) {
      for(int p=v.begin;p<v.end;p++) {int j=ids[p];if(j==i)continue;double dd=dist(i,j);
        if((dd>0||j<i) && z[i]>=z[j]+L*dd)return true;
      }
      return false;
    }
    int l=v.left,r=v.right;
    if(nodes[l].zlo+L*lower(i,l)>nodes[r].zlo+L*lower(i,r))std::swap(l,r);
    return dominated(i,l,L)||dominated(i,r,L);
  }
  // Heap stores (-violation,j), so the worst of the top 20 is first.
  void cuts(int i,int id,double L,double tol,int k,Heap &h)const {
    const Node &v=nodes[id];if(v.imax<=i)return;
    double bound=std::max(z[i]-v.zlo,v.zhi-z[i])-L*lower(i,id);
    double threshold=h.size()==(size_t)k?std::max(tol,-h.top().first):tol;
    if(bound<threshold-1e-12)return;
    if(v.left<0) {
      for(int p=v.begin;p<v.end;p++) {
        int j=ids[p];if(j<=i)continue;
        double violation=std::abs(z[i]-z[j])-L*dist(i,j);
        if(violation<=tol)continue;
        std::pair<double,int>a(-violation,j);
        if((int)h.size()<k)h.push(a);else if(a<h.top()){h.pop();h.push(a);}
      }
    } else {cuts(i,v.left,L,tol,k,h);cuts(i,v.right,L,tol,k,h);}
  }
};
// [[Rcpp::export]]
DataFrame exact_l1_knn_cpp(NumericMatrix x,int k) {
  Tree tree(x,NumericVector(0));int n=x.nrow();k=std::min(k,n-1);
  std::vector<std::pair<int,int>> edges;edges.reserve((size_t)n*k);
  for(int i=0;i<n;i++) {
    if(i%256==0)checkUserInterrupt();Tree::Heap h;tree.knn(i,0,k,h);
    while(!h.empty()) {int j=h.top().second;h.pop();edges.emplace_back(std::min(i,j),std::max(i,j));}
  }
  std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
  IntegerVector a(edges.size()),b(edges.size());NumericVector dd(edges.size());
  for(size_t p=0;p<edges.size();p++){a[p]=edges[p].first+1;b[p]=edges[p].second+1;dd[p]=tree.dist(a[p]-1,b[p]-1);}
  return DataFrame::create(_["i"]=a,_["j"]=b,_["distance"]=dd);
}
// [[Rcpp::export]]
LogicalVector exact_l1_dominated_cpp(NumericMatrix x,NumericVector z,double L) {
  Tree tree(x,z);LogicalVector remove(x.nrow());
  for(int i=0;i<x.nrow();i++){if(i%256==0)checkUserInterrupt();remove[i]=z[i]>=1||tree.dominated(i,0,L);}
  return remove;
}
// [[Rcpp::export]]
List exact_l1_cuts_cpp(NumericMatrix x,NumericVector z,double L,double tolerance) {
  Tree tree(x,z);std::vector<int> a,b;std::vector<double> dd;double maximum=tolerance;
  for(int i=0;i<x.nrow();i++) {
    if(i%256==0)checkUserInterrupt();Tree::Heap h;tree.cuts(i,0,L,tolerance,20,h);
    std::vector<std::pair<double,int>> selected;
    while(!h.empty()){selected.push_back(h.top());h.pop();}
    std::sort(selected.begin(),selected.end());
    for(auto v:selected){a.push_back(i+1);b.push_back(v.second+1);dd.push_back(tree.dist(i,v.second));maximum=std::max(maximum,-v.first);}
  }
  return List::create(_["edges"]=DataFrame::create(_["i"]=a,_["j"]=b,_["distance"]=dd),
    _["certified_max_violation_upper_bound"]=maximum);
}

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif
#ifdef __linux__
#include <fstream>
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>
#endif
// Process lifetime high-water RSS, not a per-job isolated peak.
// [[Rcpp::export]]
double mab_process_peak_rss_mb() {
#if defined(__APPLE__) || defined(__linux__)
  struct rusage r;
  if (getrusage(RUSAGE_SELF, &r) != 0) return NA_REAL;
#ifdef __APPLE__
  return ((double)r.ru_maxrss)/(1024.0*1024.0);
#else
  return ((double)r.ru_maxrss)/1024.0;
#endif
#else
  return NA_REAL;
#endif
}
// [[Rcpp::export]]
double mab_current_rss_mb() {
#ifdef __APPLE__
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
  if(task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count)==KERN_SUCCESS)
    return ((double)info.resident_size)/(1024.0*1024.0);
#elif defined(__linux__)
  std::ifstream input("/proc/self/statm");
  long total=0, resident=0;
  if(input >> total >> resident)
    return ((double)resident)*sysconf(_SC_PAGESIZE)/(1024.0*1024.0);
#endif
  return NA_REAL;
}
)---"

mab_geometry_init <- function() {
  if (!.mab_geometry_loaded) {
    if (!requireNamespace("Rcpp",quietly=TRUE)) stop("Install Rcpp for exact geometry.")
    Rcpp::sourceCpp(code=.mab_geometry_source,
      env=environment(mab_geometry_init),
      cacheDir=getOption("mab.rcpp_cache",file.path(tempdir(),"mab_rcpp")),showOutput=FALSE)
    .mab_geometry_loaded <<- TRUE
  }
  invisible(TRUE)
}

sweep_dominated_blocked <- function(X,z,L,block=64L) {
  mab_geometry_init(); exact_l1_dominated_cpp(as.matrix(X),as.numeric(z),L)
}

solve_empirical_lipschitz_lp <- function(X_obj,theta,weights,L,nu=1,
    lp_neighbor_k=20L,zero_obj_tol=0,initial_edges=NULL) {
  check_required_packages(); X_obj<-as.matrix(X_obj);n<-nrow(X_obj);d<-ncol(X_obj)
  stopifnot(n>0,length(theta)==d+1L,length(weights)==n,all(is.finite(c(X_obj,theta,weights,L))),
    all(weights>=0),sum(weights)>0,L>=0,nu==1)
  weights<-weights/sum(weights);adv<-drop(add_intercept(X_obj)%*%theta);obj<-weights*adv
  initial_edge_count<-0L
  answer<-function(z,E,rounds,maximum) list(z_ref=z,p_obj=z,objective_value=sum(obj*z),
    n_constraints=2*nrow(E),n_lipschitz_edges=nrow(E),lp_constraint_mode="full_history_exact_separation",
    all_pairs_certified=TRUE,max_pair_violation=maximum,initial_max_pair_violation=NA_real_,
    constraint_generation_rounds=rounds,n_initial_edges=initial_edge_count,adv=adv)
  empty<-data.frame(i=integer(),j=integer(),distance=numeric())
  if (L==0 || n==1 || all(obj>=0) || all(obj<=0)) {
    # Global sign objective or L=0 has an exact constant optimizer.
    z<-rep(if(sum(obj)>0)1 else if(sum(obj)<0)0 else .5,n)
    return(answer(z,empty,0L,0))
  }
  if(d==1L) {
    z<-mab_d1_optimizer(X_obj[,1],theta,obj,L)
    ord<-order(X_obj[,1]);v<-max(0,abs(diff(z[ord]))-L*diff(X_obj[ord,1]))
    # Canonical clipped affine form is globally L-Lipschitz analytically.
    if(v>MAB_NUMERICAL_TOLERANCE)stop("One-dimensional optimizer failed its fairness certificate")
    ans<-answer(z,empty,0L,v);ans$lp_constraint_mode<-"exact_d1_breakpoint_sweep";return(ans)
  }
  mab_geometry_init()
  if(!is.null(initial_edges)) E<-initial_edges else if(d==1L) {
    ord<-order(X_obj[,1],seq_len(n));i<-head(ord,-1);j<-tail(ord,-1)
    E<-data.frame(i=pmin(i,j),j=pmax(i,j),distance=abs(X_obj[i,1]-X_obj[j,1]))
  } else if(is.null(lp_neighbor_k)||identical(lp_neighbor_k,"full")||is.na(lp_neighbor_k)||is.infinite(lp_neighbor_k)) {
    if(n>10000)stop("Dense reference LP refused above 10000 points; use exact separation.")
    ij<-which(upper.tri(matrix(FALSE,n,n)),arr.ind=TRUE)
    E<-data.frame(i=ij[,1],j=ij[,2],distance=rowSums(abs(X_obj[ij[,1],,drop=FALSE]-X_obj[ij[,2],,drop=FALSE])))
  } else E<-exact_l1_knn_cpp(X_obj,as.integer(min(n-1,max(1,lp_neighbor_k))))
  initial_edge_count<-nrow(E)
  tolerance<-MAB_NUMERICAL_TOLERANCE;last_n<-0L
  for(rounds in 0:100) {
    m<-nrow(E);if(m>getOption("mab.max_lp_edges",2000000L))stop("Exact LP edge memory budget exceeded; history was not truncated")
    rr<-seq_len(m)
    # Scale short-distance constraints so GLPK's absolute row tolerance cannot
    # erase the fairness constraint. Scale objective without changing argmax.
    scale<-1/pmax(L*E$distance,1e-7)
    A<-slam::simple_triplet_matrix(i=c(rr,rr,m+rr,m+rr),j=c(E$i,E$j,E$i,E$j),
      v=c(scale,-scale,-scale,scale),nrow=2L*m,ncol=n)
    sol<-Rglpk::Rglpk_solve_LP(obj/max(abs(obj)),A,rep("<=",2*m),rep(L*E$distance*scale,2),
      bounds=list(lower=list(ind=seq_len(n),val=rep(0,n)),upper=list(ind=seq_len(n),val=rep(1,n))),
      max=TRUE,control=list(presolve=TRUE,tm_limit=getOption("mab.lp_time_ms",0L)))
    if(length(sol$status)!=1L || is.na(sol$status) || sol$status!=0L)
      stop("Finite Lipschitz LP did not certify optimality: status ",sol$status)
    z<-as.numeric(sol$solution)
    if(length(z)!=n || any(!is.finite(z)))stop("LP returned invalid probabilities")
    if(any(z< -tolerance|z>1+tolerance))stop("LP box feasibility failed")
    z<-clip01(z)
    sep<-exact_l1_cuts_cpp(X_obj,z,L,tolerance)
    maximum<-sep$certified_max_violation_upper_bound
    if(maximum<=tolerance)return(answer(z,E,rounds,maximum))
    E<-unique(rbind(E,sep$edges))
    if(nrow(E)==last_n)stop("All-pair feasibility not certified despite scaled rows: ",maximum)
    last_n<-nrow(E)
  }
  stop("Constraint generation did not certify all-pair feasibility after 100 rounds")
}

mab_policy_from_lp <- function(theta,X,z,L,nu=1,meta=list()) {
  d<-ncol(X);stopifnot(nu==1)
  lo<-theta[1]+sum(pmin(theta[-1],0));hi<-theta[1]+sum(pmax(theta[-1],0))
  if(all(theta==0))p<-make_constant_policy(.5,d) else if(lo>=0||hi<=0)
    p<-make_constant_policy(as.numeric(lo>=0),d) else if(d==1L) {
      H<-matrix(-theta[1]/theta[2],ncol=1)
      q<-mcshane_upper_extension_batch(H,X,z,L,nu)
      p<-make_d1_boundary_policy(theta,H,q,L,nu)
    } else if(d==2L)p<-make_exact_d2_l1_boundary_policy(theta,X,z,L,meta=meta) else
      p<-make_l1_boundary_tight_policy(theta,X,z,L,meta=meta,
        time_limit_ms=getOption("mab.boundary_time_ms",0L))
  p$meta<-utils::modifyList(p$meta,meta);p
}

mab_d1_optimizer <- function(x,theta,obj,L) {
  h<- -theta[1]/theta[2];offset<-sign(theta[2])*L*(x-h)
  q<-c(-offset,1-offset);delta<-c(obj,-obj)
  active<-q>0 & q<1;q<-q[active];delta<-delta[active]
  events<-if(length(q))tapply(delta,q,sum) else numeric()
  knots<-if(length(events))as.numeric(names(events)) else numeric()
  ord<-order(knots);knots<-knots[ord];events<-events[ord]
  slope<-sum(obj[offset>=0 & offset<1]);value<-sum(obj*clip01(offset))
  best<-value;bestq<-0;prev<-0
  for(j in seq_along(knots)) {
    value<-value+slope*(knots[j]-prev)
    if(value>best){best<-value;bestq<-knots[j]}
    slope<-slope+events[j];prev<-knots[j]
  }
  value<-value+slope*(1-prev)
  if(value>best)bestq<-1
  clip01(bestq+offset)
}

make_l1_boundary_tight_policy <- function(theta, X_ref, z_ref, L,
    label="continuous_reference", meta=list(), time_limit_ms=getOption("mab.boundary_time_ms",0L)) {
  X_ref<-as.matrix(X_ref);z_ref<-as.numeric(z_ref);d<-ncol(X_ref)
  n_before<-nrow(X_ref)
  stopifnot(length(theta)==d+1L,length(z_ref)==n_before,
    all(is.finite(c(theta,X_ref,z_ref,L))),all(X_ref>=0 & X_ref<=1),
    all(z_ref>=0 & z_ref<=1),L>=0)
  if(L==0)return(make_constant_policy(min(z_ref),d,label))
  stopifnot(hyperplane_intersects_unit_box(theta),max(abs(theta[-1]))>0)
  # Remove reference points dominated by another Lipschitz cone.
  remove<-sweep_dominated_blocked(X_ref,z_ref,L)
  keep<-which(!remove);if(!length(keep))keep<-1L
  X_ref<-X_ref[keep,,drop=FALSE];z_ref<-z_ref[keep];n<-length(keep)
  lo<-rep(0,d);hi<-rep(1,d);normal<-theta[-1]
  for(j in which(normal!=0)) {
    other<-normal[-j]
    limits<-sort((-theta[1]-c(sum(pmax(other,0)),sum(pmin(other,0))))/normal[j])
    # Widen, never shrink, to cover floating-point rounding in the division.
    margin<-1e-12*max(1,abs(limits))
    lo[j]<-max(0,limits[1]-margin);hi[j]<-min(1,limits[2]+margin)
  }
  stopifnot(all(lo<=hi))
  cache<-new.env(hash=TRUE,parent=emptyenv())
  stats<-new.env();stats$solves<-0L;stats$queries<-0L;stats$max_gap<-0
  stats$max_active<-0L;stats$total_active<-0L;stats$early_zero<-0L
  stats$numerical_retries<-0L;stats$certified_fallbacks<-0L
  model<-function(active) {
    key<-paste(active,collapse=",")
    # R limits environment-name length. Long keys bypass both lookup and storage.
    cacheable<-nchar(key,type="bytes")<=MAB_CACHE_KEY_MAX_BYTES
    if(cacheable && exists(key,cache,inherits=FALSE))return(get(key,cache))
    r_col<-d+1L;w_cols<-d+1L+seq_len(d);k<-2L*d+1L
    lower<-c(lo,0,rep(0,d));upper<-c(hi,1,rep(1,d));types<-rep("C",k)
    ii<-jj<-integer();vv<-rhs<-numeric();directions<-character()
    add<-function(cols,vals,bound,dir="<=") {
      row<-length(rhs)+1L;ii<<-c(ii,rep(row,length(cols)));jj<<-c(jj,cols)
      vv<<-c(vv,vals);rhs<<-c(rhs,bound);directions<<-c(directions,dir);row
    }
    variable<-function(up,type="C") {
      k<<-k+1L;lower<<-c(lower,0);upper<<-c(upper,up);types<<-c(types,type);k
    }
    sc<-max(abs(normal));add(seq_len(d),normal/sc,-theta[1]/sc,"==")
    qr<-matrix(0L,d,2)
    for(j in seq_len(d)) {
      qr[j,1]<-add(c(j,w_cols[j]),c(1,-1),0)
      qr[j,2]<-add(c(j,w_cols[j]),c(-1,-1),0)
    }
    for(i in active) {
      us<-integer(d)
      for(j in seq_len(d)) {
        y<-X_ref[i,j];u<-variable(max(abs(c(lo[j],hi[j])-y)));us[j]<-u
        if(y<=lo[j])add(c(u,j),c(1,-1),-y,"==")
        else if(y>=hi[j])add(c(u,j),c(1,1),y,"==")
        else {
          b<-variable(1,"B")
          add(c(j,u),c(1,-1),y);add(c(j,u),c(-1,-1),-y)
          add(c(u,j,b),c(1,-1,2*(y-lo[j])),y-2*lo[j])
          add(c(u,j,b),c(1,1,-2*(hi[j]-y)),y)
        }
      }
      add(c(r_col,us),c(1,rep(-L,d)),z_ref[i])
    }
    objective<-numeric(k);objective[r_col]<-1;objective[w_cols]<--L
    ans<-list(A=slam::simple_triplet_matrix(ii,jj,vv,nrow=length(rhs),ncol=k),
      rhs=rhs,directions=directions,objective=objective,types=types,
      qr=qr,w_cols=w_cols,bounds=list(lower=list(ind=seq_len(k),val=lower),
      upper=list(ind=seq_len(k),val=upper)))
    if(cacheable && length(ls(cache,all.names=TRUE))<128L)assign(key,ans,cache)
    ans
  }
  eval_one<-function(x) {
    x<-as.numeric(x);stopifnot(length(x)==d,all(is.finite(x)),all(x>=0 & x<=1))
    stats$queries<-stats$queries+1L
    if(sum(c(1,x)*theta)>=0)
      return(min(1,min(z_ref+L*boundary_l1_two_hop_distances(x,X_ref,theta))))
    values<-z_ref+L*rowSums(abs(sweep(X_ref,2L,x)))
    active<-sort(order(values)[seq_len(min(4L,n))]);use_presolve<-TRUE
    repeat {
      m<-model(active);rn<-m$rhs;rn[m$qr[,1]]<-x;rn[m$qr[,2]]<--x
      bounds<-m$bounds
      bounds$upper$val[m$w_cols]<-pmax(abs(lo-x),abs(hi-x))
      sol<-tryCatch(Rglpk::Rglpk_solve_LP(m$objective,m$A,m$directions,rn,
        bounds=bounds,types=m$types,max=TRUE,
        control=list(presolve=use_presolve,tm_limit=as.integer(time_limit_ms))),
        error=function(err)list(status=99L,numerical_error=conditionMessage(err)))
      stats$solves<-stats$solves+1L;stats$max_active<-max(stats$max_active,length(active))
      stats$total_active<-stats$total_active+length(active)
      if(is.null(sol$status)||sol$status!=0L) {
        if(use_presolve) {use_presolve<-FALSE;stats$numerical_retries<-stats$numerical_retries+1L;next}
        unseen<-setdiff(seq_len(n),active)
        if(length(unseen)) {active<-sort(c(active,head(unseen[order(values[unseen])],4L)));use_presolve<-TRUE;next}
        stats$certified_fallbacks<-stats$certified_fallbacks+1L
        return(mab_boundary_box_bound(x,X_ref,z_ref,theta,L,lo,hi))
      }
      if(length(sol$solution)!=length(m$objective) ||
         any(!is.finite(c(sol$solution,sol$optimum))) || length(sol$optimum)!=1L)
        stop("Boundary solver returned nonfinite or incomplete output")
      h<-sol$solution[seq_len(d)]
      if(abs(sum(c(1,h)*theta))/max(abs(normal))>MAB_NUMERICAL_TOLERANCE ||
         any(h<lo-MAB_NUMERICAL_TOLERANCE | h>hi+MAB_NUMERICAL_TOLERANCE)) {
        stats$certified_fallbacks<-stats$certified_fallbacks+1L
        return(mab_boundary_box_bound(x,X_ref,z_ref,theta,L,lo,hi))
      }
      upper<-clip01(sol$optimum)
      if(upper==0) {stats$early_zero<-stats$early_zero+1L;return(0)}
      all_values<-z_ref+L*rowSums(abs(sweep(X_ref,2L,h)))
      lower<-clip01(min(1,all_values)-L*sum(abs(h-x)));gap<-upper-lower
      if(abs(gap)<=MAB_NUMERICAL_TOLERANCE) {stats$max_gap<-max(stats$max_gap,abs(gap));return(lower)}
      unseen<-setdiff(seq_len(n),active)
      violating<-unseen[all_values[unseen]<sol$solution[d+1L]-MAB_NUMERICAL_TOLERANCE]
      if(!length(violating)) {
        if(use_presolve) {use_presolve<-FALSE;stats$numerical_retries<-stats$numerical_retries+1L;next}
        if(length(unseen)) {active<-sort(c(active,head(unseen[order(all_values[unseen])],4L)));use_presolve<-TRUE;next}
        stats$certified_fallbacks<-stats$certified_fallbacks+1L
        return(mab_boundary_box_bound(x,X_ref,z_ref,theta,L,lo,hi))
      }
      active<-sort(c(active,head(violating[order(all_values[violating])],1L)))
    }
  }
  # Outside the entire H-box slab, |h_j-x_j| is affine for every h in H.
  # Translating AWAY from the slab adds L*distance on the positive side and
  # subtracts it on the negative side. Clipping commutes in that direction.
  # Thus repeated tangential locations reuse one exact boundary-slab query.
  axis<-which.max(abs(normal));slab_cache<-new.env(hash=TRUE,parent=emptyenv())
  stats$slab_cache_hits<-0L;stats$slab_cache_entries<-0L
  eval_batch<-function(Xnew) {
    Xnew<-as.matrix(Xnew);stopifnot(ncol(Xnew)==d)
    vapply(seq_len(nrow(Xnew)),function(i) {
      mab_check_deadline()
      x<-Xnew[i,]
      side<-if(x[axis]<=lo[axis])-1L else if(x[axis]>=hi[axis])1L else 0L
      if(side==0L)return(eval_one(x))
      key<-paste(c(side,sprintf("%.17g",x[-axis])),collapse=",")
      cacheable<-nchar(key,type="bytes")<=MAB_CACHE_KEY_MAX_BYTES
      positive<-side*sign(normal[axis])>0
      base<-if(side<0)lo[axis] else hi[axis]
      if(cacheable && exists(key,slab_cache,inherits=FALSE)) {
        v<-get(key,slab_cache);stats$slab_cache_hits<-stats$slab_cache_hits+1L
      } else {
        y<-x;y[axis]<-base;v<-eval_one(y)
        if(cacheable && stats$slab_cache_entries<250000L) {
          assign(key,v,slab_cache);stats$slab_cache_entries<-stats$slab_cache_entries+1L
        }
      }
      shift<-L*abs(x[axis]-base)
      if(positive)min(1,v+shift) else max(0,v-shift)
    },numeric(1))
  }
  structure(list(type="boundary_canonical_l1_reference",d=d,
    params=list(theta=theta,X_ref=X_ref,z_ref=z_ref,L=L,nu=1),
    meta=utils::modifyList(meta,list(label=label,solver_acceleration="exact_hyperplane_bounds_lazy_cone_models_v2",
      n_cones_before=n_before,n_cones_after=n,canonical_eval="continuous_l1_milp_reference",
      fairness_guaranteed=TRUE,uses_boundary_samples=FALSE,n_boundary_samples=0L,
      probability_certificate_tolerance=MAB_NUMERICAL_TOLERANCE,solver_time_limit_ms=time_limit_ms)),
    solver_stats=stats,eval_one=eval_one,eval_batch=eval_batch),class="estimated_policy")
}

mab_boundary_box_bound <- function(x,X,z,theta,L,lo,hi,
    tolerance=MAB_NUMERICAL_TOLERANCE,max_nodes=20000L) {
  stopifnot(length(tolerance)==1L,is.finite(tolerance),tolerance>0,
    tolerance<=MAB_NUMERICAL_TOLERANCE)
  n<-nrow(X);d<-ncol(X);normal<-theta[-1];scale<-max(abs(normal))
  best<-0;queue<-list(list(lo=lo,hi=hi,upper=1));visited<-0L
  bound<-function(a,b) {
    # Skip boxes that cannot intersect the hyperplane.
    if(theta[1]+sum(pmin(normal*a,normal*b))>1e-12 ||
       theta[1]+sum(pmax(normal*a,normal*b))< -1e-12)return(NULL)
    width<-b-a;A<-C<-matrix(0,n,d)
    for(j in seq_len(d)) {
      if(width[j]<1e-14) {A[,j]<-0;C[,j]<-abs(a[j]-X[,j])} else {
        A[,j]<-(abs(b[j]-X[,j])-abs(a[j]-X[,j]))/width[j]
        C[,j]<-abs(a[j]-X[,j])-A[,j]*a[j]
      }
    }
    mat<-matrix(0,1+2*d+n,2*d+1);rhs<-numeric(nrow(mat));dir<-rep("<=",nrow(mat))
    mat[1,1:d]<-normal/scale;rhs[1]<- -theta[1]/scale;dir[1]<-"=="
    for(j in seq_len(d)) {
      mat[1+2*j-1,c(j,d+1+j)]<-c(1,-1);rhs[1+2*j-1]<-x[j]
      mat[1+2*j,c(j,d+1+j)]<-c(-1,-1);rhs[1+2*j]<- -x[j]
    }
    rr<-1+2*d+seq_len(n);mat[rr,1:d]<- -L*A;mat[rr,d+1]<-1;rhs[rr]<-z+L*rowSums(C)
    objective<-c(rep(0,d),1,rep(-L,d))
    sol<-Rglpk::Rglpk_solve_LP(objective,mat,dir,rhs,
      bounds=list(lower=list(ind=1:(2*d+1),val=c(a,0,rep(0,d))),
        upper=list(ind=1:(2*d+1),val=c(b,1,rep(1,d)))),max=TRUE,
      control=list(presolve=TRUE,canonicalize_status=FALSE,tm_limit=getOption("mab.boundary_time_ms",60000L)))
    if(sol$status==4L)return(NULL) # GLP_NOFEAS, not a time-limit result
    if(sol$status!=5L)stop("Box-bound LP did not certify optimality: ",sol$status)
    h<-sol$solution[1:d]
    if(any(!is.finite(c(h,sol$optimum))))stop("Box-bound LP returned nonfinite output")
    if(abs(sum(c(1,h)*theta))/scale>MAB_NUMERICAL_TOLERANCE ||
       any(h<a-MAB_NUMERICAL_TOLERANCE | h>b+MAB_NUMERICAL_TOLERANCE))
      stop("Box-bound LP witness is infeasible")
    values<-z+L*rowSums(abs(sweep(X,2,h)))
    lower<-clip01(min(1,values)-L*sum(abs(h-x)))
    upper<-clip01(sol$optimum)
    if(lower>upper+tolerance)stop("Inconsistent box-bound certificate")
    list(lo=a,hi=b,upper=upper,lower=lower,h=h,values=values,A=A,C=C)
  }
  while(length(queue)) {
    mab_check_deadline();visited<-visited+1L
    if(visited>max_nodes)stop("Certified boundary fallback exceeded subdivision budget")
    index<-which.max(vapply(queue,`[[`,numeric(1),"upper"));node<-queue[[index]];queue[[index]]<-NULL
    if(node$upper<=best+tolerance)next
    node<-bound(node$lo,node$hi);if(is.null(node))next
    best<-max(best,node$lower)
    if(node$upper<=best+tolerance)next
    cone<-which.min(node$values)
    error<-node$A[cone,]*node$h+node$C[cone,]-abs(node$h-X[cone,])
    eligible<-X[cone,]>node$lo+1e-13 & X[cone,]<node$hi-1e-13
    error[!eligible]<--Inf;j<-which.max(error)
    if(!any(eligible)||!is.finite(error[j])||error[j]<=0)
      stop("Boundary fallback cannot close floating-point certificate gap")
    cut<-X[cone,j];left<-node$hi;left[j]<-cut;right<-node$lo;right[j]<-cut
    queue[[length(queue)+1L]]<-list(lo=node$lo,hi=left,upper=node$upper)
    queue[[length(queue)+1L]]<-list(lo=right,hi=node$hi,upper=node$upper)
  }
  best
}

mab_alpha <- function(T,d,lambda,sigma,eta,C_beta) {
  stopifnot(T>=1,lambda>0,sigma>=0,eta>0,eta<1,C_beta>0)
  sqrt(2*sigma^2*(log(16/eta)+(d+1)/2*log1p(T/lambda)))+C_beta*sqrt(lambda)
}

mab_vertices <- function(V,beta,alpha) {
  a<-eigen(V,symmetric=TRUE);if(any(a$values<=0))stop("Ridge Gram matrix not positive definite")
  inverse_root<-a$vectors%*%diag(1/sqrt(a$values),nrow(V))%*%t(a$vectors)
  deviations<-sqrt(length(beta))*alpha*inverse_root
  # Fixed order: +e1,-e1,+e2,-e2,...; deterministic first maximizing pair.
  out<-matrix(0,2*length(beta),length(beta))
  for(j in seq_along(beta)) {out[2*j-1,]<-beta+deviations[,j];out[2*j,]<-beta-deviations[,j]}
  out
}

mab_vertex_update <- function(X,V,beta,alpha,L,nu=1,lp_neighbor_k=20L) {
  lp_neighbor_k<-mab_validate_neighbor_k(lp_neighbor_k)
  vertices<-lapply(1:2,function(k)mab_vertices(V[[k]],beta[[k]],alpha))
  best<- -Inf;log<-list();winner<-NULL;idx<-0L;xm<-colMeans(add_intercept(X));n<-nrow(X)
  mab_geometry_init()
  edges<-if(ncol(X)>1 && n>1)exact_l1_knn_cpp(X,as.integer(min(lp_neighbor_k,n-1))) else NULL
  for(i in seq_len(nrow(vertices[[1]])))for(j in seq_len(nrow(vertices[[2]]))) {
    mab_check_deadline()
    theta<-vertices[[1]][i,]-vertices[[2]][j,]
    fit<-solve_empirical_lipschitz_lp(X,theta,rep(1/n,n),L,nu,
      lp_neighbor_k=lp_neighbor_k,initial_edges=edges)
    score<-sum(xm*vertices[[2]][j,])+fit$objective_value;idx<-idx+1L
    log[[idx]]<-data.frame(vertex1=i,vertex2=j,score=score,n_history=n,
      lp_neighbor_k=lp_neighbor_k,n_initial_edges=fit$n_initial_edges,
      constraint_generation_rounds=fit$constraint_generation_rounds,
      max_violation=fit$max_pair_violation,all_pairs_certified=fit$all_pairs_certified,
      boundary_intersects=hyperplane_intersects_unit_box(theta))
    if(score>best) {
      best<-score;winner<-list(theta=theta,z=fit$z_ref,beta1=vertices[[1]][i,],
        beta2=vertices[[2]][j,],vertex1=i,vertex2=j,fit=fit)
    }
  }
  meta<-list(n_opt=n,n_lp_points=n,lp_neighbor_k=lp_neighbor_k,
    n_initial_edges=winner$fit$n_initial_edges,
    constraint_generation_rounds=winner$fit$constraint_generation_rounds,
    all_pairs_certified=winner$fit$all_pairs_certified,
    max_pair_violation=winner$fit$max_pair_violation,optimistic_score=best,
    vertex1=winner$vertex1,vertex2=winner$vertex2,vertex_pairs=idx)
  list(policy=mab_policy_from_lp(winner$theta,X,winner$z,L,nu,meta),
    selected_beta1=winner$beta1,selected_beta2=winner$beta2,vertex_log=do.call(rbind,log))
}

mab_pack_policy <- function(p) list(type=p$type,d=p$d,params=p$params,meta=p$meta)

mab_unpack_policy <- function(p) {
  a <- p$params
  ans <- switch(p$type,
    constant = make_constant_policy(a$prob, p$d),
    boundary_canonical = make_d1_boundary_policy(a$theta, a$H_ref, a$q_ref, a$L, a$nu),
    boundary_canonical_continuous_d2_l1 = make_d2_l1_trace_policy(a$theta, a$trace_origin, a$trace_direction, a$trace_knots, a$q_ref, a$L),
    boundary_canonical_l1_reference = make_l1_boundary_tight_policy(a$theta, a$X_ref, a$z_ref, a$L, meta = p$meta),
    standard_ucb = mab_standard_policy(a$beta, a$Vinv, a$alpha),
    stop("Unsupported saved policy: ", p$type))
  ans$meta <- utils::modifyList(ans$meta, p$meta)
  if(p$type=="boundary_canonical_l1_reference")
    ans$meta$probability_certificate_tolerance <- MAB_NUMERICAL_TOLERANCE
  ans
}

mab_standard_policy <- function(beta,Vinv,alpha) {
  evaluate<-function(X)standard_ucb_policy_matrix(X,beta[[1]],beta[[2]],Vinv[[1]],Vinv[[2]],alpha)
  structure(list(type="standard_ucb",d=length(beta[[1]])-1L,params=list(beta=beta,Vinv=Vinv,alpha=alpha),
    meta=list(fairness_guaranteed=FALSE),eval_batch=evaluate,eval_one=function(x)evaluate(matrix(x,nrow=1))[1]),
    class="estimated_policy")
}

mab_check_deadline <- function() {
  if(isTRUE(.mab_geometry_loaded)) {
    rss<-mab_current_rss_mb()
    if(is.finite(rss)&&rss>getOption("mab.max_process_mb",Inf))stop("Configured process RSS budget reached")
  }
  deadline<-getOption("mab.deadline",Inf)
  if(as.numeric(Sys.time())>deadline)stop("Configured wall-clock budget reached")
}

mab_train <- function(X,y1,y2,algorithm=c("original","polyhedral_ucb","standard_ucb"),
    L=.5,nu=1,lambda=1,sigma=.3,eta=.05,C_beta=2,action_uniforms,
    exploration_multiplier=1,progress_file=NULL,lp_neighbor_k=20L,
    checkpoint_file=NULL,checkpoint_every=128L) {
  lp_neighbor_k<-mab_validate_neighbor_k(lp_neighbor_k)
  algorithm<-match.arg(algorithm);X<-as.matrix(X);T<-nrow(X);d<-ncol(X);p<-d+1L
  stopifnot(length(y1)==T,length(y2)==T,length(action_uniforms)==T,
    all(is.finite(c(X,y1,y2,action_uniforms))),all(action_uniforms>=0 & action_uniforms<=1),
    nu==1,lambda>0,L>=0,exploration_multiplier>0)
  To<-if(algorithm=="original")ceiling(sqrt(T*log(4*T))) else if(algorithm=="polyhedral_ucb")ceiling(sqrt(T)) else 0L
  To<-ceiling(To*exploration_multiplier)
  if(2*To>=T)stop("Algorithm requires 2*T_o<T; increase T or explicitly change exploration_multiplier in sensitivity analysis.")
  V<-replicate(2,diag(lambda,p),simplify=FALSE);Vinv<-replicate(2,diag(1/lambda,p),simplify=FALSE)
  b<-replicate(2,numeric(p),simplify=FALSE);ell<-rep(p*log(lambda),2);snap<-ell
  alpha<-mab_alpha(T,d,lambda,sigma,eta,C_beta)
  current<-make_constant_policy(.5,d);history<-list(list(start=1L,policy=mab_pack_policy(current),
    beta1_hat=rep(0,p),beta2_hat=rep(0,p),selected_beta1=rep(0,p),selected_beta2=rep(0,p)))
  vertex_logs<-list();update_logs<-list();s<-2*To+1L
  prob<-reward<-numeric(T);arm<-integer(T);Xt<-add_intercept(X)
  input<-list(X=X,y1=y1,y2=y2,action_uniforms=action_uniforms,algorithm=algorithm,
    L=L,nu=nu,lambda=lambda,sigma=sigma,eta=eta,C_beta=C_beta,
    exploration_multiplier=exploration_multiplier,lp_neighbor_k=lp_neighbor_k,
    numerical_tolerance=MAB_NUMERICAL_TOLERANCE)
  next_t<-1L;policy_ready<-FALSE
  stopifnot(length(checkpoint_every)==1L,is.finite(checkpoint_every),checkpoint_every>=1)
  save_checkpoint<-function() {
    if(is.null(checkpoint_file))return(invisible(NULL))
    state<-list(input=input,next_t=next_t,policy_ready=policy_ready,V=V,Vinv=Vinv,
      b=b,ell=ell,snap=snap,current=mab_pack_policy(current),history=history,
      vertex_logs=vertex_logs,update_logs=update_logs,s=s,prob=prob,reward=reward,arm=arm)
    mab_save_rds(state,checkpoint_file)
    if(!is.null(progress_file))mab_save_rds(list(t=next_t-1L,T=T,
      updates=length(update_logs),algorithm=algorithm),progress_file)
  }
  saved<-if(is.null(checkpoint_file))NULL else mab_read_rds(checkpoint_file)
  if(!is.null(saved)) {
    if(!identical(input,saved$input))stop("Checkpoint inputs differ; use a new output directory.")
    next_t<-saved$next_t;policy_ready<-saved$policy_ready
    V<-saved$V;Vinv<-saved$Vinv;b<-saved$b;ell<-saved$ell;snap<-saved$snap
    current<-mab_unpack_policy(saved$current);history<-saved$history
    vertex_logs<-saved$vertex_logs;update_logs<-saved$update_logs;s<-saved$s
    prob<-saved$prob;reward<-saved$reward;arm<-saved$arm
  } else save_checkpoint()
  for(t in if(next_t<=T)seq.int(next_t,T) else integer()) {
    mab_check_deadline()
    determinant_trigger<-algorithm=="polyhedral_ucb" && any(ell-snap>=log(2))
    update<-t>2*To && (t==2*To+1L || t>=2*s || determinant_trigger || algorithm=="standard_ucb")
    if(update && !policy_ready) {
      if(algorithm!="standard_ucb")save_checkpoint()
      s<-t;snap<-ell
      # Greedy uses the updated inverse; UCB solves the ridge system directly.
      beta<-if(algorithm=="original")lapply(1:2,function(k)drop(Vinv[[k]]%*%b[[k]])) else
        lapply(1:2,function(k)drop(solve(V[[k]],b[[k]])))
      selected<-beta;tm<-proc.time()[[3]]
      if(algorithm=="polyhedral_ucb") {
        optimistic<-mab_vertex_update(X[seq_len(t-1),,drop=FALSE],V,beta,alpha,L,nu,lp_neighbor_k)
        current<-optimistic$policy;selected<-list(optimistic$selected_beta1,optimistic$selected_beta2)
        optimistic$vertex_log$update_time<-t;vertex_logs[[length(vertex_logs)+1L]]<-optimistic$vertex_log
      } else if(algorithm=="original") {
        xx<-X[seq_len(t-1),,drop=FALSE];theta<-beta[[1]]-beta[[2]]
        fit<-solve_empirical_lipschitz_lp(xx,theta,rep(1/(t-1),t-1),L,nu,lp_neighbor_k=lp_neighbor_k)
        current<-mab_policy_from_lp(theta,xx,fit$z_ref,L,nu,list(n_opt=t-1,
          lp_neighbor_k=lp_neighbor_k,n_initial_edges=fit$n_initial_edges,
          constraint_generation_rounds=fit$constraint_generation_rounds,
          all_pairs_certified=fit$all_pairs_certified,max_pair_violation=fit$max_pair_violation))
      } else current<-mab_standard_policy(beta,Vinv,alpha)
      rec<-list(start=t,policy=mab_pack_policy(current),beta1_hat=beta[[1]],beta2_hat=beta[[2]],
        selected_beta1=selected[[1]],selected_beta2=selected[[2]])
      if(t==1)history[[1]]<-rec else history[[length(history)+1L]]<-rec
      update_logs[[length(update_logs)+1L]]<-data.frame(time=t,n_history=t-1,
        lp_neighbor_k=if(algorithm=="standard_ucb")NA_integer_ else lp_neighbor_k,
        determinant_trigger=determinant_trigger,logdet1=ell[1],logdet2=ell[2],
        fit_seconds=proc.time()[[3]]-tm,full_history=TRUE)
      policy_ready<-TRUE
      if(algorithm!="standard_ucb")save_checkpoint()
    }
    prob[t]<-if(t<=2*To).5 else current$eval_one(X[t,])
    if(!is.finite(prob[t]) || prob[t]< -MAB_NUMERICAL_TOLERANCE ||
       prob[t]>1+MAB_NUMERICAL_TOLERANCE)stop("Invalid action probability")
    prob[t]<-clip01(prob[t])
    arm[t]<-if(action_uniforms[t]<prob[t])1L else 2L
    k<-arm[t];reward[t]<-if(k==1)y1[t] else y2[t]
    upd<-rank1_update(Vinv[[k]],b[[k]],ell[k],Xt[t,],reward[t])
    Vinv[[k]]<-upd$Vinv;b[[k]]<-upd$b;ell[k]<-upd$ell
    V[[k]]<-V[[k]]+tcrossprod(Xt[t,])
    if(algorithm=="polyhedral_ucb")ell[k]<-as.numeric(determinant(V[[k]],logarithm=TRUE)$modulus)
    next_t<-t+1L;policy_ready<-FALSE
    if(t%%checkpoint_every==0L || t==T)save_checkpoint()
  }
  list(algorithm=algorithm,T=T,d=d,T_o=To,initial_rounds=2L*To,alpha=alpha,full_history=TRUE,
    lp_neighbor_k=if(algorithm=="standard_ucb")NA_integer_ else lp_neighbor_k,
    history=history,actions=data.frame(t=seq_len(T),prob=prob,arm=arm,reward=reward),
    updates=if(length(update_logs))do.call(rbind,update_logs) else data.frame(),
    vertices=if(length(vertex_logs))do.call(rbind,vertex_logs) else data.frame())
}

# Complete-file replacement keeps checkpoints readable after interruption.
mab_save_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary<-paste0(path,".tmp-",Sys.getpid())
  on.exit(unlink(temporary),add=TRUE)
  saveRDS(x,temporary,compress=FALSE)
  if(.Platform$OS.type=="windows" && file.exists(path)) {
    backup<-paste0(path,".previous")
    if(file.exists(backup))unlink(backup)
    if(!file.rename(path,backup))stop("Cannot preserve previous checkpoint: ",path)
    if(!file.rename(temporary,path)) {
      file.rename(backup,path)
      stop("Cannot save checkpoint: ",path)
    }
    unlink(backup)
  } else if(!file.rename(temporary,path))stop("Cannot save checkpoint: ",path)
  invisible(path)
}

mab_read_rds <- function(path,default=NULL) {
  if(!file.exists(path) && file.exists(paste0(path,".previous")))
    return(readRDS(paste0(path,".previous")))
  if(file.exists(path))readRDS(path) else default
}

mab_release_options <- function(study) {
  args<-commandArgs(trailingOnly=TRUE)
  allowed<-args=="--plot-only" |
    grepl("^--(output|worker-id|worker-count|max-tasks|algorithm|oracle)=",args)
  if(any(!allowed))stop("Unknown argument: ",paste(args[!allowed],collapse=", "))
  value<-function(name,default) {
    flag<-paste0("--",gsub("_","-",name),"=")
    hits<-args[startsWith(args,flag)]
    if(length(hits))substring(tail(hits,1L),nchar(flag)+1L) else
      getOption(paste0("mab.release.",name),default)
  }
  plot_only<-"--plot-only"%in%args || isTRUE(getOption("mab.release.plot_only",FALSE))
  root<-dirname(MAB_SOURCE_FILE)
  output_dir<-value("output",file.path(root,"results",study))
  worker_id<-as.integer(value("worker_id",1L));worker_count<-as.integer(value("worker_count",1L))
  max_tasks<-as.numeric(value("max_tasks",Inf))
  if(anyNA(c(worker_id,worker_count,max_tasks)) || worker_id<1L ||
     worker_count<worker_id || max_tasks<0)stop("Invalid worker or task-count option.")
  list(plot_only=plot_only,output_dir=output_dir,
    worker_id=worker_id,worker_count=worker_count,max_tasks=max_tasks)
}
