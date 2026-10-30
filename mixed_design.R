# Mixed discrete--continuous experiment: six horizons, 50 replications each.
# Rscript mixed_design.R                 runs/resumes the paper experiment
# Rscript mixed_design.R --plot-only     redraws completed results
# Set options(mab.release.autorun=FALSE) before source() to load functions only.

.mx_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  path<-if(length(paths))tail(paths,1)[[1]] else if(length(arg))arg[1] else "mixed_design.R"
  if(!file.exists(path))path<-gsub("~+~"," ",path,fixed=TRUE)
  normalizePath(path,mustWork=TRUE)
})
.mx_root <- dirname(.mx_file)
source(file.path(.mx_root, "algorithms.R"), local=TRUE)

mx_with_seed <- function(seed, fun) {
  old_kind <- RNGkind()
  old_seed <- if(exists(".Random.seed", .GlobalEnv, inherits=FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if(is.null(old_seed)) {
      if(exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
    } else assign(".Random.seed", old_seed, .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fun()
}

mx_settings <- function(algorithm="original") {
  stopifnot(algorithm %in% c("original", "polyhedral_ucb"))
  list(version="mixed_d4_release_1", d=4L, categories=6L, m=4L,
    seed=202610111L, environment_seed=202610111L, rho=.5,
    continuous_reward_coefficient=.03, L=.5, nu=1, sigma=.3,
    lambda=1, eta=.05, C_beta=2, exploration_multiplier=1,
    lp_neighbor_k=20L, horizons=as.integer(2^(10:15)),
    repetitions=50L, algorithm=algorithm,
    qmc_points=32768L, halton_skip=128L,
    oracle_grid=64L, numerical_tolerance=MAB_NUMERICAL_TOLERANCE)
}

mx_environment <- function(cfg) {
  mx_with_seed(cfg$environment_seed, function() sample(c(-1,1), cfg$m^3, replace=TRUE))
}

mx_H <- function(U, signs, m=4L) {
  cell <- pmin(floor(U*m), m-1L)
  local <- U*m-cell
  index <- as.integer(1 + drop(cell %*% m^(0:2)))
  inside <- rowSums(local >= .25 & local < .75) == 3L
  signs[index] * ifelse(inside, ifelse(local[,1] < .5, 1, -1), 0)
}

mx_category_probabilities <- function(U, signs, cfg) {
  (1 + outer(mx_H(U, signs, cfg$m)*cfg$rho, c(-1,-1,-1,1,1,1)))/6
}

mx_stream <- function(T, replication, cfg, signs) {
  seeds <- setNames(lapply(1:5, function(kind)
    as.integer((as.double(cfg$seed) + 4*10000019 + replication*100003 + T*31 + kind*7919) %% 2147483646 + 1)),
    c("environment", "context", "noise1", "noise2", "action"))
  seeds$environment <- cfg$environment_seed
  V <- mx_with_seed(seeds$context, function() matrix(runif(T*4L), ncol=4L))
  U <- V[,2:4,drop=FALSE]
  probability <- mx_category_probabilities(U, signs, cfg)
  cdf <- t(apply(probability,1,cumsum))
  category <- 1L + rowSums(V[,1] > cdf[,1:5,drop=FALSE])
  X <- cbind((category-1)/5, U)
  gap <- X[,1]-.5 + cfg$continuous_reward_coefficient*(rowSums(U)-1.5)
  list(X=X,
    y1=gap + mx_with_seed(seeds$noise1, function() rnorm(T,sd=cfg$sigma)),
    y2=mx_with_seed(seeds$noise2, function() rnorm(T,sd=cfg$sigma)),
    uniforms=mx_with_seed(seeds$action, function() runif(T)), seeds=seeds)
}

mx_halton <- function(n, skip=128L) {
  radical <- function(base) {
    index <- seq.int(skip+1L, length.out=n)
    out <- numeric(n); factor <- 1/base
    while(any(index>0L)) {
      out <- out + factor*(index %% base)
      index <- index %/% base
      factor <- factor/base
    }
    out
  }
  cbind(radical(2L), radical(3L), radical(5L))
}

mx_oracle_environment <- function() {
  list(seed=202610111L, m=4L, rho=.5, L=.5,
    continuous_reward_coefficient=.03,
    signs=c(1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,1,1,1,1,-1,1,1,1,-1))
}

mx_oracle_python <- function() r"---(
"""Build the six-category oracle with trilinear policies on a uniform grid."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import highspy
import numpy as np

TOL = 1e-6


def h_value(u, environment):
    m = int(environment["m"])
    cell = np.minimum((u * m).astype(int), m - 1)
    local = u * m - cell
    inside = ((local >= .25) & (local < .75)).all(axis=1)
    index = cell[:, 0] + m * cell[:, 1] + m * m * cell[:, 2]
    signs = np.asarray(environment["signs"])
    return signs[index] * inside * np.where(local[:, 0] < .5, 1., -1.)


def lattice(q):
    index = np.arange(q ** 3)
    return np.column_stack((index % q, (index // q) % q, index // q ** 2))


def coefficients(n, environment):
    # Split at grid boundaries and at every discontinuity of H. On each piece,
    # two-point tensor Gauss quadrature integrates the objective exactly.
    q = int(np.lcm(n, 4 * int(environment["m"])))
    cells = lattice(q)
    cost = np.zeros((6, (n + 1) ** 3))
    absolute_gap, continuum_error = 0., 0.
    mass = np.zeros(6)
    for gauss in itertools.product((-1., 1.), repeat=3):
        u = (cells + .5 + np.asarray(gauss) / (2 * np.sqrt(3))) / q
        h = h_value(u, environment)
        v = u * n
        lo = np.minimum(v.astype(int), n - 1)
        t = v - lo
        distance = (2 * t * (1 - t)).sum(axis=1) / n
        for g in range(6):
            weight = (1 + environment["rho"] * (-1 if g < 3 else 1) * h) / (48 * q ** 3)
            gap = g / 5 - .5 + environment["continuous_reward_coefficient"] * (u.sum(axis=1) - 1.5)
            obj = weight * gap
            mass[g] += weight.sum()
            absolute_gap += np.sum(weight * np.abs(gap))
            continuum_error += environment["L"] * np.sum(weight * np.abs(gap) * distance)
            for corner in itertools.product((0, 1), repeat=3):
                c = np.asarray(corner)
                index = (lo[:, 0] + c[0]) + (n + 1) * (lo[:, 1] + c[1]) + (n + 1) ** 2 * (lo[:, 2] + c[2])
                basis = np.prod(np.where(c, t, 1 - t), axis=1)
                np.add.at(cost[g], index, obj * basis)
    if np.max(np.abs(mass - 1 / 6)) > TOL:
        raise RuntimeError("Invalid category masses in oracle objective")
    return cost.ravel(), continuum_error, absolute_gap, mass


def edges(n):
    # Coordinate 1 varies fastest inside each category block.
    ids = np.arange(6 * (n + 1) ** 3, dtype=np.int32).reshape((6, n + 1, n + 1, n + 1))
    ii, jj, distances = [], [], []
    for axis in range(4):
        a, b = [slice(None)] * 4, [slice(None)] * 4
        a[axis], b[axis] = slice(None, -1), slice(1, None)
        i, j = ids[tuple(a)].ravel(), ids[tuple(b)].ravel()
        ii.append(i)
        jj.append(j)
        distances.append(np.full(len(i), 1 / 5 if axis == 0 else 1 / n))
    return np.concatenate(ii), np.concatenate(jj), np.concatenate(distances)


def solve(n, environment, output, threads=1):
    output.mkdir(parents=True, exist_ok=True)
    objective, error, absolute_gap, mass = coefficients(n, environment)
    i, j, distance = edges(n)
    scale = 1 / distance
    factor = max(np.max(np.abs(objective)), 1e-15)
    cost = -objective / factor
    L = environment["L"]
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = len(cost), len(i)
    lp.col_cost_, lp.col_lower_, lp.col_upper_ = cost, np.zeros(len(cost)), np.ones(len(cost))
    lp.row_lower_, lp.row_upper_ = np.full(len(i), -L), np.full(len(i), L)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.arange(len(i) + 1, dtype=np.int64) * 2
    lp.a_matrix_.index_ = np.column_stack((i, j)).ravel().astype(np.int32)
    lp.a_matrix_.value_ = np.column_stack((scale, -scale)).ravel()
    solver = highspy.Highs()
    options = dict(threads=threads, parallel="on", solver="ipm", presolve="on",
                   run_crossover="off", primal_feasibility_tolerance=TOL,
                   dual_feasibility_tolerance=TOL, ipm_optimality_tolerance=TOL,
                   random_seed=0, output_flag=True, log_file=str(output / "solver.log"))
    for name, value in options.items():
        if solver.setOptionValue(name, value) != highspy.HighsStatus.kOk:
            raise RuntimeError(f"Could not set HiGHS option {name}")
    if solver.passModel(lp) != highspy.HighsStatus.kOk:
        raise RuntimeError("HiGHS rejected the oracle model")
    started = time.monotonic()
    solver.run()
    if solver.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        raise RuntimeError(solver.modelStatusToString(solver.getModelStatus()))
    solution = solver.getSolution()
    raw = np.asarray(solution.col_value)
    if not np.isfinite(raw).all() or raw.min() < -TOL or raw.max() > 1 + TOL:
        raise RuntimeError("Oracle probability bounds failed")
    z = np.clip(raw, 0, 1)
    slope = float(np.max(np.abs(z[i] - z[j]) / distance))
    if slope > L + TOL:
        raise RuntimeError("Oracle Lipschitz bound failed")
    # Remove accepted floating-point edge excess without changing L.
    contraction = min(1., L / (slope + 1e-12))
    z = .5 + contraction * (z - .5)
    y = np.asarray(solution.row_dual)
    at_y = np.zeros(len(z))
    np.add.at(at_y, i, scale * y)
    np.add.at(at_y, j, -scale * y)
    reduced = cost - at_y
    dual_min = (np.sum(-L * np.abs(y), dtype=np.longdouble)
                + np.sum(np.minimum(reduced, 0), dtype=np.longdouble))
    value = float(np.sum(objective * z, dtype=np.longdouble))
    upper = float(-dual_min * factor)
    if not -TOL <= upper - value <= TOL:
        raise RuntimeError("Oracle finite-grid objective gap exceeds tolerance")
    path = output / f"mixed_d4_grid{n}.bin"
    temp = path.with_suffix(".tmp")
    z.astype("<f8").tofile(temp)
    temp.replace(path)
    record = dict(grid=n, categories=6, continuous_dimension=3,
                  layout="little-endian float64; U1 varies fastest, then U2, U3, category",
                  environment=environment,
                  finite_objective=value, finite_dual_upper=upper, finite_gap=upper-value,
                  raw_max_slope=slope, contraction=contraction, category_masses=mass.tolist(),
                  expected_absolute_gap=absolute_gap, continuum_discretization_bound=error,
                  continuum_optimality_certified=False, numerical_tolerance=TOL,
                  seconds=time.monotonic()-started, solver_version=solver.version(),
                  probabilities_md5=hashlib.md5(path.read_bytes()).hexdigest(),
                  probabilities_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    path.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    args = parser.parse_args()
    if args.grid < 1 or args.threads < 1:
        parser.error("grid and threads must be positive")
    environment = json.loads(args.environment.read_text())
    solve(args.grid, environment, args.output, args.threads)
)---"

mx_build_oracle <- function(grid=64L, threads=1L,
    output_dir=file.path(.mx_root,"results","oracle"),
    python=Sys.getenv("MAB_PYTHON",unset="python3")) {
  if(length(grid)!=1L || !is.finite(grid) || grid<1 || grid!=as.integer(grid) ||
     length(threads)!=1L || !is.finite(threads) || threads<1 || threads!=as.integer(threads))
    stop("grid and threads must be positive integers")
  for(package in c("jsonlite","filelock"))
    if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," to construct the numerical oracle")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  path <- file.path(output_dir,sprintf("mixed_d4_grid%d.bin",as.integer(grid)))
  lock <- filelock::lock(paste0(path,".lock"),timeout=0)
  if(is.null(lock)) {
    message("Waiting for the numerical oracle being calculated by another worker...")
    lock <- filelock::lock(paste0(path,".lock"),timeout=Inf)
  }
  if(is.null(lock)) stop("Could not lock the numerical oracle cache")
  on.exit(filelock::unlock(lock),add=TRUE)
  if(file.exists(path)) {
    cached <- tryCatch(mx_oracle(path=path,expected_grid=as.integer(grid)),error=function(e) NULL)
    if(!is.null(cached)) return(invisible(path))
    message("The cached oracle is incomplete or incompatible; recalculating it.")
  }
  executable <- if(file.exists(python)) normalizePath(python) else unname(Sys.which(python))
  if(!nzchar(executable)) stop("Python was not found; set MAB_PYTHON to its executable path")
  stage <- tempfile(sprintf(".grid%d_build_",as.integer(grid)),tmpdir=output_dir)
  dir.create(stage,recursive=TRUE,showWarnings=FALSE)
  script <- tempfile("mixed_oracle_",fileext=".py")
  environment <- tempfile("mixed_environment_",fileext=".json")
  on.exit(unlink(c(script,environment)),add=TRUE)
  writeLines(mx_oracle_python(),script,useBytes=TRUE)
  jsonlite::write_json(mx_oracle_environment(),environment,auto_unbox=TRUE,pretty=TRUE,digits=NA)
  message(sprintf("Calculating the mixed numerical oracle on the %d x %d x %d grid; progress follows below.",
    as.integer(grid)+1L,as.integer(grid)+1L,as.integer(grid)+1L))
  status <- system2(executable,c(shQuote(script),"--grid",as.integer(grid),
    "--threads",as.integer(threads),"--output",shQuote(stage),
    "--environment",shQuote(environment)))
  if(status!=0L) stop("Numerical oracle construction failed; Python requires numpy and highspy. Solver output: ",stage)
  built <- file.path(stage,basename(path))
  mx_oracle(path=built,expected_grid=as.integer(grid))
  for(suffix in c(".bin",".json")) {
    target <- sub("\\.bin$",suffix,path)
    source <- sub("\\.bin$",suffix,built)
    if(file.exists(target) && unlink(target)!=0L) stop("Could not replace incomplete oracle file: ",target)
    if(!file.rename(source,target)) stop("Could not save numerical oracle: ",target)
  }
  if(file.exists(file.path(stage,"solver.log")))
    file.copy(file.path(stage,"solver.log"),file.path(output_dir,sprintf("grid%d_solver.log",as.integer(grid))),overwrite=TRUE)
  unlink(stage,recursive=TRUE)
  message("Saved numerical oracle: ",path)
  invisible(path)
}

mx_oracle <- function(root=.mx_root, path=NULL, expected_grid=64L) {
  if(is.null(path)) path <- mx_build_oracle(grid=expected_grid,output_dir=file.path(root,"results","oracle"))
  if(!file.exists(path)) stop("Numerical oracle file does not exist: ",path)
  if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite to read the numerical oracle")
  meta_path <- sub("\\.bin$",".json",path)
  if(identical(meta_path,path) || !file.exists(meta_path)) stop("The numerical oracle requires its matching .json metadata")
  meta <- jsonlite::read_json(meta_path,simplifyVector=TRUE)
  actual_md5 <- unname(tools::md5sum(path))
  scalar <- function(x) is.numeric(x) && length(x)==1L && is.finite(x)
  if(!identical(meta$probabilities_md5,actual_md5) ||
     !isTRUE(all.equal(meta$environment,mx_oracle_environment(),tolerance=0,check.attributes=FALSE)) ||
     !isTRUE(meta$grid==expected_grid) || !isTRUE(meta$categories==6) ||
     !isTRUE(meta$continuous_dimension==3) ||
     !scalar(meta$finite_gap) || abs(meta$finite_gap)>MAB_NUMERICAL_TOLERANCE ||
     !scalar(meta$finite_objective) || !scalar(meta$finite_dual_upper) ||
     abs(meta$finite_dual_upper-meta$finite_objective-meta$finite_gap)>MAB_NUMERICAL_TOLERANCE ||
     !isTRUE(meta$numerical_tolerance==MAB_NUMERICAL_TOLERANCE))
    stop("Numerical oracle metadata does not match this mixed experiment")
  grid <- as.integer(meta$grid)
  expected_length <- 6L*(grid+1L)^3
  probabilities <- readBin(path,"double",n=expected_length,endian="little")
  if(file.info(path)$size!=8*expected_length || length(probabilities)!=expected_length ||
     !all(is.finite(probabilities)) ||
     !all(probabilities >= -MAB_NUMERICAL_TOLERANCE & probabilities <= 1+MAB_NUMERICAL_TOLERANCE))
    stop("Numerical oracle probabilities are invalid")
  probabilities <- clip01(probabilities)
  a <- array(probabilities,dim=c(rep(grid+1L,3L),6L))
  slope <- max(max(abs(a[-1,,,]-a[-(grid+1L),,,]))*grid,
    max(abs(a[,-1,,]-a[,-(grid+1L),,]))*grid,
    max(abs(a[,,-1,]-a[,,-(grid+1L),]))*grid,
    max(abs(a[,,,-1]-a[,,,-6]))*5)
  if(slope>mx_oracle_environment()$L+MAB_NUMERICAL_TOLERANCE)
    stop("Numerical oracle violates the fairness bound")
  list(grid=grid,probabilities=probabilities,md5=actual_md5)
}

mx_oracle_predict <- function(oracle, U, category) {
  n <- oracle$grid; m <- n+1L
  v <- U*n; lo <- pmin(floor(v),n-1L); t <- v-lo
  out <- numeric(nrow(U))
  for(a in 0:1) for(b in 0:1) for(c in 0:1) {
    index <- 1+lo[,1]+a + m*(lo[,2]+b) + m^2*(lo[,3]+c) + (category-1)*m^3
    w <- (if(a) t[,1] else 1-t[,1]) * (if(b) t[,2] else 1-t[,2]) * (if(c) t[,3] else 1-t[,3])
    out <- out + oracle$probabilities[index]*w
  }
  out
}

mx_integration <- function(cfg, signs, oracle) {
  U <- mx_halton(cfg$qmc_points, cfg$halton_skip)
  weights <- mx_category_probabilities(U, signs, cfg)/nrow(U)
  gap <- sapply(1:6, function(g) (g-1)/5-.5 + cfg$continuous_reward_coefficient*(rowSums(U)-1.5))
  oracle_probability <- sapply(1:6, function(g) mx_oracle_predict(oracle,U,g))
  list(U=U, weights=weights, gap=gap, oracle_probability=oracle_probability,
    oracle_reward=sum(vapply(1:6, function(g) sum(weights[,g]*gap[,g]*oracle_probability[,g]), numeric(1))))
}

mx_evaluate <- function(training, integration, cfg, folder) {
  # Every completed block is saved. An interrupted block is recalculated on
  # restart using the same policy, points and weights; completed blocks survive.
  starts <- vapply(training$history, `[[`, integer(1), "start")
  duration <- diff(c(starts, training$T+1L))
  n <- nrow(integration$U)
  block_size <- 512L
  block_starts <- seq.int(1L,n,by=block_size)
  values <- rates <- numeric(length(starts))
  progress_path <- file.path(folder,"evaluation_progress.rds")
  dir.create(file.path(folder,"evaluation"), recursive=TRUE, showWarnings=FALSE)
  completed <- 0L; total <- length(starts)*6L*length(block_starts)
  for(j in seq_along(starts)) {
    message(sprintf("  evaluating policy %d/%d on %d shared Halton points",j,length(starts),n))
    policy <- mab_unpack_policy(training$history[[j]]$policy)
    for(g in 1:6) {
      p <- numeric(n)
      for(first in block_starts) {
        index <- first:min(n,first+block_size-1L)
        cache <- file.path(folder,"evaluation",sprintf("policy_%03d_category_%d_block_%06d.rds",j,g,first))
        block <- mab_read_rds(cache)
        if(is.null(block)) {
          X <- cbind(rep((g-1)/5,length(index)), integration$U[index,,drop=FALSE])
          block <- policy$eval_batch(X)
          if(length(block)==length(index) && all(is.finite(block)) &&
             all(block >= -MAB_NUMERICAL_TOLERANCE & block <= 1+MAB_NUMERICAL_TOLERANCE)) {
            block <- clip01(block)
            mab_save_rds(block,cache)
          }
        }
        if(length(block)!=length(index) || any(!is.finite(block)) ||
           any(block < -MAB_NUMERICAL_TOLERANCE | block > 1+MAB_NUMERICAL_TOLERANCE))
          stop("Invalid mixed-design evaluation probabilities")
        p[index] <- block
        completed <- completed+1L
        mab_save_rds(list(stage="evaluation",blocks=completed,total_blocks=total,
          policy=j,category=g,updated=as.character(Sys.time())),progress_path)
      }
      wgap <- integration$weights[,g]*integration$gap[,g]
      values[j] <- values[j]+sum(wgap*p)
      rates[j] <- rates[j]+sum(wgap*(integration$oracle_probability[,g]-p))
    }
  }
  initial_duration <- pmax(0,pmin(starts+duration-1L,training$initial_rounds)-starts+1L)
  total_regret <- sum(duration*rates)
  list(total_regret=total_regret, initial_regret=sum(initial_duration*rates),
    post_regret=sum((duration-initial_duration)*rates),
    relative_regret=total_regret/(training$T*integration$oracle_reward),
    oracle_value=integration$oracle_reward,
    decomposition=data.frame(start=starts,duration=duration,policy_reward=values,regret_rate=rates))
}

mx_tasks <- function(cfg) {
  tasks <- expand.grid(T=cfg$horizons,rep=seq_len(cfg$repetitions))
  tasks$stage <- findInterval(tasks$rep-1L,c(0L,5L,10L,20L))
  tasks <- tasks[order(tasks$stage,tasks$rep,match(tasks$T,cfg$horizons)),]
  tasks$index <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  tasks
}

mx_summary <- function(runs) {
  do.call(rbind,lapply(split(runs,runs$T),function(x) data.frame(T=x$T[1],n=nrow(x),
    mean_regret=mean(x$total_regret),sd_regret=if(nrow(x)>1) sd(x$total_regret) else NA_real_,
    mean_relative=mean(x$relative_regret),sd_relative=if(nrow(x)>1) sd(x$relative_regret) else NA_real_)))
}

mx_plot_theme <- function() {
  ggplot2::theme_bw(base_size=12,base_family="Helvetica") + ggplot2::theme(
    text=ggplot2::element_text(colour="black"),
    axis.title=ggplot2::element_text(size=13,colour="black"),
    axis.text=ggplot2::element_text(size=11,colour="black"),
    axis.ticks=ggplot2::element_line(colour="black",linewidth=.45),axis.ticks.length=grid::unit(3,"pt"),
    panel.grid.major=ggplot2::element_line(colour="#DCDCDC",linewidth=.35),
    panel.grid.minor=ggplot2::element_blank(),
    panel.border=ggplot2::element_rect(colour="black",fill=NA,linewidth=.5),
    panel.spacing=grid::unit(12.7,"pt"),
    strip.background=ggplot2::element_rect(fill="#EBEBEB",colour="black",linewidth=.5),
    strip.text=ggplot2::element_text(size=12,colour="black",face="plain",margin=ggplot2::margin(5.5,5.5,5.5,5.5)),
    axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=3)),
    axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=6)),
    legend.position="bottom",legend.direction="horizontal",legend.text=ggplot2::element_text(size=11),
    legend.title=ggplot2::element_text(size=11),legend.key.width=grid::unit(16,"pt"),
    legend.key.height=grid::unit(14,"pt"),legend.spacing.x=grid::unit(5,"pt"),
    legend.box.spacing=grid::unit(22,"pt"),legend.margin=ggplot2::margin(0,0,0,0),
    plot.margin=ggplot2::margin(57,56,10.5,38,"pt"))
}

mx_plot <- function(runs, output_dir) {
  if(!nrow(runs)) return(invisible(NULL))
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 to create figures")
  s <- mx_summary(runs); s <- s[order(s$T),]
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(runs,file.path(output_dir,"mixed_runs.csv"),row.names=FALSE)
  write.csv(s,file.path(output_dir,"mixed_summary.csv"),row.names=FALSE)
  if(nrow(s)>1 && all(s$mean_regret>0)) {
    fit <- lm(log(mean_regret)~log(T),data=s)
    write.csv(data.frame(slope=unname(coef(fit)[2]),intercept=unname(coef(fit)[1]),
      n_horizons=nrow(s)),file.path(output_dir,"mixed_loglog_fit.csv"),row.names=FALSE)
  }
  for(relative in c(FALSE,TRUE)) {
    d <- data.frame(T=s$T,value=if(relative) s$mean_relative else s$mean_regret,
      sd=if(relative) s$sd_relative else s$sd_regret,dimension=factor("4"))
    d$lower <- d$value-d$sd; d$upper <- d$value+d$sd
    # A nonpositive lower SD endpoint cannot be displayed on a logarithmic axis.
    if(!relative) d$lower[d$lower<=0] <- NA_real_
    p <- ggplot2::ggplot(d,ggplot2::aes(x=T,y=value,colour=dimension,fill=dimension)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=.15,colour=NA,na.rm=TRUE) +
      ggplot2::geom_line(data=if(nrow(d)>1) d else d[FALSE,],linewidth=.75) +
      ggplot2::geom_point(shape=18,size=1.7,stroke=.5) +
      ggplot2::scale_colour_manual(values=c("4"="#CC79A7"),name="d") +
      ggplot2::scale_fill_manual(values=c("4"="#CC79A7"),guide="none") +
      mx_plot_theme()
    if(relative) {
      p <- p + ggplot2::scale_x_continuous(breaks=if(max(d$T)>=1024L) 8192*(0:4) else sort(unique(d$T)),
        labels=scales::label_comma(accuracy=1),limits=c(0,max(d$T)),expand=ggplot2::expansion(mult=c(.02,.04))) +
        ggplot2::scale_y_continuous(labels=scales::label_number(accuracy=.01),expand=ggplot2::expansion(mult=c(.05,.07))) +
        ggplot2::labs(x="T",y="Relative regret")
    } else {
      loglabels <- function(x) sub("\\.?0+$","",sprintf("%.2f",log10(x)))
      p <- p + ggplot2::scale_x_log10(breaks=sort(unique(d$T)),labels=loglabels,expand=ggplot2::expansion(mult=c(.06,.08))) +
        ggplot2::scale_y_log10(breaks=c(1,3,10,30,100,300,1000),labels=loglabels,expand=ggplot2::expansion(mult=c(.05,.07))) +
        ggplot2::labs(x=expression(log[10](T)),y=expression(log[10](plain("Cumulative regret"))))
    }
    stem <- if(relative) "mixed_d4_relative_regret" else "mixed_d4_cumulative_loglog"
    ggplot2::ggsave(file.path(output_dir,paste0(stem,".pdf")),p,width=8,height=5.2,device="pdf",bg="white")
    ggplot2::ggsave(file.path(output_dir,paste0(stem,".png")),p,width=8,height=5.2,dpi=180,bg="white")
  }
  invisible(s)
}

mx_collect <- function(output_dir) {
  files <- list.files(file.path(output_dir,"tasks"),pattern="^result\\.rds$",recursive=TRUE,full.names=TRUE)
  if(!length(files)) return(data.frame())
  rows <- lapply(files,function(path) {
    r <- readRDS(path)
    data.frame(d=4L,T=r$T,rep=r$rep,total_regret=r$evaluation$total_regret,
      initial_regret=r$evaluation$initial_regret,post_regret=r$evaluation$post_regret,
      oracle_value=r$evaluation$oracle_value,relative_regret=r$evaluation$relative_regret)
  })
  out <- do.call(rbind,rows); out[order(out$T,out$rep),]
}

mixed_design_main <- function() {
  opts <- mab_release_options("mixed")
  argument <- grep("^--algorithm=",commandArgs(trailingOnly=TRUE),value=TRUE)
  algorithm <- if(length(argument)) sub("^--algorithm=","",tail(argument,1)) else getOption("mab.release.algorithm","original")
  cfg <- mx_settings(algorithm)
  output_dir <- opts$output_dir
  if(algorithm!="original") output_dir <- file.path(output_dir,algorithm)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if(opts$plot_only) {
    if(!requireNamespace("filelock",quietly=TRUE)) stop("Install filelock to update figures")
    plot_lock <- filelock::lock(file.path(output_dir,"plot.lock"),timeout=10000)
    if(is.null(plot_lock)) stop("Another process is updating these figures; rerun after it finishes")
    on.exit(filelock::unlock(plot_lock),add=TRUE)
    runs <- mx_collect(output_dir)
    if(!nrow(runs)) {
      saved_runs <- file.path(output_dir,"mixed_runs.csv")
      if(!file.exists(saved_runs)) stop("No completed mixed results are available in ",output_dir,
        ". Run mixed_design.R first; figures are updated after each completed replication.")
      runs <- read.csv(saved_runs)
      required <- c("T","total_regret","relative_regret")
      if(!nrow(runs) || !all(required %in% names(runs))) stop("No valid completed mixed results to plot")
    }
    return(invisible(mx_plot(runs,output_dir)))
  }
  oracle_argument <- grep("^--oracle=",commandArgs(trailingOnly=TRUE),value=TRUE)
  oracle_path <- if(length(oracle_argument)) sub("^--oracle=","",tail(oracle_argument,1)) else getOption("mab.release.oracle",NULL)
  oracle <- mx_oracle(path=oracle_path,expected_grid=cfg$oracle_grid)
  cfg$oracle_md5 <- oracle$md5
  settings_path <- file.path(output_dir,"settings.rds")
  if(!requireNamespace("filelock",quietly=TRUE)) stop("Install filelock to run resumable experiments")
  settings_lock <- filelock::lock(file.path(output_dir,"settings.lock"),timeout=10000)
  if(is.null(settings_lock)) stop("Another process is configuring this output directory; retry shortly")
  tryCatch({
    prior <- mab_read_rds(settings_path)
    extra <- if(is.null(prior)) list() else prior[setdiff(names(prior),names(cfg))]
    compatible <- is.null(prior) || (identical(prior[names(cfg)],cfg) &&
      length(extra)<=1L && all(vapply(extra,identical,logical(1),FALSE)))
    if(!compatible) stop("Settings changed: choose a separate --output directory")
    if(is.null(prior)) mab_save_rds(cfg,settings_path)
  },finally=filelock::unlock(settings_lock))
  signs <- mx_environment(cfg)
  integration <- mx_integration(cfg,signs,oracle)
  tasks <- mx_tasks(cfg)
  tasks <- tasks[(tasks$index-1L) %% opts$worker_count == opts$worker_id-1L,,drop=FALSE]
  completed <- 0L
  for(i in seq_len(nrow(tasks))) {
    task <- tasks[i,]
    folder <- file.path(output_dir,"tasks",sprintf("T%05d_rep%02d",task$T,task$rep))
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    result_path <- file.path(folder,"result.rds")
    if(file.exists(result_path)) next
    if(completed >= opts$max_tasks) break
    # File locks prevent two terminals from calculating the same task.
    if(!requireNamespace("filelock",quietly=TRUE)) stop("Install filelock to run resumable experiments")
    lock <- filelock::lock(file.path(folder,"task.lock"),timeout=0)
    if(is.null(lock)) next
    tryCatch({
      if(!file.exists(result_path)) {
        message(sprintf("Mixed d=4 %s: T=%d, replication=%d/%d",algorithm,task$T,task$rep,cfg$repetitions))
        stream <- mx_stream(task$T,task$rep,cfg,signs)
        training_path <- file.path(folder,"training.rds")
        training <- mab_read_rds(training_path)
        if(is.null(training)) {
          training <- mab_train(stream$X,stream$y1,stream$y2,algorithm,
            L=cfg$L,nu=cfg$nu,lambda=cfg$lambda,sigma=cfg$sigma,eta=cfg$eta,C_beta=cfg$C_beta,
            action_uniforms=stream$uniforms,exploration_multiplier=cfg$exploration_multiplier,
            lp_neighbor_k=cfg$lp_neighbor_k,progress_file=file.path(folder,"training_progress.rds"),
            checkpoint_file=file.path(folder,"checkpoint.rds"),checkpoint_every=128L)
          mab_save_rds(training,training_path)
        }
        evaluation <- mx_evaluate(training,integration,cfg,folder)
        mab_save_rds(list(T=task$T,rep=task$rep,settings=cfg,seeds=stream$seeds,
          evaluation=evaluation),result_path)
        completed <- completed+1L
        message(sprintf("Saved: regret %.6f; relative regret %.6f",evaluation$total_regret,evaluation$relative_regret))
      }
    },finally=filelock::unlock(lock))
    plot_lock <- filelock::lock(file.path(output_dir,"plot.lock"),timeout=0)
    if(!is.null(plot_lock)) tryCatch(mx_plot(mx_collect(output_dir),output_dir),finally=filelock::unlock(plot_lock))
  }
  invisible(mx_collect(output_dir))
}

if(!isFALSE(getOption("mab.release.autorun",TRUE))) mixed_design_main()
