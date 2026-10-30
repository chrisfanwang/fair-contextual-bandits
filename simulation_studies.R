# Ordinary simulations, fairness-parameter comparison, and policy illustration.
# Default: 50 replications at each of the six paper horizons.
# Run using source(file.choose()) or Rscript simulation_studies.R.
.sim_entry_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  path <- if(length(paths))tail(paths,1L)[[1L]] else args[[1L]]
  if(!file.exists(path))path<-gsub("~+~"," ",path,fixed=TRUE)
  normalizePath(path,mustWork=TRUE)
})
.sim_release_dir <- dirname(.sim_entry_file)
local({
  old <- getOption("mab.release.autorun")
  on.exit(options(mab.release.autorun = old))
  options(mab.release.autorun = FALSE)
  source(file.path(.sim_release_dir, "cost_of_fairness.R"), local = FALSE)
})

# Independent streams for contexts, the two potential-reward noises, and actions.
# The seed formula is shared across L and algorithms. T enters the formula, so
# each horizon is a separate experiment, rather than a prefix of a longer run.
mab_simulation_seeds <- function(d, T, replication, seed = 20260924L) {
  setNames(lapply(1:5, function(k) as.integer((as.double(seed) + d*10000019 +
    replication*100003 + T*31 + k*7919) %% 2147483646 + 1)),
    c("environment", "context", "noise1", "noise2", "action"))
}

mab_simulation_with_seed <- function(seed, fun) {
  old_kind <- RNGkind()
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old)) {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fun()
}

mab_simulation_stream <- function(task, cfg = mab_simulation_settings()) {
  seeds <- mab_simulation_seeds(task$d, task$T, task$rep, cfg$seed)
  X <- mab_simulation_with_seed(seeds$context,
    function() generate_contexts(task$T, task$d, task$covariate, cfg$ordinary$sigma_X))
  beta <- get_beta_setting(task$regime, task$d)
  mu1 <- drop(add_intercept(X) %*% beta$beta1)
  mu2 <- drop(add_intercept(X) %*% beta$beta2)
  noise1 <- mab_simulation_with_seed(seeds$noise1, function() rnorm(task$T, sd = cfg$sigma))
  noise2 <- mab_simulation_with_seed(seeds$noise2, function() rnorm(task$T, sd = cfg$sigma))
  uniforms <- mab_simulation_with_seed(seeds$action, function() runif(task$T))
  list(X = X, beta = beta, mu1 = mu1, mu2 = mu2, y1 = mu1 + noise1,
    y2 = mu2 + noise2, noise1 = noise1, noise2 = noise2, uniforms = uniforms, seeds = seeds)
}

mab_simulation_tasks <- function() {
  T <- as.integer(2^(10:15))
  reps <- 1:50
  baseline <- expand.grid(d = 1:4,
    covariate = c("uniform", "truncated_gaussian"),
    regime = c("homogeneous", "heterogeneous"), L = .5,
    T = T, rep = reps, algorithm = "original", stringsAsFactors = FALSE)
  vary <- expand.grid(d = 2L, covariate = "uniform", regime = "homogeneous",
    L = c(.25, 1, 2, 4), T = T, rep = reps, algorithm = "original", stringsAsFactors = FALSE)
  illustration <- expand.grid(d = 2L, covariate = "uniform",
    regime = "heterogeneous_d2_arm1_small", L = .5,
    T = c(8192L, 16384L), rep = 15L,
    algorithm = c("original", "standard_ucb"), stringsAsFactors = FALSE)
  tasks <- rbind(baseline, vary, illustration)
  tasks <- tasks[order(tasks$rep, tasks$T, tasks$d, tasks$covariate,
    tasks$regime, tasks$L, tasks$algorithm), ]
  rownames(tasks) <- NULL
  tasks$id <- vapply(seq_len(nrow(tasks)), function(i) mab_simulation_id(as.list(tasks[i, ])), character(1))
  tasks
}

mab_simulation_id <- function(task) paste(task$covariate, task$regime, paste0("d", task$d),
  paste0("L", format(task$L, trim = TRUE, scientific = FALSE)), paste0("T", task$T),
  sprintf("rep%02d", task$rep), task$algorithm, sep = "_")

# Oracle and learner are evaluated by the same reward operator. In d = 1,2
# affine policy regions are integrated analytically. In d = 3,4 both use the
# fixed 16,384 / 32,768 Halton points, skipping 128, with the same weights.
mab_simulation_evaluate <- function(training, task, cfg, checkpoint_file) {
  beta <- get_beta_setting(task$regime, task$d)
  oracle <- mab_ordinary_oracle(task$d, task$regime, task$L, task$covariate)
  saved <- mab_read_rds(checkpoint_file)
  if (is.null(saved)) saved <- list(task = task, values = list(), pending = NULL)
  if (!identical(saved$task, task)) stop("Evaluation checkpoint belongs to another experiment")
  starts <- vapply(training$history, `[[`, integer(1), "start")
  duration <- diff(c(starts, task$T + 1L))
  if (task$d == 1L) {
    density <- mab_d1_density(task$covariate, cfg$ordinary$sigma_X)
    oracle_value <- mab_d1_policy_value(oracle, beta, density)
    evaluate <- function(policy) c(value = mab_d1_policy_value(policy, beta, density),
      regret = mab_d1_regret_rate(oracle, policy, beta, density))
    method <- "analytic_d1_moments"
  } else if (task$d == 2L) {
    density <- mab_d2_density(task$covariate, cfg$ordinary$sigma_X)
    oracle_value <- mab_d2_policy_value(oracle, beta, density)
    oracle_gap_value <- mab_d2_gap_value(mab_d2_policy_regions(oracle), beta$beta1 - beta$beta2, density)
    evaluate <- function(policy) {
      regions <- mab_d2_policy_regions(policy)
      if (is.null(regions)) stop("Nonlinear policy is not used in the fair-regret experiment")
      c(value = mab_d2_policy_value(policy, beta, density, regions),
        regret = oracle_gap_value - mab_d2_gap_value(regions, beta$beta1 - beta$beta2, density))
    }
    method <- "analytic_d2_polygon_moments"
  } else {
    n <- cfg$evaluation[[paste0("qmc_n_d", task$d)]]
    X <- unit_to_context(halton_sequence(n, task$d, cfg$evaluation$halton_skip),
      task$covariate, cfg$ordinary$sigma_X)
    mu2 <- drop(add_intercept(X) %*% beta$beta2)
    gap <- drop(add_intercept(X) %*% (beta$beta1 - beta$beta2))
    op <- oracle$eval_batch(X)
    weights <- rep(1/n, n)
    oracle_value <- sum(weights * (mu2 + gap*op))
    method <- paste0("halton_d", task$d)
  }
  for (i in seq_along(training$history)) {
    if (length(saved$values) >= i && !is.null(saved$values[[i]])) next
    message(sprintf("  evaluating policy %d/%d", i, length(training$history)))
    policy <- mab_unpack_policy(training$history[[i]]$policy)
    if (task$d <= 2L) value <- evaluate(policy) else {
      pending <- saved$pending
      if (is.null(pending)) pending <- list(policy = i, probability = rep(NA_real_, n), next_node = 1L)
      stopifnot(pending$policy == i)
      while (pending$next_node <= n) {
        ii <- seq.int(pending$next_node, min(n, pending$next_node + 255L))
        pp <- policy$eval_batch(X[ii, , drop = FALSE])
        if (any(!is.finite(pp)) || any(pp < -MAB_NUMERICAL_TOLERANCE |
            pp > 1 + MAB_NUMERICAL_TOLERANCE)) stop("Invalid evaluation probabilities")
        pending$probability[ii] <- clip01(pp)
        pending$next_node <- max(ii) + 1L
        saved$pending <- pending
        mab_save_rds(saved, checkpoint_file)
      }
      # Retain node order and sum only after the entire probability vector exists.
      pp <- pending$probability
      value <- c(value = sum(weights * (mu2 + gap*pp)), regret = sum(weights * gap * (op - pp)))
    }
    saved$values[[i]] <- value
    saved$pending <- NULL
    mab_save_rds(saved, checkpoint_file)
  }
  values <- do.call(rbind, saved$values)
  total <- sum(duration * values[, "regret"])
  initial_duration <- pmax(0, pmin(starts + duration - 1L, training$initial_rounds) - starts + 1L)
  list(total_regret = total, relative_regret = total / (task$T * oracle_value),
    oracle_value = oracle_value, initial_regret = sum(initial_duration * values[, "regret"]),
    final_reward_gap = tail(values[, "regret"], 1L), method = method,
    eval_n = if (task$d >= 3L) n else 0L,
    decomposition = data.frame(start = starts, duration = duration,
      oracle_value = oracle_value, policy_value = values[, "value"], regret_rate = values[, "regret"]))
}

mab_simulation_collect <- function(output_dir) {
  paths <- list.files(file.path(output_dir, "runs"), "^result\\.rds$", recursive = TRUE, full.names = TRUE)
  rows <- lapply(paths, function(path) {
    result <- readRDS(path)
    task <- as.data.frame(result$task, stringsAsFactors = FALSE)
    if (result$task$algorithm != "original") return(NULL)
    cbind(task, data.frame(total_regret = result$evaluation$total_regret,
      relative_regret = result$evaluation$relative_regret,
      oracle_value = result$evaluation$oracle_value, final_reward_gap = result$evaluation$final_reward_gap))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

mab_simulation_summary <- function(runs) {
  if (!nrow(runs)) return(data.frame())
  keys <- c("d", "covariate", "regime", "L", "T")
  groups <- split(seq_len(nrow(runs)), interaction(runs[keys], drop = TRUE))
  do.call(rbind, lapply(groups, function(ii) {
    z <- runs[ii, , drop = FALSE]
    cbind(z[1L, keys, drop = FALSE], data.frame(n = nrow(z),
      mean_cumulative = mean(z$total_regret), sd_cumulative = sd(z$total_regret),
      mean_relative = mean(z$relative_regret), sd_relative = sd(z$relative_regret)))
  }))
}

# Plot means with shading for one sample standard deviation.
mab_release_regret_plot <- function(data, metric = c("cumulative", "relative"), series = "d", facet = TRUE) {
  metric <- match.arg(metric)
  stopifnot(series %in% c("d", "L"))
  levels <- if (series == "d") as.character(1:4) else as.character(c(.25, .5, 1, 2, 4))
  colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")[seq_along(levels)], levels)
  shapes <- setNames((if (series == "d") c(16, 24, 22, 18) else c(16, 17, 15, 18, 8)), levels)
  sizes <- setNames(if (series == "d") c(1.7, 1.15, 1.15, 1.7) else rep(1.7, 5), levels)
  data$series <- factor(as.character(data[[series]]), levels = levels)
  data$value <- data[[paste0("mean_", metric)]]
  data$lower_sd <- data$value - data[[paste0("sd_", metric)]]
  data$upper_sd <- data$value + data[[paste0("sd_", metric)]]
  if (facet) {
    data$covariate <- factor(data$covariate, levels = c("uniform", "truncated_gaussian"))
    data$regime <- factor(data$regime, levels = c("homogeneous", "heterogeneous"))
  }
  if (metric == "cumulative") data$lower_sd[data$lower_sd <= 0] <- NA_real_
  groups <- interaction(data$series, data$covariate, data$regime, drop = TRUE)
  # Never connect across an uncomputed intermediate horizon.
  horizons <- sort(unique(c(as.integer(2^(10:15)), data$T)))
  data$line_id <- NA_character_
  for (g in unique(groups)) {
    ii <- which(groups == g); ii <- ii[order(data$T[ii])]
    data$line_id[ii] <- paste(g, cumsum(c(TRUE, diff(match(data$T[ii], horizons)) != 1L)), sep = "|")
  }
  markers <- do.call(rbind, lapply(split(data, groups), function(z) {
    z <- z[order(z$T), ]
    if (series == "d" && z$d[1] == 2L) z <- z[seq.int(1L, nrow(z), by = 2L), , drop = FALSE]
    if (series == "d" && z$d[1] == 3L && nrow(z) > 1L) z <- z[seq.int(2L, nrow(z), by = 2L), , drop = FALSE]
    z
  }))
  markers$outline <- if (series == "d") ifelse(markers$d %in% c(2L, 3L), .35, .5) else .5
  lines <- data[as.integer(table(data$line_id)[data$line_id]) > 1L, , drop = FALSE]
  p <- ggplot2::ggplot(data, ggplot2::aes(T, value, colour = series, shape = series, fill = series)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_sd, ymax = upper_sd, group = line_id),
      colour = NA, alpha = .15, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(data = lines, ggplot2::aes(group = line_id), linewidth = .75) +
    ggplot2::geom_point(data = markers, ggplot2::aes(size = series, stroke = I(outline))) +
    ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(values = shapes) + ggplot2::scale_size_manual(values = sizes, guide = "none") +
    ggplot2::labs(colour = series, shape = series, fill = series) + mab_release_plot_theme()
  if (facet) p <- p + ggplot2::facet_grid(covariate ~ regime, scales = "free_y", labeller = ggplot2::labeller(
    covariate = c(uniform = "Uniform", truncated_gaussian = "Truncated Gaussian"),
    regime = c(homogeneous = "Homogeneous", heterogeneous = "Heterogeneous")))
  loglabels <- function(values) sub("\\.?0+$", "", sprintf("%.2f", log10(values)))
  if (metric == "cumulative") p <- p +
    ggplot2::scale_x_log10(breaks = sort(unique(data$T)), labels = loglabels,
      expand = ggplot2::expansion(mult = c(.06, .08))) +
    ggplot2::scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000), labels = loglabels,
      expand = ggplot2::expansion(mult = c(.05, .07))) +
    ggplot2::labs(x = expression(log[10](T)), y = expression(log[10](plain("Cumulative regret")))) else p <- p +
    ggplot2::scale_x_continuous(breaks = if (max(data$T) >= 1024L) 8192*(0:4) else sort(unique(data$T)),
      labels = scales::label_comma(accuracy = 1), limits = c(0, max(data$T)),
      expand = ggplot2::expansion(mult = c(.02, .04))) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = .01),
      expand = ggplot2::expansion(mult = c(.05, .07))) +
    ggplot2::labs(x = "T", y = "Relative regret")
  present <- levels[levels %in% as.character(unique(data$series))]
  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(
    size = unname(sizes[present]), fill = unname(colors[present]),
    stroke = if (series == "d") ifelse(present %in% c("2", "3"), .35, .5) else .5)))
  p
}

mab_simulation_export <- function(plot, stem, width = 12, height = 12*1230/1788) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 180, bg = "white")
}

mab_policy_figure <- function(data, output_dir) {
  kinds <- c("Oracle fair policy", "Our algorithm", "Standard UCB")
  colors <- setNames(c("#0072B2", "#D55E00", "#009E73"), kinds)
  times <- sort(unique(data$T))
  draw <- function() {
    graphics::par(mfrow = c(length(times), 3), mar = c(2.1, 1.8, 2.7, 1.2),
      oma = c(.3, .2, 1.2, .2), mgp = c(1.8, .5, 0), family = "Helvetica")
    for (T in times) for (kind in kinds) {
      q <- data[data$T == T & data$policy == kind, , drop = FALSE]
      if (!nrow(q)) { graphics::plot.new(); next }
      x <- sort(unique(q$x1)); y <- sort(unique(q$x2))
      values <- gap <- matrix(NA_real_, length(x), length(y))
      ii <- cbind(match(q$x1, x), match(q$x2, y))
      values[ii] <- q$probability; gap[ii] <- q$true_gap
      tr <- graphics::persp(x, y, values, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1),
        theta = 35, phi = 25, expand = .8, scale = FALSE, col = colors[[kind]],
        border = grDevices::adjustcolor("#333333", .20), shade = .25, ltheta = 120, lphi = 35,
        axes = FALSE, xlab = "", ylab = "", zlab = "", box = TRUE)
      ends <- rbind(c(0, .82), c(1, .55))
      plane <- grDevices::trans3d(c(ends[, 1], rev(ends[, 1])),
        c(ends[, 2], rev(ends[, 2])), c(0, 0, 1, 1), tr)
      graphics::polygon(plane$x, plane$y, col = grDevices::adjustcolor("#D62728", .22),
        border = "#C61C30", lwd = .9)
      span <- diff(graphics::par("usr")[1:2])
      label_at <- function(x, y, z, label, rotation = 0, dx = 0) {
        a <- grDevices::trans3d(x, y, z, tr)
        graphics::text(a$x + dx, a$y, labels = label, srt = rotation, cex = 1.05, xpd = NA)
      }
      for (tick in seq(.2, .8, by = .2)) {
        a <- grDevices::trans3d(tick, 0, 0, tr); b <- grDevices::trans3d(tick, -.035, 0, tr)
        graphics::segments(a$x, a$y, b$x, b$y)
        a <- grDevices::trans3d(tick, -.09, 0, tr)
        graphics::text(a$x, a$y, sprintf("%.1f", tick), cex = .82, xpd = NA)
        a <- grDevices::trans3d(1, tick, 0, tr); b <- grDevices::trans3d(1.035, tick, 0, tr)
        graphics::segments(a$x, a$y, b$x, b$y)
        a <- grDevices::trans3d(1.085, tick, 0, tr)
        graphics::text(a$x, a$y, sprintf("%.1f", tick), cex = .82, xpd = NA)
      }
      for (tick in seq(0, 1, by = .2)) {
        a <- grDevices::trans3d(0, 0, tick, tr)
        graphics::segments(a$x, a$y, a$x - .024*span, a$y)
        graphics::text(a$x - .034*span, a$y, sprintf("%.1f", tick), adj = 1, cex = .82, xpd = NA)
      }
      label_at(.50, -.24, 0, expression(x[1]), -28)
      label_at(1.23, .50, 0, expression(x[2]), 52)
      label_at(0, 0, .50, expression(p[1](x)), 90, -.12*span)
      graphics::title(main = kind, font.main = 2, cex.main = 1.12, line = .8)
      if (kind == "Our algorithm") graphics::mtext(paste0("T = ", format(T, big.mark = ",")),
        side = 3, line = 2.25, cex = 1.05)
    }
  }
  grDevices::pdf(file.path(output_dir, "policy.pdf"), width = 15, height = 5.4*length(times), family = "Helvetica")
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(file.path(output_dir, "policy.png"), width = 15, height = 5.4*length(times), units = "in", res = 180)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(NULL)
}

mab_simulation_report <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("filelock", quietly = TRUE)) stop("Install filelock")
  lock <- filelock::lock(file.path(output_dir, ".report.lock"), timeout = 0)
  if (is.null(lock)) return(invisible(NULL))
  on.exit(filelock::unlock(lock), add = TRUE)
  runs <- mab_simulation_collect(output_dir)
  summary_file <- file.path(output_dir, "summary.csv")
  if (nrow(runs)) {
    write.csv(runs, file.path(output_dir, "replication_results.csv"), row.names = FALSE)
    summary <- mab_simulation_summary(runs)
    write.csv(summary, summary_file, row.names = FALSE)
  } else if (file.exists(summary_file)) {
    summary <- read.csv(summary_file)
  } else return(invisible(NULL))
  grids <- list.files(file.path(output_dir, "runs"), "^policy_grid\\.rds$", recursive = TRUE, full.names = TRUE)
  grid_file <- file.path(output_dir, "policy_grid.csv")
  policy_grid <- if (length(grids)) unique(do.call(rbind, lapply(grids, readRDS))) else
    if (file.exists(grid_file)) read.csv(grid_file) else data.frame()
  if (nrow(policy_grid)) write.csv(policy_grid, grid_file, row.names = FALSE)
  baseline <- summary[summary$L == .5 & summary$regime %in% c("homogeneous", "heterogeneous"), ]
  vary <- summary[summary$d == 2L & summary$covariate == "uniform" & summary$regime == "homogeneous", ]
  for (metric in c("cumulative", "relative")) {
    if (nrow(baseline)) mab_simulation_export(mab_release_regret_plot(baseline, metric),
      file.path(output_dir, paste0("ordinary_", metric)))
    if (nrow(vary)) mab_simulation_export(mab_release_regret_plot(vary, metric, series = "L", facet = FALSE),
      file.path(output_dir, paste0("vary_L_", metric)), 8, 5.2)
  }
  slopes <- lapply(split(baseline, interaction(baseline$d, baseline$covariate, baseline$regime, drop = TRUE)), function(z) {
    z <- z[is.finite(z$mean_cumulative) & z$mean_cumulative > 0, ]
    if (nrow(z) < 2L) return(NULL)
    cbind(z[1L, c("d", "covariate", "regime")], data.frame(
      slope = unname(coef(lm(log(mean_cumulative) ~ log(T), data = z))[2]),
      horizons = nrow(z), smallest_n = min(z$n), largest_n = max(z$n)))
  })
  slopes <- Filter(Negate(is.null), slopes)
  if (length(slopes)) write.csv(do.call(rbind, slopes), file.path(output_dir, "loglog_slopes.csv"), row.names = FALSE)
  if (nrow(policy_grid)) mab_policy_figure(policy_grid, output_dir)
  invisible(summary)
}

mab_simulation_save_policy_grid <- function(task, training, output_dir) {
  if (task$regime != "heterogeneous_d2_arm1_small") return(invisible(NULL))
  axis <- seq(0, 1, length.out = 101L)
  grid <- as.matrix(expand.grid(x1 = axis, x2 = axis))
  oracle <- mab_ordinary_oracle(2L, task$regime, task$L)
  learner <- mab_unpack_policy(tail(training$history, 1L)[[1L]]$policy)
  kinds <- list()
  kinds[[if (task$algorithm == "original") "Our algorithm" else "Standard UCB"]] <- learner
  kinds[["Oracle fair policy"]] <- oracle
  gap <- grid[, 2] + .27*grid[, 1] - .82
  rows <- lapply(names(kinds), function(kind) data.frame(policy = kind, T = task$T,
    rep = task$rep, L = task$L, x1 = grid[, 1], x2 = grid[, 2],
    probability = kinds[[kind]]$eval_batch(grid), true_gap = gap))
  mab_save_rds(do.call(rbind, rows), file.path(output_dir, "runs", mab_simulation_id(task), "policy_grid.rds"))

}

mab_simulation_settings_file <- function(output_dir, cfg) {
  lock <- filelock::lock(file.path(output_dir, ".settings.lock"))
  on.exit(filelock::unlock(lock))
  path <- file.path(output_dir, "settings.rds")
  current <- list(version = 1L, settings = cfg,
    numerical_tolerance = MAB_NUMERICAL_TOLERANCE)
  previous <- mab_read_rds(path)
  extra <- if (is.null(previous)) list() else previous[setdiff(names(previous), names(current))]
  compatible <- is.null(previous) || (identical(previous[names(current)], current) &&
    length(extra) <= 1L && all(vapply(extra, identical, logical(1), FALSE)))
  if (!compatible)
    stop("This output directory contains different simulation settings; select a new --output directory")
  if (is.null(previous)) mab_save_rds(current, path)
}

run_simulation_studies <- function(options = mab_release_options("simulation")) {
  output_dir <- options$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (options$plot_only) {
    has_runs <- nrow(mab_simulation_collect(output_dir)) > 0L
    if (!has_runs && !file.exists(file.path(output_dir, "summary.csv")))
      stop("No completed simulation results; first run without --plot-only")
    return(invisible(mab_simulation_report(output_dir)))
  }
  cfg <- mab_simulation_settings()
  if (!requireNamespace("filelock", quietly = TRUE)) stop("Install filelock")
  mab_simulation_settings_file(output_dir, cfg)
  tasks <- mab_simulation_tasks()
  write.csv(tasks, file.path(output_dir, "experiment_plan.csv"), row.names = FALSE)
  chosen <- which((seq_len(nrow(tasks)) - 1L) %% options$worker_count == options$worker_id - 1L)
  completed_now <- 0L
  on.exit(try(mab_simulation_report(output_dir), silent = TRUE), add = TRUE)
  for (i in chosen) {
    task <- as.list(tasks[i, ]); id <- task$id; task$id <- NULL
    run_dir <- file.path(output_dir, "runs", id)
    result_file <- file.path(run_dir, "result.rds")
    if (file.exists(result_file)) next
    if (completed_now >= options$max_tasks) break
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    task_lock <- filelock::lock(file.path(run_dir, ".task.lock"), timeout = 0)
    if (is.null(task_lock)) next
    tryCatch({
    if (file.exists(result_file)) next
    message(sprintf("Simulation %d/%d: %s", i, nrow(tasks), id))
    stream <- mab_simulation_stream(task, cfg)
    training_file <- file.path(run_dir, "training.rds")
    training <- mab_read_rds(training_file)
    if (is.null(training)) {
      training <- mab_train(stream$X, stream$y1, stream$y2, algorithm = task$algorithm,
        L = task$L, nu = 1, lambda = cfg$lambda, sigma = cfg$sigma,
        eta = cfg$eta, C_beta = cfg$C_beta, action_uniforms = stream$uniforms,
        checkpoint_file = file.path(run_dir, "checkpoint.rds"), checkpoint_every = 128L)
      mab_save_rds(training, training_file)
    }
    evaluation <- if (task$algorithm == "original") mab_simulation_evaluate(training, task, cfg,
      file.path(run_dir, "evaluation_checkpoint.rds")) else NULL
    mab_simulation_save_policy_grid(task, training, output_dir)
    mab_save_rds(list(task = task, seeds = stream$seeds, settings = cfg,
      numerical_tolerance = MAB_NUMERICAL_TOLERANCE, evaluation = evaluation), result_file)
    completed_now <- completed_now + 1L
    mab_simulation_report(output_dir)
    message(sprintf("Saved %s; %d tasks completed in this invocation", id, completed_now))
    }, finally = filelock::unlock(task_lock))
  }
  invisible(mab_simulation_report(output_dir))
}

if (!identical(getOption("mab.release.autorun", TRUE), FALSE)) run_simulation_studies()
