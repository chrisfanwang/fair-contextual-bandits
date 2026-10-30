# Adult and Jobs: fixed preprocessing, 50 replications, and online accuracy.
# Run: Rscript real_data_analysis.R
# Plot completed runs: Rscript real_data_analysis.R --plot-only

.real_entry_file <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  path<-if(length(paths))tail(paths,1)[[1]] else arg[[1]]
  if(!file.exists(path))path<-gsub("~+~"," ",path,fixed=TRUE)
  normalizePath(path,mustWork=TRUE)
})
source(file.path(dirname(.real_entry_file), "algorithms.R"))

real_with_seed <- function(seed, fun) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if(had_seed)old_seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if(had_seed)assign(".Random.seed",old_seed,envir=.GlobalEnv) else
      if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fun()
}

real_seed <- function(dataset, replication, kind, base_seed=20260924L) {
  ds <- match(dataset,c("adult","jobs"))
  k <- match(kind,c("evaluation","permutation","action"))
  stopifnot(!is.na(ds),!is.na(k))
  as.integer((as.double(base_seed)+ds*10000019+replication*100003+k*7919)%%2147483646+1)
}

.real_preprocessing_python <- r"---("""Prepare Adult and Jobs using the paper's fixed, label-free PCA split."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import urllib.request
import zipfile

import numpy as np
import pandas as pd

SOURCES = json.loads(r'''{
  "verified_on": "2026-10-30",
  "adult": {
    "name": "Adult",
    "creator": "Barry Becker and Ronny Kohavi",
    "year": 1996,
    "url": "https://archive.ics.uci.edu/dataset/2/adult",
    "filename": "adult.zip",
    "sha256": "7537312dd56c2b98035880805ce99e68183a30ee468aa5329d6df0fbb3cc21bb",
    "md5": "1252efa3d4fd9b7a70b3460adb3deb28",
    "license": "CC-BY-4.0",
    "bundled_records": false,
    "reference_array_md5": {
      "label_target.npy": "82dc37f2c3544038d63d818ec7eb2dbf",
      "raw_id_target.npy": "e81cc7000880b633b162d52efe56dfbf",
      "preprocess_source_2.npy": "ba6b555a4a615bd04eb82ef0106dceb0",
      "raw_id_source_6.npy": "2e437d74f2b9ec64098403015ed02d80",
      "preprocess_target.npy": "ca8d48dc919a930b85c2fb48323be4dd",
      "preprocess_source_3.npy": "7c94f1c61f202abed5cbddfbfc17b611",
      "preprocess_source_1.npy": "f3665106526ef4e83791783380c6069a",
      "raw_id_source_4.npy": "006f6bc98bdfaea64a233e3c27ebd999",
      "raw_id_source_5.npy": "c25c4d1f48adc75bd6df2ded957c6bd2",
      "preprocess_source_0.npy": "0cbe5aa2051564350e44001c00cacdd1",
      "preprocess_source_4.npy": "b2db26c27d0fb95fbc2899508db6c601",
      "raw_id_source_1.npy": "a81124f2643bc8b5491adcdceff2e5c7",
      "raw_id_source_0.npy": "6feec11a7ab1a925aa6c160d5a0091a1",
      "preprocess_source_5.npy": "b15d238646ee7995f0a20e3a379c88e5",
      "raw_id_source_2.npy": "ee3454e063da620aed6fd0276cc3b0d9",
      "raw_id_source_3.npy": "8bbb0f3c5da472aa03ec34c288f8f1a3",
      "preprocess_source_6.npy": "dae90a68b544fe4edf799385718c4a7f",
      "X_source_6.npy": "45316b0c2f554c5456ec26a1169f6bee",
      "label_source_5.npy": "b652c47db7199f5e1d6cd2a5e3c7d4ad",
      "label_source_4.npy": "9bd64fee7f51c7f29c4e1675272a43c4",
      "X_source_5.npy": "ac0c955b9e470cdba05b9591ae083d45",
      "label_source_6.npy": "21d27edc7c36b9d607a5140c3743b25d",
      "X_source_4.npy": "ed25461986f2e4063aaa8b6f683e70f0",
      "X_source_0.npy": "5b06a67c932b7f6eeb114830235b4b08",
      "label_source_3.npy": "ef2b3786971428d57898960ec295cd49",
      "label_source_2.npy": "69ebb09683481c07d3962ab868035f3e",
      "X_source_1.npy": "a836142795ef756485a89683e31ed706",
      "X_source_3.npy": "598f527da1906fabdc40137b240347fc",
      "label_source_0.npy": "eaeb27ea62ebf13a2e2cca8d6d4b3f23",
      "X_target.npy": "e6f45497f299756dbdb9c3c40e2d8f39",
      "label_source_1.npy": "23a7e7616e6321b99c6d946bf3b73b6c",
      "X_source_2.npy": "a5a0e878a56f703092a45387818958ab"
    },
    "download_url": "https://archive.ics.uci.edu/static/public/2/adult.zip",
    "doi": "10.24432/C5XW20",
    "license_url": "https://creativecommons.org/licenses/by/4.0/"
  },
  "jobs": {
    "name": "70k+ Job Applicants Data (Human Resource)",
    "creator": "Ayush Tankha",
    "year": 2023,
    "url": "https://www.kaggle.com/datasets/ayushtankha/70k-job-applicants-data-human-resource",
    "filename": "stackoverflow_full.csv",
    "sha256": "5947c0cd2f527716700c2f71c3ac32b155a63b5d7ba88adc4887b8d9f6c6cf38",
    "md5": "8f305383183bbe4fb9bec935f2944a8e",
    "license": "Unknown (as reported by the original Kaggle dataset metadata)",
    "bundled_records": false,
    "reference_array_md5": {
      "label_target.npy": "4ba7b5ec25b83d08ae76562246de68b2",
      "raw_id_target.npy": "6dfb3e0b93a2a831f2a03ea43cf6f972",
      "preprocess_target.npy": "1fe5286d6b7fc2083035eb8d4881821a",
      "preprocess_source_0.npy": "43d0e7613f4826a6e8a8d363c485ad4c",
      "raw_id_source_0.npy": "a4585c00c28ea3efe143bf20d9961c8c",
      "X_source_0.npy": "7db901ef5b16b3d4a58b7626feeec8ba",
      "label_source_0.npy": "0074ea6c99b70b6c0424c5ba59cf73c0",
      "X_target.npy": "f3da599f3cac9f48e5cbd144bb6a1b77"
    },
    "access": "Automatic public download; manual original-page download is a fallback if Kaggle requires sign-in.",
    "license_metadata_url": "https://www.kaggle.com/api/v1/datasets/view/ayushtankha/70k-job-applicants-data-human-resource",
    "download_url": "https://www.kaggle.com/api/v1/datasets/download/ayushtankha/70k-job-applicants-data-human-resource?datasetVersionNumber=1",
    "dataset_version": 1
  }
}
''')
COUNTRIES = ["United-States", "Mexico", "Philippines", "Germany", "Puerto-Rico",
             "Canada", "El-Salvador", "India"]
BASE_SEED = 20260720
DATASETS = ("adult", "jobs")


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_adult(path):
    columns = ["age", "workclass", "fnlwgt", "education", "education-num",
               "marital-status", "occupation", "relationship", "race", "sex",
               "capital-gain", "capital-loss", "hours-per-week", "native-country", "income"]
    with zipfile.ZipFile(path) as archive:
        frame = pd.concat([pd.read_csv(archive.open(name), names=columns,
                                      skipinitialspace=True, comment="|", na_values="?")
                           for name in ("adult.data", "adult.test")], ignore_index=True)
    frame = frame.dropna().copy()
    frame = frame[frame["native-country"].isin(COUNTRIES)].copy()
    labels = frame.income.str.rstrip(".").map({"<=50K": 0, ">50K": 1})
    groups = frame["native-country"].map({country: j for j, country in enumerate(COUNTRIES)})
    return frame.drop(columns=["native-country", "income", "education"]), labels, groups


def load_jobs(path):
    frame = pd.read_csv(path).drop(columns=["Unnamed: 0", "HaveWorkedWith"], errors="ignore")
    frame = frame.dropna().copy()
    labels = pd.to_numeric(frame["Employed"])
    # Group non-US records before US records for the seeded sampling order.
    groups = pd.Series(np.where(frame.Country == "United States of America", 1, 0), index=frame.index)
    return frame.drop(columns=["Country", "Employed"]), labels, groups


def fit_schema(frame):
    numeric = list(frame.select_dtypes(include="number").columns)
    categories = {name: sorted(frame[name].astype(str).unique().tolist())
                  for name in frame.columns if name not in numeric}
    return {"numeric": numeric, "categories": categories}


def encode(frame, schema):
    pieces = [frame[schema["numeric"]].to_numpy(dtype=float)]
    for name, levels in schema["categories"].items():
        values = frame[name].astype(str).to_numpy()
        pieces.append(np.column_stack([values == level for level in levels]))
    return np.column_stack(pieces)


def npy_md5(array):
    stream = io.BytesIO()
    np.save(stream, array)
    return hashlib.md5(stream.getvalue()).hexdigest()


def prepare(dataset, raw_dir, output_dir, download=True):
    source = SOURCES[dataset]
    raw = raw_dir / source["filename"]
    if not raw.exists() and download:
        raw.parent.mkdir(parents=True, exist_ok=True)
        temporary = raw.with_suffix(raw.suffix + ".download")
        try:
            urllib.request.urlretrieve(source["download_url"], temporary)
            if dataset == "jobs":
                with zipfile.ZipFile(temporary) as archive:
                    names = [name for name in archive.namelist()
                             if Path(name).name == source["filename"]]
                    if len(names) != 1:
                        raise ValueError("The Jobs archive must contain exactly one stackoverflow_full.csv.")
                    content = archive.read(names[0])
                temporary.write_bytes(content)
            temporary.replace(raw)
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f"Could not download {dataset}: {exc}. "
                               f"Download {source['filename']} manually from {source['url']} "
                               f"and save it as {raw}.") from exc
    if not raw.exists():
        raise FileNotFoundError(f"Missing {raw}. Obtain {source['filename']} from "
                                f"{source['url']}; see README.md.")
    raw_hash = digest(raw)
    if raw_hash != source["sha256"]:
        raise ValueError(f"{raw.name} differs from the data version used in the paper. "
                         f"Expected SHA256 {source['sha256']}; obtained {raw_hash}.")
    destination = output_dir / dataset
    manifest_path = destination / "manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old.get("raw_sha256") == raw_hash and old.get("script_sha256") == digest(__file__):
            if all((destination / name).exists() and digest(destination / name) == value
                   for name, value in old["output_sha256"].items()):
                print(f"{dataset}: prepared data already present", flush=True)
                return old
        raise FileExistsError(f"{destination} contains another preprocessing build; "
                              "use --output to choose a new directory.")

    frame, labels, groups = (load_adult if dataset == "adult" else load_jobs)(raw)
    if labels.isna().any() or set(labels.unique()) != {0, 1}:
        raise ValueError("The response must contain exactly the labels 0 and 1.")
    raw_ids = frame.index.to_numpy(dtype=np.float64)
    frame = frame.reset_index(drop=True)
    labels = labels.to_numpy(dtype=np.float64)
    groups = groups.to_numpy()
    n = len(frame)
    n_pre = int(np.floor(.1 * n))
    actual_seed = BASE_SEED + DATASETS.index(dataset) + 1
    rng = np.random.default_rng(actual_seed)
    indices = np.sort(rng.choice(n, n_pre, replace=False))
    mask = np.zeros(n, dtype=np.float64)
    mask[indices] = 1
    training = frame.iloc[indices]
    schema = fit_schema(training)
    encoded = encode(training, schema)
    center = encoded.mean(axis=0)
    scale = encoded.std(axis=0, ddof=1)
    scale[scale == 0] = 1
    standardized = (encoded - center) / scale
    covariance = standardized.T @ standardized / (n_pre - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    rotation = eigenvectors[:, order[:3]]
    rotation *= np.where(rotation[np.argmax(abs(rotation), axis=0), np.arange(3)] < 0, -1, 1)
    training_scores = standardized @ rotation
    lower = training_scores.min(axis=0)
    width = np.ptp(training_scores, axis=0)
    width[width == 0] = 1
    scores = np.empty((n, 3), dtype=np.float64)
    clipped = 0
    for first in range(0, n, 65536):
        batch = encode(frame.iloc[first:first + 65536], schema)
        transformed = (((batch - center) / scale) @ rotation - lower) / width
        clipped += np.count_nonzero((transformed < 0) | (transformed > 1))
        scores[first:first + len(transformed)] = np.clip(transformed, 0, 1)

    # Group ordering affects which records are drawn by R's seeded sample.int.
    pooled_order = np.concatenate([np.flatnonzero(groups == g) for g in sorted(set(groups))])
    records = pd.DataFrame({"label": labels[pooled_order].astype(int)})
    records["preprocessing"] = mask[pooled_order].astype(int)
    records["raw_id"] = raw_ids[pooled_order].astype(int)
    destination.mkdir(parents=True, exist_ok=True)
    # Store column-major little-endian float64 for R's readBin.
    scores[pooled_order].T.astype("<f8").tofile(destination / "covariates.bin")
    records.to_csv(destination / "records.csv.gz", index=False, float_format="%.17g",
                   compression={"method": "gzip", "mtime": 0})
    np.savez(destination / "preprocessing_model.npz", center=center, scale=scale,
             rotation=rotation, lower=lower, width=width, numeric_limits=np.array([]))
    generated_md5 = {}
    for group in sorted(set(groups)):
        part = "target" if group == 0 else f"source_{group - 1}"
        take = groups == group
        for field, values in [("X", scores[take]), ("label", labels[take]),
                              ("preprocess", mask[take]), ("raw_id", raw_ids[take])]:
            generated_md5[f"{field}_{part}.npy"] = npy_md5(values)
    manifest = {
        "dataset": dataset, "n_total": n, "n_preprocessing": n_pre, "dimension": 3,
        "preprocess_fraction": .1, "preprocess_seed": BASE_SEED,
        "numpy_split_seed": actual_seed, "raw_sha256": raw_hash,
        "script_sha256": digest(__file__), "feature_columns": list(frame.columns),
        "feature_schema": schema, "encoded_dimension": encoded.shape[1],
        "explained_variance": (eigenvalues[order[:3]] / eigenvalues.sum()).tolist(),
        "clipped_fraction": clipped / (3 * n),
        "reference_array_md5": source["reference_array_md5"],
        "generated_array_md5": generated_md5,
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
        "output_sha256": {name: digest(destination / name)
                          for name in ("records.csv.gz", "covariates.bin", "preprocessing_model.npz")}}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{dataset}: {n:,} records; {n_pre:,} preprocessing records; three covariates", flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", choices=DATASETS, default=list(DATASETS))
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-download", action="store_true")
    args = parser.parse_args()
    for name in args.datasets:
        prepare(name, args.raw, args.output, download=not args.no_download)
)---"

real_prepare_data <- function(datasets=c("adult","jobs"), root=dirname(MAB_SOURCE_FILE),
                              raw_dir=getOption("mab.release.raw",file.path(root,"data","raw")),
                              output=getOption("mab.release.data",file.path(root,"data","processed")),
                              download=TRUE) {
  if(!length(datasets) || any(!datasets%in%c("adult","jobs")))
    stop("Choose adult, jobs, or both datasets.")
  datasets <- unique(datasets)
  if(!requireNamespace("filelock",quietly=TRUE))stop("Install filelock for data preparation.")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  dir.create(raw_dir,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output,mustWork=TRUE)
  raw_dir <- normalizePath(raw_dir,mustWork=TRUE)
  lock_dir <- file.path(root,"data")
  dir.create(lock_dir,recursive=TRUE,showWarnings=FALSE)
  lock <- filelock::lock(file.path(lock_dir,"preparation.lock"))
  on.exit(filelock::unlock(lock))
  missing <- datasets[!file.exists(file.path(output,datasets,"records.csv.gz")) |
                      !file.exists(file.path(output,datasets,"covariates.bin")) |
                      !file.exists(file.path(output,datasets,"manifest.json"))]
  if(length(missing)) {
    python <- Sys.getenv("MAB_PYTHON",unset="python3")
    script <- tempfile("prepare-real-data-",fileext=".py")
    writeLines(.real_preprocessing_python,script,useBytes=TRUE)
    on.exit(unlink(script),add=TRUE)
    args <- c(shQuote(script),"--datasets",missing,"--raw",shQuote(raw_dir),
              "--output",shQuote(output),if(!download)"--no-download")
    result <- system2(python,args)
    if(result!=0L)stop("Data preparation failed. See README.md for Python dependencies and data sources.")
  }
  output
}

real_load_data <- function(dataset, data_dir) {
  path <- file.path(data_dir,dataset,"records.csv.gz")
  frame <- utils::read.csv(path,check.names=FALSE)
  stopifnot(all(c("label","preprocessing","raw_id")%in%names(frame)),
            all(frame$label%in%c(0,1)),all(frame$preprocessing%in%c(0,1)))
  covariate_path <- file.path(data_dir,dataset,"covariates.bin")
  if(file.info(covariate_path)$size!=nrow(frame)*3L*8L)stop("Covariate file has the wrong size.")
  X <- matrix(readBin(covariate_path,what="double",n=nrow(frame)*3L,size=8L,
                      endian="little"),nrow=nrow(frame),ncol=3L)
  stopifnot(all(is.finite(X)),all(X>=0 & X<=1),!anyDuplicated(frame$raw_id))
  pre <- which(frame$preprocessing==1L)
  pool <- setdiff(seq_len(nrow(frame)),pre)
  holdout_seed <- real_seed(dataset,0L,"evaluation")
  selected <- real_with_seed(holdout_seed,function()sample.int(length(pool),512L))
  holdout <- pool[selected]
  online <- pool[-selected]
  stopifnot(!any(holdout%in%online),!any(online%in%pre),!any(holdout%in%pre))
  list(dataset=dataset,X=X[online,,drop=FALSE],label=frame$label[online],
       raw_id=frame$raw_id[online],heldout_X=X[holdout,,drop=FALSE],
       heldout_label=frame$label[holdout],heldout_raw_id=frame$raw_id[holdout],
       preprocessing_raw_id=frame$raw_id[pre],n_total=nrow(frame),
       n_preprocessing=length(pre),n_heldout=length(holdout),n_online=length(online),
       data_md5=paste(unname(tools::md5sum(c(path,covariate_path))),collapse=":"),holdout_seed=holdout_seed)
}

real_stream <- function(data, replication, T) {
  if(T>nrow(data$X))stop("The requested horizon exceeds the online sample size.")
  order <- real_with_seed(real_seed(data$dataset,replication,"permutation"),
                         function()sample.int(nrow(data$X),replace=FALSE))[seq_len(T)]
  uniforms <- real_with_seed(real_seed(data$dataset,replication,"action"),function()runif(T))
  list(X=data$X[order,,drop=FALSE],label=data$label[order],
       raw_id=data$raw_id[order],uniforms=uniforms)
}

# Optional evaluation of a saved fit on the separate 512-record sample.
# These expected policy accuracies are distinct from realised online accuracy.
real_heldout_accuracy <- function(fit,data) {
  accuracy <- vapply(fit$history,function(epoch) {
    p <- mab_unpack_policy(epoch$policy)$eval_batch(data$heldout_X)
    mean(ifelse(data$heldout_label==1L,p,1-p))
  },numeric(1))
  starts <- vapply(fit$history,function(epoch)epoch$start,integer(1))
  durations <- diff(c(starts,fit$T+1L))
  c(final=tail(accuracy,1L),time_average=sum(durations*accuracy)/fit$T)
}

real_summarise <- function(records) {
  if(!length(records))return(data.frame())
  rows <- do.call(rbind,lapply(records,function(x)data.frame(
    dataset=x$dataset,T=x$T,L=x$L,rep=x$rep,accuracy=x$accuracy)))
  key <- paste(rows$dataset,rows$T,rows$L,rows$rep,sep=":")
  if(anyDuplicated(key))stop("A replication is present more than once.")
  groups <- split(rows,interaction(rows$dataset,rows$T,rows$L,drop=TRUE))
  answer <- do.call(rbind,lapply(groups,function(x)data.frame(dataset=x$dataset[1],
    L=x$L[1],T=x$T[1],metric="online_accuracy",n=nrow(x),planned_n=50L,
    mean=mean(x$accuracy),SD=if(nrow(x)>1L)stats::sd(x$accuracy) else NA_real_,
    SE=if(nrow(x)>1L)stats::sd(x$accuracy)/sqrt(nrow(x)) else NA_real_,
    reps=paste(sort(x$rep),collapse=";"))))
  rownames(answer) <- NULL
  answer[order(answer$dataset,answer$L,answer$T),]
}

real_accuracy_plot <- function(summary, output_dir) {
  if(!requireNamespace("ggplot2",quietly=TRUE))stop("Install the R package ggplot2 to draw figures.")
  summary <- summary[summary$dataset%in%c("adult","jobs") & is.finite(summary$mean),]
  if(!nrow(summary))return(invisible(NULL))
  palette <- c("0.25"="#0072B2","0.5"="#D55E00","1"="#009E73","2"="#CC79A7","4"="#8E44AD")
  summary$dataset <- factor(summary$dataset,c("adult","jobs"),c("Adult","Jobs"))
  summary$L <- factor(as.character(summary$L),levels=names(palette))
  line_data <- summary[ave(summary$T,interaction(summary$dataset,summary$L),FUN=length)>1L,]
  p <- ggplot2::ggplot(summary,ggplot2::aes(x=T,y=mean,colour=L,fill=L,group=L)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=mean-SD,ymax=mean+SD),alpha=.13,colour=NA,na.rm=TRUE) +
    ggplot2::geom_line(data=line_data,linewidth=.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-SD,ymax=mean+SD),width=.055,linewidth=.35,alpha=.55,na.rm=TRUE) +
    ggplot2::geom_point(size=2.1,stroke=.3) +
    ggplot2::facet_wrap(~dataset,nrow=1,scales="free_y") +
    ggplot2::scale_colour_manual(values=palette,name="L",drop=FALSE) +
    ggplot2::scale_fill_manual(values=palette,name="L",drop=FALSE) +
    ggplot2::scale_x_continuous(trans="log2",breaks=sort(unique(summary$T)),
      labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)) +
    ggplot2::scale_y_continuous(labels=function(x)sprintf("%.0f%%",100*x)) +
    ggplot2::labs(x="Time horizon T",y="Cumulative online accuracy") +
    ggplot2::theme_bw(base_size=12,base_family="Helvetica") +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(colour="#DCDCDC",linewidth=.35),
      panel.border=ggplot2::element_rect(colour="black",fill=NA,linewidth=.5),
      strip.background=ggplot2::element_rect(fill="#EBEBEB",colour="black",linewidth=.5),
      strip.text=ggplot2::element_text(size=13),axis.title=ggplot2::element_text(size=13,colour="black"),
      axis.text=ggplot2::element_text(size=11,colour="black"),
      axis.text.x=ggplot2::element_text(size=10,angle=35,hjust=1),
      legend.position="right",legend.justification=c(0,1),
      legend.title=ggplot2::element_text(size=11),legend.text=ggplot2::element_text(size=11))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(file.path(output_dir,"real_accuracy.pdf"),p,width=10,height=4.5)
  ggplot2::ggsave(file.path(output_dir,"real_accuracy.png"),p,width=10,height=4.5,dpi=240)
  invisible(p)
}

real_update_report <- function(opts, expected_data=NULL, allow_saved_summary=FALSE) {
  lock <- filelock::lock(file.path(opts$output_dir,"report.lock"),timeout=0)
  if(is.null(lock))return(invisible(NULL))
  on.exit(filelock::unlock(lock))
  T_grid <- c(1024L,2048L,4096L,8192L,16384L,32768L)
  L_grid <- c(.25,.5,1,2,4)
  repetitions <- 50L
  files <- list.files(file.path(opts$output_dir,"runs"),pattern="^result[.]rds$",full.names=TRUE,recursive=TRUE)
  records <- lapply(files,readRDS)
  records <- Filter(function(x)x$dataset%in%c("adult","jobs") && x$T%in%T_grid &&
    x$L%in%L_grid && x$rep%in%seq_len(repetitions),records)
  for(dataset in unique(vapply(records,function(x)x$dataset,character(1)))) {
    same_data <- Filter(function(x)x$dataset==dataset,records)
    fingerprints <- unique(vapply(same_data,function(x)x$data_md5,character(1)))
    if(length(fingerprints)!=1L)stop("Cannot combine results from different preprocessing builds.")
    if(!is.null(expected_data) && !identical(fingerprints,expected_data[[dataset]]$data_md5))
      stop("An existing result was generated with different prepared data.")
  }
  settings <- list(algorithm="original",nu=1,lambda=1,sigma=.5,eta=.05,C_beta=2,
    exploration_multiplier=1,lp_neighbor_k=20L,numerical_tolerance=MAB_NUMERICAL_TOLERANCE,
    base_seed=20260924L,preprocess_seed=20260720L)
  for(record in records) {
    if(!identical(record$numerical_tolerance,MAB_NUMERICAL_TOLERANCE) ||
       !identical(record$identity[names(settings)],settings))
      stop("Cannot combine results generated with different experimental settings.")
  }
  summary <- real_summarise(records)
  path <- file.path(opts$output_dir,"online_accuracy_summary.csv")
  if(nrow(summary)) {
    summary$planned_n <- repetitions
    utils::write.csv(summary,path,row.names=FALSE)
  } else if(allow_saved_summary) {
    if(!file.exists(path))stop("No completed real-data results are available. Run real_data_analysis.R first.")
    summary <- utils::read.csv(path)
  }
  if(nrow(summary))real_accuracy_plot(summary,opts$output_dir)
  invisible(summary)
}

run_real_data_analysis <- function() {
  opts <- mab_release_options("real")
  dir.create(opts$output_dir,recursive=TRUE,showWarnings=FALSE)
  if(!requireNamespace("filelock",quietly=TRUE))stop("Install filelock for resumable experiments.")
  if(opts$plot_only) {
    return(real_update_report(opts,allow_saved_summary=TRUE))
  }
  T_grid <- c(1024L,2048L,4096L,8192L,16384L,32768L)
  L_grid <- c(.25,.5,1,2,4)
  repetitions <- 50L
  datasets <- c("adult","jobs")
  data_dir <- real_prepare_data(datasets)
  data <- setNames(lapply(datasets,real_load_data,data_dir=data_dir),datasets)
  sizes <- do.call(rbind,lapply(data,function(x)data.frame(dataset=x$dataset,
    eligible=x$n_total,preprocessing=x$n_preprocessing,heldout=x$n_heldout,online_pool=x$n_online)))
  size_lock <- filelock::lock(file.path(opts$output_dir,"report.lock"),timeout=0)
  if(!is.null(size_lock))tryCatch(
    utils::write.csv(sizes,file.path(opts$output_dir,"data_sizes.csv"),row.names=FALSE),
    finally=filelock::unlock(size_lock))
  tasks <- expand.grid(dataset=datasets,T=T_grid,L=L_grid,rep=seq_len(repetitions),stringsAsFactors=FALSE)
  tasks <- tasks[order(tasks$rep,tasks$T,tasks$L,tasks$dataset),]
  tasks$id <- seq_len(nrow(tasks))
  todo <- tasks[(tasks$id-1L)%%opts$worker_count==opts$worker_id-1L,]
  count <- 0L
  for(i in seq_len(nrow(todo))) {
    task <- todo[i,]
    directory <- file.path(opts$output_dir,"runs",sprintf("%s_T%d_L%s_rep%02d",task$dataset,
      task$T,format(task$L,scientific=FALSE,trim=TRUE),task$rep))
    result_path <- file.path(directory,"result.rds")
    identity <- list(dataset=task$dataset,T=task$T,L=task$L,rep=task$rep,
      data_md5=data[[task$dataset]]$data_md5,algorithm="original",nu=1,lambda=1,
      sigma=.5,eta=.05,C_beta=2,exploration_multiplier=1,lp_neighbor_k=20L,
      numerical_tolerance=MAB_NUMERICAL_TOLERANCE,base_seed=20260924L,preprocess_seed=20260720L)
    dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    lock <- filelock::lock(file.path(directory,"task.lock"),timeout=0)
    if(is.null(lock))next
    state <- tryCatch({
      saved <- mab_read_rds(result_path)
      if(!is.null(saved)) {
        if(!identical(saved$identity,identity))stop("Data or experimental settings changed for an existing result.")
        "existing"
      } else if(count>=opts$max_tasks) "limit" else {
        stream <- real_stream(data[[task$dataset]],task$rep,task$T)
        message(sprintf("%s: T=%d, L=%s, replication %d/%d",task$dataset,task$T,task$L,task$rep,repetitions))
        fit <- mab_train(stream$X,stream$label,1-stream$label,algorithm="original",L=task$L,
          nu=1,lambda=1,sigma=.5,eta=.05,C_beta=2,action_uniforms=stream$uniforms,
          exploration_multiplier=1,lp_neighbor_k=20L,checkpoint_file=file.path(directory,"checkpoint.rds"),
          progress_file=file.path(directory,"progress.rds"))
        answer <- list(dataset=task$dataset,T=task$T,L=task$L,rep=task$rep,
          accuracy=mean(fit$actions$reward),data_md5=data[[task$dataset]]$data_md5,identity=identity,
          numerical_tolerance=MAB_NUMERICAL_TOLERANCE,
          seeds=list(preprocessing=20260720L,evaluation=data[[task$dataset]]$holdout_seed,
            permutation=real_seed(task$dataset,task$rep,"permutation"),action=real_seed(task$dataset,task$rep,"action")))
        if(isTRUE(getOption("mab.release.heldout",FALSE)))answer$heldout <- real_heldout_accuracy(fit,data[[task$dataset]])
        mab_save_rds(answer,result_path)
        "completed"
      }
    },finally=filelock::unlock(lock))
    if(state=="limit")break
    if(state=="completed") {
      count <- count+1L
      real_update_report(opts,expected_data=data)
    }
  }
  files <- list.files(file.path(opts$output_dir,"runs"),pattern="^result[.]rds$",full.names=TRUE,recursive=TRUE)
  summary <- real_update_report(opts,expected_data=data)
  message(sprintf("Saved %d completed real-data experiments; figures use their actual replication counts.",length(files)))
  invisible(summary)
}

if(!isFALSE(getOption("mab.release.autorun",TRUE)))run_real_data_analysis()
