# Small dependency-free CLI. Source this before loading analysis packages.
parse_cli <- function(description, value_options = character(), flag_options = character()) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat(description, "\n\nOptions (values use separate arguments):\n",
        "  --project-root DIR  (or QBC_PROJECT_ROOT)\n",
        paste0("  --", c("input-dir", "output-dir", value_options), " VALUE", collapse = "\n"),
        "\n", paste0("  --", flag_options, collapse = "\n"), "\n", sep = "")
    quit(status = 0)
  }
  out <- list("project-root" = Sys.getenv("QBC_PROJECT_ROOT", unset = ""))
  allowed <- c("project-root", "input-dir", "output-dir", value_options)
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--")) stop("Expected --option, got: ", args[i], call. = FALSE)
    if (key %in% flag_options) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      if (!key %in% allowed) stop("Unknown option: ", args[i], call. = FALSE)
      if (i == length(args) || startsWith(args[i+1], "--"))
        stop("Missing value for: ", args[i], call. = FALSE)
      out[[key]] <- args[i+1]
      i <- i + 2
    }
  }
  out
}

option_path <- function(opt, key, subdir = NULL, input = TRUE) {
  path <- opt[[key]]
  if (is.null(path) && !is.null(subdir) && nzchar(opt[["project-root"]]))
    path <- file.path(opt[["project-root"]], subdir)
  if (is.null(path) || !nzchar(path))
    stop("Required: --", key, " (or --project-root/QBC_PROJECT_ROOT when supported)", call. = FALSE)
  if (input && !file.exists(path)) stop("Input does not exist: ", path, call. = FALSE)
  if (!input) dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

species_options <- function(opt, allow_both = TRUE) {
  sp <- opt[["species"]]
  if (is.null(sp)) sp <- if (allow_both) "both" else ""
  if (allow_both && sp == "both") return(c("human", "rhesus"))
  if (!sp %in% c("human", "rhesus")) stop("--species must be human or rhesus", call. = FALSE)
  sp
}

neutral_donors <- function(values, species, modality) {
  if (anyNA(values)) stop("Missing donor values")
  values <- as.character(values)
  prefix <- paste0(species, "_", modality, "_")
  neutral <- grepl(paste0("^", prefix, "[0-9]{2,}$"), values)
  if (all(neutral)) return(values)
  if (any(neutral)) stop("Mixed neutral and original donor labels")
  old_locale <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_locale), add = TRUE)
  Sys.setlocale("LC_COLLATE", "C")
  lev <- sort(unique(values))
  unname(setNames(sprintf("%s%02d", prefix, seq_along(lev)), lev)[values])
}
