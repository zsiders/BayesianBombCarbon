#' Load stan model
#'
#' @description Loads a specific stan model for model fitting in `BayesianBombCarbon`.
#'
#' @param model_metainfo A `list` containing meta information about the specified model to be fitted. The required stan model is automatically inferred from the `model_metainfo`.
#' @param model_filename File name of a specific stan model to load. This is an alternative to supplying `model_metainfo`.
#' @param model_folder Path to the folder containing the stan models for `BayesianBombCarbon`.
#' @param profile Should profiling be run during model fitting? Default is `TRUE`. If `FALSE`, will remove all profiling statements from the model before fitting. This can decrease runtime in some cases.
#' @param force_recompile If `FALSE` (default), the model is only recompiled if changes to the model code are detected. However, as the change detection is not fully reliable, it is sometimes necessary to force recompilation after having made changes to the stan code.
#' @param package Name of the package in which to search for stan files (defaults to BayesianBombCarbon). If NULL, will search in the normal working directory.
#'
#' @return A `list` containing the model definition and a link to the compiled stan model.
#' @export
#' @keywords internal
get_stan_model <- function(
    model_name = NULL,
    model_filename = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE,
    package = "BayesianBombCarbon") {
  model_stan <- list()

  if (is.null(model_filename)) {
    if (is.null(model_name)) {
      cli::cli_abort(
        c(
          "Please either provide ",
          "`model_filename` and `model_folder` (i.e. path to a stan model) or",
          "`model_name` (to infer a suitable stan model)."
        )
      )
    }
    if (model_name %in% c("ref-only")) {
      model_filename <- "delta14C_spline.stan"
    } else if (model_name == "integrated") {
      model_filename <- "delta14C_bspline_shift_wt.stan"
    } else {
      cli::cli_abort(
        paste(
          "No suitable model available,",
          "please supply a stan model yourself",
          "in `model_stan_opts`, or open an issue."
        )
      )
    }
  }

  model_stan[["model_filename"]] <- model_filename
  model_stan[["model_folder"]] <- model_folder
  model_stan[["force_recompile"]] <- force_recompile
  model_stan[["threads"]] <- threads
  model_stan[["profile"]] <- profile
  model_stan[["package"]] <- package

  model_stan <- update_compiled_stanmodel(model_stan, force_recompile)

  return(model_stan)
}

update_compiled_stanmodel <- function(model_stan, force_recompile = FALSE) {
  n_models <- length(model_stan$model_filename)

  model_path_list <- lapply(1:n_models, function(i) {
    if (is.null(model_stan$package)) {
      file.path(
        model_stan$model_folder,
        model_stan$model_filename[[i]]
      )
    } else {
      system.file(
        model_stan$model_folder,
        model_stan$model_filename[[i]],
        package = model_stan$package
      )
    }

  })
  include_paths_list <- lapply(1:n_models, function(i) {
    if (is.null(model_stan$package)) {
      file.path(
        model_stan$model_folder
      )
    } else {
      system.file(
        model_stan$model_folder,
        package = model_stan$package
      )
    }
  }) # identical

  if (!model_stan$profile) {
    stan_no_profiling_list <- lapply(1:n_models, function(i) {
      write_stan_files_no_profiling(
        model_path_list[[i]],
        include_paths_list[[i]]
      )
    })
    model_path_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$model
    })
    include_paths_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$include_paths
    })
  }

  cpp_options <- list()
  if (model_stan$threads) {
    cpp_options[["stan_threads"]] <- TRUE
  }

  stanmodel_list <- lapply(1:n_models, function(i) {
    (cmdstanr::cmdstan_model(model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      dir = dirname(model_path_list[[i]]),
      cpp_options = cpp_options,
      force_recompile = force_recompile
    ))
  })

  model_stan[["load_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_stan)
}

#' Compile BayesianBombCarbon STAN models
#'
#' @description The stan models used by BayesianBombCarbon need to be compiled for your device. This is only necessary once, after installing or updating the package. This function compiles all models from the package.
#'
#' @param model A character vector with a specific model to compile. If NULL (default), all models are compiled.
#' @param force_recompile If TRUE, then models will be recompiled even if they have already been successfully compiled.
#' @param verbose If TRUE, warnings and detailed errors from the compilation are  printed. This can help to diagnose compilation issues.
#'
#' @details If one or several models are not successfully compiled, please ensure that `cmdstan` is properly set up and try updating it to a newer version using [cmdstanr::install_cmdstan()]. If the problem persists, please run [bombcarbon_compile(verbose = TRUE)] and post the output in a new issue on GitHub, along with your [cmdstanr::cmdstan_version()].
#'
#' @export
bombcarbon_compile <- function(model = NULL, force_recompile = FALSE, verbose = FALSE) {
  all_models <- c("delta14C_spline.stan",
                  "delta14C_bspline_shift_wt.stan")
  if (!is.null(model)) {
    if (model %in% c("ref-only", "delta14C_spline")) {
      model <- "delta14C_spline.stan"
    } else if (model %in% c("integrated", "delta14C_bspline_shift")) {
      model <- "delta14C_bspline_shift_wt.stan"
    }
    if (!(model %in% all_models)) {
      cli::cli_abort(paste(
          "The model", model, "is not available."
      ))
    }
    all_models <- model
  }
  comp_success <- NULL

  for (i in 1:length(all_models)) {
    model_name <- all_models[i]
    cat("\r                                                      \r", sep = "");
    cat(sprintf("\r| Compiling model %d/%d ", i, length(all_models)), sep = "")
    #flush.console()
    success <- tryCatch(
      {
        if (verbose == FALSE) {
          suppressMessages(get_stan_model(
            model_filename = model_name,
            force_recompile = force_recompile
          ))
        } else {
          get_stan_model(
            model_filename = model_name,
            force_recompile = force_recompile
          )
        }
        TRUE
      },
      error = function(e) { FALSE }
    )
    cat(paste(sprintf(
      "\r| Compiling model %d/%d",
      i, length(all_models)), ifelse(success, "(success) ", "(failed) ")
      ), sep = "")
    Sys.sleep(0.5)
    comp_success <- c(comp_success, success)
  }
  cat("\r                                                        \r", sep = "");
  if(!all(comp_success)) {
    cli::cli_warn(
      paste(
        "The following models could not be compiled:",
        paste(all_models[!comp_success], collapse = ", ")
      )
    )
  } else {

    cli::cli_alert("All models compiled successfully.")
  }
  return(invisible(NULL))
}