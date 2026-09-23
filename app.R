# =============================================================================
# gasampling Shiny app
#
# Workflow (three steps, one tab each):
#   1. Upload data          -> read_survey_data()
#   2. Build sampling frame -> prepare_sampling_frame()
#   3. Run GA allocation    -> allocate_twostage_ga_joint()
#
# Sources exactly these files, in this order (dependencies before the
# functions that use them), from an R/ subfolder next to app.R if one
# exists, otherwise from app.R's own folder:
#   variance_twostage.R, cost_twostage.R, check_constraints.R,
#   repair_joint.R, fitness_ch_kh.R, feasible_joint_start.R,
#   allocate_twostage_joint.R, read_survey_data.R, sampling_frame.R,
#   frame_c.R
#
# Only top-level `name <- function(...) {...}` definitions are pulled out
# of each file -- any other top-level code (library() calls, loops,
# cat()/print() statements, a pipeline/driver script) is skipped. This
# makes it safe to keep these files in a working directory that also has
# other project scripts (tier-runners, R2BEAT comparisons, etc.) without
# those running every time the app starts.
# =============================================================================

library(shiny)
library(DT)
library(GA)
library(ggplot2)

REQUIRED_FILES <- c(
  "variance_twostage.R",
  "cost_twostage.R",
  "check_constraints.R",
  "repair_joint.R",
  "fitness_ch_kh.R",
  "feasible_joint_start.R",
  "allocate_twostage_joint.R",
  "read_survey_data.R",
  "sampling_frame.R",
  "frame_c.R"
)

# app_dir: folder app.R itself lives in, so this works regardless of the
# working directory runApp() was called from.
app_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
if (is.null(app_dir) || !dir.exists(app_dir)) app_dir <- getwd()

lookup_dir <- if (dir.exists(file.path(app_dir, "R"))) file.path(app_dir, "R") else app_dir

# Only evaluate top-level expressions of the form `name <- function(...) ...`
# (or `name = function(...) ...` / `assign("name", function(...) ...)`).
# Everything else in the file is left untouched.
.load_functions_only <- function(path, envir = globalenv()) {
  exprs <- tryCatch(parse(path), error = function(e) NULL)
  if (is.null(exprs)) {
    warning("Could not parse ", path, " -- skipped.", call. = FALSE)
    return(character(0))
  }
  loaded <- character(0)
  for (e in exprs) {
    is_assign <- is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=", "assign")
    if (!is_assign) next
    rhs <- e[[3]]
    if (!(is.call(rhs) && identical(rhs[[1]], as.name("function")))) next
    fn_name <- as.character(e[[2]])
    eval(e, envir = envir)
    loaded <- c(loaded, fn_name)
  }
  loaded
}

missing_files <- REQUIRED_FILES[!file.exists(file.path(lookup_dir, REQUIRED_FILES))]
loaded_by_file <- list()

if (length(missing_files) > 0) {
  missing_fn_msg <- paste0(
    "Missing file(s) in ", lookup_dir, ": ", paste(missing_files, collapse = ", "),
    ". Place them next to app.R (or in an R/ subfolder) and restart the app."
  )
} else {
  missing_fn_msg <- NULL
  for (fname in REQUIRED_FILES) {
    loaded_by_file[[fname]] <- .load_functions_only(file.path(lookup_dir, fname), envir = globalenv())
  }
  message("gasampling app: sourced all required files --")
  for (fname in names(loaded_by_file)) {
    message("  ", fname, ": ", paste(loaded_by_file[[fname]], collapse = ", "))
  }
}

# Sanity-check the three functions the app calls directly.
entry_point_fns <- c("read_survey_data", "prepare_sampling_frame", "allocate_twostage_ga_joint")
missing_fns <- entry_point_fns[!vapply(entry_point_fns, exists, logical(1))]
if (is.null(missing_fn_msg) && length(missing_fns) > 0) {
  missing_fn_msg <- paste0(
    "All required files were found, but these functions still aren't defined: ",
    paste(missing_fns, collapse = ", "),
    ". Check the function names inside those files match what the app expects."
  )
}
options(shiny.maxRequestSize = 100 * 1024^2)

# =============================================================================
# UI
# =============================================================================

ui <- navbarPage(
  title = "gasampling: Two-Stage Stratified Cluster Sampling",
  id = "main_nav",
  
  # ---------------------------------------------------------------------
  tabPanel(
    "1. Upload Data",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        if (!is.null(missing_fn_msg)) div(class = "alert alert-danger", missing_fn_msg),
        fileInput("file_upload", "Upload SSU-level CSV", accept = ".csv"),
        uiOutput("column_pickers"),
        radioButtons(
          "na_action", "Missing y_var handling",
          choices = c("warn", "stop", "keep"), selected = "warn", inline = TRUE
        ),
        actionButton("btn_read", "Read & standardise data", class = "btn-primary"),
        hr(),
        helpText("Output feeds directly into Step 2.")
      ),
      mainPanel(
        width = 8,
        h4("Raw data preview"),
        DTOutput("raw_preview"),
        hr(),
        h4("read_survey_data() summary"),
        verbatimTextOutput("read_summary"),
        h4("Standardised data preview"),
        DTOutput("clean_preview")
      )
    )
  ),
  
  # ---------------------------------------------------------------------
  tabPanel(
    "2. Sampling Frame",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        numericInput("total_c", "total_c (clusters to allocate)", value = 20, min = 1, step = 1),
        numericInput("cost_cluster", "cost_cluster (per stratum, or single value)", value = 1, min = 0),
        numericInput("cost_unit", "cost_unit (per stratum, or single value)", value = 1, min = 0),
        actionButton("btn_frame", "Build sampling frame", class = "btn-primary"),
        hr(),
        helpText("Requires Step 1 to be completed first. Output feeds directly into Step 3.")
      ),
      mainPanel(
        width = 8,
        h4("prepare_sampling_frame() summary"),
        verbatimTextOutput("frame_summary"),
        h4("Strata table"),
        DTOutput("strata_table"),
        h4("Suggested initial c_h"),
        verbatimTextOutput("ch_suggested")
      )
    )
  ),
  
  # ---------------------------------------------------------------------
  tabPanel(
    "3. Allocation",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        numericInput("total_n", "total_n (target sample size)", value = 500, min = 1, step = 1),
        numericInput("max_cost", "max_cost", value = 10000, min = 0),
        numericInput("popSize", "GA popSize", value = 100, min = 10, step = 10),
        numericInput("maxiter", "GA maxiter", value = 1000, min = 10, step = 10),
        numericInput("run", "GA run (early stop)", value = 50, min = 5, step = 5),
        numericInput("seed", "seed (blank = random)", value = NA),
        actionButton("btn_ga", "Run GA allocation", class = "btn-primary"),
        hr(),
        helpText("Requires Step 2 to be completed first. May take a while for large popSize/maxiter.")
      ),
      mainPanel(
        width = 8,
        uiOutput("feasibility_badge"),
        fluidRow(
          column(4, wellPanel(h5("Total variance"), textOutput("out_variance"))),
          column(4, wellPanel(h5("Total cost"), textOutput("out_cost"))),
          column(4, wellPanel(h5("Total n"), textOutput("out_n")))
        ),
        h4("Allocation by stratum"),
        DTOutput("allocation_table"),
        h4("c_h and k_h by stratum"),
        plotOutput("allocation_plot", height = "350px")
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(raw = NULL, clean = NULL, frame = NULL, alloc = NULL)
  
  # --- Step 1: upload & read ------------------------------------------------
  
  observeEvent(input$file_upload, {
    rv$raw <- tryCatch(
      read.csv(input$file_upload$datapath, stringsAsFactors = FALSE),
      error = function(e) {
        showNotification(paste("Could not read CSV:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    rv$clean <- NULL
    rv$frame <- NULL
    rv$alloc <- NULL
  })
  
  output$column_pickers <- renderUI({
    req(rv$raw)
    cols <- names(rv$raw)
    tagList(
      selectInput("stratum_id", "Stratum column", choices = cols),
      selectInput("cluster_id", "Cluster column", choices = cols),
      selectInput("y_var", "Analysis variable (y_var)", choices = cols)
    )
  })
  
  output$raw_preview <- renderDT({
    req(rv$raw)
    datatable(head(rv$raw, 20), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  observeEvent(input$btn_read, {
    req(rv$raw, input$stratum_id, input$cluster_id, input$y_var)
    
    summary_txt <- capture.output({
      result <- tryCatch(
        withCallingHandlers(
          read_survey_data(
            data       = rv$raw,
            stratum_id = input$stratum_id,
            cluster_id = input$cluster_id,
            y_var      = input$y_var,
            na_action  = input$na_action
          ),
          warning = function(w) {
            showNotification(conditionMessage(w), type = "warning", duration = 8)
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          showNotification(paste("read_survey_data() failed:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      rv$clean <- result
    })
    
    output$read_summary <- renderText(paste(summary_txt, collapse = "\n"))
    
    if (!is.null(rv$clean)) {
      showNotification("Data standardised. Continue to Step 2.", type = "message")
      updateNavbarPage(session, "main_nav", selected = "2. Sampling Frame")
    }
  })
  
  output$clean_preview <- renderDT({
    req(rv$clean)
    datatable(head(rv$clean, 20), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  # --- Step 2: sampling frame -------------------------------------------
  
  observeEvent(input$btn_frame, {
    req(rv$clean)
    y_var_name <- setdiff(names(rv$clean), c("strata", "cluster"))[1]
    
    summary_txt <- capture.output({
      result <- tryCatch(
        prepare_sampling_frame(
          data         = rv$clean,
          y_var        = y_var_name,
          total_c      = input$total_c,
          cost_cluster = input$cost_cluster,
          cost_unit    = input$cost_unit
        ),
        error = function(e) {
          showNotification(paste("prepare_sampling_frame() failed:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      rv$frame <- result
    })
    
    output$frame_summary <- renderText(paste(summary_txt, collapse = "\n"))
    
    if (!is.null(rv$frame)) {
      showNotification("Sampling frame built. Continue to Step 3.", type = "message")
      updateNavbarPage(session, "main_nav", selected = "3. Allocation")
    }
  })
  
  output$strata_table <- renderDT({
    req(rv$frame)
    datatable(rv$frame$strata, options = list(scrollX = TRUE, pageLength = 10)) |>
      formatRound(columns = c("W_h", "sigma2_bh", "sigma2_wh"), digits = 4)
  })
  
  output$ch_suggested <- renderText({
    req(rv$frame)
    paste0(
      "c_h = [", paste(rv$frame$c_h, collapse = ", "), "]",
      "  (total = ", sum(rv$frame$c_h), ")"
    )
  })
  
  # --- Step 3: GA allocation ----------------------------------------------
  
  observeEvent(input$btn_ga, {
    req(rv$frame)
    if (!is.null(missing_fn_msg)) {
      showNotification(missing_fn_msg, type = "error", duration = NULL)
      return(invisible(NULL))
    }
    
    seed_val <- if (is.na(input$seed)) NULL else input$seed
    
    withProgress(message = "Running GA allocation...", value = 0.3, {
      result <- tryCatch(
        allocate_twostage_ga_joint(
          strata   = rv$frame$strata,
          total_n  = input$total_n,
          max_cost = input$max_cost,
          popSize  = input$popSize,
          maxiter  = input$maxiter,
          run      = input$run,
          seed     = seed_val,
          suggestion = NULL,
          monitor  = FALSE
        ),
        error = function(e) {
          showNotification(paste("allocate_twostage_ga_joint() failed:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      rv$alloc <- result
      incProgress(0.7)
    })
    
    if (!is.null(rv$alloc)) {
      showNotification("GA allocation complete.", type = "message")
    }
  })
  
  output$feasibility_badge <- renderUI({
    req(rv$alloc)
    if (isTRUE(rv$alloc$feasible)) {
      div(class = "alert alert-success", strong("Feasible solution found."))
    } else {
      div(class = "alert alert-warning", strong("Solution is NOT fully feasible -- check constraints below."))
    }
  })
  
  output$out_variance <- renderText({
    req(rv$alloc)
    format(rv$alloc$variance, digits = 6)
  })
  
  output$out_cost <- renderText({
    req(rv$alloc)
    scales::label_number(prefix = "R", big.mark = ",")(rv$alloc$total_cost)
  })
  
  output$out_n <- renderText({
    req(rv$alloc)
    format(rv$alloc$total_n, big.mark = ",")
  })
  
  output$allocation_table <- renderDT({
    req(rv$alloc)
    datatable(rv$alloc$allocation, options = list(scrollX = TRUE, pageLength = 10)) |>
      formatRound(columns = c("var_between", "var_within", "var_total"), digits = 6)
  })
  
  output$allocation_plot <- renderPlot({
    req(rv$alloc)
    df <- rv$alloc$allocation
    df_long <- rbind(
      data.frame(stratum = df$stratum, metric = "c_h (clusters)", value = df$c_h),
      data.frame(stratum = df$stratum, metric = "k_h (units per cluster)", value = df$k_h)
    )
    ggplot(df_long, aes(x = stratum, y = value, fill = metric)) +
      geom_col(position = "dodge") +
      labs(x = "Stratum", y = "Value", fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  })
}

shinyApp(ui, server)