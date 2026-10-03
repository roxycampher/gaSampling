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

# Shiny's default upload limit is 5 MB, too small for full census-scale
# SSU CSVs -- raise it here (adjust the multiplier if your files are bigger).
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

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

# =============================================================================
# UI
# =============================================================================

app_theme_head <- tags$head(
  tags$title("gaSampling"),
  tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
  tags$link(
    href = paste0(
      "https://fonts.googleapis.com/css2?",
      "family=Inter:wght@400;500;600&",
      "family=Plus+Jakarta+Sans:wght@600;700&display=swap"
    ),
    rel = "stylesheet"
  ),
  tags$style(HTML("
    :root{
      --bg:#E6E3DC;
      --panel:#FBF9F5;
      --panel-border:#DAD4C7;
      --ink:#26231F;
      --ink-soft:#8A8276;
      --accent:#33544A;
      --accent-soft:#DCE6E1;
      --danger:#9A3B32;
    }
    html, body{
      background:var(--bg) !important;
      color:var(--ink);
      font-family:'Inter', -apple-system, sans-serif;
    }
    .page-wrap{
      max-width:1400px;
      margin:0 auto;
      padding:36px 32px 60px;
    }
    .main-frame{
      background:var(--panel);
      border:1px solid var(--panel-border);
      border-radius:16px;
      box-shadow:0 2px 8px rgba(38,35,31,0.06);
      padding:48px 56px 56px;
    }
    .app-header{ text-align:center; margin-bottom:8px; }
    .wordmark{
      font-family:'Plus Jakarta Sans', sans-serif;
      font-weight:700;
      font-size:2.5rem;
      letter-spacing:-0.02em;
      margin:0;
      color:var(--ink);
    }
    .app-subtitle{ color:var(--ink-soft); margin-top:4px; font-size:1rem; }

    .stepper{
      display:flex; justify-content:center; align-items:center;
      gap:0; margin:32px 0 8px; flex-wrap:wrap;
    }
    .step{
      display:flex; align-items:center; gap:10px;
      background:none; border:none; padding:4px 6px; cursor:pointer;
      opacity:0.5;
    }
    .step.is-locked{ cursor:not-allowed; }
    .step.is-active{ opacity:1; }
    .step.is-done{ opacity:0.9; }
    .step-circle{
      width:30px; height:30px; border-radius:50%; flex:0 0 auto;
      display:flex; align-items:center; justify-content:center;
      border:1.5px solid var(--ink-soft);
      font-weight:600; font-size:0.85rem; color:var(--ink-soft);
      background:transparent;
    }
    .step.is-active .step-circle{ background:var(--accent); border-color:var(--accent); color:#fff; }
    .step.is-done .step-circle{ background:var(--accent-soft); border-color:var(--accent); color:var(--accent); }
    .step-label{ font-size:0.95rem; color:var(--ink-soft); white-space:nowrap; }
    .step.is-active .step-label{ color:var(--ink); font-weight:600; }
    .step-connector{ width:48px; height:1.5px; background:var(--panel-border); margin:0 4px; }

    .frame-divider{ border:none; border-top:1px solid var(--panel-border); margin:28px 0 30px; }

    .stage-title{
      font-family:'Plus Jakarta Sans', sans-serif; font-weight:700; font-size:1.4rem;
      margin:0 0 4px;
    }
    .stage-help{ color:var(--ink-soft); font-size:0.92rem; margin-bottom:22px; }

    /* Blend Bootstrap's sidebar/well chrome into the single off-white frame */
    .main-frame .well{ background:transparent; border:none; box-shadow:none; padding:0; }
    .main-frame .col-sm-4{ border-right:1px solid var(--panel-border); padding-right:26px; }
    .main-frame .col-sm-8{ padding-left:30px; }
    .main-frame h4{ font-family:'Plus Jakarta Sans', sans-serif; font-weight:600; color:var(--ink); margin-top:20px; }
    .main-frame h4:first-child{ margin-top:0; }
    .main-frame label{ font-weight:500; color:var(--ink); }
    .main-frame hr{ border-top:1px solid var(--panel-border); }

    .btn-primary{ background:var(--accent); border-color:var(--accent); }
    .btn-primary:hover, .btn-primary:focus{ background:#284339; border-color:#284339; }

    .metric-box{
      background:var(--accent-soft); border:1px solid var(--panel-border);
      border-radius:8px; padding:14px 16px; text-align:center;
    }
    .metric-box h5{ margin:0 0 4px; color:var(--ink-soft); font-size:0.82rem; font-weight:500; }
    .metric-box .metric-value{ font-family:'Plus Jakarta Sans', sans-serif; font-size:1.4rem; font-weight:700; color:var(--ink); }

    /* Custom Styling for Output Tabs */
    .nav-tabs {
      border-bottom: 1px solid var(--panel-border);
      margin-bottom: 20px;
    }
    .nav-tabs > li > a {
      color: var(--ink-soft);
      font-weight: 500;
      border: none;
      border-bottom: 2px solid transparent;
      background: transparent !important;
    }
    .nav-tabs > li.active > a, 
    .nav-tabs > li.active > a:hover, 
    .nav-tabs > li.active > a:focus {
      color: var(--ink) !important;
      font-weight: 600;
      border: none;
      border-bottom: 2px solid var(--accent);
      background: transparent !important;
    }

    @media (max-width: 767px){
      .main-frame .col-sm-4{ border-right:none; border-bottom:1px solid var(--panel-border); padding-right:0; padding-bottom:20px; margin-bottom:20px; }
      .main-frame .col-sm-8{ padding-left:0; }
    }
  "))
)

ui <- fluidPage(
  title = "gaSampling",
  app_theme_head,
  
  div(
    class = "page-wrap",
    div(
      class = "main-frame",
      
      div(
        class = "app-header",
        h1("gaSampling", class = "wordmark"),
        p("Two-stage stratified cluster sampling", class = "app-subtitle")
      ),
      
      uiOutput("stepper"),
      
      if (!is.null(missing_fn_msg)) div(class = "alert alert-danger", missing_fn_msg),
      
      hr(class = "frame-divider"),
      
      tabsetPanel(
        id = "wizard",
        type = "hidden",
        
        # ---------------------------------------------------------------
        tabPanelBody(
          "step1",
          h3("Upload data", class = "stage-title"),
          p("Input census data (or previous study raw data)", class = "stage-help"),
          sidebarLayout(
            sidebarPanel(
              width = 4,
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
              tabsetPanel(
                id = "step1_tabs",
                tabPanel("Raw Data Preview", DTOutput("raw_preview")),
                tabPanel("Read Summary", verbatimTextOutput("read_summary")),
                tabPanel("Standardised Data Preview", DTOutput("clean_preview"))
              )
            )
          )
        ),
        
        # ---------------------------------------------------------------
        tabPanelBody(
          "step2",
          h3("Sampling frame", class = "stage-title"),
          p("Set cluster and cost targets to build the per-stratum frame.", class = "stage-help"),
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
              tabsetPanel(
                id = "step2_tabs",
                tabPanel("Frame Summary", verbatimTextOutput("frame_summary")),
                tabPanel("Strata Table", DTOutput("strata_table")),
                tabPanel("Suggested Initial c_h", verbatimTextOutput("ch_suggested"))
              )
            )
          )
        ),
        
        # ---------------------------------------------------------------
        tabPanelBody(
          "step3",
          h3("Allocation", class = "stage-title"),
          p("Run the GA and review the resulting allocation.", class = "stage-help"),
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
                column(4, div(class = "metric-box", h5("Total variance"), div(class = "metric-value", textOutput("out_variance")))),
                column(4, div(class = "metric-box", h5("Total cost"), div(class = "metric-value", textOutput("out_cost")))),
                column(4, div(class = "metric-box", h5("Total n"), div(class = "metric-value", textOutput("out_n"))))
              ),
              h4("Allocation by stratum"),
              DTOutput("allocation_table"),
              h4("c_h and k_h by stratum"),
              plotOutput("allocation_plot", height = "350px")
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(raw = NULL, clean = NULL, frame = NULL, alloc = NULL)
  
  # --- Stepper: reflects real progress, and is how you navigate back ----
  
  output$stepper <- renderUI({
    current <- input$wizard %||% "step1"
    done1 <- !is.null(rv$clean)
    done2 <- !is.null(rv$frame)
    
    step_div <- function(value, num, label, is_done, is_locked) {
      cls <- c("step")
      if (identical(current, value)) cls <- c(cls, "is-active")
      if (is_done) cls <- c(cls, "is-done")
      if (is_locked) cls <- c(cls, "is-locked")
      actionLink(
        inputId = paste0("goto_", value),
        label = tagList(
          span(class = "step-circle", if (is_done) "\u2713" else num),
          span(class = "step-label", label)
        ),
        class = paste(cls, collapse = " ")
      )
    }
    
    tagList(
      div(
        class = "stepper",
        step_div("step1", "1", "Upload data", done1, FALSE),
        div(class = "step-connector"),
        step_div("step2", "2", "Sampling frame", done2, !done1),
        div(class = "step-connector"),
        step_div("step3", "3", "Allocation", FALSE, !done2)
      )
    )
  })
  
  observeEvent(input$goto_step1, updateTabsetPanel(session, "wizard", selected = "step1"))
  
  observeEvent(input$goto_step2, {
    if (is.null(rv$clean)) {
      showNotification("Finish Step 1 (upload & read data) first.", type = "warning")
    } else {
      updateTabsetPanel(session, "wizard", selected = "step2")
    }
  })
  
  observeEvent(input$goto_step3, {
    if (is.null(rv$frame)) {
      showNotification("Finish Step 2 (build the sampling frame) first.", type = "warning")
    } else {
      updateTabsetPanel(session, "wizard", selected = "step3")
    }
  })
  
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
      updateTabsetPanel(session, "step1_tabs", selected = "Standardised Data Preview")
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
      updateTabsetPanel(session, "step2_tabs", selected = "Strata Table")
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
      div(class = "alert alert-warning", strong("Solution is NOT fully feasible --> check constraints below."))
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