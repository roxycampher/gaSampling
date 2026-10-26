
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


EXAMPLE_FILE <- "ecd_subset.csv"
example_hits <- c(
  # installed package (works when the app ships inside gaSampling, e.g. inst/app)
  system.file("extdata", EXAMPLE_FILE, package = "gaSampling"),
  # running from the package source folder
  file.path(app_dir, "inst", "extdata", EXAMPLE_FILE),
  file.path(app_dir, "..", "inst", "extdata", EXAMPLE_FILE)
)
example_hits <- example_hits[file.exists(example_hits)]
example_path <- if (length(example_hits) > 0) example_hits[1] else NA_character_


EXAMPLE <- list(
  stratum = "province", cluster = "id_ward", y_var = "count_children_present",
  total_c = 50, cost_cluster = 1000, cost_unit = 100,
  total_n = 500, max_cost = 1000000,
  popSize = 50, maxiter = 100, run = 15, seed = 2026
)
PAPER_GA <- list(popSize = 200, maxiter = 500, run = 25, seed = 123)

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
entry_point_fns <- c("read_survey_data", "prepare_sampling_frame",
                     "allocate_twostage_ga_joint", "find_feasible_joint_start")
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

# Re-typeset maths after DT tables (re)draw, e.g. on sorting
typeset_js <- DT::JS("function(){ if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]); }")

app_theme_head <- tags$head(
  tags$title("gaSampling"),
  tags$script(
    src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.9/MathJax.js?config=TeX-AMS-MML_HTMLorMML",
    type = "text/javascript"
  ),
  tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function(){
      var typeset = function(){
        if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
      };
      if (window.jQuery) {
        jQuery(document).on('shiny:visualchange shown.bs.tab', function(){ setTimeout(typeset, 60); });
      }
    });
  ")),
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

    table.dataTable thead th{ font-weight:600; color:var(--ink); }
    table.dataTable tfoot th{
      font-weight:600; text-align:right;
      border-top:2px solid var(--panel-border);
    }
    table.dataTable tfoot th:first-child{ text-align:left; }

    .or-divider{
      display:flex; align-items:center; gap:10px;
      margin:14px 0 12px; color:var(--ink-soft); font-size:0.85rem;
    }
    .or-divider::before, .or-divider::after{
      content:''; flex:1; border-top:1px solid var(--panel-border);
    }

    .hint-box{
      background:var(--accent-soft); border:1px solid var(--panel-border);
      border-radius:8px; padding:12px 14px; margin-top:14px;
      font-size:0.85rem; color:var(--ink);
    }
    .hint-box p{ margin:0 0 6px; }
    .hint-box p:last-child{ margin-bottom:0; }
    .hint-box .btn{ margin-top:6px; }
    .hint-bad{ background:#F3E1DE; border-color:#D9B4AE; }

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
              actionButton("btn_example", "Load example", class = "btn-default btn-sm"),
              helpText("ECDC 2021 data -> Click standardise next "),
              div(class = "or-divider", span("or")),
              fileInput("file_upload", "Upload your own dataset (SSU-level CSV)", accept = ".csv"),
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
              div(
                class = "hint-box",
                p(strong("Choosing total_c")),
                p("Keep it at or below the clusters you can afford:",
                  "(max_cost - cost_unit x total_n) / cost_cluster."),
                p("Paper case study (ECDC): 360 clusters at R2,500 per ward and R100 per centre."),
                actionButton("btn_example_frame", "Use example values", class = "btn-default btn-sm")
              ),
              br(),
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
                tabPanel("Suggested initial \\(c_h\\)", uiOutput("ch_suggested"))
              )
            )
          )
        ),
        
        # ---------------------------------------------------------------
        tabPanelBody(
          "step3",
          h3("Allocation", class = "stage-title"),
          p("A feasible Neyman starting point is built first, then the GA refines it.", class = "stage-help"),
          sidebarLayout(
            sidebarPanel(
              width = 4,
              numericInput("total_n", "total_n (target sample size)", value = 5000, min = 1, step = 1),
              numericInput("max_cost", "max_cost", value = 1400000, min = 0),
              checkboxInput("use_neyman", "Start the GA from a feasible Neyman suggestion (faster)", value = TRUE),
              numericInput("popSize", "GA popSize", value = 100, min = 10, step = 10),
              numericInput("maxiter", "GA maxiter", value = 500, min = 10, step = 10),
              numericInput("run", "GA run (early stop)", value = 30, min = 5, step = 5),
              numericInput("seed", "seed (blank = random)", value = NA),
              uiOutput("feasibility_hint"),
              div(
                class = "hint-box",
                p(strong("Getting a feasible run")),
                p("The budget must cover total_n x cost_unit first; what is left buys clusters.",
                  "total_n also cannot exceed the SSUs in the frame."),
                p("Paper case study: n = 5,000 and budget R1,400,000.",
                  "The quick settings (popSize 50, maxiter 100, run 15) finish faster than the",
                  "paper's (200, 500, 25) and may land on a slightly different allocation."),
                actionButton("btn_example_ga", "Use quick example values", class = "btn-default btn-sm"),
                actionButton("btn_example_paper", "Use paper settings (slower)", class = "btn-default btn-sm")
              ),
              br(),
              actionButton("btn_ga", "Run allocation", class = "btn-primary"),
              hr(),
              helpText("Requires Step 2 to be completed first. May take a while for large popSize/maxiter.")
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                id = "step3_tabs",
                tabPanel(
                  "Allocation",
                  uiOutput("feasibility_badge"),
                  fluidRow(
                    column(4, div(class = "metric-box", h5("Total variance"), div(class = "metric-value", textOutput("out_variance")))),
                    column(4, div(class = "metric-box", h5("Total cost"), div(class = "metric-value", textOutput("out_cost")))),
                    column(4, div(class = "metric-box", h5("Total n"), div(class = "metric-value", textOutput("out_n"))))
                  ),
                  h4("Allocation by stratum"),
                  DTOutput("allocation_table"),
                  h4("\\(c_h\\) and \\(k_h\\) by stratum"),
                  plotOutput("allocation_plot", height = "350px")
                ),
                tabPanel(
                  "Starting Suggestion",
                  uiOutput("start_summary"),
                  DTOutput("start_table")
                )
              )
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
  
  rv <- reactiveValues(raw = NULL, clean = NULL, frame = NULL, alloc = NULL,
                       map = NULL, costs = NULL, start = NULL)
  
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
    rv$start <- NULL
    rv$map <- NULL
  })
  
  # --- Example data: ECDC 2021, with the column mapping pre-filled -----------
  
  observeEvent(input$btn_example, {
    if (is.na(example_path)) {
      showNotification(
        paste0("Example file '", EXAMPLE_FILE, "' not found next to app.R (or in a data/ folder)."),
        type = "error", duration = NULL
      )
      return(invisible(NULL))
    }
    dat <- tryCatch(
      read.csv(example_path, stringsAsFactors = FALSE),
      error = function(e) {
        showNotification(paste("Could not read example:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    req(dat)
    # count_children_present contains codes such as '.a'; parse to numeric (NA where not a number)
    if (EXAMPLE$y_var %in% names(dat)) {
      raw_y <- as.character(dat[[EXAMPLE$y_var]])
      dat[[EXAMPLE$y_var]] <- if (requireNamespace("readr", quietly = TRUE)) {
        suppressWarnings(readr::parse_number(raw_y))
      } else {
        suppressWarnings(as.numeric(raw_y))
      }
    }
    rv$raw <- dat
    rv$clean <- NULL
    rv$frame <- NULL
    rv$alloc <- NULL
    rv$start <- NULL
    rv$map <- list(stratum = EXAMPLE$stratum, cluster = EXAMPLE$cluster, y_var = EXAMPLE$y_var)
    showNotification("ECDC example loaded. Click 'Read & standardise data' to continue.", type = "message")
    updateTabsetPanel(session, "step1_tabs", selected = "Raw Data Preview")
  })
  
  observeEvent(input$btn_example_frame, {
    updateNumericInput(session, "total_c", value = EXAMPLE$total_c)
    updateNumericInput(session, "cost_cluster", value = EXAMPLE$cost_cluster)
    updateNumericInput(session, "cost_unit", value = EXAMPLE$cost_unit)
  })
  
  observeEvent(input$btn_example_ga, {
    updateNumericInput(session, "total_n", value = EXAMPLE$total_n)
    updateNumericInput(session, "max_cost", value = EXAMPLE$max_cost)
    updateNumericInput(session, "popSize", value = EXAMPLE$popSize)
    updateNumericInput(session, "maxiter", value = EXAMPLE$maxiter)
    updateNumericInput(session, "run", value = EXAMPLE$run)
    updateNumericInput(session, "seed", value = EXAMPLE$seed)
    updateCheckboxInput(session, "use_neyman", value = TRUE)
  })
  
  observeEvent(input$btn_example_paper, {
    updateNumericInput(session, "total_n", value = EXAMPLE$total_n)
    updateNumericInput(session, "max_cost", value = EXAMPLE$max_cost)
    updateNumericInput(session, "popSize", value = PAPER_GA$popSize)
    updateNumericInput(session, "maxiter", value = PAPER_GA$maxiter)
    updateNumericInput(session, "run", value = PAPER_GA$run)
    updateNumericInput(session, "seed", value = PAPER_GA$seed)
    updateCheckboxInput(session, "use_neyman", value = TRUE)
  })
  
  output$column_pickers <- renderUI({
    req(rv$raw)
    cols <- names(rv$raw)
    pick <- function(key) {
      if (!is.null(rv$map) && rv$map[[key]] %in% cols) rv$map[[key]] else cols[1]
    }
    tagList(
      selectInput("stratum_id", "Stratum column", choices = cols, selected = pick("stratum")),
      selectInput("cluster_id", "Cluster column", choices = cols, selected = pick("cluster")),
      selectInput("y_var", "Analysis variable (y_var)", choices = cols, selected = pick("y_var"))
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
    rv$costs <- list(cluster = input$cost_cluster, unit = input$cost_unit)
    rv$alloc <- NULL
    rv$start <- NULL
    
    output$frame_summary <- renderText(paste(summary_txt, collapse = "\n"))
    
    if (!is.null(rv$frame)) {
      showNotification("Sampling frame built. Continue to Step 3.", type = "message")
      updateTabsetPanel(session, "step2_tabs", selected = "Strata Table")
    }
  })
  
  output$strata_table <- renderDT({
    req(rv$frame)
    st <- rv$frame$strata
    math_labels <- c(
      C_h = "\\(C_h\\)", N_h = "\\(N_h\\)", max_k_h = "\\(k_h\\)",
      W_h = "\\(W_h\\)", sigma2_bh = "\\(\\sigma^2_{bh}\\)", sigma2_wh = "\\(\\sigma^2_{wh}\\)"
    )
    hit <- intersect(names(math_labels), names(st))
    dt <- datatable(
      st,
      rownames = FALSE,
      colnames = if (length(hit) > 0) stats::setNames(hit, unname(math_labels[hit])) else names(st),
      options = list(scrollX = TRUE, pageLength = 10, drawCallback = typeset_js)
    )
    round_idx <- which(names(st) %in% c("W_h", "sigma2_bh", "sigma2_wh"))
    if (length(round_idx) > 0) dt <- formatRound(dt, columns = round_idx, digits = 4)
    dt
  })
  
  output$ch_suggested <- renderUI({
    req(rv$frame)
    ch <- rv$frame$c_h
    div(
      style = "overflow-x:auto;",
      p(sprintf("\\[c_h = (%s)\\]", paste(ch, collapse = ",\\ "))),
      p(sprintf("Total clusters: %s", format(sum(ch), big.mark = ",")))
    )
  })
  
  # --- Step 3: GA allocation ----------------------------------------------
  
  observeEvent(input$btn_ga, {
    req(rv$frame, rv$costs)
    if (!is.null(missing_fn_msg)) {
      showNotification(missing_fn_msg, type = "error", duration = NULL)
      return(invisible(NULL))
    }
    
    strata   <- rv$frame$strata
    cc       <- rv$costs$cluster
    cu       <- rv$costs$unit
    seed_val <- if (is.na(input$seed)) NULL else input$seed
    
    if (cu * input$total_n >= input$max_cost) {
      showNotification(
        "Unit costs alone use up the whole budget. Raise max_cost or lower total_n.",
        type = "error", duration = NULL
      )
      return(invisible(NULL))
    }
    
    # Clusters we can both afford and actually have (same rule as the R2BEAT tiers)
    available_c  <- sum(strata$C_h)
    affordable_c <- floor((input$max_cost - cu * input$total_n) / cc)
    max_c        <- min(available_c, max(1, affordable_c))
    
    withProgress(message = "Running allocation", value = 0, {
      
      # Step A: feasible Neyman starting suggestion, so the GA starts close to a solution
      ga_suggestion <- NULL
      rv$start <- NULL
      if (isTRUE(input$use_neyman)) {
        setProgress(0.1, detail = "Building feasible starting suggestion...")
        ns <- tryCatch(
          find_feasible_joint_start(
            strata                      = strata,
            total_n                     = input$total_n,
            max_cost                    = input$max_cost,
            max_candidate_c             = max_c,
            default_starting_suggestion = "3. Neyman proportional",
            require_all_starts          = FALSE
          ),
          error = function(e) {
            showNotification(
              paste0("Starting suggestion failed (", conditionMessage(e),
                     "). The GA will build its own start instead."),
              type = "warning", duration = 12
            )
            NULL
          }
        )
        if (!is.null(ns)) {
          ga_suggestion <- ns$suggestion
          rv$start <- ns
        }
      }
      
      # Step B: the GA, seeded with that suggestion
      setProgress(0.3, detail = "Running genetic algorithm...")
      result <- tryCatch(
        allocate_twostage_ga_joint(
          strata     = strata,
          total_n    = input$total_n,
          max_cost   = input$max_cost,
          popSize    = input$popSize,
          maxiter    = input$maxiter,
          run        = input$run,
          seed       = seed_val,
          suggestion = ga_suggestion,
          monitor    = FALSE
        ),
        error = function(e) {
          showNotification(paste("allocate_twostage_ga_joint() failed:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      rv$alloc <- result
      setProgress(1)
    })
    
    if (!is.null(rv$alloc)) {
      showNotification("Allocation complete.", type = "message")
      updateTabsetPanel(session, "step3_tabs", selected = "Allocation")
    }
  })
  
  # Live check so people know before running whether the inputs can work
  output$feasibility_hint <- renderUI({
    req(rv$frame, rv$costs, input$total_n, input$max_cost)
    strata <- rv$frame$strata
    cc <- rv$costs$cluster
    cu <- rv$costs$unit
    money <- scales::label_number(prefix = "R", big.mark = ",")
    
    unit_total <- cu * input$total_n
    if (unit_total >= input$max_cost) {
      return(div(class = "hint-box hint-bad", p(sprintf(
        "Unit costs alone (%s) use up the whole budget. Raise max_cost or lower total_n.",
        money(unit_total)
      ))))
    }
    if (input$total_n > sum(strata$N_h)) {
      return(div(class = "hint-box hint-bad", p(sprintf(
        "total_n (%s) is larger than the %s SSUs in the frame.",
        format(input$total_n, big.mark = ","), format(sum(strata$N_h), big.mark = ",")
      ))))
    }
    affordable <- floor((input$max_cost - unit_total) / cc)
    if (affordable < nrow(strata)) {
      return(div(class = "hint-box hint-bad", p(sprintf(
        "After unit costs the budget buys only %d cluster(s), fewer than the %d strata. Raise max_cost.",
        affordable, nrow(strata)
      ))))
    }
    div(class = "hint-box", p(sprintf(
      "Looks workable: unit costs take %s, leaving room for up to %s clusters (%s available).",
      money(unit_total), format(affordable, big.mark = ","), format(sum(strata$C_h), big.mark = ",")
    )))
  })
  
  output$start_summary <- renderUI({
    if (is.null(rv$start)) {
      return(helpText("No starting suggestion yet. Run the allocation with the Neyman start option ticked."))
    }
    ns <- rv$start
    div(
      class = "hint-box",
      p(strong("Selected start: "), as.character(ns$starting_suggestion %||% "n/a")),
      p(strong("Starting clusters: "), as.character(ns$total_c %||% "n/a"))
    )
  })
  
  output$start_table <- renderDT({
    req(rv$start, rv$start$diagnostics)
    datatable(rv$start$diagnostics, options = list(scrollX = TRUE, pageLength = 10))
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
  
  # The table is the `allocation` data frame exactly as returned by
  # allocate_twostage_ga_joint(); only the display is formatted. Any extra
  # columns the function adds later show up automatically.
  output$allocation_table <- renderDT({
    req(rv$alloc)
    al <- rv$alloc$allocation
    
    labels <- c(
      stratum = "Stratum", c_h = "\\(c_h\\) (clusters)", k_h = "\\(k_h\\) (units per cluster)",
      n_h = "\\(n_h\\) (total units)", var_between = "Between-cluster var",
      var_within = "Within-cluster var", var_total = "Total var"
    )
    header_names <- ifelse(names(al) %in% names(labels), labels[names(al)], names(al))
    
    is_num   <- vapply(al, is.numeric, logical(1))
    is_whole <- vapply(al, function(x) is.numeric(x) && all(abs(x - round(x)) < 1e-9, na.rm = TRUE), logical(1))
    var_cols   <- names(al)[is_num & grepl("^var_", names(al))]
    count_cols <- names(al)[is_whole & !(names(al) %in% var_cols)]
    
    # Footer totals (sums only for additive columns; k_h is left blank)
    foot_cell <- function(nm, i) {
      if (i == 1) return(tags$th("Total"))
      if (nm %in% var_cols) return(tags$th(formatC(sum(al[[nm]]), format = "f", digits = 6, big.mark = ",")))
      if (nm %in% count_cols && nm != "k_h") return(tags$th(formatC(sum(al[[nm]]), format = "f", digits = 0, big.mark = ",")))
      tags$th("")
    }
    container <- tags$table(
      class = "display",
      DT::tableHeader(header_names),
      tags$tfoot(tags$tr(lapply(seq_along(al), function(i) foot_cell(names(al)[i], i))))
    )
    
    dt <- datatable(
      al,
      container = container,
      rownames = FALSE,
      class = "row-border hover",
      options = list(dom = "t", paging = FALSE, scrollX = TRUE, scrollY = "420px", scrollCollapse = TRUE,
                     drawCallback = typeset_js)
    )
    if (length(count_cols) > 0) dt <- formatRound(dt, columns = count_cols, digits = 0, mark = ",")
    if (length(var_cols) > 0)   dt <- formatRound(dt, columns = var_cols, digits = 6)
    dt
  })
  
  output$allocation_plot <- renderPlot({
    req(rv$alloc)
    df <- rv$alloc$allocation
    df_long <- rbind(
      data.frame(stratum = df$stratum, metric = "ch", value = df$c_h),
      data.frame(stratum = df$stratum, metric = "kh", value = df$k_h)
    )
    df_long$metric <- factor(df_long$metric, levels = c("ch", "kh"))
    ggplot(df_long, aes(x = stratum, y = value, fill = metric)) +
      geom_col(position = "dodge") +
      scale_fill_discrete(
        labels = expression(c[h] ~ "(clusters)", k[h] ~ "(units per cluster)")
      ) +
      labs(x = "Stratum", y = "Value", fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  })
}

shinyApp(ui, server)