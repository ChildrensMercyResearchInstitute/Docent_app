# Load libraries ----
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(shinycssloaders)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(mongolite)
  library(glue)
  library(ggplot2)
  library(parallel)
  library(yaml)
  library(RSQLite)
})

log_message <- function(...) {
  message(Sys.time(), " ", paste0(..., collapse=""))
}

# Load configuration -----
log_message("Reading configuration...")

config_file <- "config.yaml"

if(!file.exists(config_file)) {
  stop("Configuration file not found: ", config.file)
}

app_config <- read_yaml(config_file)

datapath <- app_config$datapath

if(!dir.exists(datapath)) {
  stop("Datapath not found: ", datapath)
}

# Utility functions -----

cores <- function() {
  return(app_config$cores)
}

pn <- function(number) {
  prettyNum(number, big.mark=",")
}

# use glue() with <> instead of {}
glue2 <- function(..., .envir=parent.frame()) {
  glue(..., .open="<", .close=">", .envir=.envir)
}

# Individual remapping -----

individual_remapping_table <- readRDS(app_config$individual_remapping)

# Mongo -----

mongo_connections <- list()

mongo_connection <- function(collection) {
  if(is.null(mongo_connections[[collection]])) {
    mongo_connection_url <- glue("mongodb://{user}:{password}@{host}/{db}",
                                 user     = app_config$mongo$user,
                                 password = app_config$mongo$password,
                                 host     = app_config$mongo$host,
                                 db       = app_config$mongo$db)
    
    mongo_connections[[collection]] <<- mongo(url=mongo_connection_url, collection=collection)
  }
  return(mongo_connections[[collection]])
}

# Helpers -----

gene_id <- function(gene_names) {
  genes_df$gene_id[match(gene_names, genes_df$gene)]
}

load_gene_data <- function(gene_name) {
  gid <- gene_id(gene_name)
  gene_rds_file <- file.path(datapath, "expression", gid %% 100, glue("gene_{gid}.rds"))
  if(!file.exists(gene_rds_file)) {
    log_message("Not found: ", gene_rds_file)
    stop()
  }
  
  secs <- system.time(gene_data_df <- readRDS(gene_rds_file))["elapsed"]
  log_message(glue("Loaded {gene_name} gene data file in {secs} seconds", 
                   gene_name=gene_name,
                   secs=round(secs, 1)))
  
  gene_data_df %<>% mutate(gene = gene_name) %>%
    select(gene, barcode_id, exp)
  
  gene_data_df
}

load_cluster_method <- function(method_name) {
  log_message("Loading clustering: ", method_name)
  readRDS(file.path(datapath, "clustering", glue("{method_name}.rds"))) %>%
    mutate(cluster = as.character(cluster))
}

log_message("Determining available clustering methods...")
cluster_options <- list.files(file.path(datapath, "clustering"), full.names=TRUE) %>%
  basename %>%
  gsub("\\.rds$", "", .)

stopifnot(length(cluster_options) > 0)

log_message("Reading genes...")
genes_df <- readRDS(file.path(datapath, "genes.rds"))

log_message("Reading individuals...")
individuals_df <- readRDS(file.path(datapath, "individuals.rds"))

indiv_metadata <- colnames(individuals_df)
indiv_metadata <- indiv_metadata[! indiv_metadata %in% c("individual", "individual_id")]

# Shiny UI -----------

plot_base_text_size <- 15

ui <- fluidPage(
  useShinyjs(),
  fluidRow(
    column(12, titlePanel(glue("Docent: {basename(datapath)}")))
  ),
  fluidRow(
    column(2, 
      tabsetPanel(
      tabPanel("Settings",
        inputPanel(  
          selectInput(
              inputId = "genes",
              label = "Select gene(s):",
              choices = NULL,
              multiple = TRUE
            ),
            selectInput(
              inputId = "clusterMethod",
              label = "Clustering method",
              choices = cluster_options,
              multiple = FALSE
            ),
            selectInput(
              inputId = "cluster",
              label = "Select cluster:",
              choices = c("Cluster"=""),
              multiple = TRUE
            )
        ),
        inputPanel(
          radioButtons(
            inputId = "includeCounts",
            label = "Include counts:",
            choices = list("None", "Individuals", "Cells")
          )
        ),
        inputPanel(
          radioButtons(
            inputId = "groupBy",
            label = "Group by:",
            choices = list("No grouping", "Individual metadata", "Genotype")
          ),
          selectInput(
            inputId = "groupIndividualsBy",
            label = "Metadata:",
            choices = indiv_metadata
          ),
          textInput(
            inputId = "rsID",
            label = "rsID:"
          ),
          textOutput("rsid_validation")
        ),
        inputPanel(
          pickerInput(
            inputId = "individuals", 
            label = "Include individuals:", 
            choices = sort(individuals_df$individual), 
            selected = sort(individuals_df$individual),
            options = list(
              `actions-box` = TRUE, 
              size = 10,
              `selected-text-format` = "count > 3"
            ), 
            multiple = TRUE
          )
        )
      ),
      tabPanel("Summary",
        tabsetPanel(type="pills",
          tabPanel("Cells", withSpinner(plotOutput("summary_plot_cells"))),
          tabPanel("Individuals", withSpinner(plotOutput("summary_plot_individuals")))
        ),
        textOutput("database_info")
      ))
    ),
    column(10, 
      withSpinner(plotOutput("by_cluster_plot", click="clusterPlot_click")),
      tabsetPanel(
        tabPanel("Individuals",
          withSpinner(plotOutput("cluster_individual_plot"))
        ),
        tabPanel("Grouped", 
          fluidRow(
            column(8, withSpinner(plotOutput("cluster_grouped_plot"))),
            column(4,  
              h4("Wilcoxon Rank Sum test"),
              withSpinner(tableOutput("grouped_stat_tests")))
          )
        )
      )
    )
  )
)

# Shiny server -----

server <- function(input, output, session) {
  
  # server-side autocomplete for genes
  updateSelectizeInput(session, inputId='genes', choices=sort(genes_df$gene), server=TRUE)

  observe({
    log_message("Updating clusters")
    updateSelectInput(session=session, inputId="cluster", choices=c("Cluster"="", current_cluster_options()))
  })
  
  observe({
    if(input$groupBy == "No grouping") {
      shinyjs::hide("groupIndividualsBy")
      shinyjs::hide("rsID")
      shinyjs::hide("affected_gene")
    }

    if(input$groupBy == "Individual metadata") {
      shinyjs::show("groupIndividualsBy")
      shinyjs::hide("rsID")
      shinyjs::hide("affected_gene")
    }

    if(input$groupBy == "Genotype") {
      shinyjs::show("rsID")
      shinyjs::hide("groupIndividualsBy")
      shinyjs::hide("affected_gene")
    }
  })
  
  observeEvent(input$clusterPlot_click, {
    clusters_df <- by_cluster_data()
    clusters <- levels(clusters_df$cluster)
    
    x_index <- as.integer(round(input$clusterPlot_click$x))
    log_message(glue("Cluster {cluster_id} toggled", cluster_id=clusters[x_index]))
    if(x_index %in% seq_along(clusters)) {
    
      current_clusters <- unlist(input$cluster)
      if(clusters[x_index] %in% input$cluster) {
        new_clusters <- input$cluster[!input$cluster %in% clusters[x_index]]
      } else {
        new_clusters <- c(input$cluster, clusters[x_index]) %>% sort
      }
      updateSelectInput(
        session,
        inputId="cluster",
        choices=as.list(clusters),
        selected=new_clusters
      )
      log_message(glue("Selected cluster(s): {clusters}", clusters=paste0(new_clusters, collapse=", ")))
    }
  })
  
  ## Reactives ----
  
  ### selected_individual_ids ----
  selected_individual_ids <- reactive({
    req(input$individuals)
    filter(individuals_df, individual %in% input$individuals)$individual_id
  })
  
  current_cluster_options <- reactive({
    cell_data()$cluster %>% unique
  })
  
  cell_data <- reactive({
    req(input$clusterMethod)
    load_cluster_method(input$clusterMethod)
  })
  
  genotypes_by_rsid <- reactive({
    if(input$rsID == "") return(NULL)
    log_message("Updating genotypes...")
    
    # genos <- mongo_connection("genotypes")$find(glue2('{ "rsid": "<input$rsID>" }'))
    #if(nrow(genos) == 0) {
    #  log_message("No genotypes for rsID ", input$rsID)?
    #  return(NULL)
    #}

    #write.table(genos, file = "genos.txt", sep = "\t",row.names = TRUE, col.names = NA)
    #log_message(i)

    conn <- dbConnect(RSQLite::SQLite(), "genotype_v2.db")
    dbListTables(conn)
    select_statement<-paste0("SELECT * FROM genotypes where rsid = '", input$rsID, "'")
    write.csv(select_statement, file = "statement.txt")
    genos<-dbGetQuery(conn, select_statement)
    if(nrow(genos) == 0) {
      log_message("No genotypes for rsID ", input$rsID)?
      return(NULL)
    }
    write.table(genos, file = "genos.txt", sep = "\t",row.names = TRUE, col.names = NA)
    
    genos %<>% select(rsid, individual=i, genotype=geno)
    log_message("Found ", nrow(genos), " genotypes for rsID ", input$rsID)
    
    # fix genotypes here
    
    genos_good <- filter(genos, individual %in% individuals_df$individual)
    genos_bad  <- filter(genos, !individual %in% individuals_df$individual)
    
    if(nrow(genos_bad) > 0) {
      genos_bad$individual <- individual_remapping_table$individual[match(genos_bad$individual, individual_remapping_table$actual)]
      genos_bad %<>% filter(!is.na(individual))
    }
    
    genos <- bind_rows(genos_good, genos_bad)
    
    return(genos)
  })
  
  
  unfiltered_gene_data <- reactive({
    req(input$genes)
    
    gene_data <- input$genes %>%
      mclapply(load_gene_data, mc.cores=cores()) %>%
      bind_rows
    gene_data
  })
  
  current_dataset <- reactive({

    current_cell_data <- cell_data() %>%
      select(barcode_id, cluster, individual_id) %>%
      inner_join(select(individuals_df, individual, individual_id), by="individual_id")
    
    current_individual_ids <- selected_individual_ids()
    
    dataset <- unfiltered_gene_data() %>%
      inner_join(current_cell_data, by="barcode_id") %>%
      filter(individual_id %in% current_individual_ids)
    
    dataset
  })
  
  by_cluster_data <- reactive({
    
    cluster_data <- current_dataset() %>%
      group_by(gene, cluster) %>%
      summarize(mean_exp = mean(exp),
                cells = n(),
                individuals = length(unique(individual_id)))

    cluster_data$cluster %<>% factor(levels=current_cluster_options())
    cluster_data
  })
  
  cluster_individual_data <- reactive({
    req(input$cluster)
    
    cluster_data <- current_dataset() %>%
            filter(cluster %in% unlist(input$cluster)) %>%
            group_by(gene, individual_id) %>%
            summarize(mean_exp = mean(exp),
                      cells = n()) %>%
            inner_join(individuals_df, by="individual_id") %>%
            ungroup

    cluster_data
  })
  
  byClusterPlot <- reactive({

    req(input$genes)
    
    data_df <- by_cluster_data()   
    max_exp <- max(data_df$mean_exp)
    
    if(length(input$cluster) == 0) {
      g <- ggplot(data_df, aes(x=cluster, y=mean_exp)) +
        geom_bar(stat="identity") +
        facet_wrap(~ gene, ncol=1, scales="free_y") +
        scale_y_continuous(limits=c(0, max_exp*1.15), expand=c(0, 0)) +
        labs(x="Cluster", 
             y="Mean expression",
             title="Click on a bar to view individual gene expression within a cluster") +
        theme_bw(plot_base_text_size) +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
      
    } else {
      data_df %<>% mutate(cluster_selected = ifelse(cluster %in% unlist(input$cluster), "selected", "unselected"))
      
      g <- ggplot(data_df, aes(x=cluster, y=mean_exp, fill=cluster_selected)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=c("selected"="blue", "unselected"="gray90"), guide="none") +
        facet_wrap(~ gene, ncol=1, scales="free_y") +
        scale_y_continuous(limits=c(0, max_exp*1.15), expand=c(0, 0)) +
        labs(x="Cluster", 
             y="Mean expression",
             title="Click on a bar to select a cluster or clusters") +
        theme_bw(plot_base_text_size) +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
      
    }  
    
    if(input$includeCounts == "Individuals") {
      g <- g + geom_text(aes(label=individuals, y=Inf), vjust=1.5) 
    }
    
    if(input$includeCounts == "Cells") {
      g <- g + geom_text(aes(label=pn(cells), y=Inf), vjust=1.5) 
    }
    
    g
  })
  
  selected_clusters_label <- reactive({
    req(input$cluster)
    paste0(input$cluster, collapse=", ")
  })
  
  selected_clusters_data <- reactive({
    req(input$cluster)
    
    current_dataset() %>%
      filter(cluster %in% input$cluster) 
  })
  
  selected_clusters_data_with_grouping <- reactive({
    req(input$groupBy)
    
    validate(need(input$groupBy != "No grouping", "Select a grouping criteria in the Group By panel"))
    
    data_df <- selected_clusters_data()
    
    if(input$groupBy == "Genotype") {
      genos <- genotypes_by_rsid()
      validate(need(!is.null(genos), "Enter a valid rsID"))
      
      data_df %<>% left_join(genos, by="individual") %>%
        mutate(metadata = ifelse(is.na(genotype), "unknown", genotype))
    }
    
    if(input$groupBy == "Individual metadata") {
      keep_columns <- c("individual", "individual_id", input$groupIndividualsBy)
      metadata_df <- individuals_df[, keep_columns]
      names(metadata_df)[3] <- "metadata"
      data_df %<>% left_join(metadata_df, by="individual")
    }
    
    data_df
  })
  
  clusterIndividualPlot <- reactive({
    req(input$cluster)
    req(input$groupBy)
    
    if(input$groupBy != "No grouping") {
      data_df <- selected_clusters_data_with_grouping()
    } else {
      data_df <- selected_clusters_data() %>%
        mutate(metadata = "none")
    }
    
    data_df %<>% group_by(gene, individual, metadata) %>%
      summarize(mean_exp = mean(exp),
                cells = n())

    if(nrow(data_df) == 0) return(NULL)
    
    max_exp <- max(data_df$mean_exp)
        
    g <- ggplot(data_df, aes(x=individual, y=mean_exp)) +
      geom_bar(stat="identity", fill="darkgreen") +
      scale_y_continuous(limits=c(0, max_exp*1.15), expand=c(0, 0)) +
      labs(x="Individual", 
           y="Mean expression",
           title=glue("Cluster {selected_clusters_label()} - average expression by individual")) +
      theme_bw(plot_base_text_size) +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    
    if(input$includeCounts == "Cells") {
      g <- g + geom_text(aes(label=pn(cells), y=Inf), vjust=1.5) 
    }
    
    if(input$groupBy != "No grouping") {
      g <- g + facet_grid(gene ~ metadata, scales="free_x", space="free_x")
    } else {
      g <- g + facet_wrap(~ gene, ncol=1, scales="free_y") 
    }
      
    g
  })
  
  clusterGroupedPlot <- reactive({

    data_df <- selected_clusters_data_with_grouping()
    
    max_exp <- max(data_df$exp)
    
    g <- ggplot(data_df, aes(x=metadata, y=exp)) +
      geom_violin(fill="gray80") +
      geom_point(size=0.5, alpha=0.5, position="jitter", color="black") +
      facet_wrap(~ gene, scales="free_y", ncol=1) +
      scale_y_continuous(limits=c(NA, max_exp*1.1)) +
      labs(x="Group", 
           y="Expression",
           title=glue("Cluster {selected_clusters_label()} - cell expression by selected grouping")) +
      theme_bw(plot_base_text_size) +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    
    if(input$includeCounts == "Cells") {
      cell_counts <- data_df %>%
        group_by(gene, metadata) %>%
        summarize(cells = n())
      
      g <- g + geom_text(data=cell_counts, inherit.aes=FALSE, 
                         aes(x=metadata, y=Inf, label=pn(cells)), vjust=1.5) 
    }
    
    g
  })
  
  grouped_stat_tests_table <- reactive({
    data_df <- selected_clusters_data_with_grouping() %>%
      filter(metadata != "unknown")
    
    genes <- data_df$gene %>% unique

    results_df <- genes %>%
      lapply(function(gene_name) {
        gene_data <- filter(data_df, .data$gene == gene_name)
        
        pairwise.wilcox.test(gene_data$exp, gene_data$metadata, p.adj="BH")$p.value %>% 
          as.data.frame %>% 
          tibble::rownames_to_column("group1") %>% 
          gather(group2, adj_pvalue, -group1) %>%
          filter(group1 != group2) %>%
          mutate(gene = gene_name,
                 neg_log10pv = -log10(adj_pvalue))
      }) %>%
      bind_rows %>%
      select(gene, everything())
    
    results_df %>% arrange(adj_pvalue)
  })
  
  output$grouped_stat_tests <- renderTable({
    grouped_stat_tests_table()
  })
  
  output$rsid_validation <- renderText({
    req(input$rsID)
    if(is.null(genotypes_by_rsid())) {
      return("No genotypes for that rsID")
    } else {
      return("Valid rsID")
    }
  })
  
  output$genotypes_table <- renderTable({
    genotypes_by_rsid()
  })
  
  output$plot_table <- renderTable({
    by_cluster_data()
  })
  
  output$by_cluster_plot <- renderPlot({
    byClusterPlot()
  })
  
  output$cluster_individual_plot <- renderPlot({
    clusterIndividualPlot()
  })
  
  output$cluster_grouped_plot <- renderPlot({
    clusterGroupedPlot()
  })
  
  output$database_info <- renderText({
    glue("{cell_count} cells and {records} billion records. Database: {db}", 
         cell_count = pn(nrow(cell_data())),
         records    = round(as.numeric(nrow(cell_data())) * nrow(genes_df) / 1e9, 2),
         db         = basename(datapath))
  })
  
  output$summary_plot_cells <- renderPlot({
    summary_df <- cell_data() %>%
      mutate(cluster = as.character(cluster)) %>%
      group_by(cluster) %>%
      summarize(cells = n(),
                individuals = length(unique(individual_id))) %>%
      arrange(cells)
    
    summary_df$cluster %<>% factor(levels=summary_df$cluster)
    
    max_cells <- max(summary_df$cells/1000)
    
    g <- ggplot(summary_df, aes(x=cluster, y=cells/1000)) +
      geom_bar(stat="identity") +
      scale_y_continuous(expand=c(0, 0), limits=c(0, max_cells*1.05)) +
      labs(x="Cluster", 
           y="Cells (x 1,000)",
           title="Cells per cluster") +
      theme_bw(plot_base_text_size) +
      coord_flip()
    
    g
  })

  output$summary_plot_individuals <- renderPlot({
    summary_df <- cell_data() %>%
      mutate(cluster = as.character(cluster)) %>%
      group_by(cluster) %>%
      summarize(cells = n(),
                individuals = length(unique(individual_id))) %>%
      arrange(cells)
    
    summary_df$cluster %<>% factor(levels=summary_df$cluster)
    
    max_indiv <- max(summary_df$individuals)
    
    g <- ggplot(summary_df, aes(x=cluster, y=individuals)) +
      geom_bar(stat="identity", fill="darkgreen") +
      scale_y_continuous(expand=c(0, 0), limits=c(0, max_indiv*1.05)) +
      labs(x="Cluster", 
           y="Individuals",
           title="Individuals per cluster") +
      theme_bw(plot_base_text_size) +
      coord_flip()
    
    g
  })
  
}

# Launch -----

log_message("Starting Shiny...")
shinyApp(ui = ui, server = server)

