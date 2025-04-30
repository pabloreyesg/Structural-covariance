library(shiny)
library(igraph)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stats)
library(visNetwork)
library(BayesFactor)

ui <- fluidPage(
  titlePanel("Structural Covariance Network (SCN) Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV", accept = ".csv"),
      textInput("group_col", "Group column name", value = "Group"),
      textInput("covariate_col", "Covariate column name (optional)", value = "Edad"),
      numericInput("threshold", "Correlation Threshold", value = 0.3, min = 0, max = 1, step = 0.05),
      selectInput("group_select", "Group to Display", choices = NULL),
      selectInput("group1", "Group 1 for Comparison", choices = NULL),
      selectInput("group2", "Group 2 for Comparison", choices = NULL),
      selectInput("test_type", "Statistical Test", 
                  choices = c("t-test (parametric)" = "ttest", 
                              "Wilcoxon (non-parametric)" = "wilcox", 
                              "Bayesian (permutation)" = "bayesperm")),
      selectInput("layout_type", "Graph Layout", 
                  choices = c("circle", "layout_nicely", "layout_with_fr", "layout_with_kk", "layout_as_tree"), 
                  selected = "layout_nicely"),
      actionButton("run", "Run Analysis")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SCN Heatmap", plotOutput("heatmap")),
        tabPanel("Graph", visNetworkOutput("graph_plot")),
        tabPanel("Graph Metrics", tableOutput("metrics")),
        tabPanel("Statistical Comparison", tableOutput("comparison")),
        tabPanel("Difference Heatmap", plotOutput("diff_heatmap"))
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveVal()
  adjusted_data <- reactiveVal()
  graphs <- reactiveVal(list())
  scns <- reactiveVal(list())
  
  observeEvent(input$file, {
    req(input$file)
    df <- read.csv(input$file$datapath)
    data(df)
    group_levels <- unique(df[[input$group_col]])
    updateSelectInput(session, "group_select", choices = group_levels)
    updateSelectInput(session, "group1", choices = group_levels)
    updateSelectInput(session, "group2", choices = group_levels, selected = group_levels[2])
  })
  
  observeEvent(input$run, {
    req(data())
    df <- data()
    region_cols <- setdiff(names(df), c(input$group_col, input$covariate_col))
    
    if (input$covariate_col %in% names(df)) {
      adjusted <- df
      for (col in region_cols) {
        model <- lm(df[[col]] ~ df[[input$covariate_col]])
        adjusted[[col]] <- residuals(model)
      }
    } else {
      adjusted <- df
    }
    
    adjusted_data(adjusted)
    
    scn_list <- list()
    group_levels <- unique(adjusted[[input$group_col]])
    for (grp in group_levels) {
      subset <- adjusted %>% filter(.[[input$group_col]] == grp) %>% select(all_of(region_cols))
      cor_mat <- cor(subset)
      scn_list[[grp]] <- cor_mat
    }
    scns(scn_list)
    
    graph_list <- list()
    for (grp in names(scn_list)) {
      mat <- scn_list[[grp]]
      mat[abs(mat) < input$threshold] <- 0
      diag(mat) <- 0
      G <- graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE, diag = FALSE)
      G <- delete_vertices(G, which(degree(G) == 0))
      graph_list[[grp]] <- G
    }
    graphs(graph_list)
  })
  
  output$heatmap <- renderPlot({
    req(scns())
    mat <- scns()[[input$group_select]]
    if (is.null(mat) || length(mat) == 0) return(NULL)
    melted_mat <- melt(mat)
    ggplot(melted_mat, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste("SCN Heatmap -", input$group_select))
  })
  
  output$graph_plot <- renderVisNetwork({
    req(graphs())
    G <- graphs()[[input$group_select]]
    
    nodes <- data.frame(id = V(G)$name, label = V(G)$name)
    edges <- igraph::as_data_frame(G, what = "edges")
    colnames(edges)[1:2] <- c("from", "to")
    
    vis <- visNetwork(nodes, edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visPhysics(enabled = FALSE)
    
    if (input$layout_type == "circle") {
      vis <- vis %>% visLayout(improvedLayout = TRUE) %>% visIgraphLayout(layout = "layout_in_circle")
    } else {
      vis <- vis %>% visIgraphLayout(layout = input$layout_type)
    }
    
    vis
  })
  
  output$metrics <- renderTable({
    req(graphs())
    G <- graphs()[[input$group_select]]
    positive_weights <- E(G)$weight
    if (any(positive_weights < 0)) {
      positive_weights <- abs(positive_weights)
      E(G)$weight <- positive_weights
    }
    data.frame(
      Density = edge_density(G),
      Clustering_Coefficient = transitivity(G, type = "average"),
      Average_Path_Length = ifelse(is_connected(G), mean_distance(G), NA)
    )
  })
  
  output$comparison <- renderTable({
    req(scns())
    mat1 <- scns()[[input$group1]]
    mat2 <- scns()[[input$group2]]
    
    upper_tri <- upper.tri(mat1)
    diff_vec <- mat1[upper_tri] - mat2[upper_tri]
    
    if (input$test_type == "ttest") {
      test_res <- t.test(diff_vec)
      stat <- test_res$statistic
      pval <- test_res$p.value
      method <- "t-test"
    } else if (input$test_type == "wilcox") {
      test_res <- wilcox.test(diff_vec)
      stat <- test_res$statistic
      pval <- test_res$p.value
      method <- "Wilcoxon test"
    } else {
      bf_res <- ttestBF(x = diff_vec, rscale = "medium")
      stat <- as.vector(extractBF(bf_res)$bf)
      pval <- NA
      method <- "Bayes Factor (t)"
    }
    
    data.frame(
      Group1 = input$group1,
      Group2 = input$group2,
      Test = method,
      Mean_Diff = mean(diff_vec),
      Statistic = stat,
      p_value = pval
    )
  })
  
  output$diff_heatmap <- renderPlot({
    req(scns())
    mat1 <- scns()[[input$group1]]
    mat2 <- scns()[[input$group2]]
    diff_mat <- mat1 - mat2
    if (is.null(diff_mat) || length(diff_mat) == 0) return(NULL)
    melted_diff <- melt(diff_mat)
    ggplot(melted_diff, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste("SCN Difference Heatmap -", input$group1, "vs", input$group2))
  })
}

shinyApp(ui, server)
