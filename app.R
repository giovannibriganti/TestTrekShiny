#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(car)


recommend_test_independent <- function(data, var1, var2) {
  
  # Determine the nature of the two variables
  nature_var1 <- ifelse(is.numeric(data[[var1]]), "continuous", "nominal")
  nature_var2 <- ifelse(is.numeric(data[[var2]]), "continuous", "nominal")
  
  # Initialize variables to store recommendations
  test <- ""
  effect_size <- ""
  posthoc <- NULL
  
  # Check for normality
  p_val_norm_var1 <- if (nature_var1 == "continuous") shapiro.test(data[[var1]])$p.value else 1
  p_val_norm_var2 <- if (nature_var2 == "continuous") shapiro.test(data[[var2]])$p.value else 1
  
  is_normal_var1 <- p_val_norm_var1 > 0.05
  is_normal_var2 <- p_val_norm_var2 > 0.05
  
  # Recommend test based on variable type for independent samples
  if (nature_var1 == "continuous" & nature_var2 == "continuous") {
    if (is_normal_var1 && is_normal_var2) {
      test <- "Pearson's correlation"
      effect_size <- "r"
    } else {
      test <- "Spearman's rho or Kendall's tau"
      effect_size <- "rho or tau"
    }
  }
  else if ((nature_var1 == "nominal" && nature_var2 == "continuous") |
           (nature_var1 == "continuous" && nature_var2 == "nominal")) {
    continuous_var <- ifelse(nature_var1 == "continuous", var1, var2)
    nominal_var <- ifelse(nature_var1 == "nominal", var1, var2)
    num_groups <- length(unique(data[[nominal_var]]))
    
    # Check for heteroskedasticity
    p_val_levene <- leveneTest(data[[continuous_var]], as.factor(data[[nominal_var]]), center = median)$"Pr(>F)"[1]
    is_homogeneous_variance <- p_val_levene > 0.05
    
    if (num_groups == 2) {
      if (is_normal_var1 && is_normal_var2) {
        if (is_homogeneous_variance) {
          test <- "Independent samples t-test"
          effect_size <- "Cohen's d"
        } else {
          test <- "Welch's t-test"
          effect_size <- "Cohen's d (adjusted)"
        }
      } else {
        test <- "Mann-Whitney U test"
        effect_size <- "r (effect size for non-parametric tests)"
      }
    }
    else if (num_groups > 2) {
      if (is_normal_var1 && is_normal_var2) {
        if (is_homogeneous_variance) {
          test <- "ANOVA"
          posthoc <- "Tukey's HSD"
          effect_size <- "eta squared"
        } else {
          test <- "Welch's ANOVA"
          posthoc <- "Games-Howell posthoc test"
          effect_size <- "omega squared"
        }
      } else {
        test <- "Kruskal-Wallis H test"
        posthoc <- "Dunn's posthoc test"
        effect_size <- "epsilon squared"
      }
    }
  }
  else if (nature_var1 == "nominal" & nature_var2 == "nominal") {
    if (length(unique(data[[var1]])) == 2 && length(unique(data[[var2]])) == 2) {
      test <- "Chi-squared test or Fisher's Exact Test (for small sample sizes)"
    } else {
      test <- "Chi-squared test"
    }
    effect_size <- "Phi or Cramer's V"
  }
  
  # Build the recommendation string
  recommendation <- paste("Recommended Test:", test, "\n",
                          "Effect Size:", effect_size, "\n")
  
  if (!is.null(posthoc)) {
    recommendation <- paste(recommendation, "Post hoc test:", posthoc, "\n")
  }
  
  return(recommendation)
}


#' Recommend Appropriate Statistical Test for Dependent Samples
#'
#' This function analyses two variables from a dataset and recommends
#' the appropriate statistical test and effect size for dependent samples.
#'
#' @param data A data frame containing the variables to analyze.
#' @param var1 The name of the first variable.
#' @param var2 The name of the second variable.
#'
#' @return Prints the recommended statistical test and effect size.
#'
#' @examples
#' \dontrun{
#' recommend_test_dependent(my_data, "var1", "var2")
#' }
#'
#' @export

recommend_test_dependent <- function(data, var1, var2) {
  # Determine the nature of the two variables
  nature_var1 <- ifelse(is.numeric(data[[var1]]), "continuous", "nominal")
  nature_var2 <- ifelse(is.numeric(data[[var2]]), "continuous", "nominal")
  
  # Initialize variables to store recommendations
  test <- ""
  effect_size <- ""
  posthoc <- NULL
  
  # Check for normality
  p_val_norm_var1 <- if (nature_var1 == "continuous") shapiro.test(data[[var1]])$p.value else 1
  p_val_norm_var2 <- if (nature_var2 == "continuous") shapiro.test(data[[var2]])$p.value else 1
  
  is_normal_var1 <- p_val_norm_var1 > 0.05
  is_normal_var2 <- p_val_norm_var2 > 0.05
  
  # Recommendations for dependent samples
  if (nature_var1 == "continuous" && nature_var2 == "continuous") {
    if (is_normal_var1 && is_normal_var2) {
      test <- "Paired samples t-test"
      effect_size <- "d (paired)"
    } else {
      test <- "Wilcoxon signed-rank test"
      effect_size <- "r (effect size for non-parametric tests)"
    }
  }
  else if ((nature_var1 == "nominal" && nature_var2 == "continuous") ||
           (nature_var1 == "continuous" && nature_var2 == "nominal")) {
    continuous_var <- ifelse(nature_var1 == "continuous", var1, var2)
    nominal_var <- ifelse(nature_var1 == "nominal", var1, var2)
    num_groups <- length(unique(data[[nominal_var]]))
    
    if (num_groups == 2) {
      test <- "Repeated measures ANOVA"
      effect_size <- "Partial eta squared"
      
      # Checking sphericity with Mauchly's test
      p_val_mauchly <- mauchly.test(data[, c(var1, var2)])$p.value
      if (p_val_mauchly <= 0.05) {
        cat("Mauchly's test indicates that the assumption of sphericity has been violated.\n")
        cat("Consider using Greenhouse-Geisser or Huynh-Feldt corrections.\n")
      }
    }
    else {
      test <- "Friedman test"
      effect_size <- "Kendall's W"
    }
  }
  else if (nature_var1 == "nominal" && nature_var2 == "nominal") {
    test <- "McNemar's test"
    effect_size <- "Phi coefficient"
  }
  
  # Build the recommendation string
  recommendation <- paste("Recommended Test:", test, "\n",
                          "Effect Size:", effect_size, "\n")
  
  if (!is.null(posthoc)) {
    recommendation <- paste(recommendation, "Post hoc test:", posthoc, "\n")
  }
  
  return(recommendation)
}






# Rest of the Shiny app (UI, server, shinyApp call)

ui <- fluidPage(
  
  # Application title
  titlePanel("TestTrek: Statistical Test Recommender"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Inputs
    sidebarPanel(
      fileInput("data", "Upload your data (CSV or Excel)", accept = c(".csv", ".xlsx")),
      uiOutput("var1_ui"),
      uiOutput("var2_ui"),
      radioButtons("sample_type", "Sample type:", choices = c("Independent", "Dependent"), selected = "Independent"),
      actionButton("analyze", "Analyze")
    ),
    
    # Outputs
    mainPanel(
      tableOutput("data_head"),
      textOutput("recommendation")
    )
  )
)

server <- function(input, output) {
  
  # Reactive expression to read data
  dataset <- reactive({
    req(input$data)
    if(grepl(".csv$", input$data$name)) {
      # Check the separator used in the CSV file
      first_line <- readLines(input$data$datapath, n = 1)
      sep_char <- ifelse(grepl(";", first_line), ";", ",")
      
      read.csv(input$data$datapath, sep = sep_char)
    } else {
      readxl::read_excel(input$data$datapath)
    }
  })
  
  # Show the first few rows of the uploaded data
  output$data_head <- renderTable({
    head(dataset())
  })
  
  # Dynamic UI for variable selection
  output$var1_ui <- renderUI({
    req(dataset())
    selectInput("var1", "Select first variable:", names(dataset()))
  })
  
  output$var2_ui <- renderUI({
    req(dataset())
    selectInput("var2", "Select second variable:", names(dataset()))
  })
  
# Recommendation 
  observeEvent(input$analyze, {
    data <- dataset()
    var1 <- input$var1
    var2 <- input$var2
    if(input$sample_type == "Independent") {
      recommendation <- recommend_test_independent(data, var1, var2)
    } else {
      recommendation <- recommend_test_dependent(data, var1, var2)
    }
    output$recommendation <- renderText({
      recommendation
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
