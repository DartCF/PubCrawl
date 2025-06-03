#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(shinycssloaders)
library(bslib)
library(bsicons)
library(shinyjs)

# Custom theme
my_theme <- bs_theme(
  bg = "white", 
  fg = "black", 
  primary = "#205493",
  base_font = font_google("Roboto"),
  heading_font = font_google("Roboto")
) %>%  bs_add_rules(
  ".custom-hover {
    --bs-btn-hover-bg: #205493 !important;
  }"
)

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  theme = my_theme,  
  useShinyjs(),
  titlePanel(title = div(img(src="PubCrawl_Icon_Cropped.png", width=500, height=165), 
                         h5("Search PubMed by keywords and molecular symbols"),style="color:#205493"),
             windowTitle = "PubCrawl"),
  fluidRow(
    column(width=3, 
           card(
             card_header("Search Terms", class="bg-primary"),
             textInput("search_keyword", label = tooltip(
               trigger = list(
                 "Enter search keyword",
                 bs_icon("info-square")
               ),
               "Keyword is used in PubMed search query, can be any word or phrase that you want to find in Title/Abstract of PubMed entry."
             ), value = ""),
             fileInput("search_symbols", label = tooltip(
               trigger = list(
                 "Upload CSV file of molecule identifiers",
                 bs_icon("info-square")
               ),
               "A CSV file with identifiers in the first column (no header). Can be any molecule identifier (gene, protein, cytokine etc.) you want to find in Title/Abstract of PubMed entry."
             ), accept = "csv"),
             selectInput("retmax", "Number of results to return (most recent)", choices = c(5,10,15,20,25), multiple = F, selected = 15),
             textInput("search_apikey", label = tooltip(
               trigger = list(
                 HTML('Enter an optional <a href="https://support.nlm.nih.gov/kbArticle/?pn=KA-05316">NCBI API key</a>'),
                 bs_icon("info-square")
               ),
               "An NCBI API increases the number of requests/second a user can send to the NCBI server, improving the execution time of the search."
             ), value=""),
             actionButton("search_pubmed", "Search", class="custom-hover"),
             actionButton("reset", "Reset", class="custom-hover")  
           )
    ),
    column(width=9, 
           card(height = "100vh",
                card_header("Results Table", class="bg-primary"),
                uiOutput("error_message", style="color:red"),
                uiOutput("results_message", style="color:red"),
                downloadButton("download_results", label = "CSV Download", class = "custom-hover"),
                withSpinner(DT::dataTableOutput("results_table"), type = 8, color = "#205493")
           )
    )
  )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  source("./Lookup_Functions.R", local=T)
  # Reactive Input Monitor
  rv <- reactiveValues(keyword=NULL, 
                       symbols=NULL, 
                       retmax=NULL,
                       apikey=NULL, 
                       failure=NULL, 
                       no_hits=NULL,
                       search_triggered=FALSE)
  
  # Trigger PubMed Search
  observeEvent(input$search_pubmed,{
    if(is.null(input$search_symbols)){
      rv$failure = "You must upload a CSV file of symbols"
      return()
    }
    rv$symbols = read.csv(input$search_symbols$datapath, header=F)[,1]
    if(length(rv$symbols) > 20){
      rv$failure = "Too many gene symbols, maximum is 20"
      return()
    }
    rv$failure = NULL
    rv$keyword = input$search_keyword
    rv$retmax = input$retmax
    rv$apikey= input$search_apikey
    
    # build in buffer with 1.5 multiplier
    estimate <- ifelse(rv$apikey=="",
                       as.numeric(rv$retmax)*0.34*length(rv$symbols)*1.5,
                       as.numeric(rv$retmax)*0.1*length(rv$symbols)*1.5)
    showNotification(paste0("Search time estimate: ", estimate, " seconds"),
                     type="message",
                     duration=5)
    
    rv$search_triggered=TRUE
    
  })
  
  # Search PubMed
  SearchPubMed <- reactive({
    if(rv$search_triggered){
      results <- NULL
      keyword <- isolate(rv$keyword)
      apikey <- isolate(rv$apikey)
      retmax <- isolate(rv$retmax)
      symbols <- isolate(rv$symbols)
      # set API key if one was provided
      if (!(apikey=="")){
        status <- set_ncbi_key(apikey)
      }else{
        apikey=NULL
      }
      result.df <- data.frame()
      tryCatch({
        result.list <- lapply(symbols,function(x){
          fetch_pubmed_data(keyword, x, max_results=retmax, apikey)})
        result.df <- do.call(rbind, result.list)
      }, error = function(e){
        rv$failure = e$message
      })
      no_hits <- setdiff(symbols, result.df$Symbol)
      if (length(no_hits) > 0){
        rv$no_hits <- no_hits
      }
      return(result.df)
    }
  })
  
  # Render results data table
  output$results_table <- DT::renderDataTable({
    req(rv$search_triggered == TRUE)
    datatable(SearchPubMed(),
              options = list(scrollX=TRUE, 
                             dom = 'frtip',
                             autoWidth=TRUE),
              escape = FALSE)
    
  })
  
  # Input listeners (reset search_triggered, failure and no hits reactive variables)
  observeEvent(input$search_keyword, {
    rv$search_triggered = FALSE
    rv$failure = NULL
    rv$no_hits = NULL
  })
  
  observeEvent(input$search_symbols, {
    rv$search_triggered = FALSE
    rv$failure = NULL
    rv$no_hits = NULL
  })
  
  observeEvent(input$retmax, {
    rv$search_triggered = FALSE
    rv$failure = NULL
    rv$no_hits = NULL
  })
  
  observeEvent(input$search_apikey, {
    rv$search_triggered = FALSE
    rv$failure = NULL
    rv$no_hits = NULL
  })
  
  # Download search results table
  output$download_results <- downloadHandler(
    filename = function(){
      date <- gsub("-","",Sys.Date())
      time <- gsub("[:\\.]","",str_split_i(Sys.time(), " ", 2))
      paste0("pubmed_search_results_",date,"_",time,".csv")
    },
    content = function(file){
      results <- SearchPubMed()
      results[,"PubMed_Link"] <- str_extract(results[,"PubMed_Link"],"https\\:\\/\\/pubmed\\.ncbi\\.nlm\\.nih.\\gov\\/[0-9]+")
      write.csv(results, file, row.names = F)
    }
  )
  
  # Render error message if PubMed search failed
  output$error_message <- renderUI({
    if(!is.null(rv$failure)){
      renderText({paste0("Error: ",rv$failure)})
    }
  })
  
  # Render results message if search didn't return any hits
  output$results_message <- renderUI({
    if(!is.null(rv$no_hits)){
      renderText({paste0("No results for: ", rv$no_hits)})
    }
  })
  # Reset inputs and UI elements
  observeEvent(input$reset, {
    # reset reactive values
    rv$keyword = rv$symbols = rv$apikey = rv$retmax = rv$failure = rv$no_hits = NULL
    rv$search_triggered = FALSE
    # clear text input elements
    updateTextInput(session, "search_keyword",value="")
    shinyjs::reset("search_symbols")
    updateTextInput(session, "search_apikey",value="")
    updateSelectInput(session, "retmax", choices = c(5,10,15,20,25), selected = 15)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
