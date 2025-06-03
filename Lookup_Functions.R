# Functions used in Shiny app to retrieve pubmed search results for given gene 
# symbole and search keyword

library(rentrez)
library(tidyverse) 

set_ncbi_key <- function(api_key){
  val = set_entrez_key(api_key)
  return(val)
}

fetch_pubmed_data <- function(base_term, symbol, max_results, apikey) {
  # Construct the search term
  
  search_term <- paste0("(", base_term, "[Title/Abstract]) AND (", symbol, "[Title/Abstract])")
  
  # Perform the Entrez search
  search_results <- entrez_search(db = "pubmed", term = search_term, retmax = max_results, sort="pub_date")
  
  if (length(search_results$ids) == 0) {
    return(NULL)
  }
  
  # Initialize an empty data frame to store results
  results_df <- data.frame(
    Symbol = character(0),
    Base_Term = character(0),
    PubMed_ID = character(0),
    Author = character(0),
    Date = character(0),
    Title = character(0),
    PubMed_Link = character(0),
    stringsAsFactors = FALSE
  )
  
  # Fetch summary for each PMID
  for (pmid in search_results$ids) {
    link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid)
    hyperlink <- paste('<a href="',link,'" target="_blank">', link, '</a>', collapse="</br>")
    tryCatch({
      esummary <- entrez_summary(db = "pubmed", id = pmid)
      if (length(esummary) > 0 && !is.null(esummary)) {
        title <- esummary$title
        pubdate <- esummary$pubdate
        first_author <- ifelse(length(esummary$authors) > 0, esummary$authors[1,"name"], NA_character_)
        results_df <- rbind(results_df, data.frame(Symbol = symbol, Base_Term = base_term, PubMed_ID = pmid, Author = first_author, Date = pubdate, Title = title, PubMed_Link = hyperlink, stringsAsFactors = FALSE))
      } else {
        results_df <- rbind(results_df, data.frame(Symbol = symbol, Base_Term = base_term, PubMed_ID = pmid, Author = NA_character_, Date = NA_character_, Title = NA_character_, PubMed_Link = hyperlink, stringsAsFactors = FALSE))
      }
    }, error = function(e) {
      message(paste("Error fetching summary for PMID", pmid, ":", e$message))
      results_df <- rbind(results_df, data.frame(Symbol = symbol, Base_Term = base_term, PubMed_ID = pmid, Author = "Error", Date = "Error", Title = "Error", PubMed_Link = hyperlink, stringsAsFactors = FALSE))
    })
    if (!is.null(apikey)){
      Sys.sleep(0.1) # API key rate limit is 10 queries / second
    }else{
      Sys.sleep(0.34) # no key rate limit is 3 queries / second
    }
    
  }
  
  return(results_df)
}