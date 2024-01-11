fetchGHdata <- function(gh_account, repo, path, output) {
    
    # Dependencies
    require(tidyverse)
    require(httr) 
    require(rlist)
    require(jsonlite)
    
    # First you have to authenticate.
    # Store a personal access token in .Renviron
    # See https://blog.exploratory.io/extract-data-from-private-github-repository-with-rest-api-db804fa43d84
    auth <- authenticate(Sys.getenv("GITHUB_PAT"), "")
    accpt <- accept("application/vnd.github.v3.raw")
    
    # Seperate the filename from the directory
    match <- regexpr("^(.*[\\/])", path)
    if (match[1] > 0) {
        dir <- path %>% substring(match[1], attr(match, "match.length"))
        file <- path %>% substring(attr(match, "match.length") + 1, nchar(path))
    } else {
        dir <- ""
        file <- path
    }
    
    # To handle files larger than 1MB, use this trick:
    # https://medium.com/@caludio/how-to-download-large-files-from-github-4863a2dbba3b
    req_meta <- 
        content(
            GET(
                paste("https://api.github.com/repos", gh_account, repo, "contents", dir, sep="/"), 
                auth
            )
        )
    
    entry <- req_meta %>% list.filter(name == file)
    sha <- entry[[1]]$sha
    
    # Download file
    system(
        paste("curl -H", dQuote(paste("Authorization: token", Sys.getenv("GITHUB_PAT")),q = "C"),
              "-H", dQuote("Accept: application/vnd.github.v3.raw+json",q = "C"),
              "-o",dQuote(output,q = "C"),
              "-L", dQuote(paste("https://api.github.com/repos", gh_account, repo, "git/blobs", sha, sep="/"),q = "C"))
    )
    
}

