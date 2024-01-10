# spnames = c("Alchemilla arvensis", 
#             "Arabis bellidifolia", 
#             "Elymus caput-medusae",
#             "Perezia microcephala",
#             "Polystachya estrellensis",
#             "Tachigali rubiginosa",
#             "Oxalis rhombeo ovata",
#             "Axonopus canescens")
# 
# clean_sp_names <- Find_Sci_Names(spnames)
# clean_sp_names

Find_Sci_Names <- function(spnames, 
                           db = c("gbif","itis","ncbi","iucn"),
                           priority = "gbif",
                           suggest_names = FALSE, 
                           return_accepted_only = FALSE){
    
    require("bdc")
    require("taxadb")
    require("traitdataform")
    require("parallel")
    require("dplyr")
    require("data.table")
    require("pbapply")
    
    # Create empty table 
    
    tofind <- data.frame(matrix(nrow = length(spnames), ncol = 8))
    names(tofind) = c("scientificName", "kingdom", "phylum", "class", "order", "family", "db", "db_code")
    tofind <- data.frame(species = spnames, tofind)
    
    for(i in 1:length(db)){ 
        
        cat("\n\n", db[i],"\n\n")
        
        if(db[i] == "gbif"){
            
            # retrieve sp names
            togo <- tofind[which(is.na(tofind$scientificName)),]
            
            cl <- makeCluster(detectCores()-2)
            clusterExport(cl, c("togo","standardize_taxa"))
            
            new_names <- pblapply(togo$species, function(x){
                try(
                    traitdataform::standardize_taxa(
                        data.frame(verbatimScientificName = x), 
                        fuzzy = FALSE,
                        silent = TRUE)
                    )
            }, cl = cl)
            
            stopCluster(cl)
            
            rem <- sapply(new_names, class)
            rem <- which(rem=="try-error")
            if(any(rem)){
                new_names <- new_names[-rem]
            }
            new_names <- rbindlist(new_names)
            
            new_names <- new_names[-which(is.na(new_names$scientificName)),]
            new_names <- new_names[which(new_names$taxonRank=='species'),]
            
            # remove duplicates
            if(any(duplicated(new_names$verbatimScientificName))){
                new_names <- new_names[-which(duplicated(new_names$verbatimScientificName)),]
            }
            
            cat("--- Summary ---\n",
                "N taxa:",nrow(togo),"\n",
                "N taxa found:",nrow(new_names), "\n",
                "N taxa not found:", nrow(togo)-nrow(new_names))
            
            if(!nrow(new_names) == 0){
                new_names <- new_names[,c("verbatimScientificName","scientificName","kingdom","phylum","class","order","family","taxonID")]
                names(new_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
                new_names$db <- "gbif"
                new_names$db_code <- gsub("http://www.gbif.org/species/","",new_names$db_code)
                new_names$db_code <- paste("GBIF:",new_names$db_code,sep = "")
                
                # Feed
                tofind <- tofind %>% 
                    rows_patch(new_names, 
                               by = "species")
            }
            
        } 
        if(db[i]=="iucn"){
            
            # retrieve sp names
            togo <- tofind[which(is.na(tofind$scientificName)),]
            
            new_names<-bdc_query_names_taxadb(togo$species, 
                                              db = "iucn",
                                              suggest_names = FALSE) # using exact names
            new_names$scientificName <- paste(new_names$genus,new_names$specificEpithet) 
            
            # remove duplicates
            if(any(duplicated(new_names$original_search))){
                new_names <- new_names[-which(duplicated(new_names$original_search)),]
            }
            
            # select only accepted names
            new_names <- new_names[which(new_names$taxonomicStatus == "accepted"),]
            
            cat("--- Summary ---\n",
                "N taxa:",nrow(togo),"\n",
                "N taxa found:",length(which(new_names$taxonomicStatus == "accepted")), "\n",
                "N taxa not found:", nrow(togo)-length(which(new_names$taxonomicStatus == "accepted")))
            
            if(!nrow(new_names) == 0){
                new_names <- new_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
                names(new_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
                new_names$db <- "iucn"
                
                # Feed
                tofind <- tofind %>% 
                    rows_patch(new_names, 
                               by = "species")
            }
            
        } 
        if(db[i]=="itis" | db[i]=="ncbi"){
            # retrieve sp names
            togo <- tofind[which(is.na(tofind$scientificName)),]
            
            # retrieve sp names
            new_names <- Find_Names(spnames = togo$species,
                                    db = "itis",
                                    suggest_names = FALSE, # using exact names
                                    return_accepted_only = TRUE) 
            
            # remove duplicates
            if(any(duplicated(new_names$original_search))){
                new_names <- new_names[-which(duplicated(new_names$original_search)),]
            }
            
            cat("\n--- Summary ---\n",
                "N taxa:",nrow(togo),"\n",
                "N taxa found:",length(which(new_names$taxonomicStatus == "accepted")), "\n",
                "N taxa not found:", nrow(togo)-length(which(new_names$taxonomicStatus == "accepted")))
            
            if(!nrow(new_names) == 0){
                
                new_names <- new_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
                names(new_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
                new_names$db <- "itis"
                
                # Feed
                tofind <- tofind %>% 
                    rows_patch(new_names, 
                               by = "species")
            }
            
        }
    }
    
    #################################
    # Fix issues with scientific names
    if(any(!is.na(tofind$scientificName))){
        new <- sapply(tofind$scientificName, function(x){
            tmp <- strsplit(x," ")[[1]]
            if(any(duplicated(tmp))){
                paste(tmp[-1], collapse = " ")
            } else {
                x
            }
        })
        
        tofind$scientificName <- new
    }
    
    #################################
    # Priority
    cat("\nPrioriotizing", priority)
    
    # retrieve sp names
    if(priority == "gbif"){
        
        # retrieve sp names
        togo <- tofind[which(!tofind$db=="gbif"),]
        
        cl <- makeCluster(detectCores()-2)
        clusterExport(cl, c("togo","standardize_taxa"))
        
        new_names <- pblapply(togo$species, function(x){
            try(standardize_taxa(data.frame(verbatimScientificName = x), 
                                 fuzzy = FALSE,
                                 silent = TRUE))
        }, cl = cl)
        
        stopCluster(cl)
        
        rem <- sapply(new_names, class)
        rem <- which(rem=="try-error")
        if(any(rem)){
            new_names <- new_names[-rem]
        }
        new_names <- rbindlist(new_names)
        
        new_names <- new_names[-which(is.na(new_names$scientificName)),]
        new_names <- new_names[which(new_names$taxonRank=='species'),]
        
        # remove duplicates
        if(any(duplicated(new_names$verbatimScientificName))){
            new_names <- new_names[-which(duplicated(new_names$verbatimScientificName)),]
        }
        
        cat("--- Summary ---\n",
            "N taxa:",nrow(togo),"\n",
            "N taxa found:",nrow(new_names), "\n",
            "N taxa not found:", nrow(togo)-nrow(new_names))
        
        if(!nrow(new_names) == 0){
            new_names <- new_names[,c("verbatimScientificName","scientificName","kingdom","phylum","class","order","family","taxonID")]
            names(new_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
            new_names$db <- "gbif"
            new_names$db_code <- gsub("http://www.gbif.org/species/","",new_names$db_code)
            new_names$db_code <- paste("GBIF:",new_names$db_code,sep = "")
            
            # Feed
            tofind <- tofind %>% 
                rows_patch(new_names, 
                           by = "species")
        }
        
        
    }
    
    return(tofind)
}



# x="Arabis bellidifolia"

Clean_Names <- function(names, show.progress = FALSE, return_gen_sps = TRUE){
    res <- c()
    for(x in 1:length(names)){
        if(show.progress){
            cat("\r",x,"from",length(names))
        }
        if(is.na(names[x])){
            # do nothing
        }
        else{
            # replace artifacts from encoding issue
            tmp. <- gsub('<U\\+FB01>', 'fi', names[x])
            tmp. <- gsub('ÿ', '', tmp.)
            # dot within name
            tmp. <- strsplit(tmp.,"[.]")[[1]]
            tmp. = paste(tmp., collapse = " ")
            # reduces repeated white space
            tmp. <- str_squish(tmp.)
            # remove letters with accent
            tmp. <- iconv(tmp., "latin1", "ASCII", "")
            # remove parentheses
            tmp. <- strsplit(tmp., " ")[[1]]
            if(any(grepl("[(]", tmp.))){
                tmp. = tmp.[-which(grepl("[(]", tmp.))]
            }
            # remove bar
            if(any(grepl("[/]", tmp.))){
                pos <- grep("[/]", tmp.)
                tmp.[pos] = gsub("[/]"," ",tmp.[pos])
                tmp. = paste(tmp., collapse = " ")
                tmp. <- strsplit(tmp., " ")[[1]]
            }
            # Genus uppercase and species lowercase
            tmp.1 = str_to_title(tmp.[1])
            if(length(tmp.)>1){
                tmp.2 <- sapply(tmp.[-1], function(x) str_to_lower(x))
                tmp.2 = paste(tmp.1, paste(tmp.2, collapse = " "), collapse = " ")
                tmp. <- strsplit(tmp.2, " ")[[1]]
            } else {
                tmp. <- tmp.1
            }
            # fix sp. var. subspecies
            tmp.[which(tmp. == "sp")] <- "sp."
            tmp.[which(tmp. == "spp")] <- "sp."
            tmp.[which(tmp. == "spp.")] <- "sp."
            tmp.[which(tmp. == "var")] <- "var."
            tmp.[which(tmp. == "cf")] <- "cf."
            tmp.[which(tmp. == "ssp")] <- "subsp."
            tmp.[which(tmp. == "ssp.")] <- "subsp."
            tmp.[which(tmp. == "subsp")] <- "subsp."
            tmp. <- gsub("[..]",".",tmp.)
            if(return_gen_sps) {
                # remove var subsp
                if(any(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")){
                    tmp. = tmp.[-which(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")]
                }
                if(length(tmp.)>2){
                    tmp. <- tmp.[1:2]
                }
            }
            if(length(tmp.)==1){
                tmp. = paste(tmp., "sp.", collapse = " ")
            } else {
                tmp. = paste(tmp., collapse = " ")
            }
        }
        
        res = c(res,tmp.)
    }
    return(res)
}


####################

Find_Names <- function(spnames, 
                       db = "gbif", 
                       suggest_names = FALSE, 
                       return_accepted_only = FALSE,
                       parallel = F,
                       ncores = detectCores()-2){
    cat("Searching names for",length(spnames),"taxa\n")
    
    sps_names<-query_names(sci_name = spnames, 
                           db = db,
                           suggest_names = suggest_names,
                           parallel = parallel,
                           ncores = ncores) 
    sps_names$scientificName <- paste(sps_names$genus,sps_names$specificEpithet)
    sps_names$scientificName[which(sps_names$scientificName=="NA NA")] <- NA
    if(any(duplicated(sps_names$original_search))){
        # keep the one with info
        dups <- unique(sps_names$original_search[which(duplicated(sps_names$original_search))])
        dups <- which(sps_names$original_search %in% dups)
        tokeep <- sps_names[dups,]
        tokeep <- tokeep %>%
            arrange(original_search, is.na(taxonID)) %>% 
            distinct(original_search, .keep_all = TRUE)
        sps_names <- sps_names[-dups,]
        sps_names <- rbind(sps_names,tokeep)
    }
    # all(spnames %in% sps_names$original_search)
    
    found <- length(unique(sps_names$original_search[
        grep("accepted",sps_names$taxonomicStatus)]))
    cat("Found accepted names for",found,"out of",length(spnames),"taxa\n")
    
    ## find accepted names for synonyms from accepted id
    cat("Searching accepted names for synonims from accepted id\n")
    ## Method 1 - from Accepted ID
    pos <- which(!sps_names$taxonomicStatus == "accepted" & !is.na(sps_names$acceptedNameUsageID))
    if(any(pos)){
        cat("Method 1 - Get names from accepted ID\n")
        cat("Searching for",length(pos),"taxa\n")
        togosyno <- sps_names[pos,]
        ## find names from ID
        togosyno_accepted <- find_names_from_id(id = togosyno$acceptedNameUsageID,
                                                name = togosyno$original_search,
                                                db = db)
        # feed
        togosyno_accepted <- togosyno_accepted[,which(colnames(togosyno_accepted) %in% colnames(sps_names))]
        # 1) rem species found
        rem <- which(sps_names$original_search %in% togosyno_accepted$original_search)
        if(any(rem)){
            sps_names <- sps_names[-rem,]
        }
        # 2) add species found
        sps_names <- bind_rows(sps_names,
                               togosyno_accepted)
        
        found <- length(unique(togosyno_accepted$original_search[
            grep("accepted",togosyno_accepted$taxonomicStatus)]))
        
        cat("Found accepted names for",found,"taxa\n")
    }
    ## Method 2 - using function taxadb::filter_name
    pos <- which(!sps_names$notes == "accepted")
    if(any(pos)){
        cat("Method 2 - Filter names \n")
        cat("Searching for",length(pos),"taxa\n")
        togosyno <- sps_names[pos,]
        togosyno_accepted <- find_names_filter(x = togosyno$original_search,
                                               db = db)
        togosyno_accepted <- togosyno_accepted[,which(colnames(togosyno_accepted) %in% colnames(sps_names))]
        # if multiple names suggested for the same species, get the accepted name based on taxa id
        dups <- which(duplicated(togosyno_accepted$original_search))
        if(length(dups)>0){
            dups <- unique(togosyno_accepted$original_search[dups])
            dups <- togosyno_accepted[which(togosyno_accepted$original_search %in% dups),]
            tmp <- find_names_from_id(id = dups$acceptedNameUsageID,
                                      name = dups$original_search,
                                      db = db)
            # 1) rem species found
            rem <- which(togosyno_accepted$original_search %in% tmp$original_search)
            if(any(rem)){
                togosyno_accepted <- togosyno_accepted[-rem,]
            }
            # 2) add species found
            togosyno_accepted <- bind_rows(togosyno_accepted,
                                           tmp)
        }
        # feed
        # 1) rem species found
        rem <- which(sps_names$original_search %in% togosyno_accepted$original_search)
        if(any(rem)){
            sps_names <- sps_names[-rem,]
        }
        # 2) add species found
        sps_names <- bind_rows(sps_names,
                               togosyno_accepted)
        
        found <- length(unique(togosyno_accepted$original_search[
            grep("accepted",togosyno_accepted$taxonomicStatus)]))
        
        cat("Found accepted names for",found,"taxa\n")
    }
    
    cat("--- Summary ---\n",
        "N taxa:",nrow(sps_names),"\n",
        "N taxa found:",length(grep("accepted",sps_names$taxonomicStatus)), "\n",
        "N taxa not found:", nrow(sps_names)-length(grep("accepted",sps_names$taxonomicStatus)))
    if(return_accepted_only){
        # select only accepted names
        sps_names <- sps_names[grep("accepted",sps_names$taxonomicStatus),]
    } 
    return(sps_names)
}

####################
find_names_from_id <- function(id, name, db = "gbif"){
    togosyno_accepted <- taxadb::filter_id(id, 
                                           provider = db, 
                                           type = "acceptedNameUsageID")
    name = data.frame(original_search = name, input = id)
    togosyno_accepted <- merge(togosyno_accepted, name, by = "input")
    
    if(any(duplicated(togosyno_accepted$acceptedNameUsageID))){
        # 1) get the accepted name
        dups <- unique(togosyno_accepted$acceptedNameUsageID[which(duplicated(togosyno_accepted$acceptedNameUsageID))])
        found <- lapply(dups, function(x){
            tmp <- togosyno_accepted[which(togosyno_accepted$acceptedNameUsageID == x),]
            tmp[which(tmp$taxonomicStatus=="accepted"),]
        })
        found <- do.call(rbind,found)
        
        # feed
        # 1) rem species found
        rem <- which(togosyno_accepted$acceptedNameUsageID %in% dups)
        if(any(rem)){
            togosyno_accepted <- togosyno_accepted[-rem,]
        }
        # 2) add species found
        togosyno_accepted <- rbind(togosyno_accepted,
                                   found)
        
    }
    
    if(any(duplicated(togosyno_accepted$acceptedNameUsageID))){
        # 2) use the most similar name
        dups <- unique(togosyno_accepted$acceptedNameUsageID[which(duplicated(togosyno_accepted$acceptedNameUsageID))])
        found <- lapply(dups, function(x){
            tmp <- togosyno_accepted[which(togosyno_accepted$acceptedNameUsageID == x),]
            tmp[order(adist(unique(tmp$original_search),tmp$scientificName))[1],]
        })
        found <- do.call(rbind,found)
        
        # feed
        # 1) rem species found
        rem <- which(togosyno_accepted$acceptedNameUsageID %in% dups)
        if(any(rem)){
            togosyno_accepted <- togosyno_accepted[-rem,]
        }
        # 2) add species found
        togosyno_accepted <- bind_rows(togosyno_accepted,
                                       found)
        
    }
    
    if(any(duplicated(togosyno_accepted$original_search))){
        # 3) get the one with more occurrence records
        dups <- unique(togosyno_accepted$original_search[which(duplicated(togosyno_accepted$original_search))])
        found <- data.frame()
        for(i in 1:length(dups)){
            dup_tmp <- togosyno_accepted[which(togosyno_accepted$original_search %in% dups[1]),]
            n_tmp <- sapply(dup_tmp$acceptedNameUsageID, function(x){
                rgbif::occ_count(taxonKey = readr::parse_number(x))
            })
            n_tmp <- names(sort(n_tmp, decreasing = TRUE)[1])
            dup_tmp <- dup_tmp %>% filter(acceptedNameUsageID %in% n_tmp)
            found <- rbind(found, dup_tmp)
        }
        
        # feed
        # 1) rem species found
        rem <- which(togosyno_accepted$original_search %in% dups)
        if(any(rem)){
            togosyno_accepted <- togosyno_accepted[-rem,]
        }
        # 2) add species found
        togosyno_accepted <- bind_rows(togosyno_accepted,
                                       found)
        
    }
    
    return(togosyno_accepted)
}

####################
query_names <- function(sci_name,
                        replace_synonyms = TRUE,
                        suggest_names = TRUE,
                        suggestion_distance = 0.9,
                        db = "gbif",
                        rank_name = NULL,
                        rank = NULL,
                        parallel = FALSE,
                        ncores = 2,
                        export_accepted = FALSE){
    
    ori_names <- sci_name # save original names
    if(db == 'gbif'){ # gbif can find anything with a '-'
        sci_name = gsub("-"," ",sci_name)
    }
    ori_names <- data.frame(original_search = ori_names, sci_name)
    
    sps_names<-bdc_query_names_taxadb(sci_name = ori_names$sci_name,
                                      replace_synonyms = replace_synonyms,
                                      suggest_names = suggest_names,
                                      suggestion_distance = suggestion_distance,
                                      db = db,
                                      parallel = parallel,
                                      ncores = ncores,
                                      export_accepted = export_accepted) 
    sps_names$scientificName <- paste(sps_names$genus,sps_names$specificEpithet)
    sps_names$scientificName[which(sps_names$scientificName=="NA NA")] <- NA
    if(db == 'gbif'){ # return original names to data
        sps_names$id <- sps_names$original_search
        sps_names <- sps_names[,-which(colnames(sps_names)=="original_search")]
        sps_names <- merge(sps_names, ori_names, by.x = "id",by.y = 'sci_name',all.x = T)
        sps_names <- data.frame(sps_names)
        sps_names <- sps_names[,-which(colnames(sps_names)=="id")]
    } 
    return(sps_names)
}

####################
# if(db == 'gbif'){
#     # select the one with the greatest N occurrence
#     ids <- gsub("GBIF:","",togosyno_accepted$acceptedNameUsageID)
#     nocc <- get_n_occur_gbif(gbif_code = ids)
#     dups_names = data.frame(original_search=togosyno_accepted$original_search,
#                             nocc)
#     tolook <- unique(dups_names$original_search)
#     found <- lapply(tolook, function(x){
#         tmp <- dups_names[which(dups_names == x),]
#         tmp[which(tmp$nocc == max(tmp$nocc)),]
#     })
#     togosyno_accepted <- do.call(rbind,found)
# }

####################
# x=togosyno$original_search
# db="itis"

# "Alchemilla arvensis"
# "Perezia microcephala"

find_names_filter <- function(x, db = db){
    togosyno_accepted <- taxadb::filter_name(x, 
                                             provider = db)
    togosyno_accepted$original_search <- togosyno_accepted$input
    if(any(is.na(togosyno_accepted$acceptedNameUsageID))){
        togosyno_accepted <- togosyno_accepted[-which(is.na(togosyno_accepted$acceptedNameUsageID)),]
    }
    return(togosyno_accepted)
}

