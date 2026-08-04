
# For genome info gathering function to search NCBI
# Testing Update which incorporates 1. GCA in assembly, BioProject codes, and searches assembly/genomes/bioproject deferring to assembly but returning others as partial data if all that is avaliable






conda_lib <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_lib)) {
  conda_r_lib <- file.path(conda_lib, "lib", "R", "library")
  .libPaths(c(conda_r_lib, .libPaths()))
}


# install packages if they are missing
if (!requireNamespace("dplyr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("dplyr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("getopt", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("getopt", repos = "https://cloud.r-project.org", quiet = TRUE)))
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("stringr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}

if (!requireNamespace("rentrez", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("rentrez", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("readr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("readr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
# data.table intentionally NOT loaded; optional via --use_data_table / BRIOCHE_USE_DT
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyverse", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("stringi", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("stringi", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("tidyselect", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyselect", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("magrittr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  suppressWarnings(suppressMessages(install.packages("tidyr", repos = "https://cloud.r-project.org", quiet = TRUE)))
}


suppressPackageStartupMessages({
  library(getopt)
  library(dplyr)
  library(stringr)
  library(rentrez)
  library(readr)
  # DO NOT: library(data.table)
})



args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--", a)) {
      a2 <- sub("^--", "", a)
      if (grepl("=", a2)) {
        kv <- strsplit(a2, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- kv[2]
      } else {
        if (i < length(args) && !grepl("^--", args[i + 1])) {
          out[[a2]] <- args[i + 1]
          i <- i + 1L
        } else {
          out[[a2]] <- "true" 
        }
      }
    }
    i <- i + 1L
  }
  out
}

opts <- parse_args(args)

get_opt <- function(opts, key, default = NULL) {
  if (!is.null(opts[[key]])) opts[[key]] else default
}


Accession <- get_opt(opts, "Accession", "")
outputfile <- get_opt(opts, "outputfile", "")
message(sprintf("Accession      : %s", Accession))
message(sprintf("outputfile      : %s", outputfile))









`%||%` <- function(a, b) {
  if (!is.null(a) &&
      length(a) > 0 &&
      !all(is.na(a)) &&
      any(nzchar(as.character(a)))) {
    a
  } else {
    b
  }
}



fetch_ncbi_meta <- function(accession = "", genomename = "",
                            search_genome_db = TRUE,
                            search_bioproject_db = TRUE,
                            retmax = 10) {
  out <- list(ok = FALSE)
  
  empty_to_na <- function(x) {
    x <- as.character(x %||% NA_character_)
    x <- trimws(x)
    
    x_norm <- tolower(x)
    x_norm <- gsub("[[:space:]]+", " ", x_norm)
    x_norm <- gsub("^[[:punct:]]+|[[:punct:]]+$", "", x_norm)
    
    missing_tokens <- c(
      "",
      "na",
      "n/a",
      "n.a",
      "nan",
      "null",
      "none",
      "nil",
      "missing",
      "unknown",
      "unk",
      "not known",
      "not applicable",
      "not available",
      "unavailable",
      "unspecified",
      "undefined",
      "undetermined",
      "unassigned",
      "not assigned",
      "no data",
      "no information",
      "no value"
    )
    
    x[x_norm %in% missing_tokens] <- NA_character_
    x
  }
  
  first_nonempty <- function(...) {
    vals <- list(...)
    for (v in vals) {
      v <- empty_to_na(v)
      if (length(v) > 0 && any(!is.na(v) & nzchar(v))) {
        return(v[which(!is.na(v) & nzchar(v))[1]])
      }
    }
    NA_character_
  }
  
  nz <- function(x) {
    x <- empty_to_na(x)
    length(x) == 1 && !is.na(x) && nzchar(x)
  }
  
  any_nz <- function(x) {
    x <- empty_to_na(x)
    any(!is.na(x) & nzchar(x))
  }
  
  strip_asm_version <- function(x) {
    sub("\\.[0-9]+$", "", toupper(trimws(x)))
  }
  
  is_assembly_accession <- function(x) {
    grepl("^(GCA|GCF)_", trimws(x), ignore.case = TRUE)
  }
  
  is_bioproject_accession <- function(x) {
    grepl("^PRJ[A-Z]{2,}[0-9]+$", trimws(x), ignore.case = TRUE)
  }
  
  esc_query <- function(x) {
    gsub('"', '\\"', x)
  }
  
  normalize_summaries <- function(sums) {
    if (is.null(sums)) {
      list()
    } else if (is.list(sums) && !is.null(sums$uid)) {
      list(sums)
    } else if (is.list(sums)) {
      sums
    } else {
      list(sums)
    }
  }
  
  safe_search <- function(db, term, retmax = 10) {
    tryCatch(
      rentrez::entrez_search(db = db, term = term, retmax = retmax),
      error = function(e) list(ids = character(0), count = 0)
    )
  }
  
  safe_summary <- function(db, ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    if (!length(ids)) return(list())
    
    tryCatch(
      normalize_summaries(rentrez::entrez_summary(db = db, id = ids)),
      error = function(e) list()
    )
  }
  
  get_asm_acc <- function(x) {
    empty_to_na(x[["assemblyaccession"]] %||% x[["accn"]] %||% NA_character_)
  }
  
  get_asm_name <- function(x) {
    empty_to_na(x[["assemblyname"]] %||% NA_character_)
  }
  
  get_ftp <- function(x, prefer = c("auto", "refseq", "genbank")) {
    prefer <- match.arg(prefer)
    
    ftp_ref <- empty_to_na(x[["ftppath_refseq"]] %||% NA_character_)
    ftp_gb  <- empty_to_na(x[["ftppath_genbank"]] %||% NA_character_)
    
    if (prefer == "refseq") {
      if (!is.na(ftp_ref)) ftp_ref else ftp_gb
    } else if (prefer == "genbank") {
      if (!is.na(ftp_gb)) ftp_gb else ftp_ref
    } else {
      acc <- get_asm_acc(x)
      
      if (!is.na(acc) && grepl("^GCA_", acc, ignore.case = TRUE)) {
        if (!is.na(ftp_gb)) ftp_gb else ftp_ref
      } else {
        if (!is.na(ftp_ref)) ftp_ref else ftp_gb
      }
    }
  }
  
  unique_or_null <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    if (length(ids) == 1) ids else NULL
  }
  
  choose_assembly_by_accession <- function(ids, acc_in) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    if (!length(ids)) return(NULL)
    
    sums <- safe_summary("assembly", ids)
    if (!length(sums)) return(NULL)
    
    acc_in_full  <- toupper(trimws(acc_in))
    acc_in_short <- strip_asm_version(acc_in)
    
    asm_accs <- vapply(sums, get_asm_acc, character(1))
    asm_accs <- empty_to_na(asm_accs)
    
    full_match <- !is.na(asm_accs) & toupper(asm_accs) == acc_in_full
    if (sum(full_match) == 1) {
      return(ids[which(full_match)])
    }
    
    short_match <- !is.na(asm_accs) & strip_asm_version(asm_accs) == acc_in_short
    if (sum(short_match) == 1) {
      return(ids[which(short_match)])
    }
    
    NULL
  }
  
  choose_assembly_by_name <- function(ids, name_in) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    if (!length(ids)) return(NULL)
    
    sums <- safe_summary("assembly", ids)
    if (!length(sums)) return(NULL)
    
    asm_names <- vapply(sums, get_asm_name, character(1))
    asm_names <- empty_to_na(asm_names)
    
    exact <- !is.na(asm_names) & tolower(asm_names) == tolower(trimws(name_in))
    
    if (sum(exact) == 1) {
      return(ids[which(exact)])
    }
    
    NULL
  }
  
  link_database_to_assembly <- function(ids, dbfrom) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    
    if (!length(ids)) return(character(0))
    
    linked <- character(0)
    
    for (id_i in ids) {
      lk <- tryCatch(
        rentrez::entrez_link(dbfrom = dbfrom, db = "assembly", id = id_i),
        error = function(e) NULL
      )
      
      if (!is.null(lk) && !is.null(lk$links)) {
        linked <- unique(c(
          linked,
          unlist(lk$links, use.names = FALSE)
        ))
      }
    }
    
    linked <- unique(as.character(linked))
    linked[!is.na(linked) & nzchar(linked)]
  }
  
  link_genome_to_assembly <- function(genome_ids) {
    link_database_to_assembly(genome_ids, dbfrom = "genome")
  }
  
  link_bioproject_to_assembly <- function(bioproject_ids) {
    link_database_to_assembly(bioproject_ids, dbfrom = "bioproject")
  }
  
  search_genome_ids <- function(query_text, mode = c("accession", "name")) {
    mode <- match.arg(mode)
    
    query_text <- trimws(query_text)
    if (!nzchar(query_text)) return(character(0))
    
    if (mode == "accession") {
      g_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('%s[All Fields]', query_text),
        if (grepl("\\.[0-9]+$", query_text)) {
          sprintf('"%s"[All Fields]', strip_asm_version(query_text))
        } else {
          character(0)
        }
      ))
    } else {
      g_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('"%s"[Title]', esc_query(query_text))
      ))
    }
    
    genome_ids <- character(0)
    
    for (q in g_queries) {
      gs <- safe_search("genome", q, retmax = retmax)
      genome_ids <- unique(c(genome_ids, gs$ids))
    }
    
    genome_ids <- unique(as.character(genome_ids))
    genome_ids[!is.na(genome_ids) & nzchar(genome_ids)]
  }
  
  search_bioproject_ids <- function(query_text, mode = c("accession", "name")) {
    mode <- match.arg(mode)
    
    query_text <- trimws(query_text)
    if (!nzchar(query_text)) return(character(0))
    
    if (mode == "accession") {
      bp_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('%s[All Fields]', query_text),
        sprintf('"%s"[Project Accession]', esc_query(query_text)),
        sprintf('"%s"[BioProject Accession]', esc_query(query_text))
      ))
    } else {
      bp_queries <- unique(c(
        sprintf('"%s"[All Fields]', esc_query(query_text)),
        sprintf('"%s"[Title]', esc_query(query_text)),
        sprintf('"%s"[Project Name]', esc_query(query_text))
      ))
    }
    
    bioproject_ids <- character(0)
    
    for (q in bp_queries) {
      bp <- safe_search("bioproject", q, retmax = retmax)
      bioproject_ids <- unique(c(bioproject_ids, bp$ids))
    }
    
    bioproject_ids <- unique(as.character(bioproject_ids))
    bioproject_ids[!is.na(bioproject_ids) & nzchar(bioproject_ids)]
  }
  
  search_genome_for_assemblies <- function(query_text, mode = c("accession", "name")) {
    genome_ids <- search_genome_ids(query_text, mode = match.arg(mode))
    link_genome_to_assembly(genome_ids)
  }
  
  search_bioproject_for_assemblies <- function(query_text, mode = c("accession", "name")) {
    bioproject_ids <- search_bioproject_ids(query_text, mode = match.arg(mode))
    link_bioproject_to_assembly(bioproject_ids)
  }
  
  extract_genome_partial <- function(genome_id, query_text = "") {
    sums <- safe_summary("genome", genome_id)
    if (!length(sums)) return(NULL)
    
    s <- sums[[1]]
    
    title <- first_nonempty(
      s[["title"]],
      s[["name"]],
      s[["genome"]],
      s[["description"]],
      query_text
    )
    
    organism <- first_nonempty(
      s[["organism"]],
      s[["organismname"]],
      s[["speciesname"]],
      s[["scientificname"]],
      s[["taxname"]]
    )
    
    taxid <- first_nonempty(
      s[["taxid"]],
      s[["taxonomyid"]]
    )
    
    acc <- first_nonempty(
      s[["accession"]],
      s[["assemblyaccession"]],
      s[["accn"]],
      query_text
    )
    
    url <- paste0("https://www.ncbi.nlm.nih.gov/genome/", genome_id)
    
    has_any <- any_nz(c(title, organism, taxid, acc))
    if (!has_any) return(NULL)
    
    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "genome",
      source_uid        = as.character(genome_id),
      source_accession  = acc,
      accession         = acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(genome_id),
      lookup_used       = "partial",
      lookup_via        = "genome",
      partial_reason    = "Genome record found but no unique linked Assembly record was resolved",
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }
  
  extract_bioproject_partial <- function(bioproject_id, query_text = "") {
    sums <- safe_summary("bioproject", bioproject_id)
    if (!length(sums)) return(NULL)
    
    s <- sums[[1]]
    
    bp_acc <- first_nonempty(
      s[["project_acc"]],
      s[["projectacc"]],
      s[["projectaccession"]],
      s[["bioprojectaccn"]],
      s[["accession"]],
      query_text
    )
    
    title <- first_nonempty(
      s[["project_name"]],
      s[["projectname"]],
      s[["title"]],
      s[["name"]],
      s[["description"]],
      query_text
    )
    
    organism <- first_nonempty(
      s[["organism"]],
      s[["organismname"]],
      s[["speciesname"]],
      s[["scientificname"]],
      s[["taxname"]]
    )
    
    taxid <- first_nonempty(
      s[["taxid"]],
      s[["taxonomyid"]]
    )
    
    url <- if (!is.na(bp_acc) && nzchar(bp_acc)) {
      paste0("https://www.ncbi.nlm.nih.gov/bioproject/", bp_acc)
    } else {
      paste0("https://www.ncbi.nlm.nih.gov/bioproject/", bioproject_id)
    }
    
    has_any <- any_nz(c(bp_acc, title, organism, taxid))
    if (!has_any) return(NULL)
    
    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "bioproject",
      source_uid        = as.character(bioproject_id),
      source_accession  = bp_acc,
      accession         = bp_acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(bioproject_id),
      lookup_used       = "partial",
      lookup_via        = "bioproject",
      partial_reason    = "BioProject record found but no unique linked Assembly record was resolved",
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }
  
  extract_assembly_summary_partial <- function(asm_sum, asm_uid, reason = "") {
    acc <- first_nonempty(
      asm_sum[["assemblyaccession"]],
      asm_sum[["accn"]]
    )
    
    title <- first_nonempty(
      asm_sum[["assemblyname"]],
      asm_sum[["title"]],
      asm_sum[["name"]]
    )
    
    organism <- first_nonempty(
      asm_sum[["organism"]],
      asm_sum[["speciesname"]],
      asm_sum[["scientificname"]],
      asm_sum[["taxname"]]
    )
    
    taxid <- first_nonempty(
      asm_sum[["taxid"]],
      asm_sum[["taxonomyid"]]
    )
    
    url <- if (!is.na(acc) && nzchar(acc)) {
      paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/", acc)
    } else {
      paste0("https://www.ncbi.nlm.nih.gov/assembly/", asm_uid)
    }
    
    has_any <- any_nz(c(acc, title, organism, taxid))
    if (!has_any) return(NULL)
    
    list(
      ok                = TRUE,
      partial           = TRUE,
      assembly_resolved = FALSE,
      source_db         = "assembly",
      source_uid        = as.character(asm_uid),
      source_accession  = acc,
      accession         = acc,
      asm_name          = title,
      organism          = organism,
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = url,
      genome_url        = url,
      ftp_path          = NA_character_,
      ncbi_uid          = as.character(asm_uid),
      lookup_used       = "partial",
      lookup_via        = "assembly",
      partial_reason    = reason,
      chr_table         = data.frame(),
      len_by_seqname    = setNames(numeric(0), character(0)),
      md5_by_seqname    = setNames(character(0), character(0))
    )
  }
  
  find_partial_meta <- function(acc_in, name_in) {
    if (nzchar(acc_in) && isTRUE(search_bioproject_db) && is_bioproject_accession(acc_in)) {
      bp_ids <- search_bioproject_ids(acc_in, mode = "accession")
      bp_id <- unique_or_null(bp_ids)
      
      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }
    
    if (nzchar(acc_in) && isTRUE(search_genome_db)) {
      genome_ids <- search_genome_ids(acc_in, mode = "accession")
      genome_id <- unique_or_null(genome_ids)
      
      if (!is.null(genome_id)) {
        partial <- extract_genome_partial(genome_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }
    
    if (nzchar(acc_in) && isTRUE(search_bioproject_db)) {
      bp_ids <- search_bioproject_ids(acc_in, mode = "accession")
      bp_id <- unique_or_null(bp_ids)
      
      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = acc_in)
        if (!is.null(partial)) return(partial)
      }
    }
    
    if (nzchar(name_in) && isTRUE(search_genome_db)) {
      genome_ids <- search_genome_ids(name_in, mode = "name")
      genome_id <- unique_or_null(genome_ids)
      
      if (!is.null(genome_id)) {
        partial <- extract_genome_partial(genome_id, query_text = name_in)
        if (!is.null(partial)) return(partial)
      }
    }
    
    if (nzchar(name_in) && isTRUE(search_bioproject_db)) {
      bp_ids <- search_bioproject_ids(name_in, mode = "name")
      bp_id <- unique_or_null(bp_ids)
      
      if (!is.null(bp_id)) {
        partial <- extract_bioproject_partial(bp_id, query_text = name_in)
        if (!is.null(partial)) return(partial)
      }
    }
    
    NULL
  }
  
  find_assembly_id <- function(acc_in, name_in) {
    if (nzchar(acc_in)) {
      acc_in <- trimws(acc_in)
      
      queries <- unique(c(
        sprintf('%s[Assembly Accession]', acc_in),
        sprintf('"%s"[Assembly Accession]', esc_query(acc_in)),
        sprintf('"%s"[All Fields]', esc_query(acc_in)),
        if (grepl("\\.[0-9]+$", acc_in)) {
          sprintf('"%s"[All Fields]', strip_asm_version(acc_in))
        } else {
          character(0)
        }
      ))
      
      ids <- character(0)
      
      for (q in queries) {
        s <- safe_search("assembly", q, retmax = retmax)
        ids <- unique(c(ids, s$ids))
      }
      
      if (is_assembly_accession(acc_in)) {
        pick_id <- choose_assembly_by_accession(ids, acc_in)
      } else {
        pick_id <- unique_or_null(ids)
      }
      
      if (!is.null(pick_id)) {
        return(list(id = pick_id, used = "acc", via = "assembly"))
      }
      
      if (isTRUE(search_genome_db)) {
        asm_ids <- search_genome_for_assemblies(acc_in, mode = "accession")
        
        if (is_assembly_accession(acc_in)) {
          pick_id <- choose_assembly_by_accession(asm_ids, acc_in)
        } else {
          pick_id <- unique_or_null(asm_ids)
        }
        
        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "acc", via = "genome_link"))
        }
      }
      
      if (isTRUE(search_bioproject_db)) {
        asm_ids <- search_bioproject_for_assemblies(acc_in, mode = "accession")
        
        if (is_assembly_accession(acc_in)) {
          pick_id <- choose_assembly_by_accession(asm_ids, acc_in)
        } else {
          pick_id <- unique_or_null(asm_ids)
        }
        
        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "acc", via = "bioproject_link"))
        }
      }
    }
    
    if (nzchar(name_in)) {
      name_in <- trimws(name_in)
      
      queries <- unique(c(
        sprintf('"%s"[Assembly Name]', esc_query(name_in)),
        sprintf('"%s"[All Fields]', esc_query(name_in))
      ))
      
      ids <- character(0)
      
      for (q in queries) {
        s <- safe_search("assembly", q, retmax = retmax)
        ids <- unique(c(ids, s$ids))
      }
      
      pick_id <- choose_assembly_by_name(ids, name_in)
      
      if (!is.null(pick_id)) {
        return(list(id = pick_id, used = "name", via = "assembly"))
      }
      
      if (isTRUE(search_genome_db)) {
        asm_ids <- search_genome_for_assemblies(name_in, mode = "name")
        pick_id <- choose_assembly_by_name(asm_ids, name_in)
        
        if (is.null(pick_id)) {
          pick_id <- unique_or_null(asm_ids)
        }
        
        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "name", via = "genome_link"))
        }
      }
      
      if (isTRUE(search_bioproject_db)) {
        asm_ids <- search_bioproject_for_assemblies(name_in, mode = "name")
        pick_id <- choose_assembly_by_name(asm_ids, name_in)
        
        if (is.null(pick_id)) {
          pick_id <- unique_or_null(asm_ids)
        }
        
        if (!is.null(pick_id)) {
          return(list(id = pick_id, used = "name", via = "bioproject_link"))
        }
      }
    }
    
    list(id = NULL, used = "", via = "")
  }
  
  tryCatch({
    acc_in  <- trimws(accession %||% "")
    name_in <- trimws(genomename %||% "")
    
    found <- find_assembly_id(acc_in, name_in)
    
    if (is.null(found$id)) {
      partial <- find_partial_meta(acc_in, name_in)
      
      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("No unique Assembly record resolved; returning partial metadata from ",
                partial$source_db)
        return(partial)
      }
      
      message("Accession or genome name were nonspecific unable to populate genomic information")
      return(out)
    }
    
    sum <- rentrez::entrez_summary(db = "assembly", id = found$id)
    
    if (found$used == "acc") {
      acc_sum <- get_asm_acc(sum)
      
      if (is_assembly_accession(acc_in)) {
        if (is.na(acc_sum) ||
            strip_asm_version(acc_sum) != strip_asm_version(acc_in)) {
          message("Accession or genome name were nonspecific unable to populate genomic information")
          return(out)
        }
      }
    } else if (found$used == "name") {
      aname <- get_asm_name(sum)
      
      if (found$via == "assembly") {
        if (is.na(aname) || tolower(aname) != tolower(name_in)) {
          message("Accession or genome name were nonspecific unable to populate genomic information")
          return(out)
        }
      }
    }
    
    acc_sum <- get_asm_acc(sum)
    aname   <- get_asm_name(sum)
    org     <- empty_to_na(sum[["organism"]] %||% sum[["speciesname"]] %||% NA_character_)
    taxid   <- sum[["taxid"]] %||% NA
    
    ftp <- get_ftp(sum, prefer = "auto")
    
    if (is.na(acc_sum) || is.na(ftp)) {
      partial <- extract_assembly_summary_partial(
        sum,
        found$id,
        reason = "Assembly summary record found but FTP path or assembly report was unavailable"
      )
      
      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("Assembly summary resolved but assembly report/FTP path was unavailable; returning partial metadata from assembly")
        return(partial)
      }
      
      partial <- find_partial_meta(acc_in, name_in)
      
      if (!is.null(partial) && isTRUE(partial$ok)) {
        message("Assembly record resolved but assembly report/FTP path was unavailable; returning partial metadata from ",
                partial$source_db)
        return(partial)
      }
      
      message("Accession or genome name were nonspecific unable to populate genomic information")
      return(out)
    }
    
    ftp_https <- sub("^ftp://", "https://", ftp)
    base      <- basename(ftp_https)
    rpt_url   <- paste0(ftp_https, "/", base, "_assembly_report.txt")
    asm_url   <- paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/", acc_sum)
    
    parse_assembly_report <- function(rpt_url) {
      lines <- readr::read_lines(rpt_url, progress = FALSE)
      
      hdr_i <- which(grepl("^#\\s*Sequence-Name\\t", lines))[1]
      
      if (is.na(hdr_i)) {
        stop("assembly_report header line not found")
      }
      
      header <- sub("^#\\s*", "", lines[hdr_i])
      data   <- lines[(hdr_i + 1L):length(lines)]
      data   <- data[!grepl("^#", data)]
      
      txt <- paste(c(header, data), collapse = "\n")
      
      readr::read_tsv(
        I(txt),
        show_col_types = FALSE,
        progress = FALSE,
        col_types = readr::cols(.default = "c")
      )
    }
    
    ar <- parse_assembly_report(rpt_url)
    
    col_len      <- intersect(c("Sequence-Length", "sequence-length", "length"), names(ar))
    col_assigned <- intersect(c("Assigned-Molecule", "assigned-molecule"), names(ar))
    col_md5      <- intersect(c("MD5", "md5"), names(ar))
    col_seqname  <- intersect(c("Sequence-Name", "sequence-name"), names(ar))
    col_role     <- intersect(c("Sequence-Role", "sequence-role"), names(ar))
    col_refseq   <- intersect(c("RefSeq-Accn", "refseq-accn", "RefSeq Accn"), names(ar))
    col_genbank  <- intersect(c("GenBank-Accn", "genbank-accn", "GenBank Accn"), names(ar))
    
    if (!length(col_len)) {
      stop("assembly_report lacks Sequence-Length column")
    }
    
    key_col <- if (length(col_seqname)) {
      col_seqname[1]
    } else if (length(col_assigned)) {
      col_assigned[1]
    } else {
      NULL
    }
    
    if (is.null(key_col)) {
      stop("assembly_report lacks Sequence-Name / Assigned-Molecule columns")
    }
    
    seq_name <- empty_to_na(ar[[key_col]])
    
    assigned <- if (length(col_assigned)) {
      empty_to_na(ar[[col_assigned[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }
    
    len_val <- suppressWarnings(as.numeric(ar[[col_len[1]]]))
    
    md5_val <- if (length(col_md5)) {
      empty_to_na(ar[[col_md5[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }
    
    role_val <- if (length(col_role)) {
      empty_to_na(ar[[col_role[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }
    
    ref_acc <- if (length(col_refseq)) {
      empty_to_na(ar[[col_refseq[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }
    
    gb_acc <- if (length(col_genbank)) {
      empty_to_na(ar[[col_genbank[1]]])
    } else {
      rep(NA_character_, nrow(ar))
    }
    
    extract_num <- function(x) {
      s <- tolower(trimws(x))
      
      m <- stringr::str_match(
        s,
        "^.*(?:chromosome|chrom|chr|contig|scaffold|segment|seg)?[ _-]*0*([0-9]+)([a-z])?(?![a-z0-9]).*"
      )
      
      suppressWarnings(as.integer(m[, 2]))
    }
    
    chrom_num <- extract_num(seq_name)
    
    chrom_num <- ifelse(
      is.na(chrom_num),
      extract_num(assigned),
      chrom_num
    )
    
    preferred <- if (grepl("^GCA_", acc_sum, ignore.case = TRUE)) {
      ifelse(!is.na(gb_acc) & nzchar(gb_acc), gb_acc, ref_acc)
    } else {
      ifelse(!is.na(ref_acc) & nzchar(ref_acc), ref_acc, gb_acc)
    }
    
    nuccore <- ifelse(
      !is.na(preferred) & nzchar(preferred),
      paste0("https://www.ncbi.nlm.nih.gov/nuccore/", preferred),
      NA_character_
    )
    
    keep <- !is.na(len_val) & !is.na(seq_name) & nzchar(seq_name)
    
    chr_table <- data.frame(
      chrom_num         = chrom_num,
      seq_name          = seq_name,
      assigned_molecule = assigned,
      sequence_role     = role_val,
      refseq_accn       = ref_acc,
      genbank_accn      = gb_acc,
      preferred_accn    = preferred,
      length            = len_val,
      md5               = md5_val,
      nuccore_url       = nuccore,
      stringsAsFactors  = FALSE
    )[keep, , drop = FALSE]
    
    key_seq <- tolower(chr_table$seq_name)
    
    len_by_seqname <- stats::setNames(chr_table$length, key_seq)
    md5_by_seqname <- stats::setNames(chr_table$md5, key_seq)
    
    list(
      ok                = TRUE,
      partial           = FALSE,
      assembly_resolved = TRUE,
      accession         = acc_sum,
      asm_name          = aname %||% "",
      organism          = org %||% "",
      taxid             = suppressWarnings(as.integer(taxid)),
      asm_url           = asm_url,
      genome_url        = asm_url,
      ftp_path          = ftp_https,
      ncbi_uid          = found$id,
      lookup_used       = found$used,
      lookup_via        = found$via,
      source_db         = "assembly",
      source_uid        = found$id,
      source_accession  = acc_sum,
      partial_reason    = NA_character_,
      chr_table         = chr_table,
      len_by_seqname    = len_by_seqname,
      md5_by_seqname    = md5_by_seqname
    )
  }, error = function(e) {
    partial <- find_partial_meta(
      trimws(accession %||% ""),
      trimws(genomename %||% "")
    )
    
    if (!is.null(partial) && isTRUE(partial$ok)) {
      message("Assembly metadata extraction failed; returning partial metadata from ",
              partial$source_db)
      return(partial)
    }
    
    message("Accession or genome name were nonspecific unable to populate genomic information")
    out
  })
}


build_extra_meta <- function(ncbi_meta, genomename,
                             override_species = NULL,
                             override_url = NULL) {
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(trimws(x))] <- "unknown"
    gsub('"', '\\\\"', x, fixed = TRUE)
  }
  
  clean <- function(x, missing = "unknown") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(missing)
    }
    
    x <- as.character(x[1])
    x <- trimws(x)
    
    if (is.na(x) || !nzchar(x)) {
      missing
    } else {
      x
    }
  }
  
  clean_taxid <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return("unknown")
    }
    
    tx <- suppressWarnings(as.integer(x[1]))
    
    if (is.na(tx)) {
      "unknown"
    } else {
      as.character(tx)
    }
  }
  
  clean_url <- function(...) {
    vals <- list(...)
    
    for (v in vals) {
      u <- clean(v, missing = "")
      if (nzchar(u)) {
        return(u)
      }
    }
    
    "unknown"
  }
  
  gn <- clean(genomename, missing = "unknown")
  
  is_ok <- isTRUE(ncbi_meta$ok)
  
  is_partial <- is_ok && (
    isTRUE(ncbi_meta$partial) ||
      identical(ncbi_meta$assembly_resolved, FALSE)
  )
  
  species_fallback <- function() {
    sp <- clean(override_species, missing = "")
    
    if (!nzchar(sp)) {
      sp <- clean(gn, missing = "")
    }
    
    if (!nzchar(sp)) {
      sp <- "unknown"
    }
    
    sp
  }
  
  if (is_ok) {
    sp <- clean(ncbi_meta$organism, missing = "")
    
    if (!nzchar(sp)) {
      sp <- species_fallback()
    }
    
    acc <- clean(ncbi_meta$accession, missing = "")
    
    if (!nzchar(acc)) {
      acc <- clean(ncbi_meta$source_accession, missing = "")
    }
    
    if (!nzchar(acc)) {
      acc <- gn
    }
    
    if (!nzchar(acc) || is.na(acc)) {
      acc <- "unknown"
    }
    
    asm_name <- clean(ncbi_meta$asm_name, missing = "")
    
    if (!nzchar(asm_name)) {
      asm_name <- gn
    }
    
    if (!nzchar(asm_name) || is.na(asm_name)) {
      asm_name <- "unknown"
    }
    
    taxid <- clean_taxid(ncbi_meta$taxid)
    
    meta_url <- clean_url(
      ncbi_meta$asm_url,
      ncbi_meta$genome_url,
      override_url
    )
    
    source_db <- clean(ncbi_meta$source_db, missing = "")
    
    if (!nzchar(source_db)) {
      source_db <- clean(ncbi_meta$lookup_via, missing = "")
    }
    
    if (!nzchar(source_db)) {
      source_db <- "unknown"
    }
    
    automated_lookup <- if (is_partial) {
      "partial metadata found only"
    } else {
      "resolved by automated process"
    }
    
    c(
      sprintf("##reference=%s", acc),
      sprintf(
        '##assembly=<accession=%s,name="%s",species="%s",taxonomy=%s,url="%s",source="%s",automated_lookup="%s">',
        acc,
        esc(asm_name),
        esc(sp),
        taxid,
        esc(meta_url),
        esc(source_db),
        esc(automated_lookup)
      ),
      sprintf("##genome_url=%s", meta_url)
    )
  } else {
    sp <- species_fallback()
    
    url_out <- clean_url(override_url)
    
    ref_val <- if (!is.na(gn) && nzchar(gn) && gn != "unknown") {
      gn
    } else if (!is.na(sp) && nzchar(sp) && sp != "unknown") {
      sp
    } else {
      "unknown"
    }
    
    c(
      sprintf("##reference=%s", ref_val),
      sprintf(
        '##assembly=<accession=%s,name="%s",species="%s",taxonomy=%s,url="%s",source="%s",automated_lookup="%s">',
        "unknown",
        esc(gn),
        esc(sp),
        "unknown",
        esc(url_out),
        "unknown",
        "not resolved by automated process"
      ),
      sprintf("##genome_url=%s", url_out)
    )
  }
}


build_contig_lines <- function(ids, meta, override_species = NULL) {
  clean <- function(x, missing = "unknown") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(missing)
    }
    
    x <- as.character(x[1])
    x <- trimws(x)
    
    if (is.na(x) || !nzchar(x)) {
      missing
    } else {
      x
    }
  }
  
  clean_num <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(NA_real_)
    }
    
    y <- suppressWarnings(as.numeric(x[1]))
    
    if (is.na(y)) {
      NA_real_
    } else {
      y
    }
  }
  
  clean_taxid <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return("unknown")
    }
    
    tx <- suppressWarnings(as.integer(x[1]))
    
    if (is.na(tx)) {
      "unknown"
    } else {
      as.character(tx)
    }
  }
  
  has_value <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(FALSE)
    }
    
    x <- as.character(x[1])
    !is.na(x) && nzchar(trimws(x))
  }
  
  esc <- function(x) {
    x <- clean(x, missing = "unknown")
    gsub('"', '\\\\"', x, fixed = TRUE)
  }
  
  strip_seq_version <- function(x) {
    sub("\\.[0-9]+$", "", trimws(as.character(x)))
  }
  
  make_keys <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    
    unique(c(
      tolower(x),
      tolower(strip_seq_version(x))
    ))
  }
  
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  
  if (!length(ids)) {
    return(character(0))
  }
  
  contig_species <- function() {
    sp <- clean(meta$organism, missing = "")
    
    if (!nzchar(sp)) {
      sp <- clean(override_species, missing = "")
    }
    
    if (!nzchar(sp)) {
      sp <- "unknown"
    }
    
    sp
  }
  
  tax <- clean_taxid(meta$taxid)
  
  mk_attr <- function(id,
                      len = NA_real_,
                      acc = NA_character_,
                      alias = NA_character_,
                      tax = "unknown") {
    id    <- clean(id, missing = "unknown")
    acc   <- clean(acc, missing = "unknown")
    alias <- clean(alias, missing = "unknown")
    
    parts <- c(sprintf("ID=%s", id))
    
    # Important:
    # Do NOT write length=unknown.
    # VCF/BCF parsers expect contig length, if present, to be numeric.
    if (!is.na(len)) {
      parts <- c(parts, sprintf("length=%d", as.integer(len)))
    }
    
    parts <- c(parts, sprintf('species="%s"', esc(contig_species())))
    parts <- c(parts, sprintf("taxonomy=%s", clean(tax, missing = "unknown")))
    parts <- c(parts, sprintf("Accession=%s", acc))
    parts <- c(parts, sprintf('Alias="%s"', esc(alias)))
    
    sprintf("##contig=<%s>", paste(parts, collapse = ","))
  }
  
  has_resolved_assembly <- isTRUE(meta$ok) &&
    !isTRUE(meta$partial) &&
    !identical(meta$assembly_resolved, FALSE)
  
  has_chr_table <- has_resolved_assembly &&
    !is.null(meta$chr_table) &&
    NROW(meta$chr_table) > 0
  
  if (!has_chr_table) {
    return(vapply(
      ids,
      function(id) {
        mk_attr(
          id    = id,
          len   = NA_real_,
          acc   = NA_character_,
          alias = NA_character_,
          tax   = tax
        )
      },
      FUN.VALUE = character(1)
    ))
  }
  
  ct <- meta$chr_table
  
  get_col <- function(df, col, default = NA_character_) {
    if (col %in% names(df)) {
      df[[col]]
    } else {
      rep(default, nrow(df))
    }
  }
  
  seq_col       <- get_col(ct, "seq_name")
  assn_col      <- get_col(ct, "assigned_molecule")
  len_col       <- suppressWarnings(as.numeric(get_col(ct, "length", NA_real_)))
  preferred_col <- get_col(ct, "preferred_accn")
  genbank_col   <- get_col(ct, "genbank_accn")
  refseq_col    <- get_col(ct, "refseq_accn")
  
  
  seq_chr  <- trimws(as.character(seq_col))
  assn_chr <- trimws(as.character(assn_col))
  
  has_seq  <- !is.na(seq_chr)  & nzchar(seq_chr)
  has_assn <- !is.na(assn_chr) & nzchar(assn_chr)
  
  alias_col <- ifelse(
    has_seq & has_assn,
    paste0(seq_chr, ";", assn_chr),
    ifelse(
      has_assn,
      assn_chr,
      seq_chr
    )
  )
  
  accession_col <- ifelse(
    !is.na(preferred_col) & nzchar(trimws(as.character(preferred_col))),
    as.character(preferred_col),
    ifelse(
      !is.na(genbank_col) & nzchar(trimws(as.character(genbank_col))),
      as.character(genbank_col),
      as.character(refseq_col)
    )
  )
  
  lookup <- new.env(parent = emptyenv())
  
  add_lookup <- function(keys, row_i) {
    keys <- make_keys(keys)
    
    for (k in keys) {
      if (!exists(k, envir = lookup, inherits = FALSE)) {
        assign(k, row_i, envir = lookup)
      }
    }
    
    invisible(NULL)
  }
  
  for (i in seq_len(nrow(ct))) {
    add_lookup(seq_col[i], i)
    add_lookup(assn_col[i], i)
    add_lookup(preferred_col[i], i)
    add_lookup(genbank_col[i], i)
    add_lookup(refseq_col[i], i)
  }
  
  vapply(
    ids,
    function(id) {
      id_keys <- make_keys(id)
      
      matched_row <- NA_integer_
      
      for (k in id_keys) {
        if (exists(k, envir = lookup, inherits = FALSE)) {
          matched_row <- get(k, envir = lookup, inherits = FALSE)
          break
        }
      }
      
      if (!is.na(matched_row)) {
        len <- clean_num(len_col[matched_row])
        
        acc <- if (has_value(accession_col[matched_row])) {
          accession_col[matched_row]
        } else {
          NA_character_
        }
        
        alias <- if (has_value(alias_col[matched_row])) {
          alias_col[matched_row]
        } else {
          NA_character_
        }
        
        mk_attr(
          id    = id,
          len   = len,
          acc   = acc,
          alias = alias,
          tax   = tax
        )
      } else {
        mk_attr(
          id    = id,
          len   = NA_real_,
          acc   = NA_character_,
          alias = NA_character_,
          tax   = tax
        )
      }
    },
    FUN.VALUE = character(1)
  )
}




user_species=""
genomename=""
user_genome_url=""

ncbi_meta <- list(ok = FALSE)
if (nzchar(Accession) || nzchar(genomename)) {
  ncbi_meta <- fetch_ncbi_meta(Accession, genomename)
  if (isTRUE(ncbi_meta$ok)) {
    message(sprintf("NCBI assembly  : %s (%s) taxid=%s",
                    ncbi_meta$accession, ncbi_meta$asm_name, ncbi_meta$taxid))
    message(sprintf("NCBI URL       : %s", ncbi_meta$asm_url))
  } else {
    message("NCBI metadata   : none found (proceeding with user/unknown fallbacks)")
  }
}

# if NCBI failed and user gave species
if (!isTRUE(ncbi_meta$ok) && nzchar(user_species)) {
  ncbi_meta$organism <- user_species
}

# Header meta lines (uses NCBI URL when present; otherwise uses user_genome_url)
extra_meta <- build_extra_meta(
  ncbi_meta,
  genomename,
  override_species = if (nzchar(user_species)) user_species else NULL,
  override_url     = if (nzchar(user_genome_url)) user_genome_url else NULL
)


VCFchroms <- ncbi_meta$chr_table$genbank_accn

contig_lines <- build_contig_lines(
  ids   = VCFchroms,
  meta  = ncbi_meta,
  override_species = if (nzchar(user_species)) user_species else NULL
)


savemetaenome <- extra_meta[2:3]


metadataallref <- c(savemetaenome, contig_lines)

write.table(metadataallref,outputfile,col.names = FALSE, row.names = FALSE,quote = FALSE)



