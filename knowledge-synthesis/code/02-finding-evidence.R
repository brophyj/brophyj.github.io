# R code for chapter: 02-finding-evidence.qmd
# From 'Knowledge Synthesis in Health Research' (J.M. Brophy).
# Run from the book's folder so that data/ is found.
# Lines starting with '## ' are not run by default (e.g. live database queries).

## -----------------------------------------------------------------------------
#| label: setup-ch02
#| include: false
# Shared setup, sourced at the top of every chapter --------------------------
knitr::opts_chunk$set(
  echo      = TRUE,      # code is part of the text: show it
  message   = FALSE,
  warning   = FALSE,
  fig.align = "center",
  fig.retina = 2,
  out.width = "85%"
)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
})
theme_set(theme_minimal(base_size = 13))
options(dplyr.summarise.inform = FALSE)
suppressPackageStartupMessages({
  library(xml2)
  library(stringr)
  library(kableExtra)
})


## flowchart TD
##   A["Records from databases<br/>CENTRAL = 841<br/>MEDLINE = 1923<br/>Embase = 2417"] --> B
##   A2["Other sources<br/>Registers = 142<br/>Grey literature = 37"] --> B
##   B["After deduplication<br/>n = 3961"] --> C
##   C["Titles and abstracts screened<br/>n = 3961"] --> D
##   C --> E["Excluded<br/>n = 3541"]
##   D["Full texts assessed<br/>n = 420"] --> F
##   D --> G["Excluded, with reasons<br/>Wrong population: 189<br/>Wrong intervention: 97<br/>Wrong outcome: 62<br/>Duplicate study: 30<br/>n = 378"]
##   F["Studies included in the review<br/>n = 42"] --> H
##   H["Studies included in the meta-analysis<br/>n = 38"]

## -----------------------------------------------------------------------------
#| label: prisma2020
#| eval: false
# library(PRISMA2020)
# 
# # Start from the package's template and overwrite the counts
# csv_file   <- system.file("extdata", "PRISMA.csv", package = "PRISMA2020")
# flow_input <- read.csv(csv_file)
# set_n <- function(df, field, n) { df$n[df$data == field] <- n; df }
# 
# flow_input <- flow_input |>
#   set_n("database_results",         4823) |>
#   set_n("register_results",          142) |>
#   set_n("duplicates",               1004) |>
#   set_n("records_screened",         3961) |>
#   set_n("records_excluded",         3541) |>
#   set_n("dbr_sought_reports",        420) |>
#   set_n("dbr_notretrieved_reports",    0) |>
#   set_n("dbr_assessed",              420) |>
#   set_n("dbr_excluded",              378) |>
#   set_n("new_studies",                42) |>
#   set_n("new_reports",                42) |>
#   set_n("total_studies",              42) |>
#   set_n("total_reports",              42) |>
#   set_n("total_studies_ma",           38) |>
#   set_n("total_reports_ma",           38)
# 
# flow_data <- PRISMA_data(flow_input)
# PRISMA_flowdiagram(flow_data,
#                    interactive   = TRUE,
#                    previous      = FALSE,   # not an update of an earlier review
#                    other         = FALSE,   # no separate "other methods" arm
#                    meta_analysis = TRUE)    # show the meta-analysis box


## -----------------------------------------------------------------------------
#| label: kappa-hand
reviewer1 <- c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1)
reviewer2 <- c(1, 0, 1, 0, 0, 1, 0, 1, 1, 1)

p_o <- mean(reviewer1 == reviewer2)                  # observed agreement
p_e <- mean(reviewer1) * mean(reviewer2) +           # both include by chance
       (1 - mean(reviewer1)) * (1 - mean(reviewer2)) # both exclude by chance
kappa <- (p_o - p_e) / (1 - p_e)
round(c(observed = p_o, chance = p_e, kappa = kappa), 2)


## -----------------------------------------------------------------------------
#| label: kappa-irr
#| eval: !expr requireNamespace("irr", quietly = TRUE)
irr::kappa2(cbind(reviewer1, reviewer2))


## -----------------------------------------------------------------------------
#| label: rentrez-steps
#| eval: false
# library(rentrez)
# library(xml2)
# 
# query <- '("heart failure"[MeSH Terms] OR "atrial fibrillation"[MeSH Terms])
#           AND "randomized controlled trial"[Publication Type]
#           AND "humans"[MeSH Terms]
#           AND (2022/01/01:2024/12/31[PDAT])'
# 
# # 1. Count the matching records (retmax = 0 returns the count only)
# entrez_search(db = "pubmed", term = query, retmax = 0)$count
# 
# # 2. Retrieve the PMIDs
# pmids <- entrez_search(db = "pubmed", term = query,
#                        retmax = 500, sort = "pub date")$ids
# 
# # 3. Fetch the records as XML (at most 200 per request)
# xml_raw <- entrez_fetch(db = "pubmed", id = pmids[1:50],
#                         rettype = "xml", parsed = FALSE)
# records <- xml_find_all(read_xml(xml_raw), ".//PubmedArticle")
# length(records)


## -----------------------------------------------------------------------------
#| label: parse-functions
# Extract the text of the first node matching an XPath, or NA if absent
extract_text <- function(node, xpath) {
  result <- xml_find_first(node, xpath)
  if (inherits(result, "xml_missing")) return(NA_character_)
  xml_text(result)
}

# One PubMed record -> one row
parse_article <- function(node) {
  year <- extract_text(node, ".//PubDate/Year")
  if (is.na(year)) year <- extract_text(node, ".//PubDate/MedlineDate")
  data.frame(
    PMID          = extract_text(node, ".//PMID"),
    Title         = extract_text(node, ".//ArticleTitle"),
    Journal       = extract_text(node, ".//Journal/Title"),
    JournalAbbrev = extract_text(node, ".//Journal/ISOAbbreviation"),
    Year          = year,
    Abstract      = paste(xml_text(xml_find_all(node, ".//Abstract/AbstractText")),
                          collapse = " "),
    stringsAsFactors = FALSE
  )
}


## -----------------------------------------------------------------------------
#| label: fetch-rcts
fetch_rcts <- function(mesh_terms, start, end, journals = NULL,
                       retmax = 5000, out_csv = NULL) {
  # 1. Build the query
  pop_clause <- paste0("(", paste(sprintf('"%s"[MeSH Terms]', mesh_terms),
                                  collapse = " OR "), ")")
  journal_clause <- if (is.null(journals)) "" else
    paste0(" AND (", paste(sprintf('"%s"[TA]', journals), collapse = " OR "), ")")
  query <- paste0(
    pop_clause, journal_clause,
    ' AND "randomized controlled trial"[Publication Type]',
    ' AND "humans"[MeSH Terms]',
    sprintf(" AND (%s:%s[PDAT])", start, end),
    ' NOT ("case reports"[Publication Type] OR "letter"[Publication Type])'
  )

  # 2. Retrieve the PMIDs
  n_total <- rentrez::entrez_search(db = "pubmed", term = query, retmax = 0)$count
  ids <- rentrez::entrez_search(db = "pubmed", term = query,
                                retmax = min(n_total, retmax),
                                sort = "pub date")$ids

  # 3. Fetch and parse in batches of 200
  batches <- split(ids, ceiling(seq_along(ids) / 200))
  results <- lapply(batches, function(batch) {
    xml_raw <- rentrez::entrez_fetch(db = "pubmed", id = batch,
                                     rettype = "xml", parsed = FALSE)
    Sys.sleep(0.34)                       # at most 3 requests per second
    recs <- xml_find_all(read_xml(xml_raw), ".//PubmedArticle")
    do.call(rbind, lapply(recs, parse_article))
  })
  df <- do.call(rbind, results)

  # 4. Save a dated copy
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    write.csv(df, out_csv, row.names = FALSE)
  }
  df
}


## -----------------------------------------------------------------------------
#| label: run-search
if (!file.exists("data/cv_rcts_demo.csv")) {
  cv_rcts <- fetch_rcts(
    mesh_terms = c("heart failure", "atrial fibrillation",
                   "coronary artery disease", "hypertension", "stroke"),
    journals   = c("N Engl J Med", "Lancet", "JAMA", "BMJ",
                   "JAMA Cardiol", "Circulation"),
    start      = "2023/01/01",
    end        = "2024/12/31",
    out_csv    = "data/cv_rcts_demo.csv")
}
cv_rcts <- read.csv("data/cv_rcts_demo.csv", stringsAsFactors = FALSE)
nrow(cv_rcts)


## -----------------------------------------------------------------------------
#| label: tbl-cv-rcts
#| tbl-cap: "The first four records retrieved."
cv_rcts |>
  select(PMID, Title, JournalAbbrev, Year) |>
  head(4) |>
  mutate(Title = str_trunc(Title, 70)) |>
  kable(col.names = c("PMID", "Title", "Journal", "Year"))


## -----------------------------------------------------------------------------
#| label: fig-cv-rcts
#| fig-cap: "Randomised trials retrieved, by journal and year."
#| fig-width: 7
#| fig-height: 4
cv_rcts |>
  count(JournalAbbrev, Year) |>
  ggplot(aes(x = factor(Year), y = n, fill = JournalAbbrev)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Year", y = "Number of trials", fill = NULL) +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: epmc
if (!file.exists("data/epmc_hf_trials.csv")) {
  epmc_hf <- europepmc::epmc_search(
    query  = 'MESH:"Heart Failure" AND PUB_TYPE:"Clinical Trial" AND FIRST_PDATE:[2022 TO 2024]',
    limit  = 500,
    output = "parsed")
  write.csv(epmc_hf, "data/epmc_hf_trials.csv", row.names = FALSE)
}
epmc_hf <- read.csv("data/epmc_hf_trials.csv", stringsAsFactors = FALSE)
epmc_hf |>
  select(title, journalTitle, pubYear, pmid) |>
  head(3) |>
  mutate(title = str_trunc(title, 70)) |>
  kable(col.names = c("Title", "Journal", "Year", "PMID"))


## -----------------------------------------------------------------------------
#| label: counts
# Stand-ins for objects your own pipeline would create
central_search <- list(count = 841)       # e.g. from a CENTRAL export
medline_search <- list(count = 1923)      # e.g. rentrez::entrez_search()
embase_search  <- list(count = 2417)
registers_n    <- 142
greylit_n      <- 37
deduped_df     <- data.frame(id = 1:3961)  # after deduplication
passed_ta_df   <- data.frame(id = 1:420)   # passed title/abstract screening
included_df    <- data.frame(id = 1:42,
                             in_meta_analysis = c(rep(TRUE, 38), rep(FALSE, 4)))

# One object holding every number in the flow diagram
counts <- list(
  central      = central_search$count,
  medline      = medline_search$count,
  embase       = embase_search$count,
  other        = registers_n + greylit_n,
  no_dupes     = nrow(deduped_df),
  ta_excluded  = nrow(deduped_df) - nrow(passed_ta_df),
  full_text    = nrow(passed_ta_df),
  ft_excluded  = nrow(passed_ta_df) - nrow(included_df),
  qualitative  = nrow(included_df),
  quantitative = sum(included_df$in_meta_analysis))
str(counts)


## -----------------------------------------------------------------------------
#| label: update-search
#| eval: false
# last_run <- as.Date("2024-12-01")
# today    <- Sys.Date()
# new_results <- fetch_rcts(
#   mesh_terms = "heart failure",
#   start      = format(last_run, "%Y/%m/%d"),
#   end        = format(today, "%Y/%m/%d"),
#   out_csv    = paste0("data/update_", today, ".csv"))
# cat(nrow(new_results), "new articles since", format(last_run), "\n")

