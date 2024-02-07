library(rentrez)
library(XML)
library(tidytext)
library(ggplot2)
library(ggrepel)
library(stopwords)
library(tidyverse)
library(widyr)
library(tidytext)
library(igraph)
library(ggraph)

# query ----
# Define groups of search terms. Each group's terms will be combined using "OR".
# The final query will require each group to be present by combining them with "AND".
term_groups <- list(
  omics_terms = c("omics", "omic"),
  # cohort_terms = c("paediatric", "paediatrics", "pediatric", "pediatrics"), # comment out to remove pediatric-only
  disease_terms = c("sepsis", "septic")
)

# Function to construct a query fragment for a single group of terms.
# This will join terms within a group using "OR".
construct_group_query <- function(terms) {
  if (length(terms) == 0) {
    # If a group is empty, it's excluded from the final query.
    return(NULL)
  }
  # Combine terms within a group using "OR" and add field restriction "[Abstract]".
  grouped_terms <- paste0(terms, "[Abstract]", collapse = " OR ")
  return(sprintf("(%s)", grouped_terms))
}

# Construct query fragments for each group and exclude any empty groups.
query_fragments <- lapply(term_groups, construct_group_query) %>% 
  Filter(Negate(is.null), .)

# Combine the query fragments of non-empty groups using "AND".
query <- paste(query_fragments, collapse = " AND ")

# Print the final PubMed query
query

print(query)

# Specify the number of records you want to fetch
num_records_to_fetch <- 300  # default number is 20 articles

# Search for publications in PubMed
papers <-
  entrez_search(
    db = "pubmed",
    term = query,
    use_history = TRUE,
    retmax = num_records_to_fetch
  )

length(papers$ids)

# Fetch records in "abstract" format
xmlData <-
  entrez_fetch(
    db = "pubmed",
    id = papers$ids,
    rettype = "abstract",
    retmode = "xml"
  )

# Parse XML
doc <- xmlParse(xmlData)

# Extract abstract nodes
abstractNodes <- getNodeSet(doc, "//AbstractText")

# Extract the content of abstract nodes
abstracts <- sapply(abstractNodes, xmlValue)

# Combine all abstracts into a single string, separated by two newline characters
all_abstracts <- paste(abstracts, collapse = "\n\n")

# Write all abstracts to a single file
writeLines(all_abstracts, "../data/all_abstracts.txt")

# Plot term frequency data ----

# Create a tibble
abstracts_tibble <- tibble(abstract = all_abstracts)

# Tokenize the abstracts
tokenized_abstracts <- abstracts_tibble %>%
  unnest_tokens(word, abstract)

# Remove common "stop words"
data("stop_words")
tokenized_abstracts <- tokenized_abstracts %>%
  anti_join(stop_words)

# Highlight words ending with "omic" or "omics" and count the frequency of each word
word_counts <- tokenized_abstracts %>%
  mutate(TermOfInterest = ifelse(grepl("omics?$", word), "Yes", "No")) %>%
  count(word, TermOfInterest, sort = TRUE)

# Make the plot, showing the most common words:
title1 <-
  "Most common words in \nabstracts related to terms of interest (quantile 0.95)"

plot1 <- word_counts %>%
  # mutate(TermOfInterest = ifelse(word %in% terms, "Yes", "No")) %>%
  filter(n > quantile(n, 0.95)) %>%  # Optional: Only show the most common words
  ggplot(aes(reorder(word, n), n, fill = TermOfInterest)) +
  geom_col() +
  geom_text(aes(label = ifelse(
    TermOfInterest == "Yes", as.character(word), ""
  )), hjust = -0.1) +
  scale_fill_manual(values = c("Yes" = "#fa7e1e", "No" = "#4f5bd5")) +
  coord_flip() +
  theme_bw() +
  labs(x = "Word",
       y = "Frequency",
       title = title1,
       fill = "Term of interest") +
  theme_bw()

plot1

# Save it as a PDF
ggsave(filename = "../data/figures/abstracts_word_frequency.pdf",
       plot = plot1)

# Plot TF-IDF ----
# Split the single string of all_abstracts into individual abstracts
abstracts_tibble <-
  tibble(abstract = strsplit(all_abstracts, "\n\n")[[1]])

# Add a unique identifier for each abstract:
abstracts_tibble <- abstracts_tibble %>%
  mutate(abstract_id = row_number())

# Tokenize the abstracts with the unique identifier:
tokenized_abstracts <- abstracts_tibble %>%
  unnest_tokens(word, abstract)

# Remove common "stop words" (like "and", "the", "of", etc.) as they are usually not informative:
data(stop_words)
tokenized_abstracts <- tokenized_abstracts %>%
  anti_join(stop_words)

# Calculate TF-IDF:
tf_idf <- tokenized_abstracts %>%
  count(abstract_id, word) %>%
  bind_tf_idf(word, abstract_id, n) %>%
  arrange(desc(tf_idf))  %>%
  # Highlight words ending with "omic" or "omics"
  mutate(TermOfInterest = ifelse(grepl("omics?$", word), "Yes", "No"))

# plot
max_terms <- 5000

title2 <- paste0(
  "TF-IDF for terms of interest in \nabstracts (top ", max_terms, " terms)")

plot2 <- tf_idf %>%
  arrange(desc(tf_idf)) %>%
  mutate(Rank = row_number()) %>%  # Add a new 'Rank' column
  filter(Rank < max_terms) %>%  # Optional: Only show the most common words
  ggplot(aes(x = Rank, y = tf_idf, fill = TermOfInterest)) +
  geom_col() +
  scale_fill_manual(values = c("Yes" = "#fa7e1e", "No" = "#4f5bd5")) +
  geom_text_repel(
    aes(label = ifelse(
      TermOfInterest == "Yes", as.character(word), ""
    )),
    max.overlaps = 10,
    alpha = .8,
    angle = 90,
    hjust = -20,
    vjust = 4,
    size = 2
  ) +
  theme_bw() +
  labs(x = "Rank",
       y = "TF-IDF",
       title = title2,
       fill = "Term of interest") +
  theme_bw()

plot2

# Save it  as a PDF
ggsave(filename = "../data/figures/abstracts_tf_idf.pdf",
       plot = plot2,
       # units = "cm", height = 10, width = 10
       )


# Plot correlation ----
# Compute a co-occurrence matrix
co_occur <- tokenized_abstracts %>%
  pairwise_count(word, abstract_id, sort = TRUE, upper = FALSE)

# Create an igraph object
g <- graph_from_data_frame(co_occur)

# Calculate the edge weights based on the frequency of co-occurrence
E(g)$weight <- E(g)$n

# Prune the graph to only keep edges with weights above a certain threshold
g <-
  delete_edges(g, E(g)[E(g)$weight < quantile(E(g)$weight, 0.95)])

# Calculate the layout of the graph
layout <- layout_with_fr(g)

# Add a node attribute based on degree
V(g)$high_degree <- degree(g) > quantile(degree(g), 0.98)

# Plot
title3 <-
  "Network plot of edge weights based \non the frequency of co-occurrence in\nabstracts related to terms of interest (quantile  0.98)"

plot3 <- ggraph(g, layout = layout) +
  geom_edge_link(aes(edge_alpha = weight),
                 edge_width = 0.5,
                 show.legend = FALSE) +
  geom_node_point(size = .3, alpha = 0.5) +
  geom_node_text(
    aes(label = ifelse(high_degree, name, "")),
    color = "#4f5bd5",
    repel = TRUE,
    max.overlaps = 50,
    hjust = -0.6,
    size = 4
  ) +
  theme_void()
# labs(title = title3, size=0.6)

# Save it  as a PDF
# ggsave(
#   filename = "../data/figures/network_plot.pdf",
#   plot = plot3,
#   width = 6,
#   height = 6,
#   units = "in"
# )

# Plot heatmap ----
library(ggplot2)
library(reshape2)

# First, we need to reshape the data to a long format
# # Filter to include only word pairs that co-occur more than a certain number of times
threshold <- 50  # Change this to suit your needs
# threshold <- 10

co_occur_long <- melt(co_occur) %>%
  filter(value > threshold)

# Plot
title4 <-
  paste0(
    "Heatmap of term co-occurrence in \nabstracts related to terms of interest\n(co-occurrence threshold >",
    threshold,
    ")"
  )

plot4 <-
  ggplot(co_occur_long, aes(x = item1, y = item2, fill = value)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Word 1",
       y = "Word 2",
       fill = "Co-occurrence",
       title = title4) +
  theme_bw()


plot4

# Save it  as a PDF
ggsave(filename = "../data/figures/abstracts_heatmap.pdf",
       plot = plot4)

# omics count ----
# Extract only "omics"-related terms from your terms list
# Assuming 'omics' and similar terms are in terms_group1
omics_terms <- term_groups[1]

# Tokenize the abstracts
tokenized_abstracts <- abstracts_tibble %>%
  unnest_tokens(word, abstract) %>%
  mutate(word = tolower(word)) # Ensure matching is case-insensitive

# Flag occurrences of "omics"-related terms in abstracts
# Here, filtering directly for words that end with 'omic' or 'omics' from the entire abstract set
flagged_abstracts <- tokenized_abstracts %>%
  mutate(TermMatch = ifelse(grepl("omics?$", word), word, NA)) %>%
  filter(!is.na(TermMatch))

# Count how many abstracts contain each "omics"-related term
articles_per_term <- flagged_abstracts %>%
  group_by(TermMatch) %>%
  summarise(ArticleCount = n_distinct(abstract_id)) %>%
  ungroup()

# Visualize the number of articles mentioning each specific "omics"-related term
plot5 <-
  ggplot(articles_per_term,
         aes(x = TermMatch, y = ArticleCount, fill = TermMatch)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "Number of articles mentioning\neach 'omics'-related term",
       x = "'Omics'-related term",
       y = "Number of articles") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none")

plot5

# Save it  as a PDF
ggsave(filename = "../data/figures/omics_count.pdf",
       plot = plot5)

# co-occurance omics ----
# Assuming 'flagged_abstracts' contains the flagged occurrences of "omics"-related terms

# Step 1: Group by article and summarize unique "omics" terms
articles_omics_cooccurrence <- flagged_abstracts %>%
  group_by(abstract_id) %>%
  summarise(UniqueOmicsTerms = n_distinct(TermMatch)) %>%
  ungroup()

# Step 2: Summarize how many articles have each count of unique "omics" terms
cooccurrence_counts <- articles_omics_cooccurrence %>%
  count(UniqueOmicsTerms) %>%
  arrange(UniqueOmicsTerms)

# Step 3: Visualize the co-occurrence distribution
plot_cooccurrence <-
  ggplot(cooccurrence_counts, aes(
    x = as.factor(UniqueOmicsTerms),
    y = n,
    fill = as.factor(UniqueOmicsTerms)
  )) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(
    begin = 0.3,
    end = 0.9,
    direction = 1,
    name = "Unique 'Omics' Terms"
  ) +
  labs(title = "Co-occurrence of 'omics' terms in articles",
       x = "Number of unique 'omics' terms",
       y = "Number of articles") +
  theme_bw() +
  theme(legend.position = "none")

# Display the plot
plot_cooccurrence

# Save the plot
ggsave(filename = "../data/figures/omics_cooccurrence.pdf",
       plot = plot_cooccurrence)


# co-occurance breakdown ----
# Function to generate all combinations of terms of all lengths
generate_all_combinations <- function(terms) {
  if (length(terms) < 2) {
    return(NULL) # Return NULL for articles with fewer than 2 terms
  }
  # Generate combinations for each length
  all_combinations <- lapply(2:length(terms), function(n) {
    combn(terms, n, simplify = FALSE)
  })
  # Flatten the list of combinations
  do.call(c, all_combinations)
}

# Apply the function to each article
combinations_per_article <- flagged_abstracts %>%
  select(abstract_id, TermMatch) %>%
  distinct() %>%
  filter(!TermMatch %in% c("omic", "omics")) %>% # Exclude general terms
  group_by(abstract_id) %>%
  summarise(Combinations = list(generate_all_combinations(unique(TermMatch)))) %>%
  filter(lengths(Combinations) > 0) %>% # Filter out articles without combinations
  unnest(Combinations) %>%
  mutate(Combination = map_chr(Combinations, ~ paste(sort(.x), collapse = " & "))) %>%
  select(-Combinations)

combination_counts <- combinations_per_article %>%
  count(Combination, sort = TRUE)

# set the minimum number of cobinations to plot
min_combinations = 2

plot_combination_counts <-
  combination_counts %>%
  filter(n > min_combinations) %>%
  ggplot(aes(
    x = reorder(Combination, n),
    y = n,
    fill = Combination
  )) +
  geom_col() +
  coord_flip() +
  labs(title = paste0("Frequency of 'omics' \nterm combinations (n > ", min_combinations, " )"),
       x = "Combination of 'omics' terms",
       y = "Number of articles") +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none")

plot_combination_counts

# Save the plot
ggsave(filename = "../data/figures/omics_combination_counts.pdf",
       plot = plot_combination_counts
       # units = "cm",
       # height = 20,
       # width = 20
       )


# patchwork ----
library(patchwork)
# patch0 <- (plot1 | plot2) / (plot3 | plot4) + plot_annotation(tag_levels = 'A')
patch0 <-
  (plot1 | plot2) / (plot4) + plot_annotation(tag_levels = 'A')
patch1 <- (plot1 | plot2) + plot_annotation(tag_levels = 'A')
# patch2 <- (plot3 | plot4) + plot_annotation(tag_levels = 'A')

# # Save it  as a PDF
# ggsave(
#   filename = "../data/figures/plot_patched_0.pdf",
#   plot = patch0,
#   width = 10,
#   height = 12,
#   units = "in"
# )
# ggsave(
#   filename = "../data/figures/plot_patched_1.pdf",
#   plot = patch1,
#   width = 16,
#   height = 8,
#   units = "in"
# )
# ggsave(
#   filename = "../data/figures/plot_patched_2.pdf",
#   plot = patch2,
#   width = 12,
#   height = 8,
#   units = "in"
# )

# Figure legends ----
title1
title2
title3
title4

# Metrics ----
length(papers$ids)
papers$ids
papers[["QueryTranslation"]]

# Pubmed URLs ----
# Concatenate the main URL with each PMID
pmid_urls <-
  paste0("https://pubmed.ncbi.nlm.nih.gov/", papers$ids)

# Save PubMed URLs to a text file in one step
writeLines(pmid_urls, "../data/pmid_urls.txt")

# Or get pubmed IDs to paste to Zotero
# Concatenate PubMed IDs with commas
pmids_for_zotero <- paste(papers$ids, collapse = ", ")

# Save the formatted PubMed IDs to a text file
writeLines(pmids_for_zotero, "../data/pmids_for_zotero.txt")
