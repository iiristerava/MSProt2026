
# A custom function to read FASTA files without stripping headers
read_fasta_full <- function(file_path) {
  lines <- readLines(file_path)
  is_header <- grepl("^>", lines)
  
  headers <- lines[is_header]
  
  # Find start and end indices for each sequence block
  starts <- which(is_header) + 1
  ends <- c(which(is_header)[-1] - 1, length(lines))
  
  # Collapse multi-line sequences into a single string per header
  seqs <- mapply(function(s, e) paste(lines[s:e], collapse=""), starts, ends)
  names(seqs) <- headers
  return(seqs)
}

# Read the files using the new robust function
contaminants <- read_fasta_full("../capstoneProject/contaminants.fasta")
f26v2 <- read_fasta_full("../capstoneProject/uniprotkb_AND_model_organism_9606_AND_r_2026_05_18.fasta")

# Extract core Accession ID for Contaminants (Text between '>' and first space)
# Example: ">P00761 SWISS-PROT..." -> "P00761"
acc_contam <- sub("^>([^ ]+).*", "\\1", names(contaminants))

# Extract core Accession ID for Main Database
# Example: ">sp|A0A087X1C5|CP2D7_HUMAN..." -> "A0A087X1C5"
acc_f26v2 <- sapply(strsplit(names(f26v2), "\\|"), function(x) {
  if (length(x) >= 2) x[2] else x[1]
})
acc_f26v2 <- sub("\\.[0-9]+$", "", acc_f26v2) # Remove version .1 if present

# Filter out duplicate human proteins
unique_contam_indices <- which(!(acc_contam %in% acc_f26v2))
unique_contaminants <- contaminants[unique_contam_indices]

# Function to standardize contaminant headers to match UniProt format
standardize_header <- function(header) {
  # Remove the leading ">"
  clean_header <- sub("^>", "", header)
  
  # Split the header by spaces
  parts <- strsplit(clean_header, " ")[[1]]
  
  # The first part is the Accession ID (e.g., "P00761" or "Q32MB2")
  acc_id <- parts[1]
  
  # The rest is the description
  description <- paste(parts[-1], collapse = " ")
  
  # Create a new UniProt-style header: >sp|Accession|CON_Accession Description
  # Adding "CON_" helps you identify contaminants easily in Spectronaut results
  new_header <- paste0(">sp|", acc_id, "|CON_", acc_id, " ", description)
  
  return(new_header)
}

# Apply the function to reformat the headers of your unique contaminants
names(unique_contaminants) <- sapply(names(unique_contaminants), standardize_header)

# Let's check the result to ensure they match f26v2!
print(names(unique_contaminants)[1:3])

# NOW you can safely combine them!
# Combine databases
final_combined_db <- c(unique_contaminants, f26v2)

cat("--- Merging Report ---\n")
cat("Total human proteins:    ", length(acc_f26v2), "\n")
cat("Total Instructor's contaminants:    ", length(contaminants), "\n")
cat("Duplicate human proteins removed:", length(contaminants) - length(unique_contaminants), "\n")
cat("Unique contaminants appended:       ", length(unique_contaminants), "\n")
cat("Total proteins in final DB:         ", length(final_combined_db), "\n")

# Write to a new FASTA file with FULL headers
output_path <- "./Comb_uniprotkb_2026_05_18_Contaminants.fasta"

# Function to save FASTA with specific line length (default 60)
save_wrapped_fasta <- function(db, file_path, line_length = 60) {
  # Open a connection to write to the file
  con <- file(file_path, "w")
  
  for (i in seq_along(db)) {
    # Write the header
    writeLines(names(db)[i], con)
    
    # Get the full sequence string
    seq_str <- as.character(db[i])
    
    # Calculate start and end positions for 60-character chunks
    starts <- seq(1, nchar(seq_str), by = line_length)
    ends <- pmin(starts + line_length - 1, nchar(seq_str))
    
    # Extract the chunks and write them to the file
    chunks <- substring(seq_str, starts, ends)
    writeLines(chunks, con)
  }
  
  # Close the file connection
  close(con)
}

# Save the final combined database with 60 amino acids per line
save_wrapped_fasta(final_combined_db, output_path, line_length = 60)

cat("\nSuccess! File saved at:\n", output_path, "\n")
