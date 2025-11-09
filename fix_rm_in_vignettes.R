# ==============================================================================
# Script to fix {\rm } in vignettes
# ==============================================================================
# This script replaces {\rm X} with \text{X} in all Rmd files
# Run this from the package root directory
# ==============================================================================

# Get all Rmd files in vignettes/
vignette_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)

cat("Found", length(vignette_files), "vignette files:\n")
print(basename(vignette_files))

# Function to replace {\rm X} with \text{X}
fix_rm <- function(file_path) {
  cat("\nProcessing:", basename(file_path), "\n")

  # Read the file
  content <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  original_content <- content

  # Replace {\rm ...} with \text{...}
  # This pattern matches {\rm followed by space and any text until }
  content <- gsub("\\{\\\\rm\\s+([^}]+)\\}", "\\\\text{\\1}", content)

  # Count changes
  n_changes <- sum(content != original_content)

  if (n_changes > 0) {
    cat("  →", n_changes, "lines modified\n")

    # Backup original file
    backup_file <- paste0(file_path, ".backup")
    writeLines(original_content, backup_file, useBytes = TRUE)
    cat("  → Backup saved:", basename(backup_file), "\n")

    # Write modified content
    writeLines(content, file_path, useBytes = TRUE)
    cat("  → File updated\n")

    return(TRUE)
  } else {
    cat("  → No changes needed\n")
    return(FALSE)
  }
}

# Process all vignette files
cat("\n=== Starting Conversion ===\n")
results <- sapply(vignette_files, fix_rm)

# Summary
cat("\n=== Summary ===\n")
cat("Total files processed:", length(vignette_files), "\n")
cat("Files modified:", sum(results), "\n")
cat("Files unchanged:", sum(!results), "\n")

cat("\n=== Next Steps ===\n")
cat("1. Review the changes (backup files saved with .backup extension)\n")
cat("2. Delete docs/ folder\n")
cat("3. Run: pkgdown::build_site()\n")
cat("4. Check that math renders correctly\n")
cat("5. If satisfied, delete .backup files:\n")
cat("   unlink('vignettes/*.backup')\n")

cat("\n=== Done ===\n")
