# ==============================================================================
# pkgdown Website Build Script
# ==============================================================================
# This script builds the pkgdown website for the twoCoprimary package.
# Run this script from the package root directory in RStudio.
#
# Usage:
#   1. Open this file in RStudio
#   2. Source the entire script (Ctrl/Cmd + Shift + S)
#   3. Or run sections individually
# ==============================================================================

# Install pkgdown if not already installed
if (!requireNamespace("pkgdown", quietly = TRUE)) {
  install.packages("pkgdown")
}

# Load required packages
library(pkgdown)

# ==============================================================================
# STEP 1: Verify Setup
# ==============================================================================

cat("\n=== Checking Package Structure ===\n")

# Check if _pkgdown.yml exists
if (file.exists("_pkgdown.yml")) {
  cat("✓ _pkgdown.yml found\n")
} else {
  stop("✗ _pkgdown.yml not found. Please create it first.")
}

# Check if NEWS.md exists
if (file.exists("NEWS.md")) {
  cat("✓ NEWS.md found\n")
} else {
  warning("⚠ NEWS.md not found. Creating a basic one...")
  writeLines("# twoCoprimary 0.1.0\n\nInitial release.", "NEWS.md")
}

# Check if DESCRIPTION exists
if (file.exists("DESCRIPTION")) {
  cat("✓ DESCRIPTION found\n")
} else {
  stop("✗ DESCRIPTION not found. Are you in the package root directory?")
}

# Check if man/figures directory exists (for logo)
if (!dir.exists("man/figures")) {
  cat("⚠ man/figures directory does not exist. Creating it...\n")
  dir.create("man/figures", recursive = TRUE)
  cat("✓ man/figures directory created\n")
}

# ==============================================================================
# STEP 2: Clean Previous Build (Optional)
# ==============================================================================

cat("\n=== Cleaning Previous Build ===\n")

if (dir.exists("docs")) {
  cat("Removing old docs/ directory...\n")
  unlink("docs", recursive = TRUE)
  cat("✓ Old docs/ directory removed\n")
} else {
  cat("✓ No previous docs/ directory to clean\n")
}

# ==============================================================================
# STEP 3: Build the Website
# ==============================================================================

cat("\n=== Building pkgdown Website ===\n")
cat("This may take a few minutes...\n\n")

# Build the complete site
build_site(
  pkg = ".",
  examples = TRUE,
  run_dont_run = FALSE,
  seed = 12345,
  lazy = FALSE,
  override = list(destination = "docs"),
  preview = TRUE,  # Automatically open in browser/viewer
  devel = FALSE,
  new_process = TRUE
)

cat("\n✓ Website build complete!\n")

# ==============================================================================
# STEP 4: Verify Build
# ==============================================================================

cat("\n=== Verifying Build ===\n")

# Check if docs directory was created
if (dir.exists("docs")) {
  cat("✓ docs/ directory created\n")

  # Check key files
  key_files <- c(
    "docs/index.html",
    "docs/reference/index.html",
    "docs/articles/index.html"
  )

  for (file in key_files) {
    if (file.exists(file)) {
      cat("✓", file, "created\n")
    } else {
      warning("⚠", file, "not found")
    }
  }
} else {
  stop("✗ docs/ directory was not created. Check for errors above.")
}

# ==============================================================================
# STEP 5: Open Website
# ==============================================================================

cat("\n=== Opening Website ===\n")

# Open the website in browser
if (file.exists("docs/index.html")) {
  browseURL("docs/index.html")
  cat("✓ Website opened in browser\n")
} else {
  warning("⚠ Could not open website. Please open docs/index.html manually.")
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== Build Summary ===\n")
cat("Website successfully built in: docs/\n")
cat("\nNext steps:\n")
cat("1. Review the website in your browser\n")
cat("2. Check that all vignettes are included\n")
cat("3. Verify that function documentation is complete\n")
cat("4. (Optional) Add a logo to man/figures/logo.png\n")
cat("5. Push to GitHub and enable GitHub Pages\n")
cat("\nTo rebuild specific sections:\n")
cat("  - Home page:   pkgdown::build_home()\n")
cat("  - Reference:   pkgdown::build_reference()\n")
cat("  - Articles:    pkgdown::build_articles()\n")
cat("  - News:        pkgdown::build_news()\n")
cat("\n=== Done! ===\n")
