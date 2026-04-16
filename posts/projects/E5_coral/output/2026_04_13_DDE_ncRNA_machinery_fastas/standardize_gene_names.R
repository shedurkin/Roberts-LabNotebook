#!/usr/bin/env Rscript
# standardize_gene_names.R
# ─────────────────────────────────────────────────────────────────────
# Adds a "CGN=" (Corrected Gene Name) field to each FASTA header
# based on a curated mapping of UniProt protein descriptions to
# standard HGNC/community-accepted gene symbols.
#
# Usage:
#   Rscript standardize_gene_names.R input.fasta output.fasta
#
# The mapping table is defined inline below. For new protein families
# not yet in the table, the script flags them as "UNMAPPED" so you
# can curate them manually.
# ─────────────────────────────────────────────────────────────────────

library(tidyverse)

# ── 1. Curated description → standard gene symbol mapping ───────────
#
# Each regex pattern is matched (case-insensitive) against the protein
# description field extracted from the FASTA header.  Order matters:
# more specific patterns should come before general ones.
#
# Format: list(pattern = "GENE_SYMBOL")

desc_to_gene <- tribble(
  ~pattern,                                                    ~gene_symbol,
  # ── Cap-binding complex ──
  "Nuclear cap-binding protein subunit 1",                     "NCBP1",
  "Nuclear cap-binding protein subunit 2",                     "NCBP2",
  
  # ── Heat shock ──
  "Heat shock protein 90",                                     "HSP90",
  
  # ── SWI/SNF chromatin remodeling ──
  "SWI/SNF complex subunit SMARCC2",                           "SMARCC2",
  "SWI/SNF-related.*SMARCA2",                                  "SMARCA2",
  "Transcription activator BRG1|SMARCA4",                      "SMARCA4",
  
  # ── RNA helicases ──
  "RNA helicase.*DDX6",                                        "DDX6",
  "RNA helicase.*DDX17",                                       "DDX17",
  "ATP-dependent RNA helicase A|DHX9",                         "DHX9",
  "Exosome RNA helicase MTR4",                                 "MTREX",
  
  # ── mRNA decapping ──
  "mRNA-decapping enzyme 1B",                                  "DCP1B",
  "mRNA-decapping enzyme 1A",                                  "DCP1A",
  "mRNA-decapping enzyme C-terminal",                          "DCP1",
  "Enhancer of mRNA-decapping protein 3",                      "EDC3",
  "Enhancer of mRNA-decapping protein 4|Enhancer of mRNA.decapping 4", "EDC4",
  
  # ── DIS3 family exonucleases ──
  "DIS3-like exonuclease 2",                                   "DIS3L2",
  "DIS3-like exonuclease 1",                                   "DIS3L",
  "DIS3 mitotic control",                                      "DIS3L",
  "Protein DIS3 homolog",                                      "DIS3",
  
  # ── RNA-binding / processing ──
  "RNA-binding protein EWS",                                   "EWSR1",
  "RNA-binding.*FUS",                                           "FUS",
  "RNA-binding protein 3",                                     "RBM3",
  "RNA-binding protein 7",                                     "RBM7",
  "Heterogeneous nuclear ribonucleoprotein A1",                "HNRNPA1",
  "Heterogeneous nuclear ribonucleoproteins A2",               "HNRNPA2B1",
  "Heterogeneous nuclear ribonucleoprotein K",                 "HNRNPK",
  "HYL1-like",                                                 "HYL1",
  "TAR DNA-binding protein 43",                                "TARDBP",
  "Nuclease-sensitive element-binding protein 1|Y-box-binding protein 1", "YBX1",
  "Non-POU domain-containing octamer-binding",                 "NONO",
  "Splicing factor, proline- and glutamine-rich",              "SFPQ",
  "Paraspeckle component 1",                                   "PSPC1",
  "Far upstream element-binding protein 2|KHSRP",              "KHSRP",
  "Polyadenylate-binding protein",                             "PABPC1",
  
  # ── ADAR editing ──
  "Double-stranded RNA-specific adenosine deaminase",          "ADAR",
  "Double-stranded RNA-specific editase 1",                    "ADARB1",
  
  # ── Bromodomain proteins ──
  "Bromodomain-containing protein 4",                          "BRD4",
  "Bromodomain-containing protein 3",                          "BRD3",
  
  # ── CTCF ──
  "Transcriptional repressor CTCF",                            "CTCF",
  
  # ── Germ cell / DAZL ──
  "Deleted in azoospermia-like",                               "DAZL",
  
  # ── DNA methyltransferases ──
  "DNA \\(cytosine-5\\)-methyltransferase 3A",                "DNMT3A",
  "DNA \\(cytosine-5\\)-methyltransferase 3B",                "DNMT3B",
  "DNA \\(cytosine-5\\)-methyltransferase 1|DNA \\(Cytosine-5\\)-methyltransferase 1", "DNMT1",
  "DNA \\(cytosine-5\\)-methyltransferase|DNA \\(cytosine-5-\\)-methyltransferase", "DNMT1",
  
  # ── Exosome complex ──
  "Exosome complex component RRP45|EXOSC9",                   "EXOSC9",
  "Exosome complex component RRP46|Exosome component 5",      "EXOSC5",
  "Exosome complex component RRP4(?!0|1)",                     "EXOSC2",
  "Exosome complex component RRP41",                           "EXOSC4",
  "Exosome complex component RRP40",                           "EXOSC3",
  "Exosome complex component MTR3",                            "EXOSC6",
  "Ribosomal RNA-processing protein 42|EXOSC7",               "EXOSC7",
  "Ribosomal RNA-processing protein 43|EXOSC8",               "EXOSC8",
  "Exosome complex component 10|EXOSC10",                     "EXOSC10",
  "Exosome complex component CSL4|EXOSC1",                    "EXOSC1",
  
  # ── PRC2 / Polycomb ──
  "Histone-lysine N-methyltransferase EZH1",                  "EZH1",
  "\\[histone H3\\]-lysine\\(27\\) N-trimethyltransferase|EZH2", "EZH2",
  "Protein Jumonji|JARID2",                                    "JARID2",
  "Polycomb protein EED|Embryonic ectoderm development",       "EED",
  "Polycomb protein suz12|Polycomb protein SUZ12",             "SUZ12",
  "PHD finger protein 1|phf1",                                 "PHF1",
  "Zinc finger protein AEBP2|Zinc finger protein aebp2",      "AEBP2",
  "Histone-binding protein RBBP4",                             "RBBP4",
  "WD repeat-containing protein 5|WDR5",                       "WDR5",
  
  # ── Histone deacetylases ──
  # Note: many cnidarian HDACs lack isoform-level annotation.
  # Where GN field specifies HDAC1/2/3, we use that; otherwise "HDAC"
  "histone deacetylase.*HDAC1|histone deacetylase OS=Stylophora", "HDAC1",
  "histone deacetylase.*HDAC2|Histone deacetylase.*HDAC2",    "HDAC2",
  "Histone deacetylase 3|histone deacetylase.*Hdac3|HDAC3",   "HDAC3",
  "[Hh]istone deacetylase",                                   "HDAC",  # generic fallback
  
  # ── Integrator complex ──
  "Integrator complex subunit 1(?![0-9])",                     "INTS1",
  "Integrator complex subunit 2(?![0-9])",                     "INTS2",
  "Integrator complex subunit 3(?![0-9])",                     "INTS3",
  "Integrator complex subunit 4(?![0-9])",                     "INTS4",
  "Integrator complex subunit 5(?![0-9])",                     "INTS5",
  "Integrator complex subunit 6-like|Integrator complex subunit 6.*ints6l", "INTS6L",
  "Integrator complex subunit 6(?![0-9])",                     "INTS6",
  "Integrator complex subunit 7(?![0-9])",                     "INTS7",
  "Integrator complex subunit 8(?![0-9])",                     "INTS8",
  "Integrator complex subunit 9(?![0-9])",                     "INTS9",
  "Integrator complex subunit 10",                             "INTS10",
  "Integrator complex subunit 11",                             "INTS11",
  "Integrator complex subunit 12",                             "INTS12",
  "Integrator complex subunit 13",                             "INTS13",
  "Integrator complex subunit 14",                             "INTS14",
  "Integrator complex subunit 15",                             "INTS15",
  
  # ── miRNA biogenesis / RNAi ──
  "Endoribonuclease Dicer|Dicer-like protein 1|Dicer 1",      "DICER1",
  "Dicer 2",                                                   "DICER2",
  "Ribonuclease 3|Drosha",                                     "DROSHA",
  "Pasha/DGCR8|DGCR8",                                        "DGCR8",
  "Protein argonaute-1",                                       "AGO1",
  "Protein argonaute-2|Argonaute-2|Argonaute 2|Putative argonaute-2", "AGO2",
  "GW182/TNRC6",                                               "TNRC6",
  "Protein Loquacious",                                        "LOQS",
  "SERRATE|Serrate RNA effector|Ars2/serrate",                 "SRRT",
  "Interferon-inducible.*protein kinase activator A",          "PRKRA",
  "Small RNA 2'-O-methyltransferase",                          "HENMT1",
  "RNA-dependent RNA polymerase",                              "RDRP",
  
  # ── CCR4-NOT deadenylation complex ──
  "CCR4-NOT transcription complex subunit 1(?![0-9])",         "CNOT1",
  "CCR4-NOT transcription complex subunit 2(?![0-9])",         "CNOT2",
  "CCR4-NOT transcription complex subunit 3(?![0-9])",         "CNOT3",
  "CCR4-NOT transcription complex subunit 4(?![0-9])",         "CNOT4",
  "CCR4-NOT transcription complex subunit 9",                  "CNOT9",
  "CCR4-NOT transcription complex subunit 10",                 "CNOT10",
  "CCR4-NOT transcription complex subunit 11",                 "CNOT11",
  "CCR4-Not complex component Not",                            "CNOT",
  
  # ── PAN2-PAN3 deadenylation ──
  "PAN2-PAN3.*subunit PAN3",                                   "PAN3",
  "PAN2-PAN3.*catalytic subunit PAN2",                          "PAN2",
  
  # ── Poly(A) machinery ──
  "polynucleotide adenylyltransferase|PAPD5",                  "PAPD5",
  "Poly\\(A\\)-specific ribonuclease PARN",                    "PARN",
  "polyribonucleotide nucleotidyltransferase",                 "PNPT1",
  
  # ── TUTases / LIN28 ──
  "Terminal uridylyltransferase 4",                             "TUT4",
  "Terminal uridylyltransferase 7",                             "TUT7",
  "Protein lin-28",                                            "LIN28A",
  
  # ── Nuclear export ──
  "Exportin-1",                                                "XPO1",
  "Exportin.5",                                                "XPO5",
  "Phosphorylated adapter RNA export protein",                 "PHAX",
  "GTP-binding nuclear protein",                               "RAN",
  
  # ── Nucleolar snoRNP ──
  "Nucleolar protein 56",                                      "NOP56",
  "Nucleolar protein 58",                                      "NOP58",
  
  # ── Transcriptional repressors / other ──
  "Transcriptional repressor NF-X1|NF-X1-type",               "NFX1",
  
  # ── Cohesin ──
  "Double-strand-break repair.*rad21|Rad21/Rec8",              "RAD21",
  "Cohesin subunit SA-1",                                      "STAG1",
  "Cohesin subunit SA-2",                                      "STAG2",
  "Structural maintenance of chromosomes protein.*Smc3|SMC3",  "SMC3",
  
  # ── Nonsense-mediated decay ──
  "Regulator of nonsense transcripts 1|UPF1|Upf1|helicase NAM7", "UPF1",
  "Protein SMG5",                                              "SMG5",
  "non-specific serine/threonine protein kinase.*Smg1",        "SMG1",
  
  # ── 5'-3' exoribonucleases ──
  "5'-3' exoribonuclease 1|Xrn1",                              "XRN1",
  "5'-3' exoribonuclease",                                     "XRN",
  
  # ── Ribonuclease P ──
  "Ribonuclease P protein subunit p20",                        "POP7",
  
  # ── Endoribonuclease ZC3H12A / Regnase ──
  "Ribonuclease ZC3H12A|Endoribonuclease ZC3H12A",           "ZC3H12A",
  "Zinc finger C3H1 domain-containing",                        "ZFC3H1",
  
  # ── tRNA methyltransferase ──
  "tRNA.*methyltransferase TARBP1",                            "TARBP1",
  
  # ── Polycomb / EED extra entries ──
  "Histone deacetylase OS=Hydractinia.*Hdac2",                 "HDAC2"
)


# ── 2. Parse FASTA headers ──────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript standardize_gene_names.R <input.fasta> <output.fasta>\n")
  cat("       Optionally add a third arg 'mapping' to also write a TSV mapping table.\n")
  quit(status = 1)
}

infile  <- args[1]
outfile <- args[2]
write_mapping <- length(args) >= 3 && args[3] == "mapping"

lines <- readLines(infile)

# Identify header lines
is_header <- str_starts(lines, ">")

# Extract protein description from UniProt FASTA format:
#   >db|accession|entry_name Description OS=... OX=... GN=... PE=... SV=...
parse_header <- function(header) {
  # Remove leading ">"
  h <- str_remove(header, "^>")
  
  # Extract accession (second field between pipes)
  acc <- str_extract(h, "(?<=\\|)[^|]+(?=\\|)")
  
  # Extract description: everything between entry_name and " OS="
  desc <- str_match(h, "\\|\\S+\\s+(.+?)\\s+OS=")[, 2]
  
  # Extract existing GN= field if present
  gn <- str_match(h, "GN=(\\S+)")[, 2]
  
  # Extract organism
  
  org <- str_match(h, "OS=(.+?)\\s+OX=")[, 2]
  
  tibble(
    raw_header  = header,
    accession   = acc,
    description = desc,
    existing_gn = gn,
    organism    = org
  )
}

headers <- map_dfr(lines[is_header], parse_header)


# ── 3. Map descriptions to standard gene symbols ────────────────────

assign_gene <- function(desc, existing_gn) {
  if (is.na(desc)) return("UNMAPPED")
  
  for (i in seq_len(nrow(desc_to_gene))) {
    if (str_detect(desc, regex(desc_to_gene$pattern[i], ignore_case = TRUE))) {
      
      gene <- desc_to_gene$gene_symbol[i]
      
      # Special handling: for generic "HDAC", try to resolve from GN field
      if (gene == "HDAC" && !is.na(existing_gn)) {
        gn_upper <- toupper(existing_gn)
        if (str_detect(gn_upper, "HDAC1|^HDAC1$")) return("HDAC1")
        if (str_detect(gn_upper, "HDAC2|^HDAC2$")) return("HDAC2")
        if (str_detect(gn_upper, "HDAC3|^HDAC3$")) return("HDAC3")
      }
      
      # For generic DCP1, try to resolve from GN field
      if (gene == "DCP1" && !is.na(existing_gn)) {
        gn_upper <- toupper(existing_gn)
        if (str_detect(gn_upper, "DCP1B")) return("DCP1B")
        if (str_detect(gn_upper, "DCP1A")) return("DCP1A")
      }
      
      # For generic DNMT1, check if already matched a specific one
      return(gene)
    }
  }
  
  # If no pattern matched, check if GN field itself is already a
  
  # recognizable standard symbol (some reviewed sp| entries)
  if (!is.na(existing_gn)) {
    gn_clean <- toupper(existing_gn)
    # Return as-is if it looks like a real gene symbol (all caps, short)
    if (str_detect(gn_clean, "^[A-Z][A-Z0-9]{1,10}$")) {
      return(gn_clean)
    }
  }
  
  return("UNMAPPED")
}

headers <- headers %>%
  mutate(
    std_gene = map2_chr(description, existing_gn, assign_gene)
  )


# ── 4. Report mapping stats ─────────────────────────────────────────

n_total    <- nrow(headers)
n_mapped   <- sum(headers$std_gene != "UNMAPPED")
n_unmapped <- n_total - n_mapped

cat(sprintf("\n── Mapping summary ──\n"))
cat(sprintf("  Total headers:  %d\n", n_total))
cat(sprintf("  Mapped:         %d (%.1f%%)\n", n_mapped, 100 * n_mapped / n_total))
cat(sprintf("  Unmapped:       %d\n", n_unmapped))

if (n_unmapped > 0) {
  cat("\n── Unmapped entries ──\n")
  headers %>%
    filter(std_gene == "UNMAPPED") %>%
    select(accession, description, existing_gn) %>%
    print(n = Inf)
}

# Show unique gene assignments
cat("\n── Unique gene symbol assignments ──\n")
headers %>%
  count(std_gene, sort = TRUE) %>%
  print(n = Inf)


# ── 5. Write corrected FASTA ────────────────────────────────────────

# Insert CGN= field into each header, right before the OS= field
output_lines <- lines
header_idx   <- which(is_header)

for (i in seq_along(header_idx)) {
  li <- header_idx[i]
  gene <- headers$std_gene[i]
  # Insert CGN=SYMBOL before " OS="
  output_lines[li] <- str_replace(
    lines[li],
    " OS=",
    paste0(" CGN=", gene, " OS=")
  )
}

writeLines(output_lines, outfile)
cat(sprintf("\n── Output written to: %s ──\n", outfile))


# ── 6. Optionally write mapping table ───────────────────────────────

if (write_mapping) {
  mapping_file <- str_replace(outfile, "\\.fasta$", "_mapping.tsv")
  headers %>%
    select(accession, description, existing_gn, std_gene, organism) %>%
    write_tsv(mapping_file)
  cat(sprintf("── Mapping table written to: %s ──\n", mapping_file))
}

cat("\nDone.\n")