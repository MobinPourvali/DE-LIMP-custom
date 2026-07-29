# Tests for the proteogenomics helper modules
# (helpers_proteogenomics, helpers_proteog_assembly, helpers_slims,
#  helpers_rnaseq, helpers_proteog_qc).
#
# All synthetic-data tests run anywhere. One integration test against the
# May 20 validation FASTA is skipped unless the file is locally available
# (typically only on Brett's dev machine or Hive).
#
# Run manually:
#   testthat::test_file("tests/testthat/test-helpers_proteogenomics.R")

# =============================================================================
# count_proteog_classes() — composition counting
# =============================================================================

test_that("count_proteog_classes handles all 5 source classes in synthetic input", {
  test_fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(test_fasta), add = TRUE)
  writeLines(c(
    ">sp|P12345|GAPDH_HUMAN OS=Homo sapiens",                            # UNIPROT
    "MGKVKVGVNG",
    ">sp|ENSMUST00000000001.5.p1|Gnai3_MM39TEST source=REF ORF_type=complete strand=+ len=354 coords=ENSMUST00000000001.5:142-1206(+) parent_gene=ENSMUSG00000000001.5 transcript=ENSMUST00000000001.5",
    "MGCTLSAEDK",
    ">sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST source=NOVEL_GENE ORF_type=5prime_partial strand=- len=112 coords=MSTRG.10029.5:344-682(-) parent_gene=MSTRG.10029 transcript=MSTRG.10029.5",
    "MAGNREALY",
    ">sp|MSTRG.10075.2.p1|Trim25_MM39TEST source=NOVEL_ISOFORM ORF_type=complete strand=+ len=285 coords=MSTRG.10075.2:39-896(+) parent_gene=ENSMUSG00000005951 transcript=MSTRG.10075.2",
    "MASEHFVCK",
    ">INDEL_ENSP00000354813_81:CAAAAAAAACTC_CAAAAAAACTC",                # VARIANT
    "MAQGCATKL"
  ), test_fasta)

  counts <- count_proteog_classes(test_fasta)
  expect_equal(counts$total, 5)
  expect_equal(counts$UNIPROT, 1)
  expect_equal(counts$REF, 1)
  expect_equal(counts$NOVEL_GENE, 1)
  expect_equal(counts$NOVEL_ISOFORM, 1)
  expect_equal(counts$VARIANT, 1)
  expect_equal(counts$UNPARSED, 0)
})

test_that("count_proteog_classes handles UNPARSED entries", {
  test_fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(test_fasta), add = TRUE)
  writeLines(c(
    ">sp|X|Y source=UNPARSED reason=test",
    "AAAA"
  ), test_fasta)
  counts <- count_proteog_classes(test_fasta)
  expect_equal(counts$UNPARSED, 1)
  expect_equal(counts$total, 1)
})

test_that("count_proteog_classes against May 20 validation FASTA (integration)", {
  skip_if(
    !file.exists("/tmp/mm39_test_92bp_relaxed_RRNAfilt.fasta"),
    "May 20 validation FASTA not staged locally — skip integration test"
  )
  counts <- count_proteog_classes("/tmp/mm39_test_92bp_relaxed_RRNAfilt.fasta")
  expect_equal(counts$total, 67386)
  expect_equal(counts$REF, 66046)
  expect_equal(counts$NOVEL_GENE, 1340)
  expect_equal(counts$UNPARSED, 0)
})

# =============================================================================
# classify_proteins() — DIA-NN report → classification data.frame
# =============================================================================

# The `source=` / `ORF_type=` / `parent_gene=` tags live in the FASTA headers, not in
# the DIA-NN report — so classify_proteins() reads them from `fasta_path`, exactly as
# production does (server_data.R passes values$proteog_active_fasta). Write a temp
# FASTA carrying the same tags rather than expecting the report to carry them.
.write_test_fasta <- function(headers) {
  p <- tempfile(fileext = ".fasta")
  writeLines(as.vector(rbind(paste0(">", headers), "MPEPTIDEK")), p)
  p
}

test_that("classify_proteins extracts source/orf_type/parent_gene from the FASTA", {
  fasta <- .write_test_fasta(c(
    "sp|P12345|GAPDH_HUMAN GAPDH_HUMAN OS=Homo sapiens",
    "sp|ENSMUST00000000001.5.p1|Gnai3_MM39TEST source=REF ORF_type=complete strand=+ len=354 parent_gene=ENSMUSG00000000001.5 transcript=ENSMUST00000000001.5",
    "sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST source=NOVEL_GENE ORF_type=5prime_partial parent_gene=MSTRG.10029 transcript=MSTRG.10029.5",
    "sp|MSTRG.10075.2.p1|Trim25_MM39TEST source=NOVEL_ISOFORM ORF_type=complete parent_gene=ENSMUSG00000005951 transcript=MSTRG.10075.2"
  ))
  on.exit(unlink(fasta), add = TRUE)

  fake_genes <- data.frame(
    Protein.Group = c(
      "sp|P12345|GAPDH_HUMAN",
      "sp|ENSMUST00000000001.5.p1|Gnai3_MM39TEST",
      "sp|MSTRG.10029.5.p2|MSTRG.10029_MM39TEST",
      "sp|MSTRG.10075.2.p1|Trim25_MM39TEST",
      "INDEL_ENSP00000354813_81:CAAAAAAAACTC_CAAAAAAACTC"
    ),
    Protein.Group.Description = c(
      "GAPDH_HUMAN OS=Homo sapiens",
      "Gnai3_MM39TEST source=REF ORF_type=complete strand=+ len=354 coords=X parent_gene=ENSMUSG00000000001.5 transcript=ENSMUST00000000001.5",
      "MSTRG.10029_MM39TEST source=NOVEL_GENE ORF_type=5prime_partial parent_gene=MSTRG.10029 transcript=MSTRG.10029.5",
      "Trim25_MM39TEST source=NOVEL_ISOFORM ORF_type=complete parent_gene=ENSMUSG00000005951 transcript=MSTRG.10075.2",
      ""
    ),
    stringsAsFactors = FALSE
  )
  classification <- classify_proteins(fake_genes, fasta_path = fasta)
  expect_equal(nrow(classification), 5)
  expect_equal(classification$source[1], "UNIPROT")
  expect_equal(classification$source[2], "REF")
  expect_equal(classification$source[3], "NOVEL_GENE")
  # The whole point of trusting the FASTA: this accession is MSTRG.*.pN, so
  # structure alone says NOVEL_GENE. The FASTA knows it is a novel ISOFORM of a
  # known gene (parent_gene=ENSMUSG...) and must win.
  expect_equal(classification$source[4], "NOVEL_ISOFORM")
  expect_equal(classification$parent_gene[4], "ENSMUSG00000005951")
  expect_equal(classification$source[5], "VARIANT")
  expect_equal(classification$orf_type[2], "complete")
  expect_equal(classification$orf_type[3], "5prime_partial")
  expect_true(is.na(classification$orf_type[1]))
  expect_equal(classification$parent_gene[3], "MSTRG.10029")
})

test_that("classify_proteins falls back to accession structure with no FASTA", {
  # No fasta_path (a non-proteogenomics run): source still classifies structurally,
  # and the FASTA-only columns are NA rather than invented.
  cls <- classify_proteins(data.frame(
    Protein.Group = c("sp|P12345|GAPDH_HUMAN", "sp|MSTRG.10075.2.p1|Trim25_MM39TEST"),
    stringsAsFactors = FALSE))
  expect_equal(cls$source, c("UNIPROT", "NOVEL_GENE"))
  expect_true(all(is.na(cls$orf_type)))
  expect_true(all(is.na(cls$parent_gene)))
})

test_that("classify_proteins keeps source and parent_gene on the same member", {
  # A ;-joined group: the class must come from the winning member and parent_gene
  # must describe THAT member, not whichever happened to be tagged first.
  fasta <- .write_test_fasta(c(
    "sp|ENSMUST00000000001.5.p1|A source=REF ORF_type=complete parent_gene=ENSMUSG_REF",
    "sp|MSTRG.10075.2.p1|B source=NOVEL_ISOFORM ORF_type=complete parent_gene=ENSMUSG_ISO"
  ))
  on.exit(unlink(fasta), add = TRUE)
  cls <- classify_proteins(data.frame(
    Protein.Group = "sp|ENSMUST00000000001.5.p1|A;sp|MSTRG.10075.2.p1|B",
    stringsAsFactors = FALSE), fasta_path = fasta)
  expect_equal(cls$source, "NOVEL_ISOFORM")
  expect_equal(cls$parent_gene, "ENSMUSG_ISO")
})

test_that("classify_proteins accepts EList-like input with $genes slot", {
  # This is about the input SHAPE ($genes slot). The class comes from the FASTA tag,
  # as it does in production — bare accession "A" has no structure to classify.
  fasta <- .write_test_fasta("sp|A|B B source=REF ORF_type=complete parent_gene=G")
  on.exit(unlink(fasta), add = TRUE)
  elist <- list(genes = data.frame(
    Protein.Group = c("sp|A|B"),
    stringsAsFactors = FALSE
  ))
  classification <- classify_proteins(elist, fasta_path = fasta)
  expect_equal(classification$source, "REF")
  expect_equal(classification$parent_gene, "G")
})

test_that("classify_proteins returns empty df on NULL or no Protein.Group column", {
  expect_equal(nrow(classify_proteins(NULL)), 0)
  expect_equal(nrow(classify_proteins(data.frame(other = 1))), 0)
})

test_that("is_proteogenomic_session detects proteogenomic classes", {
  pg <- data.frame(source = c("REF"), stringsAsFactors = FALSE)
  expect_true(is_proteogenomic_session(pg))
  uniprot_only <- data.frame(source = c("UNIPROT"), stringsAsFactors = FALSE)
  expect_false(is_proteogenomic_session(uniprot_only))
  expect_false(is_proteogenomic_session(NULL))
})

# =============================================================================
# build_proteog_*() — Claude export prompt helpers
# =============================================================================

values_proteog <- list(
  is_proteogenomics = TRUE,
  protein_classification = data.frame(
    Protein.Group = c("a", "b", "c", "d", "e"),
    source        = c("UNIPROT", "REF", "NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT"),
    orf_type      = c(NA, "complete", "complete", "complete", NA),
    parent_gene   = c(NA, "G1", "G2", "G3", NA),
    stringsAsFactors = FALSE
  )
)

test_that("build_proteog_note formats counts correctly", {
  note <- build_proteog_note(values_proteog)
  expect_true(grepl("PROTEOGENOMICS-EXPANDED", note))
  expect_true(grepl("1 canonical UniProt", note))
  expect_true(grepl("1 reference", note))
  expect_true(grepl("1 novel genes", note))
  expect_true(grepl("1 novel isoforms", note))
  expect_true(grepl("1 variant proteoforms", note))
})

test_that("all build_proteog_* helpers return empty string when not proteogenomic", {
  off <- list(is_proteogenomics = FALSE)
  expect_equal(build_proteog_note(off), "")
  expect_equal(build_proteog_file_note(off), "")
  expect_equal(build_proteog_section(off, "brief"), "")
  expect_equal(build_proteog_section(off, "full"), "")
  expect_equal(build_proteog_section(off, "manuscript"), "")
  expect_equal(build_biosynth_proteog_note(off), "")
})

test_that("build_proteog_section returns non-empty text for each template", {
  for (tmpl in c("brief", "full", "manuscript")) {
    s <- build_proteog_section(values_proteog, tmpl)
    expect_gt(nchar(s), 100)
  }
})

test_that("build_proteog_file_note mentions Proteogenomics_Glossary.txt", {
  expect_true(grepl("Proteogenomics_Glossary.txt",
                    build_proteog_file_note(values_proteog)))
})

test_that("build_biosynth_proteog_note mentions follow-up experiment classes", {
  bnote <- build_biosynth_proteog_note(values_proteog)
  expect_true(grepl("NOVEL_GENE", bnote))
  expect_true(grepl("NOVEL_ISOFORM", bnote))
})

# =============================================================================
# helpers_slims: URL validation, project name sanitization
# =============================================================================

test_that("is_slims_url accepts valid SLIMS URLs", {
  expect_true(is_slims_url("http://slimsdata.genomecenter.ucdavis.edu/Data/abc123/Unaligned/"))
  expect_true(is_slims_url("https://slimsdata.genomecenter.ucdavis.edu/Data/i5om268pkp/Unaligned/"))
  expect_false(is_slims_url("https://example.com/Data/abc/"))
  expect_false(is_slims_url(""))
  expect_false(is_slims_url(NA))
})

test_that("sanitize_project_name strips disallowed characters", {
  expect_equal(sanitize_project_name("simple_name"), "simple_name")
  expect_equal(sanitize_project_name("with spaces"), "with_spaces")
  expect_equal(sanitize_project_name("weird/chars$here"), "weird_chars_here")
  expect_equal(sanitize_project_name("dash-and.dot"), "dash-and.dot")
  expect_error(sanitize_project_name(""))
})

test_that("load_reference_registry handles missing file gracefully", {
  old <- Sys.getenv("DELIMP_REFERENCE_REGISTRY", unset = "")
  Sys.setenv(DELIMP_REFERENCE_REGISTRY = "/nonexistent/registry.json")
  on.exit({
    if (nzchar(old)) Sys.setenv(DELIMP_REFERENCE_REGISTRY = old)
    else Sys.unsetenv("DELIMP_REFERENCE_REGISTRY")
  })
  expect_equal(length(load_reference_registry()), 0)
})

# =============================================================================
# helpers_rnaseq: STAR tier selection + sbatch generators
# =============================================================================

test_that("select_star_params returns correct tier for each read length", {
  expect_equal(select_star_params(150)$tier, "default")
  expect_equal(select_star_params(150)$qc_gate_unique_pct, 60)
  expect_equal(select_star_params(100)$tier, "mildly_relaxed")
  expect_equal(select_star_params(100)$qc_gate_unique_pct, 45)
  expect_equal(select_star_params(92)$tier, "significantly_relaxed")
  expect_equal(select_star_params(92)$qc_gate_unique_pct, 25)
  expect_equal(select_star_params(50)$tier, "refuse")

  # Boundaries
  expect_equal(select_star_params(130)$tier, "default")
  expect_equal(select_star_params(129)$tier, "mildly_relaxed")
  expect_equal(select_star_params(99)$tier, "significantly_relaxed")
  expect_equal(select_star_params(60)$tier, "significantly_relaxed")
  expect_equal(select_star_params(59)$tier, "refuse")
})

test_that("detect_read_length samples median from gzipped FASTQ", {
  tmp_fq <- tempfile(fileext = ".fastq.gz")
  on.exit(unlink(tmp_fq), add = TRUE)
  con <- gzfile(tmp_fq, "wt")
  for (i in 1:50) {
    writeLines(c(
      sprintf("@read_%d", i),
      paste(rep("A", 92), collapse = ""),
      "+",
      paste(rep("I", 92), collapse = "")
    ), con)
  }
  close(con)
  expect_equal(detect_read_length(tmp_fq, n_reads = 50L), 92)
})

test_that("generate_fastp_sbatch produces correct array directive and module", {
  script <- generate_fastp_sbatch("/tmp/proj_test", c("S1", "S2", "S3"))
  expect_true(grepl("#!/bin/bash -l", script))
  expect_true(grepl("#SBATCH --array=1-3", script))
  expect_true(grepl("module load fastp/0.23.4", script))
  expect_true(grepl("--detect_adapter_for_pe", script))
})

test_that("generate_rrna_sbatch uses --very-sensitive-local + --un-conc-gz", {
  script <- generate_rrna_sbatch(
    "/tmp/proj_test", c("S1", "S2"),
    rrna_index_prefix = "/q/p/d/references/rrna_index/mm39/rrna_mm39"
  )
  expect_true(grepl("module load bowtie2/2.5.2", script))
  expect_true(grepl("--very-sensitive-local", script))
  expect_true(grepl("--un-conc-gz", script))
})

test_that("generate_star_sbatch propagates tier flags + strand + multi attributes", {
  tp <- select_star_params(92)
  script <- generate_star_sbatch("/tmp/proj_test", c("S1"), "/q/star_index/mm39", tp)
  expect_true(grepl("--outFilterScoreMinOverLread 0.30", script))
  expect_true(grepl("--outFilterMultimapNmax 20", script))
  expect_true(grepl("--outSAMstrandField intronMotif", script))
  expect_true(grepl("--outSAMattributes Standard XS NH", script))
})

test_that("generate_star_sbatch refuses on refuse tier", {
  expect_error(
    generate_star_sbatch("/tmp/proj", "S1", "/idx", select_star_params(50))
  )
})

test_that("generate_qc_gate_sbatch surfaces all 4 likely-cause categories", {
  script <- generate_qc_gate_sbatch("/tmp/proj_test", c("S1", "S2"),
                                     qc_gate_unique_pct = 25)
  expect_true(grepl("THRESHOLD=25", script))
  expect_true(grepl("exit 1", script))
  expect_true(grepl("Wrong reference genome", script))
  expect_true(grepl("Heavy non-rRNA contamination", script))
  expect_true(grepl("Library type unsuited", script))
  expect_true(grepl("sample degradation", script))
})

test_that("generate_stringtie_sbatch uses conda env + validates strand flag", {
  script <- generate_stringtie_sbatch(
    "/tmp/proj_test", c("S1"), "/q/gtf/mm39.gtf", strand_flag = "--rf"
  )
  expect_true(grepl("source activate /quobyte/proteomics-grp/de-limp/envs/proteog_helpers",
                    script))
  expect_true(grepl("module load conda", script))
  expect_true(grepl("--rf", script))
  expect_error(
    generate_stringtie_sbatch("/tmp/proj_test", c("S1"), "/q/gtf", strand_flag = "--bogus")
  )
})

test_that("generate_merge_sbatch + generate_gffcompare_sbatch produce stringtie/gffcompare commands", {
  ms <- generate_merge_sbatch("/tmp/proj_test", "/q/gtf/mm39.gtf", c("S1", "S2"))
  expect_true(grepl("stringtie --merge", ms))
  gs <- generate_gffcompare_sbatch("/tmp/proj_test", "/q/gtf/mm39.gtf")
  expect_true(grepl("gffcompare", gs))
  expect_true(grepl("-o gffcmp", gs))
})

test_that("generate_transdecoder_sbatch handles optional DIAMOND db", {
  td <- generate_transdecoder_sbatch("/tmp/proj_test")
  expect_true(grepl("TransDecoder.LongOrfs", td))
  expect_true(grepl("--single_best_only", td))

  td_d <- generate_transdecoder_sbatch("/tmp/proj_test", diamond_db = "/q/uniprot.dmnd")
  expect_true(grepl("module load diamond/2.1.7", td_d))
  expect_true(grepl("diamond blastp", td_d))
})

test_that("generate_rewrite_sbatch uses v1.0 canonical signature only", {
  rs <- generate_rewrite_sbatch("/tmp/proj_test", "MYPROJ")
  expect_true(grepl("--transdecoder", rs))
  expect_true(grepl("--merged-gtf", rs))
  expect_true(grepl("--gffcompare-tmap", rs))
  expect_true(grepl("--project-tag", rs))
  expect_true(grepl("--output", rs))
  # v0.1 args must NOT appear
  expect_false(grepl("--gtf \"", rs, fixed = TRUE))
  expect_false(grepl("--pep \"", rs, fixed = TRUE))

  # Reject malformed project_tag
  expect_error(generate_rewrite_sbatch("/tmp/proj_test", "MY PROJ WITH SPACES"))

  # Reject removed parameter — no silent acceptance
  expect_error(
    generate_rewrite_sbatch("/tmp/proj_test", "MYPROJ", use_gffcompare_args = TRUE)
  )
})

# =============================================================================
# helpers_proteog_qc: log parsers + gate enforcement
# =============================================================================

.write_synthetic_star_log <- function(unique_pct = 26.31, multi_pct = 73.69) {
  log_file <- tempfile()
  writeLines(c(
    sprintf("                          Number of input reads |\t4380211"),
    sprintf("                      Average input read length |\t183"),
    sprintf("                   Uniquely mapped reads number |\t1152374"),
    sprintf("                        Uniquely mapped reads %% |\t%.2f%%", unique_pct),
    sprintf("                          Average mapped length |\t93.18"),
    sprintf("                       Number of splices: Total |\t651848"),
    sprintf("            Number of splices: Annotated (sjdb) |\t650173"),
    sprintf("                      Mismatch rate per base, %% |\t0.15%%"),
    sprintf("        Number of reads mapped to multiple loci |\t3227584"),
    sprintf("             %% of reads mapped to multiple loci |\t%.2f%%", multi_pct),
    sprintf("        Number of reads mapped to too many loci |\t148"),
    sprintf("             %% of reads mapped to too many loci |\t0.00%%"),
    sprintf("  Number of reads unmapped: too many mismatches |\t0"),
    sprintf("       %% of reads unmapped: too many mismatches |\t0.00%%"),
    sprintf("            Number of reads unmapped: too short |\t4"),
    sprintf("                 %% of reads unmapped: too short |\t0.00%%"),
    sprintf("                Number of reads unmapped: other |\t101"),
    sprintf("                     %% of reads unmapped: other |\t0.00%%")
  ), log_file)
  log_file
}

test_that("parse_star_log extracts metrics from synthetic STAR log", {
  log_file <- .write_synthetic_star_log()
  on.exit(unlink(log_file), add = TRUE)
  parsed <- parse_star_log(log_file)
  expect_equal(parsed$unique_pct, 26.31)
  expect_equal(parsed$multi_pct, 73.69)
  expect_equal(parsed$n_input_reads, 4380211)
  expect_equal(parsed$avg_input_read_length, 183)
  expect_equal(parsed$mismatch_rate_pct, 0.15)
  expect_lt(abs(parsed$pct_annotated_splices - 99.74), 0.1)
})

test_that("check_alignment_quality applies tier gate correctly", {
  log_file <- .write_synthetic_star_log(unique_pct = 26.31)
  on.exit(unlink(log_file), add = TRUE)

  # 92bp tier (gate=25) should PASS at 26.31% unique
  gate_pass <- check_alignment_quality(log_file, select_star_params(92))
  expect_true(gate_pass$pass)

  # Default tier (gate=60) should FAIL at 26.31% unique
  gate_fail <- check_alignment_quality(log_file, select_star_params(150))
  expect_false(gate_fail$pass)
  expect_true(grepl("Alignment quality below threshold", gate_fail$message))
})

test_that("render_gate_failure includes all 4 candidate causes", {
  log_file <- .write_synthetic_star_log(unique_pct = 26.31)
  on.exit(unlink(log_file), add = TRUE)
  gate_fail <- check_alignment_quality(log_file, select_star_params(150))
  rendered <- render_gate_failure(gate_fail)
  expect_true(grepl("QC GATE FAILED", rendered))
  expect_true(grepl("Wrong reference genome", rendered))
  expect_true(grepl("Heavy non-rRNA contamination", rendered))
})

test_that("parse_rrna_log extracts overall alignment rate from bowtie2 output", {
  rrna_log <- tempfile()
  on.exit(unlink(rrna_log), add = TRUE)
  writeLines(c(
    "4396135 reads; of these:",
    "  4396135 (100.00%) were paired; of these:",
    "    4380211 (99.64%) aligned concordantly 0 times",
    "    0 (0.00%) aligned concordantly exactly 1 time",
    "    15924 (0.36%) aligned concordantly >1 times",
    "10.57% overall alignment rate"
  ), rrna_log)
  r <- parse_rrna_log(rrna_log)
  expect_equal(r$overall_alignment_pct, 10.57)
  expect_equal(r$n_input_pairs, 4396135)
  expect_equal(r$concordant_pairs_aligned, 15924)
})

test_that("parse_fastp_json extracts read-length + pct_passed", {
  fastp_json <- tempfile(fileext = ".json")
  on.exit(unlink(fastp_json), add = TRUE)
  writeLines(
    '{"summary":{"before_filtering":{"total_reads":5000000,"read1_mean_length":150},"after_filtering":{"total_reads":4900000,"read1_mean_length":148}}}',
    fastp_json
  )
  fp <- parse_fastp_json(fastp_json)
  expect_equal(fp$total_reads_pre, 5000000)
  expect_equal(fp$total_reads_post, 4900000)
  expect_equal(fp$median_read_length_pre, 150)
  expect_equal(fp$median_read_length_post, 148)
  expect_lt(abs(fp$pct_passed - 98), 0.5)
})

test_that("check_pipeline_gates aggregates per-sample pass/fail", {
  log_pass <- tempfile()
  log_fail <- tempfile()
  on.exit(unlink(c(log_pass, log_fail)), add = TRUE)
  writeLines("                        Uniquely mapped reads % |\t80.0%", log_pass)
  writeLines("                        Uniquely mapped reads % |\t40.0%", log_fail)
  gates <- check_pipeline_gates(c(log_pass, log_fail), select_star_params(150))
  expect_false(gates$pass)  # one of two failed
  expect_true(gates$sample_results[[1]]$pass)
  expect_false(gates$sample_results[[2]]$pass)
})
