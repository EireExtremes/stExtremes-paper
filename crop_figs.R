##==============================================================================
## Trim the white margins from the generated figures
##
## Every figure is saved at a fixed canvas size, but the map panels use
## coord_sf and so hold their own aspect ratio. The panel is letterboxed
## inside the canvas and the leftover white is baked into the file, which
## is what puts a visible gap between a figure and its caption in the
## rendered paper. Measured on the first run, the worst offenders were
## mostly white:
##
##   compare_qalpha_mean          59%   3900x3300 -> 3856x1365
##   compare_field_sd             58%
##   compare_field_mean           57%
##   compare_rlevel_map           45%
##   compare_qalpha_diff          27%
##   compare_rlevel_diff          27%
##   pp_rlevel_ciwidth_map        25%
##   pp_rlevel_median_map         23%
##   pp_rlevel_quantile_map_50yr  23%
##   data_glyph_map               11%
##   map_study_area                8%
##   map_mesh                      7%
##
## The three-panel comparison maps are a row of wide maps on a nearly
## square canvas, which is why they were two thirds air.
##
## Doing it at source would mean solving for the canvas size that makes
## each panel exactly fill it, which depends on the projection, the
## legend and the axis labels, and would have to be redone whenever any
## of those change. Trimming afterwards is the robust version of the same
## thing, and it is what pdfcrop exists for.
##
## PNG (what the manuscript actually includes) is trimmed with magick;
## PDF (the vector twin, for submission) with pdfcrop. Both then get a
## small uniform margin back, so nothing touches the caption.
##
## IDEMPOTENT BY CONSTRUCTION: trimming removes the uniform border,
## including the margin a previous run added, and the same margin is then
## added back. Running this twice is a no-op, so it is safe to put at the
## end of any pipeline. Verified, not assumed; see the check at the foot.
##
## Usage:
##   Rscript crop_figs.R          # every figure in figs
##   Rscript crop_figs.R map_mesh # just the ones matching
##==============================================================================

here::i_am("crop_figs.R")

suppressPackageStartupMessages({
    library(tidyverse)
    library(magick)
})

figs_dir <- here::here("figs")

## Margin in device pixels for PNG, in PostScript points for PDF. Both
## are about 2 mm at the sizes used here.
png_margin <- 24
pdf_margin <- 6

## Anti-aliased edges mean the outer border is not exactly one colour, so
## the trim needs a little tolerance. 2% is enough to catch near-white
## without eating a pale panel background.
trim_fuzz <- 2

pat <- commandArgs(trailingOnly = TRUE)
keep <- function(f) length(pat) == 0 || any(str_detect(f, fixed(pat)))

pngs <- list.files(figs_dir, "\\.png$", full.names = TRUE)
pngs <- pngs[map_lgl(basename(pngs), keep)]
pdfs <- list.files(figs_dir, "\\.pdf$", full.names = TRUE)
pdfs <- pdfs[map_lgl(basename(pdfs), keep)]

##==============================================================================
## PNG

crop_png <- function(f) {
    img <- image_read(f)
    before <- image_info(img)
    out <- img |>
        image_trim(fuzz = trim_fuzz) |>
        image_border(color = "white",
            geometry = sprintf("%dx%d", png_margin, png_margin))
    image_write(out, f)
    after <- image_info(out)
    tibble(file = basename(f),
        w0 = before$width, h0 = before$height,
        w1 = after$width, h1 = after$height,
        saved = 1 - (after$width * after$height) /
            (before$width * before$height))
}

png_res <- map(pngs, crop_png) |> list_rbind()

##==============================================================================
## PDF
##
## pdfcrop writes to a separate file, so crop to a temporary and move it
## back only on success. A failure must leave the original alone rather
## than truncate it.

## The temporary goes NEXT TO the target, not in tempdir(). file.rename()
## cannot move across filesystems, and /tmp is usually its own, so a
## tempfile() there fails with "Invalid cross-device link" AFTER pdfcrop
## has already reported success. That combination reports a crop that
## never happened, which is worse than failing outright.
page_size <- function(f) {
    x <- suppressWarnings(system2("pdfinfo", shQuote(f), stdout = TRUE,
        stderr = FALSE))
    str_trim(str_remove(x[str_detect(x, "^Page size:")], "^Page size:"))
}

crop_pdf <- function(f) {
    tmp <- file.path(dirname(f), paste0(".crop_", basename(f)))
    on.exit(unlink(tmp), add = TRUE)
    before <- page_size(f)
    m <- paste(rep(pdf_margin, 4), collapse = " ")
    st <- suppressWarnings(system2("pdfcrop",
        c("--margins", shQuote(m), shQuote(f), shQuote(tmp)),
        stdout = FALSE, stderr = FALSE))
    ok <- st == 0 && file.exists(tmp) && file.size(tmp) > 0
    if(ok) ok <- file.rename(tmp, f)
    tibble(file = basename(f), ok = ok, before = before,
        after = if(ok) page_size(f) else NA_character_)
}

pdf_res <- if(nzchar(Sys.which("pdfcrop"))) {
    map(pdfs, crop_pdf) |> list_rbind()
} else {
    warning("pdfcrop not found: the PDF figures were left uncropped. ",
        "It ships with TeX Live.")
    tibble(file = character(), ok = logical(), before = character(),
        after = character())
}

## Report a crop only if the page size actually moved, so a silent no-op
## cannot be counted as a success.
pdf_res <- mutate(pdf_res, changed = ok & !is.na(after) & after != before)

##==============================================================================
## Report

png_res |>
    arrange(desc(saved)) |>
    mutate(saved = sprintf("%.0f%%", 100 * saved)) |>
    print(n = Inf)

message(sprintf("## %d png trimmed, %d of %d pdf cropped (%d unchanged)",
    nrow(png_res), sum(pdf_res$changed), nrow(pdf_res),
    sum(!pdf_res$changed)))
if(nrow(pdf_res) > 0 && any(!pdf_res$changed)) {
    print(filter(pdf_res, !changed))
}

## Idempotence check on the PNGs: a second trim of an already-trimmed
## file must find nothing but the margin we just added.
if(nrow(png_res) > 0) {
    f <- file.path(figs_dir, png_res$file[which.max(png_res$saved)])
    a <- image_info(image_read(f))
    b <- image_info(image_border(image_trim(image_read(f),
        fuzz = trim_fuzz), "white",
        sprintf("%dx%d", png_margin, png_margin)))
    stopifnot(a$width == b$width, a$height == b$height)
    message("## Idempotence checked on ", basename(f))
}
