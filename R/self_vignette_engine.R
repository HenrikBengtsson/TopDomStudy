register_vignette_engine <- function(pkgname) {
  tools::vignetteEngine("self", package = pkgname, pattern = "[.]md$",
    weave = function(file, ...) {
      output <- sprintf("%s.html", tools::file_path_sans_ext(basename(file)))
      md <- readLines(file)

      title <- grep("%\\VignetteIndexEntry{", md, fixed = TRUE, value = TRUE)
      title <- gsub(".*[{](.*)[}].*", "\\1", title)

      md <- grep("%\\\\Vignette", md, invert = TRUE, value = TRUE)

      ## Inject vignette title
      md <- c(sprintf("# %s\n\n", title), md)

      html <- commonmark::markdown_html(md, smart = TRUE, extensions = TRUE,
                                        normalize = FALSE)
      css <- readLines("incl/clean.css", warn = FALSE)

      ## Inject HTML environment
      html <- c("<!DOCTYPE html>",
                "<html lang=\"en\">",
                "<head>",
                sprintf("<title>%s</title>", title),
                "<style>", css, "</style>",
                "</head>",
                "<body>", html,	"</body>",
                "</html>")
      html <- paste(html, collapse = "\n")		
      html <- enc2utf8(html)
      cat(html, file = output)
      cat(html, file = "~/foo.html")
      output
    },

    tangle = function(file, ...) {
      ## As of R 3.3.2, vignette engines must produce tangled output, but as
      ## long as it contains all comments then 'R CMD build' will drop it.
      output <- sprintf("%s.R", tools::file_path_sans_ext(basename(file)))
      cat(sprintf("### This is an R script tangled from %s\n",
                  sQuote(basename(file))), file = output)
      output
    }
  )
}
