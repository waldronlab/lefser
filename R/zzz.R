.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      "\n",
      "Please cite our software 😃 \n \n",
      "Asya Khleborodova, Samuel D Gamboa-Tuz, Marcel Ramos, Nicola Segata, Levi Waldron, Sehyun Oh,",
      "Lefser: Implementation of metagenomic biomarker discovery tool, LEfSe, in R, Bioinformatics, 2024;, btae707 \n \n",
      "DOI:  https://doi.org/10.1093/bioinformatics/btae707 \n \n",
      "Citation Options: https://github.com/waldronlab/lefser \n \n"
    )
  )
  invisible()
}