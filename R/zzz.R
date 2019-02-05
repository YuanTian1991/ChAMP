# ==========================================================================
# ChAMP package initialization
# ==========================================================================

.onAttach <- function(libname, pkgname) {
    mes <- "       ___ _      _   __  __ ___ 
      / __| |_   /_\\ |  \\/  | _ \\
     | (__| ' \\ / _ \\| |\\/| |  _/
      \\___|_||_/_/ \\_\\_|  |_|_|  
      ------------------------------
"
    mes2 <- "    If you have any question or suggestion about ChAMP, please email to champ450k@gmail.com or tianyuan1991hit@gmail.com.\n    ChAMP paper just get published on Bioconductor! ^_^\n\n    Yuan Tian, Tiffany J Morris, Amy P Webster, Zhen Yang, Stephan Beck, Andrew Feber, Andrew E Teschendorff; ChAMP: updated methylation analysis pipeline for Illumina BeadChips, Bioinformatics, btx513, https://doi.org/10.1093/bioinformatics/btx513\n     --------------------------\n"
    packageStartupMessage(">> Package version 2.12.4 loaded <<\n",mes,mes2)
}

