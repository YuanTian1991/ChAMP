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
    mes2 <- "ChAMP provides comprehensive integrated analysis pipeline for DNA methylation HumanMethylation Beadchip.\n\n You may use vignette(\"ChAMP\") to view html version guildbook.\n\n If you have any question or suggestion about ChAMP, please email to champ450k@gmail.com."
    packageStartupMessage("Package loaded\n",mes,mes2)
}

