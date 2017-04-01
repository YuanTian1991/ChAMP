\name{champnews}
\title{ChAMP News}
\encoding{UTF-8}
\section{Version 2.6.2}{
\itemize{
    \item DMRcate pacakge get updated, Error like "Error in if (nsig == 0) { : missing value where TRUE/FALSE needed" has been solved.
    \item In champ.load(), instead of replacing all 0 and negative value into 0.0001, we relplace them as smallest positive now.
    \item Fixed warnings() in GUI() functions.
    \item In champ.runCombat() function, removed restriction on factors like Sample_Group. Also, added "variable" parameter so that user may assign other variables other then "Sample_Group".
    \item Modified champ.DMR() function, for ProbeLasso, there is no need to input myDMP anymore, ProbeLasso function would calculate inside the function.
    }
}
