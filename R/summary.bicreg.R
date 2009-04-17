`summary.bicreg` <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, display.dropped = FALSE, ...) 
{
    x <- object
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (display.dropped & x$reduced) {
        cat("\nThe following variables were dropped prior to averaging:\n")
        cat(x$dropped)
        cat("\n")
    }
    n.models <- min(n.models, x$n.models)
    sel <- 1:n.models
    cat("\n ", length(x$postprob), " models were selected")
    cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(x$postprob[sel]), digits), "): \n\n")
    nms <- length(x$namesx) + 1
    r2 <- format(round(x$r2[sel]/100, 3), digits = 3)
    nvar <- rbind(rep(1, length(x$namesx) + 1)) %*% t(x$ols[sel, , drop = FALSE
        ] != 0) - 1
    modelposts <- format(round(x$postprob[sel], 3), digits = 3)
    coeffs <- t(x$ols[sel, ,drop=FALSE])
    cfbic <- rbind(x$bic[sel], coeffs)
    cfbicf <- format(cfbic, digits = digits)
    coeffsf <- cfbicf[-1, ,drop=FALSE]
    bic <- cfbicf[1, ]
    dotoffset <- round(max(nchar(coeffsf))/2)
    zerocoefstring <- paste(paste(rep(" ", times = dotoffset), 
        collapse = "", sep = ""), ".", sep = "")
    coeffsf[coeffs == 0] <- zerocoefstring
    postmeans <- format(x$postmean, digits = digits)
    postsds <- format(x$postsd, digits = digits)
    if (conditional) {
        cpostmeans <- format(x$condpostmean, digits = digits)
        cpostsds <- format(x$condpostsd, digits = digits)
    }
    varposts <- format(round(c(100, x$probne0), 1), digits = 3)
    strlength <- nchar(coeffsf[1, 1])
    decpos <- nchar(unlist(strsplit(coeffsf[1, 1], "\\."))[1])
    offset <- paste(rep(" ", times = decpos - 1), sep = "", collapse = "")
    offset2 <- paste(rep(" ", times = decpos + 1), sep = "", 
        collapse = "")
    r2 <- paste(offset, r2, sep = "")
    modelposts <- paste(offset, modelposts, sep = "")
    nvar <- paste(offset2, nvar, sep = "")
    top <- cbind(varposts, postmeans, postsds)
    if (conditional) 
        top <- cbind(top, cpostmeans, cpostsds)
    top <- cbind(top, coeffsf)
    linesep <- rep("", times = ncol(top))
    offset <- c("", "", "")
    if (conditional) 
        offset <- c(offset, c("", ""))
    bottom <- rbind(c(offset, nvar), c(offset, r2), c(offset, 
        bic), c(offset, modelposts))
    out <- rbind(top, linesep, bottom)
    row.names(out) <- c("Intercept", x$namesx, "", "nVar", "r2", 
        "BIC", "post prob")
    colnms <- c("p!=0", " EV", "SD")
    if (conditional) 
        colnms <- c(colnms, "cond EV", "cond SD")
    colnms <- c(colnms, paste("model ", 1:n.models, sep = ""))
    dimnames(out)[[2]] <- colnms
    print.default(out, print.gap = 2, quote = FALSE, ...)
}

