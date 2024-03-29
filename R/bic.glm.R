"bic.glm" <-
function (x, ...) 
UseMethod("bic.glm")

"bic.glm.data.frame" <-
function (x, y, glm.family, wt = rep(1, nrow(x)), strict = FALSE, 
    prior.param = c(rep(0.5, ncol(x))), OR = 20, maxCol = 30, 
    OR.fix = 2, nbest = 150, dispersion = NULL, factor.type = TRUE, 
    factor.prior.adjust = FALSE, occam.window = TRUE, call = NULL, 
    ...) 
{

    namx <- names(x)

    if (is.null(colnames(x))) colnames(x) <- 1:ncol(x)

# this will add a suffix ".x" to the names to prevent duplicate names when
# there are factors

    LAST <- sapply(colnames(x), function(z) substring(z,nchar(z)))
    m <- match(LAST,as.character(0:9),nomatch=0)
    pad <- ""
    if (any(m != 0)) {
      pad <- ".x"
      colnames(x) <- paste( colnames(x), ".x", sep = "")
    }

    facx <- sapply(x,is.factor)

    leaps.glm <- function(info, coef, names.arg, nbest = nbest) {
        names.arg <- names.arg
        if (is.null(names.arg)) 
            names.arg <- c(as.character(1:9), LETTERS, letters)[1:ncol(info)]
        if (length(names.arg) < ncol(info)) 
            stop("Too few names")
        bIb <- coef %*% info %*% coef
        kx <- ncol(info)
        maxreg <- nbest * kx
        if (kx < 3) 
            stop("Too few independent variables")
        imeth <- 1
        df <- kx + 1
        Ib <- info %*% coef
        rr <- cbind(info, Ib)
        rr <- rbind(rr, c(Ib, bIb))
        it <- 0
        n.cols <- kx + 1
        nv <- kx + 1
        nf <- 0
        no <- 1e+05
        ib <- 1
        mb <- nbest
        nd <- n.cols
        nc <- 4 * n.cols
        rt <- matrix(rep(0, times = nd * nc), ncol = nc)
        rt[, 1:n.cols] <- rr
        iw <- c(1:(kx + 1), rep(0, times = 4 * nd))
        nw <- length(iw)
        rw <- rep(0, times = 2 * mb * kx + 7 * nd)
        nr <- length(rw)
        t1 <- 2
        s2 <- -1
        ne <- 0
        iv <- 0
        nret <- mb * kx
        Subss <- rep(0, times = nret)
        RSS <- Subss
        ans <- .Fortran("fwleaps", as.integer(nv), as.integer(it), 
            as.integer(kx), as.integer(nf), as.integer(no), as.integer(1), 
            as.double(2), as.integer(mb), as.double(rt), as.integer(nd), 
            as.integer(nc), as.integer(iw), as.integer(nw), as.double(rw), 
            as.integer(nr), as.double(t1), as.double(s2), as.integer(ne), 
            as.integer(iv), as.double(Subss), as.double(RSS), 
            as.integer(nret), PACKAGE = "BMA")
        regid <- ans[[21]]/2
        r2 <- ans[[20]]
        nreg <- sum(regid > 0)
        regid <- regid[1:nreg]
        r2 <- r2[1:nreg]
        which <- matrix(TRUE, nreg, kx)
        z <- regid
        which <- matrix(as.logical((rep.int(z, kx)%/%rep.int(2^((kx - 
            1):0), rep.int(length(z), kx)))%%2), byrow = FALSE, 
            ncol = kx)
        size <- which %*% rep(1, kx)
        label <- character(nreg)
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        for (i in 1:nreg) label[i] <- paste(names.arg[which[i, 
            ]], collapse = sep)
        ans <- list(r2 = r2, size = size, label = label, which = which)
        return(ans)
    }

    factor.names <- function(x) {
        out <- list()
        for (i in 1:ncol(x)) if (is.factor(x[, i])) 
            out[[i]] <- levels(x[, i])
        else out <- c(out, list(NULL))
        attributes(out)$names <- names(x)
        return(out)
    }

    create.assign <- function(xx) {
        asgn <- list()
        asgn[[1]] <- 1
        cnt <- 2
        for (i in 1:ncol(xx)) {
            if (!is.factor(xx[, i])) 
                size <- 1
            else size <- length(levels(xx[, i])) - 1
            asgn[[i + 1]] <- cnt:(cnt + size - 1)
            cnt <- cnt + size
        }
        names(asgn) <- c("(Intercept)", attributes(xx)$names)
        return(asgn)
    }

    dropcols <- function(x, y, glm.family, wt, maxCols = 30) {
        vnames <- attributes(x)$names
        nvar <- length(vnames)
        isfac <- rep(FALSE, times = nvar)
        for (i in 1:nvar) isfac[i] <- is.factor(x[, i])
        nlevels <- rep(NA, times = nvar)
        for (i in 1:nvar) if (isfac[i]) 
            nlevels[i] <- length(levels(x[, i]))
        any.dropped <- FALSE
        mm <- model.matrix(terms.formula(~., data = x), data = x)
        designx <- attributes(mm)$assign
        n.designx <- length(designx)
        designx.levels <- rep(1, times = n.designx)
        for (i in 2:n.designx) if (isfac[designx[i]]) 
            designx.levels[i] <- sum(designx[1:i] == designx[i]) + 
               1
        x.df <- data.frame(x = x)
        glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
            data = x.df)
        glm.assign <- create.assign(x)
        while (length(glm.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            dropglm <- drop1(glm.out, test = "Chisq")
            lrtname <- if(!is.null(dropglm$LRT)) "LRT" else "scaled dev."
            dropped <- which.max(dropglm[[lrtname]][-1]) + 1
            if (length(dropped) == 0) 
                stop("dropped == 0")
            x.df <- x.df[, -(dropped - 1)]
            designx.levels <- designx.levels[-dropped]
            designx <- designx[-dropped]
            glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
                data = x.df)
        }
        remaining.vars <- unique(designx[-1])
        new.nvar <- length(remaining.vars)
        dropped.vars <- vnames[-remaining.vars]
        dropped.levels <- NULL
        ncol.glm <- ncol(x.df) - 1
        xx <- data.frame(matrix(rep(NA, times = new.nvar * nrow(x.df)), 
            ncol = new.nvar))
        colnames(xx) <- sapply(colnames(x.df), 
                               function(s) substring(s,3,nchar(s)))
        x.df <- x.df[-(ncol.glm + 1)]
        new.names = rep(NA, times = new.nvar)
        for (i in 1:new.nvar) {
            cvar <- remaining.vars[i]
            lvls <- designx.levels[cvar == designx]
            if (isfac[cvar]) {
                if (length(lvls) != length(levels(x[, cvar]))) {
                  newvar <- (as.matrix(x.df[, cvar == designx[-1]]) %*% 
                    cbind(lvls - 1)) + 1
                  xx[, i] <- factor(levels(x[, cvar])[newvar])
                  new.names[i] <- vnames[cvar]
                  removed.levels <- levels(x[, cvar])[-c(1, lvls)]
                  dropped.levels <- c(dropped.levels, paste(vnames[cvar], 
                    "_", removed.levels, sep = ""))
                }
                else {
                  xx[, i] <- factor(x[, cvar])
                  new.names[i] <- vnames[cvar]
                }
            }
            else {
                xx[, i] <- x[, cvar]
                new.names[i] <- vnames[cvar]
            }
        }
        dropped <- c(dropped.vars, dropped.levels)
        return(list(mm = xx, any.dropped = any.dropped, dropped = dropped, 
            var.names = new.names, remaining.vars = remaining.vars))
    }

    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    options(contrasts = c("contr.treatment", "contr.treatment"))
    prior.weight.denom <- 0.5^ncol(x)
    x <- as.data.frame(x)
    LEVELS <- lapply(x, levels)
    names.arg <- names(x)
    if (is.null(names.arg)) 
        names.arg <- paste("X", 1:ncol(x), sep = "")
    x2 <- na.omit(x)
    used <- match(row.names(x), row.names(x2))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
       wt <- wt[-omitted]
       x <- x2
       y <- y[-omitted]
       warning(paste("There were ", length(omitted), 
                     "records deleted due to NA's"))
    }
    leaps.x <- x
    output.names <- names(x)
    fn <- factor.names(x)
    factors <- !all(unlist(lapply(fn, is.null)))
    x.df <- data.frame(x = x)
    glm.out <- glm(y ~ ., family = glm.family, weights = wt, data = x.df)
    glm.assign <- create.assign(x)
    fac.levels <- unlist(lapply(glm.assign, length)[-1])
    varNames <- names.arg
    if (factors) {
# factors get turned into vectors in leaps.x using coefficient values
        cdf <- cbind.data.frame(y = y, x)
        ncoly <- if (is.null(dim(y))) 
            1
        else ncol(y)
        mm <- model.matrix(formula(cdf), data = cdf)[, -(1:ncoly), 
            drop = FALSE]
        varNames <- colnames(mm)
        mmm <- data.frame(matrix(mm, nrow = nrow(mm), byrow = FALSE))
        names(mmm) <- dimnames(mm)[[2]]
        output.names <- names(mmm)
        if (factor.type) { 
           for (i in 1:length(names(x))) {
                if (!is.null(fn[[i]])) {
                  nx <- names(x)[i]
                  coefs <- glm.out$coef[glm.assign[[i + 1]]]
                  old.vals <- x[, i]
                  new.vals <- c(0, coefs)
                  new.vec <- as.vector(new.vals[match(old.vals, 
                    fn[[i]])])
                  leaps.x[, nx] <- new.vec
                }
            }
        }
        else {
            new.prior <- NULL
            for (i in 1:length(names(x))) {
                addprior <- prior.param[i]
                if (!is.null(fn[[i]])) {
                  k <- length(fn[[i]])
                  if (factor.prior.adjust) 
                    addprior <- rep(1 - (1 - prior.param[i])^(1/(k - 
                      1)), k - 1)
                  else addprior <- rep(prior.param[i], k - 1)
                }
                new.prior <- c(new.prior, addprior)
            }
            prior.param <- new.prior
            x <- leaps.x <- mmm
        }
    }
    xx <- data.frame()
    xx <- dropcols(leaps.x, y, glm.family, wt, maxCol)
    var.names <- xx$var.names
    remaining <- xx$remaining.vars
    leaps.x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- 0
    if (reduced) {
        dropped <- pmatch(xx$dropped, varNames, nomatch = 0)
        varNames <- varNames[-dropped]
        prior.param <- prior.param[-dropped]
    }
    nvar <- length(x[1, ])
    x <- x[, remaining, drop = FALSE]
    x <- data.frame(x)
    fac.levels <- fac.levels[remaining]
    output.names <- list()
    for (i in 1:length(var.names)) {
        if (is.factor(x[, i])) 
            output.names[[i]] <- levels(x[, i])
        else output.names[[i]] <- NA
    }
    xnames <- names(x)
    names(leaps.x) <- var.names
    x.df <- data.frame(x = leaps.x)
    glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
        data = x.df, x = TRUE)
#   glm.assign <- create.assign(leaps.x)
    glm.assign <- create.assign(x)
    if (factor.type == FALSE) 
        fac.levels <- unlist(lapply(glm.assign, length)[-1])
    famname <- glm.out$family["family"]$family
    linkinv <- glm.out$family["linkinv"]$linkinv
    if (is.null(dispersion)) {
        if (famname == "poisson" | famname == "binomial") 
            dispersion <- FALSE
        else dispersion <- TRUE
    }
    nobs <- length(y)
    resid <- resid(glm.out, "pearson")
    rdf <- glm.out$df.resid
    is.wt <- !all(wt == rep(1, nrow(x)))
    if (is.wt) {
        resid <- resid * sqrt(wt)
        excl <- wt == 0
        if (any(excl)) {
            warning(paste(sum(excl), "rows with zero wts not counted"))
            resid <- resid[!excl]
        }
    }
    phihat <- sum(resid^2)/rdf
    if (dispersion) 
        disp <- phihat
    else disp <- 1
    coef <- glm.out$coef[-1]
    p <- glm.out$rank
    R <- glm.out$R
    rinv <- diag(p)
    rinv <- backsolve(R, rinv)
    rowlen <- drop(((rinv^2) %*% rep(1, p))^0.5)
    sigx <- rowlen %o% sqrt(disp)
    correl <- rinv %*% t(rinv) * outer(1/rowlen, 1/rowlen)
    cov <- correl * sigx %*% t(sigx)
    info <- solve(cov[-1, -1])
    if (ncol(x) > 2) {
        a <- leaps.glm(info, coef, names.arg = names(leaps.x), 
            nbest = nbest)
        a$r2 <- pmin(pmax(0, a$r2), 0.999)
        a$r2 <- c(0, a$r2)
        a$size <- c(0, a$size)
        a$label <- c("NULL", a$label)
        a$which <- rbind(rep(FALSE, ncol(x)), a$which)
        nmod <- length(a$size)
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(a$which*prior.mat + (!a$which)*(1 - prior.mat), 1, prod)
        bIb <- as.numeric(coef %*% info %*% coef)
        lrt <- bIb - (a$r2 * bIb)
        bic <- lrt + (a$size) * log(nobs) - 2 * log(prior)
        occam <- bic - min(bic) < 2 * OR.fix * log(OR)
        size <- a$size[occam]
        label <- a$label[occam]
        which <- a$which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
    }
    else {
        nmod <- switch(ncol(x), 2, 4)
        bic <- label <- rep(0, nmod)
        which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)[1:(nmod*nmod/2)], 
                          nmod, nmod/2)
        size <- c(0, 1, 1, 2)[1:nmod]
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(which * prior.mat + (!which) * (1 - prior.mat), 1, prod)
        for (k in 1:nmod) {
            if (k == 1) 
                label[k] <- "NULL"
            else label[k] <- paste(names.arg[which[k, ]], collapse = sep)
        }
    }
    nmod <- length(label)
    model.fits <- as.list(rep(0, nmod))
    dev <- rep(0, nmod)
    df <- rep(0, nmod)
    xdf <- x.df
    colnames(xdf) <- sapply(colnames(xdf), function(s) substring(s,3,nchar(s)))
    for (k in 1:nmod) {
        if (sum(which[k, ]) == 0) {
            glm.out <- glm(y ~ 1, family = glm.family, weights = wt)
        }
        else {
            x.df <- data.frame(x = x[, which[k, ], drop = F])
            if (ncol(x.df) != 1) {
              colnames(x.df) <- sapply(colnames(x.df), 
                                  function(s) substring(s,3,nchar(s)))
            }
            glm.out <- glm(y ~ ., data = x.df, family = glm.family, weights = wt)
        }
        dev[k] <- glm.out$deviance
        df[k] <- glm.out$df.residual
        model.fits[[k]] <- matrix(0, nrow = length(glm.out$coef), 
            ncol = 2, dimnames = list(names(glm.out$coef), NULL))
        model.fits[[k]][, 1] <- glm.out$coef
        coef <- glm.out$coef
        p <- glm.out$rank
        R <- glm.out$R
        rinv <- diag(p)
        rinv <- backsolve(R, rinv)
        rowlen <- drop(((rinv^2) %*% rep(1, p))^0.5)
        sigx <- rowlen %o% sqrt(disp)
        correl <- rinv %*% t(rinv) * outer(1/rowlen, 1/rowlen)
        cov <- correl * sigx %*% t(sigx)
        model.fits[[k]][, 2] <- sqrt(diag(cov))
    }
    bic <- dev/disp - df * log(nobs) - 2 * log(prior)
    if (occam.window) 
        occam <- bic - min(bic) < 2 * log(OR)
    else occam = rep(TRUE, length(bic))
    dev <- dev[occam]
    df <- df[occam]
    size <- size[occam]
    label <- label[occam]
    which <- which[occam, , drop = FALSE]
    bic <- bic[occam]
    prior <- prior[occam]
    model.fits <- model.fits[occam]
    postprob <- exp(-0.5 * (bic - min(bic)))/sum(exp(-0.5 * (bic - min(bic))))
    order.bic <- order(bic, size, label)
    dev <- dev[order.bic]
    df <- df[order.bic]
    size <- size[order.bic]
    label <- label[order.bic]
    which <- which[order.bic, , drop = FALSE]
    bic <- bic[order.bic]
    prior <- prior[order.bic]
    postprob <- postprob[order.bic]
    model.fits <- model.fits[order.bic]
    nmod <- length(bic)
    if (strict & (nmod != 1)) {
        occam <- rep(TRUE, nmod)
        for (k in (2:nmod)) for (j in (1:(k - 1))) {
            which.diff <- which[k, ] - which[j, ]
            if (all(which.diff >= 0)) 
                occam[k] <- FALSE
        }
        dev <- dev[occam]
        df <- df[occam]
        size <- size[occam]
        label <- label[occam]
        which <- which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
        postprob <- postprob[occam]
        postprob <- postprob/sum(postprob)
        model.fits <- model.fits[occam]
    }
    bic <- bic + 2 * log(prior)
    probne0 <- round(100 * t(which) %*% as.matrix(postprob), 1)
    nmod <- length(bic)
    nvar <- max(unlist(glm.assign))
    Ebi <- SDbi <- rep(0, nvar)
    EbiMk <- sebiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    for (i in (1:ncol(x))) {
        whereisit <- glm.assign[[i + 1]]
        if (any(which[, i])) 
            for (k in (1:nmod)) if (which[k, i] == TRUE) {
                spot <- sum(which[k, (1:i)])
                posMk <- (c(0, cumsum(fac.levels[which[k, ]])) + 1)[spot]
                posMk <- posMk:(posMk + fac.levels[i] - 1) + 1
                EbiMk[k, whereisit] <- model.fits[[k]][posMk, 1]
                sebiMk[k, whereisit] <- model.fits[[k]][posMk, 2]
            }
    }
    for (k in 1:nmod) {
        EbiMk[k, 1] <- model.fits[[k]][1, 1]
        sebiMk[k, 1] <- model.fits[[k]][1, 2]
    }
    Ebi <- postprob %*% EbiMk
    Ebimat <- matrix(rep(Ebi, nmod), nrow = nmod, byrow = TRUE)
    SDbi <- sqrt(postprob %*% (sebiMk^2) + postprob %*% ((EbiMk - Ebimat)^2))
    CEbi <- CSDbi <- rep(0, nvar)
    for (i in (1:ncol(x))) {
        sel <- which[, i]
        if (sum(sel) > 0) {
            cpp <- rbind(postprob[sel]/sum(postprob[sel]))
            CEbi[glm.assign[[i + 1]]] <- as.numeric(cpp %*% EbiMk[sel, 
                glm.assign[[i + 1]]])
            CSDbi[glm.assign[[i + 1]]] <- sqrt(cpp %*% (sebiMk[sel, 
                glm.assign[[i + 1]]]^2) + cpp %*% ((EbiMk[sel, 
                glm.assign[[i + 1]]] - CEbi[glm.assign[[i + 1]]])^2))
        }
    }
    CSDbi[1] <- SDbi[1]
    CEbi[1] <- Ebi[1]
    names(output.names) <- var.names
    postmean <- as.vector(Ebi)
    varNames <- gsub("`", "", varNames)
    if (pad != "") varNames <- sapply( varNames, function(z) substring(z,1,nchar(z)-2))
    names(varNames) <- NULL
    colnames(EbiMk) <- names(postmean) <- c("(Intercept)", varNames)
    names(probne0) <- if (factor.type) {
        if (!all(dropped == 0)) 
            namx[-dropped]
        else namx
    }
    else varNames
    
    result <- list(postprob = postprob, label = label, deviance = dev, 
        size = size, bic = bic, prior.param = prior.param, 
        prior.model.weights = prior/prior.weight.denom, 
        family = famname, linkinv = linkinv, levels = LEVELS, 
        disp = disp, which = which, probne0 = c(probne0), postmean = postmean, 
        postsd = as.vector(SDbi), condpostmean = CEbi, condpostsd = CSDbi, 
        mle = EbiMk, se = sebiMk, namesx = var.names, reduced = reduced, 
        dropped = dropped, call = cl, n.models = length(postprob), 
        n.vars = length(probne0), nests = length(Ebi), 
        output.names = output.names, 
        assign = glm.assign, factor.type = factor.type, design = leaps.x, 
        x = x, y = y)
    class(result) <- "bic.glm"
    result
}

"bic.glm.formula" <-
function (f, data, glm.family, wt = rep(1, nrow(data)), strict = FALSE, 
    prior.param = c(rep(0.5, ncol(x))), OR = 20, maxCol = 30, 
    OR.fix = 2, nbest = 150, dispersion = NULL, factor.type = TRUE, 
    factor.prior.adjust = FALSE, occam.window = TRUE, na.action = na.omit, ...) 
{
    cl <- match.call()
    tms <- terms(f, data = data)
    fmatrix <- attr(tms, "factors")
    tms.order <- attr(tms, "order")
    tms.labels <- attr(tms, "term.labels")
    attr(data, "na.action") <- na.action
    if(!is.null(na.action)) data <- na.action(data)
    mm <- model.matrix(tms, data = data)
############################################################################
## change to facilitate predict 10/2011 CF
#   tms.labels <- colnames(mm) 
#   if (tms.labels[1] == "(Intercept)") tms.labels <- tms.labels[-1]
############################################################################
    assn <- attr(mm, "assign")
    nterms <- max(assn)
    datalist <- eval(attr(tms, "variables"), envir = data)
    nvar <- nrow(fmatrix) - 1
    isvarfac <- rep(NA, times = nvar)
    for (i in 1:nvar) isvarfac[i] <- is.factor(datalist[[i + 
        1]])
    istermfac <- rep(NA, times = nterms)
    for (i in 1:nterms) {
        cterms <- fmatrix[-1, i] == 1
        istermfac[i] <- sum(isvarfac[cterms] == FALSE) == 0
    }
    resp.name <- all.vars(f)[1]
    moddata <- data.frame(rep(NA, times = dim(mm)[1]))
    cnames <- NULL
    for (i in 1:nterms) {
        if (istermfac[i]) {
            if (tms.order[i] == 1) {
                moddata <- cbind(moddata, datalist[[i + 1]])
                cnames <- c(cnames, tms.labels[i])
            }
            else {
                sel <- assn == i
                nlev <- sum(sel)
                newfac.index <- (mm[, sel] %*% cbind(1:nlev)) + 
                  1
                facnames <- c("ref", colnames(mm)[sel])
                newfac <- facnames[newfac.index]
                newfac <- factor(newfac)
                moddata <- cbind(moddata, newfac)
                cnames <- c(cnames, paste(tms.labels[i], "..", 
                  sep = ""))
            }
        }
        else {
            sel <- assn == i
            moddata <- cbind(moddata, mm[, sel])
            cnames <- c(cnames, colnames(mm)[sel])
        }
    }
    moddata <- moddata[, -1, drop = FALSE]
########################################################################
## deleted to facilitate predict 10/2011 CF
##  cnames <- gsub(":", ".", cnames)
##  moddata <- moddata
##  colnames(moddata) <- c(cnames)
    colnames(moddata) <- cnames
########################################################################
# caution: at least x needs to be assigned before the call (don't know why) CF
    y <- datalist[[1]]
    if (!is.null(dim(y))) {
      ncases <- apply( y, 1, sum)
      index <- rep( 1:nrow(y), ncases)
      x <- moddata[index,]
      wt <- wt[index]
      y <- as.vector(apply( y, 1, function(n) c(rep(0,n[1]),rep(1,n[2]))))
    }
    else {
      x <- moddata
    }
    result <- bic.glm(x=x, y=y, 
        glm.family, wt = wt, strict = FALSE, 
        prior.param = prior.param, 
        OR = OR, maxCol = maxCol, OR.fix = OR.fix, nbest = nbest, 
        dispersion = dispersion, factor.type = factor.type, 
        factor.prior.adjust = factor.prior.adjust, 
        occam.window = occam.window, call = cl)
########################################################################
## added to facilitate predict 10/2011 CF
########################################################################
    result$formula <- f
    result$x <- moddata
    result$y <- datalist[[1]]
    result$na.action <- length(attr(data, "na.action"))
    result

}

bic.glm.matrix <-
function (x, y, glm.family, wt = rep(1, nrow(x)), strict = FALSE, 
    prior.param = c(rep(0.5, ncol(x))), OR = 20, maxCol = 30, 
    OR.fix = 2, nbest = 150, dispersion = NULL, factor.type = TRUE, 
    factor.prior.adjust = FALSE, occam.window = TRUE, call = NULL, 
    ...) 
{
    leaps.glm <- function(info, coef, names.arg, nbest = nbest) {
        names.arg <- names.arg
        if (is.null(names.arg)) 
            names.arg <- c(as.character(1:9), LETTERS, letters)[1:ncol(info)]
        if (length(names.arg) < ncol(info)) 
            stop("Too few names")
        bIb <- coef %*% info %*% coef
        kx <- ncol(info)
        maxreg <- nbest * kx
        if (kx < 3) 
            stop("Too few independent variables")
        imeth <- 1
        df <- kx + 1
        Ib <- info %*% coef
        rr <- cbind(info, Ib)
        rr <- rbind(rr, c(Ib, bIb))
        it <- 0
        n.cols <- kx + 1
        nv <- kx + 1
        nf <- 0
        no <- 1e+05
        ib <- 1
        mb <- nbest
        nd <- n.cols
        nc <- 4 * n.cols
        rt <- matrix(rep(0, times = nd * nc), ncol = nc)
        rt[, 1:n.cols] <- rr
        iw <- c(1:(kx + 1), rep(0, times = 4 * nd))
        nw <- length(iw)
        rw <- rep(0, times = 2 * mb * kx + 7 * nd)
        nr <- length(rw)
        t1 <- 2
        s2 <- -1
        ne <- 0
        iv <- 0
        nret <- mb * kx
        Subss <- rep(0, times = nret)
        RSS <- Subss
        ans <- .Fortran("fwleaps", as.integer(nv), as.integer(it), 
            as.integer(kx), as.integer(nf), as.integer(no), as.integer(1), 
            as.double(2), as.integer(mb), as.double(rt), as.integer(nd), 
            as.integer(nc), as.integer(iw), as.integer(nw), as.double(rw), 
            as.integer(nr), as.double(t1), as.double(s2), as.integer(ne), 
            as.integer(iv), as.double(Subss), as.double(RSS), 
            as.integer(nret), PACKAGE = "BMA")
        regid <- ans[[21]]/2
        r2 <- ans[[20]]
        nreg <- sum(regid > 0)
        regid <- regid[1:nreg]
        r2 <- r2[1:nreg]
        which <- matrix(TRUE, nreg, kx)
        z <- regid
        which <- matrix(as.logical((rep.int(z, kx)%/%rep.int(2^((kx - 
            1):0), rep.int(length(z), kx)))%%2), byrow = FALSE, 
            ncol = kx)
        size <- which %*% rep(1, kx)
        label <- character(nreg)
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        for (i in 1:nreg) label[i] <- paste(names.arg[which[i, 
            ]], collapse = sep)
        ans <- list(r2 = r2, size = size, label = label, which = which)
        return(ans)
    }
    factor.names <- function(x) {
        out <- list()
        for (i in 1:ncol(x)) if (is.factor(x[, i])) 
            out[[i]] <- levels(x[, i])
        else out <- c(out, list(NULL))
        attributes(out)$names <- names(x)
        return(out)
    }
    create.assign <- function(xx) {
        asgn <- list()
        asgn[[1]] <- 1
        cnt <- 2
        for (i in 1:ncol(xx)) {
            if (!is.factor(xx[, i])) 
                size <- 1
            else size <- length(levels(xx[, i])) - 1
            asgn[[i + 1]] <- cnt:(cnt + size - 1)
            cnt <- cnt + size
        }
        names(asgn) <- c("(Intercept)", attributes(xx)$names)
        return(asgn)
    }
    dropcols <- function(x, y, glm.family, wt, maxCols = 30) {
        vnames <- attributes(x)$names
        nvar <- length(vnames)
        isfac <- rep(FALSE, times = nvar)
        for (i in 1:nvar) isfac[i] <- is.factor(x[, i])
        nlevels <- rep(NA, times = nvar)
        for (i in 1:nvar) if (isfac[i]) 
            nlevels[i] <- length(levels(x[, i]))
        any.dropped <- FALSE
        mm <- model.matrix(terms.formula(~., data = x), data = x)
        designx <- attributes(mm)$assign
        n.designx <- length(designx)
        designx.levels <- rep(1, times = n.designx)
        for (i in 2:n.designx) if (isfac[designx[i]]) 
            designx.levels[i] <- sum(designx[1:i] == designx[i]) + 
                1
        x.df <- data.frame(x = x)
        glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
            data = x.df)
        glm.assign <- create.assign(x)
        while (length(glm.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            dropglm <- drop1(glm.out, test = "Chisq")
#           dropped <- which.max(dropglm$"Pr(Chi)"[-1]) + 1
            dropped <- which.max(dropglm$LRT[-1]) + 1
            if (length(dropped) == 0) stop("dropped == 0")
            x.df <- x.df[, -(dropped - 1)]
            designx.levels <- designx.levels[-dropped]
            designx <- designx[-dropped]
            glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
                data = x.df)
        }
        remaining.vars <- unique(designx[-1])
        new.nvar <- length(remaining.vars)
        dropped.vars <- vnames[-remaining.vars]
        dropped.levels <- NULL
        ncol.glm <- ncol(x.df) - 1
        x.df <- x.df[-(ncol.glm + 1)]
        xx <- data.frame(matrix(rep(NA, times = new.nvar * nrow(x.df)), 
            ncol = new.nvar))
        new.names = rep(NA, times = new.nvar)
        for (i in 1:new.nvar) {
            cvar <- remaining.vars[i]
            lvls <- designx.levels[cvar == designx]
            if (isfac[cvar]) {
                if (length(lvls) != length(levels(x[, cvar]))) {
                  newvar <- (as.matrix(x.df[, cvar == designx[-1]]) %*% 
                    cbind(lvls - 1)) + 1
                  xx[, i] <- factor(levels(x[, cvar])[newvar])
                  new.names[i] <- vnames[cvar]
                  removed.levels <- levels(x[, cvar])[-c(1, lvls)]
                  dropped.levels <- c(dropped.levels, paste(vnames[cvar], 
                    "_", removed.levels, sep = ""))
                }
                else {
                  xx[, i] <- factor(x[, cvar])
                  new.names[i] <- vnames[cvar]
                }
            }
            else {
                xx[, i] <- x[, cvar]
                new.names[i] <- vnames[cvar]
            }
        }
        dropped <- c(dropped.vars, dropped.levels)
        return(list(mm = xx, any.dropped = any.dropped, dropped = dropped, 
            var.names = new.names, remaining.vars = remaining.vars))
    }
    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    options(contrasts = c("contr.treatment", "contr.treatment"))
    prior.weight.denom <- 0.5^ncol(x)
    x <- data.frame(x)
    names.arg <- names(x)
    if (is.null(names.arg)) 
        names.arg <- paste("X", 1:ncol(x), sep = "")
    x2 <- na.omit(x)
    used <- match(row.names(x), row.names(x2))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
        wt <- wt[-omitted]
        x <- x2
        y <- y[-omitted]
        warning(paste("There were ", length(omitted), "records deleted due to NA's"))
    }
    leaps.x <- x
    output.names <- names(x)
    fn <- factor.names(x)
    factors <- !all(unlist(lapply(fn, is.null)))
    x.df <- data.frame(x = x)
    glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
        data = x.df)
    glm.assign <- create.assign(x)
    fac.levels <- unlist(lapply(glm.assign, length)[-1])
    if (factors) {
        cdf <- cbind.data.frame(y = y, x)
        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
        mmm <- data.frame(matrix(mm, nrow = nrow(mm), byrow = FALSE))
        names(mmm) <- dimnames(mm)[[2]]
        output.names <- names(mmm)
        if (factor.type) {
            for (i in 1:length(names(x))) {
                if (!is.null(fn[[i]])) {
                  nx <- names(x)[i]
                  coefs <- glm.out$coef[glm.assign[[i + 1]]]
                  old.vals <- x[, i]
                  new.vals <- c(0, coefs)
                  new.vec <- as.vector(new.vals[match(old.vals, 
                    fn[[i]])])
                  leaps.x[, nx] <- new.vec
                }
            }
        }
        else {
            new.prior <- NULL
            for (i in 1:length(names(x))) {
                addprior <- prior.param[i]
                if (!is.null(fn[[i]])) {
                  k <- length(fn[[i]])
                  if (factor.prior.adjust) 
                    addprior <- rep(1 - (1 - prior.param[i])^(1/(k - 
                      1)), k - 1)
                  else addprior <- rep(prior.param[i], k - 1)
                }
                new.prior <- c(new.prior, addprior)
            }
            prior.param <- new.prior
            x <- leaps.x <- mmm
        }
    }
    xx <- data.frame()
    xx <- dropcols(leaps.x, y, glm.family, wt, maxCol)
    var.names <- xx$var.names
    remaining <- xx$remaining.vars
    leaps.x <- xx$mm
    reduced <- xx$any.dropped
#   dropped <- NULL
#   if (reduced) 
#       dropped <- xx$dropped
    dropped <- 0
    varNames <- var.names
    if (reduced) {
       dropped <- match(xx$dropped,varNames,nomatch=0)
       varNames <- varNames[-dropped]
    }
    nvar <- length(x[1, ])
    x <- x[, remaining, drop = FALSE]
    x <- data.frame(x)
    fac.levels <- fac.levels[remaining]
    output.names <- list()
    for (i in 1:length(var.names)) {
        if (is.factor(x[, i])) 
            output.names[[i]] <- levels(x[, i])
        else output.names[[i]] <- NA
    }
    xnames <- names(x)
    names(leaps.x) <- var.names
    x.df <- data.frame(x = leaps.x)
    glm.out <- glm(y ~ ., family = glm.family, weights = wt, 
        data = x.df, x = TRUE)
    glm.assign <- create.assign(leaps.x)
    if (factor.type == FALSE) 
        fac.levels <- unlist(lapply(glm.assign, length)[-1])
    famname <- glm.out$family["family"]$family
    linkinv <- glm.out$family["linkinv"]$linkinv
    if (is.null(dispersion)) {
        if (famname == "poisson" | famname == "binomial") 
            dispersion <- FALSE
        else dispersion <- TRUE
    }
    nobs <- length(y)
    resid <- resid(glm.out, "pearson")
    rdf <- glm.out$df.resid
    is.wt <- !all(wt == rep(1, nrow(x)))
    if (is.wt) {
        resid <- resid * sqrt(wt)
        excl <- wt == 0
        if (any(excl)) {
            warning(paste(sum(excl), "rows with zero wts not counted"))
            resid <- resid[!excl]
        }
    }
    phihat <- sum(resid^2)/rdf
    if (dispersion) 
        disp <- phihat
    else disp <- 1
    coef <- glm.out$coef[-1]
    p <- glm.out$rank
    R <- glm.out$R
    rinv <- diag(p)
    rinv <- backsolve(R, rinv)
    rowlen <- drop(((rinv^2) %*% rep(1, p))^0.5)
    sigx <- rowlen %o% sqrt(disp)
    correl <- rinv %*% t(rinv) * outer(1/rowlen, 1/rowlen)
    cov <- correl * sigx %*% t(sigx)
    info <- solve(cov[-1, -1])
    if (ncol(x) > 2) {
        a <- leaps.glm(info, coef, names.arg = names(leaps.x), 
            nbest = nbest)
        a$r2 <- pmin(pmax(0, a$r2), 0.999)
        a$r2 <- c(0, a$r2)
        a$size <- c(0, a$size)
        a$label <- c("NULL", a$label)
        a$which <- rbind(rep(FALSE, ncol(x)), a$which)
        nmod <- length(a$size)
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(leaps.x), 
            byrow = TRUE)
        prior <- apply(a$which * prior.mat + (!a$which) * (1 - 
            prior.mat), 1, prod)
        bIb <- as.numeric(coef %*% info %*% coef)
        lrt <- bIb - (a$r2 * bIb)
        bic <- lrt + (a$size) * log(nobs) - 2 * log(prior)
        occam <- bic - min(bic) < 2 * OR.fix * log(OR)
        size <- a$size[occam]
        label <- a$label[occam]
        which <- a$which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
    }
    else {
        nmod <- switch(ncol(x), 2, 4)
        bic <- label <- rep(0, nmod)
        which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
            TRUE, TRUE)[1:(nmod*nmod/2)], nmod, nmod/2)
        size <- c(0, 1, 1, 2)[1:nmod]
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(which * prior.mat + (!which) * (1 - prior.mat), 
            1, prod)
        for (k in 1:nmod) {
            if (k == 1) 
                label[k] <- "NULL"
            else label[k] <- paste(names.arg[which[k, ]], collapse = sep)
        }
    }
    nmod <- length(label)
    model.fits <- as.list(rep(0, nmod))
    dev <- rep(0, nmod)
    df <- rep(0, nmod)
    for (k in 1:nmod) {
        if (sum(which[k, ]) == 0) {
            glm.out <- glm(y ~ 1, family = glm.family, weights = wt)
        }
        else {
            x.df <- data.frame(x = x[, which[k, ]])
            glm.out <- glm(y ~ ., data = x.df, family = glm.family, 
                weights = wt)
        }
        dev[k] <- glm.out$deviance
        df[k] <- glm.out$df.residual
        model.fits[[k]] <- matrix(0, nrow = length(glm.out$coef), 
            ncol = 2)
        model.fits[[k]][, 1] <- glm.out$coef
        coef <- glm.out$coef
        p <- glm.out$rank
        R <- glm.out$R
        rinv <- diag(p)
        rinv <- backsolve(R, rinv)
        rowlen <- drop(((rinv^2) %*% rep(1, p))^0.5)
        sigx <- rowlen %o% sqrt(disp)
        correl <- rinv %*% t(rinv) * outer(1/rowlen, 1/rowlen)
        cov <- correl * sigx %*% t(sigx)
        model.fits[[k]][, 2] <- sqrt(diag(cov))
    }
    bic <- dev/disp - df * log(nobs) - 2 * log(prior)
    if (occam.window) 
        occam <- bic - min(bic) < 2 * log(OR)
    else occam = rep(TRUE, length(bic))
    dev <- dev[occam]
    df <- df[occam]
    size <- size[occam]
    label <- label[occam]
    which <- which[occam, , drop = FALSE]
    bic <- bic[occam]
    prior <- prior[occam]
    model.fits <- model.fits[occam]
    postprob <- exp(-0.5 * (bic - min(bic)))/sum(exp(-0.5 * (bic - 
        min(bic))))
    order.bic <- order(bic, size, label)
    dev <- dev[order.bic]
    df <- df[order.bic]
    size <- size[order.bic]
    label <- label[order.bic]
    which <- which[order.bic, , drop = FALSE]
    bic <- bic[order.bic]
    prior <- prior[order.bic]
    postprob <- postprob[order.bic]
    model.fits <- model.fits[order.bic]
    nmod <- length(bic)
    if (strict & (nmod != 1)) {
        occam <- rep(TRUE, nmod)
        for (k in (2:nmod)) for (j in (1:(k - 1))) {
            which.diff <- which[k, ] - which[j, ]
            if (all(which.diff >= 0)) 
                occam[k] <- FALSE
        }
        dev <- dev[occam]
        df <- df[occam]
        size <- size[occam]
        label <- label[occam]
        which <- which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
        postprob <- postprob[occam]
        postprob <- postprob/sum(postprob)
        model.fits <- model.fits[occam]
    }
    bic <- bic + 2 * log(prior)
    probne0 <- round(100 * t(which) %*% as.matrix(postprob), 
        1)
    nmod <- length(bic)
    nvar <- max(unlist(glm.assign))
    Ebi <- rep(0, nvar)
    SDbi <- rep(0, nvar)
    EbiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    sebiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    for (i in (1:ncol(x))) {
        whereisit <- glm.assign[[i + 1]]
        if (any(which[, i])) 
            for (k in (1:nmod)) if (which[k, i] == TRUE) {
                spot <- sum(which[k, (1:i)])
                posMk <- (c(0, cumsum(fac.levels[which[k, ]])) + 
                  1)[spot]
                posMk <- posMk:(posMk + fac.levels[i] - 1) + 
                  1
                EbiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  1]
                sebiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  2]
            }
    }
    for (k in 1:nmod) {
        EbiMk[k, 1] <- model.fits[[k]][1, 1]
        sebiMk[k, 1] <- model.fits[[k]][1, 2]
    }
    Ebi <- postprob %*% EbiMk
    Ebimat <- matrix(rep(Ebi, nmod), nrow = nmod, byrow = TRUE)
    SDbi <- sqrt(postprob %*% (sebiMk^2) + postprob %*% ((EbiMk - 
        Ebimat)^2))
    CSDbi <- rep(0, nvar)
    CEbi <- CSDbi
    for (i in (1:ncol(x))) {
        sel <- which[, i]
        if (sum(sel) > 0) {
            cpp <- rbind(postprob[sel]/sum(postprob[sel]))
            CEbi[glm.assign[[i + 1]]] <- as.numeric(cpp %*% EbiMk[sel, 
                glm.assign[[i + 1]]])
            CSDbi[glm.assign[[i + 1]]] <- sqrt(cpp %*% (sebiMk[sel, 
                glm.assign[[i + 1]]]^2) + cpp %*% ((EbiMk[sel, 
                glm.assign[[i + 1]]] - CEbi[glm.assign[[i + 1]]])^2))
        }
    }
    CSDbi[1] <- SDbi[1]
    CEbi[1] <- Ebi[1]
    names(output.names) <- var.names
    result <- list(postprob = postprob, label = label, deviance = dev, 
        size = size, bic = bic, prior.param = prior.param, prior.model.weights
                   = prior/prior.weight.denom, 
        family = famname, linkinv = linkinv, 
        disp = disp, which = which, probne0 = c(probne0), 
        postmean = as.vector(Ebi), postsd = as.vector(SDbi), 
        condpostmean = CEbi, condpostsd = CSDbi, mle = EbiMk, 
        se = sebiMk, namesx = var.names, reduced = reduced, dropped = dropped, 
        call = cl, n.models = length(postprob), n.vars = length(probne0), 
        nests = length(Ebi), output.names = output.names, assign = glm.assign, 
        factor.type = factor.type, design = leaps.x, x = x, y = y)
    class(result) <- "bic.glm"
    result
}

