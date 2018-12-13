## Modified by Thomas Lumley 26 Apr 97
## Added boundary checks and step halving
## Modified detection of fitted 0/1 in binomial
## Updated by KH as suggested by BDR on 1998/06/16
nobs   = NROW(y)

weights = rep(1, nobs)
start = NULL
etastart = NULL
mustart = NULL
offset = rep(0, nobs)
family=poisson(link = "identity")
control = glm.control()
intercept = TRUE
x = cbind(cntAll, PC15, x)
y = yd

glm.fit <- 
    function (x, y, weights = rep(1, nobs), start = NULL,
	      etastart = NULL, mustart = NULL, offset = rep(0, nobs),
	      family = gaussian(), control = glm.control(), intercept = TRUE)
{
  x      = as.matrix(x)
  xnames = dimnames(x)[[2L]]
  ynames = if(is.matrix(y)) rownames(y) else names(y)
  conv   = FALSE
  nobs   = NROW(y)
  nvars  = ncol(x)
  EMPTY  = nvars == 0
  
  ## define weights and offset if needed
  if (is.null(weights))
    weights = rep.int(1, nobs)
  if (is.null(offset))
    offset = rep.int(0, nobs)

  ## get family functions:
  variance   = family$variance
  dev.resids = family$dev.resids
  aic     = family$aic
  linkinv = family$linkinv
  mu.eta  = family$mu.eta
  
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object")
  
  valideta = family$valideta
  if (is.null(valideta))
    valideta = function(eta) TRUE
  
  validmu = family$validmu
  if (is.null(validmu))
    validmu = function(mu) TRUE
  
  if(is.null(mustart)) {
      ## calculates mustart and may change y and weights and set n (!)
      eval(family$initialize)
  } else {
      mukeep = mustart
      eval(family$initialize)
      mustart = mukeep
  }
  
  if(EMPTY) {
    eta = rep.int(0, nobs) + offset
    if (!valideta(eta))
        stop("invalid linear predictor values in empty model")
    mu = linkinv(eta)
    ## calculate initial deviance and coefficient
    if (!validmu(mu))
        stop("invalid fitted means in empty model")
    dev = sum(dev.resids(y, mu, weights))
    w = ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals = (y - mu)/mu.eta(eta)
    good = rep(TRUE, length(residuals))
    boundary = conv = TRUE
    coef = numeric(0L)
    iter = 0L
  } else {
    coefold = NULL
    
    if(!is.null(etastart)){ 
      eta = etastart
    }else if(!is.null(start)){
      if (length(start) != nvars){
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", nvars, paste(deparse(xnames), collapse=", ")),
                     domain = NA)
      }else {
                coefold = start
                eta = offset + as.vector(if (NCOL(x) == 1) x * start else x %*% start)
            }
    }else {
      eta = family$linkfun(mustart)
    }
    mu = linkinv(eta)
    
    if (!(validmu(mu) && valideta(eta)))
        stop("cannot find valid starting values: please specify some")
    ## calculate initial deviance and coefficient
    devold = sum(dev.resids(y, mu, weights))
    boundary = conv = FALSE

    ##------------- THE Iteratively Reweighting L.S. iteration -----------
    for (iter in 1L:control$maxit) {
        good = weights > 0
        varmu = variance(mu)[good]
        if (any(is.na(varmu)))
            stop("NAs in V(mu)")
        if (any(varmu == 0))
            stop("0s in V(mu)")
        mu.eta.val = mu.eta(eta)
        if (any(is.na(mu.eta.val[good])))
            stop("NAs in d(mu)/d(eta)")
        ## drop observations for which w will be zero
        good = (weights > 0) & (mu.eta.val != 0)

        if (all(!good)) {
            conv = FALSE
            warning("no observations informative at iteration ", iter)
            break
        }
        z = (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w = sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        ngoodobs = as.integer(nobs - sum(!good))
        ## call Fortran code
        fit = .Fortran("dqrls",
                        qr = x[good, ] * w, n = ngoodobs,
                        p = nvars, y = w * z, ny = 1L,
                        tol = min(1e-7, control$epsilon/1000),
                        coefficients = double(nvars),
                        residuals = double(ngoodobs),
                        effects = double(ngoodobs),
                        rank = integer(1L),
                        pivot = 1L:nvars, qraux = double(nvars),
                        work = double(2 * nvars),
                        PACKAGE = "base")
        if (any(!is.finite(fit$coefficients))) {
            conv = FALSE
            warning("non-finite coefficients at iteration ", iter)
            break
        }
        ## stop if not enough parameters
        if (nobs < fit$rank)
            stop(gettextf("X matrix has rank %d, but only %d observations",
                          fit$rank, nobs), domain = NA)
        ## calculate updated values of eta and mu with the new coef:
        start[fit$pivot] = fit$coefficients
        eta = drop(x %*% start)
        mu  = linkinv(eta = eta + offset)
        dev = sum(dev.resids(y, mu, weights))
        
        if (control$trace)
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        ## check for divergence
        boundary = FALSE
        if (!is.finite(dev)) {
            if(is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated due to divergence", call. = FALSE)
            ii = 1
            while (!is.finite(dev)) {
                if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size")
                ii = ii + 1
                start = (start + coefold)/2
                eta = drop(x %*% start)
                mu = linkinv(eta = eta + offset)
                dev = sum(dev.resids(y, mu, weights))
            }
            boundary = TRUE
            if (control$trace)
                cat("Step halved: new deviance =", dev, "\n")
        }
        ## check for fitted values outside domain.
        if (!(valideta(eta) && validmu(mu))) {
            if(is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
            warning("step size truncated: out of bounds", call. = FALSE)
            ii = 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size")
                ii = ii + 1
                start = (start + coefold)/2
                eta = drop(x %*% start)
                mu = linkinv(eta = eta + offset)
            }
            boundary = TRUE
            dev = sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Step halved: new deviance =", dev, "\n")
        }
        ## check for convergence
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv = TRUE
            coef = start
            break
        } else {
            devold = dev
            coef = coefold = start
        }
    } 
    ##-------------- end IRLS iteration -------------------------------

    if (!conv) warning("algorithm did not converge")
    if (boundary) warning("algorithm stopped at boundary value")
    eps = 10*.Machine$double.eps
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps))
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
        if (any(mu < eps))
            warning("fitted rates numerically 0 occurred")
    }
    ## If X matrix was not full rank then columns were pivoted,
    ## hence we need to re-label the names ...
    ## Original code changed as suggested by BDR---give NA rather
    ## than 0 for non-estimable parameters
    if (fit$rank < nvars) coef[fit$pivot][seq.int(fit$rank+1, nvars)] = NA
    xxnames = xnames[fit$pivot]
    
    ## update by accurate calculation, including 0-weight cases.
    residuals =  (y - mu)/mu.eta(eta)
    
    ## residuals = rep.int(NA, nobs)
    ## residuals[good] = z - (eta - offset)[good] # z does not have offset in.
    fit$qr = as.matrix(fit$qr)
    nr = min(sum(good), nvars)
    if (nr < nvars) {
        Rmat = diag(nvars)
        Rmat[1L:nr, 1L:nvars] = fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat = fit$qr[1L:nvars, 1L:nvars]
    Rmat = as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] = 0
    names(coef) = xnames
    colnames(fit$qr) = xxnames
    dimnames(Rmat) = list(xxnames, xxnames)
  }
  
  names(residuals) = ynames
  names(mu) = ynames
  names(eta) = ynames
  
  # for compatibility with lm, which has a full-length weights vector
  wt = rep.int(0, nobs)
  wt[good] = w^2
  names(wt) = ynames
  names(weights) = ynames
  names(y) = ynames
  if(!EMPTY)
      names(fit$effects) =
          c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
  ## calculate null deviance -- corrected in glm() if offset and intercept
  wtdmu =
    if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev = sum(dev.resids(y, wtdmu, weights))
  ## calculate df
  n.ok = nobs - sum(weights==0)
  nulldf = n.ok - as.integer(intercept)
  rank = if(EMPTY) 0 else fit$rank
  resdf  = n.ok - rank
  
  ## calculate AIC
  aic.model = aic(y, n, mu, weights, dev) + 2*rank
  
  ##     ^^ is only initialize()d for "binomial" [yuck!]
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if(!EMPTY) fit$effects, R = if(!EMPTY) Rmat, rank = rank,
    qr = if(!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot", "tol")], 
                              class="qr"), family = family,
    linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = iter, weights = wt,
    prior.weights = weights, df.residual = resdf, df.null = nulldf,
    y = y, converged = conv, boundary = boundary)
}


# ----------------------------------------------------------------------
# summary.glm 
# ----------------------------------------------------------------------

summary.glm <- function(object, dispersion = NULL,
    correlation = FALSE, symbolic.cor = FALSE, ...)
{
  est.disp = FALSE
  df.r = object$df.residual
  if(is.null(dispersion))	# calculate dispersion if needed
	dispersion =
  if(object$family$family %in% c("poisson", "binomial"))  1
  else if(df.r > 0) {
    est.disp = TRUE
		if(any(object$weights==0))
    warning("observations with zero weight not used for calculating dispersion")
		sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
  } else {
    est.disp = TRUE
    NaN
  }
  
  ## calculate scaled and unscaled covariance matrix
  
  aliased = is.na(coef(object))  # used in print method
  p = object$rank
  if (p > 0) {
    p1 = 1L:p
    Qr = object$qr
    
    ## WATCHIT! does not this rely on pivoting not permuting 1L:p? -- it is quaranteed
    coef.p = object$coefficients[Qr$pivot[p1]]
    covmat.unscaled = chol2inv(Qr$qr[p1,p1,drop=FALSE])
    dimnames(covmat.unscaled) = list(names(coef.p),names(coef.p))
    covmat = dispersion*covmat.unscaled
    var.cf = diag(covmat)
    
    ## calculate coef table
    s.err = sqrt(var.cf)
    tvalue = coef.p/s.err
    
    dn = c("Estimate", "Std. Error")
    if(!est.disp) { # known dispersion
      pvalue = 2*pnorm(-abs(tvalue))
      coef.table = cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) = list(names(coef.p),
                                   c(dn, "z value","Pr(>|z|)"))
    } else if(df.r > 0) {
      pvalue = 2*pt(-abs(tvalue), df.r)
      coef.table = cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) = list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    } else { # df.r == 0
      coef.table = cbind(coef.p, NaN, NaN, NaN)
      dimnames(coef.table) = list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    }
    df.f = NCOL(Qr$qr)
  } else {
    coef.table = matrix(, 0L, 4L)
    dimnames(coef.table) =
    list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    covmat.unscaled = covmat = matrix(, 0L, 0L)
    df.f = length(aliased)
  }
  ## return answer
  
  ## these need not all exist, e.g. na.action.
  keep = match(c("call","terms","family","deviance", "aic",
                  "contrasts", "df.residual","null.deviance","df.null",
                  "iter", "na.action"), names(object), 0L)
  ans = c(object[keep],
           list(deviance.resid = residuals(object, type = "deviance"),
                coefficients = coef.table,
                aliased = aliased,
                dispersion = dispersion,
                df = c(object$rank, df.r, df.f),
                cov.unscaled = covmat.unscaled,
                cov.scaled = covmat))
  
  if(correlation && p > 0) {
    dd = sqrt(diag(covmat.unscaled))
    ans$correlation =
    covmat.unscaled/outer(dd,dd)
    ans$symbolic.cor = symbolic.cor
  }
  class(ans) = "summary.glm"
  return(ans)
}

# ----------------------------------------------------------------------
# anova.glm 
# ----------------------------------------------------------------------

anova.glm <- function(object, ..., dispersion=NULL, test=NULL)
{
  ## check for multiple objects
  dotargs = list(...)
  named = if (is.null(names(dotargs)))
	rep(FALSE, length(dotargs)) else (names(dotargs) != "")
  if(any(named))
	warning("the following arguments to 'anova.glm' are invalid and dropped: ",
          paste(deparse(dotargs[named]), collapse=", "))
  dotargs = dotargs[!named]
  is.glm = unlist(lapply(dotargs,function(x) inherits(x,"glm")))
  dotargs = dotargs[is.glm]
  if (length(dotargs))
	return(anova.glmlist(c(list(object), dotargs),
                       dispersion = dispersion, test=test))
  
  ## extract variables from model
  
  varlist = attr(object$terms, "variables")
  
  ## must avoid partial matching here.
  x =
	if (n = match("x", names(object), 0L))
  object[[n]]
	else model.matrix(object)
  varseq = attr(x, "assign")
  nvars = max(0, varseq)
  resdev = resdf = NULL
  
  ## if there is more than one explanatory variable then
  ## recall glm.fit to fit variables sequentially
  
  if(nvars > 1) {
    method = object$method
    if(!is.function(method))
    method = get(method, mode = "function", envir=parent.frame())
    
    ## allow for 'y = FALSE' in the call (PR#13098)
    y = object$y
    if(is.null(y)) { ## code from residuals.glm
      mu.eta = object$family$mu.eta
      eta = object$linear.predictors
      y =   object$fitted.values + object$residuals * mu.eta(eta)
    }
    for(i in 1L:(nvars-1)) {
      ## explanatory variables up to i are kept in the model
      ## use method from glm to find residual deviance
      ## and df for each sequential fit
	    fit = method(x=x[, varseq <= i, drop = FALSE],
                    y=y,
                    weights=object$prior.weights,
                    start	 =object$start,
                    offset =object$offset,
                    family =object$family,
                    control=object$control)
	    resdev = c(resdev, fit$deviance)
	    resdf = c(resdf, fit$df.residual)
    }
  }
  
  ## add values from null and full model
  
  resdf = c(object$df.null, resdf, object$df.residual)
  resdev = c(object$null.deviance, resdev, object$deviance)
  
  ## construct table and title
  
  table = data.frame(c(NA, -diff(resdf)),
                      c(NA, pmax(0, -diff(resdev))), resdf, resdev)
  tl = attr(object$terms, "term.labels")
  if (length(tl) == 0L) table = table[1,,drop=FALSE] # kludge for null model
  dimnames(table) = list(c("NULL", tl),
                          c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
  title = paste("Analysis of Deviance Table", "\n\nModel: ",
                 object$family$family, ", link: ", object$family$link,
                 "\n\nResponse: ", as.character(varlist[-1L])[1L],
                 "\n\nTerms added sequentially (first to last)\n\n", sep="")
  
  ## calculate test statistics if needed
  
  df.dispersion = Inf
  if(is.null(dispersion)) {
    dispersion = summary(object, dispersion=dispersion)$dispersion
    df.dispersion = if (dispersion == 1) Inf else object$df.residual
  }
  if(!is.null(test)) {
    if(test == "F" && df.dispersion == Inf) {
      fam = object$family$family
      if(fam == "binomial" || fam == "poisson")
      warning(gettextf("using F test with a %s family is inappropriate",
                       fam),
              domain = NA)
      else
      warning("using F test with a fixed dispersion is inappropriate")
    }
    table = stat.anova(table=table, test=test, scale=dispersion,
                        df.scale=df.dispersion, n=NROW(x))
  }
  structure(table, heading = title, class= c("anova", "data.frame"))
}


