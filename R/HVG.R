# find highly variable genes in each data set separately a la Brennecke et al

# define a convenience function for HVG definition
find_hvg <- function(dataframe, plot=FALSE, p.threshold=1e-2, return.ranks=FALSE,
                     return.p=FALSE){
  # define a set of highly variable gene for the GFP+ and GFP- separately

  # assume gene names are in the rows
  # even if they aren't this will still get
  # the input row ordering
  gene.names <- rownames(dataframe)
  means <- rowMeans(dataframe, na.rm = T)
  vars <- apply(dataframe, 1, var, na.rm=T)
  cv2 <- vars/(means^2)

  minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))

  # select genes with mean value greater than min value for fitting
  # remove values with 1/means == infinite
  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0

  useForFit <- recip.means <= 0.1

  # fit with a gamma-distributed GLM
  fit <- glmgam.fit(cbind(a0 = 1, a1tilde=recip.means[!useForFit]), cv2[!useForFit])

  # calculate % variance explained by the model fit
  resid.var <- var(fitted.values(fit) - cv2[!useForFit])
  total.var <- var(cv2[!useForFit])

  # get fitted values and mean-dispersion dependence line
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])

  xg <- seq(0, max(means[means != Inf]), length.out=100000)
  vfit <- (a1/xg) + a0

  # add confidence intervals
  d.f <- ncol(dataframe) - 1

  # rank genes by the significance of their deviation from the fit
  # to call HVGs
  a.fit <- (a1/means) + a0
  varFitRatio <- vars/(a.fit * means^2)
  varOrder <- order(varFitRatio, decreasing=T)

  oed <- dataframe[varOrder, ]

  if(plot == TRUE){
    smoothScatter(means, cv2, xlab="Mean expression", ylab=expression("CV"^2))
    lines(xg, vfit, col="black", lwd=3 )
    lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col="black")
    lines(xg, vfit * qchisq(0.025, d.f)/d.f,lty=2,col="black")
    # display the 100 most highly variable genes
    points(means[varOrder[1:100]], cv2[varOrder[1:100]], col='red')
  }

  pvals <- pchisq(varFitRatio * d.f, d.f, lower.tail = F)
  pvals[is.na(pvals)] <- 1.0
  adj.pvals <- p.adjust(pvals, method='fdr')
  HVG <- adj.pvals <= p.threshold

  if(return.ranks){
    # order p-values, then subset past a threshold
    rank.p <- adj.pvals[order(adj.pvals, decreasing=FALSE)]
    order.names <- gene.names[order(adj.pvals, decreasing=FALSE)]
    thr.p <- rank.p <= p.threshold
    HVG <- order.names[thr.p]
  }

  if(return.p){
    HVG <- cbind(HVG, adj.pvals)
  }

  return(HVG)
}
