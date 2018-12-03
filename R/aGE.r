#' aGE interaction test
#'
#' @param Y a numeric vector of phenotype values
#' @param G a matrix for all RVs in the test gene or genomic region. The order of rows must match the order of Y. Missing is imputed as 0.
#' @param cov a matrix with first column as the environmental variable to be tested. The order of rows must match the order of Y.
#' @param model "binomial" for binary traits or "gaussian" for quantitative traits.
#' @param pow Gamma set used to build a family of tests, default=c(1:6) for rare variants
#' @param n.perm number of simulation to calculate the p-values, default=1000. Can increase to higher value depending on the signficiance level.
#' @param method only have one option: "Simulation", also called Monte Carlo Method.
#' @param nonparaE "T": use cubic splines for the environmental variable to fit the model; "F": use a linear function of the environmental variable to fit the model
#' @param DF degree of freedom to use in the cubic splines, default=10. This option only works when nonparaE is set to "T"
#' @param stepwise an option to speed up the simulation procedure for large n.perm number in real-data application. Up to $n.perm=10^8$
#'
#' @return p-values
#' @export
#'
#' @examples {
#'  set.seed(12345)
#'  phenotype <- c(rep(1,50),rep(0,50))
#'  genotype <- data.frame(g1=sample(c(rep(1,10),rep(0,90))),g2=sample(c(rep(1,5), rep(0,95))))
#'  covariates <- data.frame(Envir=rnorm(100), Age=rnorm(100,60,5))
#'  exD <- list(Y=phenotype, G=genotype, X=covariates)
#'  aGE(Y=exD$Y, G=exD$G, cov=exD$X, model='binomial', nonparaE=FALSE, stepwise=FALSE)
#'  }
aGE <- function(Y, G, cov = NULL, model = c("gaussian", "binomial"), pow = c(1:6), n.perm = 1000,method='Simulation', nonparaE=F ,  DF=10, stepwise=T){

    model = match.arg(model, c('gaussian','binomial'))
  method = match.arg(method, c('Simulation'))
  G[is.na(G)] <- 0
  n <- length(Y)
  Ind <- stats::complete.cases(Y) & stats::complete.cases(cov) & stats::complete.cases(G)
  Y <- Y[Ind]
  cov <- cov[Ind,]
  G<-as.matrix(G[Ind,])
  cov <- scale(cov,center=T, scale=T)
  G <- scale(G,center=T, scale=F)
  n <- sum(Ind)
  E <- ME <- cov[,1]
  XUs <- E*G

  if (nonparaE) {ME <- spl.X(x=E,df=DF); cov <- cbind(ME, cov[,-1])}
  tdat1 <- data.frame(cbind(Y=Y,cov))
  pis<-res<-tau<-phi<-W<-NA

  nc <- ncol(cov)
  tdat1 <- data.frame(trait=Y,cov,G,id = rep(1, n))
  nullmodel <- GLMMaSPU(tdat1=tdat1,model=model,G=G, cov=cov, nc=nc, n=n)
  pis <- nullmodel$pis; tau<-nullmodel$tau; phi<-nullmodel$phi; W<-nullmodel$W
  res <- as.vector(tdat1$trait - pis)

  if (tau==0) return(pvs=rep(NA,(length(pow)+4))) else {
  Gw<-W%*%G
  D <- solve( diag(1/tau, ncol(G)) + t(G) %*%  Gw )
  Vinv <- W - Gw %*% D %*% t(Gw)

  if (all(is.na(res))) { return(pvs=rep(NA,(length(pow)+4))) } else if (method=='Simulation') {  
    s <- sample(1:10^5,1)
    set.seed(s)
    if (!stepwise) {
  result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=n.perm, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
 } else {
  if ( n.perm <= 1000 ) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1000, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
  if (min(result$pvs) < 5/1e3 &  n.perm <= 1e4 ) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1e4, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
  if (min(result$pvs) < 5/1e4 &  n.perm <= 1e5) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1e5, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
  if (min(result$pvs) < 5/1e5 &  n.perm <= 1e6) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1e6, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
  if (min(result$pvs) < 5/1e6 &  n.perm <= 1e7) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1e7, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
  if (min(result$pvs) < 5/1e7 &  n.perm <= 1e8) result <- RobustSim(n=n, cov=cov, G=G, XUs=XUs, res=res,  n.perm=1e8, pow=pow, Vinv=Vinv,  model=model, Y=Y,   phi=phi)
 }
  pvs <- result$pvs
  rm(Vinv)
  rm(tdat1)
  if (result$nperm < n.perm*0.9)  { warning("more than 10% of error in simulation! ") }
  return(pvs)
  }
  }
  }
