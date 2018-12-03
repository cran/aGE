#' aGE joint test
#'
#' @param Y a numeric vector of phenotype values
#' @param G a matrix or data frame for all RVs in the test gene or genomic region. The order of rows must match the order of Y. Missing is imputed as 0.
#' @param cov a matrix or data frane with first column as the environmental variable to be tested. The order of rows must match the order of Y.
#' @param model "binomial" for binary traits or "gaussian" for quantitative traits.
#' @param pow Gamma set used to build a family of tests, default=c(1:6) for rare variants
#' @param n.perm number of simulation to calculate the p-values, default=1000. Can increase to higher value depending on the signficiance level.
#' @param method 'Simulation': Monte Carlo Method
#' @param nonparaE  "T": use cubic splines for the environmental variable to fit the model; "F": use a linear function of the environmental variable to fit the model
#' @param DF degree of freedom to use in the cubic splines, default=10. This option only works when nonparaE is set to "T".
#'
#' @return p-values
#' @export
#' @examples {
#'  set.seed(12345)
#'  phenotype <- c(rep(1,50),rep(0,50))
#'  genotype <- data.frame(g1=sample(c(rep(1,10),rep(0,90))),g2=sample(c(rep(1,5), rep(0,95))))
#'  covariates <- data.frame(Envir=rnorm(100), Age=rnorm(100,60,5))
#'  exD <- list(Y=phenotype, G=genotype, X=covariates)
#'  aGE.joint(Y=exD$Y, G=exD$G, cov=exD$X, model='binomial') 
#'  }
aGE.joint<-function(Y, G, cov = NULL, model = c("gaussian", "binomial"), pow = c(1:6), n.perm = 1000, method=c('Simulation'), nonparaE=F, DF=10){
  model = match.arg(model)
  G[is.na(G)] <- 0
  cov<-as.matrix(cov)
  if (!is.null(cov)) cov <- scale(cov,center=T, scale=T)
  if (!is.null(G))  G <- scale(G,center=T, scale=F)
  tdat1 <- data.frame(trait=Y,cov,G)
  Ind <- stats::complete.cases(tdat1)
  tdat2 <- tdat1[Ind,]
  n <- nrow(tdat2)
  nc <- ncol(cov)
  nZ <- ncol(G)
  cov1 <- as.matrix(cov[Ind,])
  G <- as.matrix(G[Ind,])
  E <- cov1[,1]
  Y <- as.matrix(Y[Ind])
  XdE <- cov1[,-1,drop=F]
  if (nonparaE) {ME <- spl.X(x=E,df=DF); cov1=cbind(cov1,ME)}


  if (is.null(cov1)) {
   res=tdat2$trait-mean(tdat2$trait)
   U_main <- as.vector(t(G)%*%t(t(res)))
   EG<-E*G
   U_int <- as.vector(t(EG) %*% t(t(res)))
  } else {
       tdat3 <- data.frame(trait=Y, cov1)
	   fit1 <- stats::glm(trait~.,family=model,data=tdat3)
	   pis <- stats::fitted.values(fit1)
	   res <- tdat3$trait - pis
	   U_main <- as.vector(t(G)%*%t(t(res)))
	   EG <- E*G
	   U_int <- as.vector(t(EG) %*% t(t(res)))
	   }

  #calculate the observed test statistics
  Ts_main = rep(NA, length(pow))
  Ts_int = rep(NA, length(pow))
  Ts=rep(NA,length(pow))
  npow = pow
  for (j in 1:length(pow)) {
        if (pow[j] < Inf) {
            Ts_main[j] = sum(U_main^pow[j])
			Ts_int[j] = sum(U_int^pow[j])
        }
        else {
             Ts_main[j] = max(abs(U_main))
			 Ts_int[j] =max(abs(U_int))
             npow[j] = 0
        }
    }
	Ts=Ts_main+Ts_int

  if (method=='Simulation') { V11<-simulation_var(n=n,cov=cov, E=E,G=G,EG=EG, res=res, model=model, pis=pis)  }
  s <- sample(1:10^5,1)
  set.seed(s)
  result<-calculation(method=method, n.perm=n.perm, res=res, XdE=XdE, nZ=nZ, n=n,  V11=V11, pow=pow, Ts_main=Ts_main,Ts_int=Ts_int,Ts=Ts, G=G, EG=EG, pis=pis, cov1=cov1, model=model )
 	return(result)
	}

