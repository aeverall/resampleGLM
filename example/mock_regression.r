library(ggplot2)
source("../R/glmFit.R")

# Number of samples
N = 1000

# Coefficients of visible parameters
gamma = c(0.5, 0.3,0.2)
# Coefficient of hidden parameter
Qgamma = 3
# Alpha of negative binomial model
alpha = 2

# Covariates
Z = cbind(1, rnorm(N,0,1), rnorm(N,0,1))
# Hidden covariate
Q = rbinom(N, 1, 0.002)

# Draw samples
mu = exp(Z%*%gamma + Q*Qgamma)
Y = rnbinom(N, size=alpha, mu=mu)

pvalue_wald <- c()
pvalue_lr <- c()
pvalue_crt <- c()
pvalue_bts <- c()
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for (i in c(1:100)) {

    # Target drawn such that Y is independent of X
    pi = 0.1
    X = rbinom(N, 1, pi)

    # Generate dataframe
    data.df <- data.frame(tgt=X, Z)

    # Run distilled-CRT resampling
    covariates <- list(count=c("X1", "X2", "X3"))
    results_CRT <- resampleGLM(data.df, Y, 'tgt',
                                family='poisson', target_family='binomial', maxit=100, zeroinfl=F,
                                nulliter = 100, covariates=covariates,
                                resampling='dCRT')
    pvalue_wald <- c(pvalue_wald, pchisq( (results_CRT$betas['tgt']/sqrt(diag(results_CRT$covar)['tgt']))^2, df=1, lower.tail = FALSE))
    pvalue_lr <- c(pvalue_lr, pchisq(-2*(results_CRT$null_loglike - results_CRT$loglike), df=1, lower.tail = FALSE))
    pvalue_crt <- c(pvalue_crt, results_CRT$pv_score)

    # Run bootstrap resampling
    results_BTS <- resampleGLM(data.df, Y, 'tgt',
                            family='poisson', target_family='binomial', maxit=100, zeroinfl=F,
                            nulliter = 100, covariates=covariates,
                                   resampling='bootstrap')

    pvalue_bts <- c(pvalue_bts, results_BTS$pv_score)

    setTxtProgressBar(pb, i)

}

#### QQ plot of results ####

n <- length(pvalue_bts)
k <- c(1:n)

data.df <- data.frame(x=sort(-log10(qbeta(0.5, k, n+1-k))),
                      lower=sort(-log10(qbeta(0.05, k, n+1-k))),
                      upper=sort(-log10(qbeta(0.95, k, n+1-k))),
                      Wald=sort(-log10(pvalue_wald)),
                      LikelihoodRatio=sort(-log10(pvalue_lr)),
                      dCRT=sort(-log10(pvalue_crt)),
                      Bootstrap=sort(-log10(pvalue_bts)))
data.df.long <- data.df %>% pivot_longer(!c(x,lower,upper), names_to="pvalue")
p <- ggplot(data.df, aes(x=x))+
    geom_line(aes(x=x, y=x)) +
    geom_ribbon(aes(x=x, ymin=lower, ymax=upper), alpha=0.2) +
    geom_point(data=data.df.long, aes(x=x, y=value, color=pvalue)) +
    theme_bw() +
    scale_x_continuous(limits=c(0,-log10(qbeta(0.5, k, n+1-k))[1]), expand=c(0,0.01)) +
    scale_y_continuous(limits=c(0,max(data.df.long$value)), expand=c(0,0)) +
    xlab("log10 pvalue expected") +
    ylab("log10 pvalue observed")
ggsave("qqplot.pdf", width = 20, height = 20, units = "cm", dpi=200)
