champ.ComBat <-
function (dat, batch, mod, numCovs = NULL, par.prior = TRUE,
prior.plots = FALSE)
{
	Beta.NA<-NULL
	int.eprior<-NULL
	
    mod = cbind(mod, batch)
    check = apply(mod, 2, function(x) all(x == 1))
    mod = as.matrix(mod[, !check])
    colnames(mod)[ncol(mod)] = "Batch"
    
    if(sum(check) > 0 & !is.null(numCovs))
    	numCovs = numCovs - 1
    design <- design.mat(mod, numCov = numCovs)
    batches <- list.batch(mod)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs = any(is.na(dat))
    
    if(NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"),
        sep = " ")
    }
    cat("Standardizing Data across genes\n")

    tryCatch(
    
    if (!NAs) {                
        B.hat <- solve(t(design) %*% design) %*% t(design) %*%
        t(as.matrix(dat))
    }
    else {
        B.hat = apply(dat, 1, Beta.NA, design)
    },
    error=function(e)
    {
        message("Combat failed...Your slides may be confounded with batch or with each other. Analysis will proceed without batch correction\n")
        
    })
    if(exists("B.hat"))
    {
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array,
        n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var,
        na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,
    n.array)))
    cat("Fitting L/S model and finding priors\n")
    batch.design <- design[, 1:n.batch]
    if (!NAs) {
        gamma.hat <- solve(t(batch.design) %*% batch.design) %*%
        t(batch.design) %*% t(as.matrix(s.data))
    }
    else {
        gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
    }
    delta.hat <- NULL
    for (i in batches) {
        delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var,
        na.rm = T))
    }
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)
    if (prior.plots & par.prior) {
        par(mfrow = c(2, 2))
        tmp <- density(gamma.hat[1, ])
        plot(tmp, type = "l", main = "Density Plot")
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
        qqnorm(gamma.hat[1, ])
        qqline(gamma.hat[1, ], col = 2)
        tmp <- density(delta.hat[1, ])
        invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
        tmp1 <- density(invgam)
        plot(tmp, typ = "l", main = "Density Plot", ylim = c(0,
        max(tmp$y, tmp1$y)))
        lines(tmp1, col = 2)
        qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles",
        ylab = "Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
        title("Q-Q Plot")
    }
    gamma.star <- delta.star <- NULL
    if (par.prior) {
        cat("Finding parametric adjustments\n")
        for (i in 1:n.batch) {
            temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i,
            ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
            b.prior[i])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    else {
        cat("Finding nonparametric adjustments\n")
        for (i in 1:n.batch) {
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
            gamma.hat[i, ], delta.hat[i, ])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,
        ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1,
        n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
    n.array)))) + stand.mean
    return(bayesdata)
    }
    else(return())

}
