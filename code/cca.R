library(CCP)
library(CCA)

docc <- function(dev, dem,dir, suf){
    cc1 <- cc(dev, dem)


    write.table(cc1$scores$corr.Y.yscores,file=paste(dir,'corr.Y.yscores.',suf,'.txt',sep=''),sep='\t')
    write.table(cc1$scores$corr.X.xscores,file=paste(dir,'corr.X.xscores.',suf,'.txt',sep=''),sep='\t')
    write.table(cc1$scores$xscores,file=paste(dir,'xscores.',suf,'.txt',sep=''),sep='\t')
    write.table(cc1$scores$yscores,file=paste(dir, 'yscores.',suf,'.txt',sep=''),sep='\t')

    write.table(cc1$ycoef,file=paste(dir,'ycoef.',suf,'.txt',sep=''),sep='\t')

    s1 = diag(sqrt(diag(cov(dev))))
    x = s1 %*% cc1$xcoef
    rownames(x)= rownames(cc1$xcoef)
    write.table(x,file=paste(dir, 'xcoef_std.',suf,'.txt',sep=''),sep='\t')

    s1 = diag(sqrt(diag(cov(dem))))
    x = s1 %*% cc1$ycoef
    rownames(x)= rownames(cc1$ycoef)
    write.table(x,file=paste(dir, 'ycoef_std.',suf,'.txt',sep=''),sep='\t')

    p = length(dev)
    q = length(dem)
    n = nrow(dev)
    ### below, "function" with cc1, n (nsamp) , p (dim1 ), q (dim2)
    ev <- (1 - cc1$cor^2)
    k <- min(p, q)
    m <- n - 3/2 - (p + q)/2

    w <- rev(cumprod(rev(ev)))

                                        # initialize
    d1 <- d2 <- f <- vector("numeric", k)

    for (i in 1:k) {
        s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
        si <- 1/s
        d1[i] <- p * q
        d2[i] <- m * s - p * q/2 + 1
        r <- (1 - w[i]^si)/w[i]^si
        f[i] <- r * d2[i]/d1[i]
        p <- p - 1
        q <- q - 1
    }

    pv <- pf(f, d1, d2, lower.tail = FALSE)
    dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv)
    return(list(cc=cc1, wilks=dmat))
}

ccFromFiles <- function(){
## 28 x 2837
demp = read.csv('results/demographics_cca.txt',sep='\t',row.names=1)

## 195 x 2837
dev = read.csv('results/deviances_cca.txt',sep='\t',row.names=1)
r = docc(dev, demp, 'results/cca/','')

res1 = p.asym(r$cc$cor,nrow(dev),ncol(dev),ncol(demp),tstat = "Wilks")
write.table(data.frame(res1),'results/cca_wilks_table.txt',sep='\t')
}
