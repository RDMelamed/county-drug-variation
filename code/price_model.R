library(MCMCglmm)

getBrandPreference <- function(cxd, druginfo){
    countylist = rownames(cxd)
    rownames(cxd) = paste("X", countylist, sep="")
    randomclass = colnames(druginfo)
    randomclass = randomclass[randomclass != "brandred"]
    cxd = cbind(druginfo, t(cxd))
    
    coefdat = data.frame(matrix(0,nrow=length(countylist),ncol=length(randomclass)), row.names=countylist)
    for (cty in countylist){
        coefs = c()    
        for (r in randomclass){
            
            mc = MCMCglmm(as.formula(paste("X", cty, " ~ 1 + brandred",sep="")),
                          random=as.formula(paste("~", r)), data=cxd,verbose=F)
            coefs = c(coefs, mean(mc$Sol[,'brandred']))
        }
        coefdat[cty,] = coefs
    }
    coefdat
}

projVbrand <- function(dat){
    res = data.frame(matrix(0,8,3),row.names=colnames(dat)[1:8])
    for (v in colnames(dat)[1:8]){
        mc = MCMCglmm(as.formula(paste(v,' ~ 1 + brandred')), random=~classf, data=dat,verbose=F)
        res[v,] = quantile(mc$Sol[,'brandred'],probs=c(.05,.5,.95)) #mean(mc$Sol[,'brandred'])    
    }
    colnames(res)=c(5,50,95)
    res
}    

run <- function(){
    
    d2c = read.table('brand_drugclass.txt',sep='\t',header=T,row.names=1)
    cxd = read.table('county_deviance_norm.txt',sep='\t',row.names=1,colClasses=list('place'='character'),header=T)
    sxd = read.table('states_deviance_norm.txt',sep='\t',row.names=1,colClasses=list('place'='character'),header=T)

    cat('Getting state coefficients... takes less than 15 minutes\n')
    res = getBrandPreference(sxd, d2c)
    write.table(res, 'state_coef.txt',sep='\t',quote=F)

    cat('Getting county coefficients... takes more than one hour \n'    )
    res = getBrandPreference(cxd, d2c)
    write.table(res, 'county_coef.txt',sep='\t',quote=F)
    
    drugproj = read.table('drug_projection.txt',sep='\t',header=T,row.names=1)
    drugproj$classf = as.factor(d2c$X0)
    drugproj$brandred = d2c$brandred
    res = projVbrand(drugproj)
    write.table(res, 'component_price_effect_estimates.txt',sep='\t',quote=F)    
}
