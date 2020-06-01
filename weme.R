library(ggplot2)
library(RColorBrewer)
library(doMC)
library(reshape2)

P=100
F=100

list_to_vec <- function(tc){
    out = mat.or.vec(length(tc),1)
    for (i in 1:length(tc)){
        out[i] = tc[[i]][1]
    }
    return(out)
}

get_data <- function(fns){
    ncs = list()
    alldf = data.frame()
    tm=mat.or.vec(P*2,0)
    ms = list()
    for (fn in fns){
        d = read.csv(fn,sep="\t")
        d=d[order(d$proportion,decreasing=TRUE),]
        if (length(d$n_ssms)==0){
            d$n_ssms = d$n_variants
        }
        d$r_ssms = d$n_ssms/sum(d$n_ssms)
        outv = mat.or.vec(P*2,2)
        for(i in 1:P){
            ind = sum(cumsum(d$r_ssms)<i/P)+1
            if (ind > nrow(d)){ind=nrow(d)}
            outv[i*2,1] = d$proportion[ind]
            outv[i*2-1,1] = d$proportion[ind]
            outv[i*2,2] = i
            outv[i*2-1,2] = i-1
        }
        m = strsplit(fn,"/")[[1]][1]
        outdf = data.frame(pssm=outv[,2],phi=outv[,1],method=m)
        ms = c(ms,m)
        alldf = rbind(alldf,outdf)
        ncs = c(ncs,nrow(d))
        tm=cbind(tm,outv[,1])
    }
    mphi = apply(tm,1,median)
    if (length(ncs) >= 4){
        if (length(ncs)>=8){
            TOR = 3
        }
        else if (length(ncs) >= 5){
            TOR = 2
        }
        else{
            TOR=1
        }
        for (i in 1:TOR){
            fits = mat.or.vec(length(ncs),1)
            for (j in 1:length(ncs)){
                fits[j] = compute_fit(tm[,j],mphi)
            }
            to_r = which(fits==min(fits))[1]
            tm = tm[,-to_r]
            alldf = alldf[alldf$method != ms[[to_r]],]
            ncs = ncs[-to_r]
            mphi = apply(tm,1,median)
            ms = ms[-to_r]
        }
    }
    outdf = data.frame(pssm=outv[,2],phi=apply(tm,1,median),method="median")
    alldf = rbind(alldf,outdf)
    return(list(median(list_to_vec(ncs)),alldf))
}


get_data_perms <- function(t_fns){
    ncs = list()
    alldf = data.frame()
    tm=mat.or.vec(P*2,0)
    ms = list()
    for (fn in t_fns){
        d = read.csv(fn,sep="\t")
        d=d[order(d$proportion,decreasing=TRUE),]
        if (length(d$n_ssms)==0){
            d$n_ssms = d$n_variants
        }
        d$r_ssms = d$n_ssms/sum(d$n_ssms)
        outv = mat.or.vec(P*2,2)
        for(i in 1:P){
            ind = sum(cumsum(d$r_ssms)<i/P)+1
            if (ind > nrow(d)){ind=nrow(d)}
            outv[i*2,1] = d$proportion[ind]
            outv[i*2-1,1] = d$proportion[ind]
            outv[i*2,2] = i
            outv[i*2-1,2] = i-1
        }
        m = strsplit(fn,"/")[[1]][1]
        outdf = data.frame(pssm=outv[,2],phi=outv[,1],method=m)
        ms = c(ms,m)
        alldf = rbind(alldf,outdf)
        ncs = c(ncs,nrow(d))
        tm=cbind(tm,outv[,1])
    }
    mphi = apply(tm,1,median)
    if (length(ncs) >= 4){
        if (length(ncs)>=8){
            TOR = 3
        }
        else if (length(ncs) >= 5){
            TOR = 2
        }
        else{
            TOR=1
        }
        for (i in 1:TOR){
            fits = mat.or.vec(length(ncs),1)
            for (j in 1:length(ncs)){
                fits[j] = compute_fit(tm[,j],mphi)
            }
            to_r = which(fits==min(fits))[1]
            tm = tm[,-to_r]
            alldf = alldf[alldf$method != ms[[to_r]],]
            ncs = ncs[-to_r]
            mphi = apply(tm,1,median)
            ms = ms[-to_r]
        }
    }
    outdf = data.frame(pssm=outv[,2],phi=apply(tm,1,median),method="median")
    alldf = rbind(alldf,outdf)
    return(list(median(list_to_vec(ncs)),alldf))
}

pssms_to_vec <- function(pssms,med){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=mean(med$phi)
        }
        else if (i==1){
            d$proportion[i] = mean(med$phi[1:round(P*2*pssms[i])])
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = mean(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,1)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        if (ind > nrow(d)){ind=nrow(d)}
        outv[i*2] = d$proportion[ind]
        outv[i*2-1] = d$proportion[ind]
    }
    return(outv)
}

pssms_to_df <- function(pssms,med){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=mean(med$phi)
        }
        else if (i==1){
            d$proportion[i] = mean(med$phi[1:round(P*2*pssms[i])])
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = mean(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,2)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        outv[i*2,1] = d$proportion[ind]
        outv[i*2-1,1] = d$proportion[ind]
        outv[i*2,2] = i
        outv[i*2-1,2] = i-1
    }
    outdf = data.frame(pssm=outv[,2],phi=outv[,1],method="consensus")
}

pssms_to_vec_median_purity <- function(pssms,med,purity){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=purity
        }
        else if (i==1){
            d$proportion[i] = purity
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,1)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        if (ind > nrow(d)){ind=nrow(d)}
        outv[i*2] = d$proportion[ind]
        outv[i*2-1] = d$proportion[ind]
    }
    return(outv)
}

pssms_to_vec_median <- function(pssms,med){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=median(med$phi)
        }
        else if (i==1){
            d$proportion[i] = median(med$phi[1:round(P*2*pssms[i])])
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,1)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        if (ind > nrow(d)){ind=nrow(d)}
        outv[i*2] = d$proportion[ind]
        outv[i*2-1] = d$proportion[ind]
    }
    return(outv)
}

pssms_to_df_median <- function(pssms,med){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=median(med$phi)
        }
        else if (i==1){
            d$proportion[i] = median(med$phi[1:round(P*2*pssms[i])])
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,2)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        outv[i*2,1] = d$proportion[ind]
        outv[i*2-1,1] = d$proportion[ind]
        outv[i*2,2] = i
        outv[i*2-1,2] = i-1
    }
    outdf = data.frame(pssm=outv[,2],phi=outv[,1],method="consensus")
}

pssms_to_df_median_purity <- function(pssms,med,purity){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=purity
        }
        else if (i==1){
            d$proportion[i] = purity
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    outv = mat.or.vec(P*2,2)
    for(i in 1:P){
        ind = sum(cumsum(pssms)<i/P)+1
        outv[i*2,1] = d$proportion[ind]
        outv[i*2-1,1] = d$proportion[ind]
        outv[i*2,2] = i
        outv[i*2-1,2] = i-1
    }
    outdf = data.frame(pssm=outv[,2],phi=outv[,1],method="consensus")
}

compute_fit <- function(v1,v2){
    return(1-sum(abs(v1-v2))/length(v1))
}

nm_search <- function(nbps,mdf){
    best_pssms = rep(1/nbps,nbps)
    if (nbps==1){
        return(best_pssms)
    }
    score <- function(pssms){
        cv = pssms_to_vec_median(pssms,mdf)
        return(compute_fit(cv,mdf$phi))
    }
    solnp(best_pssms, #starting values (random - obviously need to be positive and sum to 15)
      score, #function to optimise
      eqfun=sum, #equality function
      eqB=1.0,   #the equality constraint
      LB=rep(0,nbps), #lower bound for parameters i.e. greater than zero
      UB=rep(1,nbps))
}


grid_search <- function(nbps,mdf){
    nattempt=0
    mv = mdf$phi
    best_pssms = mat.or.vec(nbps,1)
    best_score=-1
    for (i in 1:(F^(nbps-1)))
    {
        pssms = gen_pssms(nbps,i)
        if (sum(pssms <= 0) > 0){
            next
        }
        nattempt = nattempt+1
        cv = pssms_to_vec_median(pssms,mdf)
        score = compute_fit(cv,mv)
        if (score > best_score){
            best_score = score
            best_pssms = pssms
        }
    }
    print(best_pssms)
    if (nbps==1){
        best_pssms[1]=1
    }
    return(best_pssms)
}

grid_search_purity <- function(nbps,mdf,purity){
    nattempt=0
    mv = mdf$phi
    best_pssms = mat.or.vec(nbps,1)
    best_score=-1
    for (i in 1:(F^(nbps-1)))
    {
        pssms = gen_pssms(nbps,i)
        if (sum(pssms <= 0) > 0){
            next
        }
        nattempt = nattempt+1
        cv = pssms_to_vec_median_purity(pssms,mdf,purity)
        score = compute_fit(cv,mv)
        if (score > best_score){
            best_score = score
            best_pssms = pssms
        }
    }
    print(best_pssms)
    if (nbps==1){
        best_pssms[1]=1
    }
    return(best_pssms)
}

gen_pssms <- function(nbps,n){
    fac=1
    pssms = mat.or.vec(nbps,1)
    for (i in 1:nbps){
        if (i==1){
            pssms[i] = (n %% F)/F
        }
        else if(i==nbps){
            pssms[i] = 1-cumsum(pssms)[i]
        }
        else{
            pssms[i] = (floor(n/fac) %% F)/F
        }
        fac = fac*F
    }
    return(pssms)
}

make_plot <- function(filename, result){
    png(file=filename,height=800,width=1000)
    p=ggplot(result,
             aes(x=pssm,y=phi,color=method,alpha=method,linetype=method)) +
        geom_line(size=1) +
        scale_alpha_manual(values=c(rep(0.45,length(levels(factor(result$method)))-2),1,1)) +
        scale_linetype_manual(values=c(rep("solid",length(levels(factor(result$method)))-2),"11","33")) +
        scale_color_manual(values=c(brewer.pal(length(levels(factor(result$method)))-2,"Set2")[1:(length(levels(factor(result$method)))-2)],
                                    "red","black")) +
        coord_flip()
    print(p)
    dev.off()
}

make_subclonal_structure <- function(sid, pssms, med){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=median(med$phi)
        }
        else if (i==1){
            d$proportion[i] = median(med$phi[1:round(P*2*pssms[i])])
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    sc = data.frame(cluster=1:length(pssms),n_ssms=pssms*100,proportion=d$proportion)
    write.table(file=paste(sid,"_subclonal_structure.txt",sep=""),sc,sep="\t",quote=FALSE,row.names=FALSE)
}

make_subclonal_structure_purity <- function(sid, pssms, med, purity){
    d = data.frame(r_ssms=pssms)
    d$proportion = rep(0,length(pssms))
    for (i in 1:length(pssms)){
        if (length(pssms)==1){
            d$proportion[i]=purity
        }
        else if (i==1){
            d$proportion[i] = purity
        }
        else{
            p1 = cumsum(pssms)[i-1]
            p2 = cumsum(pssms)[i]
            d$proportion[i] = median(med$phi[round(P*2*p1):round(P*2*p2)])
        }
    }
    sc = data.frame(cluster=1:length(pssms),n_ssms=pssms*100,proportion=d$proportion)
    write.table(file=paste(sid,"_subclonal_structure.txt",sep=""),sc,sep="\t",quote=FALSE,row.names=FALSE)
}

genconsensus <- function(sids, fns, ncores=NULL, rounddown=TRUE, purities=NULL)
{
    if (is.null(ncores)){
        registerDoMC(1)
    }
    else{
        registerDoMC(ncores)
    }
    if (is.null(purities)){
        res = foreach (i=1:length(sids)) %dopar% {
            sid = sids[i]
            r = get_data(fns[i])
            if (rounddown){
                best_pssms = grid_search(floor(r[[1]]),r[[2]][r[[2]]$method=="median",])
            }
            else{
                best_pssms = grid_search(ceiling(r[[1]]),r[[2]][r[[2]]$method=="median",])
                make_subclonal_structure(sid,best_pssms,r[[2]][r[[2]]$method=="median",])
                condf = pssms_to_df_median(best_pssms,r[[2]][r[[2]]$method=="median",])
                tpdf = rbind(r[[2]],condf)
                make_plot(paste(sids[i],".png",sep=""),tpdf)
            }
        }
    }
    else{
        res = foreach (i=1:length(sids)) %dopar% {
            sid = sids[i]
            purity = purities[i]
            r = get_data(fns[i])
            if (rounddown){
                best_pssms = grid_search_purity(floor(r[[1]]),r[[2]][r[[2]]$method=="median",],purity)
            }
            else{
                best_pssms = grid_search_purity(ceiling(r[[1]]),r[[2]][r[[2]]$method=="median",],purity)
            }
            make_subclonal_structure_purity(sid,best_pssms,r[[2]][r[[2]]$method=="median",],purity)
            condf = pssms_to_df_median_purity(best_pssms,r[[2]][r[[2]]$method=="median",],purity)
            tpdf = rbind(r[[2]],condf)
            make_plot(paste(sids[i],".png",sep=""),tpdf)
        }
    }
}


