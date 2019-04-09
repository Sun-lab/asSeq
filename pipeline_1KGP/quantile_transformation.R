normscore = function(vec) {
    len  = length(na.omit(vec))+1
    rank = rank(na.omit(vec))
    ties = (rank - floor(rank)) > 0
    new.vec = vec[!is.na(vec)] 
    new.vec[!ties]=qnorm(rank[!ties]/len)
    new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
    vec[!is.na(vec)] = new.vec
    vec
}
