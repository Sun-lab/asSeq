fishBB = function(n, y, rho, pi1){
n
y
rho
pi1

  nind = length(n)
  pi2 = 1-pi1
  y2 = n-y
  out2 = out = matrix(0, nrow=2, ncol=2)
  x1a = x2a = x3a = 0    

  for(j in 1:nind){
  #j = 1
    if(n[j]>0){
    for(r in 1:n[j]){
      div1 = ((1-rho)*pi1 + (r-1)*rho)^2
      div2 = ((1-rho)*pi2 + (r-1)*rho)^2
      div3 = ((1-rho)     + (r-1)*rho)^2
      
      x1a = x1a + pbetabinom(q=n[j]-r, size=n[j], prob=pi2, rho=rho)/div1 +
                  pbetabinom(q=n[j]-r, size=n[j], prob=pi1, rho=rho)/div2

      x2a = x2a + pbetabinom(q=n[j]-r, size=n[j], prob=pi2, rho=rho)/div1*pi1 -
                  pbetabinom(q=n[j]-r, size=n[j], prob=pi1, rho=rho)/div2*pi2
      x3a = x3a + pbetabinom(q=n[j]-r, size=n[j], prob=pi2, rho=rho)/div1*pi1^2 +
                  pbetabinom(q=n[j]-r, size=n[j], prob=pi1, rho=rho)/div2*pi2^2 -
                  1/div3
    }    
  }
  }
  out2[2,2] = x1a*(1-rho)^2
  out2[1,2] = out2[2,1] = x2a*(rho-1)/rho
  out2[1,1] = x3a/rho^2
  out2
}



nllBB = function(pars, n, y){
  -llBB(rho=pars[1], pi1=pars[2], n=n, y=y)
}
llBB = function(n, y, rho, pi1){
  nind = length(n)
  pi2 = 1-pi1
  y2 = n-y
  out = 0
  x1 = 0

  for(j in 1:nind){
  #j = 1
    if(y[j]>0){  
    for(r in 1:y[j]){
      out = out + log((1-rho)*pi1 + (r-1)*rho)
    }
    }
    if(y2[j]>0){
    for(r in 1:y2[j]){
      out = out + log((1-rho)*pi2 + (r-1)*rho)
    }
    }
    if(n[j]>0){
    for(r in 1:n[j]){
      out = out - log((1-rho) + (r-1)*rho)
    }
    out = out + lchoose(n[j],y[j])
    }
  }    
  out
}





