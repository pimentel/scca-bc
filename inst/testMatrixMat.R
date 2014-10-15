a <- c(3, 8, 4)
b <- c(2, 7, 10)
X <- matrix(rnorm(10*15, mean = 5, sd = 2), nrow = 10)
Y <- matrix(rnorm(10*15), nrow = 10)
Dx <- diag(runif(15))

standardMethod <- (t(a) %*% t(X%*% Dx)) %*% (Y %*% Dx) %*% b


sccaProd <- a %*% t(X) %*% Y %*% b
diag(Dx)

# (t(a) %*% t(X%*% Dx)) %*% (Y %*% Dx) %*% b
# from Haiyan's notes
(a %*% X %*% Dx) %*% t( (b %*% Y %*% Dx) )

Xtilde <- t(X %*% Dx)
Ytilde <- t(Y %*% Dx)

(a %*% t(Xtilde)) %*% (Ytilde %*% b)

# incorrect dim
t((Dx %*% X %*% a)) %*% (Dx %*% Y %*% b)


t(a %*% X %*% Dx) %*% (b %*% Y %*% Dx) 

t(X %*% a) %*% (Y %*% b) 

a %*% t(X) %*% (Y %*% b) 

(X %*% a) %*% t( (Y %*% b) )



t(diag(Dx)) %*% ( t(a) %*% X %*% t(t(b) %*% Y) ) %*% diag(Dx)
diag(Dx) %*% ( t(a) %*% X %*% t(t(b) %*% Y) ) %*% diag(Dx)


diag(Dx) %*% diag(as.numeric((a %*% X) * (b %*% Y))) %*% diag(Dx)


Q <- - 2 * diag(as.numeric((a %*% X) * (b %*% Y)))
d <- diag(Dx)

library(kernlab)
ipop(c = rep(0, length(d)), 
     H = Q, 
     A = diag(length(d)), 
     b = rep(0, length(d)),
     l = rep(0, length(d)),
     u = rep(1, length(d)),
     r = rep(1, length(d)),
             )

