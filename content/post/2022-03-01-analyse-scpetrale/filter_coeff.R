##ESTIMATED FILTER FUNCTION
##FROM ASTSA PACKAGE

function(series, L, M, max.freq, min.freq){
  L <- 2*floor((L-1)/2)+1
  if (max.freq < 0 || max.freq > 0.5) 
    stop("max.freq must be between 0 and 1/2")
  chek <- FALSE
  if (max.freq < (1/M)) {
    M <- 2.5 * (1/max.freq)
    chek <- TRUE
  }
  M <- 2*floor(M/2)
  if (chek == TRUE) {
    cat("WARNING: must have max.freq > 1/M -- M changed to", 
        M, "\n")
  }
  tspar <- tsp(series) #start end and freq 
  series <- ts(series, frequency=1)
  spectra <- spec.pgram(series, spans=L, plot=FALSE) #smooth periodogramme
  A <- function(nu) {
    qwe <- ifelse((nu > min.freq && nu < max.freq), 1, 0) #freq max à capturer
    qwe
  }
  N <- 2*length(spectra$freq)
  sampled.indices <- (N/M) * (1:(M/2))
  fr.N <- spectra$freq
  fr.M <- fr.N[sampled.indices]
  A.desired <- vector(length=length(fr.M))
  for (k in 1:length(fr.M)) {
    A.desired[k]  <- A(fr.M[k])
  }
  delta <- 1/M
  Omega <- seq(from = 1/M, to = 0.5, length = M/2)
  aa <- function(s) {
    2 * delta * sum(exp((0+2i) * pi * Omega * s) * A.desired)
  }
  S <- ((-M/2+1):(M/2-1))
  a <- vector(length=length(S))
  for (k in 1:length(S)) {
    a[k] <- aa(S[k])
  }
  a <- Re(a)
  h <- 0.5*(1 + cos(2*pi*S/length(S))) #cosinus taper
  a <- a*h 
  qwe <- cbind(S, a)
  colnames(qwe) <- c("s", "a(s)")
  return(qwe)
}