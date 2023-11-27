##SIGNALE FILTERED EXTRACTED
##FROM ASTSA PACKAGE

function (series, L = c(3, 3), M = 50, max.freq = 0.05, min.freq) 
{
  ts = stats::ts
  tsp = stats::tsp
  par = graphics::par
  plot = graphics::plot
  dev.new = grDevices::dev.new
  abline = graphics::abline
  ccf = stats::ccf
  ts.intersect = stats::ts.intersect
  ts.plot = stats::ts.plot
  na.omit = stats::na.omit
  lines = graphics::lines
  L = 2 * floor((L - 1)/2) + 1
  if (max.freq < 0 || max.freq > 0.5) 
    stop("max.freq must be between 0 and 1/2")
  chek = FALSE
  if (max.freq < (1/M)) {
    M = 2.5 * (1/max.freq)
    chek = TRUE
  }
  M = 2 * floor(M/2)
  if (chek == TRUE) {
    cat("WARNING: must have max.freq > 1/M -- M changed to", 
        M, "\n")
  }
  tspar = tsp(series)
  series = ts(series, frequency = 1)
  spectra = stats::spec.pgram(series, spans = L, plot = FALSE)
  A <- function(nu) {
    qwe = ifelse((nu > min.freq && nu < max.freq), 1, 0)
    qwe
  }
  N = 2 * length(spectra$freq)
  sampled.indices = (N/M) * (1:(M/2))
  fr.N = spectra$freq
  fr.M = fr.N[sampled.indices]
  spec.N = spectra$spec
  spec.M = spec.N[sampled.indices]
  A.desired = vector(length = length(fr.M))
  for (k in 1:length(fr.M)) A.desired[k] = A(fr.M[k])
  delta = 1/M
  Omega = seq(from = 1/M, to = 0.5, length = M/2)
  aa = function(s) 2 * delta * sum(exp((0+2i) * pi * Omega * 
                                         s) * A.desired)
  S = ((-M/2 + 1):(M/2 - 1))
  a = vector(length = length(S))
  for (k in 1:length(S)) a[k] = aa(S[k])
  a = Re(a)
  h = 0.5 * (1 + cos(2 * pi * S/length(S)))
  a = a * h
  A.M = function(nu) Re(sum(exp(-(0+2i) * pi * nu * S) * a))
  A.attained = vector(length = length(fr.N))
  A.theoretical = vector(length = length(fr.N))
  for (k in 1:length(fr.N)) {
    A.attained[k] = A.M(fr.N[k])
    A.theoretical[k] = A(fr.N[k])
  }
  series.filt = stats::filter(series, a, sides = 2)
  series = ts(series, start = tspar[1], frequency = tspar[3])
  series.filt = ts(series.filt, start = tspar[1], frequency = tspar[3])
  return(invisible(series.filt))
}