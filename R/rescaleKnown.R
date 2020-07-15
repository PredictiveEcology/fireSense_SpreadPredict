rescaleKnown <- function(x, minNew, maxNew, minOrig, maxOrig) {
  a1 <- x - minOrig # brings min to zero
  a2 <- a1 * maxNew/max(a1)
  a2
}