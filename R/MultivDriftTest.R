#'Multivariate genetic drift test for 2 populations
#'
#'This function estimates populations evolving through drift from an ancestral 
#'population, given an effective population size, number of generations separating
#'them and the ancestral G-matrix. It calculates the magnitude of morphological 
#'divergence expected and compare it to the observed magnitude of morphological 
#'change. 
#'
#'@param population1 data.frame with original measurements for the ancestral population
#'@param population2 data.frame with original measurements for the derived population
#'@param G ancestral G matrix 
#'@param Ne effective population size estimated for the populations
#'@param generations time in generations separating both populations
#'@param iterations number of simulations to perform
#'@note Each trait is estimated independently. 
#'@return list with the 95% IC magnitude of morphological change expected through 
#'drift and the range of the observed magnitude of morphological change
#'@export
#'@examples
#'\donttest{
#'data(dentus)
#'A <- dentus[dentus$species== "A",-5]
#'B <- dentus[dentus$species== "B",-5]
#'G <- cov(A)
#'MultivDriftTest(A, B, G, Ne = 1000, generations = 250)
#'}
#'@references Hohenlohe, P.A ; Arnold, S.J. (2008). MIPod: a hypothesis testing framework for 
#'microevolutionary inference from patterns of divergence. 
#'American Naturalist, 171(3),
#' 366-385. doi: 10.1086/527498
#'@author Ana Paula Assis
MultivDriftTest <- function(population1, population2, G , Ne, generations, iterations = 1000){
  populations <- list(as.matrix(population1), as.matrix(population2))
  std.error <- function(x) sd(x)/sqrt(length(x))
  means <- lapply(populations, function(x) apply(x, 2, mean)) 
  se <- lapply(populations, function(x) apply(x, 2, std.error))
  D <- as.matrix((G * generations) / Ne)
  n.traits <- ncol(populations[[1]])
  population1.se <- matrix(NA, nrow = iterations, ncol = n.traits)
  population2.se <- matrix(NA, nrow = iterations, ncol = ncol(populations[[2]]))
  for ( i in seq_len(n.traits)) {
    population1.se[,i] <- runif(iterations, min= means[[1]][i]-se[[1]][i],
                                max=means[[1]][i]+se[[1]][i] ) 
    population2.se[,i] <- runif(iterations, min= means[[2]][i]-se[[2]][i],
                                max=means[[2]][i]+se[[2]][i] ) }
  
  drift.estimate <- population1.se + .rmvnorm(iterations, sigma = D)
  
  drift.delta.z <- drift.estimate - population1.se
  drift.magnitude <- sqrt(rowSums(drift.delta.z^2))
  empirical.delta.z <- population2.se - population1.se
  observed.magnitude <- sqrt(rowSums(empirical.delta.z^2))
  drift.quantile <- quantile(drift.magnitude, c(0.025,0.975))
  observed.magnitude <- range(observed.magnitude)
  lista <- list("magnitude of morphological change expected through drift" = drift.quantile, 
                "observed magnitude" = observed.magnitude)
  return(lista)
}
