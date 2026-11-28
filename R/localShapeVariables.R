interpol.TPS <- function (x, y)
{
  distance <- Norm (x - y)
  if (length (x) == 2)
  {
    if (distance < .Machine $ double.eps)
      return (0)
    else
      return (- log (distance ^ 2) * (distance ^ 2))
  }
  if (length (x) == 3)
    return (distance)
}

interpol.L1.TPS <- function (y, X)
  return (.tps_kernel (sqrt (colSums ((t (X) - y)^2)), ncol (X)))

interpol.L2.TPS <- function (X, Y)
  return (.tps_kernel (as.matrix (dist (Y)), ncol (Y)))

# TPS kernel as a function of the distance: -r^2 log r^2 (0 at 0) in 2D, r in 3D
.tps_kernel <- function (d, m)
{
  if (m == 2)
    ifelse (d < .Machine $ double.eps, 0, - log (d ^ 2) * (d ^ 2))
  else if (m == 3)
    d
  else NULL
}

# Closed form jacobian of the TPS function theta.p at point p
.tps_jacobian <- function (p, Q, W, A)
{
  delta <- t (p - t (Q))                      # k x m, rows p - Q_i
  r2 <- rowSums (delta ^ 2)
  if (ncol (Q) == 2) {
    # d/dp of -r2 log(r2) is -(2 log(r2) + 2) (p - Q_i)
    coef <- ifelse (r2 < .Machine $ double.eps ^ 2, 0, - (2 * log (r2) + 2))
  } else {
    # d/dp of r is (p - Q_i) / r
    r <- sqrt (r2)
    coef <- ifelse (r < .Machine $ double.eps, 0, 1 / r)
  }
  A + t (W) %*% (coef * delta)
}

#' TPS transform
#'
#' Calculates the Thin Plate Spline transform between a reference shape and a target shape
#'
#' @param target.shape Target shape
#' @param reference.shape Reference shape
#'
#' @returns A list with the transformation parameters and a function that gives
#' the value of the TPS function at each point for numerical differentiation
#'
#' @author Guilherme Garcia
#' @export
#' @importFrom stats spline
TPS <- function (target.shape, reference.shape)
{
  ## calculates thin plate splines deformations between reference and target shapes
  theta.p <- function (p)
    return (b + (A %*% p) + t (W) %*% interpol.L1.TPS (p, Q))
  Q <- reference.shape
  P <- target.shape
  k <- nrow (Q)
  m <- ncol (Q)
  R <- interpol.L2.TPS (Q, Q)
  # print (R)
  one.k <- array (1, c (k, 1))
  zero.11 <- array (0, c (1, 1))
  zero.m1 <- array (0, c (m, 1))
  zero.mm <- array (0, c (m, m))
  L <- rbind (R, t (one.k), t (Q))
  L <- cbind (L, rbind (one.k, zero.11, zero.m1))
  L <- cbind (L, rbind (Q, t (zero.m1), zero.mm))
  P.zero <- rbind (P, t (zero.m1), zero.mm)
  spline <- solve (L, P.zero)
  W <- spline [1:k, ]
  b <- t (t (spline [1+k, ]))
  A <- t (spline [(2:(m+1))+k, ])
  return (list ('W' = W, 'A' = A,
                'b' = b, 'Q' = Q, 'P' = P, 'theta.p' = theta.p))
}


#' Midline rotate
#'
#' Returns the rotation matrix that aligns a specimen sagital line
#' to plane y = 0 (2D) or z = 0 (3D)
#'
#' @param X shape array
#' @param midline rows for the midline landmarks
#'
#' @returns Rotation matrix
#'
#' @author Guilherme Garcia
#' @export
Rotate2MidlineMatrix <- function (X, midline)
{
  ncl <- ncol (X)
  Xm <- na.omit (X [midline, ])
  Mm <- matrix (apply (Xm, 2, mean), byrow = TRUE, nrow = nrow (X), ncol = ncl)
  Xc <- X - Mm
  W <- na.omit (Xc [midline, ])
  RM <-svd (var (W))$v
  return (RM)
}

#' Local Jacobian calculation
#'
#' Calculates jacobians for a given interpolation in a set of points
#' determined from tesselation (as centroids of each tetrahedron defined, for now...)
#'
#' @param spline Thin plate spline calculated by the TPS function
#' @param tesselation matrix of landmarks.
#' @param ... Additional arguments to some function
#' @note Jacobians are calculated on the row centroids of the tesselation matrix.
#'
#' @returns array of jacobians calculated at the centroids
#'
#' @author Guilherme Garcia
#' @export
JacobianArray <- function (spline, tesselation, ...)
{
  if (length (list (...)) > 0)
    warning ("Arguments in '...' were used for numerical differentiation and are now ignored: ",
             "the jacobian is calculated analytically.")
  with (spline,
        {
          Q.tetra <- Q [tesselation, ]
          dim (Q.tetra) <- c (dim (tesselation), ncol (Q))
          Q.centroids <- apply (Q.tetra, 1, colMeans)
          if (!is.matrix (Q.centroids)) Q.centroids <- matrix (Q.centroids, nrow = ncol (Q))
          jacobs <- apply (Q.centroids, 2, function (p) as.vector (.tps_jacobian (p, Q, W, A)))
          dim (jacobs) <- c (ncol (Q), ncol (Q), ncol (Q.centroids))
          return (jacobs)
        })
}

#' Centered jacobian residuals
#'
#' Centers the log determinants of a set of jacobian matrices
#' describing a local aspect of shape deformation for a given set of volumes,
#' returning log determinants of deviations from the mean jacobian (Woods, 2003).
#' The mean of the log determinants is the log determinant of the Karcher
#' (log-Euclidean) mean jacobian, so the result is
#' \code{logdet(J_i) - mean(logdet(J))}.
#'
#' @param jacobArray Arrays of Jacobian calculated in the JacobianArray function
#'
#' @returns vector of centered log determinants of the jacobians, with mean zero
#'
#' @author Guilherme Garcia
#' @author Diogo Melo
#' @references Woods, Roger P. 2003. “Characterizing Volume and Surface
#' Deformations in an Atlas Framework: Theory, Applications, and Implementation.” NeuroImage 18 (3):769-88.
#' @export
Center2MeanJacobianFast <- function (jacobArray)
{
  log.det <- vapply (seq_len (dim (jacobArray) [3]),
                     function (i) as.numeric (determinant (jacobArray [, , i], logarithm = TRUE) $ modulus),
                     numeric (1))
  return (log.det - mean (log.det))
}

#' Local Shape Variables
#'
#' Calculates the local shape variables of a set of landmarks using the sequence:
#' - TPS transform between all shapes and the mean shape
#' - Jacobian of the TPS transforms at the centroid of rows of the landmarks in the tesselation argument
#' - Take the log determinant of the jacobians
#' - Mean center the log determinants (equivalent to centering the jacobians
#'   on their Karcher mean)
#'
#' @param gpa Procustes aligned landmarks.
#' @param cs Centoid sizes
#' @param landmarks unaligned landmarks. Ignored if gpa is passed. If both gpa
#' and cs are NULL, a generalized Procrustes analysis is run on the landmarks
#' (requires the Morpho package).
#' @param tesselation matrix of rows of the landmarks. The centroid of each row
#' is used to mark the position of the jacobians
#' @param run_parallel Logical. If computation should be paralleled. Use with
#' caution, can make things worse. Requires that at parallel back-end like doParallel
#' be registered
#' @param verbose Logical. If TRUE, progress messages are shown.
#' @returns List with TPS functions, jacobian matrices, local shape variables, mean shape, centroid sizes and individual IDs
#' @author Guilherme Garcia
#' @author Diogo Melo
#' @export
LocalShapeVariables <- function (gpa = NULL, cs = NULL, landmarks = NULL, tesselation,
                         run_parallel = FALSE, verbose = TRUE)
{
  if(is.null(gpa)){
    if(!is.null(cs))
      stop("gpa must be supplied when cs is supplied. Pass both gpa and cs, or neither (with landmarks).")
    if(is.null(landmarks))
      stop("Supply either gpa and cs, or the unaligned landmarks.")
    if(!requireNamespace("Morpho", quietly = TRUE))
      stop("Running the Procrustes alignment requires the Morpho package: install Morpho or supply gpa and cs.")
    if(verbose) message("Running GPA")
    gpa <- Morpho::ProcGPA(landmarks, CSinit = TRUE)
    cs <- apply(landmarks, 3, .csize)
    gpa <- gpa$rotated
    dimnames(gpa)[[3]] = dimnames(landmarks)[[3]]
  } else if(is.null(cs)){
    if(!is.null(landmarks)) {
      cs <- apply(landmarks, 3, .csize)
    } else {
      warning("cs not supplied, computed from gpa; if gpa was scaled to unit centroid size these are all 1.")
      cs <- apply(gpa, 3, .csize)
    }
  }
  mshape <- rowMeans(gpa, dims = 2)
  dimnames(mshape) <- dimnames(gpa)[1:2]
  ids <- dimnames(gpa)[[3]]
  n.ind <- dim(gpa)[3]
  tps <- .evq_lapply(lapply(seq_len(n.ind), function(i) gpa[, , i]),
                     TPS, reference.shape = mshape, parallel = run_parallel)
  names(tps) <- ids
  names(tps) <- as.character(seq_len(n.ind))
  if(verbose) message("tps done")
  jacobs <- .evq_lapply(tps, JacobianArray, tesselation = tesselation,
                        parallel = run_parallel)
  n.dim <- ncol(mshape)
  n.vol <- dim(jacobs[[1]])[3]
  jacobs <- array(unlist(jacobs, use.names = FALSE), c(n.dim, n.dim, n.vol, n.ind))
  jacobs <- aperm(jacobs, c(1, 2, 4, 3))
  dimnames(jacobs) <- list(as.character(seq_len(n.dim)), as.character(seq_len(n.dim)),
                           X1 = ids, as.character(seq_len(n.vol)))
  if(verbose) message("jacobs done")
  local <- vapply(seq_len(n.vol), function(v) Center2MeanJacobianFast(jacobs[, , , v]),
                  numeric(n.ind))
  local <- matrix(local, n.ind, n.vol,
                  dimnames = list(as.character(seq_len(n.ind)), as.character(seq_len(n.vol))))
  return (list ('tps' = tps,
                'jacobians' = jacobs,
                'local' = local,
                'reference' = mshape,
                'cs' = cs,
                'info' = ids))
}

# Centroid size
.csize <- function (x) sqrt (sum (sweep (x, 2, colMeans (x)) ^ 2))
