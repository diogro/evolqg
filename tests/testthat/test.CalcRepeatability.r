test_that("CalcRepeatability returns correct values for balanced pairs",
          {
            local_seed(1)
            num.ind <- length(iris[,1])
            ID <- as.factor(rep(1:num.ind, 2))
            ind.data <- rbind(iris[,1:4], iris[,1:4] + array(rnorm(num.ind*4, 0, 0.1), dim(iris[,1:4])))
            models.list <- apply(ind.data, 2, function(vec) lm(vec ~ ID))
            models.list <- lapply(models.list, anova)
            rep.itself <- function(lm.model){
                msq <- lm.model$'Mean Sq' ## 1 entre, 2 dentro
                s2a <- (msq[1] - msq[2])/2
                s2a / (s2a + msq[2])
            }
            out <- sapply(models.list, rep.itself)
            expect_equal(CalcRepeatability(ID, ind.data), out)
          }
)

test_that("CalcRepeatability uses n0 for unbalanced designs",
          {
            # 3 individuals measured 3, 2 and 1 times (a = 3, N = 6)
            ID <- c(1, 1, 1, 2, 2, 3)
            y <- c(1, 2, 3, 6, 8, 10)
            # one-way ANOVA by hand
            grand <- mean(y)
            group.means <- tapply(y, ID, mean)
            n.i <- c(3, 2, 1)
            ss.among <- sum(n.i * (group.means - grand)^2)
            ss.within <- sum((y - group.means[as.character(ID)])^2)
            ms.among <- ss.among/2
            ms.within <- ss.within/3
            n0 <- (6 - sum(n.i^2)/6)/2
            expect_equal(n0, 11/6)
            s2a <- (ms.among - ms.within)/n0
            expected <- s2a/(s2a + ms.within)
            expect_equal(unname(CalcRepeatability(ID, cbind(y))), expected)
          })

test_that("CalcRepeatability accepts vectors and keeps column names",
          {
            ID <- rep(1:5, each = 2)
            y <- c(1, 1.1, 2, 2.2, 3, 2.9, 4, 4.3, 5, 5.1)
            expect_equal(CalcRepeatability(ID, y), CalcRepeatability(ID, cbind(y)) |> unname())
            out <- CalcRepeatability(ID, cbind(a = y, b = rev(y)))
            expect_named(out, c("a", "b"))
          })
