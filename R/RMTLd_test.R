#' @name RMTLd_test
#' @title Comparing restricted mean survival time
#' @description Performs two-sample comparisons using the restricted mean time lost(RMTL) as a summary measure of the cumulative incidence function under competing risks
#' @usage  RMTLd_test(time, status, group, alpha = 0.05, tau = NULL)
#' @importFrom survival survfit survdiff Surv
#' @importFrom cmprsk cuminc timepoints
#' @importFrom stats pchisq pnorm qnorm
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 1 = event of interest, 2 = competing event and 0 = right censored.
#' @param group The group indicator for comparison. The elements of this vector take either 1 or 0. Normally, 0 = control group, 1 = active treatment group.
#' @param tau  A scaler value to specify the truncation time point for the RMST calculation.
#' The defaultis the minimum of the maximum follow-up time of two groups
#' @param alpha The default is 0.05. (1-\code{alpha}) confidence intervals are reported.
#' @return an object of class RMTLd_test.
#' @return \item{tau}{the truncation time used in the analyses}
#' @return \item{Note}{a note regarding the truncation time}
#' @return \item{Base}{RMTL results in two groups, including RMTL and corresponding variance in group 0 and 1}
#' @return \item{Estimation}{RMSTd results between groups, including RMTLd and corresponding (1-\code{alpha}) confidence intervals.}
#' @return \item{Test}{Results of the log-rank test, Gray test and RMTLd test analyses.}
#' @export
#' @examples
#' library(crRMTL)
#' data(Simdata)
#' # The time point is set as the minimum of the maximum follow-up time of two groups when tau = NULL.
#' RMTLd_test(Simdata$time, Simdata$status, Simdata$group, alpha = 0.05, tau = NULL)
#' # The time point is set as tau = 3.5.
#' RMTLd_test(Simdata$time, Simdata$status, Simdata$group, alpha = 0.05, tau = 3.5)
#'
RMTLd_test <- function(time, status, group, alpha = 0.05,tau = NULL){

  ddd <- table(time, status, group)
  judge <- as.data.frame(table(status, group))

  if (dim(ddd)[3] != 2)
    stop("RMTLd test is for two groups")
  if (dim(ddd)[2] > 3)
    stop("All competing risks should be grouped under code 2")
  if (dim(ddd)[2] == 1  )
    stop("Either all observations are censored or \there is only one type of
         event and no censor observations")
  if (dim(ddd)[2] == 2 & sum(judge$status == 0) != 0)
    stop("There is no competing event or interest event")
  if (judge[judge$status == 2, 3][1] == 0 | judge[judge$status == 2, 3][2] == 0)
    stop("There is no competing event at least one group")
  if (judge[judge$status == 1, 3][1] == 0 | judge[judge$status == 1, 3][2] == 0)
    stop("There is no interest event at least one group")

  s0 <- 1 * (status == 1 | status == 2)
  s1 <- 1 * (status == 1)
  s2 <- 1 * (status == 2)
  d <- data.frame(group, time, status, s0, s1, s2)
  d0 <- d[d$group == 0,]
  d1 <- d[d$group == 1,]
  n1 <- table(group)[[1]]
  n2 <- table(group)[[2]]

  tau_max <- min(max(d0$time[d0$status != 2]), max(d1$time[d1$status != 2]))
  if(!is.null(tau)){
    if(tau <= tau_max){
      NOTE <- paste("The truncation time: tau =", tau, " was specified.")
    }
    if(tau > tau_max){
      stop(paste("The truncation time, tau, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(tau_max, digits=3)))
    }
  }
  if(is.null(tau)){
    tau <- tau_max
    NOTE <-(paste("The truncation time, tau, was not specified. Thus, the default tau (the minimum of the largest observed time on each of the two groups)", round(tau_max, digits=3)," is used."))
  }

  inner <- function(data,tau){

    sur_all <- survfit(Surv(time, s0) ~ 1, data = data)
    sur_int <- survfit(Surv(time, s1) ~ 1, data = data)

    index1 <- grep("TRUE", (sur_all[["n.event"]] != 0), value = F)
    index2 <- (sur_all[["time"]][index1] < tau)

    point <- sur_all[["time"]][index1][index2]
    t_fin <- c(point, tau)

    cif <- timepoints(cuminc(d$time, d$status), t_fin)$est[1,]
    cif1 <- timepoints(cuminc(data$time, data$status), t_fin)$est[1,]
    cif2 <- timepoints(cuminc(data$time, data$status), t_fin)$est[2,]

    e <- sur_int[["n.event"]][index1][index2]
    e_all <- sur_all[["n.event"]][index1][index2]

    r <- sur_all[["n.risk"]][index1][index2]
    s <- sur_all[["surv"]][index1][index2]
    N <- cumsum(sur_all[["n.event"]])[index1][index2]
    N1 <- cumsum(e)

    rmtl <- diff(t_fin)*cif1[-length(cif1)]
    R_L <- sum(rmtl)

    var.bk <- diff(c(0, cif1[-length(cif1)])) * ((tau - point) * (1 - cif2[-length(cif2)])- rev(cumsum(rev(rmtl)))) ^ 2 / s / r +
      diff(c(0, cif2[-length(cif2)])) * ((tau - point) * cif1[-length(cif1)] - rev(cumsum(rev(rmtl)))) ^ 2 / s / r
    var.bk[length(var.bk)] <- ifelse(s[length(s)] == 0, 0, var.bk[length(var.bk)])
    var.BK <- sum(var.bk)

    output <- c(tau, R_L, var.BK)
    output
  }
  G1 <- inner(d0, tau)
  G2 <- inner(d1, tau)

  out1 <- matrix(0, 2, 3)
  out1[1, ] <- G1
  out1[2, ] <- G2
  rownames(out1) <- c("group=0", "group=1")
  colnames(out1) <- c("tau", "RMTL", "Var")

  RMTLd <- (out1[2,2] - out1[1,2])
  out2 <- matrix(0, 1, 3)
  out2[1,] <- c(RMTLd, RMTLd - qnorm(1 - alpha / 2) * sqrt(out1[1, 3] + out1[2, 3]),
                RMTLd + qnorm(1 - alpha / 2) * sqrt(out1[1, 3] + out1[2, 3]))
  rownames(out2) <- c("RMTLd")
  colnames(out2) <- c("Statistics", "(1-alpha)%L", "(1-alpha)%U")


  z1 <- survdiff(Surv(d$time, d$s1) ~ group)$chisq
  z2 <- cuminc(d$time, d$status, d$group, cencode = 0)$Tests[1, ]
  z3 <- (out1[2, 2] - out1[1, 2]) / sqrt(out1[1, 3] + out1[2, 3])
  out3 <- matrix(0, 3, 2)
  out3[1, ] <- c(z1, 1 - pchisq(z1, 1))
  out3[2, ] <- c(z2[1], z2[2])
  out3[3, ] <- c(z3, 2 * (1 - pnorm(abs(z3))))
  rownames(out3) <- c("Logrank", "Gray", "RMTLd")
  colnames(out3) <- c("Z","P-value")
  z <- list()
  z$Note <- NOTE
  z$Base <- out1
  z$Estimation <- out2
  z$Test <- out3
  z
}
