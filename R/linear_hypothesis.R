#' Testiranje linearnih hipoteza
#'
#' Funkcija za testiranje hipoteza oblika C*beta = gamma
#'
#' @param model linearni model za koji se testira
#' @param C matrica C
#' @param gamma vektor gama
#' @return Ispisuje vrednost F statistike, odgovarajucih stepeni slobode, kao i
#'   p vrednost testa i vraca listu sa pomenutim vrednostima.
#' @export
linear_hypothesis <- function(model, C, gamma) {
  # Prvo racunamo potrebne vrednosti za statistiku Q
  beta <- model$coefficients
  X <- model.matrix(model)
  cbg <- C %*% beta - gamma
  XtXi <- solve(t(X) %*% X)
  Q <- t(cbg) %*% solve(C %*% XtXi %*% t(C)) %*% cbg

  # Racunamo SSE
  SSE <- sum(model$residuals^2)
  n <- nrow(X)
  p <- ncol(X) - 1
  m <- nrow(C)

  # Test statistika je sledeca...
  Fstat <- (Q / m) / (SSE / (n-p-1))

  # Racunamo p vrednost na osnovu F raspodele
  pval <- 1 - pf(Fstat, m, n-p-1)

  # belezimo stepene slobode
  df <- c(m, n-p-1)

  # Ispisujemo rezultat testiranja
  cat(sprintf("F-statistic: %.4f on %d and %d DF,  p-value: %.4f\n",
              Fstat, df[1], df[2], pval))

  ret_val <- list("statistic" = Fstat,
                  "p.value" = pval,
                  "df" = df)
  invisible(ret_val) # vracamo listu ret_val, ali je ne prikazujemo odmah
}
