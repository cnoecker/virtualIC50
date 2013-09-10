# Mehmet Gonen (mehmet.gonen@aalto.fi)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science

logdet <- function(Sigma) {
    2 * sum(log(diag(chol(Sigma))))
}

repmat <- function(M, row, column) {
    kronecker(matrix(1, row, column), M)
}
