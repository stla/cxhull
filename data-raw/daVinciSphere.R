# da Vinci's 72-sided sphere (septuaginta)

n <- 3L
p <- 4L * n
polyg1 <- vapply(seq_len(p), function(i) {
  c(cos(2*pi*i/p), sin(2*pi*i/p), 0)
}, numeric(3L))

R <- rgl::rotationMatrix(angle = 2*pi/p, x = 0, y = 1, z = 0)[-4L, -4L]

polyg2 <- R %*% polyg1[, -c(3L, 9L)]
polyg3 <- R %*% polyg2
polyg4 <- R %*% polyg3
polyg5 <- R %*% polyg4
polyg6 <- R %*% polyg5

daVinciSphere <- t(cbind(polyg1, polyg2, polyg3, polyg4, polyg5, polyg6))
