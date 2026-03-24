
####################
load("cox_single.RData")

i <- 99
plot(density(rwmh_chain[[1]][, i]), main = "Posterior distribution of beta_100", xlab = "beta_1")
lines(density(bf_chain[[1]][, i]), col = "blue")

mean(ess(bf_chain[[1]]))
min(ess(bf_chain[[1]]))/bf_time[3]
mean(bf_chain[[2]])
multiESS(bf_chain[[1]], r = 1)




rwmh_chain[[2]]
mean(ess(rwmh_chain[[1]]))
min(ess(rwmh_chain[[1]]))/rwmh_time[3]
mean(rwmh_chain[[2]])
multiESS(rwmh_chain[[1]], r = 1)

plot(density(rwmh_chain[[1]][, 100]), col = "red", main = "Posterior distribution of beta_100", xlab = "beta_1")
lines(density(bf_chain[[1]][, 100]), col = "blue")
# plot(density(bf_chain[[1]][, 99]), col = "blue", main = "Posterior distribution of beta_100", xlab = "beta_1")
# lines(density(rwmh_chain[[1]][, 100]), col = "red")
# legend("topright", legend = c("RWMH", "Exact proposal"), col = c("blue", "red"), lty = 1)
