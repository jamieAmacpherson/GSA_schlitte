cat = read.table("time_entropies.dat")
cat1 = as.data.frame(cbind(as.data.frame(c(1: nrow(cat)) * 0.05), (cat$V1)))
names(cat1) = c("timestep", "entropy")

pdf("plot_entropies.pdf")
plot(cat1$timestep, cat1$entropy, xlab="Time (ns)", ylab = "S' (J/mol K)", panel.first=grid())

dev.off()
