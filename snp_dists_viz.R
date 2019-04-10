snp_dist = scan('1240k_snp_dists.txt', what = integer(), sep = '\n')

log2_d = log2(snp_dist + 1)
log2_pos = seq(0, 25, 5)
real_pos = (2^log2_pos) - 1
axis(side=1, at = log2_pos, labels = real_pos)

real_pos2 = c(0, 10, 100, 400, 1250, 32800)

log2_pos2 = log2(real_pos2 +1)

pdf('tmp.pdf', 6, 6)
plot(1:5, xlim = c(0,15), ylim = c(0,0.20), type = 'n', 
     xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', main = '', 
     axes = F)
grid(ny=10, nx = 0)
hist(log2_d, col = 'darkolivegreen', xaxt = 'n', xlab = '', main= '',
     xlim = c(0,15), br = 30, prob=T, add=T)
axis(side=1, at = log2_pos2, labels = real_pos2, las=3)
dev.off()




dens = density(log2_d)

plot(dens, col = 'darkolivegreen', xaxt = 'n', xlab = '', main= '',
     xlim = c(0,15))
x1 = min(dens$x)
x2 = max(dens$x)

with(density(log2_d), polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gray"))
axis(side=1, at = log2_pos2, labels = real_pos2, las=4)


q75 = quantile(log2_d, .75)
q95 = quantile(log2_d, .95)
dd <- with(dens,data.frame(x,y))
library(ggplot2)
pdf('tmp1.pdf', 6, 6)
qplot(x,y,data=dd,geom="line") +
  geom_ribbon(aes(ymax=y),ymin=0,fill="darkolivegreen",colour=NA,alpha=0.9) +
  scale_x_continuous(limits = c(0, 15), labels = real_pos2, breaks = log2_pos2) + 
  labs(x = 'Distance', y = 'Proportion', title = 'Distance between consecutive SNPs') +
  theme(title = element_text(color = 'gray4', face = 'bold', size = 15)
dev.off()
