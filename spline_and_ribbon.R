library(ggplot2)
library(splines)
library(quantreg)
library(reshape2)
library(dplyr)


num_points = 100
num_samples = 50

x = seq(0, 4*pi, length.out=num_points)

sim <- lapply(1:num_samples, function(f) {
    x = runif(num_points, 0, 4*pi)
    y = sin(x) + rnorm(num_points, 0, 0.4)
    data.frame(x=x, y=y)
})

sim.df <- do.call(rbind, sim)
actual = data.frame(x=x, y=sin(x))

ggplot(sim.df, aes(x=x, y=y)) +
geom_point(alpha=0.7) +
geom_line(data=actual, colour='blue', size=1.5)

nq = 50 # Numbre of quantiles
qq = seq(0,1, length.out=nq)

m1 = rq(y ~ ns(x,10), data=sim.df, tau=qq)

xvals = seq(min(sim.df$x), max(sim.df$x), length.out=100)
rqs = data.frame(x=xvals, predict(m1, newdata=data.frame(x=xvals)))
names(rqs) = c("x", paste0("p",100*qq))

dat1 = rqs[, -length(rqs)]
names(dat1)[-1] = paste0(names(dat1)[-1])
dat2 = rqs[, -2]
names(dat2)[-1] = paste0(names(dat1)[-1])

dat1 = melt(dat1, id.var="x")
names(dat1) = c("x","group","min")
dat2 = melt(dat2, id.var="x")
names(dat2) = c("x","group1","max")

dat = bind_cols(dat1, dat2)

haha <- ggplot() +
  geom_point(data=sim.df, aes(x,y), alpha=0.1, size=0.5, colour="red") +
  geom_ribbon(data=dat, aes(x=x, ymin=min, ymax=max, group=group, alpha=group),
          fill="blue", lwd=0, show.legend=FALSE) +
  theme_bw() +
  scale_alpha_manual(values=c(seq(0.05,0.9,length.out=floor(0.5*length(qq))),
                              seq(0.9,0.05,length.out=floor(0.5*length(qq)))))
