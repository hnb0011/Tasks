setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_11')
rnorm(100, mean=5,  sd = 2)
x <- rnorm(100, mean=5,  sd = 2)
y <- (((5*x)+2)+runif(100, min=0, max=0.1))
plot(x,y)
lm(x~y)
Q: Using rise over run, my slope~1/5. My y-int is at ~ -0.4. The mean and standard deviation help determine these numbers in the rnorm function. They set up parameters that the numbers cannot stray too far from.
{
z[[var]] <- mean(y*as.numeric(x[100]))
z <- runif(1)
y = (((x*z)+2) + runif(100, min =0, max =0.1))
}
print({
z[[var]] <- mean(y*as.numeric(x[100]))
z <- runif(1)
y = (((x*z)+2) + runif(100, min =0, max =0.1))
})
zz <- print({
z[[var]] <- mean(y*as.numeric(x[100]))
z <- runif(1)
y = (((x*z)+2) + runif(100, min =0, max =0.1))
})
plot(zz)

install.packages("meme")
library("meme")
vignette("meme", package="meme")
u <- system.file("angry8.jpg", package="meme")
meme(u, "357 lines of code complete", "3 missing parentheses")
