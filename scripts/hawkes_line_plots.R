par(mfrow=c(2,1), xpd=TRUE)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(-10, -10, xlim=c(0,1000), ylim=c(0,1), axes=FALSE,
     ylab="", xlab="", main="")
text(1070,.975, "May 17")
clip(0,960,-.04,1)

abline(v=which(mml[1,]==1), lwd=1)
axis(1, at = seq(0,1000,120), labels = seq(0,500,60))
legend(x=950, y=.975, 
       legend=c(
                expression(paste(mu, " = ", .008)),
                expression(paste(gamma, " = ", .049)),
                expression(paste(alpha, " = ", .018)),
                expression(paste(beta, " = ", "0.020"))
                
       ), ncol=2, bty="n"

       )


plot(-10, -10, xlim=c(0,1000), ylim=c(0,1), axes=FALSE,
     ylab="", xlab="Time (seconds)", main="")
text(1070,.975, "June 22")

clip(0,960,-.04,1)

abline(v=which(babcock622[7,]==1), lwd=1)
axis(1, at = seq(0,1000,120), labels = seq(0,500,60))
legend(x=950, y=.975, 
       legend=c(
                expression(paste(mu, " = ", .011)),
                expression(paste(gamma, " = ", "0.010")),
                expression(paste(alpha, " = ", .002)),
                expression(paste(beta, " = ", "0.002"))
                
       ), ncol=2, bty="n"
)