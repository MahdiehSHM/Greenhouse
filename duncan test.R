
#https://rcompanion.org/handbook/G_09.html


library(agricolae)

# follow and repeat this example for all of your columns
#DATA INPUT

GH.DATA<-read.csv("Cucumis.csv",header = TRUE)
?duncan.test()
#do not report this ANOVA results! this is just for the duncan text! 
# you should report the ANOVA table with the interactions that you already have 
anova.1<-aov(Lshoot~treatment.code,data=GH.DATA)

output <- duncan.test(anova.1,"treatment.code", group=TRUE)
# in the following plot groups are show with different colors
# same colors means that those treatments are not sttatistically different fom eachother
plot(output,horiz=TRUE,las=1)
print(out$groups)
# end