##########################################################################################
################################### Codes for model mechanism ############################
##########################################################################################



library(diagram)
library(Cairo);
library(grid)
library(ggsci)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#setwd('/home/spatial-r/E盘/Project/Mumps/写作模板/figures')

size.state <- 1.5
size.true <- 1.2

cairo_pdf("Figures/flowchart.pdf",width = 10, height = 8, 
          family = "serif",bg="transparent",pointsize = 15)

openplotmat()

par(mai =c(0,0,0,0),mar =c(0,0,0,0))

elpos <- coordinates(c(4)); elpos[,2] <- elpos[,2] -0.2
elpos[2:3,1] <- elpos[2:3,1]-0.03
elpos[1,1] <- elpos[1,1]+ 0.02
elpos[4,1] <- elpos[4,1]- 0.05

elpos_add <- matrix(c(elpos[3,1],elpos[3,1],0.1,0.5),nrow = 2)
elpos <- rbind(elpos,elpos_add)

for (i in 1:(3)){
  straightarrow(from = c(elpos[i, 1],elpos[i, 2]),
                to  = c(elpos[i + 1, 1],elpos[i + 1, 2]),
                lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)
}

straightarrow(from = c(elpos[2, 1]+0.082,elpos[2, 2]),to  =  c(elpos[5, 1]-0.08,elpos[5, 2]),
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)
straightarrow(from = c(elpos[2, 1]+0.082,elpos[2, 2]),to  = c(elpos[6, 1]-0.08,elpos[6, 2]),
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)

straightarrow(from = c(elpos[5, 1]+0.082,elpos[5, 2]),to  =  c(elpos[4, 1]-0.08,elpos[4, 2]),
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)
straightarrow(from = c(elpos[6, 1]+0.082,elpos[6, 2]),to  = c(elpos[4, 1]-0.08,elpos[4, 2]),#lty = 2,
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)


# curvedarrow(from = c(elpos[1, 1],elpos[1, 2]),
#             to  = c(elpos[4, 1],elpos[4, 2]), lwd = 1, 
#             arr.pos = 0.5, curve = -0.3, dr = 0.01,arr.length = 0.2,arr.width = 0.15)

labels <- c("S","E","IH","R","IU","IQ")

for (i in c(3,5,6)){
  textellipse(elpos[i,], 0.06, lab = labels[i], box.col = add.alpha("red",0.5),
              shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.08,radx = 0.08)
}

textellipse(elpos[1,], 0.08, lab = labels[1], box.col = add.alpha("red",0.5),
            shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.08,radx = 0.08)
textellipse(elpos[2,], 0.08, lab = labels[2], box.col = add.alpha("red",0.5),
            shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.08,radx = 0.08)
textellipse(elpos[4,], 0.08, lab = labels[4], box.col = add.alpha("red",0.5),
            shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.08,radx = 0.08)


curvedarrow(from = c(0.05,elpos[1,2]), to = c(0.95,elpos[1,2]),lty = 1,
            lcol = add.alpha("red",0.5),
            arr.pos = 0.5, curve = -0.35, dr = 0.01,arr.length = 0,arr.width = 0)
curvedarrow(from = c(0.05,elpos[1,2]), to = c(0.95,elpos[1,2]),lty = 1,
            lcol = add.alpha("red",0.5),
            arr.pos = 0.5, curve = 0.35, dr = 0.01,arr.length = 0,arr.width = 0)


#curvedarrow(from = c(0.15,0.78), to = c(0.88,0.78),lty = 2,
#            arr.pos = 0.5, curve = -0, dr = 0.01,arr.length = 0,arr.width = 0)
#curvedarrow(from = c(0.15,0.63), to = c(0.88,0.63),lty = 2,
#            arr.pos = 0.5, curve = -0, dr = 0.01,arr.length = 0,arr.width = 0)


textellipse(c(0.2,0.9), 0.08, lab = "Region B", box.col = add.alpha("dodgerblue",0.3),
            lcol = add.alpha("dodgerblue",0.3),
            shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.1,radx = 0.1)

straightarrow(from = c(0.22,0.75),to  =  c(0.25,0.64),lcol = add.alpha("dodgerblue",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)
straightarrow(from = c(0.28,0.67),to  =  c(0.25,0.77),lcol = add.alpha("red",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)


straightarrow(from = c(0.35,0.85),to  =  c(0.65,0.85),lcol = add.alpha("dodgerblue",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)

straightarrow(from = c(0.65,0.92),to  =  c(0.35,0.92),lcol = add.alpha("darkmagenta",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)


#curvedarrow(from = c(0.15,0.92), to = c(0.88,0.92),lty = 2,lcol = add.alpha("dodgerblue",0.3),
#            arr.pos = 0.5, curve = -0.2, dr = 0.01,arr.length = 0,arr.width = 0)

#curvedarrow(from = c(0.15,0.63), to = c(0.88,0.63),lty = 2,lcol = add.alpha("darkmagenta",0.3),
#            arr.pos = 0.5, curve = -0, dr = 0.01,arr.length = 0,arr.width = 0)

textellipse(c(0.8,0.9), 0.08, lab = "Region C,D,...", box.col = add.alpha("darkmagenta",0.3),
            lcol = add.alpha("darkmagenta",0.3),
            shadow.col = "white", shadow.size = 0.005, cex = 1,rady = 0.10,radx = 0.10)

straightarrow(from = c(0.74,0.77),to  =  c(0.7,0.67),lcol = add.alpha("darkmagenta",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)
straightarrow(from = c(0.74,0.64),to  =  c(0.78,0.75),lcol = add.alpha("red",0.5),
              lwd = 4, arr.pos = 0.8, arr.length = 0.5,arr.width = 0.5)

textempty(c((elpos[1,1]+elpos[2,1])/2+0.1,elpos[1,2]+0.2),lab = "Region A")

textempty(c(0.12,0.71),lab = "Population \n movement",srt = 0)
textempty(c(0.5,0.98),lab = "Population movement",srt = 0)

# textplain(c(f.text.pos[1,1],f.text.pos[1,2] + 0.24), 0.08, lab = expression(lambda(t)), cex = size.true)
# textplain(c((elpos[2,1] + elpos[3,1])/2,f.text.pos[1,2] + 0.24), 0.08, lab = expression(sigma * (1-omega) * q(t)), cex = size.true)
# textplain(c((elpos[2,1] + elpos[3,1])/2,f.text.pos[1,2] + 0.34), 0.08, lab = expression(sigma * (1-omega) * (1-q(t))), cex = size.true)
# textplain(c((elpos[2,1] + elpos[3,1])/2,f.text.pos[1,2] + 0.13), 0.08, lab = expression(sigma * omega), cex = size.true)
# textplain(c((elpos[3,1] + elpos[4,1])/2,f.text.pos[1,2] + 0.22), 0.08, lab = expression(gamma), cex = size.true)
# textplain(c((elpos[3,1] + elpos[4,1])/2,f.text.pos[1,2] + 0.34), 0.08, lab = expression(gamma), cex = size.true)
# textplain(c((elpos[3,1] + elpos[4,1])/2,f.text.pos[1,2] + 0.13), 0.08, lab = expression(gamma), cex = size.true)


dev.off()

