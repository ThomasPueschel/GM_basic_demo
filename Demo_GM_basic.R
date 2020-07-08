#Basic GM demo 

#load packages
library(jpeg)
library(png)
library(scatterplot3d)
library(Morpho)
library(geomorph)


#2D
#reading an image and interacting with the image
oo<-readPNG("pusheen.png")
plot(c(1, dim(oo)[2]), c(1, dim(oo)[1]), type = "n", xlab = "", ylab = "", asp=1)
a<-rasterImage(oo, 0, 0, dim(oo)[2], dim(oo)[1])
#getting coordinates from the image
coord<-locator(5,type="p",bg=5, pch=21)
coord

# basic data interaction
dat<-scan("TURTLE.DTA", skip=8, what="character")
turt1<-matrix(as.numeric(dat[1:(74*3)]), 74,3, byrow=T)
#use scatterplot3d to draw the configuration
scatterplot3d(turt1, scale.y=1)
plot(turt1[,1], turt1[,2], asp=1)
#another way to visualize x, y and z dimensions
pairs(turt1, asp=1)
link<-c((1:25)*2-1,1,NA,(1:25)*2,2,NA,1,2,NA,3,4,NA,5,6,NA,7,8,NA,9,10,NA,11,12,NA,13,14,NA,15,16,NA,17,18,NA,19,20,NA,21,22,NA,23,24,NA,25,26,NA,27,28,NA,29,30,NA,31, 32,NA,33,34,NA,35,36,NA,37,38,NA,39,40,NA,41,42,NA,43,44,NA,45,46,NA,47,48,NA,49,50, NA,51,61:74,60,51,NA,61,74,NA,63,72,NA,70,65,NA,67,24,NA,28,68,NA,66,54,NA,64,53,NA,62, 52,NA,73,59,NA,71,58,NA,69,57)
plot(turt1[link,1:2],type="b", asp=1)   
scatterplot3d(turt1[link,1], turt1[link,2], turt1[link,3], type="b", angle=120)
#the rgl library
plot3d(turt1[,1], turt1[,2], turt1[,3], aspect="iso")
#draw spheres instead of points
plot3d(turt1[,1], turt1[,2], turt1[,3], aspect="iso",type="s",col="purple", box=F, axes=F)
#adjust the size of spheres
plot3d(turt1[,1], turt1[,2], turt1[,3], aspect="iso",type="s",size=0.6,col="red", box=F, axes=F)
#add links
lines3d(turt1[link,1], turt1[link,2], turt1[link,3],col="blue", lwd=3)
#identify a landmark in 3d
identify3d(turt1, buttons="right", n=1, adj=2)
# Another example
#read PLY files and plot mesh files
#load the geomorph library
adapis<-read.ply("Adapis_parisiensis_left_ear.ply")
#clear the rgl scene
rgl.clear()
#use two colors for vertices
shade3d(adapis, color=c("magenta", "cyan",'yellow'))
# another example
data(nose)
shade3d(shortnose.mesh, color="pink")

#
#digitise 3d landmarks  #https://rdrr.io/cran/geomorph/man/digitsurface.html
## .ply has to be in ASCII format not binary
x <-read.ply("surface.for.specimen1.ply", ShowSpecimen=TRUE,addNormals = F) #read ply file and generate individual data
colors='green'
x$material=list(color=colors)
digit.fixed(x, fixed=5, index = FALSE, ptsize =1, center = TRUE) # landmark collection


# some examples
#getting familiar with GM data
data(plethodon)
#
oo<-readPNG("plethodon.png")
plot(c(1, dim(oo)[2]), c(1, dim(oo)[1]), type = "n", xlab = "", ylab = "", asp=1)
a<-rasterImage(oo, 0, 0, dim(oo)[2], dim(oo)[1])
# check coordinates 1st individual
plethodon$land[,,1] 
# check covariates
plethodon$species
# Procrustes
Y <- gpagen(plethodon$land)
plot(Y)
plotAllSpecimens(Y$coords,links=plethodon$links)
Y$coords #shape variables
Csize<-Y$Csize # Centroid size
#
#VISUALIZATION
ref <- mshape(Y$coords) #creates a reference from the average of the GPA specimens
#Principal Components Analysis
PCA<-plotTangentSpace(Y$coords, groups=as.factor(paste(plethodon$species, plethodon$site)))
#
#plot shape differences between a reference and target specimen
plotRefToTarget(ref, Y$coords[,,39], method="TPS")
#can magnify the differences seen using mag=#
plotRefToTarget(ref, Y$coords[,,39], mag=3, method="TPS") #magnifies by 3
#by changing the method to "vector" instead of "TPS" can see vector displacements rather than a thin plate spline
#a similar change is induced by switching the method to "points" to show reference and landmark points
plotRefToTarget(ref, Y$coords[,,39], mag=3, method="points", links=plethodon$links) 
#

#3d
#for a "fun diversion", see how to use a 3D dataset with semilandmarks
#load in the scallop data and check out the landmarks and semilandmarks. This dataset has semilandmarks on both curves and surfaces, differing from the demo dataset using salamanders (only landmarks)
data(scallops)
scallops$curvslide # Matrix defining which points are semilandmarks (middle column) and in which directions they slide (columns 1 [before] vs. 3 [after])
scallops$surfslide # Matrix (1 column) defining which points are semilandmarks to slide over the surface
#Procrustes Distance for the semilandmarks
Sc <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
#thin plate spline
scallopref <- mshape(Sc$coords)
plotRefToTarget(scallopref, Sc$coords[,,1], method="TPS", mag=3)
#plot 3D specimen with landmarks and semilandmarks
data(scallopPLY)
ply <- scallopPLY$ply
digitdat <- scallopPLY$coords
plotspec(spec=ply, digitspec=digitdat, fixed=16, centered=T)
#PCA
PCA_Sc<-plotTangentSpace(Sc$coords)
# plot PC extremes
PC1max<-plotRefToTarget(scallopref,PCA_Sc$pc.shapes$PC1max, mesh = ply, method="surface", mag = 3)
PC1min<-plotRefToTarget(scallopref,PCA_Sc$pc.shapes$PC1min, mesh = ply, method="surface", mag = 3)


############################################################
#Collecting 2D data
############################################################
# import image showing the landmark map for reference in another window

# landmarks 
filelist <- list.files(pattern = ".jpg") # open files to be digitised

# LANDMARKS: Data collection
# ----------------------------------------------------------------------------------------------------------

# set scale # e.g. 1 cm = scale=1
# collect lmks creating .tps file called 2D_data.tps
# we'll collect 3 landmarks
digitize2d(filelist, nlandmarks=3, scale=1, tpsfile = "2D_data.tps",verbose = T)

# ANALYSIS
# ----------------------------------------------------------------------------------------------------------

#import lmks ## 
data <- readland.tps("2D_data.tps", specID = "ID")
#Procrustes superimposition; this step removes the differences due to translation, rotation and scale
data.gpa <- gpagen(data)
plot(data.gpa)
two.d<-two.d.array(data.gpa$coords)
#Save the Procrustes residuals for further analysis
data.Y <- data.gpa$coords 
data.PCA <- plotTangentSpace(data.Y,label=row.names(two.d)) # PCA
#
data.links <- rbind(c(1,2),c(2,3))
# mean shape
ref <- mshape(data.Y) 
layout(1)
gp1 <- gridPar(pt.bg=c(rep("magenta",3)),txt.cex = c(rep(1, 3)))
#plot mean shape
plotRefToTarget(ref, ref, method="points",links=data.links, gridPars=gp1, label = T)
## Shape changes 
layout(mat=matrix(c(1:4), ncol = 2,byrow = T))
# sets the margins 
#c(bottom, left, top, right) 
par(mar=c(1, 1, 1, 1))  
# TPS grids
gp <- gridPar(grid.col = "lightblue", tar.pt.bg=c(rep("magenta",3)))
# shape changes 
# change mag to increase or decrease the change factor for visualization purposes
# default mag=1
plotRefToTarget(ref, data.PCA$pc.shapes$PC1min, links=data.links ,gridPars=gp, mag=1); title("PC1min")
plotRefToTarget(ref, data.PCA$pc.shapes$PC1max, links=data.links,gridPars=gp, mag=1); title("PC1max")
plotRefToTarget(ref, data.PCA$pc.shapes$PC2min, links=data.links,gridPars=gp, mag=1); title("PC2min")
plotRefToTarget(ref, data.PCA$pc.shapes$PC2max, links=data.links,gridPars=gp, mag=1); title("PC2max")




