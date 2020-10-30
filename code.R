#Nina Mortensen take home R exam 

#1. 

v1 <- c(-1.4,1.4,-0.5, -0.4)
v2 <- c(0,0.5, -1.7, 1.6) 
v3 <- c(-2.8, 2.3, 0.7, -2.4) 
v4 <- c(1, -2.9, -1.6, 3.0)
v5 <- c(2.4, -4.8, 0.6, 1.8)
v6 <- c(0.9, -1.9, 1.8, -1.0)
v7 <- c(0.8, -0.4, 3.5, -3.4) 

v <- cbind(v1,v2,v3,v4,v5,v6,v7)

v[1,] <- v[1,]/-1.4; v[1,]
v[2,] <- v[2,] - 1.4*v[1,]
v[3,] <- v[3,] + .5*v[1,]
v[4,] <- v[4,] + .4*v[1,]

v[2,] <- v[2,]*2
v[3,] <- v[3,] + 1.7*v[2,]
v[4,] <- v[4,] - 1.6*v[2,]

v[3,]<- v[3,]/-8.4171429
v[4,]<- v[4,] - 8.7942857*v[3,]

v[4,] <- v[4,]/-.0646639

v[2,] = v[2,]+3.8*v[3,]
v[1,] = v[1,]+.7143*v[3,]

v[1,] = v[1,] + .4798 * v[4,]
v[2,] = v[2,] + 1.1326 * v[4,]
v[3,] = v[3,] - .2283 * v[4,]

#final rref 
round(v, digits=2)

#the basis of the image are composed of the 1st, 2nd, 4th, 6th cols of the original matrix 
#aka v1, v2, v4, v6

# v1 <- c(-1.4,1.4,-0.5, -0.4)
# v2 <- c(0,0.5, -1.7, 1.6) 
# 
# v4 <- c(1, -2.9, -1.6, 3.0)
# 
# v6 <- c(0.9, -1.9, 1.8, -1.0)


#basis of kernel 
#To find a basis for the kernel, use the three nonpivotal columns.
#Vectors of the form (a1 b1 1 c1 0 d1 0) and (a2 b2 0 c2 1 d2 0) and (a3 b3 0 c3 0 d3 1) must be independent,
#because no linear combination can have 0 as its 3rd and 5th and 7th component.

#The top row of the row reduced matrix applied to (a1 b1 1 c1 0 d1 0) says that a1 = -2
#The second row of the row reduced matrix applied to (a1 b1 1 c1 0 d1 0) says that b1 = 1
#the third row c1 = 0 
#fourth, d1 = 0


k1 <- c(-2, 1, 1, 0, 0, 0, 0)  
#first basis vector for the kernel
round(v%*%k1, digits=3)     #yes, it is in the kernel

#The top row of the row reduced matrix applied to (a2 b2 0 c2 1 d2 0) says that a1 = 1
#The second row of the row reduced matrix applied to (a2 b2 0 c2 1 d2 0) says that b1 = 1
#the third row c1 = -1
#fourth, d1 = 0

k2 <- c(1, 1, 0, -1, 1, 0, 0)  
#first basis vector for the kernel
round(v%*%k2, digits=3)     #yes, it is in the kernel

#The top row of the row reduced matrix applied to (a3 b3 0 c3 0 d3 1) says that a1 = 0
#The second row of the row reduced matrix applied to (a3 b3 0 c3 0 d3 1) says that b1 = -1
#the third row c1 = 1
#fourth, d1 = -2

k3 <- c(0, -1, 0, 1, 0, -2, 1)  
#first basis vector for the kernel
round(v%*%k3, digits=3)     #yes, it is in the kernel


# the basis of the image spans all of R4 & therefore the matrix is onto; it however is not one-to-one
#since the kernel isn't empty 





#2. 
#part a) 

x1 <- rep(1,9); x1
x2 <- c(0,0.2,0.5, 0.9, 1.1,1.2, 1.6, 1.8, 2.1); x2
x3 <- x2*x2; x3
x4 <- x2*x3; x4

#Step 1: make the first unit vector by normalization
x_1 <- x1/Norm(x1); x_1; Norm(x_1)
#Step 2: convert w2 to a vector that is orthogonal to v1.
x <- x2 - (x2%.%x_1)*x_1; round(x%.% x_1)
#Then convert x to a unit vector
x_2 <- x/Norm(x); x_2
#Step 3: convert w3 to a vector that is orthogonal to both v1 and v2.
x_ <- x3 - (x3%.%x_1)*x_1 - (x3%.%x_2)*x_2; x_%.% x_1; x_%.% x_2
#Then convert x to a unit vector
x_3 <- x_/Norm(x_); x_3
#Step 4: convert w4 to a vector that is orthogonal to v1,v2,v3.
x__ <- x4 - (x4%.%x_1)*x_1 - (x4%.%x_2)*x_2 - (x4%.%x_3)*x_3; x__%.% x_1; x__%.% x_2; x__%.% x_3
#Then convert x to a unit vector
x_4 <- x__/Norm(x__); x_4

C <- cbind(x_1,x_2, x_3, x_4) ;C  #basis vectors are the columns
round(t(C)%*%C, digits = 6) 

#part b) 

y <- c(-0.4, -3.1, .3, .8, .9, .7, 3.0, 3.1, 5.9)
y1 <- (y%.%x_1)*x_1 + (y%.%x_2)*x_2 + (y%.%x_3)*x_3 + (y%.%x_4)*x_4; y1
y2 <- y - y1; y2
y1+y2 #final check 




#3.

c1 <- c(-0.69, 0.3, 1.58, -0.31, 1.1, -0.55, -1.43); c1
c2 <- c(0.96, 2.05, 0.03, -1.6, -1.67, 0.43, -0.18); c2
c3 <- c(0.8, -1.79, 0.35, 0.08, 1.17, -1.59, 0.98); c3
c4 <- c(0.66, 0.52, 0.47, -1.51, -2.35, 0.77, 1.44); c4
B<- cbind(c1,c2,c3,c4)
C<- t(B)%*%B; C

w <- c(1,0,0,0)
T <- cbind(w, C%*%w, C%*%C%*%w, C%*%C%*%C%*%w, C%*%C%*%C%*%C%*%w); T
rref(T)
#So A^4 - 37.6336A^3 +396.8666A^2-1413.2727A + 1342.87
p <- function(t) t^4 - 37.6336*t^3 +396.8666*t^2-1413.2727*t + 1342.87
curve(p(x), from = 1, to = 25); abline(h=0, col = "red")

f1 <- function(x) x^4 - 37.6336*x^3 +396.8666*x^2-1413.2727*x + 1342.87
ev_1 <- uniroot(f1, c(10,25))
ev_2 <- uniroot(f1, c(5,10))
ev_3 <- uniroot(f1, c(2,5))
ev_4 <- uniroot(f1, c(0,2))

#Since C has four distinct real roots 22.87052, 8.777701, 4.498348, 1.487048, we will get four eigenvectors.
I <- diag(c(1,1,1,1))
v0 <- (C - 8.777701*I)%*%(C - 4.498348*I)%*%(C - 1.487048*I)%*%w; v0
C%*%v0; 22.87052*v0      #eigenvector for eigenvalue 22.87052
v1 <- (C - 22.87052*I)%*%(C - 4.498348*I)%*%(C - 1.487048*I)%*%w; v1
C%*%v1; 8.777701*v1      #eigenvector for eigenvalue 8.777701
v2 <- (C - 8.777701*I)%*%(C - 22.87052*I)%*%(C - 1.487048*I)%*%w; v2
C%*%v2; 4.498348*v2      #eigenvector for eigenvalue 4.498348
v3 <- (C - 8.777701*I)%*%(C - 22.87052*I)%*%(C - 4.498348*I)%*%w; v3
C%*%v3; 1.487048*v3      #eigenvector for eigenvalue 1.487048

v_0 <- v0/Norm(v0); v_0
v_1 <- v1/Norm(v1); v_1
v_2 <- v2/Norm(v2); v_2
v_3 <- v3/Norm(v3); v_3

#check PDP_inv
P<- cbind(v_0, v_1, v_2, v_3)
P_<- solve(P)
D <- diag(c(22.87052, 8.777701, 4.498348, 1.487048))
P%*%D%*%P_
C







