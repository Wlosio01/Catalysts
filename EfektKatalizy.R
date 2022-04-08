
# Algorithm
HellwigForTwoX <- function(Y, listX)
{
  H <- 0 
  for (j in 1:2){
    # We count rj
    corXY <- cor(listX[[j]],Y)
    sum <- 0 
    for(k in 1:2){
      #we count sum of rij
      sum = sum + abs(cor(listX[[j]], listX[[k]]))
    }
    # Here we count H by add all h in loop
    H = H + (corXY^2)/sum
  }
  return(H)
}

X1 <- mtcars[[1]]
X2 <- mtcars[[2]]
X3 <- mtcars[[3]]
X4 <- mtcars[[4]]
X5 <- mtcars[[5]]
X6 <- mtcars[[7]]
Y <- mtcars[[6]]
# We put X to list in order
xlist <- list(X1,X2,X3,X4, X5, X6)

Catalysis <- function(Y,xlist){
  n <- length(xlist)
  # We create empty vector for X corelations with Y
  R0 <- c()
  # Vector for information about direction of corelation
  Rsign <- c()
  for(i in 1:n){
    if(cor(Y,xlist[[i]]) >= 0){
      Rsign <- c(Rsign,1)
    }
    else{
      Rsign <- c(Rsign,0)
    }
    R0 <- c(R0,abs(cor(Y,xlist[[i]])))
  }
  R1 <- order(R0)
  # We order also corelations of X 
  Rsign <- Rsign[R1]
  # If direction is positive we have 1, if negative -1
  Rsign[Rsign != 0] = 1
  Rsign[Rsign == 0] = -1
  # We sort xlist in order of R0
  xlist <- xlist[R1]
  R0 <- sort(R0)
  R <- matrix(nrow = n, ncol = n)
  # List with i, j index of catalysts
  catalystsIndex <- list()
  for(i in 1:n){
    for(j in 1:n){
      # We multiply cor by our vector with direction of corelations X
      if(i == j){
      R[i,j] <- 1}
      else{
        R[i,j] <- cor(xlist[[i]]*Rsign[[i]],xlist[[j]]*Rsign[[j]])
      }
      # We find catalysts
      if(j > i && (R[i,j] > (R0[i]/R0[j]) || R[i,j] < 0)){
        n1 <- length(catalystsIndex)
        catalystsIndex[[n1 +1]] <- c(i,j)
      }
    }
  }
  if(length(catalystsIndex) != 0){
  for(i in 1:length(catalystsIndex)){
    R0prim <- c(R0[catalystsIndex[[i]][1]],R0[catalystsIndex[[i]][2]])
    Rprim <- c(1,R[catalystsIndex[[i]][1],catalystsIndex[[i]][2]],R[catalystsIndex[[i]][1],catalystsIndex[[i]][2]],1)
    R0prim <- matrix(R0prim, nrow = 2)
    Rprim <- matrix(Rprim, ncol = 2, nrow = 2)
    R2 = (t(R0prim)%*%solve(Rprim))%*%R0prim;
    catalystsIndex[[i]][3] <- R2
    catalystsIndex[[i]][4] <- HellwigForTwoX(Y, list(xlist[[catalystsIndex[[i]][1]]],xlist[[catalystsIndex[[i]][2]]]))
    catalystsIndex[[i]][5] <- catalystsIndex[[i]][3] -catalystsIndex[[i]][4]
    # We change index of X to basic in first list
    catalystsIndex[[i]][1] <- R1[catalystsIndex[[i]][1]]
    catalystsIndex[[i]][2] <- R1[catalystsIndex[[i]][2]]
  }
  return(catalystsIndex)}
  else{
    print("Brak katalizatorów")
    return(NULL)
  }
}
# We put our pairs catalysts into list
catalysts <- Catalysis(Y,xlist)
#catalysts[[i]][1] and catalysts[[i]][2] shows indexes of X
#catalysts[[i]][3] shows R2
#catalysts[[i]][4] shows H
#catalysts[[i]][5] shows R2 - H

