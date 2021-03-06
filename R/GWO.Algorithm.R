#############################################################################
#
#  This file is a part of the R package "metaheuristicOpt".
#
#  Author: Iip
#  Co-author: -
#  Supervisors: Lala Septem Riza, Eddy Prasetyo Nugroho
#  Copyright (c) Department of Computer Science Education, Universitas Pendidikan Indonesia.
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This is the internal function that implements Grey Wolf Optimizer 
#' Algorithm. It is used to solve continuous optimization tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by Mirjalili in 2014, inspired by the behaviour of grey
#' wolf (Canis lupus). The GWO algorithm mimics the leadership hierarchy and hunting 
#' mechanism of grey wolves in nature. Four types of grey wolves such as alpha, beta, 
#' delta, and omega are employed for simulating the leadership hierarchy. 
#' In addition, the three main steps of hunting, searching for prey, encircling prey, 
#' and attacking prey, are implemented.
#' 
#' @title Optimization using Grey Wolf Optimizer
#'
#' @param FUN an objective function or cost function,
#'
#' @param optimType a string value that represent the type of optimization.
#'        There are two option for this arguments: \code{"MIN"} and \code{"MAX"}.
#'        The default value is \code{"MIN"}, which the function will do minimization. 
#'        Otherwise, you can use \code{"MAX"} for maximization problem.
#'
#' @param numVar a positive integer to determine the number variable.
#'
#' @param numPopulation a positive integer to determine the number population.
#'
#' @param maxIter a positive integer to determine the maximum number of iteration.
#'
#' @param rangeVar a matrix (\eqn{2 \times n}) containing the range of variables, 
#'        where \eqn{n} is the number of variables, and first and second rows
#'        are the lower bound (minimum) and upper bound (maximum) values, respectively. 
#'        If all variable have equal upper bound, you can define \code{rangeVar} as 
#'        matrix (\eqn{2 \times 1}).
#'
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{metaOpt}}
#' 
#' @examples
#' ################################## 
#' ## Optimizing the sphere function
#' 
#' # define sphere function as objective function
#' sphere <- function(X){
#'     return(sum(X^2))
#' }
#' 
#' ## Define parameter 
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#' 
#' ## calculate the optimum solution using Ant Lion Optimizer 
#' best.variable <- GWO(sphere, optimType="MIN", numVar, numPopulation=20, 
#'                  maxIter=100, rangeVar)
#' 
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(best.variable)
#' 
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable 
#'         and \code{vn} is value of \code{n-th} variable.
#' 
#' @references
#' Seyedali Mirjalili, Seyed Mohammad Mirjalili, Andrew Lewis, Grey Wolf Optimizer, 
#' Advances in Engineering Software, Volume 69, 2014, Pages 46-61, ISSN 0965-9978, 
#' http://dx.doi.org/10.1016/j.advengsoft.2013.12.007.
#' @export

GWO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
	# calculate the dimension of problem if not specified by user
	dimension <- ncol(rangeVar)

	# parsing rangeVar to lowerBound and upperBound
	lowerBound <- rangeVar[1,]
	upperBound <- rangeVar[2,]
	
	# if user define the same upper bound and lower bound for each dimension
	if(dimension==1){
		dimension <- numVar
	}

	## convert optimType to numerical form
	## 1 for minimization and -1 for maximization
	if(optimType == "MAX") optimType <- -1 else optimType <- 1

	# generate initial population of wolf
	wolf <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
	
	# find the best position
	bestPos <- engineGWO(FUN, optimType, maxIter, lowerBound, upperBound, wolf)
	
	return(bestPos)
}

## support function for calculating best position with GWO algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param wolf population of wolf

engineGWO <- function(FUN, optimType, maxIter, lowerBound, upperBound, wolf){
	# calculate the wolf fitness
	wolfFitness <- calcFitness(FUN, optimType, wolf)

	# sort wolf location based on fitness value
	index <- order(wolfFitness)
	wolfFitness <- sort(wolfFitness)
	wolf <- wolf[index,]

	# set the current alpha, beta, and delta position
	alpha <- wolf[1,]
	Falpha <- wolfFitness[1]

	beta <- wolf[2,]
	Fbeta <- wolfFitness[2]

	delta <- wolf[3,]
	Fdelta <- wolfFitness[3]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# value a decreased linearly from 2 to 0
		a <- 2-t*((2)/maxIter)

		for (i in 1:nrow(wolf)){
			for (j in 1:ncol(wolf)) {
				# generate random number [0,1]
				r1 <- runif(1)
				r2 <- runif(1)

				A1 <- 2*a*r1-a
	            C1 <- 2*r2
	            
	            D_alpha <- abs(C1*alpha[j]-wolf[i,j])
	            X1 <- alpha[j]-A1*D_alpha
	                       
	            r1 <- runif(1)
	            r2 <- runif(1)
	            
	            A2 <- 2*a*r1-a
	            C2 <- 2*r2
	            
	            D_beta <- abs(C2*beta[j]-wolf[i,j])
	            X2 <- beta[j]-A2*D_beta
	            
	            r1 <- runif(1)
	            r2 <- runif(1)
	            
	            A3 <- 2*a*r1-a
	            C3 <- 2*r2
	            
	            D_delta <- abs(C3*delta[j]-wolf[i,j])
	            X3 <- delta[j]-A3*D_delta
	            
	            wolf[i,j]=(X1+X2+X3)/3
			}

			# check boundary for each dimension
			wolf[i,] <- checkBound(wolf[i,], lowerBound, upperBound)

			fitness <- optimType*FUN(wolf[i,])

			# update alpha, beta and delta
	        if(fitness<Falpha){ 
	            Falpha <- fitness
	            alpha <- wolf[i,]
	        }
	        
	        if(fitness>Falpha & fitness<Fbeta){
	            Fbeta <- fitness
	            beta <- wolf[i,]
	        }
	        
	        if(fitness>Falpha & fitness>Fbeta & fitness<Fdelta){
	            Fdelta <- fitness
	            delta <- wolf[i,]
	        }
		}
		
		# save the best fitness for iteration t
		curve[t] <- Falpha
		
		setTxtProgressBar(progressbar, t)
	}
	
	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="GWO", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(alpha)
}
