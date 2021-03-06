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
#' This is the internal function that implements Particle Swarm Optimization 
#' Algorithm. It is used to solve continuous optimization tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by Kennedy and Eberhart in 1995, inspired by
#' the behaviour of the social animals/particles, like a flock of birds in
#' a swarm. The inertia weight that proposed by Shi and Eberhart is used to
#' increasing the performance of PSO.
#' 
#' @title Optimization using Prticle Swarm Optimization
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
#' @param Vmax a positive integer to determine the maximum particle's velocity.
#'
#' @param ci a positive integer to determine individual cognitive.
#'
#' @param cg a positive integer to determine group cognitive.
#'
#' @param w a positive integer to determine inertia weight.
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
#' Vmax <- 2
#' ci <- 1.5
#' cg <- 1.5
#' w <- 0.7
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#' 
#' ## calculate the optimum solution using Ant Lion Optimizer 
#' best.variable <- PSO(sphere, optimType="MIN", numVar, numPopulation=20, 
#'                  maxIter=100, rangeVar, Vmax, ci, cg, w)
#' 
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(best.variable)
#' 
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable 
#'         and \code{vn} is value of \code{n-th} variable.
#' 
#' @references
#' Kennedy, J. and Eberhart, R. C. Particle swarm optimization.
#' Proceedings of IEEE International Conference on Neural Networks, Piscataway, NJ. pp. 1942-1948, 1995
#' 
#' Shi, Y. and Eberhart, R. C. A modified particle swarm optimizer. 
#' Proceedings of the IEEE Congress on Evolutionary Computation (CEC 1998), Piscataway, NJ. pp. 69-73, 1998
#' @export

PSO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, Vmax=2, ci=1.49445, cg=1.49445, w=0.729){
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

	# generate initial population of particle
	particles <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# calculate the initial local best
	# local best for this step is a global best
	Gbest <- calcBest(FUN, optimType, particles)
	Lbest <- particles

	# initial velocity of each particle
	velocity <- generateRandom(numPopulation, dimension, -Vmax, Vmax)

	# find the best particle position
	bestParticle <- engine.PSO(FUN, optimType, maxIter, lowerBound, upperBound, Vmax, ci, cg, w, Gbest, Lbest, particles, velocity)

	return(bestParticle)
}

## support function for calculating best position with PSO algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param Vmax maximum velocity
# @param ci individual cognitive
# @param cg group cognitive
# @param w inertia weight
# @param Gbest initial global best
# @param Lbest initial local best
# @param particles population of particle
# @param velocity velocity for particles

engine.PSO <- function(FUN, optimType, maxIter, lowerBound, upperBound, Vmax, ci, cg, w, Gbest, Lbest, particles, velocity){
	FLbest <- calcFitness(FUN, optimType, Lbest)
	FGbest <- optimType*FUN(Gbest)
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
	for (t in 1:maxIter){
		for (i in 1:nrow(particles)){
			for (d in 1:ncol(particles)){
				# pick random rumber
				ri <- runif(1)
				rg <- runif(1)

				# update the particle velocity
				newV <- w * velocity[i,d] + ci*ri*(Lbest[i,d]-particles[i,d]) + cg*rg*(Gbest[d]-particles[i,d])
				
				# check range velocity
				if(newV < -Vmax) newV <- -Vmax
				if(newV > Vmax) newV <- Vmax
				velocity[i,d] <- newV

				newPos <- particles[i,d] + velocity[i,d]
				# check range search space
				if(length(lowerBound)==1){
					if(newPos < lowerBound) newPos <- lowerBound
					if(newPos > upperBound) newPos <- upperBound
				}else{
					if(newPos < lowerBound[d]) newPos <- lowerBound[d]
					if(newPos > upperBound[d]) newPos <- upperBound[d]
				}
				particles[i,d] <- newPos

				# check Local best and Global best
				F <- optimType*FUN(particles[i,])
				if(F < FLbest[i]){
					Lbest[i,] <- particles[i,]
					FLbest[i] <- F
					if(FLbest[i] < FGbest){
						Gbest <- Lbest[i,]
						FGbest <- FLbest[i]
					}
				}
			}
		}
		curve[t] <- FGbest
		setTxtProgressBar(progressbar, t)
	}
	close(progressbar)
	curve <- curve*optimType
	## plot(c(1:maxIter), curve, type="l", main="PSO", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  ## ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(Gbest)
}
