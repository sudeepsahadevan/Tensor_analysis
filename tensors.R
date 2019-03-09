# TODO: Add comment
# 
# Author: sahadeva
# do tensor analysis
###############################################################################

#' helper function to return an object of atensor class
#' @param ndim : a vector of (length 3) dimensions, [row,column,number of matrices]
#' @param mode: mode of the tensor (scalar)
#' @param dimnames: row and column names of the matrices (list)
#' @param data: acutal data (array object)
#' @return an object of "atensor" class
.helper.get.tensor <- function(ndim, mode, dimnames, data){
	setClass("atensor", slots = list(ndim = "numeric", mode = "numeric", dimnames = "list", data = "array"))
	tout <- new("atensor", ndim = ndim, mode = mode, dimnames = dimnames, data = data)
	return(tout)
}



#' given a list of matrices, convert the list of matrices into a multidimensional array
#' @param mat.list : list of matrices, nrow and ncol of every matrix MUST be same
#' @param mode: mode for the tensor
#' @return array
get.tensor <- function (mat.list, mode = 3) {
	if (is.null(mode)) {
		stop("ERROR! mode MUST be specified\n")
	}
	if (mode < 3) {
		stop("ERROR! for modes less than 3, look into vector, matrix represenations\n")
	}
	array.out <- array(0, c(dim(mat.list[[1]]), length(mat.list)))
	ndims <- vector("numeric")
	namedim <- list()
	for (i in 1:length(mat.list)) {
		if (i == 1) {
			ndims <- dim(mat.list[[i]])
			namedim <- dimnames(mat.list[[i]])
			if (is.null(namedim)) {
				stop("ERROR! input matrices are expected to have row and column names\n")
			}
			array.out[, , i] <- mat.list[[i]]
		}
		else {
			if (all.equal(ndims, dim(mat.list[[i]]))) {
				temp.mat <- mat.list[[i]]
				array.out[, , i] <- temp.mat[namedim[[1]], namedim[[2]]]
			}
			else {
				stop("ERROR! the row and column dimensions of matrix ", i, " do not match with the rest of the list\n")
			}
		}
	}
	setClass("atensor", slots = list(ndim = "numeric", mode = "numeric", dimnames = "list", data = "array"))
	outt <- new("atensor", ndim = dim(array.out), mode = mode, dimnames = namedim, data = array.out)
	return(outt)
}


#' given a tensor (a multidimensional array), unfold it to mode 1,2 or 3 matrix
#' @param tensor: an object of "atensor" class
#' @param mode: the mode to unfold the tensor,
#'			  : three modes are supported at present, 1 column wise, 2 row wise, 3 tube wise
#' @param verbose: boolean to control the verbosity level
#' @return matrix
tensor.unfold <- function(tensor,mode = 1,verbose = TRUE) {
	tdata <- tensor@data
	tdim <- dim(tensor@data)
	if (verbose) {cat("Tensor dimensions: ", tdim, "\n")}
	switch(mode,
			`1` = {
				if (verbose) {cat("Mode 1 unfolding dimensions: ", tdim[1], ",", tdim[2] * tdim[3], "\n")}
				outm <- tdata[, , 1]
				for (i in 2:tdim[3]) { outm <- cbind(outm, tdata[, , i])}
			}, 
			`2` = {
				if (verbose) {cat("Mode 2 unfolding dimensions: ", tdim[2], ",", tdim[1] * tdim[3], "\n")}
				outm <- t(tdata[, , 1])
				for (i in 2:tdim[3]) { outm <- cbind(outm, t(tdata[, , i]))}
			}, 
			`3` = {
				if (verbose) { cat("Mode 3 unfolding dimensions: ", tdim[3], ",", tdim[1] * tdim[2], "\n")}
				outm <- as.vector(tdata[, , 1])
				for (i in 2:tdim[3]) {outm <- rbind(outm, as.vector(tdata[, , i]))}
				rownames(outm) <- NULL
			}, 
			stop("ERROR! only modes 1,2 and 3 are supported\n"))
	return(outm)
}

#' given a matrix, the mode of the unfolded matrix and dimensions of the tensor,
#' fold the matrix into  a tensor
#' @param mat: input matrix
#' @param mode: mode of the matrix, 1 for column, 2 for row and 3 for tube
#' @param nrow: number of row elements for the tensor
#' @param ncol: number of column elements for the tensor
#' @param ntube: number of tube elements
#' @param namedim: a list of row and column ids 
#' @return object of "atensor" class
tensor.fold <- function(mat,mode = 1,nmode = 3,n.row,n.col,n.tube,namedim = NULL,verbose = TRUE) {
	matdim <- dim(mat)
	if (is.null(namedim)) {namedim = list()}
	if ((matdim[1] * matdim[2]) != (n.row * n.col * n.tube)) {
		stop("ERROR! the dimensions of the given matrix and input tensor dimensions (n.row,n.col,n.tube) are not equal\n")
	}
	array.out <- array(0, c(n.row, n.col, n.tube))
	switch(mode, 
			`1` = {
				if (matdim[1] != n.row) {
					stop("ERROR! nrow(mat) should be equal to given n.row\n")
				} else if (matdim[2] != (n.col * n.tube)) {
					stop("ERROR! ncol(mat) should be equal to n.col*n.tube\n")
				}
				if (verbose) {
					cat("Folding mode 1 matrix...\n")
				}
				colsplit <- split(c(1:ncol(mat)), f = cut(x = c(1:ncol(mat)), breaks = n.tube, labels = FALSE))
				for (i in 1:n.tube) {
					cols <- as.integer(colsplit[[i]])
					if (length(cols) == n.col) {
						array.out[c(1:n.row), c(1:n.col), i] <- mat[c(1:n.row), cols]
					} else {
						stop("ERROR! ", i, " tube folding has morethan ", n.col, " values in column\n")
					}
				}
			}, 
			`2` = {
				if (matdim[1] != n.col) {
					stop("ERROR! ncol(mat) should be equal to given n.col\n")
				} else if (matdim[2] != (n.row * n.tube)) {
					stop("ERROR! ncol(mat) should be equal to n.row*n.tube\n")
				}
				if (verbose) {
					cat("Folding mode 2 matrix...\n")
				}
				colsplit <- split(c(1:ncol(mat)), f = cut(x = c(1:ncol(mat)), breaks = n.tube, labels = FALSE))
				for (i in 1:n.tube) {
					cols <- as.integer(colsplit[[i]])
					if (length(cols) == n.row) {
						array.out[c(1:n.row), c(1:n.col), i] <- t(mat[c(1:nrow(mat)), cols])
					} else {
						stop("ERROR! ", i, " tube folding has morethan ", n.row, " values in row\n")
					}
				}
			}, 
			`3` = {
				if (matdim[1] != n.tube) {
					stop("ERROR! nrow(mat) should be equal to given n.tube\n")
				} else if (matdim[2] != (n.row * n.col)) {
					stop("ERROR! ncol(mat) should be equal to n.row*n.col\n")
				}
				if (verbose) {
					cat("Folding mode 3 matrix...\n")
				}
				for (i in 1:n.tube) {
					array.out[, , i] <- matrix(mat[i, ], nrow = n.row, ncol = n.col, byrow = FALSE)
				}
			}, 
			stop("ERROR! only modes 1,2 and 3 are supported\n"))
	return(.helper.get.tensor(ndim = dim(array.out), mode = nmode, dimnames = namedim, data = array.out))
}


#' given an object of class atensor, array, matrix or data.frame calculate the Frobenius norm
#' defined as sqrt(sum(tensor*tensor)) 
#' @param obj: an object of "atensor" class OR "array" OR "matrix" OR "data.frame"  
#' @return scalar
fnorm.tensor <- function(obj){
	cobj <- class(obj)
	if(cobj=="atensor"){
		return(sqrt(sum(obj@data*obj@data)))
	}else if((cobj=="array")||(cobj=="matrix")||(cobj=="data.frame")){
		return(sqrt(sum(obj*obj)))
	}else{
		stop("ERROR! obj has unknown class\n")
	}
}

#' does the same thing, don't know why I keep these two functions
get.fnorm <- function (obj) {
	cobj <- class(obj)
	if (cobj == "atensor") {
		return(sqrt(sum(obj@data * obj@data)))
	}
	else if ((cobj == "array") || (cobj == "matrix") || (cobj == "data.frame")) {
		return(sqrt(sum(obj * obj)))
	}
	else {
		stop("ERROR! obj has unknown class\n")
	}
}
	
#' tensor higher order singular value decomposition
#' given an object of "atensor" class perform higher order svd
#' @param tensor: an object of "atensor" class
#' @param ranks: rank for decomposing the tensor, either a scalar or a vector
#' @return atensor object
ho.svd <- function(tensor,ranks = NULL,verbose = TRUE) {
	if (is.null(ranks)) {
		ranks <- tensor@ndim
		if (verbose) {
			cat("Info: using tensor modes as ranks\n")
		}
	}
	else if (length(ranks) == 1) {
		if (verbose) {
			cat("Info: only one rank is given, using rank: ", ranks, " for all decompositions\n")
		}
		ranks <- rep(ranks, length(tensor@ndim))
	}
	else if (length(ranks) != length(tensor@ndim)) {
		stop("Error! number of ranks given does not match the tensor mode (they must be equal)\n")
	}
	ulist <- list()
	for (i in 1:length(ranks)) {
		matun <- tensor.unfold(tensor, mode = i)
		ulist[[i]] <- svd(matun, nu = ranks[i])$u
	}
	tensor.svd <- tensor.list.mult(tensor, mat.list = lapply(ulist, t), mode.list = c(1:length(ulist)))
	tensor.est <- tensor.list.mult(tensor.svd, mat.list = ulist,mode.list = c(1:length(ulist)))
	resid <- fnorm.tensor(tensor.est@data - tensor@data)
	return(list(t.svd = tensor.svd, t.est = tensor.est, residual = resid))
}

#' given a list of matrices and a tensor, unroll the tensor, multiply with the list and 
#' roll it back to  a tensor
#' @param tensor: an object of "atensor" class
#' @param mat.list: a list of matrices
#' @param mode.list: an integer vector with length(list.mode) == length(matlist) and the values 
#' 				   : denote the mode in which the tensor must be unrolled
#' @param verbose: A boolean controlling the verbosity level
#' @return object of class "atensor" object
tensor.list.mult <- function (tensor, mat.list, mode.list, verbose = TRUE) {
	tdim <- tensor@ndim
	nmode <- tensor@mode
	if (length(mat.list) != length(mode.list)) {stop("ERROR! length(mat.list) must be equal to length(mode.list)\n")}
	if (length(unique(mode.list)) > length(tdim)) {stop("ERROR! unique elements in mode.list exceeds mode of the tensor\n")}
	modtensor <- tensor
	mats <- list()
	for (i in 1:length(mode.list)) {
		tmpdim <- modtensor@ndim
		unmat <- tensor.unfold(modtensor, mode = mode.list[i], verbose = verbose)
		mmat <- mat.list[[mode.list[i]]]
		tmpdim[mode.list[i]] <- nrow(mmat)
		tmat <- mmat %*% unmat
		tensor.dim <- .helper.dimensions(dim(tmat), tmpdim, mode.list[i])
		modtensor <- tensor.fold(mat = tmat, mode = mode.list[i], nmode = nmode, n.row = tensor.dim[1], n.col = tensor.dim[2], n.tube = tensor.dim[3], namedim = tensor@dimnames, verbose = verbose)
	}
	return(modtensor)
}

#' helper function to get the tensor dimensions during tensor folding 
.helper.dimensions <- function(mat.dim,t.dim,mode){
	switch(mode, 
			`1` = {
				if ((mat.dim[2])%%3 != 0) {
					stop("ERROR! tensor is not rank 3, only upto rank 3 tensors are supported")
				}
				return(c(mat.dim[1], (mat.dim[2]/3), 3))
			}, 
			`2` = {
				if ((mat.dim[2])%%3 != 0) {
					stop("ERROR! tensor is not rank 3, only upto rank 3 tensors are supported")
				}
				return(c((mat.dim[2]/3), mat.dim[1], 3))
			}, 
			`3` = {
				if (mat.dim[2] != (t.dim[1] * t.dim[2])) {
					stop("ERROR! mode 3 (tube folding) ncol(unfolded matrix) is not a product of nrow(tensor)*ncol(tensor)\n")
				}
				return(c(t.dim[1], t.dim[2], 3))
			},
			stop("ERROR! only rank 3 tensors are supported\n"))
}

#' given a tensor do cp decomposition
decompose.cp <- function (tensor, rank = NULL, init =  "svd", tol = 1e-06, maxit = 500, verbose = TRUE) {
	if ((is.null(rank)) || rank == 0) {
		stop("ERROR! rank MUST be provided\n")
	}
	if (length(init) > 1) {
		init <- init[1]
	}
	tmodes <- length(tensor@ndim)
	initmats <- unfoldmats <- list()
	for (i in 1:tmodes) {
		umat <- tensor.unfold(tensor = tensor, i, FALSE)
		unfoldmats[[i]] <- umat
		initmats[[i]] <- .helper.init.cp(umat, rank, init)
	}
	est.tensor <- tensor
	prev.fn <- 0
	cur.fn <- orig.fn <- (fnorm.tensor(tensor))^2
	iter <- 1
	dif.se <- 0
	s.se <- 0
	converged <- true.converge <- FALSE
	cat("start iteration\n")
	while (!converged) {
		for (m in 1:tmodes) {
			temp.initmats <- rev(initmats[-m])
			temp.init <- unfoldmats[[m]] %*% (.helper.khatrirao.list(temp.initmats)) %*% (.helper.pinv.cp(temp.initmats))
			temp.lambda <- sqrt(colSums(temp.init^2))
			initmats[[m]] <- sweep(x = temp.init, MARGIN = 2, temp.lambda, "/")
			sdt <- superdiagonal.tensor(nmode = tmodes, ncomp = rank, diag.el = temp.lambda)
			est.tensor <- tensor.list.mult(sdt, initmats, c(1:tmodes), verbose = FALSE)
		}
		prev.fn <- cur.fn
		cur.fn <- (fnorm.tensor(est.tensor))^2
		if ((abs(cur.fn - prev.fn)/orig.fn) < tol) {
			converged <- true.converge <- TRUE
		}
		else if (iter >= maxit) {
			converged <- max.iteration <- TRUE
			cat("Warning! tensors did not converge, but maximum iterations completed\n")
		}
		if (!verbose) {
			if ((iter%%50 == 0) || (iter == maxit) || (converged)) {
				cat("iteration: ", iter, "\r")
			}
			flush.console()
		}
		else {
			cat("iteration: ", iter, "\n")
		}
		iter <- iter + 1
	}
	cat("\n")
	return(list(U.tensor = initmats, cp.tensor = est.tensor, tol = dif.se, convergence = true.converge, total.iteration = iter, max.iteration = maxit))
}

#' helper to compute psuedo inverse of a matrix
.helper.pinv.cp <- function(in.list){
	library(MASS)
	if (class(in.list[[1]]) != "matrix") {
		stop("ERROR! the elements of the list MUST be matrix\n")
	}
	retmat <- matrix()
	for (i in 1:length(in.list)) {
		if (i == 1) {
			retmat <- crossprod(in.list[[i]])
		}
		else {
			retmat <- retmat * crossprod(in.list[[i]])
		}
	}
	return(ginv(retmat))
}


#' function for cp initializations
.helper.init.cp <- function(mat,r,init ="svd"){
	init <- match.arg(init,c("rand","svd","unif"),several.ok = FALSE)
	switch(init, 
			rand = {
				if (r > ncol(mat)) {
					stop("ERROR! the rank value (no. of components) exceeds the dimensions of the unfolded matrix\n")
				}
				rmat <- matrix(rnorm(n = nrow(mat) * r), nrow = nrow(mat), ncol = r)
				return(rmat)
			},
			svd = {
				if (r > nrow(mat)) {
					message("WARNING! the rank value (no. of components) exceeds the no. of left singular values of the unfolded matrix, returning the maximum possible")
					r <- nrow(mat)
				}
				svdmat <- svd(x = mat, nu = r)$u
				return(svdmat)
			},
			unif = {
				if (r > ncol(mat)) {
					stop("ERROR! the rank value (no. of components) exceeds the dimensions of the unfolded matrix\n")
				}
				rmat <- matrix(runif(n = nrow(mat) * r), nrow = nrow(mat), ncol = r)
				return(rmat)
			}, 
			stop("ERROR! unknown init option, MUST be one of rand, svd or unif\n"))
}

#' khatri rao product helper
.helper.khatrirao.list <- function(in.list){
	library(Matrix)
	if (class(in.list[[1]]) != "matrix") {
		stop("ERROR! the elements of the list MUST be matrix\n")
	}
	retmat <- matrix()
	for (i in 1:length(in.list)) {
		if (i == 1) {
			tempmat <- in.list[[i]]
			retmat <- as.matrix(KhatriRao(matrix(rep(1, ncol(tempmat)), nrow = 1, ncol = ncol(tempmat)), tempmat))
		}
		else {
			retmat <- as.matrix(KhatriRao(retmat, in.list[[i]]))
		}
	}
	return(retmat)
}

#' return super diagonal tensor
superdiagonal.tensor <- function(nmode, ncomp, diag.el = 1){
	if(length(diag.el) == 1){
		cat("Info: only one diag.el found, repeating the same value for all diagonal elements of the tensor\n")
		diag.el <- rep(diag.el, ncomp)
	}
	else if (length(diag.el) != ncomp) {stop("ERROR! ncomp value: ", ncomp, " doesnot match the number of diag.el values, ", length(diag.el), "\n")}
	ndim <- rep(ncomp, nmode)
	sdt <- array(0, ndim)
	for(i in 1:length(diag.el)){
		sdt[i, i, i] <- diag.el[i]
	}
	return(.helper.get.tensor(ndim = ndim, mode = nmode, dimnames = list(), data = sdt))
}

#' compute right singular matrix for higher order svd
#' Given a list of input matrices, compute the common right singular eigen vectors according to publication
#' Ponnapalli SP, Saunders MA, Van Loan CF, Alter O (2011) 
#' A Higher-Order Generalized Singular Value Decomposition for Comparison of Global mRNA Expression from Multiple Organisms. 
#' PLoS ONE 6(12): e28072
#' NOTE: Moreover, W is not guaranteed to be non-defective and have a full set of real eigenvalues and eigenvectors. However, even in the absence of an exact common decomposition, the real part of the complex eigenvectors can be used to derive a low rank approximation of the common subspace and extract common and differential covariance structures from the data.
#' @param matlist: a list of matrices
#' @param pinv: Boolean, if TRUE, calculates Moore–Penrose pseudoinverse of the matrix, else calculates the general inverse of a matrix
#' @return a list
compute.right.singular <- function(matlist, pinv = TRUE){
	library(MASS)
	smat <- matrix()
	first <- TRUE
	invfn <- function(bool, mat) {
		if (bool){ return(ginv(mat)) }
		else { return(solve(mat)) }
	}
	for (i in 1:length(matlist)) {
		mat1 <- matlist[[i]]
		mat1inv <- invfn(pinv, mat1)
		j <- i + 1
		while (j <= length(matlist)) {
			mat2 <- matlist[[j]]
			mat2inv <- invfn(pinv, mat2)
			if (first) {
				smat <- (mat1 %*% mat2inv) + (mat2 %*% mat1inv)
				first <- FALSE
			}
			else {
				smat <- smat + ((mat1 %*% mat2inv) + (mat2 %*% mat1inv))
			}
			j <- j + 1
		}
	}
	n <- length(matlist)
	smat <- smat/(n * (n - 1))
#	return(smat)
	return(eigen(x = smat, symmetric = isSymmetric(smat)))
}

#' compute higher order singular value decomposition
#' Given a list of matrices, compute higher order generalized singular value decomposition as defined in the publication:
#' Ponnapalli SP, Saunders MA, Van Loan CF, Alter O (2011) 
#' A Higher-Order Generalized Singular Value Decomposition for Comparison of Global mRNA Expression from Multiple Organisms. 
#' PLoS ONE 6(12): e28072
#' NOTE: "Moreover, W is not guaranteed to be non-defective and have a full set of real eigenvalues and eigenvectors. However, even in the absence of an exact common decomposition, the real part of the complex eigenvectors can be used to derive a low rank approximation of the common subspace and extract common and differential covariance structures from the data."
#' @param matlist: a list of matrices
#' @param crossp: Boolean, if TRUE, calculates matrix crossproduct t(x)%*%x
#' @param pinv: Boolean, if TRUE, calculates Moore–Penrose pseudoinverse of the matrix, else calculates the general inverse of a matrix
#' @return An object of "hogsvd" class 
compute.hogsvd <- function(matlist, crossp = FALSE, pinv = TRUE){
	matlist <- lapply(matlist, as.matrix)
	tmplist <- matlist
	if (crossp) {
		tmplist <- lapply(matlist, crossprod)
	}
	eigen.smat <- compute.right.singular(tmplist, pinv)
	ulist <- list()
	dlist <- list()
	for (i in 1:length(matlist)) {
		dmat <- as.matrix(t(matlist[[i]]))
		bmat <- t(solve(eigen.smat$vectors, dmat))
		bdiag <- apply(bmat, 2, function(x) sqrt(sum(x^2)))
		bu <- bmat/bdiag
		matnam <- names(matlist[i])
		ulist[[matnam]] <- bu
		dlist[[matnam]] <- bdiag
	}
	return(.helper.get.hogsvd(ulist, dlist, eigen.smat$vectors, eigen.smat$values))
}

#' Function to return an object of hogsvd class
#' given a list of left singular matrices, 
#' a  list of diagonal values, 
#' a matrix of right singular vectors
#' and a vector of right singular diagonal values, 
#' return an object of "hogsvd" class
.helper.get.hogsvd <- function (ulist, dlist, v, vd){
	setClass("hogsvd", slots = list(Ulist = "list", Dlist = "list", V = "matrix", Vd = "vector"))
	tout <- new("hogsvd", Ulist = ulist, Dlist = dlist, V = v, Vd = vd)
	return(tout)
}

#' given a matrix and a number n, split the matrix column wise into n submatrices and return as a list
#' @param inmat: input matrix
#' @param nmat: split inmat into nmat sub matrices
#' @return list
matrix.to.list <- function (inmat, nmat) {
	if (ncol(inmat)%%nmat != 0) {
		stop("ERROR! ncol(matrix) is not a multiple of nmat, cannot split the matrix evenly\n")
	}
	colsplit <- split(c(1:ncol(inmat)), f = cut(c(1:ncol(inmat)), 
					breaks = nmat, labels = FALSE))
	matlist <- list()
	for (i in 1:length(colsplit)) {
		tmpmat <- inmat[, colsplit[[i]]]
		nam <- paste(unique(gsub(colnames(tmpmat), pattern = "\\d*(\\_|\\.)*\\d{1,}$", 
								replacement = "")), collapse = "_")
		matlist[[nam]] <- tmpmat
	}
	return(matlist)
}

