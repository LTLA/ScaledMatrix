#' @export
setClass("ScaledMatrixSeed", slots=c(.matrix="ANY", center="numeric", scale="numeric", use_center="logical", use_scale="logical", transposed="logical"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ScaledMatrix", contains="DelayedMatrix", slots=c(seed="ScaledMatrixSeed"))
