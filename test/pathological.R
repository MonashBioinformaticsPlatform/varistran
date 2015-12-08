
source("test/common.R")

cat("\n0x0 matrix\n")
varistran::vst( matrix(numeric(0),nrow=0,ncol=0) )

cat("\n10x0 matrix\n")
varistran::vst( matrix(numeric(0),nrow=10,ncol=0) )

cat("\n0x10 matrix\n")
varistran::vst( matrix(numeric(0),nrow=0,ncol=10) )

cat("\n1x10 matrix of 0\n")
varistran::vst( matrix(0,nrow=1,ncol=10) )

cat("\n1x10 matrix of 10\n")
varistran::vst( matrix(10,nrow=1,ncol=10) )