BatchEntropy <- function(dataset, batch0, L=100, M=100, k=100) {
  #entropy of batch mixing
  # L is the number bootstrapping times
  # M is the number of randomly picked biosamples
  # k is the number of nearest neighbours of biosamples (from all batches) to check

  nbatches<-length(unique(batch0))

  entropy<-matrix(0,L,1)
  set.seed(0)
  for (boot in 1:L) {
    bootsamples<-sample(1:nrow(dataset),M)
    W21<-nn2(dataset,query=dataset[bootsamples,],k)

    for (i in 1:length(bootsamples)){

      for (j in 1:nbatches) {
        xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
        entropy[boot]<-entropy[boot]+xi*log(xi)
      }
    }
  }

  return( (-1)*entropy/length(bootsamples) )
}
