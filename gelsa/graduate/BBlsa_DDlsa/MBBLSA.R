BlockResample<-function(x){
  EstimatedCoefficient<-arima(x,order=c(1,0,0),include.mean=TRUE,
            method='ML',optim.control =list(maxit = 1000))$coef[[1]]
  lopt<-max(ceiling((6^(1/2) * abs(EstimatedCoefficient) / (1 - EstimatedCoefficient^2))^(2/3) * length(x)^(1/3) - 0.5),1)
  block_number<-length(x) %/% lopt+1
  block_index<-sample.int(max(length(x)-lopt+1,1),size = block_number,replace = T)
  perm_x<-rep(0,length(x))
  if (block_number == 1){
    middle <- sample.int(length(x), block_number)
    return(c(x[middle:length(x)],x[1:(middle-1)]))
  }
  else{
    for (i in 1:block_number){
      perm_x[(1+lopt*(i-1)):(lopt*i)] <- x[block_index[i]:(block_index[i]+lopt-1)]
    }
    return(perm_x[1:length(x)])
  }
}

MovingBlockBootstrap<-function(x, y, numPermu=1000, maxDelay=3){
  x <- rankNormalization(x)
  y <- rankNormalization(y)
  originalSimilarity<- LocalSimilarity(x, y, maxDelay)
  localscore <- originalSimilarity[1]
  scoreArray<-rep(0,numPermu)
  for (idx in 1:numPermu){
    scoreArray[idx]<-LocalSimilarity(BlockResample(x), y, maxDelay)[1]
  }
  pValue <- sum(scoreArray>localscore) / numPermu
  value <- t(c(localscore, pValue, originalSimilarity[2:6]))
  colnames(value)<-c('Maxscore','p-value', 'startX', 'startY', 'delay', 'length', 'PosOrNeg')
  return(value)
}
