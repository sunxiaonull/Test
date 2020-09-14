rm(list=ls())
testa=c(1,2,3,4,5,6,7,8,9,10)
testb=c(1,2,3,4,6,5,7,8,9,10,11,12)

#返回两个变量（向量）中相同的元素个数
testTwoArray<-function(a,b) {
  c=0
  num=min(length(a),length(b))
  for(i in 1:num){
    if(a[i]==b[i]){
      c=c+1
    }
    i=i+1
  }
  return (c)
}


result=testTwoArray(testa,testb)


testc=1:200
dim(testc)=c(10,20)
testc=testc[-c((nrow(testc)-3):nrow(testc)),]
