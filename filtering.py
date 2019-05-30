 

Raw <- read.csv("/home/subin95/바탕화면/g/g.txt",header = T, sep = '\t', stringsAsFactors = F)

similar <- read.csv("/home/subin95/바탕화면/g/similar_cut.txt",header = F, sep = '\t', stringsAsFactors = F)

result <- read.csv("/home/subin95/바탕화면/g/result.txt",header = T, sep = '\t', stringsAsFactors = F)

seed1_result <- read.csv("/home/subin95/바탕화면/g/result.txt",header = T, sep = '\t', stringsAsFactors = F)

seed2_result <- read.csv("/home/subin95/바탕화면/g/result.txt",header = T, sep = '\t', stringsAsFactors = F)

 

 

fileConn<-file("output.txt")

writeLines(Raw[i,], fileConn)

 

for (i in 0:3659){

  i = i+1

  j1 = 1

  j2 = 1

  while(Raw[i,1] != similar[j1,2]){

    j1 = j1+1

    if(j1==19684){

      j1=1

      break

    }

  }

  while(Raw[i,2] != similar[j2,2]){

    j2 = j2+1

    if(j2==19684){

      j2=1

      break

    }

  }

  if(similar[j1,1] != similar[j2,1]){

    result <- rbind(result,Raw[i,])

  }

}

 

for(i in 0:2675){

  i = i + 1

  if(result[i,5] == 1){

    if(result[i,6] >= 2){

      seed1_result<-rbind(seed1_result,result[i,])

    }

  }else{

    seed2_result<-rbind(seed2_result,result[i,])

  }

}

 