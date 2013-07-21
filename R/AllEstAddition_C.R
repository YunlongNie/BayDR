## compare to IffDiscardPoints.R this part gg function changed  using Bernoulli not binomial, which will match paul's results in Augest3rd email. 


NewCom_C <- function(comb1,Head,AddBin=10,IntMean,IntCon){
  
  temp <- comb1[with(comb1, order(-w)), ]
  
  AddBox <- as.matrix(head(temp, Head))
  LeftBox <- temp[-c(1:Head),]
  AddBox[,3] = rep(1/nrow(temp),nrow(AddBox))
  LeftBox[,3] = rep(1/nrow(temp),nrow(LeftBox))
  
  add <- as.data.frame(do.call(rbind,
                               NewComb_C(Add=AddBox[,-3],weight=AddBox[,3],Addbin=AddBin,
                                         intMean=IntMean,intCon=IntCon)))[,3:5];
  names(add)=names(LeftBox)
  Com1 <- rbind(add,LeftBox) 
  return(Com1)
}





 

 
