wormSimBatchBeta<-function(a, b, reps, batches, maxCFU){
  #a and b are parameters of the beta distribution (double)
  #reps is the number of runs in the simulation (integer)
  #batches is the number of data points to generate for each run (integer)
  #maxCFU is the maximum value you want for the CFU/worm synthetic data 
	pv1<-numeric(reps)
	pv5<-numeric(reps)
	pv10<-numeric(reps)
	pv20<-numeric(reps)
	pv50<-numeric(reps)

	wpv1<-numeric(reps)
	wpv5<-numeric(reps)
	wpv10<-numeric(reps)
	wpv20<-numeric(reps)
	wpv50<-numeric(reps)

	for (j in 1:reps){
 	  tempA<-rbeta(batches, a, b)*maxCFU
 	  tempB<-rbeta(batches, a, b)*maxCFU
 	  temp5A<-numeric(batches)
 	  temp5B<-numeric(batches)
	  temp10A<-numeric(batches)
	  temp10B<-numeric(batches)
	  temp20A<-numeric(batches)
	  temp20B<-numeric(batches)
	  temp50A<-numeric(batches)
	  temp50B<-numeric(batches)

	  for(k in 1:batches){
		temp1<-rbeta(5, a, b)*maxCFU
	    temp5A[k]<-log10(mean(temp1))
	    temp2<-rbeta(5, a, b)*maxCFU
	    temp5B[k]<-log10(mean(temp2))
		temp1<-rbeta(10, a, b)*maxCFU
	    temp10A[k]<-log10(mean(temp1))
	    temp2<-rbeta(10, a, b)*maxCFU
	    temp10B[k]<-log10(mean(temp2))
		temp1<-rbeta(20, a, b)*maxCFU
	    temp20A[k]<-log10(mean(temp1))
	    temp2<-rbeta(20, a, b)*maxCFU
	    temp20B[k]<-log10(mean(temp2))
		temp1<-rbeta(50, a, b)*maxCFU
	    temp50A[k]<-log10(mean(temp1))
	    temp2<-rbeta(50, a, b)*maxCFU
	    temp50B[k]<-log10(mean(temp2))
  		}  
    
	  #t tests on log data
	  test1<-t.test(tempA, tempB)
	  test5<-t.test(temp5A, temp5B)
	  test10<-t.test(temp10A, temp10B)
	  test20<-t.test(temp20A, temp20B)
	  test50<-t.test(temp50A, temp50B)

	  #Mann-Whitney U tests
	  wtest1<-wilcox.test(tempA,tempB)
	  wtest5<-wilcox.test(temp5A, temp5B)
	  wtest10<-wilcox.test(temp10A, temp10B)
	  wtest20<-wilcox.test(temp20A, temp20B)
	  wtest50<-wilcox.test(temp50A, temp50B)
  
 	   pv1[j]<-test1$p.value
	  pv5[j]<-test5$p.value
	  pv10[j]<-test10$p.value
	  pv20[j]<-test20$p.value
	  pv50[j]<-test50$p.value

	  wpv1[j]<-wtest1$p.value
	  wpv5[j]<-wtest5$p.value
	  wpv10[j]<-wtest10$p.value
	  wpv20[j]<-wtest20$p.value
	  wpv50[j]<-wtest50$p.value
	}
	t.pvals<-c(sum(pv1<0.05)/reps, sum(pv5<0.05)/reps, sum(pv10<0.05)/reps, sum(pv20<0.05)/reps, sum(pv50<0.05)/reps)
	w.pvals<-c(sum(wpv1<0.05)/reps, sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, sum(wpv20<0.05)/reps, sum(wpv50<0.05)/reps)
	print(t.pvals)
	print(w.pvals)
	batch<-c(rep(1,batches), rep(5, batches), rep(10, batches), rep(20, batches), rep(50, batches))
	dataA<-c(tempA, temp5A, temp10A, temp20A, temp50A)
	dataB<-c(tempB, temp5B, temp10B, temp20B, temp50B)
	frameA<-data.frame(batch, data=dataA)
	frameA$set<-"A"
	frameB<-data.frame(batch, data=dataB)
	frameB$set<-"B"
	mydata<-rbind(frameA, frameB)
	return(mydata)
}