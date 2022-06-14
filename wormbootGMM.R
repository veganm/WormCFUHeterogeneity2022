# Parameters for Salmonella enterica three days data
# (mean1, var1)=(2.985569, 0.4568056) and (mean2, var2)=(4.742767, 0.1239928), containing 58% and 42% of the mass respectively

#reps<-100
#fUP<-0.58 #fraction of the total population that is "highly colonized"
#fUPb<-fUP
#meanD<-2.986
#meanU<-4.743
#mySDD<-sqrt(0.457)
#mySDU<-sqrt(0.124)

wormbootGMM<-function(reps, batches, meanD, varD, meanU, varU, fUP, fDIFF){
	# Script for generating log-scale fake data based on a Gaussian mixture model
	# Prints p-values for t-tests and Mann-Whitney U tests in order of batch size
	# Returns an example data frame of simulated data
	# where "batch" is the number of individual measurements in the batch digest
	# "Count" is CFU/worm, and "logCount" is the original log-scale simulated data
	#"fUP" should be the larger of the fractions [0,1] of high-mode individuals
	# and "fDIFF" is the difference in modal fraction across the two populations
	fUPb<-fUP-fDIFF
	mySDD<-sqrt(varD)
	mySDU<-sqrt(varU)
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
	  up<-sum(rbinom(batches,1, fUP))
	  down=batches-up
	  tempA<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	  tempA[tempA<1.3]<-0 
	  up<-sum(rbinom(batches,1, fUPb))
	  down=batches-up
	  tempB<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	  tempB[tempB<1.3]<-0 
	# batch 5
	  temp5A<-numeric(batches)
	  temp5B<-numeric(batches)
	  temp10A<-numeric(batches)
	  temp10B<-numeric(batches)
	  temp20A<-numeric(batches)
	  temp20B<-numeric(batches)
	  temp50A<-numeric(batches)
	  temp50B<-numeric(batches)

	  for(k in 1:batches){
	    up<-sum(rbinom(5, 1, fUP))
	    down=5-up
	    temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    up<-sum(rbinom(5, 1, fUPb))
	    down=5-up
	    temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    temp1[temp1<1.3]<-0
	    temp2[temp2<1.3]<-0
	    temp5A[k]<-log10(mean(10^temp1))
	    temp5B[k]<-log10(mean(10^temp2))
	#now batch 10
	    up<-sum(rbinom(10, 1, fUP))
	    down=10-up
	    temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    up<-sum(rbinom(10, 1, fUPb))
	    down=10-up
	    temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    temp1[temp1<1.3]<-0
	    temp2[temp2<1.3]<-0
	    temp10A[k]<-log10(mean(10^temp1))
	    temp10B[k]<-log10(mean(10^temp2))
	# batch 20
	    up<-sum(rbinom(20, 1, fUP))
	    down=20-up
	    temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    up<-sum(rbinom(20, 1, fUPb))
	    down=20-up
	    temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    temp1[temp1<1.3]<-0
	    temp2[temp2<1.3]<-0
	    temp20A[k]<-log10(mean(10^temp1))
	    temp20B[k]<-log10(mean(10^temp2))
	#batch50
	    up<-sum(rbinom(50, 1, fUP))
	    down=50-up
	    temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    up<-sum(rbinom(50, 1, fUPb))
	    down=50-up
	    temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
	    temp1[temp1<1.3]<-0
	    temp2[temp2<1.3]<-0
	    temp50A[k]<-log10(mean(10^temp1))
	    temp50B[k]<-log10(mean(10^temp2))
 	 }
    
  #t tests on log data
  test1<-t.test(tempA, tempB)
  test5<-t.test(temp5A, temp5B)
  test10<-t.test(temp10A, temp10B)
  test20<-t.test(temp20A, temp20B)
  test50<-t.test(temp50A, temp50B)

  #Wilcoxon rank sum tests
  wtest1<-wilcox.test(tempA,tempB)
  wtest5<-wilcox.test(temp5A,temp5B)
  wtest10<-wilcox.test(temp10A,temp10B)
  wtest20<-wilcox.test(temp20A,temp20B)
  wtest50<-wilcox.test(temp50A,temp50B)
    
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

pt<-c(sum(pv1<0.05)/reps, sum(pv5<0.05)/reps, sum(pv10<0.05)/reps, sum(pv20<0.05)/reps, sum(pv50<0.05)/reps)
pw<-c(sum(wpv1<0.05)/reps, sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, sum(wpv20<0.05)/reps, sum(wpv50<0.05)/reps)
print(pt)
print(pw)
batch<-c(rep(1,batches), rep(5, batches), rep(10, batches), rep(20, batches), rep(50, batches))
dataA<-c(tempA, temp5A, temp10A, temp20A, temp50A)
dataB<-c(tempB, temp5B, temp10B, temp20B, temp50B)
frameA<-data.frame(batch, logCount=dataA)
frameA$set<-"A"
frameB<-data.frame(batch, logCount=dataB)
frameB$set<-"B"
mydata<-rbind(frameA, frameB)
mydata$Count<-10^mydata$logCount
return(mydata)
}