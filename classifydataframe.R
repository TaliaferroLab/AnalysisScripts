myvector <- rep('none',nrow(test))
myvector2 <- rep('no', nrow(test))
myvector[which(test$kmer == 'GAGGG')] = 'HNRNPA2B1'
myvector[which(test$kmer == 'CAGGG')] = 'HNRNPA2B1'
myvector[which(test$kmer == 'AGGGC')] = 'HNRNPA2B1'
myvector[which(test$kmer == 'GGGAG')] = 'SRSF'
myvector[which(test$kmer == 'GGAGG')] = 'SRSF'
myvector[which(test$kmer == 'CAGGA')] = 'SRSF'
#myvector[which(test$kmer == 'UGCAU')] = 'RBFOX'
#myvector[which(test$kmer == 'GCAUG')] = 'RBFOX'
#myvector[which(test$kmer == 'UGCUU')] = 'MBNL'
#myvector[which(test$kmer == 'AUGCU')] = 'MBNL'
#myvector[which(test$kmer == 'UGCUG')] = 'MBNL'
#myvector[which(test$kmer == 'ACUAA')] = 'QKI'
#myvector[which(test$kmer == 'CUAAC')] = 'QKI'
#myvector[which(test$kmer == 'GUAUG')] = 'HNRNPU'
myvector2[which(test$corrected_pvalue < 0.1)] = 'yes'



test2 <- data.frame(test, protein = myvector, sigenriched = myvector2)