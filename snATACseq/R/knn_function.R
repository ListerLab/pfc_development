myknn = function(ref_embed, ref_class, test_embed, test_class, k){
    
    output = list()
    kest =knn(train = ref_embed, test = test_embed, cl = 
                  as.factor(ref_class), prob=TRUE, k =k)
    result = character(nrow(test_embed))
    
    for(i in 1:nrow(test_embed)){
        if(kest[i] == as.factor(test_class)[i]) {    
            result[i] = "correct"
        }
        else {
            result[i] = "incorrect"
        }
    }
    
    output$result = result
    correct = subset(test_class, result == "correct")
    incorrect = subset(test_class, result == "incorrect")
    output$correct = correct
    output$incorrect = incorrect
    output$accuracy = length(correct)/length(result)
    output$corprop = summary(correct)/length(correct)
    output$est = kest
    output$percents = summary(correct)/summary(test_class)
    output$celltypes = list()
    output$celltypes$Astro = summary(subset(output$est, test_class == 'Astro'))
    output$celltypes$CGE_dev = summary(subset(output$est, test_class == 'CGE_dev'))
    output$celltypes$ID2 = summary(subset(output$est, test_class == 'ID2'))
    output$celltypes$L2_3 = summary(subset(output$est, test_class == 'L2/3_CUX2'))
    output$celltypes$L4 = summary(subset(output$est, test_class == 'L4_RORB'))
    output$celltypes$L5_6 = summary(subset(output$est, test_class == 'L5/6_THEMIS_TLE4'))
    output$celltypes$MGE_dev = summary(subset(output$est, test_class == 'MGE_dev'))
    output$celltypes$Micro = summary(subset(output$est, test_class == 'Micro'))
    output$celltypes$Oligo = summary(subset(output$est, test_class == 'Oligo'))
    output$celltypes$OPC = summary(subset(output$est, test_class == 'OPC'))
    output$celltypes$PN_dev = summary(subset(output$est, test_class == 'PN_dev'))
    output$celltypes$Poor_Quality = summary(subset(output$est, test_class == 'Poor-Quality'))
    output$celltypes$PV = summary(subset(output$est, test_class == 'PV'))
    output$celltypes$SST = summary(subset(output$est, test_class == 'SST'))
    output$celltypes$Vas = summary(subset(output$est, test_class == 'Vas'))
    output$celltypes$VIP = summary(subset(output$est, test_class == 'VIP'))
    return(output)
}