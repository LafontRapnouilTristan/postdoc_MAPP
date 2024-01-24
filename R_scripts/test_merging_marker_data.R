dummy_16s <- data.frame(asv_id=sprintf("asv_%03d", 1:15),
                        seq=replicate(15,paste(sample(LETTERS[c(1,3,7,20)],50,replace=T),collapse = '')),
                        taxo_1=replicate(15,paste0("k_",sample(letters[1:3],1,replace = T),";o_",sample(letters[1:3],1,replace = T),";c_",sample(letters[1:3],1,replace = T))),
                        taxo_2=replicate(15,paste0("k_",sample(letters[1:3],1,replace = T),";o_",sample(letters[1:3],1,replace = T),";c_",sample(letters[1:3],1,replace = T)))
)
