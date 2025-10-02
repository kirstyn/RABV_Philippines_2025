#  install.packages("redcapAPI")

library(redcapAPI)

unlockREDCap(c(rcon= 'rage-redcap'),  
             keyring= "login",
             envir= globalenv(),
             url= 'https://cvr-redcap.mvls.gla.ac.uk/redcap/redcap_v15.5.13/API/')
