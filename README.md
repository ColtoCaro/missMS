# missMS
Missing Data Models for Label-Free Mass Spectrometry Proteomics

The code in this package was used to fit all of the models found in the paper 
"The effects of nonignorable missing data on label-free mass spectrometry 
proteomics experiments", which can be found at the following URL, 
https://projecteuclid.org/euclid.aoas/1542078037

The selection model can be fit by calling the smp() function.  This function 
takes a dataframe with a header similar to the one found in sampleDat.rda
Further examples of a similar header structure can be found in the 
compositionalMS package, https://github.com/ColtoCaro/compositionalMS

While the code in this repository works, no software engineering has 
been done to make the use of these models simple for the casual user and
we have no plans to work further on the development of tools for LFQ 
proteomics.  That said, if anyone has a serious interest in running the smp() 
function and would like help getting the code working, I'm more than happy to help. 
Feel free to send me a message.  
