# TODO: Generates 34 testfunctions, one for each testcase taken from the Roche Diagnostics method comparison
#       intranet module testcase collection. Results obtained from the approximative Passing-Bablok implementation
#       (PaBaLarge) are compared to results of the exact implementation (PaBa).
# 
# Author: Andre Schuetzenmeister
###############################################################################



# one could increase the number of bins (NBins) to minimize the differences, here we use the default setting of NBins=1e06

PaBaLargePrecision <- 1e-04                                                       

cat("\n\n********************************************\nmcPaBaLarge.R method comparison test cases\n********************************************")

load(".\\TestCaseCollection\\testcases.RData")

TCnames <- names(testcases)

# 

genericPaBaLargeTest <- function(Data, Name, Exception=FALSE)
{
    X <- Data[,1]
    Y <- Data[,2]
    
    if(Exception)
    {
        cat("\n\nTestcase:", Name, "\n")
        NData <- nrow(Data)
        cat("\nN(Data)         =", NData)
        NUData <- nrow(unique(Data))
        cat("\nN(unique(Data)) =", NUData)
        cat("\nN(Ties)         =", NData-NUData)
        cat("\n\n")
        checkException(mcreg(X, Y, method.reg="PaBaLarge", method.ci="analytical",NBins=1e06)@para)
        checkException(mcreg(X, Y, method.reg="PaBa",      method.ci="analytical")@para)
    }
    else
    {
        cat("\n\nTestcase: ", Name, "\n")
        NData <- nrow(Data)
        cat("\nN(Data)         =", NData)
        NUData <- nrow(unique(Data))
        cat("\nN(unique(Data)) =", NUData)
        cat("\nN(Ties)         =", NData-NUData)
        cat("\n\n")
        resPaBaL <- mcreg(X, Y, method.reg="PaBaLarge", method.ci="analytical", NBins=1e06)@para
        resExact <- mcreg(X, Y, method.reg="PaBa",      method.ci="analytical")@para
        
        checkEquals(resPaBaL, resExact, tolerance=PaBaLargePrecision)
    }    
}

# imitate call-by-value argument passing

cloneLocalArgs <- function(Data, Name, Exception)
{
    cloneData <- Data                           # generate local copies
    cloneName <- Name
    cloneExpt <- Exception
    locFunc <- function(){genericPaBaLargeTest(cloneData, cloneName, cloneExpt)}
    return(locFunc) 
}

# generate test-function for each dataset of the testcase collection

for( i in 1:length(testcases))
{
    Fname <- paste("test.PaBaLargeTestcase_", TCnames[i], sep="")

    assign( Fname, cloneLocalArgs(testcases[[i]], TCnames[i], TCnames[i] %in% c("part_1_dataset_2", "part_1_dataset_12")) )
}


test.PaBaLarge.angle1 <- function()
{
	data(creatinine)
	crea <- na.omit(creatinine)
	
	fit.PaBa.radian  		<- mcreg(crea[,1], crea[,2], method.reg="PaBa", 	 slope.measure="radian",  rng.seed=331, method.ci="analytical")
	fit.PaBa.tangent 		<- mcreg(crea[,1], crea[,2], method.reg="PaBa", 	 slope.measure="tangent", rng.seed=420, method.ci="analytical")
	                                                                                                      
	fit.PaBaLarge.radian  	<- mcreg(crea[,1], crea[,2], method.reg="PaBaLarge", slope.measure="radian",  rng.seed=331, NBins=1e8, method.ci="analytical")
	fit.PaBaLarge.tangent 	<- mcreg(crea[,1], crea[,2], method.reg="PaBaLarge", slope.measure="tangent", rng.seed=420, NBins=1e8, method.ci="analytical")
	
	checkEquals(fit.PaBa.radian@para[,"EST"],  fit.PaBaLarge.radian@para[,"EST"],  tol=1e-7)
	checkEquals(fit.PaBa.tangent@para[,"EST"], fit.PaBaLarge.tangent@para[,"EST"], tol=1e-7)
	
	# pathological example from the mcreg-Rdoc
	
	x1 <- 1:10
	y1 <- 0.5*x1
	x <- c(x1,y1)
	y <- c(y1,x1)
	
	m1 <- mcreg(x,y,method.reg="PaBa",method.ci="analytical",slope.measure="radian")
	m2 <- mcreg(x,y,method.reg="PaBa",method.ci="analytical",slope.measure="tangent")
	
	m1.2 <- mcreg(x,y,method.reg="PaBaLarge",method.ci="analytical",slope.measure="radian",  NBins=1e8)
	m2.2 <- mcreg(x,y,method.reg="PaBaLarge",method.ci="analytical",slope.measure="tangent", NBins=1e8)
	
	checkEquals(m1@para[,"EST"], m1.2@para[,"EST"], tol=1e-7)
	checkEquals(m2@para[,"EST"], m2.2@para[,"EST"], tol=1e-7)
}

