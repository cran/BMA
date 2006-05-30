





iBMA.glm<-function(x, ...)
UseMethod("iBMA.glm")

iBMA.glm.data.frame<- function(x, Y, wt = rep(1, nrow(X)), thresProbne0 = 5, glm.family, maxNvar = 30, nIter=100, verbose = FALSE, sorted = FALSE, factor.type = TRUE, ...) 
{

   printCGen<- function(printYN)
   {
        printYN<- printYN
        return(function(x) if (printYN) cat(paste(paste(x,sep="", collapse = " "),"\n", sep="")))
   }
   
        
   
  sortX<- function(Y,X, glm.family, wt)
   {
        fitvec<- rep(NA, times = ncol(X))

        nastyHack_____glm.family <<- glm.family
        nastyHack_____x.df <<- data.frame(X)
        
        glm.out <- glm(Y ~ 1 , family = nastyHack_____glm.family, weights = wt, data = nastyHack_____x.df) 
        scp<- formula(paste("~", paste(colnames(X), sep = "", collapse = " + ")))
        
        addglm <- add1(glm.out, scope = scp , test = "Chisq", data = nastyHack_____x.df)
        fitvec <- addglm$"Pr(Chi)"[-1]
        initial.order<- order(fitvec,decreasing = FALSE)
        sortedX<- X[, initial.order]
        
        return(list(sortedX = sortedX, initial.order = initial.order))
   }
     X<-x
 
    cl <- match.call()
    printC<- printCGen(verbose)
   
   # if factor.type = FALSE then individual factor levels can be dropped, and we
   # need to convert X into its design matrix before we do any iterations
   if (factor.type == FALSE)
   {
     x.df <- data.frame( X)
     X <- model.matrix(terms.formula(~., data = x.df), data = x.df)[,-1]
   }
   
   ### sort variables prior to running iterations if required
   
   if (!sorted)
   {
        printC("sorting X")
        sorted<- sortX(Y,X, glm.family, wt = wt)
        sortedX<- sorted$sortedX
        initial.order<- sorted$initial.order
   }
   else 
   {
        sortedX<- X
        initial.order<- 1:ncol(sortedX)
   }
   
   
   
   
   
   #### Iteration Initiation

   nVar<- ncol(sortedX)
   # make sure that we do not try to use more columns at a time than are present
   maxNvar <- min (maxNvar, nVar)

   stopVar <- 0
   nextVar <- maxNvar + 1
   current.probne0<- rep(0, maxNvar)
   maxProbne0<- rep(0, times = nVar)
   nTimes<- rep(0, times = nVar)
   currIter <- 0
   first.in.model<- rep(NA, times = nVar)
   new.vars<- 1:maxNvar
   first.in.model[new.vars]<- currIter + 1
   iter.dropped<- rep(NA, times = nVar)
   currentSet<- NULL

   current_state<- list(Y = Y,
                        sortedX = sortedX,
                        wt = wt,
                        call = cl,
                        initial.order = initial.order,
                        thresProbne0 = thresProbne0,
                        maxNvar = maxNvar,
                        glm.family = glm.family,
                        nIter = nIter,
                        verbose = verbose,
                        nVar = nVar,
                        currentSet = currentSet,
                        new.vars= new.vars,
                        stopVar = stopVar,
                        nextVar = nextVar,
                        current.probne0 = current.probne0,
                        maxProbne0 = maxProbne0,
                        nTimes = nTimes,
                        currIter = currIter,
                        first.in.model = first.in.model,
                        iter.dropped = iter.dropped) 
   class(current_state)<- "iBMA.intermediate.glm"
   result<- iBMA.glm.iBMA.intermediate.glm(current_state, ...)
   result

 
}







### this function does a set number of iterations of iBMA, returning an intermediate result unless it is finished,
### in which case it returns a final result
        
iBMA.glm.iBMA.intermediate.glm<- function (x, nIter = NULL, verbose = NULL, ...) 
{

   printCGen<- function(printYN)
   {
        printYN<- printYN
        return(function(x) if (printYN) cat(paste(paste(x,sep="", collapse = " "),"\n", sep="")))
   }
   
   
   cs<- x
      
   # check if nIter has been redefined
   if (!is.null(nIter)) cs$nIter<- nIter
   if (!is.null(verbose)) cs$verbose<- verbose
   printC<- printCGen(cs$verbose)
   
   
   finalIter<- cs$currIter + cs$nIter
   
### iterate until a final result is produced (cs$stopVar == 1) or nIter more iterations have been done
   while (cs$stopVar == 0 && cs$currIter < finalIter) 
   {
   
        # add in the new variables
        nextSet<- c(cs$currentSet, cs$new.vars)
        cs$currIter<- cs$currIter + 1    
        
        printC(paste("\n\n starting iteration ",cs$currIter,"   nextVar =",cs$nextVar))
        printC("applying bic.glm now")
    
        currentX<- cs$sortedX[,nextSet]
        colnames(currentX)<- colnames(cs$sortedX)[nextSet]
        ret.bic.glm <- bic.glm (x = currentX, y = cs$Y, glm.family= cs$glm.family, maxCol = cs$maxNvar + 1, ...)
        
        printC(ret.bic.glm$probne0)

        cs$maxProbne0[nextSet]<- pmax(ret.bic.glm$probne0, cs$maxProbne0[nextSet])
        cs$nTimes[nextSet]<- cs$nTimes[nextSet] + 1
        cs$rmVector <- ret.bic.glm$probne0 < cs$thresProbne0

        # adaptive threshold
        if (any(cs$rmVector) == FALSE) 
        {
            # no var to swap in!!, increase threshold
            currMin <- min (ret.bic.glm$probne0)
            printC (paste("no var to swap! Min probne0 = ", currMin, sep=""))
            newThresProbne0 <- currMin + 1
            printC(paste("new probne0 threshold = ", newThresProbne0, sep=""))
            cs$rmVector <- ret.bic.glm$probne0 < newThresProbne0
            
            # check that we do not drop everything!
            if (all(cs$rmVector))
                cs$rmVector<- c(rep(FALSE, times = length(cs$rmVector)-1), TRUE)
        }

        # drop the bad ones...
        cs$iter.dropped[nextSet[cs$rmVector]]<- cs$currIter
        cs$currentSet<- nextSet[!cs$rmVector]

        # now if there are more variables to examine add the new set of variables to the current set
        
        if ( cs$nextVar <= cs$nVar) 
        {
  
            # set up new X
            printC ("generating next set of variables")
            lastVar<- sum(cs$rmVector) + cs$nextVar - 1
            
            # add in bulk if we are not close to the end of the variables,
            if (lastVar <= cs$nVar) 
            {
                cs$new.vars<- cs$nextVar:lastVar
                cs$first.in.model[cs$new.vars]<- cs$currIter + 1
                cs$nextVar <- lastVar + 1
            } 
            # add in one by one until no variables left
            else 
            {
                cs$new.vars<- NULL
                for (i in length(cs$rmVector):1) 
                {
                   if (cs$rmVector[i] == TRUE && cs$nextVar <= cs$nVar) 
                   {
                        cs$new.vars<- c(cs$new.vars, cs$nextVar)
                        cs$first.in.model[cs$nextVar]<- cs$currIter + 1
                        cs$nextVar <- cs$nextVar + 1
                   }
                }
            }
        }
        else 
        {
            # exhausted all data
            cs$stopVar <- 1
            cs$new.vars = NULL
        }
    }
 
   # if we have finished (all variables) do some wrap-up and generate output values 
   if (cs$stopVar == 1) 
   {
   
        printC("finished iterating")
        
        
        
    
        currentX<- cs$sortedX[,cs$currentSet]
        colnames(currentX)<- colnames(cs$sortedX)[cs$currentSet]
        ret.bic.glm <- bic.glm (x = currentX, y = cs$Y, glm.family= cs$glm.family, maxCol = cs$maxNvar + 1, ...)
        
        
         
        output<- cs
        output$bma<- ret.bic.glm
        output$selected<- cs$currentSet
        output$nIterations<- cs$currIter
        class(output)<- "iBMA.glm"
   }
   
   else
   {
        output<- cs
        class(output)<- "iBMA.intermediate.glm"
   }
   output
}





iBMA.glm.matrix<- iBMA.glm.data.frame

