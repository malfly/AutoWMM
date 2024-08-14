##############################################################################
##############################################################################
####################           AutoWMM Package        ########################
#################### Automation of WMM Method on Tree ########################
####################        Structured Data           ########################
####################       By: Mallory Flynn          ########################
##############################################################################
##############################################################################

# Assuming a specific structure, create a tree with the following columns:
# from (node label), to (node label), Estimate (+ integer), Total (+ integer), 
# and Count (for terminal nodes with marginal counts). 
# 'from' and 'to' describe the edge for that row of data, where 'Estimate'
# and 'Total' are assumed to come from surveys of size 'Total' (a sample of
# the population at node 'from'), and observe 'Estimate' number of those
# individuals at 'Total' which move to the node described by 'to'.
# 'Estimate' and 'Total' columns are used for branching probabilities only.
# 'Count' column is NA for rows where 'to' nodes are not leaves; and also
# for all leaves without a marginal count.
# A Population (logical) column is not needed, but can be added if 'Estimate' 
# and 'Total' come from population numbers, rather than samples.
# A 'Description' column (string) is also possible to include if particulars 
# are desired on the tree diagram.
# 'TerminalCount' (binary) will be created for functional purposes, where
# marginal counts are included on leaves.

# Load libraries
if (!require(data.tree)){
  install.packages("data.tree")
  install.packages("DiagrammeR")
  library(data.tree)
} 

if (!require(rlang)){
  install.packages("rlang")
} 

if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}

if(!require(matlib)){
  install.packages("matlib")
  library(matlib)
}

if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

if(!require(networkD3)){
  install.packages("networkD3")
  library(networkD3)
}

if(!require(matrixStats)){
  install.packages("matrixStats")
  library(matrixStats)
}

##############################################################################
##############################################################################
# 1. Helper functions
##############################################################################
##############################################################################
# 1.1 Make tree
makeTree <- function(data){
  # check structure of data here
  if(!is.data.frame(data)){
    print('data must be dataframe type.')
  }
  # if no column x, print 'no column x'
  if(is.null(data$from)){
    print('data must have \'from\' column')
  }
  if(is.null(data$to)){
    print('data must have \'to\' column')
  }
  if(is.null(data$Estimate)){
    print('data must have \'Estimate\' column')
  }
  if(is.null(data$Total)){
    print('data must have \'Total\' column')
  }
  if(is.null(data$Count)){
    print('data must have \'Count\' column')
  }
  if(all(is.na(data$Total))){
    print('\'Total\' values cannot all be NA')
  }
  if(all(is.na(data$Estimate))){
    print('\'Estimate\' values cannot all be NA')
  }
  # if no 'Population' column exists, assume all values are samples
  if(is.null(data$Population)){
    data$Population <- rep(FALSE, times = dim(data)[1])
    print('WARNING: No \'Population\' column exists. We assume all 
          values are sample estimates')
  }
  
  # Now create tree structure
  tree <- FromDataFrameNetwork(data)
  tree$Do(function(node){
    node$probability <- round(node$Estimate/node$Total,2)
    # Create a binary value called 'TerminalCount' which is true for leaves with
    # marginal count data only
    if(node$isLeaf & !is.null(node$Count)){
      node$TerminalCount <- TRUE
    }else{
      node$TerminalCount <- FALSE
    }
  })
  return(tree)
}

##############################################################################
# 1.2 Visualize tree with descriptions and probabilities
drawTree <- function(tree, probs=TRUE, desc = TRUE){
  SetGraphStyle(tree,scale=2)
  if(probs){
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2,
                 label = function(node) round(mean(node$probability), 2))
  }else{
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2)
  }
  tree$Do(function(node){
    if(!is.null(node$Description) && desc){ # incorporate description on visual diagram
      SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                   label = function(node) node$Description) 
    }else{
      SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                   label = function(node) node$to) 
    }
  })
  plot(tree)
}


##############################################################################
# 1.3 Helper function. Method for sampling from a Beta distribution given the 
# survey estimates. 
# what happens if estimate or total is empty????
sampleBeta <- function(x,n,pop,node){ 
  result <- NA
  # if x,n are population values and not sample estimates, set probability exactly
  if(pop){                   
    result <- x/n
  }
  # else, sample a probability from a beta distribution with parameters below
  else{
    result <- rbeta(1,x + 1,n - x + 1)
  }
  return(result)
}

##############################################################################
# 1.4 Helper functions. Performs the multiplier method from a single terminal 
# node (o) and returns the root estimate given that path, and the probabilities 
# of each branch on that path.
mmEstimate <- function(o,test=FALSE){
  if(o$isLeaf){
    # first, generate an estimate of the root of the leaf
    estimate <- o$Count*(1/o$probability) 
    current_node <- o$parent # set the current node to now be that parent
    
    # while current node exists (ie the path back to the root is not yet
    # done), keep moving backwards and generating estimates for successively
    # higher nodes, until you reach the root, and generate a root estimate
    # this while loop performs the multiplicative back-calculation of one path
    # which goes from node o back to the root
    while(!is.null(current_node$parent)){ 
      if(test) print(current_node)        
      estimate <- estimate *(1/current_node$probability)
      current_node <- current_node$parent
      if(test) print(estimate)
    }
    o$targetEst <- estimate  # o$targetEst is the path-specific estimate
  }
  
  # The following assures an estimate of the root if not generated if 
  # the node with marginal count is not a leaf.
  else{
    o$targetEst <- NA
  }
}

##############################################################################
# 1.5 Helper function. Performs the closed form calculation of variances and means
# based on a "single-source sibling" tree.  Engages when single.source = TRUE 
# in wmmTree function.  See documentation for further details.
ssEstimate <- function(o,test=FALSE){
  if(o$isLeaf){
    # first, multiply the marginal count by first branch segment's contribution to product
    meanCalc <- o$Count*((o$Total + o$Estimate - 1)/(o$Estimate - 1))
    varCalc <- (o$Count)^2*(((o$Total + o$Estimate - 1)*(o$Total + o$Estimate - 2))/
                              ((o$Estimate - 1)*(o$Estimate - 2)))
    current_node <- o$parent # set the current node to now be that parent
    
    # while current node exists (ie the path back to the root is not yet
    # done), keep moving backwards and multiplying by the quotient for successively
    # higher nodes, until you reach the root, and complete the mean and var calculation
    # this while loop performs the multiplicative back-calculation of one path
    # which goes from node o back to the root
    while(!is.null(current_node$parent)){ 
      if(test) print(current_node)        
      meanCalc <- meanCalc*((current_node$Total + current_node$Estimate - 1)/(current_node$Estimate - 1))
      varCalc <- varCalc*(((current_node$Total + current_node$Estimate - 1)*
                             (current_node$Total + current_node$Estimate - 2))/
                            ((current_node$Estimate - 1)*(current_node$Estimate - 2)))
      current_node <- current_node$parent
      #if(test) print(estimate)
    }
    o$targetEst <- meanCalc  # o$targetEst is the path-specific estimate
    o$variance <- varCalc
  }
  
  # The following assures an estimate of the root if not generated if 
  # the node with marginal count is not a leaf.
  else{
    o$targetEst <- NA
    o$variance <- NA
  }
}




##############################################################################
##############################################################################
# 2. Generate weighted estimates
##############################################################################
##############################################################################
##############################################################################
# 2.1 Main function, including many cases for study data
wmmTree <- function(tree, sample_length = 10, method ='mmEstimate', 
                    int.type ='quants', single.source = FALSE){
  # choose which method to use - currently supports mmEstimate
  methodFunction <- NULL
  if(method =='mmEstimate'){
    print('using variance-weighted mean with multiplier method sampled path estimates')
    methodFunction <- mmEstimate
  }else{
    stop(paste(method,'is not a known method'))
  }
  
  tree$Set(targetEst_samples=numeric())
  tree$Set(probability_samples=numeric())
  tree$Set(imp_weights=numeric())
  
  if(single.source){
    # If all informative paths are obtained using a single, fully informative
    # source for all sibling data, closed form calculation can be used
    print('using closed-form expressions to generate estimates - be sure tree satisfies required assumptions')
    
    # Use parameters for each Dirichlet/Beta distributions from data table
    methodFunction <- ssEstimate
    
    # calculate target estimates from each leaf based on method above
    tree$Do(methodFunction, traversal = "post-order")
    
    # add targetEst to target estimates sample list
    tree$Do(function(node){
      node$targetEst_samples <- c(node$targetEst_samples,
                                  node$targetEst)
    })
    
    # set targetEst_samples to numeric() if the entire vector is NA
    # use for drawing functions
    tree$Do(function(node){
      if(all(is.na(node$targetEst_samples))){
        node$targetEst_samples <- numeric()
      }
    })
    
    # generate outputs
    # extract mean values for each path (column) (only leaves are required for
    # because that is where the estimates are stored)
    means <- tree$Get('targetEst', filterFun = function(node) node$isLeaf, 
                      traversal = 'post-order')
    vars <- tree$Get('variance', filterFun = function(node) node$isLeaf, 
                     traversal = 'post-order')
    
    ## incase the above is not in matrix form...
    if(is.list(means)){
      mat.means <- NULL
      for (i in 1:length(means)) {
        if (length(means[[i]]) > 0) {
          mat.means <- cbind(mat.means, means[[i]])
          colnames(mat.means)[dim(mat.means)[2]] <- names(means)[[i]]
        }
      }
      means <- mat.means
    }
    if(is.list(vars)){
      mat.vars <- NULL
      for (i in 1:length(vars)) {
        if (length(vars[[i]]) > 0) {
          mat.vars <- cbind(mat.vars, vars[[i]])
          colnames(mat.vars)[dim(mat.vars)[2]] <- names(vars)[[i]]
        }
      }
      vars <- mat.vars
    }
    
    # set as data frame
    means <- as.data.frame(means)
    vars <- as.data.frame(vars)
    
    ## select only leaves with marginal counts
    getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, 
                                traversal = 'post-order'))
    
    # if number of columns of x is greater than number of leaves, choose leaves only
    if(dim(means)[2]>length(getleaves)){
      means <- means %>% 
        select(all_of(getleaves))
    }
    if(dim(vars)[2]>length(getleaves)){
      vars <- vars %>% 
        select(all_of(getleaves))
    }
    
    # get mean estimate
    prec <- 1/vars
    w <- prec/sum(prec)
    rootEst <- as.matrix(means) %*% t(w) 
    rootVar <- as.matrix(vars) %*% t(w^2)
    
    # final values of the root are retained in 'estimate' and 'variance' at root
    tree$Do(function(node){
      if(isRoot(node)){node$Estimate <- round(rootEst, 2)}
    })
    tree$Do(function(node){
      if(isRoot(node)){node$variance <- rootVar}
    })
    
    m <- log(tree$Get('Estimate', filterFun = isRoot))
    
    # add 95% confidence intervals
    tree$Do(function(node) {
      node$uncertainty <- ss.confInts(node)
    })
    
  }else{
    for(m in 1:sample_length){
      # create sample probabilities among sibling bunches, so that sampled
      # branching probabilites add to 1 within a sibling group
      tree$Do(function(node){
        if(!is.null(node$children)){
          # start with sibling group of branches 
          siblist <- node$children
          k <- length(siblist)
          est.vec <- numeric(k) # create vector to save parameter values for dirichlet
          # or accepted beta draws
          total.vec <- numeric(k) # create vector to save 'Total' values for siblings
          
          # Get 'Estimate' and 'Total' from each sibling
          for(i in 1:k){
            childNode <- siblist[[i]]
            if(is.null(childNode$Estimate) | is.null(childNode$Total)){
              est.vec[i] <- NA
              total.vec[i] <- NA
            }else{
              est.vec[i] <- childNode$Estimate
              total.vec[i] <- childNode$Total
            }
          }
          
          # remove siblings from sib_list with NA 'estimate' (won't be sampled)
          # remove same entries from est.vec, total.vec
          na_index <- which(is.na(est.vec))
          est_childs <- which(!is.na(est.vec))
          sub.siblist <- siblist
          
          if(!rlang::is_empty(na_index)){
            sub.siblist <- siblist[-na_index]
            est.vec <- est.vec[-na_index]
            total.vec <- total.vec[-na_index]
          }
          
          # create sample_vec where samples will go
          sample.vec <- numeric(length(sub.siblist))
          names(sample.vec) <- names(sub.siblist)
          
          # FIRST CASE: IF any sibling 'estimate' in the group is NA,
          # then at least one branch is uninformed
          if(length(sub.siblist) < k){
            # while sib_list is non-empty:
            while(!rlang::is_empty(sub.siblist)){
              # start with first informed branch, current.branch<-sub.siblist[1]
              # find other branches informed by same study
              branch.samps <- which(total.vec == total.vec[1])
              
              # IF 'total' is the same for any other informed branches
              # && sum('estimate')<'total' for those branches
              if(length(branch.samps)>1 && sum(est.vec[branch.samps])<total.vec[1]){
                dir.params <- c(est.vec[branch.samps] + 1,
                                total.vec[1] - sum(est.vec[branch.samps]) + 1)
                probs <- rdirichlet(1, dir.params)
                
                # discard last probability (represents complement branches)
                probs <- probs[1:length(probs)-1]
              }
              
              else{
                probs <- rbeta(1, est.vec[1] + 1, total.vec[1] - est.vec[1] + 1)
              }
              
              # Now decide whether to accept or reject probs
              # IF sum(sample.vec) and new sampled probs < 1:
              if(sum(sample.vec) + sum(probs) < 1){
                # take branch number(s) off lists
                est.vec <- est.vec[-branch.samps]
                total.vec <- total.vec[-branch.samps]
                
                # set sample.vec equal to samples at the correct branch
                # positions
                sampled <- names(sub.siblist[branch.samps])
                sample.vec[names(sample.vec) %in% sampled] <- probs
                sub.siblist <- sub.siblist[-branch.samps]
              }
              # else: reject samples.  branches stay on sub.siblist and loop
              # re-runs to try sampling these branches again
            }
            
            # once we cycle through sampling for each informed sibling branch,
            # set 'probability' to be equal to the sample.vec
            t <- 1
            for(i in est_childs){
              node$children[[i]]$probability <- sample.vec[t]
              t <- t+1
            }
          }
          
          # SECOND CASE: ELSE, sub.siblist==siblist, so all branches are 
          # informed. In addition to rejection sampling, we also use  
          # importance sampling in this case
          else{
            # start with first branch (special case if Dirichlet)
            # find other branches informed by same study
            branch.samps <- which(total.vec == total.vec[1])
            
            # IF 'total' is the same for all branches && 
            # sum('estimate')='total' for those branches assume one 
            # study informs all branches. use a
            # Dir(Estimate_sibling1 + 1, Estimate_sibling2 + 1, ...) 
            # to sample those branches, save as sample.vec
            if(length(branch.samps)==length(sub.siblist) && 
               sum(est.vec[branch.samps])==total.vec[1]){
              dir.params <- c(est.vec[branch.samps] + 1)
              probs <- rdirichlet(1, dir.params)
              sample.vec[branch.samps] <- probs
              
              # set 'probability' value in tree
              for(i in 1:k){
                node$children[[i]]$probability <- sample.vec[i]
              }
            }
            
            # ELSE, a single study does not inform all branches
            # so while there is at least one branch left on sub.siblist 
            # we sample accordingly (by Dirichlet if same 'total' among 
            # some branches, with sum('estimate')< 'total' for those 
            # branches, or by Beta if not).
            else{
              # must create a process here where sampling occurs many times over
              # to initiate importance sampling.
              # set new names for sample.vec, total.vec, est.vec so they can be reset
              num.samp <- 100
              samp.siblist <- sub.siblist
              group.samples <- matrix(nrow = num.samp, ncol = length(sub.siblist))
              
              sample.vec.imp <- sample.vec
              est.vec.imp <- est.vec
              total.vec.imp <- total.vec
              
              # first generate num.samp samples of sibling group probs from the
              # 'wrong' distribution (sample each branch, deterministically
              # set the last branch of the group)
              for(w in 1:num.samp){
                while(length(samp.siblist)>1){
                  # start with first informed branch, current.branch<-sub.siblist[1]
                  # find other branches informed by same study
                  
                  # IF 'total' is the same for any other informed branches
                  # && sum('estimate')<1 for those branches
                  if(length(branch.samps)>1 && sum(est.vec[branch.samps])<total.vec.imp[1]){
                    dir.params <- c(est.vec.imp[branch.samps] + 1,
                                    total.vec.imp[1] - sum(est.vec.imp[branch.samps]) + 1)
                    probs <- rdirichlet(1, dir.params)
                    
                    # discard last probability (represents complement branches)
                    probs <- probs[1:length(probs)-1]
                    
                    # if length of probs is the same as the number of branches 
                    # remaining to sample, take last branch prob off as it will
                    # be set deterministically
                    if(length(probs)==length(samp.siblist)){
                      probs <- probs[1:length(probs)-1]
                      branch.samps <- branch.samps[1:length(branch.samps)-1]
                    }
                  }
                  
                  # ELSE, there are no other branches with same 'total' 
                  # or doesn't sum in way that suggest same study, and we
                  # sample branch as Beta(Estimate + 1, Total - Estimate + 1) 
                  else{
                    branch.samps <- branch.samps[1]
                    probs <- rbeta(1, est.vec.imp[1] + 1, 
                                   total.vec.imp[1] - est.vec.imp[1] + 1)
                  }
                  
                  # Now decide whether to accept or reject probs.  IF sum
                  # of existing sample.vec and sum of probs sampled < 1, accept
                  # samples and update est.vec, total.vec, samp.siblist so indices
                  # continue to match
                  if(sum(sample.vec.imp) + sum(probs) < 1){
                    # take branch number(s) off lists
                    est.vec.imp <- est.vec.imp[-branch.samps]
                    total.vec.imp <- total.vec.imp[-branch.samps]
                    
                    # set sample.vec equal to samples at the correct branch
                    # positions
                    sampled <- names(samp.siblist[branch.samps])
                    sample.vec.imp[names(sample.vec.imp) %in% sampled] <- probs
                    samp.siblist <- samp.siblist[-branch.samps]
                  }
                  # ELSE, we reject the samples and loop through again
                }
                
                # now there should be one branch left, to be set deterministically
                # set last branch on list equal to 1-sum(sample.vec) in sample.vec
                sample.vec.imp[length(sample.vec.imp)] <- 1-sum(sample.vec.imp)
                
                # save samples by row in group.samples matrix
                group.samples[w,] <- sample.vec.imp
                
                # reset sibling list for next iterations in for loop
                samp.siblist <- sub.siblist
                
                # reset sample.vec.imp, est.vec.imp, total.vec.imp
                sample.vec.imp <- sample.vec
                est.vec.imp <- est.vec
                total.vec.imp <- total.vec
                
              }# end of sub-sampling for importance scheme
              
              # calculate the importance weights associated with each set
              p.last <- group.samples[,dim(group.samples)[2]]
              imp.weights <- ((p.last)^(est.vec[1])) * (1-p.last)^(total.vec[1]-est.vec[1])
              
              #browser()
              #log.imp.weights1 <- (est.vec[1])*log(p.last) 
              #log.imp.weights2 <- (total.vec[1]-est.vec[1])*log(1-p.last)
              #logsum.imp.weights <- logSumExp(((p.last)^(est.vec[1])) * (1-p.last)^(total.vec[1]-est.vec[1]))
              #log.imp.weights <- log.imp.weights - logSumExp(((p.last)^(est.vec[1])) * (1-p.last)^(total.vec[1]-est.vec[1]))
              #imp.weights <- exp(log.imp.weights)
              
              # if sum of imp.weights is 0, computation error - use logged values
              if (sum(imp.weights)==0){
                log.imp.weights <- (est.vec[1])*log(p.last) + (total.vec[1]-est.vec[1])*log(1-p.last)
                log.imp.weights <- log.imp.weights - min(log.imp.weights)
                imp.weights <- exp(log.imp.weights)
              }
              imp.weights <- imp.weights/(sum(imp.weights))
              
              # if any imp.weights are NA, compute error resulting from large variance in p.last samples and extreme
              # values.  set imp.weights uniformly
              if (any(is.nan(imp.weights))){
                imp.weights <- rep(1/length(imp.weights), length(imp.weights))
              }
              
              # then create sampling scheme where we choose randomly from each of 
              # the samples we generated, weighted by the normalized importance weights
              # use that sample as the 'right' one for this pass, and move on as usual
              sample.vec <- group.samples[sample(x = 1:length(imp.weights),
                                                 size = 1, prob = imp.weights),]
              
              # set 'probability' value within tree to be equal to the chosen sample
              # from importance weighting scheme
              for(i in 1:k){
                node$children[[i]]$probability <- sample.vec[i]
                node$children[[i]]$impweight <- imp.weights
              } 
            }          
          } # end of case 2
        } # end of process for one node with children
      }) # end of node and treeDo function
      
      tree$Do(function(node){
        if(node$isRoot | is.null(node$Estimate)){
          node$probability <- NA
        }
      })
      
      # calculate target estimates from each leaf based on method above
      tree$Do(methodFunction, traversal = "post-order")
      
      # add targetEst to target estimates sample list
      tree$Do(function(node){
        node$targetEst_samples <- c(node$targetEst_samples,
                                    node$targetEst)
        node$probability_samples <- c(node$probability_samples,
                                      node$probability)
        node$imp_weights <- c(node$imp_weights,
                              node$impweight)
      })
      m <- m + 1
    } # end of for loop up to sample_length
    
    # Importance weight is the same within a sibling group, on each pass of wmm 
    # Calculate normalized imp_weights for each group of siblings, per pass
    tree$Do(function(node){
      node$imp_weights <- (node$imp_weights)/(sum(node$imp_weights))
    })
    
    # set targetEst_samples to numeric() if the entire vector is NA
    # use for drawing functions
    tree$Do(function(node){
      if(all(is.na(node$targetEst_samples))){
        node$targetEst_samples <- numeric()
      }
    })
    
    # after method is applied to each of the leaf node, targetEst_samples of each
    # node are used to calculate weights
    w <- ko.weights(tree)
    
    # weights are multiplied by mean estimates of each path to calculate a
    # final estimate of the root
    m <- logEstimates(tree)
    
    # final estimate of the root is retained in a value called 'estimate' at root
    tree$Do(function(node){
      if(isRoot(node)){node$Estimate <- round(exp(mean(m)), 2)}
    })
    
    # add 95% confidence intervals
    tree$Do(function(node) {
      node$uncertainty <- if(isRoot(node)){root.confInt(tree, int.type)}else{confInts(node$targetEst_samples)}
    })
  }
  # Finally, output the full list of Nhat estimates
  output <- list("root" = tree$Get("Estimate", filterFun = isRoot), 
                 "uncertainty" = tree$Get("uncertainty", filterFun = isRoot),
                 "estimates" = round(exp(m),2), 
                 "weights" = w)
  return(output)
}




##############################################################################
# 2.2 Visualize tree with root estimate and marginal counts, post-analysis
# also displays average of probability samples on each branch
countTree <- function(tree){
  # check if probability samples are empty everywhere - this indicates
  # weightedTree has not yet been used
  if(is.na(tree$Get('Estimate', filterFun = isRoot))){
    print('weightedTree() function has not yet been applied. conduct root estimation first.')
  }else{  
    SetGraphStyle(tree,scale=2)
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2,
                 label = function(node) round(mean(node$probability_samples), digits = 2))
    SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                 label = function(node) if(isRoot(node)){round(node$Estimate)}else{node$Count})
    plot(tree)
  }
}


##############################################################################
# 2.3 Visualize tree with root estimate given by each branch, and weighted sum 
# at the root (for post-analysis)
# Also displays average of probability samples on each branch.
estTree <- function(tree){
  # check if probability samples are empty everywhere - this indicates
  # weightedTree has not yet been used
  if(is.na(tree$Get('Estimate', filterFun = isRoot))){
    print('weightedTree() function has not yet been applied. conduct root estimation first.')
  }else{
    SetGraphStyle(tree,scale=2)
    SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2,
                 label = function(node) round(mean(node$probability_samples), digits = 2))
    SetNodeStyle(tree, fontsize=25, penwidth=3,width=1,
                 label = function(node) if(isRoot(node)){round(node$Estimate)}
                 else{round(mean(node$targetEst_samples))})
    plot(tree)
  }
}



##############################################################################
##############################################################################
# 3. Estimation helper functions
##############################################################################
##############################################################################
# 3.1 Method that takes samples and generates confidence intervals for 
# node other than the root. Assume raw data (not log), with normal distributed
# log data for confidence interval construction
confInts <- function(v){
  v <- log(v)
  uc <- quantile(v,0.975,na.rm = TRUE)
  lc <- quantile(v,0.025,na.rm = TRUE)
  int <- c(exp(lc), exp(uc))
  names(int) <- c("lower", "upper")
  
  return(int)
}

##############################################################################
# 3.2 Method for generating confidence interval for root node.  Assume unlogged
# input data and normally distributed logged data.  Completes conversion 
# internally.
root.confInt <- function(tree, int.type){
  # extract sample values for each path (column) (only leaves are required for
  # because that is where the estimates are stored)
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf, 
                traversal = 'post-order')
  
  ## incase the above is not in matrix form...
  if(is.list(x)){
    mat.x <- NULL
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        mat.x <- cbind(mat.x, x[[i]])
        colnames(mat.x)[dim(mat.x)[2]] <- names(x)[[i]]
      }
    }
    x <- mat.x
  }
  
  # set as data frame
  x <- as.data.frame(log(x))
  
  ## select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, 
                              traversal = 'post-order'))
  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>% 
      select(all_of(getleaves))
  }
  
  # get mean estimate
  m <- tree$Get('Estimate', filterFun = isRoot)
  
  # calculate interval based on int.type argument from wmmTree
  if(int.type =='var'){   #provides variance based interval
    print('using variance-weighted confidence interval - compare to quantiles')
    
    # calculate variance
    sig <- cov(log(x)) # covariance matrix
    #prec <- Inverse(sig) # precision matrix
    prec <- ginv(sig) # precision matrix
    e <- seq(1,1,length.out = dim(sig)[1])
    den.w <- t(e)%*%prec%*%e
    var <- 1/den.w[1]
    
    # generate endpoints
    lc <- log(m)-2*sqrt(var)
    uc <- log(m)+2*sqrt(var)
    lc <- exp(lc)
    uc <- exp(uc)
  }else if(int.type == 'cox'){        # provides cox interval
    print('using cox interval for uncertainty - compare to quantiles')
    
    # get weights
    weights <- as.matrix(ko.weights(tree))
    
    # get logN for each sample
    logrootEsts <- as.data.frame(as.matrix(x) %*% t(weights))
    
    # calculate sample variance of estimates
    nsamps <- len(logrootEsts)
    s2 <-sum((logrootEsts - log(m))^2)
    s2 <- s2/(nsamps-1)
    
    # generate interval endpoints
    m <- exp(log(m) + s2/2)
    lc <- m - 1.96*sqrt((s2/nsamps) + (s2^2/(2*(nsamps-1))))
    uc <- m + 1.96*sqrt((s2/nsamps) + (s2^2/(2*(nsamps-1))))
  }else{        # provides central 95% using quantiles
    # get weights
    weights <- as.matrix(ko.weights(tree))
    
    # get logN for each sample
    logrootEsts <- as.data.frame(as.matrix(x) %*% t(weights))
    
    uc <- exp(quantile(logrootEsts,0.975,na.rm = TRUE))
    lc <- exp(quantile(logrootEsts,0.025,na.rm = TRUE))
  }
  
  # return printed interval
  int <- c(round(lc),round(uc))
  names(int) <- c("lower", "upper")
  return(int)
}

##############################################################################
# 3.3 Creates a vector of mean estimate values given by each path
# where the leaf has a marginal count
logEstimates <- function(tree){
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf, 
                traversal = 'post-order')
  
  ## incase the above is not in matrix form...
  if(is.list(x)){
    mat.x <- NULL
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        mat.x <- cbind(mat.x, x[[i]])
        colnames(mat.x)[dim(mat.x)[2]] <- names(x)[[i]]
      }
    }
    x <- mat.x
  }
  
  # set as data frame
  x <- as.data.frame(log(x))
  
  ## select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, 
                              traversal = 'post-order'))
  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>% 
      select(all_of(getleaves))
  }
  
  weights <- as.matrix(ko.weights(tree)) #weights are internally calculated with log(data)
  logNhats <- as.matrix(x) %*% t(weights)
  
  return(logNhats)
}


##############################################################################
# 3.4 returns Nhat for each sample from the WMM (rather than the aggregate
# value given by the average, this calculates weights and applies the 
# weighted sum to each of the samples)
Nhats <- function(tree){
  weights <- as.matrix(ko.weights(tree))

  # extract sample values for each path with terminal node counts and
  # informative paths
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf,
                traversal = 'post-order')

  ## incase the above is not in matrix form...
  if(is.list(x)){
    mat.x <- NULL
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        mat.x <- cbind(mat.x, x[[i]])
        colnames(mat.x)[dim(mat.x)[2]] <- names(x)[[i]]
      }
    }
    x <- mat.x
  }

  # set as data frame
  x <- as.data.frame(log(x))

  # select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf,
                              traversal = 'post-order'))

  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>%
      select(all_of(getleaves))
  }

  # calculate logN for each sample
  logNhats <- as.data.frame(as.matrix(log(x)) %*% t(weights))

  # return unlogged values, Nhat
  return(round(exp(logNhats)))
}

##############################################################################
# 3.5 Calculating weights via Keller and Olkin solution
# This variant only assigns weights to leaves with population values
ko.weights <- function(tree){
  # extract sample values for each path with terminal node counts and
  # informative paths
  x <- tree$Get('targetEst_samples', filterFun = function(node) node$isLeaf, 
                traversal = 'post-order')
  
  
  ## incase the above is not in matrix form...
  if(is.list(x)){
    mat.x <- NULL
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        mat.x <- cbind(mat.x, x[[i]])
        colnames(mat.x)[dim(mat.x)[2]] <- names(x)[[i]]
      }
    }
    x <- mat.x
  }
  
  # make this a data frame
  x <- as.data.frame(x)
  
  # select only leaves with marginal counts
  getleaves <- which(tree$Get('TerminalCount', filterFun = isLeaf, 
                              traversal = 'post-order'))
  # if number of columns of x is greater than number of leaves, choose leaves only
  if(dim(x)[2]>length(getleaves)){
    x <- x %>% 
      select(all_of(getleaves))
  }
  
  # calculate weights
  sig <- cov(log(x)) # covariance matrix- use log to avoid numerical errors
  #prec <- Inverse(sig) # precision matrix 
  prec <- MASS::ginv(sig) # precision matrix 
  e <- seq(1,1,length.out = dim(sig)[1])
  num.w <- t(e)%*%prec
  den.w <- t(e)%*%prec%*%e
  w <- 1/den.w[1]*num.w
  colnames(w) <- names(getleaves)
  return(w)
}

##############################################################################
# 3.6 Method that takes samples and generates confidence intervals for 
# nodes in single source sibling tree (single.source = TRUE)
ss.confInts <- function(o, digits = 3){
  
  # get estimate, variance of root from node 
  m <- o$Estimate
  v <- o$variance
  
  # generate endpoints
  lc <- signif(m-1.96*sqrt(var), digits = digits)
  uc <- signif(m+1.96*sqrt(var), digits = digits)
  int <- c(lc, uc)
  names(int) <- c("lower", "upper")
  
  return(int)
}