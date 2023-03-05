# Script for the treatments in the article "L'homophilie sociale au college: amities et inimities entre
# eleves socialement distants dans quatre etablissements mixtes".
# Octobre 2022

# Run "setup general - w16" first to load raw data and create network objects.
# Script options : remove NAs, impute = soft.

########## 0 - Setup ##########
##### functions #####

create_metamat <- function(nets,fill_structural_zeros = NA, return_id_equivalence = T){
  
  if(class(nets)!="list"){
    print("INPUT SHOULD BE A LIST OF NETWORKS AND/OR MATRICES")
    break
  }
  
  # Convert to matrices if needed.
  for(p in 1:length(nets)){
    if(class(nets[[p]])=="network"){nets[[p]] <- as.sociomatrix(nets[[p]])}
  }
  
  # Create the empty metamat:
  mynames <- vector()
  
  for(p in 1:length(nets)){
    # add a numeric ID in case the same nodes are in different networks (e.g. if the networks are repeated observations over time)
    add_names <- paste(rownames(nets[[p]]),p,sep="_") 
    mynames <- c(mynames,add_names)
  }
  
  metamat <- matrix(nrow = length(mynames),ncol = length(mynames), data = fill_structural_zeros)
  rownames(metamat) <- mynames
  colnames(metamat) <- mynames
  
  # And fill :
  for(p in 1:length(nets)){
    add_names <- paste(rownames(nets[[p]]),p,sep="_") 
    metamat[add_names,add_names] <- nets[[p]]
  }
  
  res <- metamat
  
  # Also (optional) an id_equivalence matrix.
  if(return_id_equivalence==T){
    
    id_equivalence <- matrix(nrow = nrow(metamat),ncol = 2)
    id_equivalence[,1] <- rownames(metamat)
    w <- 0
    x <- 0
    for(p in 1:length(nets)){
      w <- x+1
      x <- x + nrow(nets[[p]])
      id_equivalence[w:x,2] <- rownames(nets[[p]])
    }
    res <- list(metamat,id_equivalence)
    names(res) <- c("metamat","id_equivalence")
  }
  
  return(res)
  
}

manual_time_qap <- function(net,indeps,id_equivalence,reps = 1000,return_all_coefs = F){
  
  # The point is to permute data with a time structure, ie where individuals can be present
  # at several time point. The function preserves that structure, by using cross-time-points
  # permutations (e.g. if id 7 becomes id 19, it does so on all time waves at once).
  
  # The input "id_equivalence" should be a matrix with 2 columns, and as many rows are there are single
  # ids in "net" (which should be a metamat created by function "create_metamat"). The first column gives the
  # time-individual unique id, and the second column gives the transwave id for that individual.
  # Rows should have the same ordering as net rownames and colnames.
  
  # NB : "indeps" can have names, which will be passed to the output (not necessary though).
  
  if(class(net)!="network" & class(net)!="matrix"){
    print("FIRST INPUT SHOULD BE A SINGLE NETWORK AND/OR MATRIX")
    break
  }
  
  if(class(indeps)!="list"){
    print("SECOND INPUT SHOULD BE A LIST OF ONE OR SEVERAL DYADIC COVARIATES")
    break
  }
  
  if(class(net)=="network"){net <- as.sociomatrix(net)} # turn networks into matrices
  
  if(length(which(rownames(net)!=id_equivalence[,1]))>0){
    print("ERROR : THE FIRST ROW OF ID_EQUIVALENCE DOES NOT CORRESPOND TO THE LABEL NAMES OF THE METAMATRIX PROVIDED")
    break
  }
  
  # Create result matrix.
  res <- matrix(nrow = reps,ncol = length(indeps)+1)
  if(is.null(names(indeps))==F){
    colnames(res) <- c("intercept",names(indeps))
  } else {colnames(res) <- c("intercept",paste("x",1:length(indeps),sep=""))}
  
  # Prepare the variables and formula for regression (only the dependent variable, i.e. the network will change).
  for(i in 1:length(indeps)){
    indeps[[i]] <- as.vector(indeps[[i]])
  }
  vars <- data.frame(indeps)
  names(vars) <- paste("x",1:length(indeps),sep="")
  
  formula <- paste(names(vars),collapse = "+")
  formula <- paste("new_net","~",formula)
  formula <- as.formula(formula)
  
  # Add a third (for now empty) row to id_equivalence.
  id_equivalence <- cbind(id_equivalence,rep("",nrow(id_equivalence)))
  
  # Time stamps.
  times <- unique(sapply(strsplit(id_equivalence[,1],"_"),"[",2))
  time_stamps <- list()
  for(t in times){
    time_stamps[[t]] <- which(sapply(strsplit(id_equivalence[,1],"_"),"[",2)==t)
  }
  rm(times)
  
  # Set up a progress bar.
  print(paste("Starting permutations with",reps,"replications..."))
  pb <- txtProgressBar(min = 0, max = reps, initial = 0,style=3,char="  ^.^' ") 
  
  # Now for the numbers of reps...
  for(it in 1:reps){
    
    setTxtProgressBar(pb,it)
    
    # Permute randomly the ids.
    # Warning : there is a little "trick", in that we want to make sure a) that permutations are done
    # only within a time-wave (don t want a wave 4 id to arrive in w1 for example) but b) that the
    # correspondence is the same across time waves.
    # Thus we start with t=1, permute, save the correspondences, then go to t=2, only permute the new ones, etc.
    
    # NB : the issue is that kids arrived at a specific wave are only permuted among themselves, and it s usually just a few of
    # them at once. But I couldn t find another solution...
    
    for(t in 1:length(time_stamps)){
      tempo_equivalence <- id_equivalence[time_stamps[[t]],]
      tempo_equivalence <- tempo_equivalence[tempo_equivalence[,3]=="",] # only keep rows that don t have an id yet
      
      if(class(tempo_equivalence)=="matrix"){ # only if there are at least 2 nodes to permute this wave.
        unique_id <- unique(tempo_equivalence[,2])
        permuted_id <- sample(unique_id,size=length(unique_id),replace=FALSE)
        
        # And the permuted ids are applied for ALL waves, not just the running one :
        for(z in 1:length(unique_id)){
          hop <- which(id_equivalence[,2]==unique_id[z])
          id_equivalence[hop,3] <- permuted_id[z]
        }
      } else { # otherwise the node stays at such.
        unique_id <- tempo_equivalence[2]
        hop <- which(id_equivalence[,2]==unique_id)
        id_equivalence[hop,3] <- unique_id
      }
    }
    
    # Now permute within each time step.
    
    new_net <- net
    for(t in 1:length(time_stamps)){
      
      rows <- time_stamps[[t]]
      tempo_net <- net[rows,rows]      
      
      # Name with the unique ids (not the time specific ones, since we are now considering only one wave).
      rownames(tempo_net) <- id_equivalence[rows,2]
      colnames(tempo_net) <- id_equivalence[rows,2]
      
      # Permute based on the news ids.
      new_order <- order(match(id_equivalence[rows,2], id_equivalence[rows,3]))
      tempo_net <- tempo_net[new_order,new_order]
      
      # And attribute in new_net.
      new_net[rows,rows] <- tempo_net
    }
    
    # So new_net is our new metamatrix, with all waves and properly permuted !
    # All that s left is to regress it on the indeps variables, and get the coefficients.
    
    # NB : we set the diagonal to 0, so that it s not used in the regression (I just use glm, after vectorizing
    # the matrix; structural zeos are NAs as well so they ll be left off).
    diag(new_net) <- NA
    
    vars$new_net <- as.vector(new_net)
    
    millesabords <- glm(formula,data = vars,family="binomial")
    res[it,] <- coef(millesabords)
    
    # Finally, we must not forget to "refresh" the 3rd column of the equivalence matrix, so the loop can start again next time.
    id_equivalence[,3] <- ""
    
  } # And done with the big loop at last.
  
  # Now we perform the true regression:
  diag(net) <- NA
  vars$new_net <- as.vector(net) 
  millesabords <- glm(formula,data = vars,family="binomial")
  
  # Output creation.
  short_res <- matrix(nrow = length(coef(millesabords)), ncol = 5)
  rownames(short_res) <- names(coef(millesabords))
  
  short_res[,1] <- coef(millesabords) # the true coefs
  short_res[,2] <- colMeans(res,na.rm=T) # average in the permutations
  
  standard_deviations <- apply(res,2,sd) # sd in the permutations
  
  short_res[,3] <- short_res[,2] - 1.96*standard_deviations
  short_res[,4] <- short_res[,2] + 1.96*standard_deviations # 95% confidence interval
  
  # And equivalent of the p-value (i.e the proportion of the distribution that is extreme in regard to the true value)
  for(x in 1:nrow(short_res)){
    if(short_res[x,1]>short_res[x,2]){ # if the observed value is superior to the mean of the sample...
      # ... count the proportion of sample values that are superior to the observed one.
      short_res[x,5] <- 2*length(which(res[,x]>short_res[x,1]))/reps
      
      # NB : we multiply the p value by 2 because this is a one-sided test.
      
    }
    if(short_res[x,1]<=short_res[x,2]){ # and the other way around.
      short_res[x,5] <- 2*length(which(res[,x]<short_res[x,1]))/reps
    }
  }
  
  colnames(short_res) <- c("observed","average simulated","95%_low","95%_high","p_value")
  
  if(return_all_coefs==T){
    full_res <- list(short_res,res)
    names(full_res) <- c("summary","coefficients")
    return(full_res)
  } else {return(short_res)}
  
  
}

##### a - create outgf networks #####

# We count all the gf ties that also correspond to a out or home tie.

outgf_rvl <- list()
outgf_jbs <- list()
outgf_pm <- list()
outgf_dmz <- list()

aa <- as.sociomatrix(myhome_rvl[[1]])+as.sociomatrix(hishome_rvl[[1]])+as.sociomatrix(outside_rvl[[1]])
aa[aa>1] <- 1
outgf_rvl[[1]] <- network(aa*as.sociomatrix(gf_rvl[[1]]))

aa <- as.sociomatrix(myhome_rvl[[2]])+as.sociomatrix(hishome_rvl[[2]])+as.sociomatrix(outside_rvl[[2]])
aa[aa>1] <- 1
outgf_rvl[[2]] <- network(aa*as.sociomatrix(gf_rvl[[3]]))

##

aa <- as.sociomatrix(myhome_jbs[[1]])+as.sociomatrix(hishome_jbs[[1]])+as.sociomatrix(outside_jbs[[1]])
aa[aa>1] <- 1
outgf_jbs[[1]] <- network(aa*as.sociomatrix(gf_jbs[[1]]))

aa <- as.sociomatrix(myhome_jbs[[2]])+as.sociomatrix(hishome_jbs[[2]])+as.sociomatrix(outside_jbs[[2]])
aa[aa>1] <- 1
outgf_jbs[[2]] <- network(aa*as.sociomatrix(gf_jbs[[3]]))

##

aa <- as.sociomatrix(myhome_pm[[1]])+as.sociomatrix(hishome_pm[[1]])+as.sociomatrix(outside_pm[[1]])
aa[aa>1] <- 1
outgf_pm[[1]] <- network(aa*as.sociomatrix(gf_pm[[1]]))

aa <- as.sociomatrix(myhome_pm[[2]])+as.sociomatrix(hishome_pm[[2]])+as.sociomatrix(outside_pm[[2]])
aa[aa>1] <- 1
outgf_pm[[2]] <- network(aa*as.sociomatrix(gf_pm[[3]]))

##

aa <- as.sociomatrix(myhome_dmz[[1]])+as.sociomatrix(hishome_dmz[[1]])+as.sociomatrix(outside_dmz[[1]])
aa[aa>1] <- 1
outgf_dmz[[1]] <- network(aa*as.sociomatrix(gf_dmz[[1]]))

aa <- as.sociomatrix(myhome_dmz[[2]])+as.sociomatrix(hishome_dmz[[2]])+as.sociomatrix(outside_dmz[[2]])
aa[aa>1] <- 1
outgf_dmz[[2]] <- network(aa*as.sociomatrix(gf_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

outgf_rvl[[1]] <- transfer_att(outgf_rvl[[1]],gf_rvl[[1]])
outgf_rvl[[2]] <- transfer_att(outgf_rvl[[2]],gf_rvl[[3]])

outgf_jbs[[1]] <- transfer_att(outgf_jbs[[1]],gf_jbs[[1]])
outgf_jbs[[2]] <- transfer_att(outgf_jbs[[2]],gf_jbs[[3]])

outgf_pm[[1]] <- transfer_att(outgf_pm[[1]],gf_pm[[1]])
outgf_pm[[2]] <- transfer_att(outgf_pm[[2]],gf_pm[[3]])

outgf_dmz[[1]] <- transfer_att(outgf_dmz[[1]],gf_dmz[[1]])
outgf_dmz[[2]] <- transfer_att(outgf_dmz[[2]],gf_dmz[[3]])

##### b - create pargf networks #####

# Friends whom parents your parents know.
# NB : we symmetrize parents declaration, weak rule (otherwise it will never be dense enough).

pargf_rvl <- list()
pargf_jbs <- list()
pargf_pm <- list()
pargf_dmz <- list()

as.sociomatrix(parents_rvl[[1]])

aa <- symmetrize(as.sociomatrix(parents_rvl[[1]]),rule="weak")
aa[aa>1] <- 1
pargf_rvl[[1]] <- network(aa*as.sociomatrix(gf_rvl[[1]]))

aa <- symmetrize(as.sociomatrix(parents_rvl[[2]]),rule="weak")
aa[aa>1] <- 1
pargf_rvl[[2]] <- network(aa*as.sociomatrix(gf_rvl[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_jbs[[1]]),rule="weak")
aa[aa>1] <- 1
pargf_jbs[[1]] <- network(aa*as.sociomatrix(gf_jbs[[1]]))

aa <- symmetrize(as.sociomatrix(parents_jbs[[2]]),rule="weak")
aa[aa>1] <- 1
pargf_jbs[[2]] <- network(aa*as.sociomatrix(gf_jbs[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_pm[[1]]),rule="weak")
aa[aa>1] <- 1
pargf_pm[[1]] <- network(aa*as.sociomatrix(gf_pm[[1]]))

aa <- symmetrize(as.sociomatrix(parents_pm[[2]]),rule="weak")
aa[aa>1] <- 1
pargf_pm[[2]] <- network(aa*as.sociomatrix(gf_pm[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_dmz[[1]]),rule="weak")
aa[aa>1] <- 1
pargf_dmz[[1]] <- network(aa*as.sociomatrix(gf_dmz[[1]]))

aa <- symmetrize(as.sociomatrix(parents_dmz[[2]]),rule="weak")
aa[aa>1] <- 1
pargf_dmz[[2]] <- network(aa*as.sociomatrix(gf_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

pargf_rvl[[1]] <- transfer_att(pargf_rvl[[1]],gf_rvl[[1]])
pargf_rvl[[2]] <- transfer_att(pargf_rvl[[2]],gf_rvl[[3]])

pargf_jbs[[1]] <- transfer_att(pargf_jbs[[1]],gf_jbs[[1]])
pargf_jbs[[2]] <- transfer_att(pargf_jbs[[2]],gf_jbs[[3]])

pargf_pm[[1]] <- transfer_att(pargf_pm[[1]],gf_pm[[1]])
pargf_pm[[2]] <- transfer_att(pargf_pm[[2]],gf_pm[[3]])

pargf_dmz[[1]] <- transfer_att(pargf_dmz[[1]],gf_dmz[[1]])
pargf_dmz[[2]] <- transfer_att(pargf_dmz[[2]],gf_dmz[[3]])

##### c - create outallf networks #####

# We count all the allf ties that also correspond to a out or home tie.

outallf_rvl <- list()
outallf_jbs <- list()
outallf_pm <- list()
outallf_dmz <- list()

aa <- as.sociomatrix(myhome_rvl[[1]])+as.sociomatrix(hishome_rvl[[1]])+as.sociomatrix(outside_rvl[[1]])
aa[aa>1] <- 1
outallf_rvl[[1]] <- network(aa*as.sociomatrix(allf_rvl[[1]]))

aa <- as.sociomatrix(myhome_rvl[[2]])+as.sociomatrix(hishome_rvl[[2]])+as.sociomatrix(outside_rvl[[2]])
aa[aa>1] <- 1
outallf_rvl[[2]] <- network(aa*as.sociomatrix(allf_rvl[[3]]))

##

aa <- as.sociomatrix(myhome_jbs[[1]])+as.sociomatrix(hishome_jbs[[1]])+as.sociomatrix(outside_jbs[[1]])
aa[aa>1] <- 1
outallf_jbs[[1]] <- network(aa*as.sociomatrix(allf_jbs[[1]]))

aa <- as.sociomatrix(myhome_jbs[[2]])+as.sociomatrix(hishome_jbs[[2]])+as.sociomatrix(outside_jbs[[2]])
aa[aa>1] <- 1
outallf_jbs[[2]] <- network(aa*as.sociomatrix(allf_jbs[[3]]))

##

aa <- as.sociomatrix(myhome_pm[[1]])+as.sociomatrix(hishome_pm[[1]])+as.sociomatrix(outside_pm[[1]])
aa[aa>1] <- 1
outallf_pm[[1]] <- network(aa*as.sociomatrix(allf_pm[[1]]))

aa <- as.sociomatrix(myhome_pm[[2]])+as.sociomatrix(hishome_pm[[2]])+as.sociomatrix(outside_pm[[2]])
aa[aa>1] <- 1
outallf_pm[[2]] <- network(aa*as.sociomatrix(allf_pm[[3]]))

##

aa <- as.sociomatrix(myhome_dmz[[1]])+as.sociomatrix(hishome_dmz[[1]])+as.sociomatrix(outside_dmz[[1]])
aa[aa>1] <- 1
outallf_dmz[[1]] <- network(aa*as.sociomatrix(allf_dmz[[1]]))

aa <- as.sociomatrix(myhome_dmz[[2]])+as.sociomatrix(hishome_dmz[[2]])+as.sociomatrix(outside_dmz[[2]])
aa[aa>1] <- 1
outallf_dmz[[2]] <- network(aa*as.sociomatrix(allf_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

outallf_rvl[[1]] <- transfer_att(outallf_rvl[[1]],allf_rvl[[1]])
outallf_rvl[[2]] <- transfer_att(outallf_rvl[[2]],allf_rvl[[3]])

outallf_jbs[[1]] <- transfer_att(outallf_jbs[[1]],allf_jbs[[1]])
outallf_jbs[[2]] <- transfer_att(outallf_jbs[[2]],allf_jbs[[3]])

outallf_pm[[1]] <- transfer_att(outallf_pm[[1]],allf_pm[[1]])
outallf_pm[[2]] <- transfer_att(outallf_pm[[2]],allf_pm[[3]])

outallf_dmz[[1]] <- transfer_att(outallf_dmz[[1]],allf_dmz[[1]])
outallf_dmz[[2]] <- transfer_att(outallf_dmz[[2]],allf_dmz[[3]])

##### d - create parallf networks #####

# Friends whom parents your parents know.
# NB : we symmetrize parents declaration, weak rule (otherwise it will never be dense enough).

parallf_rvl <- list()
parallf_jbs <- list()
parallf_pm <- list()
parallf_dmz <- list()

as.sociomatrix(parents_rvl[[1]])

aa <- symmetrize(as.sociomatrix(parents_rvl[[1]]),rule="weak")
aa[aa>1] <- 1
parallf_rvl[[1]] <- network(aa*as.sociomatrix(allf_rvl[[1]]))

aa <- symmetrize(as.sociomatrix(parents_rvl[[2]]),rule="weak")
aa[aa>1] <- 1
parallf_rvl[[2]] <- network(aa*as.sociomatrix(allf_rvl[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_jbs[[1]]),rule="weak")
aa[aa>1] <- 1
parallf_jbs[[1]] <- network(aa*as.sociomatrix(allf_jbs[[1]]))

aa <- symmetrize(as.sociomatrix(parents_jbs[[2]]),rule="weak")
aa[aa>1] <- 1
parallf_jbs[[2]] <- network(aa*as.sociomatrix(allf_jbs[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_pm[[1]]),rule="weak")
aa[aa>1] <- 1
parallf_pm[[1]] <- network(aa*as.sociomatrix(allf_pm[[1]]))

aa <- symmetrize(as.sociomatrix(parents_pm[[2]]),rule="weak")
aa[aa>1] <- 1
parallf_pm[[2]] <- network(aa*as.sociomatrix(allf_pm[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_dmz[[1]]),rule="weak")
aa[aa>1] <- 1
parallf_dmz[[1]] <- network(aa*as.sociomatrix(allf_dmz[[1]]))

aa <- symmetrize(as.sociomatrix(parents_dmz[[2]]),rule="weak")
aa[aa>1] <- 1
parallf_dmz[[2]] <- network(aa*as.sociomatrix(allf_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

parallf_rvl[[1]] <- transfer_att(parallf_rvl[[1]],allf_rvl[[1]])
parallf_rvl[[2]] <- transfer_att(parallf_rvl[[2]],allf_rvl[[3]])

parallf_jbs[[1]] <- transfer_att(parallf_jbs[[1]],allf_jbs[[1]])
parallf_jbs[[2]] <- transfer_att(parallf_jbs[[2]],allf_jbs[[3]])

parallf_pm[[1]] <- transfer_att(parallf_pm[[1]],allf_pm[[1]])
parallf_pm[[2]] <- transfer_att(parallf_pm[[2]],allf_pm[[3]])

parallf_dmz[[1]] <- transfer_att(parallf_dmz[[1]],allf_dmz[[1]])
parallf_dmz[[2]] <- transfer_att(parallf_dmz[[2]],allf_dmz[[3]])

##### e - create outf networks #####

# We count all the f ties that also correspond to a out or home tie.

outf_rvl <- list()
outf_jbs <- list()
outf_pm <- list()
outf_dmz <- list()

aa <- as.sociomatrix(myhome_rvl[[1]])+as.sociomatrix(hishome_rvl[[1]])+as.sociomatrix(outside_rvl[[1]])
aa[aa>1] <- 1
outf_rvl[[1]] <- network(aa*as.sociomatrix(f_rvl[[1]]))

aa <- as.sociomatrix(myhome_rvl[[2]])+as.sociomatrix(hishome_rvl[[2]])+as.sociomatrix(outside_rvl[[2]])
aa[aa>1] <- 1
outf_rvl[[2]] <- network(aa*as.sociomatrix(f_rvl[[3]]))

##

aa <- as.sociomatrix(myhome_jbs[[1]])+as.sociomatrix(hishome_jbs[[1]])+as.sociomatrix(outside_jbs[[1]])
aa[aa>1] <- 1
outf_jbs[[1]] <- network(aa*as.sociomatrix(f_jbs[[1]]))

aa <- as.sociomatrix(myhome_jbs[[2]])+as.sociomatrix(hishome_jbs[[2]])+as.sociomatrix(outside_jbs[[2]])
aa[aa>1] <- 1
outf_jbs[[2]] <- network(aa*as.sociomatrix(f_jbs[[3]]))

##

aa <- as.sociomatrix(myhome_pm[[1]])+as.sociomatrix(hishome_pm[[1]])+as.sociomatrix(outside_pm[[1]])
aa[aa>1] <- 1
outf_pm[[1]] <- network(aa*as.sociomatrix(f_pm[[1]]))

aa <- as.sociomatrix(myhome_pm[[2]])+as.sociomatrix(hishome_pm[[2]])+as.sociomatrix(outside_pm[[2]])
aa[aa>1] <- 1
outf_pm[[2]] <- network(aa*as.sociomatrix(f_pm[[3]]))

##

aa <- as.sociomatrix(myhome_dmz[[1]])+as.sociomatrix(hishome_dmz[[1]])+as.sociomatrix(outside_dmz[[1]])
aa[aa>1] <- 1
outf_dmz[[1]] <- network(aa*as.sociomatrix(f_dmz[[1]]))

aa <- as.sociomatrix(myhome_dmz[[2]])+as.sociomatrix(hishome_dmz[[2]])+as.sociomatrix(outside_dmz[[2]])
aa[aa>1] <- 1
outf_dmz[[2]] <- network(aa*as.sociomatrix(f_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

outf_rvl[[1]] <- transfer_att(outf_rvl[[1]],f_rvl[[1]])
outf_rvl[[2]] <- transfer_att(outf_rvl[[2]],f_rvl[[3]])

outf_jbs[[1]] <- transfer_att(outf_jbs[[1]],f_jbs[[1]])
outf_jbs[[2]] <- transfer_att(outf_jbs[[2]],f_jbs[[3]])

outf_pm[[1]] <- transfer_att(outf_pm[[1]],f_pm[[1]])
outf_pm[[2]] <- transfer_att(outf_pm[[2]],f_pm[[3]])

outf_dmz[[1]] <- transfer_att(outf_dmz[[1]],f_dmz[[1]])
outf_dmz[[2]] <- transfer_att(outf_dmz[[2]],f_dmz[[3]])

##### f - create parf networks #####

# Friends whom parents your parents know.
# NB : we symmetrize parents declaration, weak rule (otherwise it will never be dense enough).

parf_rvl <- list()
parf_jbs <- list()
parf_pm <- list()
parf_dmz <- list()

as.sociomatrix(parents_rvl[[1]])

aa <- symmetrize(as.sociomatrix(parents_rvl[[1]]),rule="weak")
aa[aa>1] <- 1
parf_rvl[[1]] <- network(aa*as.sociomatrix(f_rvl[[1]]))

aa <- symmetrize(as.sociomatrix(parents_rvl[[2]]),rule="weak")
aa[aa>1] <- 1
parf_rvl[[2]] <- network(aa*as.sociomatrix(f_rvl[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_jbs[[1]]),rule="weak")
aa[aa>1] <- 1
parf_jbs[[1]] <- network(aa*as.sociomatrix(f_jbs[[1]]))

aa <- symmetrize(as.sociomatrix(parents_jbs[[2]]),rule="weak")
aa[aa>1] <- 1
parf_jbs[[2]] <- network(aa*as.sociomatrix(f_jbs[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_pm[[1]]),rule="weak")
aa[aa>1] <- 1
parf_pm[[1]] <- network(aa*as.sociomatrix(f_pm[[1]]))

aa <- symmetrize(as.sociomatrix(parents_pm[[2]]),rule="weak")
aa[aa>1] <- 1
parf_pm[[2]] <- network(aa*as.sociomatrix(f_pm[[3]]))

##

aa <- symmetrize(as.sociomatrix(parents_dmz[[1]]),rule="weak")
aa[aa>1] <- 1
parf_dmz[[1]] <- network(aa*as.sociomatrix(f_dmz[[1]]))

aa <- symmetrize(as.sociomatrix(parents_dmz[[2]]),rule="weak")
aa[aa>1] <- 1
parf_dmz[[2]] <- network(aa*as.sociomatrix(f_dmz[[3]]))

###

# Transfer attributes :
transfer_att <- function(net1,net2){
  for(att in list.vertex.attributes(net2)){
    value <- get.vertex.attribute(net2,attrname = att)
    set.vertex.attribute(net1,att,value = value)
  }
  return(net1)
}

parf_rvl[[1]] <- transfer_att(parf_rvl[[1]],f_rvl[[1]])
parf_rvl[[2]] <- transfer_att(parf_rvl[[2]],f_rvl[[3]])

parf_jbs[[1]] <- transfer_att(parf_jbs[[1]],f_jbs[[1]])
parf_jbs[[2]] <- transfer_att(parf_jbs[[2]],f_jbs[[3]])

parf_pm[[1]] <- transfer_att(parf_pm[[1]],f_pm[[1]])
parf_pm[[2]] <- transfer_att(parf_pm[[2]],f_pm[[3]])

parf_dmz[[1]] <- transfer_att(parf_dmz[[1]],f_dmz[[1]])
parf_dmz[[2]] <- transfer_att(parf_dmz[[2]],f_dmz[[3]])

############################################################
################## 1 - Homophilie vague 1 ##################
############################################################
########## 0 - Various descriptives for the data section ##########

# Ecart moyen de factbrev entre PCS :
a <- mean(occupation$factbrev[occupation$pcs_chef_ag=="3"],na.rm=T)
b <- mean(occupation$factbrev[occupation$pcs_chef_ag=="6"],na.rm=T)
a
b
abs(a-b)

sd(occupation$factbrev,na.rm=T)

# Conversion ecart-type / notes :
sd(ind1$g_mean_T2,na.rm=T)

names(ind2)
sd(ind2$g_mean_T1,na.rm=T)
sd(ind3$g_mean_T2,na.rm=T)
sd(ind4$g_mean_T2,na.rm=T)

# Mean and median degree distributions :
aa_gf <- c(unlist(lapply(gf_rvl,degree,cmode="outdegree")),unlist(lapply(gf_jbs,degree,cmode="outdegree")),
           unlist(lapply(gf_pm,degree,cmode="outdegree")),unlist(lapply(gf_dmz,degree,cmode="outdegree")))
aa_allf <- c(unlist(lapply(allf_rvl,degree,cmode="outdegree")),unlist(lapply(allf_jbs,degree,cmode="outdegree")),
             unlist(lapply(allf_pm,degree,cmode="outdegree")),unlist(lapply(allf_dmz,degree,cmode="outdegree")))
aa_dl <- c(unlist(lapply(dl_rvl,degree,cmode="outdegree")),unlist(lapply(dl_jbs,degree,cmode="outdegree")),
           unlist(lapply(dl_pm,degree,cmode="outdegree")),unlist(lapply(dl_dmz,degree,cmode="outdegree")))
aa_outgf <- c(unlist(lapply(outgf_rvl,degree,cmode="outdegree")),unlist(lapply(outgf_jbs,degree,cmode="outdegree")),
              unlist(lapply(outgf_pm,degree,cmode="outdegree")),unlist(lapply(outgf_dmz,degree,cmode="outdegree")))
aa_pargf <- c(unlist(lapply(pargf_rvl,degree,cmode="outdegree")),unlist(lapply(pargf_jbs,degree,cmode="outdegree")),
              unlist(lapply(pargf_pm,degree,cmode="outdegree")),unlist(lapply(pargf_dmz,degree,cmode="outdegree")))
aa_bully <- c(unlist(lapply(bully_rvl,degree,cmode="outdegree")),unlist(lapply(bully_jbs,degree,cmode="outdegree")),
              unlist(lapply(bully_pm,degree,cmode="outdegree")),unlist(lapply(bully_dmz,degree,cmode="outdegree")))

mean(aa_gf)
sd(aa_gf)
median(aa_gf)

mean(aa_allf)
sd(aa_allf)
median(aa_allf)

mean(aa_dl)
sd(aa_dl)
median(aa_dl)

mean(aa_outgf)
sd(aa_outgf)
median(aa_outgf)

mean(aa_pargf)
sd(aa_pargf)
median(aa_pargf)

mean(aa_bully)
sd(aa_bully)
median(aa_bully)

# Recoupement bully*amis :
table(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(bully_rvl[[1]]))
9/(37+9)

table(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(bully_jbs[[1]]))
28/(176+28)

table(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(bully_pm[[1]]))
10/(57+10)

table(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(bully_dmz[[1]]))
12/(104+12)

#
##### prop de sous ensembles #####

myprop <- function(x,y){
  table(x,y)["1","1"]/sum(table(x,y)["1",])
}

res_prop <- matrix(nrow = 24, ncol = 6, data = NA)
rownames(res_prop) <- c("allf.p1","gf.p1","outallf.p1","parallf.p1","out.p1","par.p1",
                        "allf.p2","gf.p2","outallf.p2","parallf.p2","out.p2","par.p2",
                        "allf.s1","gf.s1","outallf.s1","parallf.s1","out.s1","par.s1",
                        "allf.s2","gf.s2","outallf.s2","parallf.s2","out.s2","par.s2")
colnames(res_prop) <- c("allf","outallf","parallf","gf","outgf","pargf")

aa <- as.sociomatrix(myhome_rvl[[1]])+as.sociomatrix(hishome_rvl[[1]])+as.sociomatrix(outside_rvl[[1]])
bb <- symmetrize(as.sociomatrix(parents_rvl[[1]]),rule="weak")

res_prop["out.p1","outallf"] <- myprop(aa,as.sociomatrix(outallf_rvl[[1]]))
res_prop["out.p1","outgf"] <- myprop(aa,as.sociomatrix(outgf_rvl[[1]]))
res_prop["par.p1","parallf"] <- myprop(bb,as.sociomatrix(parallf_rvl[[1]]))
res_prop["par.p1","pargf"] <- myprop(bb,as.sociomatrix(pargf_rvl[[1]]))
res_prop["allf.p1","outallf"] <- myprop(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(outallf_rvl[[1]]))
res_prop["allf.p1","parallf"] <- myprop(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(parallf_rvl[[1]]))
res_prop["allf.p1","gf"] <- myprop(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(gf_rvl[[1]]))
res_prop["allf.p1","outgf"] <- myprop(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(outgf_rvl[[1]]))
res_prop["allf.p1","pargf"] <- myprop(as.sociomatrix(allf_rvl[[1]]),as.sociomatrix(pargf_rvl[[1]]))
res_prop["gf.p1","outgf"] <- myprop(as.sociomatrix(gf_rvl[[1]]),as.sociomatrix(outgf_rvl[[1]]))
res_prop["gf.p1","pargf"] <- myprop(as.sociomatrix(gf_rvl[[1]]),as.sociomatrix(pargf_rvl[[1]]))
res_prop["outallf.p1","outgf"] <- myprop(as.sociomatrix(outallf_rvl[[1]]),as.sociomatrix(outgf_rvl[[1]]))
res_prop["parallf.p1","pargf"] <- myprop(as.sociomatrix(parallf_rvl[[1]]),as.sociomatrix(pargf_rvl[[1]]))

aa <- as.sociomatrix(myhome_jbs[[1]])+as.sociomatrix(hishome_jbs[[1]])+as.sociomatrix(outside_jbs[[1]])
bb <- symmetrize(as.sociomatrix(parents_jbs[[1]]),rule="weak")

res_prop["out.p2","outallf"] <- myprop(aa,as.sociomatrix(outallf_jbs[[1]]))
res_prop["out.p2","outgf"] <- myprop(aa,as.sociomatrix(outgf_jbs[[1]]))
res_prop["par.p2","parallf"] <- myprop(bb,as.sociomatrix(parallf_jbs[[1]]))
res_prop["par.p2","pargf"] <- myprop(bb,as.sociomatrix(pargf_jbs[[1]]))
res_prop["allf.p2","outallf"] <- myprop(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(outallf_jbs[[1]]))
res_prop["allf.p2","parallf"] <- myprop(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(parallf_jbs[[1]]))
res_prop["allf.p2","gf"] <- myprop(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(gf_jbs[[1]]))
res_prop["allf.p2","outgf"] <- myprop(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(outgf_jbs[[1]]))
res_prop["allf.p2","pargf"] <- myprop(as.sociomatrix(allf_jbs[[1]]),as.sociomatrix(pargf_jbs[[1]]))
res_prop["gf.p2","outgf"] <- myprop(as.sociomatrix(gf_jbs[[1]]),as.sociomatrix(outgf_jbs[[1]]))
res_prop["gf.p2","pargf"] <- myprop(as.sociomatrix(gf_jbs[[1]]),as.sociomatrix(pargf_jbs[[1]]))
res_prop["outallf.p2","outgf"] <- myprop(as.sociomatrix(outallf_jbs[[1]]),as.sociomatrix(outgf_jbs[[1]]))
res_prop["parallf.p2","pargf"] <- myprop(as.sociomatrix(parallf_jbs[[1]]),as.sociomatrix(pargf_jbs[[1]]))

aa <- as.sociomatrix(myhome_pm[[1]])+as.sociomatrix(hishome_pm[[1]])+as.sociomatrix(outside_pm[[1]])
bb <- symmetrize(as.sociomatrix(parents_pm[[1]]),rule="weak")

res_prop["out.s1","outallf"] <- myprop(aa,as.sociomatrix(outallf_pm[[1]]))
res_prop["out.s1","outgf"] <- myprop(aa,as.sociomatrix(outgf_pm[[1]]))
res_prop["par.s1","parallf"] <- myprop(bb,as.sociomatrix(parallf_pm[[1]]))
res_prop["par.s1","pargf"] <- myprop(bb,as.sociomatrix(pargf_pm[[1]]))
res_prop["allf.s1","outallf"] <- myprop(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(outallf_pm[[1]]))
res_prop["allf.s1","parallf"] <- myprop(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(parallf_pm[[1]]))
res_prop["allf.s1","gf"] <- myprop(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(gf_pm[[1]]))
res_prop["allf.s1","outgf"] <- myprop(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(outgf_pm[[1]]))
res_prop["allf.s1","pargf"] <- myprop(as.sociomatrix(allf_pm[[1]]),as.sociomatrix(pargf_pm[[1]]))
res_prop["gf.s1","outgf"] <- myprop(as.sociomatrix(gf_pm[[1]]),as.sociomatrix(outgf_pm[[1]]))
res_prop["gf.s1","pargf"] <- myprop(as.sociomatrix(gf_pm[[1]]),as.sociomatrix(pargf_pm[[1]]))
res_prop["outallf.s1","outgf"] <- myprop(as.sociomatrix(outallf_pm[[1]]),as.sociomatrix(outgf_pm[[1]]))
res_prop["parallf.s1","pargf"] <- myprop(as.sociomatrix(parallf_pm[[1]]),as.sociomatrix(pargf_pm[[1]]))

aa <- as.sociomatrix(myhome_dmz[[1]])+as.sociomatrix(hishome_dmz[[1]])+as.sociomatrix(outside_dmz[[1]])
bb <- symmetrize(as.sociomatrix(parents_dmz[[1]]),rule="weak")

res_prop["out.s2","outallf"] <- myprop(aa,as.sociomatrix(outallf_dmz[[1]]))
res_prop["out.s2","outgf"] <- myprop(aa,as.sociomatrix(outgf_dmz[[1]]))
res_prop["par.s2","parallf"] <- myprop(bb,as.sociomatrix(parallf_dmz[[1]]))
res_prop["par.s2","pargf"] <- myprop(bb,as.sociomatrix(pargf_dmz[[1]]))
res_prop["allf.s2","outallf"] <- myprop(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(outallf_dmz[[1]]))
res_prop["allf.s2","parallf"] <- myprop(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(parallf_dmz[[1]]))
res_prop["allf.s2","gf"] <- myprop(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(gf_dmz[[1]]))
res_prop["allf.s2","outgf"] <- myprop(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(outgf_dmz[[1]]))
res_prop["allf.s2","pargf"] <- myprop(as.sociomatrix(allf_dmz[[1]]),as.sociomatrix(pargf_dmz[[1]]))
res_prop["gf.s2","outgf"] <- myprop(as.sociomatrix(gf_dmz[[1]]),as.sociomatrix(outgf_dmz[[1]]))
res_prop["gf.s2","pargf"] <- myprop(as.sociomatrix(gf_dmz[[1]]),as.sociomatrix(pargf_dmz[[1]]))
res_prop["outallf.s2","outgf"] <- myprop(as.sociomatrix(outallf_dmz[[1]]),as.sociomatrix(outgf_dmz[[1]]))
res_prop["parallf.s2","pargf"] <- myprop(as.sociomatrix(parallf_dmz[[1]]),as.sociomatrix(pargf_dmz[[1]]))

res_prop # exemple : 25% des nominations de allf sont aussi des outallf a paris 1.

#
########## 1 - Create data frames ##########

net_rvl <- lapply(gf_rvl,off_NA,"factbrev")
net_jbs <- lapply(gf_jbs,off_NA,"factbrev")
net_pm <- lapply(gf_pm,off_NA,"factbrev")
net_dmz <- lapply(gf_dmz,off_NA,"factbrev")

net2_rvl <- lapply(allf_rvl,off_NA,"factbrev")
net2_jbs <- lapply(allf_jbs,off_NA,"factbrev")
net2_pm <- lapply(allf_pm,off_NA,"factbrev")
net2_dmz <- lapply(allf_dmz,off_NA,"factbrev")

net3_rvl <- lapply(outgf_rvl,off_NA,"factbrev")
net3_jbs <- lapply(outgf_jbs,off_NA,"factbrev")
net3_pm <- lapply(outgf_pm,off_NA,"factbrev")
net3_dmz <- lapply(outgf_dmz,off_NA,"factbrev")

net4_rvl <- lapply(pargf_rvl,off_NA,"factbrev")
net4_jbs <- lapply(pargf_jbs,off_NA,"factbrev")
net4_pm <- lapply(pargf_pm,off_NA,"factbrev")
net4_dmz <- lapply(pargf_dmz,off_NA,"factbrev")

net5_rvl <- lapply(dl_rvl,off_NA,"factbrev")
net5_jbs <- lapply(dl_jbs,off_NA,"factbrev")
net5_pm <- lapply(dl_pm,off_NA,"factbrev")
net5_dmz <- lapply(dl_dmz,off_NA,"factbrev")

net6_rvl <- lapply(bully_rvl,off_NA,"factbrev")
net6_jbs <- lapply(bully_jbs,off_NA,"factbrev")
net6_pm <- lapply(bully_pm,off_NA,"factbrev")
net6_dmz <- lapply(bully_dmz,off_NA,"factbrev")

###

mynet <- net_rvl[[1]]
mynet2 <- net2_rvl[[1]]
mynet3 <- net3_rvl[[1]]
mynet4 <- net4_rvl[[1]]
mynet5 <- net5_rvl[[1]]
mynet6 <- net6_rvl[[1]]
drvl_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_1) <- c("sender","receiver")
drvl_1$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_1$tie_outf <- as.vector(as.sociomatrix(mynet3))
drvl_1$tie_parf <- as.vector(as.sociomatrix(mynet4))
drvl_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_1[,i] <- as.vector(aa[[i]])
}
drvl_1 <- drvl_1[-which(drvl_1$sender==drvl_1$receiver),]
head(drvl_1)

mynet <- net_jbs[[1]]
mynet2 <- net2_jbs[[1]]
mynet3 <- net3_jbs[[1]]
mynet4 <- net4_jbs[[1]]
mynet5 <- net5_jbs[[1]]
mynet6 <- net6_jbs[[1]]
djbs_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_1) <- c("sender","receiver")
djbs_1$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_1$tie_outf <- as.vector(as.sociomatrix(mynet3))
djbs_1$tie_parf <- as.vector(as.sociomatrix(mynet4))
djbs_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_1[,i] <- as.vector(aa[[i]])
}
djbs_1 <- djbs_1[-which(djbs_1$sender==djbs_1$receiver),]
head(djbs_1)

mynet <- net_pm[[1]]
mynet2 <- net2_pm[[1]]
mynet3 <- net3_pm[[1]]
mynet4 <- net4_pm[[1]]
mynet5 <- net5_pm[[1]]
mynet6 <- net6_pm[[1]]
dpm_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_1) <- c("sender","receiver")
dpm_1$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_1$tie_outf <- as.vector(as.sociomatrix(mynet3))
dpm_1$tie_parf <- as.vector(as.sociomatrix(mynet4))
dpm_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_1[,i] <- as.vector(aa[[i]])
}
dpm_1 <- dpm_1[-which(dpm_1$sender==dpm_1$receiver),]
head(dpm_1)

mynet <- net_dmz[[1]]
mynet2 <- net2_dmz[[1]]
mynet3 <- net3_dmz[[1]]
mynet4 <- net4_dmz[[1]]
mynet5 <- net5_dmz[[1]]
mynet6 <- net6_dmz[[1]]
ddmz_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_1) <- c("sender","receiver")
ddmz_1$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_1$tie_outf <- as.vector(as.sociomatrix(mynet3))
ddmz_1$tie_parf <- as.vector(as.sociomatrix(mynet4))
ddmz_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_1[,i] <- as.vector(aa[[i]])
}
ddmz_1 <- ddmz_1[-which(ddmz_1$sender==ddmz_1$receiver),]
head(ddmz_1)


########## 2 - logistic regressions ##########

mynames <- c("factbrev_absdiff","factbrev_sender","factbrev_receiver")

m1_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m1_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m1_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m1_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m2_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m2_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m2_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m2_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m3_rvl <- glm(as.formula(paste("tie_outf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m3_jbs <- glm(as.formula(paste("tie_outf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m3_pm <- glm(as.formula(paste("tie_outf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m3_dmz <- glm(as.formula(paste("tie_outf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m4_rvl <- glm(as.formula(paste("tie_parf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m4_jbs <- glm(as.formula(paste("tie_parf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m4_pm <- glm(as.formula(paste("tie_parf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m4_dmz <- glm(as.formula(paste("tie_parf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m5_rvl <- glm(as.formula(paste("tie_dl ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m5_jbs <- glm(as.formula(paste("tie_dl ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m5_pm <- glm(as.formula(paste("tie_dl ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m5_dmz <- glm(as.formula(paste("tie_dl ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m6_rvl <- glm(as.formula(paste("tie_bully ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m6_jbs <- glm(as.formula(paste("tie_bully ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m6_pm <- glm(as.formula(paste("tie_bully ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m6_dmz <- glm(as.formula(paste("tie_bully ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

stargazer(m1_rvl,m2_rvl,m3_rvl,m4_rvl,m5_rvl,m6_rvl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_jbs,m2_jbs,m3_jbs,m4_jbs,m5_jbs,m6_jbs,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_pm,m2_pm,m3_pm,m4_pm,m5_pm,m6_pm,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_dmz,m2_dmz,m3_dmz,m4_dmz,m5_dmz,m6_dmz,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

########## 3 - clustered standard errors ##########

library(sandwich)
library(lmtest)

m1_rvl_cl <- coeftest(m1_rvl, vcov = vcovCL, cluster = ~sender)
m2_rvl_cl <- coeftest(m2_rvl, vcov = vcovCL, cluster = ~sender)
m3_rvl_cl <- coeftest(m3_rvl, vcov = vcovCL, cluster = ~sender)
m4_rvl_cl <- coeftest(m4_rvl, vcov = vcovCL, cluster = ~sender)
m5_rvl_cl <- coeftest(m5_rvl, vcov = vcovCL, cluster = ~sender)
m6_rvl_cl <- coeftest(m6_rvl, vcov = vcovCL, cluster = ~sender)

m1_jbs_cl <- coeftest(m1_jbs, vcov = vcovCL, cluster = ~sender)
m2_jbs_cl <- coeftest(m2_jbs, vcov = vcovCL, cluster = ~sender)
m3_jbs_cl <- coeftest(m3_jbs, vcov = vcovCL, cluster = ~sender)
m4_jbs_cl <- coeftest(m4_jbs, vcov = vcovCL, cluster = ~sender)
m5_jbs_cl <- coeftest(m5_jbs, vcov = vcovCL, cluster = ~sender)
m6_jbs_cl <- coeftest(m6_jbs, vcov = vcovCL, cluster = ~sender)

m1_pm_cl <- coeftest(m1_pm, vcov = vcovCL, cluster = ~sender)
m2_pm_cl <- coeftest(m2_pm, vcov = vcovCL, cluster = ~sender)
m3_pm_cl <- coeftest(m3_pm, vcov = vcovCL, cluster = ~sender)
m4_pm_cl <- coeftest(m4_pm, vcov = vcovCL, cluster = ~sender)
m5_pm_cl <- coeftest(m5_pm, vcov = vcovCL, cluster = ~sender)
m6_pm_cl <- coeftest(m6_pm, vcov = vcovCL, cluster = ~sender)

m1_dmz_cl <- coeftest(m1_dmz, vcov = vcovCL, cluster = ~sender)
m2_dmz_cl <- coeftest(m2_dmz, vcov = vcovCL, cluster = ~sender)
m3_dmz_cl <- coeftest(m3_dmz, vcov = vcovCL, cluster = ~sender)
m4_dmz_cl <- coeftest(m4_dmz, vcov = vcovCL, cluster = ~sender)
m5_dmz_cl <- coeftest(m5_dmz, vcov = vcovCL, cluster = ~sender)
m6_dmz_cl <- coeftest(m6_dmz, vcov = vcovCL, cluster = ~sender)

stargazer(m1_rvl_cl,m2_rvl_cl,m3_rvl_cl,m4_rvl_cl,m5_rvl_cl,m6_rvl_cl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_jbs_cl,m2_jbs_cl,m3_jbs_cl,m4_jbs_cl,m5_jbs_cl,m6_jbs_cl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_pm_cl,m2_pm_cl,m3_pm_cl,m4_pm_cl,m5_pm_cl,m6_pm_cl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

stargazer(m1_dmz_cl,m2_dmz_cl,m3_dmz_cl,m4_dmz_cl,m5_dmz_cl,m6_dmz_cl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

##

# stargazer(m1_rvl_cl,m2_rvl_cl,m3_rvl_cl,m4_rvl_cl,m5_rvl_cl,m6_rvl_cl,type="html",out="a1.html",style = "ajps",
#           column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)
# 
# stargazer(m1_jbs_cl,m2_jbs_cl,m3_jbs_cl,m4_jbs_cl,m5_jbs_cl,m6_jbs_cl,type="html",out="a2.html",
#           column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)
# 
# stargazer(m1_pm_cl,m2_pm_cl,m3_pm_cl,m4_pm_cl,m5_pm_cl,m6_pm_cl,type="html",out="a3.html",
#           column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)
# 
# stargazer(m1_dmz_cl,m2_dmz_cl,m3_dmz_cl,m4_dmz_cl,m5_dmz_cl,m6_dmz_cl,type="html",out="a4.html",
#           column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)


########## 4 - test de robustesse - AME ##########

# Pour verifier que le biais lie a la comparaison des coefficients log odd entre modeles ne fausse pas les resultats,
# on calcule des AME (non sensibles a ce biais, cf. Mood 2010)

library(margins)
ame1_rvl <- margins(m1_rvl)
ame2_rvl <- margins(m2_rvl)
ame3_rvl <- margins(m3_rvl)
ame4_rvl <- margins(m4_rvl)
ame5_rvl <- margins(m5_rvl)
ame6_rvl <- margins(m6_rvl)

ame1_jbs <- margins(m1_jbs)
ame2_jbs <- margins(m2_jbs)
ame3_jbs <- margins(m3_jbs)
ame4_jbs <- margins(m4_jbs)
ame5_jbs <- margins(m5_jbs)
ame6_jbs <- margins(m6_jbs)

ame1_pm <- margins(m1_pm)
ame2_pm <- margins(m2_pm)
ame3_pm <- margins(m3_pm)
ame4_pm <- margins(m4_pm)
ame5_pm <- margins(m5_pm)
ame6_pm <- margins(m6_pm)

ame1_dmz <- margins(m1_dmz)
ame2_dmz <- margins(m2_dmz)
ame3_dmz <- margins(m3_dmz)
ame4_dmz <- margins(m4_dmz)
ame5_dmz <- margins(m5_dmz)
ame6_dmz <- margins(m6_dmz)

summary(ame1_rvl)
summary(ame2_rvl)
summary(ame3_rvl)
summary(ame4_rvl)
summary(ame5_rvl)
summary(ame6_rvl)

stargazer(m1_rvl,m2_rvl,m3_rvl,m4_rvl,m5_rvl,m6_rvl,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

summary(ame1_jbs)
summary(ame2_jbs)
summary(ame3_jbs)
summary(ame4_jbs)
summary(ame5_jbs)
summary(ame6_jbs)

stargazer(m1_jbs,m2_jbs,m3_jbs,m4_jbs,m5_jbs,m6_jbs,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

summary(ame1_pm)
summary(ame2_pm)
summary(ame3_pm)
summary(ame4_pm)
summary(ame5_pm)
summary(ame6_pm)

stargazer(m1_pm,m2_pm,m3_pm,m4_pm,m5_pm,m6_pm,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

summary(ame1_dmz)
summary(ame2_dmz)
summary(ame3_dmz)
summary(ame4_dmz)
summary(ame5_dmz)
summary(ame6_dmz)

stargazer(m1_dmz,m2_dmz,m3_dmz,m4_dmz,m5_dmz,m6_dmz,type="text",
          column.labels = c("allf","gf","outf","parf","dl","bully"),dep.var.labels.include=F)

########## 5 - Exposure ##########

# Quick reminder: what is the average degree per type of tie?
ind1$outdeg_allf <- c(degree(allf_rvl[[1]],cmode="outdegree"),degree(allf_jbs[[1]],cmode="outdegree"),
                      degree(allf_pm[[1]],cmode="outdegree"),degree(allf_dmz[[1]],cmode="outdegree"))

ind1$outdeg_gf <- c(degree(gf_rvl[[1]],cmode="outdegree"),degree(gf_jbs[[1]],cmode="outdegree"),
                    degree(gf_pm[[1]],cmode="outdegree"),degree(gf_dmz[[1]],cmode="outdegree"))

ind1$outdeg_outgf <- c(degree(outgf_rvl[[1]],cmode="outdegree"),degree(outgf_jbs[[1]],cmode="outdegree"),
                       degree(outgf_pm[[1]],cmode="outdegree"),degree(outgf_dmz[[1]],cmode="outdegree"))

ind1$outdeg_pargf <- c(degree(pargf_rvl[[1]],cmode="outdegree"),degree(pargf_jbs[[1]],cmode="outdegree"),
                       degree(pargf_pm[[1]],cmode="outdegree"),degree(pargf_dmz[[1]],cmode="outdegree"))

summary(ind1$outdeg_allf)
summary(ind1$outdeg_gf)
summary(ind1$outdeg_outgf)
summary(ind1$outdeg_pargf)


# Here I want to look at quartiles of factbrev, and ask : how many kids in the lowest quartile have at least
# a friend in the highest one, and vice versa ?

# We use sended tie for this (how many people one named, regardless of reciprocation)
# We create variables in ind data frames.

# First, create the factbrev quantiles, per school.

levels(quant.cut(rvl1$factbrev,4))
levels(quant.cut(jbs1$factbrev,4))
levels(quant.cut(pm1$factbrev,4))
levels(quant.cut(dmz1$factbrev,4))

rvl1$factbrev_quartiles <- as.character(quant.cut(rvl1$factbrev,4,labels = c("inf","low","up","sup")))
jbs1$factbrev_quartiles <- as.character(quant.cut(jbs1$factbrev,4,labels = c("inf","low","up","sup")))
pm1$factbrev_quartiles <- as.character(quant.cut(pm1$factbrev,4,labels = c("inf","low","up","sup")))
dmz1$factbrev_quartiles <- as.character(quant.cut(dmz1$factbrev,4,labels = c("inf","low","up","sup")))

# And back to ind1:
ind1$factbrev_quartiles <- c(rvl1$factbrev_quartiles,jbs1$factbrev_quartiles,pm1$factbrev_quartiles,dmz1$factbrev_quartiles)
table(ind1$factbrev_quartiles,useNA="ifany")

ind1$nb_allf_inf_factbrev <- NA
ind1$nb_gf_inf_factbrev <- NA
ind1$nb_outgf_inf_factbrev <- NA
ind1$nb_pargf_inf_factbrev <- NA

ind1$nb_allf_sup_factbrev <- NA
ind1$nb_gf_sup_factbrev <- NA
ind1$nb_outgf_sup_factbrev <- NA
ind1$nb_pargf_sup_factbrev <- NA

for(aa in 1:nrow(ind1)){
  
  id_ego <- ind1$id[aa]
  
  if(strsplit(id_ego,"")[[1]][2]=="R"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_rvl[[1]])[allf_rvl[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_rvl[[1]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_rvl[[1]])[gf_rvl[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_rvl[[1]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(outgf_rvl[[1]])[outgf_rvl[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- outgf_rvl[[1]]%v%"vertex.names"
    outgf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(pargf_rvl[[1]])[pargf_rvl[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- pargf_rvl[[1]]%v%"vertex.names"
    pargf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%allf_ego]
    gf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%gf_ego]
    outgf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%outgf_ego]
    pargf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%pargf_ego]
    
    # Count number of inf and sup (nb : NAs are not counted, thanks to how "which" works by default).
    ind1$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind1$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    ind1$nb_outgf_inf_factbrev[aa] <- length(which(outgf_ego_quartiles=="inf"))
    ind1$nb_pargf_inf_factbrev[aa] <- length(which(pargf_ego_quartiles=="inf"))
    
    ind1$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind1$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
    ind1$nb_outgf_sup_factbrev[aa] <-length(which(outgf_ego_quartiles=="sup"))
    ind1$nb_pargf_sup_factbrev[aa] <-length(which(pargf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="J"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_jbs[[1]])[allf_jbs[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_jbs[[1]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_jbs[[1]])[gf_jbs[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_jbs[[1]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(outgf_jbs[[1]])[outgf_jbs[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- outgf_jbs[[1]]%v%"vertex.names"
    outgf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(pargf_jbs[[1]])[pargf_jbs[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- pargf_jbs[[1]]%v%"vertex.names"
    pargf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%allf_ego]
    gf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%gf_ego]
    outgf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%outgf_ego]
    pargf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%pargf_ego]
    
    # Count number of inf and sup.
    ind1$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind1$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    ind1$nb_outgf_inf_factbrev[aa] <- length(which(outgf_ego_quartiles=="inf"))
    ind1$nb_pargf_inf_factbrev[aa] <- length(which(pargf_ego_quartiles=="inf"))
    
    ind1$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind1$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
    ind1$nb_outgf_sup_factbrev[aa] <-length(which(outgf_ego_quartiles=="sup"))
    ind1$nb_pargf_sup_factbrev[aa] <-length(which(pargf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="P"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_pm[[1]])[allf_pm[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_pm[[1]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_pm[[1]])[gf_pm[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_pm[[1]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(outgf_pm[[1]])[outgf_pm[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- outgf_pm[[1]]%v%"vertex.names"
    outgf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(pargf_pm[[1]])[pargf_pm[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- pargf_pm[[1]]%v%"vertex.names"
    pargf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%allf_ego]
    gf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%gf_ego]
    outgf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%outgf_ego]
    pargf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%pargf_ego]
    
    # Count number of inf and sup.
    ind1$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind1$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    ind1$nb_outgf_inf_factbrev[aa] <- length(which(outgf_ego_quartiles=="inf"))
    ind1$nb_pargf_inf_factbrev[aa] <- length(which(pargf_ego_quartiles=="inf"))
    
    ind1$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind1$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
    ind1$nb_outgf_sup_factbrev[aa] <-length(which(outgf_ego_quartiles=="sup"))
    ind1$nb_pargf_sup_factbrev[aa] <-length(which(pargf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="D"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_dmz[[1]])[allf_dmz[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_dmz[[1]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_dmz[[1]])[gf_dmz[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_dmz[[1]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(outgf_dmz[[1]])[outgf_dmz[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- outgf_dmz[[1]]%v%"vertex.names"
    outgf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(pargf_dmz[[1]])[pargf_dmz[[1]]%v%"vertex.names"==id_ego,]
    names(a) <- pargf_dmz[[1]]%v%"vertex.names"
    pargf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%allf_ego]
    gf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%gf_ego]
    outgf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%outgf_ego]
    pargf_ego_quartiles <- ind1$factbrev_quartiles[ind1$id%in%pargf_ego]
    
    # Count number of inf and sup.
    ind1$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind1$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    ind1$nb_outgf_inf_factbrev[aa] <- length(which(outgf_ego_quartiles=="inf"))
    ind1$nb_pargf_inf_factbrev[aa] <- length(which(pargf_ego_quartiles=="inf"))
    
    ind1$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind1$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
    ind1$nb_outgf_sup_factbrev[aa] <-length(which(outgf_ego_quartiles=="sup"))
    ind1$nb_pargf_sup_factbrev[aa] <-length(which(pargf_ego_quartiles=="sup"))
  }
  
}

# So ! How many kids from the last quartile have at least one (symmetric) friend in the high quantile ?
round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)

# The opposite :
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)

##

# To better see homophily, we can also compare with intra-quartile nominations:

round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)
round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="sup"])),digits=2)

# The opposite :
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="inf"])),digits=2)

############################################################
############### 3 - Evolution dans le temps ################
############################################################
########## 1 - create data frames ##########

net_rvl <- lapply(gf_rvl,off_NA,"factbrev") # see custom function "off_NA" in the setup
net_jbs <- lapply(gf_jbs,off_NA,"factbrev")
net_pm <- lapply(gf_pm,off_NA,"factbrev")
net_dmz <- lapply(gf_dmz,off_NA,"factbrev")

net2_rvl <- lapply(allf_rvl,off_NA,"factbrev")
net2_jbs <- lapply(allf_jbs,off_NA,"factbrev")
net2_pm <- lapply(allf_pm,off_NA,"factbrev")
net2_dmz <- lapply(allf_dmz,off_NA,"factbrev")

net5_rvl <- lapply(dl_rvl,off_NA,"factbrev")
net5_jbs <- lapply(dl_jbs,off_NA,"factbrev")
net5_pm <- lapply(dl_pm,off_NA,"factbrev")
net5_dmz <- lapply(dl_dmz,off_NA,"factbrev")

net6_rvl <- lapply(bully_rvl,off_NA,"factbrev")
net6_jbs <- lapply(bully_jbs,off_NA,"factbrev")
net6_pm <- lapply(bully_pm,off_NA,"factbrev")
net6_dmz <- lapply(bully_dmz,off_NA,"factbrev")

###

mynet <- net_rvl[[1]]
mynet2 <- net2_rvl[[1]]
mynet5 <- net5_rvl[[1]]
mynet6 <- net6_rvl[[1]]
drvl_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_1) <- c("sender","receiver")
drvl_1$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_1[,i] <- as.vector(aa[[i]])
}
drvl_1 <- drvl_1[-which(drvl_1$sender==drvl_1$receiver),]
head(drvl_1)

mynet <- net_rvl[[2]]
mynet2 <- net2_rvl[[2]]
mynet5 <- net5_rvl[[2]]
mynet6 <- net6_rvl[[2]]
drvl_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_2) <- c("sender","receiver")
drvl_2$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_2$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_2$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_2[,i] <- as.vector(aa[[i]])
}
drvl_2 <- drvl_2[-which(drvl_2$sender==drvl_2$receiver),]
head(drvl_2)

mynet <- net_rvl[[3]]
mynet2 <- net2_rvl[[3]]
mynet5 <- net5_rvl[[3]]
mynet6 <- net6_rvl[[3]]
drvl_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_3) <- c("sender","receiver")
drvl_3$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_3$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_3$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_3[,i] <- as.vector(aa[[i]])
}
drvl_3 <- drvl_3[-which(drvl_3$sender==drvl_3$receiver),]
head(drvl_3)

mynet <- net_rvl[[4]]
mynet2 <- net2_rvl[[4]]
mynet5 <- net5_rvl[[4]]
mynet6 <- net6_rvl[[4]]
drvl_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_4) <- c("sender","receiver")
drvl_4$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_4$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_4$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_4[,i] <- as.vector(aa[[i]])
}
drvl_4 <- drvl_4[-which(drvl_4$sender==drvl_4$receiver),]
head(drvl_4)

mynet <- net_rvl[[5]]
mynet2 <- net2_rvl[[5]]
mynet5 <- net5_rvl[[5]]
mynet6 <- net6_rvl[[5]]
drvl_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_5) <- c("sender","receiver")
drvl_5$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_5$tie_dl <- as.vector(as.sociomatrix(mynet5))
drvl_5$tie_bully <- as.vector(as.sociomatrix(mynet6))
drvl_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_5[,i] <- as.vector(aa[[i]])
}
drvl_5 <- drvl_5[-which(drvl_5$sender==drvl_5$receiver),]
head(drvl_5)

###

mynet <- net_jbs[[1]]
mynet2 <- net2_jbs[[1]]
mynet5 <- net5_jbs[[1]]
mynet6 <- net6_jbs[[1]]
djbs_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_1) <- c("sender","receiver")
djbs_1$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_1[,i] <- as.vector(aa[[i]])
}
djbs_1 <- djbs_1[-which(djbs_1$sender==djbs_1$receiver),]
head(djbs_1)

mynet <- net_jbs[[2]]
mynet2 <- net2_jbs[[2]]
mynet5 <- net5_jbs[[2]]
mynet6 <- net6_jbs[[2]]
djbs_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_2) <- c("sender","receiver")
djbs_2$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_2$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_2$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
# djbs_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
# djbs_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
# djbs_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_2[,i] <- as.vector(aa[[i]])
}
djbs_2 <- djbs_2[-which(djbs_2$sender==djbs_2$receiver),]
head(djbs_2)

mynet <- net_jbs[[3]]
mynet2 <- net2_jbs[[3]]
mynet5 <- net5_jbs[[3]]
mynet6 <- net6_jbs[[3]]
djbs_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_3) <- c("sender","receiver")
djbs_3$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_3$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_3$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
# djbs_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
# djbs_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
# djbs_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_3[,i] <- as.vector(aa[[i]])
}
djbs_3 <- djbs_3[-which(djbs_3$sender==djbs_3$receiver),]
head(djbs_3)

mynet <- net_jbs[[4]]
mynet2 <- net2_jbs[[4]]
mynet5 <- net5_jbs[[4]]
mynet6 <- net6_jbs[[4]]
djbs_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_4) <- c("sender","receiver")
djbs_4$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_4$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_4$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
# djbs_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
# djbs_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
# djbs_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_4[,i] <- as.vector(aa[[i]])
}
djbs_4 <- djbs_4[-which(djbs_4$sender==djbs_4$receiver),]
head(djbs_4)

mynet <- net_jbs[[5]]
mynet2 <- net2_jbs[[5]]
mynet5 <- net5_jbs[[5]]
mynet6 <- net6_jbs[[5]]
djbs_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_5) <- c("sender","receiver")
djbs_5$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_5$tie_dl <- as.vector(as.sociomatrix(mynet5))
djbs_5$tie_bully <- as.vector(as.sociomatrix(mynet6))
djbs_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
# djbs_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
# djbs_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
# djbs_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_5[,i] <- as.vector(aa[[i]])
}
djbs_5 <- djbs_5[-which(djbs_5$sender==djbs_5$receiver),]
head(djbs_5)

###

mynet <- net_pm[[1]]
mynet2 <- net2_pm[[1]]
mynet5 <- net5_pm[[1]]
mynet6 <- net6_pm[[1]]
dpm_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_1) <- c("sender","receiver")
dpm_1$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_1[,i] <- as.vector(aa[[i]])
}
dpm_1 <- dpm_1[-which(dpm_1$sender==dpm_1$receiver),]
head(dpm_1)

mynet <- net_pm[[2]]
mynet2 <- net2_pm[[2]]
mynet5 <- net5_pm[[2]]
mynet6 <- net6_pm[[2]]
dpm_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_2) <- c("sender","receiver")
dpm_2$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_2$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_2$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_2[,i] <- as.vector(aa[[i]])
}
dpm_2 <- dpm_2[-which(dpm_2$sender==dpm_2$receiver),]
head(dpm_2)

mynet <- net_pm[[3]]
mynet2 <- net2_pm[[3]]
mynet5 <- net5_pm[[3]]
mynet6 <- net6_pm[[3]]
dpm_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_3) <- c("sender","receiver")
dpm_3$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_3$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_3$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_3[,i] <- as.vector(aa[[i]])
}
dpm_3 <- dpm_3[-which(dpm_3$sender==dpm_3$receiver),]
head(dpm_3)

mynet <- net_pm[[4]]
mynet2 <- net2_pm[[4]]
mynet5 <- net5_pm[[4]]
mynet6 <- net6_pm[[4]]
dpm_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_4) <- c("sender","receiver")
dpm_4$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_4$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_4$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_4[,i] <- as.vector(aa[[i]])
}
dpm_4 <- dpm_4[-which(dpm_4$sender==dpm_4$receiver),]
head(dpm_4)

mynet <- net_pm[[5]]
mynet2 <- net2_pm[[5]]
mynet5 <- net5_pm[[5]]
mynet6 <- net6_pm[[5]]
dpm_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_5) <- c("sender","receiver")
dpm_5$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_5$tie_dl <- as.vector(as.sociomatrix(mynet5))
dpm_5$tie_bully <- as.vector(as.sociomatrix(mynet6))
dpm_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_5[,i] <- as.vector(aa[[i]])
}
dpm_5 <- dpm_5[-which(dpm_5$sender==dpm_5$receiver),]
head(dpm_5)

###

mynet <- net_dmz[[1]]
mynet2 <- net2_dmz[[1]]
mynet5 <- net5_dmz[[1]]
mynet6 <- net6_dmz[[1]]
ddmz_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_1) <- c("sender","receiver")
ddmz_1$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_1$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_1$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_1[,i] <- as.vector(aa[[i]])
}
ddmz_1 <- ddmz_1[-which(ddmz_1$sender==ddmz_1$receiver),]
head(ddmz_1)

mynet <- net_dmz[[2]]
mynet2 <- net2_dmz[[2]]
mynet5 <- net5_dmz[[2]]
mynet6 <- net6_dmz[[2]]
ddmz_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_2) <- c("sender","receiver")
ddmz_2$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_2$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_2$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_2[,i] <- as.vector(aa[[i]])
}
ddmz_2 <- ddmz_2[-which(ddmz_2$sender==ddmz_2$receiver),]
head(ddmz_2)

mynet <- net_dmz[[3]]
mynet2 <- net2_dmz[[3]]
mynet5 <- net5_dmz[[3]]
mynet6 <- net6_dmz[[3]]
ddmz_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_3) <- c("sender","receiver")
ddmz_3$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_3$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_3$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_3[,i] <- as.vector(aa[[i]])
}
ddmz_3 <- ddmz_3[-which(ddmz_3$sender==ddmz_3$receiver),]
head(ddmz_3)

mynet <- net_dmz[[4]]
mynet2 <- net2_dmz[[4]]
mynet5 <- net5_dmz[[4]]
mynet6 <- net6_dmz[[4]]
ddmz_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_4) <- c("sender","receiver")
ddmz_4$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_4$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_4$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_4[,i] <- as.vector(aa[[i]])
}
ddmz_4 <- ddmz_4[-which(ddmz_4$sender==ddmz_4$receiver),]
head(ddmz_4)

mynet <- net_dmz[[5]]
mynet2 <- net2_dmz[[5]]
mynet5 <- net5_dmz[[5]]
mynet6 <- net6_dmz[[5]]
ddmz_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_5) <- c("sender","receiver")
ddmz_5$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_5$tie_dl <- as.vector(as.sociomatrix(mynet5))
ddmz_5$tie_bully <- as.vector(as.sociomatrix(mynet6))
ddmz_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_5[,i] <- as.vector(aa[[i]])
}
ddmz_5 <- ddmz_5[-which(ddmz_5$sender==ddmz_5$receiver),]
head(ddmz_5)

########## 2 - logistic regressions + clustered s.e. - gf ##########

mynames <- c("factbrev_absdiff","factbrev_sender","factbrev_receiver")

m1_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m1_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m1_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m1_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m2_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_2)),collapse="+"),collapse = " ")),
              data = drvl_2,family="binomial")
m2_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_2)),collapse="+"),collapse = " ")),
              data = djbs_2,family="binomial")
m2_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_2)),collapse="+"),collapse = " ")),
             data = dpm_2,family="binomial")
m2_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_2)),collapse="+"),collapse = " ")),
              data = ddmz_2,family="binomial")

m3_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_3)),collapse="+"),collapse = " ")),
              data = drvl_3,family="binomial")
m3_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_3)),collapse="+"),collapse = " ")),
              data = djbs_3,family="binomial")
m3_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_3)),collapse="+"),collapse = " ")),
             data = dpm_3,family="binomial")
m3_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_3)),collapse="+"),collapse = " ")),
              data = ddmz_3,family="binomial")

m4_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_4)),collapse="+"),collapse = " ")),
              data = drvl_4,family="binomial")
m4_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_4)),collapse="+"),collapse = " ")),
              data = djbs_4,family="binomial")
m4_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_4)),collapse="+"),collapse = " ")),
             data = dpm_4,family="binomial")
m4_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_4)),collapse="+"),collapse = " ")),
              data = ddmz_4,family="binomial")

m5_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_5)),collapse="+"),collapse = " ")),
              data = drvl_5,family="binomial")
m5_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_5)),collapse="+"),collapse = " ")),
              data = djbs_5,family="binomial")
m5_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_5)),collapse="+"),collapse = " ")),
             data = dpm_5,family="binomial")
m5_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_5)),collapse="+"),collapse = " ")),
              data = ddmz_5,family="binomial")

library(sandwich)
library(lmtest)

m1_rvl_cl <- coeftest(m1_rvl, vcov = vcovCL, cluster = ~sender)
m2_rvl_cl <- coeftest(m2_rvl, vcov = vcovCL, cluster = ~sender)
m3_rvl_cl <- coeftest(m3_rvl, vcov = vcovCL, cluster = ~sender)
m4_rvl_cl <- coeftest(m4_rvl, vcov = vcovCL, cluster = ~sender)
m5_rvl_cl <- coeftest(m5_rvl, vcov = vcovCL, cluster = ~sender)

m1_jbs_cl <- coeftest(m1_jbs, vcov = vcovCL, cluster = ~sender)
m2_jbs_cl <- coeftest(m2_jbs, vcov = vcovCL, cluster = ~sender)
m3_jbs_cl <- coeftest(m3_jbs, vcov = vcovCL, cluster = ~sender)
m4_jbs_cl <- coeftest(m4_jbs, vcov = vcovCL, cluster = ~sender)
m5_jbs_cl <- coeftest(m5_jbs, vcov = vcovCL, cluster = ~sender)

m1_pm_cl <- coeftest(m1_pm, vcov = vcovCL, cluster = ~sender)
m2_pm_cl <- coeftest(m2_pm, vcov = vcovCL, cluster = ~sender)
m3_pm_cl <- coeftest(m3_pm, vcov = vcovCL, cluster = ~sender)
m4_pm_cl <- coeftest(m4_pm, vcov = vcovCL, cluster = ~sender)
m5_pm_cl <- coeftest(m5_pm, vcov = vcovCL, cluster = ~sender)

m1_dmz_cl <- coeftest(m1_dmz, vcov = vcovCL, cluster = ~sender)
m2_dmz_cl <- coeftest(m2_dmz, vcov = vcovCL, cluster = ~sender)
m3_dmz_cl <- coeftest(m3_dmz, vcov = vcovCL, cluster = ~sender)
m4_dmz_cl <- coeftest(m4_dmz, vcov = vcovCL, cluster = ~sender)
m5_dmz_cl <- coeftest(m5_dmz, vcov = vcovCL, cluster = ~sender)

########## 3 - plots - gf ##########

res <- matrix(nrow = 20, ncol = 4)

res[1,1] <- m1_rvl_cl[2,1]
res[1,2] <- m1_rvl_cl[2,2]
res[1,3] <- "Paris 1"
res[1,4] <- "w1"
res[2,1] <- m2_rvl_cl[2,1]
res[2,2] <- m2_rvl_cl[2,2]
res[2,3] <- "Paris 1"
res[2,4] <- "w2"
res[3,1] <- m3_rvl_cl[2,1]
res[3,2] <- m3_rvl_cl[2,2]
res[3,3] <- "Paris 1"
res[3,4] <- "w3"
res[4,1] <- m4_rvl_cl[2,1]
res[4,2] <- m4_rvl_cl[2,2]
res[4,3] <- "Paris 1"
res[4,4] <- "w4"
res[5,1] <- m5_rvl_cl[2,1]
res[5,2] <- m5_rvl_cl[2,2]
res[5,3] <- "Paris 1"
res[5,4] <- "w6"

res[1+5,1] <- m1_jbs_cl[2,1]
res[1+5,2] <- m1_jbs_cl[2,2]
res[1+5,3] <- "Paris 2"
res[1+5,4] <- "w1"
res[2+5,1] <- m2_jbs_cl[2,1]
res[2+5,2] <- m2_jbs_cl[2,2]
res[2+5,3] <- "Paris 2"
res[2+5,4] <- "w2"
res[3+5,1] <- m3_jbs_cl[2,1]
res[3+5,2] <- m3_jbs_cl[2,2]
res[3+5,3] <- "Paris 2"
res[3+5,4] <- "w3"
res[4+5,1] <- m4_jbs_cl[2,1]
res[4+5,2] <- m4_jbs_cl[2,2]
res[4+5,3] <- "Paris 2"
res[4+5,4] <- "w4"
res[5+5,1] <- m5_jbs_cl[2,1]
res[5+5,2] <- m5_jbs_cl[2,2]
res[5+5,3] <- "Paris 2"
res[5+5,4] <- "w6"

res[1+10,1] <- m1_pm_cl[2,1]
res[1+10,2] <- m1_pm_cl[2,2]
res[1+10,3] <- "Savoie 1"
res[1+10,4] <- "w1"
res[2+10,1] <- m2_pm_cl[2,1]
res[2+10,2] <- m2_pm_cl[2,2]
res[2+10,3] <- "Savoie 1"
res[2+10,4] <- "w2"
res[3+10,1] <- m3_pm_cl[2,1]
res[3+10,2] <- m3_pm_cl[2,2]
res[3+10,3] <- "Savoie 1"
res[3+10,4] <- "w3"
res[4+10,1] <- m4_pm_cl[2,1]
res[4+10,2] <- m4_pm_cl[2,2]
res[4+10,3] <- "Savoie 1"
res[4+10,4] <- "w4"
res[5+10,1] <- m5_pm_cl[2,1]
res[5+10,2] <- m5_pm_cl[2,2]
res[5+10,3] <- "Savoie 1"
res[5+10,4] <- "w6"

res[1+15,1] <- m1_dmz_cl[2,1]
res[1+15,2] <- m1_dmz_cl[2,2]
res[1+15,3] <- "Savoie 2"
res[1+15,4] <- "w1"
res[2+15,1] <- m2_dmz_cl[2,1]
res[2+15,2] <- m2_dmz_cl[2,2]
res[2+15,3] <- "Savoie 2"
res[2+15,4] <- "w2"
res[3+15,1] <- m3_dmz_cl[2,1]
res[3+15,2] <- m3_dmz_cl[2,2]
res[3+15,3] <- "Savoie 2"
res[3+15,4] <- "w3"
res[4+15,1] <- m4_dmz_cl[2,1]
res[4+15,2] <- m4_dmz_cl[2,2]
res[4+15,3] <- "Savoie 2"
res[4+15,4] <- "w4"
res[5+15,1] <- m5_dmz_cl[2,1]
res[5+15,2] <- m5_dmz_cl[2,2]
res[5+15,3] <- "Savoie 2"
res[5+15,4] <- "w6"

res <- data.frame(res,stringsAsFactors = F)
names(res) <- c("coef","sd","school","wave")
res$coef <- round(as.numeric(res$coef), digits = 4)
res$sd <- round(as.numeric(res$sd), digits = 4)
res$lower_bound <- res$coef - res$sd*2
res$upper_bound <- res$coef + res$sd*2

ggplot(res,aes(x = wave,y = coef, group = school, colour = school)) +
  geom_point(aes(shape=school, color=school), size = 8) +
  #geom_ribbon(aes(ymin=lower_bound,ymax=upper_bound),alpha=0.1,linetype="solid",lwd=1) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = .2) +
  geom_line(linetype="solid",lwd=1.2) +
  geom_hline(yintercept=0, color = "grey50", size=1) +
  scale_y_reverse() +
  scale_x_discrete(limits=c("w1","w2","w3","w4","w5 (missing)","w6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n (5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Collège",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  scale_shape_manual(name = "Collège",
                     values = c(15,16,17,18),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  ggtitle("Coefficients log-odd d\'homophilie sociale, par collège et vague \n Très bons amis") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24)) +
  labs(x="vague",color="Collège")

##

ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement variance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Paris 1 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

########## 4 - logistic regressions + clustered s.e. - allf ##########

mynames <- c("factbrev_absdiff","factbrev_sender","factbrev_receiver")

m1_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m1_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m1_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m1_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")

m2_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_2)),collapse="+"),collapse = " ")),
              data = drvl_2,family="binomial")
m2_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_2)),collapse="+"),collapse = " ")),
              data = djbs_2,family="binomial")
m2_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_2)),collapse="+"),collapse = " ")),
             data = dpm_2,family="binomial")
m2_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_2)),collapse="+"),collapse = " ")),
              data = ddmz_2,family="binomial")

m3_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_3)),collapse="+"),collapse = " ")),
              data = drvl_3,family="binomial")
m3_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_3)),collapse="+"),collapse = " ")),
              data = djbs_3,family="binomial")
m3_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_3)),collapse="+"),collapse = " ")),
             data = dpm_3,family="binomial")
m3_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_3)),collapse="+"),collapse = " ")),
              data = ddmz_3,family="binomial")

m4_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_4)),collapse="+"),collapse = " ")),
              data = drvl_4,family="binomial")
m4_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_4)),collapse="+"),collapse = " ")),
              data = djbs_4,family="binomial")
m4_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_4)),collapse="+"),collapse = " ")),
             data = dpm_4,family="binomial")
m4_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_4)),collapse="+"),collapse = " ")),
              data = ddmz_4,family="binomial")

m5_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_5)),collapse="+"),collapse = " ")),
              data = drvl_5,family="binomial")
m5_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_5)),collapse="+"),collapse = " ")),
              data = djbs_5,family="binomial")
m5_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_5)),collapse="+"),collapse = " ")),
             data = dpm_5,family="binomial")
m5_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_5)),collapse="+"),collapse = " ")),
              data = ddmz_5,family="binomial")

library(sandwich)
library(lmtest)

m1_rvl_cl <- coeftest(m1_rvl, vcov = vcovCL, cluster = ~sender)
m2_rvl_cl <- coeftest(m2_rvl, vcov = vcovCL, cluster = ~sender)
m3_rvl_cl <- coeftest(m3_rvl, vcov = vcovCL, cluster = ~sender)
m4_rvl_cl <- coeftest(m4_rvl, vcov = vcovCL, cluster = ~sender)
m5_rvl_cl <- coeftest(m5_rvl, vcov = vcovCL, cluster = ~sender)

m1_jbs_cl <- coeftest(m1_jbs, vcov = vcovCL, cluster = ~sender)
m2_jbs_cl <- coeftest(m2_jbs, vcov = vcovCL, cluster = ~sender)
m3_jbs_cl <- coeftest(m3_jbs, vcov = vcovCL, cluster = ~sender)
m4_jbs_cl <- coeftest(m4_jbs, vcov = vcovCL, cluster = ~sender)
m5_jbs_cl <- coeftest(m5_jbs, vcov = vcovCL, cluster = ~sender)

m1_pm_cl <- coeftest(m1_pm, vcov = vcovCL, cluster = ~sender)
m2_pm_cl <- coeftest(m2_pm, vcov = vcovCL, cluster = ~sender)
m3_pm_cl <- coeftest(m3_pm, vcov = vcovCL, cluster = ~sender)
m4_pm_cl <- coeftest(m4_pm, vcov = vcovCL, cluster = ~sender)
m5_pm_cl <- coeftest(m5_pm, vcov = vcovCL, cluster = ~sender)

m1_dmz_cl <- coeftest(m1_dmz, vcov = vcovCL, cluster = ~sender)
m2_dmz_cl <- coeftest(m2_dmz, vcov = vcovCL, cluster = ~sender)
m3_dmz_cl <- coeftest(m3_dmz, vcov = vcovCL, cluster = ~sender)
m4_dmz_cl <- coeftest(m4_dmz, vcov = vcovCL, cluster = ~sender)
m5_dmz_cl <- coeftest(m5_dmz, vcov = vcovCL, cluster = ~sender)

########## 5 - plots - allf ##########

res <- matrix(nrow = 20, ncol = 4)

res[1,1] <- m1_rvl_cl[2,1]
res[1,2] <- m1_rvl_cl[2,2]
res[1,3] <- "Paris 1"
res[1,4] <- "w1"
res[2,1] <- m2_rvl_cl[2,1]
res[2,2] <- m2_rvl_cl[2,2]
res[2,3] <- "Paris 1"
res[2,4] <- "w2"
res[3,1] <- m3_rvl_cl[2,1]
res[3,2] <- m3_rvl_cl[2,2]
res[3,3] <- "Paris 1"
res[3,4] <- "w3"
res[4,1] <- m4_rvl_cl[2,1]
res[4,2] <- m4_rvl_cl[2,2]
res[4,3] <- "Paris 1"
res[4,4] <- "w4"
res[5,1] <- m5_rvl_cl[2,1]
res[5,2] <- m5_rvl_cl[2,2]
res[5,3] <- "Paris 1"
res[5,4] <- "w6"

res[1+5,1] <- m1_jbs_cl[2,1]
res[1+5,2] <- m1_jbs_cl[2,2]
res[1+5,3] <- "Paris 2"
res[1+5,4] <- "w1"
res[2+5,1] <- m2_jbs_cl[2,1]
res[2+5,2] <- m2_jbs_cl[2,2]
res[2+5,3] <- "Paris 2"
res[2+5,4] <- "w2"
res[3+5,1] <- m3_jbs_cl[2,1]
res[3+5,2] <- m3_jbs_cl[2,2]
res[3+5,3] <- "Paris 2"
res[3+5,4] <- "w3"
res[4+5,1] <- m4_jbs_cl[2,1]
res[4+5,2] <- m4_jbs_cl[2,2]
res[4+5,3] <- "Paris 2"
res[4+5,4] <- "w4"
res[5+5,1] <- m5_jbs_cl[2,1]
res[5+5,2] <- m5_jbs_cl[2,2]
res[5+5,3] <- "Paris 2"
res[5+5,4] <- "w6"

res[1+10,1] <- m1_pm_cl[2,1]
res[1+10,2] <- m1_pm_cl[2,2]
res[1+10,3] <- "Savoie 1"
res[1+10,4] <- "w1"
res[2+10,1] <- m2_pm_cl[2,1]
res[2+10,2] <- m2_pm_cl[2,2]
res[2+10,3] <- "Savoie 1"
res[2+10,4] <- "w2"
res[3+10,1] <- m3_pm_cl[2,1]
res[3+10,2] <- m3_pm_cl[2,2]
res[3+10,3] <- "Savoie 1"
res[3+10,4] <- "w3"
res[4+10,1] <- m4_pm_cl[2,1]
res[4+10,2] <- m4_pm_cl[2,2]
res[4+10,3] <- "Savoie 1"
res[4+10,4] <- "w4"
res[5+10,1] <- m5_pm_cl[2,1]
res[5+10,2] <- m5_pm_cl[2,2]
res[5+10,3] <- "Savoie 1"
res[5+10,4] <- "w6"

res[1+15,1] <- m1_dmz_cl[2,1]
res[1+15,2] <- m1_dmz_cl[2,2]
res[1+15,3] <- "Savoie 2"
res[1+15,4] <- "w1"
res[2+15,1] <- m2_dmz_cl[2,1]
res[2+15,2] <- m2_dmz_cl[2,2]
res[2+15,3] <- "Savoie 2"
res[2+15,4] <- "w2"
res[3+15,1] <- m3_dmz_cl[2,1]
res[3+15,2] <- m3_dmz_cl[2,2]
res[3+15,3] <- "Savoie 2"
res[3+15,4] <- "w3"
res[4+15,1] <- m4_dmz_cl[2,1]
res[4+15,2] <- m4_dmz_cl[2,2]
res[4+15,3] <- "Savoie 2"
res[4+15,4] <- "w4"
res[5+15,1] <- m5_dmz_cl[2,1]
res[5+15,2] <- m5_dmz_cl[2,2]
res[5+15,3] <- "Savoie 2"
res[5+15,4] <- "w6"

res <- data.frame(res,stringsAsFactors = F)
names(res) <- c("coef","sd","school","wave")
res$coef <- round(as.numeric(res$coef), digits = 4)
res$sd <- round(as.numeric(res$sd), digits = 4)
res$lower_bound <- res$coef - res$sd*2
res$upper_bound <- res$coef + res$sd*2

ggplot(res,aes(x = wave,y = coef, group = school, colour = school)) +
  geom_point(aes(shape=school, color=school), size = 8) +
  #geom_ribbon(aes(ymin=lower_bound,ymax=upper_bound),alpha=0.1,linetype="solid",lwd=1) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = .2) +
  geom_line(linetype="solid",lwd=1.2) +
  geom_hline(yintercept=0, color = "grey50", size=1) +
  scale_y_reverse() +
  scale_x_discrete(limits=c("w1","w2","w3","w4","w5 (missing)","w6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n (5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Collège",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  scale_shape_manual(name = "Collège",
                     values = c(15,16,17,18),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  ggtitle("Coefficients log-odd d\'homophilie sociale, par collège et vague \n Amis") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24)) +
  labs(x="vague",color="Collège")

##

ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement variance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Paris 1 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

########## 6 - tests de robustesse - AME ##########

# Le code fera gf ou allf selon que la section 2 ou 4 a ete executee avant.

library(margins)
ame1_rvl <- margins(m1_rvl)
ame2_rvl <- margins(m2_rvl)
ame3_rvl <- margins(m3_rvl)
ame4_rvl <- margins(m4_rvl)
ame5_rvl <- margins(m5_rvl)

ame1_jbs <- margins(m1_jbs)
ame2_jbs <- margins(m2_jbs)
ame3_jbs <- margins(m3_jbs)
ame4_jbs <- margins(m4_jbs)
ame5_jbs <- margins(m5_jbs)

ame1_pm <- margins(m1_pm)
ame2_pm <- margins(m2_pm)
ame3_pm <- margins(m3_pm)
ame4_pm <- margins(m4_pm)
ame5_pm <- margins(m5_pm)

ame1_dmz <- margins(m1_dmz)
ame2_dmz <- margins(m2_dmz)
ame3_dmz <- margins(m3_dmz)
ame4_dmz <- margins(m4_dmz)
ame5_dmz <- margins(m5_dmz)

###

res <- matrix(nrow = 20, ncol = 4)

res[1,1] <- summary(ame1_rvl)[1,"AME"]
res[1,2] <- summary(ame1_rvl)[1,"SE"]
res[1,3] <- "Paris 1"
res[1,4] <- "w1"
res[2,1] <- summary(ame2_rvl)[1,"AME"]
res[2,2] <- summary(ame2_rvl)[1,"SE"]
res[2,3] <- "Paris 1"
res[2,4] <- "w2"
res[3,1] <- summary(ame3_rvl)[1,"AME"]
res[3,2] <- summary(ame3_rvl)[1,"SE"]
res[3,3] <- "Paris 1"
res[3,4] <- "w3"
res[4,1] <- summary(ame4_rvl)[1,"AME"]
res[4,2] <- summary(ame4_rvl)[1,"SE"]
res[4,3] <- "Paris 1"
res[4,4] <- "w4"
res[5,1] <- summary(ame5_rvl)[1,"AME"]
res[5,2] <- summary(ame5_rvl)[1,"SE"]
res[5,3] <- "Paris 1"
res[5,4] <- "w6"

res[1+5,1] <- summary(ame1_jbs)[1,"AME"]
res[1+5,2] <- summary(ame1_jbs)[1,"SE"]
res[1+5,3] <- "Paris 2"
res[1+5,4] <- "w1"
res[2+5,1] <- summary(ame2_jbs)[1,"AME"]
res[2+5,2] <- summary(ame2_jbs)[1,"SE"]
res[2+5,3] <- "Paris 2"
res[2+5,4] <- "w2"
res[3+5,1] <- summary(ame3_jbs)[1,"AME"]
res[3+5,2] <- summary(ame3_jbs)[1,"SE"]
res[3+5,3] <- "Paris 2"
res[3+5,4] <- "w3"
res[4+5,1] <- summary(ame4_jbs)[1,"AME"]
res[4+5,2] <- summary(ame4_jbs)[1,"SE"]
res[4+5,3] <- "Paris 2"
res[4+5,4] <- "w4"
res[5+5,1] <- summary(ame5_jbs)[1,"AME"]
res[5+5,2] <- summary(ame5_jbs)[1,"SE"]
res[5+5,3] <- "Paris 2"
res[5+5,4] <- "w6"

res[1+10,1] <- summary(ame1_pm)[1,"AME"]
res[1+10,2] <- summary(ame1_pm)[1,"SE"]
res[1+10,3] <- "Savoie 1"
res[1+10,4] <- "w1"
res[2+10,1] <- summary(ame2_pm)[1,"AME"]
res[2+10,2] <- summary(ame2_pm)[1,"SE"]
res[2+10,3] <- "Savoie 1"
res[2+10,4] <- "w2"
res[3+10,1] <- summary(ame3_pm)[1,"AME"]
res[3+10,2] <- summary(ame3_pm)[1,"SE"]
res[3+10,3] <- "Savoie 1"
res[3+10,4] <- "w3"
res[4+10,1] <- summary(ame4_pm)[1,"AME"]
res[4+10,2] <- summary(ame4_pm)[1,"SE"]
res[4+10,3] <- "Savoie 1"
res[4+10,4] <- "w4"
res[5+10,1] <- summary(ame5_pm)[1,"AME"]
res[5+10,2] <- summary(ame5_pm)[1,"SE"]
res[5+10,3] <- "Savoie 1"
res[5+10,4] <- "w6"

res[1+15,1] <- summary(ame1_dmz)[1,"AME"]
res[1+15,2] <- summary(ame1_dmz)[1,"SE"]
res[1+15,3] <- "Savoie 2"
res[1+15,4] <- "w1"
res[2+15,1] <- summary(ame2_dmz)[1,"AME"]
res[2+15,2] <- summary(ame2_dmz)[1,"SE"]
res[2+15,3] <- "Savoie 2"
res[2+15,4] <- "w2"
res[3+15,1] <- summary(ame3_dmz)[1,"AME"]
res[3+15,2] <- summary(ame3_dmz)[1,"SE"]
res[3+15,3] <- "Savoie 2"
res[3+15,4] <- "w3"
res[4+15,1] <- summary(ame4_dmz)[1,"AME"]
res[4+15,2] <- summary(ame4_dmz)[1,"SE"]
res[4+15,3] <- "Savoie 2"
res[4+15,4] <- "w4"
res[5+15,1] <- summary(ame5_dmz)[1,"AME"]
res[5+15,2] <- summary(ame5_dmz)[1,"SE"]
res[5+15,3] <- "Savoie 2"
res[5+15,4] <- "w6"

###

res <- data.frame(res,stringsAsFactors = F)
names(res) <- c("coef","sd","school","wave")
res$coef <- round(as.numeric(res$coef), digits = 4)
res$sd <- round(as.numeric(res$sd), digits = 4)
res$lower_bound <- res$coef - res$sd*2
res$upper_bound <- res$coef + res$sd*2

res

ggplot(res,aes(x = wave,y = coef, group = school, colour = school)) +
  geom_point(aes(shape=school, color=school), size = 8) +
  #geom_ribbon(aes(ymin=lower_bound,ymax=upper_bound),alpha=0.1,linetype="solid",lwd=1) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = .2) +
  geom_line(linetype="solid",lwd=1.2) +
  geom_hline(yintercept=0, color = "grey50", size=1) +
  scale_y_reverse() +
  scale_x_discrete(limits=c("w1","w2","w3","w4","w5 (missing)","w6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n (5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Collège",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  scale_shape_manual(name = "Collège",
                     values = c(15,16,17,18),
                     breaks = c("Paris 1","Paris 2","Savoie 1","Savoie 2"),
                     labels = c("Paris 1","Paris 2","Savoie 1","Savoie 2")) +
  ggtitle("AME d\'homophilie sociale, par collège et vague \n Ami·es / Très bon·nes ami·es") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24)) +
  labs(x="vague",color="Collège")


############################################################
############ 4 - Comparaison entre homophilies #############
############################################################
########## 1 - logistic regressions ##########
##### data frames par college et vague #####

net_rvl <- lapply(gf_rvl,off_NA,"factbrev")
net_jbs <- lapply(gf_jbs,off_NA,"factbrev")
net_pm <- lapply(gf_pm,off_NA,"factbrev")
net_dmz <- lapply(gf_dmz,off_NA,"factbrev")

net_rvl <- lapply(net_rvl,off_NA,"ethn")
net_jbs <- lapply(net_jbs,off_NA,"ethn")
net_pm <- lapply(net_pm,off_NA,"ethn")
net_dmz <- lapply(net_dmz,off_NA,"ethn")

net_rvl <- lapply(net_rvl,off_NA,"grade")
net_jbs[[1]] <- off_NA(net_jbs[[1]],"grade")
net_pm <- lapply(net_pm,off_NA,"grade")
net_dmz <- lapply(net_dmz,off_NA,"grade")

###

net2_rvl <- lapply(allf_rvl,off_NA,"factbrev")
net2_jbs <- lapply(allf_jbs,off_NA,"factbrev")
net2_pm <- lapply(allf_pm,off_NA,"factbrev")
net2_dmz <- lapply(allf_dmz,off_NA,"factbrev")

net2_rvl <- lapply(net2_rvl,off_NA,"ethn")
net2_jbs <- lapply(net2_jbs,off_NA,"ethn")
net2_pm <- lapply(net2_pm,off_NA,"ethn")
net2_dmz <- lapply(net2_dmz,off_NA,"ethn")

net2_rvl <- lapply(net2_rvl,off_NA,"grade")
net2_jbs[[1]] <- off_NA(net2_jbs[[1]],"grade")
net2_pm <- lapply(net2_pm,off_NA,"grade")
net2_dmz <- lapply(net2_dmz,off_NA,"grade")

###

mynet <- net_rvl[[1]]
mynet2 <- net2_rvl[[1]]
drvl_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_1) <- c("sender","receiver")
drvl_1$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_1[,i] <- as.vector(aa[[i]])
}
drvl_1 <- drvl_1[-which(drvl_1$sender==drvl_1$receiver),]
head(drvl_1)

mynet <- net_rvl[[2]]
mynet2 <- net2_rvl[[2]]
drvl_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_2) <- c("sender","receiver")
drvl_2$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_2[,i] <- as.vector(aa[[i]])
}
drvl_2 <- drvl_2[-which(drvl_2$sender==drvl_2$receiver),]
head(drvl_2)

mynet <- net_rvl[[3]]
mynet2 <- net2_rvl[[3]]
drvl_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_3) <- c("sender","receiver")
drvl_3$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_3[,i] <- as.vector(aa[[i]])
}
drvl_3 <- drvl_3[-which(drvl_3$sender==drvl_3$receiver),]
head(drvl_3)

mynet <- net_rvl[[4]]
mynet2 <- net2_rvl[[4]]
drvl_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_4) <- c("sender","receiver")
drvl_4$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_4[,i] <- as.vector(aa[[i]])
}
drvl_4 <- drvl_4[-which(drvl_4$sender==drvl_4$receiver),]
head(drvl_4)

mynet <- net_rvl[[5]]
mynet2 <- net2_rvl[[5]]
drvl_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(drvl_5) <- c("sender","receiver")
drvl_5$tie_gf <- as.vector(as.sociomatrix(mynet))
drvl_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
drvl_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
drvl_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
drvl_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
drvl_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
drvl_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
drvl_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
drvl_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
drvl_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
drvl_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
drvl_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  drvl_5[,i] <- as.vector(aa[[i]])
}
drvl_5 <- drvl_5[-which(drvl_5$sender==drvl_5$receiver),]
head(drvl_5)

###

mynet <- net_jbs[[1]]
mynet2 <- net2_jbs[[1]]
djbs_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_1) <- c("sender","receiver")
djbs_1$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_1[,i] <- as.vector(aa[[i]])
}
djbs_1 <- djbs_1[-which(djbs_1$sender==djbs_1$receiver),]
head(djbs_1)

mynet <- net_jbs[[2]]
mynet2 <- net2_jbs[[2]]
djbs_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_2) <- c("sender","receiver")
djbs_2$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_2[,i] <- as.vector(aa[[i]])
}
djbs_2 <- djbs_2[-which(djbs_2$sender==djbs_2$receiver),]
head(djbs_2)

mynet <- net_jbs[[3]]
mynet2 <- net2_jbs[[3]]
djbs_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_3) <- c("sender","receiver")
djbs_3$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_3[,i] <- as.vector(aa[[i]])
}
djbs_3 <- djbs_3[-which(djbs_3$sender==djbs_3$receiver),]
head(djbs_3)

mynet <- net_jbs[[4]]
mynet2 <- net2_jbs[[4]]
djbs_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_4) <- c("sender","receiver")
djbs_4$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_4[,i] <- as.vector(aa[[i]])
}
djbs_4 <- djbs_4[-which(djbs_4$sender==djbs_4$receiver),]
head(djbs_4)

mynet <- net_jbs[[5]]
mynet2 <- net2_jbs[[5]]
djbs_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(djbs_5) <- c("sender","receiver")
djbs_5$tie_gf <- as.vector(as.sociomatrix(mynet))
djbs_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
djbs_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
djbs_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
djbs_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
djbs_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
djbs_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
djbs_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
djbs_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
djbs_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
djbs_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
djbs_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  djbs_5[,i] <- as.vector(aa[[i]])
}
djbs_5 <- djbs_5[-which(djbs_5$sender==djbs_5$receiver),]
head(djbs_5)

###

mynet <- net_pm[[1]]
mynet2 <- net2_pm[[1]]
dpm_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_1) <- c("sender","receiver")
dpm_1$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_1[,i] <- as.vector(aa[[i]])
}
dpm_1 <- dpm_1[-which(dpm_1$sender==dpm_1$receiver),]
head(dpm_1)

mynet <- net_pm[[2]]
mynet2 <- net2_pm[[2]]
dpm_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_2) <- c("sender","receiver")
dpm_2$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_2[,i] <- as.vector(aa[[i]])
}
dpm_2 <- dpm_2[-which(dpm_2$sender==dpm_2$receiver),]
head(dpm_2)

mynet <- net_pm[[3]]
mynet2 <- net2_pm[[3]]
dpm_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_3) <- c("sender","receiver")
dpm_3$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_3[,i] <- as.vector(aa[[i]])
}
dpm_3 <- dpm_3[-which(dpm_3$sender==dpm_3$receiver),]
head(dpm_3)

mynet <- net_pm[[4]]
mynet2 <- net2_pm[[4]]
dpm_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_4) <- c("sender","receiver")
dpm_4$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_4[,i] <- as.vector(aa[[i]])
}
dpm_4 <- dpm_4[-which(dpm_4$sender==dpm_4$receiver),]
head(dpm_4)

mynet <- net_pm[[5]]
mynet2 <- net2_pm[[5]]
dpm_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(dpm_5) <- c("sender","receiver")
dpm_5$tie_gf <- as.vector(as.sociomatrix(mynet))
dpm_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
dpm_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
dpm_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
dpm_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
dpm_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
dpm_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
dpm_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
dpm_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
dpm_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
dpm_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
dpm_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  dpm_5[,i] <- as.vector(aa[[i]])
}
dpm_5 <- dpm_5[-which(dpm_5$sender==dpm_5$receiver),]
head(dpm_5)

###

mynet <- net_dmz[[1]]
mynet2 <- net2_dmz[[1]]
ddmz_1 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_1) <- c("sender","receiver")
ddmz_1$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_1$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_1$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_1$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_1$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_1$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_1$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_1$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_1$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_1$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_1$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_1$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_1[,i] <- as.vector(aa[[i]])
}
ddmz_1 <- ddmz_1[-which(ddmz_1$sender==ddmz_1$receiver),]
head(ddmz_1)

mynet <- net_dmz[[2]]
mynet2 <- net2_dmz[[2]]
ddmz_2 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_2) <- c("sender","receiver")
ddmz_2$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_2$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_2$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_2$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_2$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_2$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_2$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_2$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_2$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_2$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_2$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_2$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_2[,i] <- as.vector(aa[[i]])
}
ddmz_2 <- ddmz_2[-which(ddmz_2$sender==ddmz_2$receiver),]
head(ddmz_2)

mynet <- net_dmz[[3]]
mynet2 <- net2_dmz[[3]]
ddmz_3 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_3) <- c("sender","receiver")
ddmz_3$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_3$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_3$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_3$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_3$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_3$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_3$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_3$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_3$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_3$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_3$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_3$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_3[,i] <- as.vector(aa[[i]])
}
ddmz_3 <- ddmz_3[-which(ddmz_3$sender==ddmz_3$receiver),]
head(ddmz_3)

mynet <- net_dmz[[4]]
mynet2 <- net2_dmz[[4]]
ddmz_4 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_4) <- c("sender","receiver")
ddmz_4$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_4$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_4$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_4$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_4$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_4$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_4$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_4$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_4$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_4$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_4$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_4$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_4[,i] <- as.vector(aa[[i]])
}
ddmz_4 <- ddmz_4[-which(ddmz_4$sender==ddmz_4$receiver),]
head(ddmz_4)

mynet <- net_dmz[[5]]
mynet2 <- net2_dmz[[5]]
ddmz_5 <- expand.grid(as.character(mynet%v%"vertex.names"),as.character(mynet%v%"vertex.names"))
names(ddmz_5) <- c("sender","receiver")
ddmz_5$tie_gf <- as.vector(as.sociomatrix(mynet))
ddmz_5$tie_allf <- as.vector(as.sociomatrix(mynet2))
ddmz_5$factbrev_absdiff <- as.vector(mat_att(mynet,"factbrev")$diff)
ddmz_5$factbrev_sender <- as.vector(mat_att(mynet,"factbrev")$sender)
ddmz_5$factbrev_receiver <- as.vector(mat_att(mynet,"factbrev")$receiver)
ddmz_5$grade_absdiff <- as.vector(mat_att(mynet,"grade")$diff)
ddmz_5$grade_sender <- as.vector(mat_att(mynet,"grade")$sender)
ddmz_5$grade_receiver <- as.vector(mat_att(mynet,"grade")$receiver) 
ddmz_5$sex_nodematch <- as.vector(mat_att(mynet,"sex")$match)
ddmz_5$sex_sender <- as.vector(mat_att(mynet,"sex")$M_sender)
ddmz_5$sex_receiver <- as.vector(mat_att(mynet,"sex")$M_receiver)
ddmz_5$ethn_nodematch <- as.vector(mat_att(mynet,"ethn")$match)
aa <- mat_att(mynet,"ethn")
names_ethn <- sort(names(aa)[which(unlist(lapply(strsplit(names(aa),"_"),function(x){"sender"%in%x | "receiver"%in%x})))])
for(i in names_ethn){
  ddmz_5[,i] <- as.vector(aa[[i]])
}
ddmz_5 <- ddmz_5[-which(ddmz_5$sender==ddmz_5$receiver),]
head(ddmz_5)

##### regressions - gf #####

"%out%" <- Negate("%in%")

mynames <- unique(c(names(drvl_1),names(drvl_2),names(drvl_3),names(drvl_4),names(drvl_5),
                    names(djbs_1),names(djbs_2),names(djbs_3),names(djbs_4),names(djbs_5),
                    names(dpm_1),names(dpm_2),names(dpm_3),names(dpm_4),names(dpm_5),
                    names(ddmz_1),names(ddmz_2),names(ddmz_3),names(ddmz_4),names(ddmz_5)))
mynames <-mynames[mynames%out%c("sender","receiver","tie_gf","tie_allf","auto_sender","auto_receiver")] # autochtones comme categorie de reference

m1_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
m2_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_2)),collapse="+"),collapse = " ")),
              data = drvl_2,family="binomial")
m3_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_3)),collapse="+"),collapse = " ")),
              data = drvl_3,family="binomial")
m4_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_4)),collapse="+"),collapse = " ")),
              data = drvl_4,family="binomial")
m5_rvl <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(drvl_5)),collapse="+"),collapse = " ")),
              data = drvl_5,family="binomial")

m1_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
m2_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_2)),collapse="+"),collapse = " ")),
              data = djbs_2,family="binomial")
m3_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_3)),collapse="+"),collapse = " ")),
              data = djbs_3,family="binomial")
m4_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_4)),collapse="+"),collapse = " ")),
              data = djbs_4,family="binomial")
m5_jbs <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(djbs_5)),collapse="+"),collapse = " ")),
              data = djbs_5,family="binomial")

m1_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
m2_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_2)),collapse="+"),collapse = " ")),
             data = dpm_2,family="binomial")
m3_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_3)),collapse="+"),collapse = " ")),
             data = dpm_3,family="binomial")
m4_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_4)),collapse="+"),collapse = " ")),
             data = dpm_4,family="binomial")
m5_pm <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(dpm_5)),collapse="+"),collapse = " ")),
             data = dpm_5,family="binomial")

m1_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")
m2_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_2)),collapse="+"),collapse = " ")),
              data = ddmz_2,family="binomial")
m3_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_3)),collapse="+"),collapse = " ")),
              data = ddmz_3,family="binomial")
m4_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_4)),collapse="+"),collapse = " ")),
              data = ddmz_4,family="binomial")
m5_dmz <- glm(as.formula(paste("tie_gf ~",paste(intersect(mynames,names(ddmz_5)),collapse="+"),collapse = " ")),
              data = ddmz_5,family="binomial")

##### regressions - allf #####

"%out%" <- Negate("%in%")

mynames <- unique(c(names(drvl_1),names(drvl_2),names(drvl_3),names(drvl_4),names(drvl_5),
                    names(djbs_1),names(djbs_2),names(djbs_3),names(djbs_4),names(djbs_5),
                    names(dpm_1),names(dpm_2),names(dpm_3),names(dpm_4),names(dpm_5),
                    names(ddmz_1),names(ddmz_2),names(ddmz_3),names(ddmz_4),names(ddmz_5)))
mynames <-mynames[mynames%out%c("sender","receiver","tie_gf","tie_allf","auto_sender","auto_receiver")] # autochtones comme categorie de reference

n1_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_1)),collapse="+"),collapse = " ")),
              data = drvl_1,family="binomial")
n2_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_2)),collapse="+"),collapse = " ")),
              data = drvl_2,family="binomial")
n3_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_3)),collapse="+"),collapse = " ")),
              data = drvl_3,family="binomial")
n4_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_4)),collapse="+"),collapse = " ")),
              data = drvl_4,family="binomial")
n5_rvl <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(drvl_5)),collapse="+"),collapse = " ")),
              data = drvl_5,family="binomial")

n1_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_1)),collapse="+"),collapse = " ")),
              data = djbs_1,family="binomial")
n2_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_2)),collapse="+"),collapse = " ")),
              data = djbs_2,family="binomial")
n3_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_3)),collapse="+"),collapse = " ")),
              data = djbs_3,family="binomial")
n4_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_4)),collapse="+"),collapse = " ")),
              data = djbs_4,family="binomial")
n5_jbs <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(djbs_5)),collapse="+"),collapse = " ")),
              data = djbs_5,family="binomial")

n1_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_1)),collapse="+"),collapse = " ")),
             data = dpm_1,family="binomial")
n2_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_2)),collapse="+"),collapse = " ")),
             data = dpm_2,family="binomial")
n3_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_3)),collapse="+"),collapse = " ")),
             data = dpm_3,family="binomial")
n4_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_4)),collapse="+"),collapse = " ")),
             data = dpm_4,family="binomial")
n5_pm <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(dpm_5)),collapse="+"),collapse = " ")),
             data = dpm_5,family="binomial")

n1_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_1)),collapse="+"),collapse = " ")),
              data = ddmz_1,family="binomial")
n2_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_2)),collapse="+"),collapse = " ")),
              data = ddmz_2,family="binomial")
n3_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_3)),collapse="+"),collapse = " ")),
              data = ddmz_3,family="binomial")
n4_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_4)),collapse="+"),collapse = " ")),
              data = ddmz_4,family="binomial")
n5_dmz <- glm(as.formula(paste("tie_allf ~",paste(intersect(mynames,names(ddmz_5)),collapse="+"),collapse = " ")),
              data = ddmz_5,family="binomial")

##### lecture des resultats #####

stargazer(m1_rvl,m2_rvl,m3_rvl,m4_rvl,m5_rvl,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(m1_jbs,m2_jbs,m3_jbs,m4_jbs,m5_jbs,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(m1_pm,m2_pm,m3_pm,m4_pm,m5_pm,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(m1_dmz,m2_dmz,m3_dmz,m4_dmz,m5_dmz,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
# Ca devrait etre strictement egal aux ERGMs (check = ok).

stargazer(n1_rvl,n2_rvl,n3_rvl,n4_rvl,n5_rvl,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(n1_jbs,n2_jbs,n3_jbs,n4_jbs,n5_jbs,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(n1_pm,n2_pm,n3_pm,n4_pm,n5_pm,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)
stargazer(n1_dmz,n2_dmz,n3_dmz,n4_dmz,n5_dmz,type="text",column.labels = rep("tie",5),dep.var.labels.include=F)

###

# Clustered s.e.

library(sandwich)
library(lmtest)

m1_rvl_cl <- coeftest(m1_rvl, vcov = vcovCL, cluster = ~sender)
m2_rvl_cl <- coeftest(m2_rvl, vcov = vcovCL, cluster = ~sender)
m3_rvl_cl <- coeftest(m3_rvl, vcov = vcovCL, cluster = ~sender)
m4_rvl_cl <- coeftest(m4_rvl, vcov = vcovCL, cluster = ~sender)
m5_rvl_cl <- coeftest(m5_rvl, vcov = vcovCL, cluster = ~sender)

m1_jbs_cl <- coeftest(m1_jbs, vcov = vcovCL, cluster = ~sender)
m2_jbs_cl <- coeftest(m2_jbs, vcov = vcovCL, cluster = ~sender)
m3_jbs_cl <- coeftest(m3_jbs, vcov = vcovCL, cluster = ~sender)
m4_jbs_cl <- coeftest(m4_jbs, vcov = vcovCL, cluster = ~sender)
m5_jbs_cl <- coeftest(m5_jbs, vcov = vcovCL, cluster = ~sender)

m1_pm_cl <- coeftest(m1_pm, vcov = vcovCL, cluster = ~sender)
m2_pm_cl <- coeftest(m2_pm, vcov = vcovCL, cluster = ~sender)
m3_pm_cl <- coeftest(m3_pm, vcov = vcovCL, cluster = ~sender)
m4_pm_cl <- coeftest(m4_pm, vcov = vcovCL, cluster = ~sender)
m5_pm_cl <- coeftest(m5_pm, vcov = vcovCL, cluster = ~sender)

m1_dmz_cl <- coeftest(m1_dmz, vcov = vcovCL, cluster = ~sender)
m2_dmz_cl <- coeftest(m2_dmz, vcov = vcovCL, cluster = ~sender)
m3_dmz_cl <- coeftest(m3_dmz, vcov = vcovCL, cluster = ~sender)
m4_dmz_cl <- coeftest(m4_dmz, vcov = vcovCL, cluster = ~sender)
m5_dmz_cl <- coeftest(m5_dmz, vcov = vcovCL, cluster = ~sender)

n1_rvl_cl <- coeftest(n1_rvl, vcov = vcovCL, cluster = ~sender)
n2_rvl_cl <- coeftest(n2_rvl, vcov = vcovCL, cluster = ~sender)
n3_rvl_cl <- coeftest(n3_rvl, vcov = vcovCL, cluster = ~sender)
n4_rvl_cl <- coeftest(n4_rvl, vcov = vcovCL, cluster = ~sender)
n5_rvl_cl <- coeftest(n5_rvl, vcov = vcovCL, cluster = ~sender)

n1_jbs_cl <- coeftest(n1_jbs, vcov = vcovCL, cluster = ~sender)
n2_jbs_cl <- coeftest(n2_jbs, vcov = vcovCL, cluster = ~sender)
n3_jbs_cl <- coeftest(n3_jbs, vcov = vcovCL, cluster = ~sender)
n4_jbs_cl <- coeftest(n4_jbs, vcov = vcovCL, cluster = ~sender)
n5_jbs_cl <- coeftest(n5_jbs, vcov = vcovCL, cluster = ~sender)

n1_pm_cl <- coeftest(n1_pm, vcov = vcovCL, cluster = ~sender)
n2_pm_cl <- coeftest(n2_pm, vcov = vcovCL, cluster = ~sender)
n3_pm_cl <- coeftest(n3_pm, vcov = vcovCL, cluster = ~sender)
n4_pm_cl <- coeftest(n4_pm, vcov = vcovCL, cluster = ~sender)
n5_pm_cl <- coeftest(n5_pm, vcov = vcovCL, cluster = ~sender)

n1_dmz_cl <- coeftest(n1_dmz, vcov = vcovCL, cluster = ~sender)
n2_dmz_cl <- coeftest(n2_dmz, vcov = vcovCL, cluster = ~sender)
n3_dmz_cl <- coeftest(n3_dmz, vcov = vcovCL, cluster = ~sender)
n4_dmz_cl <- coeftest(n4_dmz, vcov = vcovCL, cluster = ~sender)
n5_dmz_cl <- coeftest(n5_dmz, vcov = vcovCL, cluster = ~sender)

stargazer(m1_rvl_cl,m2_rvl_cl,m3_rvl_cl,m4_rvl_cl,m5_rvl_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(m1_jbs_cl,m2_jbs_cl,m3_jbs_cl,m4_jbs_cl,m5_jbs_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(m1_pm_cl,m2_pm_cl,m3_pm_cl,m4_pm_cl,m5_pm_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(m1_dmz_cl,m2_dmz_cl,m3_dmz_cl,m4_dmz_cl,m5_dmz_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)

stargazer(n1_rvl_cl,n2_rvl_cl,n3_rvl_cl,n4_rvl_cl,n5_rvl_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(n1_jbs_cl,n2_jbs_cl,n3_jbs_cl,n4_jbs_cl,n5_jbs_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(n1_pm_cl,n2_pm_cl,n3_pm_cl,n4_pm_cl,n5_pm_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
stargazer(n1_dmz_cl,n2_dmz_cl,n3_dmz_cl,n4_dmz_cl,n5_dmz_cl,type="text",
          column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)

###

# stargazer(m1_rvl_cl,m2_rvl_cl,m3_rvl_cl,m4_rvl_cl,m5_rvl_cl,type="html", out="a1.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(m1_jbs_cl,m2_jbs_cl,m3_jbs_cl,m4_jbs_cl,m5_jbs_cl,type="html", out="a2.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(m1_pm_cl,m2_pm_cl,m3_pm_cl,m4_pm_cl,m5_pm_cl,type="html", out="a3.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(m1_dmz_cl,m2_dmz_cl,m3_dmz_cl,m4_dmz_cl,m5_dmz_cl,type="html", out="a4.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# 
# stargazer(n1_rvl_cl,n2_rvl_cl,n3_rvl_cl,n4_rvl_cl,n5_rvl_cl,type="html", out="b1.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(n1_jbs_cl,n2_jbs_cl,n3_jbs_cl,n4_jbs_cl,n5_jbs_cl,type="html", out="b2.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(n1_pm_cl,n2_pm_cl,n3_pm_cl,n4_pm_cl,n5_pm_cl,type="html", out="b3.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)
# stargazer(n1_dmz_cl,n2_dmz_cl,n3_dmz_cl,n4_dmz_cl,n5_dmz_cl,type="html", out="b4.html",
#           column.labels = c("w1","w2","w3","w4","w6"),dep.var.labels.include=F)

##### comparaison des tailles d effets - estimation gf #####

# On estime les modeles en enlevant un parametre.

m1_rvl_factbrev <- glm(update(m1_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_1,family="binomial")
m1_rvl_grade <- glm(update(m1_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_1,family="binomial")
m1_rvl_ethn <- glm(update(m1_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_1,family="binomial")
m1_rvl_sex <- glm(update(m1_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_1,family="binomial")

m2_rvl_factbrev <- glm(update(m2_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_2,family="binomial")
m2_rvl_grade <- glm(update(m2_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_2,family="binomial")
m2_rvl_ethn <- glm(update(m2_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_2,family="binomial")
m2_rvl_sex <- glm(update(m2_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_2,family="binomial")

m3_rvl_factbrev <- glm(update(m3_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_3,family="binomial")
m3_rvl_grade <- glm(update(m3_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_3,family="binomial")
m3_rvl_ethn <- glm(update(m3_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_3,family="binomial")
m3_rvl_sex <- glm(update(m3_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_3,family="binomial")

m4_rvl_factbrev <- glm(update(m4_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_4,family="binomial")
m4_rvl_grade <- glm(update(m4_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_4,family="binomial")
m4_rvl_ethn <- glm(update(m4_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_4,family="binomial")
m4_rvl_sex <- glm(update(m4_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_4,family="binomial")

m5_rvl_factbrev <- glm(update(m5_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_5,family="binomial")
m5_rvl_grade <- glm(update(m5_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_5,family="binomial")
m5_rvl_ethn <- glm(update(m5_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_5,family="binomial")
m5_rvl_sex <- glm(update(m5_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_5,family="binomial")

res_rvl <- matrix(nrow = 4, ncol = 5)
rownames(res_rvl) <- c("factbrev","grade","ethn","sex")
colnames(res_rvl) <- c("v1","v2","v3","v4","v6")
res_rvl_aic <- res_rvl

res_rvl[1,1] <- 100*(m1_rvl_factbrev[["deviance"]] - m1_rvl[["deviance"]])/m1_rvl[["deviance"]]
res_rvl[2,1] <- 100*(m1_rvl_grade[["deviance"]] - m1_rvl[["deviance"]])/m1_rvl[["deviance"]]
res_rvl[3,1] <- 100*(m1_rvl_ethn[["deviance"]] - m1_rvl[["deviance"]])/m1_rvl[["deviance"]]
res_rvl[4,1] <- 100*(m1_rvl_sex[["deviance"]] - m1_rvl[["deviance"]])/m1_rvl[["deviance"]]

res_rvl[1,2] <- 100*(m2_rvl_factbrev[["deviance"]] - m2_rvl[["deviance"]])/m2_rvl[["deviance"]]
res_rvl[2,2] <- 100*(m2_rvl_grade[["deviance"]] - m2_rvl[["deviance"]])/m2_rvl[["deviance"]]
res_rvl[3,2] <- 100*(m2_rvl_ethn[["deviance"]] - m2_rvl[["deviance"]])/m2_rvl[["deviance"]]
res_rvl[4,2] <- 100*(m2_rvl_sex[["deviance"]] - m2_rvl[["deviance"]])/m2_rvl[["deviance"]]

res_rvl[1,3] <- 100*(m3_rvl_factbrev[["deviance"]] - m3_rvl[["deviance"]])/m3_rvl[["deviance"]]
res_rvl[2,3] <- 100*(m3_rvl_grade[["deviance"]] - m3_rvl[["deviance"]])/m3_rvl[["deviance"]]
res_rvl[3,3] <- 100*(m3_rvl_ethn[["deviance"]] - m3_rvl[["deviance"]])/m3_rvl[["deviance"]]
res_rvl[4,3] <- 100*(m3_rvl_sex[["deviance"]] - m3_rvl[["deviance"]])/m3_rvl[["deviance"]]

res_rvl[1,4] <- 100*(m4_rvl_factbrev[["deviance"]] - m4_rvl[["deviance"]])/m4_rvl[["deviance"]]
res_rvl[2,4] <- 100*(m4_rvl_grade[["deviance"]] - m4_rvl[["deviance"]])/m4_rvl[["deviance"]]
res_rvl[3,4] <- 100*(m4_rvl_ethn[["deviance"]] - m4_rvl[["deviance"]])/m4_rvl[["deviance"]]
res_rvl[4,4] <- 100*(m4_rvl_sex[["deviance"]] - m4_rvl[["deviance"]])/m4_rvl[["deviance"]]

res_rvl[1,5] <- 100*(m5_rvl_factbrev[["deviance"]] - m5_rvl[["deviance"]])/m5_rvl[["deviance"]]
res_rvl[2,5] <- 100*(m5_rvl_grade[["deviance"]] - m5_rvl[["deviance"]])/m5_rvl[["deviance"]]
res_rvl[3,5] <- 100*(m5_rvl_ethn[["deviance"]] - m5_rvl[["deviance"]])/m5_rvl[["deviance"]]
res_rvl[4,5] <- 100*(m5_rvl_sex[["deviance"]] - m5_rvl[["deviance"]])/m5_rvl[["deviance"]]

res_rvl_aic[1,1] <- 100*(m1_rvl_factbrev[["aic"]] - m1_rvl[["aic"]])/m1_rvl[["aic"]]
res_rvl_aic[2,1] <- 100*(m1_rvl_grade[["aic"]] - m1_rvl[["aic"]])/m1_rvl[["aic"]]
res_rvl_aic[3,1] <- 100*(m1_rvl_ethn[["aic"]] - m1_rvl[["aic"]])/m1_rvl[["aic"]]
res_rvl_aic[4,1] <- 100*(m1_rvl_sex[["aic"]] - m1_rvl[["aic"]])/m1_rvl[["aic"]]

res_rvl_aic[1,2] <- 100*(m2_rvl_factbrev[["aic"]] - m2_rvl[["aic"]])/m2_rvl[["aic"]]
res_rvl_aic[2,2] <- 100*(m2_rvl_grade[["aic"]] - m2_rvl[["aic"]])/m2_rvl[["aic"]]
res_rvl_aic[3,2] <- 100*(m2_rvl_ethn[["aic"]] - m2_rvl[["aic"]])/m2_rvl[["aic"]]
res_rvl_aic[4,2] <- 100*(m2_rvl_sex[["aic"]] - m2_rvl[["aic"]])/m2_rvl[["aic"]]

res_rvl_aic[1,3] <- 100*(m3_rvl_factbrev[["aic"]] - m3_rvl[["aic"]])/m3_rvl[["aic"]]
res_rvl_aic[2,3] <- 100*(m3_rvl_grade[["aic"]] - m3_rvl[["aic"]])/m3_rvl[["aic"]]
res_rvl_aic[3,3] <- 100*(m3_rvl_ethn[["aic"]] - m3_rvl[["aic"]])/m3_rvl[["aic"]]
res_rvl_aic[4,3] <- 100*(m3_rvl_sex[["aic"]] - m3_rvl[["aic"]])/m3_rvl[["aic"]]

res_rvl_aic[1,4] <- 100*(m4_rvl_factbrev[["aic"]] - m4_rvl[["aic"]])/m4_rvl[["aic"]]
res_rvl_aic[2,4] <- 100*(m4_rvl_grade[["aic"]] - m4_rvl[["aic"]])/m4_rvl[["aic"]]
res_rvl_aic[3,4] <- 100*(m4_rvl_ethn[["aic"]] - m4_rvl[["aic"]])/m4_rvl[["aic"]]
res_rvl_aic[4,4] <- 100*(m4_rvl_sex[["aic"]] - m4_rvl[["aic"]])/m4_rvl[["aic"]]

res_rvl_aic[1,5] <- 100*(m5_rvl_factbrev[["aic"]] - m5_rvl[["aic"]])/m5_rvl[["aic"]]
res_rvl_aic[2,5] <- 100*(m5_rvl_grade[["aic"]] - m5_rvl[["aic"]])/m5_rvl[["aic"]]
res_rvl_aic[3,5] <- 100*(m5_rvl_ethn[["aic"]] - m5_rvl[["aic"]])/m5_rvl[["aic"]]
res_rvl_aic[4,5] <- 100*(m5_rvl_sex[["aic"]] - m5_rvl[["aic"]])/m5_rvl[["aic"]]

###

m1_jbs_factbrev <- glm(update(m1_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_1,family="binomial")
m1_jbs_grade <- glm(update(m1_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_1,family="binomial")
m1_jbs_ethn <- glm(update(m1_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_1,family="binomial")
m1_jbs_sex <- glm(update(m1_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_1,family="binomial")

m2_jbs_factbrev <- glm(update(m2_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_2,family="binomial")
m2_jbs_grade <- glm(update(m2_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_2,family="binomial")
m2_jbs_ethn <- glm(update(m2_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_2,family="binomial")
m2_jbs_sex <- glm(update(m2_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_2,family="binomial")

m3_jbs_factbrev <- glm(update(m3_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_3,family="binomial")
m3_jbs_grade <- glm(update(m3_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_3,family="binomial")
m3_jbs_ethn <- glm(update(m3_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_3,family="binomial")
m3_jbs_sex <- glm(update(m3_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_3,family="binomial")

m4_jbs_factbrev <- glm(update(m4_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_4,family="binomial")
m4_jbs_grade <- glm(update(m4_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_4,family="binomial")
m4_jbs_ethn <- glm(update(m4_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_4,family="binomial")
m4_jbs_sex <- glm(update(m4_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_4,family="binomial")

m5_jbs_factbrev <- glm(update(m5_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_5,family="binomial")
m5_jbs_grade <- glm(update(m5_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_5,family="binomial")
m5_jbs_ethn <- glm(update(m5_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_5,family="binomial")
m5_jbs_sex <- glm(update(m5_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_5,family="binomial")

res_jbs <- matrix(nrow = 4, ncol = 5)
rownames(res_jbs) <- c("factbrev","grade","ethn","sex")
colnames(res_jbs) <- c("v1","v2","v3","v4","v6")
res_jbs_aic <- res_jbs

res_jbs[1,1] <- 100*(m1_jbs_factbrev[["deviance"]] - m1_jbs[["deviance"]])/m1_jbs[["deviance"]]
res_jbs[2,1] <- 100*(m1_jbs_grade[["deviance"]] - m1_jbs[["deviance"]])/m1_jbs[["deviance"]]
res_jbs[3,1] <- 100*(m1_jbs_ethn[["deviance"]] - m1_jbs[["deviance"]])/m1_jbs[["deviance"]]
res_jbs[4,1] <- 100*(m1_jbs_sex[["deviance"]] - m1_jbs[["deviance"]])/m1_jbs[["deviance"]]

res_jbs[1,2] <- 100*(m2_jbs_factbrev[["deviance"]] - m2_jbs[["deviance"]])/m2_jbs[["deviance"]]
res_jbs[2,2] <- 100*(m2_jbs_grade[["deviance"]] - m2_jbs[["deviance"]])/m2_jbs[["deviance"]]
res_jbs[3,2] <- 100*(m2_jbs_ethn[["deviance"]] - m2_jbs[["deviance"]])/m2_jbs[["deviance"]]
res_jbs[4,2] <- 100*(m2_jbs_sex[["deviance"]] - m2_jbs[["deviance"]])/m2_jbs[["deviance"]]

res_jbs[1,3] <- 100*(m3_jbs_factbrev[["deviance"]] - m3_jbs[["deviance"]])/m3_jbs[["deviance"]]
res_jbs[2,3] <- 100*(m3_jbs_grade[["deviance"]] - m3_jbs[["deviance"]])/m3_jbs[["deviance"]]
res_jbs[3,3] <- 100*(m3_jbs_ethn[["deviance"]] - m3_jbs[["deviance"]])/m3_jbs[["deviance"]]
res_jbs[4,3] <- 100*(m3_jbs_sex[["deviance"]] - m3_jbs[["deviance"]])/m3_jbs[["deviance"]]

res_jbs[1,4] <- 100*(m4_jbs_factbrev[["deviance"]] - m4_jbs[["deviance"]])/m4_jbs[["deviance"]]
res_jbs[2,4] <- 100*(m4_jbs_grade[["deviance"]] - m4_jbs[["deviance"]])/m4_jbs[["deviance"]]
res_jbs[3,4] <- 100*(m4_jbs_ethn[["deviance"]] - m4_jbs[["deviance"]])/m4_jbs[["deviance"]]
res_jbs[4,4] <- 100*(m4_jbs_sex[["deviance"]] - m4_jbs[["deviance"]])/m4_jbs[["deviance"]]

res_jbs[1,5] <- 100*(m5_jbs_factbrev[["deviance"]] - m5_jbs[["deviance"]])/m5_jbs[["deviance"]]
res_jbs[2,5] <- 100*(m5_jbs_grade[["deviance"]] - m5_jbs[["deviance"]])/m5_jbs[["deviance"]]
res_jbs[3,5] <- 100*(m5_jbs_ethn[["deviance"]] - m5_jbs[["deviance"]])/m5_jbs[["deviance"]]
res_jbs[4,5] <- 100*(m5_jbs_sex[["deviance"]] - m5_jbs[["deviance"]])/m5_jbs[["deviance"]]

res_jbs_aic[1,1] <- 100*(m1_jbs_factbrev[["aic"]] - m1_jbs[["aic"]])/m1_jbs[["aic"]]
res_jbs_aic[2,1] <- 100*(m1_jbs_grade[["aic"]] - m1_jbs[["aic"]])/m1_jbs[["aic"]]
res_jbs_aic[3,1] <- 100*(m1_jbs_ethn[["aic"]] - m1_jbs[["aic"]])/m1_jbs[["aic"]]
res_jbs_aic[4,1] <- 100*(m1_jbs_sex[["aic"]] - m1_jbs[["aic"]])/m1_jbs[["aic"]]

res_jbs_aic[1,2] <- 100*(m2_jbs_factbrev[["aic"]] - m2_jbs[["aic"]])/m2_jbs[["aic"]]
res_jbs_aic[2,2] <- 100*(m2_jbs_grade[["aic"]] - m2_jbs[["aic"]])/m2_jbs[["aic"]]
res_jbs_aic[3,2] <- 100*(m2_jbs_ethn[["aic"]] - m2_jbs[["aic"]])/m2_jbs[["aic"]]
res_jbs_aic[4,2] <- 100*(m2_jbs_sex[["aic"]] - m2_jbs[["aic"]])/m2_jbs[["aic"]]

res_jbs_aic[1,3] <- 100*(m3_jbs_factbrev[["aic"]] - m3_jbs[["aic"]])/m3_jbs[["aic"]]
res_jbs_aic[2,3] <- 100*(m3_jbs_grade[["aic"]] - m3_jbs[["aic"]])/m3_jbs[["aic"]]
res_jbs_aic[3,3] <- 100*(m3_jbs_ethn[["aic"]] - m3_jbs[["aic"]])/m3_jbs[["aic"]]
res_jbs_aic[4,3] <- 100*(m3_jbs_sex[["aic"]] - m3_jbs[["aic"]])/m3_jbs[["aic"]]

res_jbs_aic[1,4] <- 100*(m4_jbs_factbrev[["aic"]] - m4_jbs[["aic"]])/m4_jbs[["aic"]]
res_jbs_aic[2,4] <- 100*(m4_jbs_grade[["aic"]] - m4_jbs[["aic"]])/m4_jbs[["aic"]]
res_jbs_aic[3,4] <- 100*(m4_jbs_ethn[["aic"]] - m4_jbs[["aic"]])/m4_jbs[["aic"]]
res_jbs_aic[4,4] <- 100*(m4_jbs_sex[["aic"]] - m4_jbs[["aic"]])/m4_jbs[["aic"]]

res_jbs_aic[1,5] <- 100*(m5_jbs_factbrev[["aic"]] - m5_jbs[["aic"]])/m5_jbs[["aic"]]
res_jbs_aic[2,5] <- 100*(m5_jbs_grade[["aic"]] - m5_jbs[["aic"]])/m5_jbs[["aic"]]
res_jbs_aic[3,5] <- 100*(m5_jbs_ethn[["aic"]] - m5_jbs[["aic"]])/m5_jbs[["aic"]]
res_jbs_aic[4,5] <- 100*(m5_jbs_sex[["aic"]] - m5_jbs[["aic"]])/m5_jbs[["aic"]]

###

m1_pm_factbrev <- glm(update(m1_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_1,family="binomial")
m1_pm_grade <- glm(update(m1_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_1,family="binomial")
m1_pm_ethn <- glm(update(m1_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_1,family="binomial")
m1_pm_sex <- glm(update(m1_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_1,family="binomial")

m2_pm_factbrev <- glm(update(m2_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_2,family="binomial")
m2_pm_grade <- glm(update(m2_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_2,family="binomial")
m2_pm_ethn <- glm(update(m2_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_2,family="binomial")
m2_pm_sex <- glm(update(m2_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_2,family="binomial")

m3_pm_factbrev <- glm(update(m3_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_3,family="binomial")
m3_pm_grade <- glm(update(m3_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_3,family="binomial")
m3_pm_ethn <- glm(update(m3_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_3,family="binomial")
m3_pm_sex <- glm(update(m3_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_3,family="binomial")

m4_pm_factbrev <- glm(update(m4_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_4,family="binomial")
m4_pm_grade <- glm(update(m4_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_4,family="binomial")
m4_pm_ethn <- glm(update(m4_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_4,family="binomial")
m4_pm_sex <- glm(update(m4_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_4,family="binomial")

m5_pm_factbrev <- glm(update(m5_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_5,family="binomial")
m5_pm_grade <- glm(update(m5_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_5,family="binomial")
m5_pm_ethn <- glm(update(m5_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_5,family="binomial")
m5_pm_sex <- glm(update(m5_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_5,family="binomial")

res_pm <- matrix(nrow = 4, ncol = 5)
rownames(res_pm) <- c("factbrev","grade","ethn","sex")
colnames(res_pm) <- c("v1","v2","v3","v4","v6")
res_pm_aic <- res_pm

res_pm[1,1] <- 100*(m1_pm_factbrev[["deviance"]] - m1_pm[["deviance"]])/m1_pm[["deviance"]]
res_pm[2,1] <- 100*(m1_pm_grade[["deviance"]] - m1_pm[["deviance"]])/m1_pm[["deviance"]]
res_pm[3,1] <- 100*(m1_pm_ethn[["deviance"]] - m1_pm[["deviance"]])/m1_pm[["deviance"]]
res_pm[4,1] <- 100*(m1_pm_sex[["deviance"]] - m1_pm[["deviance"]])/m1_pm[["deviance"]]

res_pm[1,2] <- 100*(m2_pm_factbrev[["deviance"]] - m2_pm[["deviance"]])/m2_pm[["deviance"]]
res_pm[2,2] <- 100*(m2_pm_grade[["deviance"]] - m2_pm[["deviance"]])/m2_pm[["deviance"]]
res_pm[3,2] <- 100*(m2_pm_ethn[["deviance"]] - m2_pm[["deviance"]])/m2_pm[["deviance"]]
res_pm[4,2] <- 100*(m2_pm_sex[["deviance"]] - m2_pm[["deviance"]])/m2_pm[["deviance"]]

res_pm[1,3] <- 100*(m3_pm_factbrev[["deviance"]] - m3_pm[["deviance"]])/m3_pm[["deviance"]]
res_pm[2,3] <- 100*(m3_pm_grade[["deviance"]] - m3_pm[["deviance"]])/m3_pm[["deviance"]]
res_pm[3,3] <- 100*(m3_pm_ethn[["deviance"]] - m3_pm[["deviance"]])/m3_pm[["deviance"]]
res_pm[4,3] <- 100*(m3_pm_sex[["deviance"]] - m3_pm[["deviance"]])/m3_pm[["deviance"]]

res_pm[1,4] <- 100*(m4_pm_factbrev[["deviance"]] - m4_pm[["deviance"]])/m4_pm[["deviance"]]
res_pm[2,4] <- 100*(m4_pm_grade[["deviance"]] - m4_pm[["deviance"]])/m4_pm[["deviance"]]
res_pm[3,4] <- 100*(m4_pm_ethn[["deviance"]] - m4_pm[["deviance"]])/m4_pm[["deviance"]]
res_pm[4,4] <- 100*(m4_pm_sex[["deviance"]] - m4_pm[["deviance"]])/m4_pm[["deviance"]]

res_pm[1,5] <- 100*(m5_pm_factbrev[["deviance"]] - m5_pm[["deviance"]])/m5_pm[["deviance"]]
res_pm[2,5] <- 100*(m5_pm_grade[["deviance"]] - m5_pm[["deviance"]])/m5_pm[["deviance"]]
res_pm[3,5] <- 100*(m5_pm_ethn[["deviance"]] - m5_pm[["deviance"]])/m5_pm[["deviance"]]
res_pm[4,5] <- 100*(m5_pm_sex[["deviance"]] - m5_pm[["deviance"]])/m5_pm[["deviance"]]

res_pm_aic[1,1] <- 100*(m1_pm_factbrev[["aic"]] - m1_pm[["aic"]])/m1_pm[["aic"]]
res_pm_aic[2,1] <- 100*(m1_pm_grade[["aic"]] - m1_pm[["aic"]])/m1_pm[["aic"]]
res_pm_aic[3,1] <- 100*(m1_pm_ethn[["aic"]] - m1_pm[["aic"]])/m1_pm[["aic"]]
res_pm_aic[4,1] <- 100*(m1_pm_sex[["aic"]] - m1_pm[["aic"]])/m1_pm[["aic"]]

res_pm_aic[1,2] <- 100*(m2_pm_factbrev[["aic"]] - m2_pm[["aic"]])/m2_pm[["aic"]]
res_pm_aic[2,2] <- 100*(m2_pm_grade[["aic"]] - m2_pm[["aic"]])/m2_pm[["aic"]]
res_pm_aic[3,2] <- 100*(m2_pm_ethn[["aic"]] - m2_pm[["aic"]])/m2_pm[["aic"]]
res_pm_aic[4,2] <- 100*(m2_pm_sex[["aic"]] - m2_pm[["aic"]])/m2_pm[["aic"]]

res_pm_aic[1,3] <- 100*(m3_pm_factbrev[["aic"]] - m3_pm[["aic"]])/m3_pm[["aic"]]
res_pm_aic[2,3] <- 100*(m3_pm_grade[["aic"]] - m3_pm[["aic"]])/m3_pm[["aic"]]
res_pm_aic[3,3] <- 100*(m3_pm_ethn[["aic"]] - m3_pm[["aic"]])/m3_pm[["aic"]]
res_pm_aic[4,3] <- 100*(m3_pm_sex[["aic"]] - m3_pm[["aic"]])/m3_pm[["aic"]]

res_pm_aic[1,4] <- 100*(m4_pm_factbrev[["aic"]] - m4_pm[["aic"]])/m4_pm[["aic"]]
res_pm_aic[2,4] <- 100*(m4_pm_grade[["aic"]] - m4_pm[["aic"]])/m4_pm[["aic"]]
res_pm_aic[3,4] <- 100*(m4_pm_ethn[["aic"]] - m4_pm[["aic"]])/m4_pm[["aic"]]
res_pm_aic[4,4] <- 100*(m4_pm_sex[["aic"]] - m4_pm[["aic"]])/m4_pm[["aic"]]

res_pm_aic[1,5] <- 100*(m5_pm_factbrev[["aic"]] - m5_pm[["aic"]])/m5_pm[["aic"]]
res_pm_aic[2,5] <- 100*(m5_pm_grade[["aic"]] - m5_pm[["aic"]])/m5_pm[["aic"]]
res_pm_aic[3,5] <- 100*(m5_pm_ethn[["aic"]] - m5_pm[["aic"]])/m5_pm[["aic"]]
res_pm_aic[4,5] <- 100*(m5_pm_sex[["aic"]] - m5_pm[["aic"]])/m5_pm[["aic"]]

###

m1_dmz_factbrev <- glm(update(m1_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_1,family="binomial")
m1_dmz_grade <- glm(update(m1_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_1,family="binomial")
m1_dmz_ethn <- glm(update(m1_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_1,family="binomial")
m1_dmz_sex <- glm(update(m1_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_1,family="binomial")

m2_dmz_factbrev <- glm(update(m2_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_2,family="binomial")
m2_dmz_grade <- glm(update(m2_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_2,family="binomial")
m2_dmz_ethn <- glm(update(m2_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_2,family="binomial")
m2_dmz_sex <- glm(update(m2_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_2,family="binomial")

m3_dmz_factbrev <- glm(update(m3_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_3,family="binomial")
m3_dmz_grade <- glm(update(m3_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_3,family="binomial")
m3_dmz_ethn <- glm(update(m3_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_3,family="binomial")
m3_dmz_sex <- glm(update(m3_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_3,family="binomial")

m4_dmz_factbrev <- glm(update(m4_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_4,family="binomial")
m4_dmz_grade <- glm(update(m4_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_4,family="binomial")
m4_dmz_ethn <- glm(update(m4_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_4,family="binomial")
m4_dmz_sex <- glm(update(m4_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_4,family="binomial")

m5_dmz_factbrev <- glm(update(m5_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_5,family="binomial")
m5_dmz_grade <- glm(update(m5_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_5,family="binomial")
m5_dmz_ethn <- glm(update(m5_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_5,family="binomial")
m5_dmz_sex <- glm(update(m5_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_5,family="binomial")

res_dmz <- matrix(nrow = 4, ncol = 5)
rownames(res_dmz) <- c("factbrev","grade","ethn","sex")
colnames(res_dmz) <- c("v1","v2","v3","v4","v6")
res_dmz_aic <- res_dmz

res_dmz[1,1] <- 100*(m1_dmz_factbrev[["deviance"]] - m1_dmz[["deviance"]])/m1_dmz[["deviance"]]
res_dmz[2,1] <- 100*(m1_dmz_grade[["deviance"]] - m1_dmz[["deviance"]])/m1_dmz[["deviance"]]
res_dmz[3,1] <- 100*(m1_dmz_ethn[["deviance"]] - m1_dmz[["deviance"]])/m1_dmz[["deviance"]]
res_dmz[4,1] <- 100*(m1_dmz_sex[["deviance"]] - m1_dmz[["deviance"]])/m1_dmz[["deviance"]]

res_dmz[1,2] <- 100*(m2_dmz_factbrev[["deviance"]] - m2_dmz[["deviance"]])/m2_dmz[["deviance"]]
res_dmz[2,2] <- 100*(m2_dmz_grade[["deviance"]] - m2_dmz[["deviance"]])/m2_dmz[["deviance"]]
res_dmz[3,2] <- 100*(m2_dmz_ethn[["deviance"]] - m2_dmz[["deviance"]])/m2_dmz[["deviance"]]
res_dmz[4,2] <- 100*(m2_dmz_sex[["deviance"]] - m2_dmz[["deviance"]])/m2_dmz[["deviance"]]

res_dmz[1,3] <- 100*(m3_dmz_factbrev[["deviance"]] - m3_dmz[["deviance"]])/m3_dmz[["deviance"]]
res_dmz[2,3] <- 100*(m3_dmz_grade[["deviance"]] - m3_dmz[["deviance"]])/m3_dmz[["deviance"]]
res_dmz[3,3] <- 100*(m3_dmz_ethn[["deviance"]] - m3_dmz[["deviance"]])/m3_dmz[["deviance"]]
res_dmz[4,3] <- 100*(m3_dmz_sex[["deviance"]] - m3_dmz[["deviance"]])/m3_dmz[["deviance"]]

res_dmz[1,4] <- 100*(m4_dmz_factbrev[["deviance"]] - m4_dmz[["deviance"]])/m4_dmz[["deviance"]]
res_dmz[2,4] <- 100*(m4_dmz_grade[["deviance"]] - m4_dmz[["deviance"]])/m4_dmz[["deviance"]]
res_dmz[3,4] <- 100*(m4_dmz_ethn[["deviance"]] - m4_dmz[["deviance"]])/m4_dmz[["deviance"]]
res_dmz[4,4] <- 100*(m4_dmz_sex[["deviance"]] - m4_dmz[["deviance"]])/m4_dmz[["deviance"]]

res_dmz[1,5] <- 100*(m5_dmz_factbrev[["deviance"]] - m5_dmz[["deviance"]])/m5_dmz[["deviance"]]
res_dmz[2,5] <- 100*(m5_dmz_grade[["deviance"]] - m5_dmz[["deviance"]])/m5_dmz[["deviance"]]
res_dmz[3,5] <- 100*(m5_dmz_ethn[["deviance"]] - m5_dmz[["deviance"]])/m5_dmz[["deviance"]]
res_dmz[4,5] <- 100*(m5_dmz_sex[["deviance"]] - m5_dmz[["deviance"]])/m5_dmz[["deviance"]]

res_dmz_aic[1,1] <- 100*(m1_dmz_factbrev[["aic"]] - m1_dmz[["aic"]])/m1_dmz[["aic"]]
res_dmz_aic[2,1] <- 100*(m1_dmz_grade[["aic"]] - m1_dmz[["aic"]])/m1_dmz[["aic"]]
res_dmz_aic[3,1] <- 100*(m1_dmz_ethn[["aic"]] - m1_dmz[["aic"]])/m1_dmz[["aic"]]
res_dmz_aic[4,1] <- 100*(m1_dmz_sex[["aic"]] - m1_dmz[["aic"]])/m1_dmz[["aic"]]

res_dmz_aic[1,2] <- 100*(m2_dmz_factbrev[["aic"]] - m2_dmz[["aic"]])/m2_dmz[["aic"]]
res_dmz_aic[2,2] <- 100*(m2_dmz_grade[["aic"]] - m2_dmz[["aic"]])/m2_dmz[["aic"]]
res_dmz_aic[3,2] <- 100*(m2_dmz_ethn[["aic"]] - m2_dmz[["aic"]])/m2_dmz[["aic"]]
res_dmz_aic[4,2] <- 100*(m2_dmz_sex[["aic"]] - m2_dmz[["aic"]])/m2_dmz[["aic"]]

res_dmz_aic[1,3] <- 100*(m3_dmz_factbrev[["aic"]] - m3_dmz[["aic"]])/m3_dmz[["aic"]]
res_dmz_aic[2,3] <- 100*(m3_dmz_grade[["aic"]] - m3_dmz[["aic"]])/m3_dmz[["aic"]]
res_dmz_aic[3,3] <- 100*(m3_dmz_ethn[["aic"]] - m3_dmz[["aic"]])/m3_dmz[["aic"]]
res_dmz_aic[4,3] <- 100*(m3_dmz_sex[["aic"]] - m3_dmz[["aic"]])/m3_dmz[["aic"]]

res_dmz_aic[1,4] <- 100*(m4_dmz_factbrev[["aic"]] - m4_dmz[["aic"]])/m4_dmz[["aic"]]
res_dmz_aic[2,4] <- 100*(m4_dmz_grade[["aic"]] - m4_dmz[["aic"]])/m4_dmz[["aic"]]
res_dmz_aic[3,4] <- 100*(m4_dmz_ethn[["aic"]] - m4_dmz[["aic"]])/m4_dmz[["aic"]]
res_dmz_aic[4,4] <- 100*(m4_dmz_sex[["aic"]] - m4_dmz[["aic"]])/m4_dmz[["aic"]]

res_dmz_aic[1,5] <- 100*(m5_dmz_factbrev[["aic"]] - m5_dmz[["aic"]])/m5_dmz[["aic"]]
res_dmz_aic[2,5] <- 100*(m5_dmz_grade[["aic"]] - m5_dmz[["aic"]])/m5_dmz[["aic"]]
res_dmz_aic[3,5] <- 100*(m5_dmz_ethn[["aic"]] - m5_dmz[["aic"]])/m5_dmz[["aic"]]
res_dmz_aic[4,5] <- 100*(m5_dmz_sex[["aic"]] - m5_dmz[["aic"]])/m5_dmz[["aic"]]

##### comparaison des tailles d effets - estimation allf #####

n1_rvl_factbrev <- glm(update(n1_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_1,family="binomial")
n1_rvl_grade <- glm(update(n1_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_1,family="binomial")
n1_rvl_ethn <- glm(update(n1_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_1,family="binomial")
n1_rvl_sex <- glm(update(n1_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_1,family="binomial")

n2_rvl_factbrev <- glm(update(n2_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_2,family="binomial")
n2_rvl_grade <- glm(update(n2_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_2,family="binomial")
n2_rvl_ethn <- glm(update(n2_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_2,family="binomial")
n2_rvl_sex <- glm(update(n2_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_2,family="binomial")

n3_rvl_factbrev <- glm(update(n3_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_3,family="binomial")
n3_rvl_grade <- glm(update(n3_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_3,family="binomial")
n3_rvl_ethn <- glm(update(n3_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_3,family="binomial")
n3_rvl_sex <- glm(update(n3_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_3,family="binomial")

n4_rvl_factbrev <- glm(update(n4_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_4,family="binomial")
n4_rvl_grade <- glm(update(n4_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_4,family="binomial")
n4_rvl_ethn <- glm(update(n4_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_4,family="binomial")
n4_rvl_sex <- glm(update(n4_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_4,family="binomial")

n5_rvl_factbrev <- glm(update(n5_rvl$formula, . ~ .-factbrev_absdiff),
                       data = drvl_5,family="binomial")
n5_rvl_grade <- glm(update(n5_rvl$formula, . ~ .-grade_absdiff),
                    data = drvl_5,family="binomial")
n5_rvl_ethn <- glm(update(n5_rvl$formula, . ~ .-ethn_nodematch),
                   data = drvl_5,family="binomial")
n5_rvl_sex <- glm(update(n5_rvl$formula, . ~ .-sex_nodematch),
                  data = drvl_5,family="binomial")

res2_rvl <- matrix(nrow = 4, ncol = 5)
rownames(res2_rvl) <- c("factbrev","grade","ethn","sex")
colnames(res2_rvl) <- c("v1","v2","v3","v4","v6")
res2_rvl_aic <- res2_rvl

res2_rvl[1,1] <- 100*(n1_rvl_factbrev[["deviance"]] - n1_rvl[["deviance"]])/n1_rvl[["deviance"]]
res2_rvl[2,1] <- 100*(n1_rvl_grade[["deviance"]] - n1_rvl[["deviance"]])/n1_rvl[["deviance"]]
res2_rvl[3,1] <- 100*(n1_rvl_ethn[["deviance"]] - n1_rvl[["deviance"]])/n1_rvl[["deviance"]]
res2_rvl[4,1] <- 100*(n1_rvl_sex[["deviance"]] - n1_rvl[["deviance"]])/n1_rvl[["deviance"]]

res2_rvl[1,2] <- 100*(n2_rvl_factbrev[["deviance"]] - n2_rvl[["deviance"]])/n2_rvl[["deviance"]]
res2_rvl[2,2] <- 100*(n2_rvl_grade[["deviance"]] - n2_rvl[["deviance"]])/n2_rvl[["deviance"]]
res2_rvl[3,2] <- 100*(n2_rvl_ethn[["deviance"]] - n2_rvl[["deviance"]])/n2_rvl[["deviance"]]
res2_rvl[4,2] <- 100*(n2_rvl_sex[["deviance"]] - n2_rvl[["deviance"]])/n2_rvl[["deviance"]]

res2_rvl[1,3] <- 100*(n3_rvl_factbrev[["deviance"]] - n3_rvl[["deviance"]])/n3_rvl[["deviance"]]
res2_rvl[2,3] <- 100*(n3_rvl_grade[["deviance"]] - n3_rvl[["deviance"]])/n3_rvl[["deviance"]]
res2_rvl[3,3] <- 100*(n3_rvl_ethn[["deviance"]] - n3_rvl[["deviance"]])/n3_rvl[["deviance"]]
res2_rvl[4,3] <- 100*(n3_rvl_sex[["deviance"]] - n3_rvl[["deviance"]])/n3_rvl[["deviance"]]

res2_rvl[1,4] <- 100*(n4_rvl_factbrev[["deviance"]] - n4_rvl[["deviance"]])/n4_rvl[["deviance"]]
res2_rvl[2,4] <- 100*(n4_rvl_grade[["deviance"]] - n4_rvl[["deviance"]])/n4_rvl[["deviance"]]
res2_rvl[3,4] <- 100*(n4_rvl_ethn[["deviance"]] - n4_rvl[["deviance"]])/n4_rvl[["deviance"]]
res2_rvl[4,4] <- 100*(n4_rvl_sex[["deviance"]] - n4_rvl[["deviance"]])/n4_rvl[["deviance"]]

res2_rvl[1,5] <- 100*(n5_rvl_factbrev[["deviance"]] - n5_rvl[["deviance"]])/n5_rvl[["deviance"]]
res2_rvl[2,5] <- 100*(n5_rvl_grade[["deviance"]] - n5_rvl[["deviance"]])/n5_rvl[["deviance"]]
res2_rvl[3,5] <- 100*(n5_rvl_ethn[["deviance"]] - n5_rvl[["deviance"]])/n5_rvl[["deviance"]]
res2_rvl[4,5] <- 100*(n5_rvl_sex[["deviance"]] - n5_rvl[["deviance"]])/n5_rvl[["deviance"]]

res2_rvl_aic[1,1] <- 100*(n1_rvl_factbrev[["aic"]] - n1_rvl[["aic"]])/n1_rvl[["aic"]]
res2_rvl_aic[2,1] <- 100*(n1_rvl_grade[["aic"]] - n1_rvl[["aic"]])/n1_rvl[["aic"]]
res2_rvl_aic[3,1] <- 100*(n1_rvl_ethn[["aic"]] - n1_rvl[["aic"]])/n1_rvl[["aic"]]
res2_rvl_aic[4,1] <- 100*(n1_rvl_sex[["aic"]] - n1_rvl[["aic"]])/n1_rvl[["aic"]]

res2_rvl_aic[1,2] <- 100*(n2_rvl_factbrev[["aic"]] - n2_rvl[["aic"]])/n2_rvl[["aic"]]
res2_rvl_aic[2,2] <- 100*(n2_rvl_grade[["aic"]] - n2_rvl[["aic"]])/n2_rvl[["aic"]]
res2_rvl_aic[3,2] <- 100*(n2_rvl_ethn[["aic"]] - n2_rvl[["aic"]])/n2_rvl[["aic"]]
res2_rvl_aic[4,2] <- 100*(n2_rvl_sex[["aic"]] - n2_rvl[["aic"]])/n2_rvl[["aic"]]

res2_rvl_aic[1,3] <- 100*(n3_rvl_factbrev[["aic"]] - n3_rvl[["aic"]])/n3_rvl[["aic"]]
res2_rvl_aic[2,3] <- 100*(n3_rvl_grade[["aic"]] - n3_rvl[["aic"]])/n3_rvl[["aic"]]
res2_rvl_aic[3,3] <- 100*(n3_rvl_ethn[["aic"]] - n3_rvl[["aic"]])/n3_rvl[["aic"]]
res2_rvl_aic[4,3] <- 100*(n3_rvl_sex[["aic"]] - n3_rvl[["aic"]])/n3_rvl[["aic"]]

res2_rvl_aic[1,4] <- 100*(n4_rvl_factbrev[["aic"]] - n4_rvl[["aic"]])/n4_rvl[["aic"]]
res2_rvl_aic[2,4] <- 100*(n4_rvl_grade[["aic"]] - n4_rvl[["aic"]])/n4_rvl[["aic"]]
res2_rvl_aic[3,4] <- 100*(n4_rvl_ethn[["aic"]] - n4_rvl[["aic"]])/n4_rvl[["aic"]]
res2_rvl_aic[4,4] <- 100*(n4_rvl_sex[["aic"]] - n4_rvl[["aic"]])/n4_rvl[["aic"]]

res2_rvl_aic[1,5] <- 100*(n5_rvl_factbrev[["aic"]] - n5_rvl[["aic"]])/n5_rvl[["aic"]]
res2_rvl_aic[2,5] <- 100*(n5_rvl_grade[["aic"]] - n5_rvl[["aic"]])/n5_rvl[["aic"]]
res2_rvl_aic[3,5] <- 100*(n5_rvl_ethn[["aic"]] - n5_rvl[["aic"]])/n5_rvl[["aic"]]
res2_rvl_aic[4,5] <- 100*(n5_rvl_sex[["aic"]] - n5_rvl[["aic"]])/n5_rvl[["aic"]]

###

n1_jbs_factbrev <- glm(update(n1_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_1,family="binomial")
n1_jbs_grade <- glm(update(n1_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_1,family="binomial")
n1_jbs_ethn <- glm(update(n1_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_1,family="binomial")
n1_jbs_sex <- glm(update(n1_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_1,family="binomial")

n2_jbs_factbrev <- glm(update(n2_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_2,family="binomial")
n2_jbs_grade <- glm(update(n2_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_2,family="binomial")
n2_jbs_ethn <- glm(update(n2_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_2,family="binomial")
n2_jbs_sex <- glm(update(n2_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_2,family="binomial")

n3_jbs_factbrev <- glm(update(n3_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_3,family="binomial")
n3_jbs_grade <- glm(update(n3_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_3,family="binomial")
n3_jbs_ethn <- glm(update(n3_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_3,family="binomial")
n3_jbs_sex <- glm(update(n3_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_3,family="binomial")

n4_jbs_factbrev <- glm(update(n4_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_4,family="binomial")
n4_jbs_grade <- glm(update(n4_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_4,family="binomial")
n4_jbs_ethn <- glm(update(n4_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_4,family="binomial")
n4_jbs_sex <- glm(update(n4_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_4,family="binomial")

n5_jbs_factbrev <- glm(update(n5_jbs$formula, . ~ .-factbrev_absdiff),
                       data = djbs_5,family="binomial")
n5_jbs_grade <- glm(update(n5_jbs$formula, . ~ .-grade_absdiff),
                    data = djbs_5,family="binomial")
n5_jbs_ethn <- glm(update(n5_jbs$formula, . ~ .-ethn_nodematch),
                   data = djbs_5,family="binomial")
n5_jbs_sex <- glm(update(n5_jbs$formula, . ~ .-sex_nodematch),
                  data = djbs_5,family="binomial")

res2_jbs <- matrix(nrow = 4, ncol = 5)
rownames(res2_jbs) <- c("factbrev","grade","ethn","sex")
colnames(res2_jbs) <- c("v1","v2","v3","v4","v6")
res2_jbs_aic <- res2_jbs

res2_jbs[1,1] <- 100*(n1_jbs_factbrev[["deviance"]] - n1_jbs[["deviance"]])/n1_jbs[["deviance"]]
res2_jbs[2,1] <- 100*(n1_jbs_grade[["deviance"]] - n1_jbs[["deviance"]])/n1_jbs[["deviance"]]
res2_jbs[3,1] <- 100*(n1_jbs_ethn[["deviance"]] - n1_jbs[["deviance"]])/n1_jbs[["deviance"]]
res2_jbs[4,1] <- 100*(n1_jbs_sex[["deviance"]] - n1_jbs[["deviance"]])/n1_jbs[["deviance"]]

res2_jbs[1,2] <- 100*(n2_jbs_factbrev[["deviance"]] - n2_jbs[["deviance"]])/n2_jbs[["deviance"]]
res2_jbs[2,2] <- 100*(n2_jbs_grade[["deviance"]] - n2_jbs[["deviance"]])/n2_jbs[["deviance"]]
res2_jbs[3,2] <- 100*(n2_jbs_ethn[["deviance"]] - n2_jbs[["deviance"]])/n2_jbs[["deviance"]]
res2_jbs[4,2] <- 100*(n2_jbs_sex[["deviance"]] - n2_jbs[["deviance"]])/n2_jbs[["deviance"]]

res2_jbs[1,3] <- 100*(n3_jbs_factbrev[["deviance"]] - n3_jbs[["deviance"]])/n3_jbs[["deviance"]]
res2_jbs[2,3] <- 100*(n3_jbs_grade[["deviance"]] - n3_jbs[["deviance"]])/n3_jbs[["deviance"]]
res2_jbs[3,3] <- 100*(n3_jbs_ethn[["deviance"]] - n3_jbs[["deviance"]])/n3_jbs[["deviance"]]
res2_jbs[4,3] <- 100*(n3_jbs_sex[["deviance"]] - n3_jbs[["deviance"]])/n3_jbs[["deviance"]]

res2_jbs[1,4] <- 100*(n4_jbs_factbrev[["deviance"]] - n4_jbs[["deviance"]])/n4_jbs[["deviance"]]
res2_jbs[2,4] <- 100*(n4_jbs_grade[["deviance"]] - n4_jbs[["deviance"]])/n4_jbs[["deviance"]]
res2_jbs[3,4] <- 100*(n4_jbs_ethn[["deviance"]] - n4_jbs[["deviance"]])/n4_jbs[["deviance"]]
res2_jbs[4,4] <- 100*(n4_jbs_sex[["deviance"]] - n4_jbs[["deviance"]])/n4_jbs[["deviance"]]

res2_jbs[1,5] <- 100*(n5_jbs_factbrev[["deviance"]] - n5_jbs[["deviance"]])/n5_jbs[["deviance"]]
res2_jbs[2,5] <- 100*(n5_jbs_grade[["deviance"]] - n5_jbs[["deviance"]])/n5_jbs[["deviance"]]
res2_jbs[3,5] <- 100*(n5_jbs_ethn[["deviance"]] - n5_jbs[["deviance"]])/n5_jbs[["deviance"]]
res2_jbs[4,5] <- 100*(n5_jbs_sex[["deviance"]] - n5_jbs[["deviance"]])/n5_jbs[["deviance"]]

res2_jbs_aic[1,1] <- 100*(n1_jbs_factbrev[["aic"]] - n1_jbs[["aic"]])/n1_jbs[["aic"]]
res2_jbs_aic[2,1] <- 100*(n1_jbs_grade[["aic"]] - n1_jbs[["aic"]])/n1_jbs[["aic"]]
res2_jbs_aic[3,1] <- 100*(n1_jbs_ethn[["aic"]] - n1_jbs[["aic"]])/n1_jbs[["aic"]]
res2_jbs_aic[4,1] <- 100*(n1_jbs_sex[["aic"]] - n1_jbs[["aic"]])/n1_jbs[["aic"]]

res2_jbs_aic[1,2] <- 100*(n2_jbs_factbrev[["aic"]] - n2_jbs[["aic"]])/n2_jbs[["aic"]]
res2_jbs_aic[2,2] <- 100*(n2_jbs_grade[["aic"]] - n2_jbs[["aic"]])/n2_jbs[["aic"]]
res2_jbs_aic[3,2] <- 100*(n2_jbs_ethn[["aic"]] - n2_jbs[["aic"]])/n2_jbs[["aic"]]
res2_jbs_aic[4,2] <- 100*(n2_jbs_sex[["aic"]] - n2_jbs[["aic"]])/n2_jbs[["aic"]]

res2_jbs_aic[1,3] <- 100*(n3_jbs_factbrev[["aic"]] - n3_jbs[["aic"]])/n3_jbs[["aic"]]
res2_jbs_aic[2,3] <- 100*(n3_jbs_grade[["aic"]] - n3_jbs[["aic"]])/n3_jbs[["aic"]]
res2_jbs_aic[3,3] <- 100*(n3_jbs_ethn[["aic"]] - n3_jbs[["aic"]])/n3_jbs[["aic"]]
res2_jbs_aic[4,3] <- 100*(n3_jbs_sex[["aic"]] - n3_jbs[["aic"]])/n3_jbs[["aic"]]

res2_jbs_aic[1,4] <- 100*(n4_jbs_factbrev[["aic"]] - n4_jbs[["aic"]])/n4_jbs[["aic"]]
res2_jbs_aic[2,4] <- 100*(n4_jbs_grade[["aic"]] - n4_jbs[["aic"]])/n4_jbs[["aic"]]
res2_jbs_aic[3,4] <- 100*(n4_jbs_ethn[["aic"]] - n4_jbs[["aic"]])/n4_jbs[["aic"]]
res2_jbs_aic[4,4] <- 100*(n4_jbs_sex[["aic"]] - n4_jbs[["aic"]])/n4_jbs[["aic"]]

res2_jbs_aic[1,5] <- 100*(n5_jbs_factbrev[["aic"]] - n5_jbs[["aic"]])/n5_jbs[["aic"]]
res2_jbs_aic[2,5] <- 100*(n5_jbs_grade[["aic"]] - n5_jbs[["aic"]])/n5_jbs[["aic"]]
res2_jbs_aic[3,5] <- 100*(n5_jbs_ethn[["aic"]] - n5_jbs[["aic"]])/n5_jbs[["aic"]]
res2_jbs_aic[4,5] <- 100*(n5_jbs_sex[["aic"]] - n5_jbs[["aic"]])/n5_jbs[["aic"]]

###

n1_pm_factbrev <- glm(update(n1_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_1,family="binomial")
n1_pm_grade <- glm(update(n1_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_1,family="binomial")
n1_pm_ethn <- glm(update(n1_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_1,family="binomial")
n1_pm_sex <- glm(update(n1_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_1,family="binomial")

n2_pm_factbrev <- glm(update(n2_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_2,family="binomial")
n2_pm_grade <- glm(update(n2_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_2,family="binomial")
n2_pm_ethn <- glm(update(n2_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_2,family="binomial")
n2_pm_sex <- glm(update(n2_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_2,family="binomial")

n3_pm_factbrev <- glm(update(n3_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_3,family="binomial")
n3_pm_grade <- glm(update(n3_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_3,family="binomial")
n3_pm_ethn <- glm(update(n3_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_3,family="binomial")
n3_pm_sex <- glm(update(n3_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_3,family="binomial")

n4_pm_factbrev <- glm(update(n4_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_4,family="binomial")
n4_pm_grade <- glm(update(n4_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_4,family="binomial")
n4_pm_ethn <- glm(update(n4_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_4,family="binomial")
n4_pm_sex <- glm(update(n4_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_4,family="binomial")

n5_pm_factbrev <- glm(update(n5_pm$formula, . ~ .-factbrev_absdiff),
                      data = dpm_5,family="binomial")
n5_pm_grade <- glm(update(n5_pm$formula, . ~ .-grade_absdiff),
                   data = dpm_5,family="binomial")
n5_pm_ethn <- glm(update(n5_pm$formula, . ~ .-ethn_nodematch),
                  data = dpm_5,family="binomial")
n5_pm_sex <- glm(update(n5_pm$formula, . ~ .-sex_nodematch),
                 data = dpm_5,family="binomial")

res2_pm <- matrix(nrow = 4, ncol = 5)
rownames(res2_pm) <- c("factbrev","grade","ethn","sex")
colnames(res2_pm) <- c("v1","v2","v3","v4","v6")
res2_pm_aic <- res2_pm

res2_pm[1,1] <- 100*(n1_pm_factbrev[["deviance"]] - n1_pm[["deviance"]])/n1_pm[["deviance"]]
res2_pm[2,1] <- 100*(n1_pm_grade[["deviance"]] - n1_pm[["deviance"]])/n1_pm[["deviance"]]
res2_pm[3,1] <- 100*(n1_pm_ethn[["deviance"]] - n1_pm[["deviance"]])/n1_pm[["deviance"]]
res2_pm[4,1] <- 100*(n1_pm_sex[["deviance"]] - n1_pm[["deviance"]])/n1_pm[["deviance"]]

res2_pm[1,2] <- 100*(n2_pm_factbrev[["deviance"]] - n2_pm[["deviance"]])/n2_pm[["deviance"]]
res2_pm[2,2] <- 100*(n2_pm_grade[["deviance"]] - n2_pm[["deviance"]])/n2_pm[["deviance"]]
res2_pm[3,2] <- 100*(n2_pm_ethn[["deviance"]] - n2_pm[["deviance"]])/n2_pm[["deviance"]]
res2_pm[4,2] <- 100*(n2_pm_sex[["deviance"]] - n2_pm[["deviance"]])/n2_pm[["deviance"]]

res2_pm[1,3] <- 100*(n3_pm_factbrev[["deviance"]] - n3_pm[["deviance"]])/n3_pm[["deviance"]]
res2_pm[2,3] <- 100*(n3_pm_grade[["deviance"]] - n3_pm[["deviance"]])/n3_pm[["deviance"]]
res2_pm[3,3] <- 100*(n3_pm_ethn[["deviance"]] - n3_pm[["deviance"]])/n3_pm[["deviance"]]
res2_pm[4,3] <- 100*(n3_pm_sex[["deviance"]] - n3_pm[["deviance"]])/n3_pm[["deviance"]]

res2_pm[1,4] <- 100*(n4_pm_factbrev[["deviance"]] - n4_pm[["deviance"]])/n4_pm[["deviance"]]
res2_pm[2,4] <- 100*(n4_pm_grade[["deviance"]] - n4_pm[["deviance"]])/n4_pm[["deviance"]]
res2_pm[3,4] <- 100*(n4_pm_ethn[["deviance"]] - n4_pm[["deviance"]])/n4_pm[["deviance"]]
res2_pm[4,4] <- 100*(n4_pm_sex[["deviance"]] - n4_pm[["deviance"]])/n4_pm[["deviance"]]

res2_pm[1,5] <- 100*(n5_pm_factbrev[["deviance"]] - n5_pm[["deviance"]])/n5_pm[["deviance"]]
res2_pm[2,5] <- 100*(n5_pm_grade[["deviance"]] - n5_pm[["deviance"]])/n5_pm[["deviance"]]
res2_pm[3,5] <- 100*(n5_pm_ethn[["deviance"]] - n5_pm[["deviance"]])/n5_pm[["deviance"]]
res2_pm[4,5] <- 100*(n5_pm_sex[["deviance"]] - n5_pm[["deviance"]])/n5_pm[["deviance"]]

res2_pm_aic[1,1] <- 100*(n1_pm_factbrev[["aic"]] - n1_pm[["aic"]])/n1_pm[["aic"]]
res2_pm_aic[2,1] <- 100*(n1_pm_grade[["aic"]] - n1_pm[["aic"]])/n1_pm[["aic"]]
res2_pm_aic[3,1] <- 100*(n1_pm_ethn[["aic"]] - n1_pm[["aic"]])/n1_pm[["aic"]]
res2_pm_aic[4,1] <- 100*(n1_pm_sex[["aic"]] - n1_pm[["aic"]])/n1_pm[["aic"]]

res2_pm_aic[1,2] <- 100*(n2_pm_factbrev[["aic"]] - n2_pm[["aic"]])/n2_pm[["aic"]]
res2_pm_aic[2,2] <- 100*(n2_pm_grade[["aic"]] - n2_pm[["aic"]])/n2_pm[["aic"]]
res2_pm_aic[3,2] <- 100*(n2_pm_ethn[["aic"]] - n2_pm[["aic"]])/n2_pm[["aic"]]
res2_pm_aic[4,2] <- 100*(n2_pm_sex[["aic"]] - n2_pm[["aic"]])/n2_pm[["aic"]]

res2_pm_aic[1,3] <- 100*(n3_pm_factbrev[["aic"]] - n3_pm[["aic"]])/n3_pm[["aic"]]
res2_pm_aic[2,3] <- 100*(n3_pm_grade[["aic"]] - n3_pm[["aic"]])/n3_pm[["aic"]]
res2_pm_aic[3,3] <- 100*(n3_pm_ethn[["aic"]] - n3_pm[["aic"]])/n3_pm[["aic"]]
res2_pm_aic[4,3] <- 100*(n3_pm_sex[["aic"]] - n3_pm[["aic"]])/n3_pm[["aic"]]

res2_pm_aic[1,4] <- 100*(n4_pm_factbrev[["aic"]] - n4_pm[["aic"]])/n4_pm[["aic"]]
res2_pm_aic[2,4] <- 100*(n4_pm_grade[["aic"]] - n4_pm[["aic"]])/n4_pm[["aic"]]
res2_pm_aic[3,4] <- 100*(n4_pm_ethn[["aic"]] - n4_pm[["aic"]])/n4_pm[["aic"]]
res2_pm_aic[4,4] <- 100*(n4_pm_sex[["aic"]] - n4_pm[["aic"]])/n4_pm[["aic"]]

res2_pm_aic[1,5] <- 100*(n5_pm_factbrev[["aic"]] - n5_pm[["aic"]])/n5_pm[["aic"]]
res2_pm_aic[2,5] <- 100*(n5_pm_grade[["aic"]] - n5_pm[["aic"]])/n5_pm[["aic"]]
res2_pm_aic[3,5] <- 100*(n5_pm_ethn[["aic"]] - n5_pm[["aic"]])/n5_pm[["aic"]]
res2_pm_aic[4,5] <- 100*(n5_pm_sex[["aic"]] - n5_pm[["aic"]])/n5_pm[["aic"]]

###

n1_dmz_factbrev <- glm(update(n1_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_1,family="binomial")
n1_dmz_grade <- glm(update(n1_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_1,family="binomial")
n1_dmz_ethn <- glm(update(n1_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_1,family="binomial")
n1_dmz_sex <- glm(update(n1_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_1,family="binomial")

n2_dmz_factbrev <- glm(update(n2_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_2,family="binomial")
n2_dmz_grade <- glm(update(n2_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_2,family="binomial")
n2_dmz_ethn <- glm(update(n2_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_2,family="binomial")
n2_dmz_sex <- glm(update(n2_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_2,family="binomial")

n3_dmz_factbrev <- glm(update(n3_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_3,family="binomial")
n3_dmz_grade <- glm(update(n3_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_3,family="binomial")
n3_dmz_ethn <- glm(update(n3_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_3,family="binomial")
n3_dmz_sex <- glm(update(n3_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_3,family="binomial")

n4_dmz_factbrev <- glm(update(n4_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_4,family="binomial")
n4_dmz_grade <- glm(update(n4_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_4,family="binomial")
n4_dmz_ethn <- glm(update(n4_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_4,family="binomial")
n4_dmz_sex <- glm(update(n4_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_4,family="binomial")

n5_dmz_factbrev <- glm(update(n5_dmz$formula, . ~ .-factbrev_absdiff),
                       data = ddmz_5,family="binomial")
n5_dmz_grade <- glm(update(n5_dmz$formula, . ~ .-grade_absdiff),
                    data = ddmz_5,family="binomial")
n5_dmz_ethn <- glm(update(n5_dmz$formula, . ~ .-ethn_nodematch),
                   data = ddmz_5,family="binomial")
n5_dmz_sex <- glm(update(n5_dmz$formula, . ~ .-sex_nodematch),
                  data = ddmz_5,family="binomial")

res2_dmz <- matrix(nrow = 4, ncol = 5)
rownames(res2_dmz) <- c("factbrev","grade","ethn","sex")
colnames(res2_dmz) <- c("v1","v2","v3","v4","v6")
res2_dmz_aic <- res2_dmz

res2_dmz[1,1] <- 100*(n1_dmz_factbrev[["deviance"]] - n1_dmz[["deviance"]])/n1_dmz[["deviance"]]
res2_dmz[2,1] <- 100*(n1_dmz_grade[["deviance"]] - n1_dmz[["deviance"]])/n1_dmz[["deviance"]]
res2_dmz[3,1] <- 100*(n1_dmz_ethn[["deviance"]] - n1_dmz[["deviance"]])/n1_dmz[["deviance"]]
res2_dmz[4,1] <- 100*(n1_dmz_sex[["deviance"]] - n1_dmz[["deviance"]])/n1_dmz[["deviance"]]

res2_dmz[1,2] <- 100*(n2_dmz_factbrev[["deviance"]] - n2_dmz[["deviance"]])/n2_dmz[["deviance"]]
res2_dmz[2,2] <- 100*(n2_dmz_grade[["deviance"]] - n2_dmz[["deviance"]])/n2_dmz[["deviance"]]
res2_dmz[3,2] <- 100*(n2_dmz_ethn[["deviance"]] - n2_dmz[["deviance"]])/n2_dmz[["deviance"]]
res2_dmz[4,2] <- 100*(n2_dmz_sex[["deviance"]] - n2_dmz[["deviance"]])/n2_dmz[["deviance"]]

res2_dmz[1,3] <- 100*(n3_dmz_factbrev[["deviance"]] - n3_dmz[["deviance"]])/n3_dmz[["deviance"]]
res2_dmz[2,3] <- 100*(n3_dmz_grade[["deviance"]] - n3_dmz[["deviance"]])/n3_dmz[["deviance"]]
res2_dmz[3,3] <- 100*(n3_dmz_ethn[["deviance"]] - n3_dmz[["deviance"]])/n3_dmz[["deviance"]]
res2_dmz[4,3] <- 100*(n3_dmz_sex[["deviance"]] - n3_dmz[["deviance"]])/n3_dmz[["deviance"]]

res2_dmz[1,4] <- 100*(n4_dmz_factbrev[["deviance"]] - n4_dmz[["deviance"]])/n4_dmz[["deviance"]]
res2_dmz[2,4] <- 100*(n4_dmz_grade[["deviance"]] - n4_dmz[["deviance"]])/n4_dmz[["deviance"]]
res2_dmz[3,4] <- 100*(n4_dmz_ethn[["deviance"]] - n4_dmz[["deviance"]])/n4_dmz[["deviance"]]
res2_dmz[4,4] <- 100*(n4_dmz_sex[["deviance"]] - n4_dmz[["deviance"]])/n4_dmz[["deviance"]]

res2_dmz[1,5] <- 100*(n5_dmz_factbrev[["deviance"]] - n5_dmz[["deviance"]])/n5_dmz[["deviance"]]
res2_dmz[2,5] <- 100*(n5_dmz_grade[["deviance"]] - n5_dmz[["deviance"]])/n5_dmz[["deviance"]]
res2_dmz[3,5] <- 100*(n5_dmz_ethn[["deviance"]] - n5_dmz[["deviance"]])/n5_dmz[["deviance"]]
res2_dmz[4,5] <- 100*(n5_dmz_sex[["deviance"]] - n5_dmz[["deviance"]])/n5_dmz[["deviance"]]

res2_dmz_aic[1,1] <- 100*(n1_dmz_factbrev[["aic"]] - n1_dmz[["aic"]])/n1_dmz[["aic"]]
res2_dmz_aic[2,1] <- 100*(n1_dmz_grade[["aic"]] - n1_dmz[["aic"]])/n1_dmz[["aic"]]
res2_dmz_aic[3,1] <- 100*(n1_dmz_ethn[["aic"]] - n1_dmz[["aic"]])/n1_dmz[["aic"]]
res2_dmz_aic[4,1] <- 100*(n1_dmz_sex[["aic"]] - n1_dmz[["aic"]])/n1_dmz[["aic"]]

res2_dmz_aic[1,2] <- 100*(n2_dmz_factbrev[["aic"]] - n2_dmz[["aic"]])/n2_dmz[["aic"]]
res2_dmz_aic[2,2] <- 100*(n2_dmz_grade[["aic"]] - n2_dmz[["aic"]])/n2_dmz[["aic"]]
res2_dmz_aic[3,2] <- 100*(n2_dmz_ethn[["aic"]] - n2_dmz[["aic"]])/n2_dmz[["aic"]]
res2_dmz_aic[4,2] <- 100*(n2_dmz_sex[["aic"]] - n2_dmz[["aic"]])/n2_dmz[["aic"]]

res2_dmz_aic[1,3] <- 100*(n3_dmz_factbrev[["aic"]] - n3_dmz[["aic"]])/n3_dmz[["aic"]]
res2_dmz_aic[2,3] <- 100*(n3_dmz_grade[["aic"]] - n3_dmz[["aic"]])/n3_dmz[["aic"]]
res2_dmz_aic[3,3] <- 100*(n3_dmz_ethn[["aic"]] - n3_dmz[["aic"]])/n3_dmz[["aic"]]
res2_dmz_aic[4,3] <- 100*(n3_dmz_sex[["aic"]] - n3_dmz[["aic"]])/n3_dmz[["aic"]]

res2_dmz_aic[1,4] <- 100*(n4_dmz_factbrev[["aic"]] - n4_dmz[["aic"]])/n4_dmz[["aic"]]
res2_dmz_aic[2,4] <- 100*(n4_dmz_grade[["aic"]] - n4_dmz[["aic"]])/n4_dmz[["aic"]]
res2_dmz_aic[3,4] <- 100*(n4_dmz_ethn[["aic"]] - n4_dmz[["aic"]])/n4_dmz[["aic"]]
res2_dmz_aic[4,4] <- 100*(n4_dmz_sex[["aic"]] - n4_dmz[["aic"]])/n4_dmz[["aic"]]

res2_dmz_aic[1,5] <- 100*(n5_dmz_factbrev[["aic"]] - n5_dmz[["aic"]])/n5_dmz[["aic"]]
res2_dmz_aic[2,5] <- 100*(n5_dmz_grade[["aic"]] - n5_dmz[["aic"]])/n5_dmz[["aic"]]
res2_dmz_aic[3,5] <- 100*(n5_dmz_ethn[["aic"]] - n5_dmz[["aic"]])/n5_dmz[["aic"]]
res2_dmz_aic[4,5] <- 100*(n5_dmz_sex[["aic"]] - n5_dmz[["aic"]])/n5_dmz[["aic"]]

##### comparaison des tailles d effets - lecture gf #####

# Exprimes en points de pourcentage (1.2 = 1.2%)
res_rvl
res_jbs
res_pm
res_dmz

# On verifie que c est pareil avec aic a la place de deviance residuelle.
res_rvl
res_rvl_aic
res_jbs
res_jbs_aic
res_pm
res_pm_aic
res_dmz
res_dmz_aic

library(reshape2)

aa <- data.frame(res_rvl)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Paris 1 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

aa <- data.frame(res_jbs[-which(rownames(res_jbs)=="grade"),])
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00"),
                     breaks = c("sex","factbrev","ethn"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17),
                     breaks = c("sex","factbrev","ethn"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire")) +
  ggtitle("Paris 2 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))


aa <- data.frame(res_pm)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Savoie 1 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

aa <- data.frame(res_dmz)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Savoie 2 - Très Bons Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

##### comparaison des tailles d effets - lecture allf #####

# Exprimes en points de pourcentage (1.2 = 1.2%)
res2_rvl
res2_jbs
res2_pm
res2_dmz

# On verifie que c est pareil avec aic a la place de deviance residuelle.
res2_rvl
res2_rvl_aic
res2_jbs
res2_jbs_aic
res2_pm
res2_pm_aic
res2_dmz
res2_dmz_aic

library(reshape2)

aa <- data.frame(res2_rvl)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Paris 1 - Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

aa <- data.frame(res2_jbs[-which(rownames(res2_jbs)=="grade"),])
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00"),
                     breaks = c("sex","factbrev","ethn"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17),
                     breaks = c("sex","factbrev","ethn"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire")) +
  ggtitle("Paris 2 - Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

aa <- data.frame(res2_pm)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Savoie 1 - Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

aa <- data.frame(res2_dmz)
aa$type <- rownames(aa)
aa <- melt(aa,id.vars=c("type"))
ggplot(aa,aes(x = variable, y = value,group = type, color=type)) +
  geom_line(lwd=1.2) +
  geom_point(aes(shape=type, color=type), size = 8) +
  theme_minimal() +
  labs(x = "vague", y = "% changement déviance résiduelle") +
  scale_x_discrete(limits=c("v1","v2","v3","v4","v5 (manquante)","v6"),
                   labels = c("v1\n(6e)","v2\n(5e)","v3\n(5e)",
                              "v4\n(4e)","v5\n(4e)\n(manquante)","v6\n(3e)")) +
  scale_color_manual(name = "Attribut",
                     values = c("#F8766D","#00BFC4","#7CAE00","#C77CFF"),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  scale_shape_manual(name = "Attribut",
                     values = c(15,16,17,18),
                     breaks = c("sex","factbrev","ethn","grade"),
                     labels = c("Genre","Or. Sociale","Or. Migratoire","Notes")) +
  ggtitle("Savoie 2 - Amis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        plot.title = element_text(hjust = 0.5,face="bold",size=32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size=24))

##### tests de robustesse - AME #####

library(margins)
ame1_rvl <- margins(m1_rvl)
ame2_rvl <- margins(m2_rvl)
ame3_rvl <- margins(m3_rvl)
ame4_rvl <- margins(m4_rvl)
ame5_rvl <- margins(m5_rvl)

ame1_jbs <- margins(m1_jbs)
ame2_jbs <- margins(m2_jbs)
ame3_jbs <- margins(m3_jbs)
ame4_jbs <- margins(m4_jbs)
ame5_jbs <- margins(m5_jbs)

ame1_pm <- margins(m1_pm)
ame2_pm <- margins(m2_pm)
ame3_pm <- margins(m3_pm)
ame4_pm <- margins(m4_pm)
ame5_pm <- margins(m5_pm)

ame1_dmz <- margins(m1_dmz)
ame2_dmz <- margins(m2_dmz)
ame3_dmz <- margins(m3_dmz)
ame4_dmz <- margins(m4_dmz)
ame5_dmz <- margins(m5_dmz)

###

# On refait les rapports des coefficients d homophilie, pour voir que ca donne pareil
# qu avec les coefs log odds.

summary(ame1_rvl)[summary(ame1_rvl)$factor=="grade_absdiff","AME"] /
  summary(ame1_rvl)[summary(ame1_rvl)$factor=="factbrev_absdiff","AME"]
summary(ame1_rvl)[summary(ame1_rvl)$factor=="ethn_nodematch","AME"] /
  summary(ame1_rvl)[summary(ame1_rvl)$factor=="factbrev_absdiff","AME"]
summary(ame1_rvl)[summary(ame1_rvl)$factor=="sex_nodematch","AME"] /
  summary(ame1_rvl)[summary(ame1_rvl)$factor=="factbrev_absdiff","AME"] # rigoureusement egal (cf. tableau 4b dans l article)

summary(ame5_rvl)[summary(ame5_rvl)$factor=="grade_absdiff","AME"] /
  summary(ame5_rvl)[summary(ame5_rvl)$factor=="factbrev_absdiff","AME"]
summary(ame5_rvl)[summary(ame5_rvl)$factor=="ethn_nodematch","AME"] /
  summary(ame5_rvl)[summary(ame5_rvl)$factor=="factbrev_absdiff","AME"]
summary(ame5_rvl)[summary(ame5_rvl)$factor=="sex_nodematch","AME"] /
  summary(ame5_rvl)[summary(ame5_rvl)$factor=="factbrev_absdiff","AME"]

#

summary(ame1_jbs)[summary(ame1_jbs)$factor=="grade_absdiff","AME"] /
  summary(ame1_jbs)[summary(ame1_jbs)$factor=="factbrev_absdiff","AME"]
summary(ame1_jbs)[summary(ame1_jbs)$factor=="ethn_nodematch","AME"] /
  summary(ame1_jbs)[summary(ame1_jbs)$factor=="factbrev_absdiff","AME"]
summary(ame1_jbs)[summary(ame1_jbs)$factor=="sex_nodematch","AME"] /
  summary(ame1_jbs)[summary(ame1_jbs)$factor=="factbrev_absdiff","AME"]

# summary(ame5_jbs)[summary(ame5_jbs)$factor=="grade_absdiff","AME"] /
#   summary(ame5_jbs)[summary(ame5_jbs)$factor=="factbrev_absdiff","AME"]
summary(ame5_jbs)[summary(ame5_jbs)$factor=="ethn_nodematch","AME"] /
  summary(ame5_jbs)[summary(ame5_jbs)$factor=="factbrev_absdiff","AME"]
summary(ame5_jbs)[summary(ame5_jbs)$factor=="sex_nodematch","AME"] /
  summary(ame5_jbs)[summary(ame5_jbs)$factor=="factbrev_absdiff","AME"]

#

summary(ame1_pm)[summary(ame1_pm)$factor=="grade_absdiff","AME"] /
  summary(ame1_pm)[summary(ame1_pm)$factor=="factbrev_absdiff","AME"]
summary(ame1_pm)[summary(ame1_pm)$factor=="ethn_nodematch","AME"] /
  summary(ame1_pm)[summary(ame1_pm)$factor=="factbrev_absdiff","AME"]
summary(ame1_pm)[summary(ame1_pm)$factor=="sex_nodematch","AME"] /
  summary(ame1_pm)[summary(ame1_pm)$factor=="factbrev_absdiff","AME"]

summary(ame5_pm)[summary(ame5_pm)$factor=="grade_absdiff","AME"] /
  summary(ame5_pm)[summary(ame5_pm)$factor=="factbrev_absdiff","AME"]
summary(ame5_pm)[summary(ame5_pm)$factor=="ethn_nodematch","AME"] /
  summary(ame5_pm)[summary(ame5_pm)$factor=="factbrev_absdiff","AME"]
summary(ame5_pm)[summary(ame5_pm)$factor=="sex_nodematch","AME"] /
  summary(ame5_pm)[summary(ame5_pm)$factor=="factbrev_absdiff","AME"]

#

summary(ame1_dmz)[summary(ame1_dmz)$factor=="grade_absdiff","AME"] /
  summary(ame1_dmz)[summary(ame1_dmz)$factor=="factbrev_absdiff","AME"]
summary(ame1_dmz)[summary(ame1_dmz)$factor=="ethn_nodematch","AME"] /
  summary(ame1_dmz)[summary(ame1_dmz)$factor=="factbrev_absdiff","AME"]
summary(ame1_dmz)[summary(ame1_dmz)$factor=="sex_nodematch","AME"] /
  summary(ame1_dmz)[summary(ame1_dmz)$factor=="factbrev_absdiff","AME"]

summary(ame5_dmz)[summary(ame5_dmz)$factor=="grade_absdiff","AME"] /
  summary(ame5_dmz)[summary(ame5_dmz)$factor=="factbrev_absdiff","AME"]
summary(ame5_dmz)[summary(ame5_dmz)$factor=="ethn_nodematch","AME"] /
  summary(ame5_dmz)[summary(ame5_dmz)$factor=="factbrev_absdiff","AME"]
summary(ame5_dmz)[summary(ame5_dmz)$factor=="sex_nodematch","AME"] /
  summary(ame5_dmz)[summary(ame5_dmz)$factor=="factbrev_absdiff","AME"]

############################################################
####################### 5 - Annexes ######################## 
############################################################
########## Annexe 1 ########## 

## Valeurs manquantes:
jonction$school <- sapply(strsplit(jonction$id,""),"[",2)

table(jonction$school,jonction$status_w1)
1 - table(jonction$school,jonction$status_w1)[,"missing"]/
  (table(jonction$school,jonction$status_w1)[,"missing"]+table(jonction$school,jonction$status_w1)[,"here"])

table(jonction$school,jonction$status_w2)
1- table(jonction$school,jonction$status_w2)[,"missing"]/
  (table(jonction$school,jonction$status_w2)[,"missing"]+table(jonction$school,jonction$status_w2)[,"here"])

table(jonction$school,jonction$status_w3)
1- table(jonction$school,jonction$status_w3)[,"missing"]/
  (table(jonction$school,jonction$status_w3)[,"missing"]+table(jonction$school,jonction$status_w3)[,"here"])

table(jonction$school,jonction$status_w4)
1-table(jonction$school,jonction$status_w4)[,"missing"]/
  (table(jonction$school,jonction$status_w4)[,"missing"]+table(jonction$school,jonction$status_w4)[,"here"])

table(jonction$school,jonction$status_w6)
1-table(jonction$school,jonction$status_w6)[,"missing"]/
  (table(jonction$school,jonction$status_w6)[,"missing"]+table(jonction$school,jonction$status_w6)[,"here"])

## PCS chef:
occupation$school <- sapply(strsplit(occupation$id,""),"[",2)
aa <- table(occupation$pcs_chef_ag,occupation$school,useNA="ifany")[,c("R","J","P","D")]
aa
colSums(aa)
cprop(aa)

## PCS chef * factbrev:
mean(occupation$factbrev[occupation$pcs_chef_ag=="1"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="2"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="3"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="4"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="5"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="6"],na.rm=T)
mean(occupation$factbrev[occupation$pcs_chef_ag=="8"],na.rm=T)

## factbrev:
summary(occupation$factbrev[occupation$school=="R"])
summary(occupation$factbrev[occupation$school=="J"])
summary(occupation$factbrev[occupation$school=="P"])
summary(occupation$factbrev[occupation$school=="D"])

sd(occupation$factbrev[occupation$school=="R"],na.rm=T)
sd(occupation$factbrev[occupation$school=="J"],na.rm=T)
sd(occupation$factbrev[occupation$school=="P"],na.rm=T)
sd(occupation$factbrev[occupation$school=="D"],na.rm=T)

## % d eleves qui restent dans la cohorte d une vague a l autre

head(jonction)
id1 <- jonction$id[jonction$status_w1%in%c("here","missing")]
id2 <- jonction$id[jonction$status_w2%in%c("here","missing")]
id3 <- jonction$id[jonction$status_w3%in%c("here","missing")]
id4 <- jonction$id[jonction$status_w4%in%c("here","missing")]
id5 <- jonction$id[jonction$status_w6%in%c("here","missing")]

# % de perte par rapport a la vague precedente
length(setdiff(id1,id2))/length(id1)
length(setdiff(id1,id3))/length(id2)
length(setdiff(id3,id4))/length(id3)
length(setdiff(id4,id5))/length(id4)
length(setdiff(id1,id5))/length(id1)

# % de nouveaux dans la vague en cours
length(setdiff(id2,id1))/length(id2)
length(setdiff(id3,id2))/length(id3)
length(setdiff(id4,id3))/length(id4)
length(setdiff(id5,id4))/length(id5)
length(setdiff(id5,id1))/length(id5)

########## Annexe 3 ##########

# Gradation des tailles d effet, methode rigoureuse.

##### allf vs gf #####

# Function that draws a random sample from net_big, of the same density as net_small.

net_big <- allf_rvl[[1]]
net_small <- gf_rvl[[1]]

net_draw <- function(net_big,net_small,nsample = 1000){
  res <- list()
  for(zz in 1:nsample){
    
    nties <- table(as.sociomatrix(net_small))["1"]
    pool_of_ties <- which(as.sociomatrix(net_big)==1,arr.ind = T) # each row is one of the ties of big net.
    
    # Draw within:
    sample_of_ties <- sample(nrow(pool_of_ties),nties,replace=F)
    sample_of_ties <- pool_of_ties[sample_of_ties,]
    
    # Create the new_net object and place 1s on the ties that were randomly drawn.
    new_net <- as.sociomatrix(net_small)
    new_net[,] <- 0
    new_net[sample_of_ties] <- 1
    
    # Stock and return:
    res[[zz]] <- new_net
  }
  
  return(res)
}

### rvl

net1_rvl <- lapply(allf_rvl,off_NA,"factbrev")
net2_rvl <- lapply(gf_rvl,off_NA,"factbrev")

diff_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$diff})
sender_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$sender})
receiver_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$receiver})

a1_rvl <- net_draw(net_big = net1_rvl[[1]],net_small = net2_rvl[[1]])
res1_rvl <- lapply(a1_rvl,function(x){
  a <- ergm(x ~ edges + edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]]))
  return(coef(a)[2])
})
res1_rvl <- unlist(res1_rvl)

coef1 <- coef(ergm(net1_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]
coef2 <- coef(ergm(net2_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]

aa <- data.frame(res1_rvl)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 105), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 105), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 1 - Amis (grand r?seau) vs Tr?s bons amis (petit r?seau)")

### jbs

net1_jbs <- lapply(allf_jbs,off_NA,"factbrev")
net2_jbs <- lapply(gf_jbs,off_NA,"factbrev")

diff_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$diff})
sender_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$sender})
receiver_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$receiver})

a1_jbs <- net_draw(net_big = net1_jbs[[1]],net_small = net2_jbs[[1]])
res1_jbs <- lapply(a1_jbs,function(x){
  a <- ergm(x ~ edges + edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]]))
  return(coef(a)[2])
})
res1_jbs <- unlist(res1_jbs)

coef1 <- coef(ergm(net1_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]
coef2 <- coef(ergm(net2_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]

aa <- data.frame(res1_jbs)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 110), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 110), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 2 - Amis (grand r?seau) vs Tr?s bons amis (petit r?seau)")

### pm

net1_pm <- lapply(allf_pm,off_NA,"factbrev")
net2_pm <- lapply(gf_pm,off_NA,"factbrev")

diff_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$diff})
sender_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$sender})
receiver_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$receiver})

a1_pm <- net_draw(net_big = net1_pm[[1]],net_small = net2_pm[[1]])
res1_pm <- lapply(a1_pm,function(x){
  a <- ergm(x ~ edges + edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]]))
  return(coef(a)[2])
})
res1_pm <- unlist(res1_pm)

coef1 <- coef(ergm(net1_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]
coef2 <- coef(ergm(net2_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]

aa <- data.frame(res1_pm)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 105), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 105), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 1 - Amis (grand r?seau) vs Tr?s bons amis (petit r?seau)")

### dmz

net1_dmz <- lapply(allf_dmz,off_NA,"factbrev")
net2_dmz <- lapply(gf_dmz,off_NA,"factbrev")

diff_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$diff})
sender_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$sender})
receiver_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$receiver})

a1_dmz <- net_draw(net_big = net1_dmz[[1]],net_small = net2_dmz[[1]])
res1_dmz <- lapply(a1_dmz,function(x){
  a <- ergm(x ~ edges + edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]]))
  return(coef(a)[2])
})
res1_dmz <- unlist(res1_dmz)

coef1 <- coef(ergm(net1_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]
coef2 <- coef(ergm(net2_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]

aa <- data.frame(res1_dmz)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 140), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 140), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 2 - Amis (grand r?seau) vs Tr?s bons amis (petit r?seau)")

##### gf vs outgf #####

# Function that draws a random sample from net_big, of the same density as net_small.

net_draw <- function(net_big,net_small,nsample = 1000){
  res <- list()
  for(zz in 1:nsample){
    
    nties <- table(as.sociomatrix(net_small))["1"]
    pool_of_ties <- which(as.sociomatrix(net_big)==1,arr.ind = T) # each row is one of the ties of big net.
    
    # Draw within:
    sample_of_ties <- sample(nrow(pool_of_ties),nties,replace=F)
    sample_of_ties <- pool_of_ties[sample_of_ties,]
    
    # Create the new_net object and place 1s on the ties that were randomly drawn.
    new_net <- as.sociomatrix(net_small)
    new_net[,] <- 0
    new_net[sample_of_ties] <- 1
    
    # Stock and return:
    res[[zz]] <- new_net
  }
  
  return(res)
}

### rvl

net1_rvl <- lapply(gf_rvl,off_NA,"factbrev")
net2_rvl <- lapply(outgf_rvl,off_NA,"factbrev")

diff_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$diff})
sender_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$sender})
receiver_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$receiver})

a1_rvl <- net_draw(net_big = net1_rvl[[1]],net_small = net2_rvl[[1]])
res1_rvl <- lapply(a1_rvl,function(x){
  a <- ergm(x ~ edges + edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]]))
  return(coef(a)[2])
})
res1_rvl <- unlist(res1_rvl)

coef1 <- coef(ergm(net1_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]
coef2 <- coef(ergm(net2_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]

aa <- data.frame(res1_rvl)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 90), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 90), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 1 - Tr?s bons amis (grand r?seau) vs Hors du coll?ge (petit r?seau)")

### jbs

net1_jbs <- lapply(gf_jbs,off_NA,"factbrev")
net2_jbs <- lapply(outgf_jbs,off_NA,"factbrev")

diff_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$diff})
sender_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$sender})
receiver_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$receiver})

a1_jbs <- net_draw(net_big = net1_jbs[[1]],net_small = net2_jbs[[1]])
res1_jbs <- lapply(a1_jbs,function(x){
  a <- ergm(x ~ edges + edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]]))
  return(coef(a)[2])
})
res1_jbs <- unlist(res1_jbs)

coef1 <- coef(ergm(net1_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]
coef2 <- coef(ergm(net2_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]

hist(res1_jbs,xlim = c(min(c(coef1,coef2,res1_jbs))-0.05,max(c(coef1,coef2,res1_jbs))+0.05),
     main = "Paris 2 - tres bons amis vs vus hors du college",
     cex.axis = 2,cex.main=2,cex.lab=2)
abline(v = coef1,col="blue")
abline(v = coef2,col="red")

# Update 15/07/2022 (v4 pour publication de l article) : version ggplot du graphique,
# et en noir et blanc?
aa <- data.frame(res1_jbs)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 110), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 110), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 2 - Tr?s bons amis (grand r?seau) vs Hors du coll?ge (petit r?seau)")

### pm

net1_pm <- lapply(gf_pm,off_NA,"factbrev")
net2_pm <- lapply(outgf_pm,off_NA,"factbrev")

diff_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$diff})
sender_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$sender})
receiver_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$receiver})

a1_pm <- net_draw(net_big = net1_pm[[1]],net_small = net2_pm[[1]])
res1_pm <- lapply(a1_pm,function(x){
  a <- ergm(x ~ edges + edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]]))
  return(coef(a)[2])
})
res1_pm <- unlist(res1_pm)

coef1 <- coef(ergm(net1_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]
coef2 <- coef(ergm(net2_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]

aa <- data.frame(res1_pm)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 105), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 105), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 1 - Tr?s bons amis (grand r?seau) vs Hors du coll?ge (petit r?seau)")

### dmz

net1_dmz <- lapply(gf_dmz,off_NA,"factbrev")
net2_dmz <- lapply(outgf_dmz,off_NA,"factbrev")

diff_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$diff})
sender_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$sender})
receiver_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$receiver})

a1_dmz <- net_draw(net_big = net1_dmz[[1]],net_small = net2_dmz[[1]])
res1_dmz <- lapply(a1_dmz,function(x){
  a <- ergm(x ~ edges + edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]]))
  return(coef(a)[2])
})
res1_dmz <- unlist(res1_dmz)

coef1 <- coef(ergm(net1_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]
coef2 <- coef(ergm(net2_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]

aa <- data.frame(res1_dmz)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 100), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 100), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 2 - Tr?s bons amis (grand r?seau) vs Hors du coll?ge (petit r?seau)")

##### gf vs pargf #####

# Function that draws a random sample from net_big, of the same density as net_small.

net_draw <- function(net_big,net_small,nsample = 1000){
  res <- list()
  for(zz in 1:nsample){
    
    nties <- table(as.sociomatrix(net_small))["1"]
    pool_of_ties <- which(as.sociomatrix(net_big)==1,arr.ind = T) # each row is one of the ties of big net.
    
    # Draw within:
    sample_of_ties <- sample(nrow(pool_of_ties),nties,replace=F)
    sample_of_ties <- pool_of_ties[sample_of_ties,]
    
    # Create the new_net object and place 1s on the ties that were randomly drawn.
    new_net <- as.sociomatrix(net_small)
    new_net[,] <- 0
    new_net[sample_of_ties] <- 1
    
    # Stock and return:
    res[[zz]] <- new_net
  }
  
  return(res)
}

### rvl

net1_rvl <- lapply(gf_rvl,off_NA,"factbrev")
net2_rvl <- lapply(pargf_rvl,off_NA,"factbrev")

diff_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$diff})
sender_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$sender})
receiver_rvl <- lapply(net1_rvl,function(x){mat_att(x,"factbrev")$receiver})

a1_rvl <- net_draw(net_big = net1_rvl[[1]],net_small = net2_rvl[[1]])
res1_rvl <- lapply(a1_rvl,function(x){
  a <- ergm(x ~ edges + edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]]))
  return(coef(a)[2])
})
res1_rvl <- unlist(res1_rvl)

coef1 <- coef(ergm(net1_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]
coef2 <- coef(ergm(net2_rvl[[1]]~edges+edgecov(diff_rvl[[1]])+edgecov(sender_rvl[[1]])+edgecov(receiver_rvl[[1]])))[2]

aa <- data.frame(res1_rvl)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 120), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 120), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 1 - Tr?s bons amis (grand r?seau) vs Connaissance des parents (petit r?seau)")

### jbs

net1_jbs <- lapply(gf_jbs,off_NA,"factbrev")
net2_jbs <- lapply(pargf_jbs,off_NA,"factbrev")

diff_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$diff})
sender_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$sender})
receiver_jbs <- lapply(net1_jbs,function(x){mat_att(x,"factbrev")$receiver})

a1_jbs <- net_draw(net_big = net1_jbs[[1]],net_small = net2_jbs[[1]])
res1_jbs <- lapply(a1_jbs,function(x){
  a <- ergm(x ~ edges + edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]]))
  return(coef(a)[2])
})
res1_jbs <- unlist(res1_jbs)

coef1 <- coef(ergm(net1_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]
coef2 <- coef(ergm(net2_jbs[[1]]~edges+edgecov(diff_jbs[[1]])+edgecov(sender_jbs[[1]])+edgecov(receiver_jbs[[1]])))[2]

aa <- data.frame(res1_jbs)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 110), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 110), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Paris 2 - Tr?s bons amis (grand r?seau) vs Connaissance des parents (petit r?seau)")

### pm

net1_pm <- lapply(gf_pm,off_NA,"factbrev")
net2_pm <- lapply(pargf_pm,off_NA,"factbrev")

diff_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$diff})
sender_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$sender})
receiver_pm <- lapply(net1_pm,function(x){mat_att(x,"factbrev")$receiver})

a1_pm <- net_draw(net_big = net1_pm[[1]],net_small = net2_pm[[1]])
res1_pm <- lapply(a1_pm,function(x){
  a <- ergm(x ~ edges + edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]]))
  return(coef(a)[2])
})
res1_pm <- unlist(res1_pm)

coef1 <- coef(ergm(net1_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]
coef2 <- coef(ergm(net2_pm[[1]]~edges+edgecov(diff_pm[[1]])+edgecov(sender_pm[[1]])+edgecov(receiver_pm[[1]])))[2]

aa <- data.frame(res1_pm)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 110), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 110), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 1 - Tr?s bons amis (grand r?seau) vs Connaissance des parents (petit r?seau)")

### dmz

net1_dmz <- lapply(gf_dmz,off_NA,"factbrev")
net2_dmz <- lapply(pargf_dmz,off_NA,"factbrev")

diff_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$diff})
sender_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$sender})
receiver_dmz <- lapply(net1_dmz,function(x){mat_att(x,"factbrev")$receiver})

a1_dmz <- net_draw(net_big = net1_dmz[[1]],net_small = net2_dmz[[1]])
res1_dmz <- lapply(a1_dmz,function(x){
  a <- ergm(x ~ edges + edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]]))
  return(coef(a)[2])
})
res1_dmz <- unlist(res1_dmz)

coef1 <- coef(ergm(net1_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]
coef2 <- coef(ergm(net2_dmz[[1]]~edges+edgecov(diff_dmz[[1]])+edgecov(sender_dmz[[1]])+edgecov(receiver_dmz[[1]])))[2]

aa <- data.frame(res1_dmz)
names(aa) <- "coef"

ggplot(aa, aes(x=coef)) + geom_histogram(color="black", fill="white") +
  geom_segment(aes(x = coef1, xend = coef1, y = 0, yend = 110), linetype = "dashed",size=1.5,colour="grey50",show.legend = T) +
  geom_segment(aes(x = coef2, xend = coef2, y = 0, yend = 110), linetype = "solid",size=1.5,colour="grey50",show.legend = T) +
  # scale_discrete_manual(name='Coefficients d\'homophilie observ?s',
  #                       labels = c("\"Grand\" r?seau","\"Petit\" r?seau"),
  #                       values = c("dashed","solid"),
  #                       aesthetics = c("linetype","linetype")) +
  theme(legend.key.size = unit(3,"line"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="plain")) +
  ggtitle("Savoie 2 - Tr?s bons amis (grand r?seau) vs Connaissance des parents (petit r?seau)")

########## Annexe 4 ##########

# Look at exposure per school.
# Run 2.2 first.

round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="rvl"])),digits=2)
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="rvl"])),digits=2)

round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="rvl"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="rvl"])),digits=2)

round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="rvl"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="rvl"])),digits=2)

round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="rvl"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="rvl"])),digits=2)

##

round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="jbs"])),digits=2)
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="jbs"])),digits=2)

round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="jbs"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="jbs"])),digits=2)

round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="jbs"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="jbs"])),digits=2)

round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="jbs"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="jbs"])),digits=2)

##

round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="pm"])),digits=2)
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="pm"])),digits=2)

round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="pm"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="pm"])),digits=2)

round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="pm"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="pm"])),digits=2)

round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="pm"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="pm"])),digits=2)

##

round(prop.table(table(ind1$nb_allf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="dmz"])),digits=2)
round(prop.table(table(ind1$nb_allf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="dmz"])),digits=2)

round(prop.table(table(ind1$nb_gf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="dmz"])),digits=2)
round(prop.table(table(ind1$nb_gf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="dmz"])),digits=2)

round(prop.table(table(ind1$nb_outgf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="dmz"])),digits=2)
round(prop.table(table(ind1$nb_outgf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="dmz"])),digits=2)

round(prop.table(table(ind1$nb_pargf_sup_factbrev[ind1$factbrev_quartiles=="inf" & ind1$school=="dmz"])),digits=2)
round(prop.table(table(ind1$nb_pargf_inf_factbrev[ind1$factbrev_quartiles=="sup" & ind1$school=="dmz"])),digits=2)

########## Annexe 5 - exposure w6 ##########

# Quick reminder: what is the average degree per type of tie?
ind6$outdeg_allf <- c(degree(allf_rvl[[5]],cmode="outdegree"),degree(allf_jbs[[5]],cmode="outdegree"),
                      degree(allf_pm[[5]],cmode="outdegree"),degree(allf_dmz[[5]],cmode="outdegree"))

ind6$outdeg_gf <- c(degree(gf_rvl[[5]],cmode="outdegree"),degree(gf_jbs[[5]],cmode="outdegree"),
                    degree(gf_pm[[5]],cmode="outdegree"),degree(gf_dmz[[5]],cmode="outdegree"))

summary(ind6$outdeg_allf)
summary(ind6$outdeg_gf)

# Here I want to look at quartiles of factbrev, and ask : how many kids in the lowest quartile have at least
# a friend in the highest one, and vice versa ?

# We use sended tie for this (how many people one named, regardless of reciprocation)
# We create variables in ind data frames.

# First, create the factbrev quantiles, per school.

levels(quant.cut(rvl6$factbrev,4))
levels(quant.cut(jbs6$factbrev,4))
levels(quant.cut(pm6$factbrev,4))
levels(quant.cut(dmz6$factbrev,4))

rvl6$factbrev_quartiles <- as.character(quant.cut(rvl6$factbrev,4,labels = c("inf","low","up","sup")))
jbs6$factbrev_quartiles <- as.character(quant.cut(jbs6$factbrev,4,labels = c("inf","low","up","sup")))
pm6$factbrev_quartiles <- as.character(quant.cut(pm6$factbrev,4,labels = c("inf","low","up","sup")))
dmz6$factbrev_quartiles <- as.character(quant.cut(dmz6$factbrev,4,labels = c("inf","low","up","sup")))

# And back to ind6:
ind6$factbrev_quartiles <- c(rvl6$factbrev_quartiles,jbs6$factbrev_quartiles,pm6$factbrev_quartiles,dmz6$factbrev_quartiles)
table(ind6$factbrev_quartiles,useNA="ifany")

ind6$nb_allf_inf_factbrev <- NA
ind6$nb_gf_inf_factbrev <- NA

ind6$nb_allf_sup_factbrev <- NA
ind6$nb_gf_sup_factbrev <- NA

for(aa in 1:nrow(ind6)){
  
  id_ego <- ind6$id[aa]
  
  if(strsplit(id_ego,"")[[1]][2]=="R"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_rvl[[5]])[allf_rvl[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_rvl[[5]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_rvl[[5]])[gf_rvl[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_rvl[[5]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%allf_ego]
    gf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%gf_ego]
    
    # Count number of inf and sup (nb : NAs are not counted, thanks to how "which" works by default).
    ind6$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind6$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    
    ind6$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind6$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="J"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_jbs[[5]])[allf_jbs[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_jbs[[5]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_jbs[[5]])[gf_jbs[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_jbs[[5]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%allf_ego]
    gf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%gf_ego]
    
    # Count number of inf and sup.
    ind6$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind6$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    
    ind6$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind6$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="P"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_pm[[5]])[allf_pm[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_pm[[5]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_pm[[5]])[gf_pm[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_pm[[5]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%allf_ego]
    gf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%gf_ego]
    
    # Count number of inf and sup.
    ind6$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind6$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    
    ind6$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind6$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
  }
  
  if(strsplit(id_ego,"")[[1]][2]=="D"){
    
    # Get the ids of ego's friends in the different networks.
    a <- symmetrize(allf_dmz[[5]])[allf_dmz[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- allf_dmz[[5]]%v%"vertex.names"
    allf_ego <- names(a[a==1 & is.na(a)==F])
    
    a <- symmetrize(gf_dmz[[5]])[gf_dmz[[5]]%v%"vertex.names"==id_ego,]
    names(a) <- gf_dmz[[5]]%v%"vertex.names"
    gf_ego <- names(a[a==1 & is.na(a)==F])
    
    # Get their factbrev quantile.
    allf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%allf_ego]
    gf_ego_quartiles <- ind6$factbrev_quartiles[ind6$id%in%gf_ego]
    
    # Count number of inf and sup.
    ind6$nb_allf_inf_factbrev[aa] <- length(which(allf_ego_quartiles=="inf"))
    ind6$nb_gf_inf_factbrev[aa] <- length(which(gf_ego_quartiles=="inf"))
    
    ind6$nb_allf_sup_factbrev[aa] <- length(which(allf_ego_quartiles=="sup"))
    ind6$nb_gf_sup_factbrev[aa] <- length(which(gf_ego_quartiles=="sup"))
  }
}


round(prop.table(table(ind6$nb_allf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="rvl"])),digits=2)
round(prop.table(table(ind6$nb_allf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="rvl"])),digits=2)

round(prop.table(table(ind6$nb_gf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="rvl"])),digits=2)
round(prop.table(table(ind6$nb_gf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="rvl"])),digits=2)

##

round(prop.table(table(ind6$nb_allf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="jbs"])),digits=2)
round(prop.table(table(ind6$nb_allf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="jbs"])),digits=2)

round(prop.table(table(ind6$nb_gf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="jbs"])),digits=2)
round(prop.table(table(ind6$nb_gf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="jbs"])),digits=2)

##

round(prop.table(table(ind6$nb_allf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="pm"])),digits=2)
round(prop.table(table(ind6$nb_allf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="pm"])),digits=2)

round(prop.table(table(ind6$nb_gf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="pm"])),digits=2)
round(prop.table(table(ind6$nb_gf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="pm"])),digits=2)

##

round(prop.table(table(ind6$nb_allf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="dmz"])),digits=2)
round(prop.table(table(ind6$nb_allf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="dmz"])),digits=2)

round(prop.table(table(ind6$nb_gf_sup_factbrev[ind6$factbrev_quartiles=="inf" & ind6$school=="dmz"])),digits=2)
round(prop.table(table(ind6$nb_gf_inf_factbrev[ind6$factbrev_quartiles=="sup" & ind6$school=="dmz"])),digits=2)

