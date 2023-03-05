#### Comments #####

# Current version : October 2022

# Script meant to setup and pre-process data from Timothee Chabot's PhD study (relationship networks
# of a cohort of ~860 students followed over 3 years ; 4 schools, 5 waves of data collection)

# The point is to have a generic setup script, that prepares the data for most analyses.
# Need to be run BEFORE any script for treaments (particularly for the article "L homophilie sociale au college")

# There is a number of options that can be specified before the setup, regarding NA management and the
# types of data that are prepared (because we don t want the exact same setup in different situation, but it
# would be quite excessive to re-write the whole script for minor changes) (see below).

# NB : see details in "3. NA management" : I can choose to exclude all kids that missed a questionnaire for a time wave,
# or leave them as long as they answered at least once at some point (which will lead to some NAs in the networks at the
# times where they did not answer).
# Note that kids that were absent for a wave (i.e. not in school) are taken off regardless.
# And those that have never answered (no parental authorization) are also always taken off.

# Use option "NA_impute" for the imputed values on individual attributes.
# 3 options : "" (no imputation), "soft" (only values for which we have reasonably good predictors)
# or "hard" (everything we can, even if it s a little barbaric).
# See script "impute missing data - ind.R" for details (can be obtained from the author on demand) .

rm(list=ls())
time1 <- Sys.time()

########################################## SCRIPT OPTIONS ##########################################
### 1. NA management options                                                                       
 NA_management <- "remove_all_nas"                                                                #
# NA_management <- "leave_time_specific_nas"                                                       #
# NA_management <- "load_imputed_networks"                                                         #
### 2. NA imputation for ind attributes (NOT networks) ("soft", "hard" or "")
NA_impute <- "soft"
### 3. Out-of-school networks                                                                      #
# Should we load and prepare them ? Pick TRUE or FALSE.                                            #                                                                     #
process_out_net <- T                                                                               #
### 4. What data to keep at the end? (to lighten working environment)                              #
# (note that is is still created and processed by the script, but it is removed at the end if F)   #
# (the point is to fasten the treatments, not that script in particular which is okay)             #
keep_dat_gf <- T                                                                                   #
keep_dat_f <- T # including allf                                                                   #
keep_dat_top5 <- T                                                                                 #
keep_dat_dl <- T                                                                                   #
keep_dat_popularities <- T # all popularity and status measures (including bullying)               #
keep_dat_psycho <- T                                                                               #
keep_dat_culture <- T                                                                              #
keep_dat_ind_all <- T                                                                              #
keep_dat_ind_per_school <- T                                                                       #
keep_osm_addresses <- T                                                                            #
keep_lockdown <- T                                                                                 #
####################################### SCRIPT OPTIONS (end) #######################################
#################### SETUP ####################
########## 1. Load packages and data ##########

.libPaths()
#setwd("D:/Sauvegarde Tim V5/Documents/These/Terrain/donnees these")
setwd("~/Documents/Sauvegarde Tim V5/Documents/These/Terrain/donnees these")

#library(xergm)
library(beepr) # make a sound when done
library(statnet)
library(R2HTML)
library(questionr)
library(stargazer)
library(ggplot2)
library(ggrepel)

## Jonction, fxd, occupation, culture and psycho :

jonction <- read.csv("./data/5 - final data/special/jonction.csv",as.is=T,fileEncoding = "UTF-8")
fxd <- read.csv("./data/5 - final data/special/fxd.csv",as.is=T,fileEncoding = "UTF-8")
fxd <- fxd[,-which(colnames(fxd)=="X")]
occupation <- read.csv("./data/5 - final data/special/occupation.csv",as.is=T,fileEncoding = "UTF-8")
culture <- read.csv("./data/5 - final data/special/culture_shortv.csv",as.is=T,fileEncoding = "UTF-8")
psycho <- read.csv("./data/5 - final data/special/psycho.csv",as.is=T,fileEncoding = "UTF-8")

## Link data sets :

load(file = "./data/5 - final data/net/gf_rvl.RData")
load(file = "./data/5 - final data/net/gf_jbs.RData")
load(file = "./data/5 - final data/net/gf_pm.RData")
load(file = "./data/5 - final data/net/gf_dmz.RData")

load(file = "./data/5 - final data/net/f_rvl.RData")
load(file = "./data/5 - final data/net/f_jbs.RData")
load(file = "./data/5 - final data/net/f_pm.RData")
load(file = "./data/5 - final data/net/f_dmz.RData")

load(file = "./data/5 - final data/net/allf_rvl.RData")
load(file = "./data/5 - final data/net/allf_jbs.RData")
load(file = "./data/5 - final data/net/allf_pm.RData")
load(file = "./data/5 - final data/net/allf_dmz.RData")

load(file = "./data/5 - final data/net/dl_rvl.RData")
load(file = "./data/5 - final data/net/dl_jbs.RData")
load(file = "./data/5 - final data/net/dl_pm.RData")
load(file = "./data/5 - final data/net/dl_dmz.RData")

load(file = "./data/5 - final data/net/top5_rvl.RData")
load(file = "./data/5 - final data/net/top5_jbs.RData")
load(file = "./data/5 - final data/net/top5_pm.RData")
load(file = "./data/5 - final data/net/top5_dmz.RData")

## Individual data sets (with short names) : 

rvl1 <- read.csv(file = "./data/5 - final data/ind/ind_c1v1_rvl.csv",as.is=T,fileEncoding = "UTF-8")
jbs1 <- read.csv(file = "./data/5 - final data/ind/ind_c1v1_jbs.csv",as.is=T, fileEncoding = "UTF-8")
pm1 <- read.csv(file = "./data/5 - final data/ind/ind_c1v1_pm.csv",as.is=T, fileEncoding = "UTF-8")
dmz1 <- read.csv(file = "./data/5 - final data/ind/ind_c1v1_dmz.csv",as.is=T, fileEncoding = "UTF-8")

rvl2 <- read.csv(file = "./data/5 - final data/ind/ind_c1v2_rvl.csv",as.is=T, fileEncoding = "UTF-8")
jbs2 <- read.csv(file = "./data/5 - final data/ind/ind_c1v2_jbs.csv",as.is=T, fileEncoding = "UTF-8")
pm2 <- read.csv(file = "./data/5 - final data/ind/ind_c1v2_pm.csv",as.is=T, fileEncoding = "UTF-8")
dmz2 <- read.csv(file = "./data/5 - final data/ind/ind_c1v2_dmz.csv",as.is=T, fileEncoding = "UTF-8")

rvl3 <- read.csv(file = "./data/5 - final data/ind/ind_c1v3_rvl.csv",as.is=T, fileEncoding = "UTF-8")
jbs3 <- read.csv(file = "./data/5 - final data/ind/ind_c1v3_jbs.csv",as.is=T, fileEncoding = "UTF-8")
pm3 <- read.csv(file = "./data/5 - final data/ind/ind_c1v3_pm.csv",as.is=T, fileEncoding = "UTF-8")
dmz3 <- read.csv(file = "./data/5 - final data/ind/ind_c1v3_dmz.csv",as.is=T, fileEncoding = "UTF-8")

rvl4 <- read.csv(file = "./data/5 - final data/ind/ind_c1v4_rvl.csv",as.is=T, fileEncoding = "UTF-8")
jbs4 <- read.csv(file = "./data/5 - final data/ind/ind_c1v4_jbs.csv",as.is=T, fileEncoding = "UTF-8")
pm4 <- read.csv(file = "./data/5 - final data/ind/ind_c1v4_pm.csv",as.is=T, fileEncoding = "UTF-8")
dmz4 <- read.csv(file = "./data/5 - final data/ind/ind_c1v4_dmz.csv",as.is=T, fileEncoding = "UTF-8")

rvl6 <- read.csv(file = "./data/5 - final data/ind/ind_c1v6_rvl.csv",as.is=T, fileEncoding = "UTF-8")
jbs6 <- read.csv(file = "./data/5 - final data/ind/ind_c1v6_jbs.csv",as.is=T, fileEncoding = "UTF-8")
pm6 <- read.csv(file = "./data/5 - final data/ind/ind_c1v6_pm.csv",as.is=T, fileEncoding = "UTF-8")
dmz6 <- read.csv(file = "./data/5 - final data/ind/ind_c1v6_dmz.csv",as.is=T, fileEncoding = "UTF-8")

## Peer status and bullying

load(file = "./data/5 - final data/net/bully_rvl.RData")
load(file = "./data/5 - final data/net/bully_jbs.RData")
load(file = "./data/5 - final data/net/bully_pm.RData")
load(file = "./data/5 - final data/net/bully_dmz.RData")

load(file = "./data/5 - final data/net/cool_rvl.RData")
load(file = "./data/5 - final data/net/cool_jbs.RData")
load(file = "./data/5 - final data/net/cool_pm.RData")
load(file = "./data/5 - final data/net/cool_dmz.RData")

load(file = "./data/5 - final data/net/liked_rvl.RData")
load(file = "./data/5 - final data/net/liked_jbs.RData")
load(file = "./data/5 - final data/net/liked_pm.RData")
load(file = "./data/5 - final data/net/liked_dmz.RData")

load(file = "./data/5 - final data/net/fear_rvl.RData")
load(file = "./data/5 - final data/net/fear_jbs.RData")
load(file = "./data/5 - final data/net/fear_pm.RData")
load(file = "./data/5 - final data/net/fear_dmz.RData")

load(file = "./data/5 - final data/net/known_rvl.RData")
load(file = "./data/5 - final data/net/known_jbs.RData")
load(file = "./data/5 - final data/net/known_pm.RData")
load(file = "./data/5 - final data/net/known_dmz.RData")

## w6 only - lockdown

load(file = "./data/5 - final data/net/lock_pers_rvl.RData")
load(file = "./data/5 - final data/net/lock_pers_jbs.RData")
load(file = "./data/5 - final data/net/lock_pers_pm.RData")
load(file = "./data/5 - final data/net/lock_pers_dmz.RData")

load(file = "./data/5 - final data/net/lock_vocal_rvl.RData")
load(file = "./data/5 - final data/net/lock_vocal_jbs.RData")
load(file = "./data/5 - final data/net/lock_vocal_pm.RData")
load(file = "./data/5 - final data/net/lock_vocal_dmz.RData")

load(file = "./data/5 - final data/net/lock_write_rvl.RData")
load(file = "./data/5 - final data/net/lock_write_jbs.RData")
load(file = "./data/5 - final data/net/lock_write_pm.RData")
load(file = "./data/5 - final data/net/lock_write_dmz.RData")

load(file = "./data/5 - final data/net/lock_work_rvl.RData")
load(file = "./data/5 - final data/net/lock_work_jbs.RData")
load(file = "./data/5 - final data/net/lock_work_pm.RData")
load(file = "./data/5 - final data/net/lock_work_dmz.RData")

jonction$school <- sapply(strsplit(jonction$id,""),"[",2)
jonction[jonction$school=="P",c("fullname","id","id_w6","status_w6","out_of_study")]

########## 2. Custom Functions ##########

attribute_in_network <- function(net,data,id_data,var_names,new_names = NA){
  
  ### HOW IT WORKS ###
  
  # "net" is a network object, "data" the data frame with nodal information, "id_data" is a variable
  # in data that corresponds to the network vertex.names, and "var_names" is a character vector giving
  # the names of the variables in data that we want to pass in the network.
  
  # NB : "data" may have more rows than there are nodes in the network, BUT all network nodes need to have a unique
  # corresponding row in data (based on id_data).
  
  # Nodal information is transfered as such (i.e as character or numeric, depending on what they are in the data frame).
  
  # "new_names" is an optional character vector that is the same lenght as var_names and give the ordered names of the 
  # attributes in the network (by default, these are the same than var_names)
  
  ### ###
  
  if(is.na(new_names[1])){
    new_names <- var_names
  }
  
  # Extract the rows in the data frames corresponding to the individuals in the network.
  new_data <- data[data[,id_data]%in%as.character(net%v%"vertex.names") & is.na(data[,id_data])==F,]
  
  # Name the rows of that data frames with the id variable :
  rownames(new_data) <- new_data[,id_data]
  
  # Use this to order new_data based on the id order in the network :
  new_data <- new_data[as.character(net%v%"vertex.names"),]
  
  # And transfer attributes :
  for(uu in 1:length(var_names)){
    var <- var_names[uu]
    name <- new_names[uu]
    set.vertex.attribute(net,name,new_data[,var])
  }
  
  return(net)
}

off_NA <- function(net, attrname){
  
  ### Function to create a new network in which all the individuals with attribute "attrname" in NA have been removed ###
  
  # Select only students whose isei is NOT NA.
  valid <- which(is.na(get.vertex.attribute(net,attrname))==F)
  
  # Built a new network from those.
  new_net <- network(as.sociomatrix(net)[valid,valid])
  
  # Transfer the attributes from old to new network.
  for(att in list.vertex.attributes(net)){
    set.vertex.attribute(new_net,att,get.vertex.attribute(net,att)[valid])
  }
  
  return(new_net)
  
}

mat_att <- function(net,att){
  
  # The function automatically recognizes whether the attribute is numeric or character. The output is a list
  # with matrices of the same size than the network.
  
  ## If numeric : matrices are "diff", "sum", "sender" and "receiver", giving for each cell [i,j], respectively :
  # abs(i-j)
  # i+j
  # i
  # j
  
  ## If character, "match" is a binary matrix that says whether i and j have the same attribute or not.
  # Then there are as many "sender" and "receiver" matrices as there are categories of the attribute, saying,
  # for each one, whether i (sender matrices) or j (emitter matrices) is of that particular category.
  
  
  if(class(get.vertex.attribute(net,att))=="numeric"){
    grid <- expand.grid(get.vertex.attribute(net,att),get.vertex.attribute(net,att))
    
    grid$sum <- grid$Var1 + grid$Var2
    grid$diff <- abs(grid$Var1 - grid$Var2)
    
    mat_sum <- 
      matrix(grid$sum,nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
    
    mat_diff <- 
      matrix(grid$diff,nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
    
    mat_sender <- 
      matrix(grid$Var1,nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
    
    mat_receiver <- 
      matrix(grid$Var2,nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
    
    rownames(mat_sum) <- as.character(net%v%"vertex.names")
    colnames(mat_sum) <- as.character(net%v%"vertex.names")
    rownames(mat_diff) <- as.character(net%v%"vertex.names")
    colnames(mat_diff) <- as.character(net%v%"vertex.names")
    rownames(mat_sender) <- as.character(net%v%"vertex.names")
    colnames(mat_sender) <- as.character(net%v%"vertex.names")
    rownames(mat_receiver) <- as.character(net%v%"vertex.names")
    colnames(mat_receiver) <- as.character(net%v%"vertex.names")
    
    res <- list()
    res[["sum"]] <- mat_sum
    res[["diff"]] <- mat_diff
    res[["sender"]] <- mat_sender
    res[["receiver"]] <- mat_receiver
    
    return(res)
    
  }
  
  if(class(get.vertex.attribute(net,att))=="character"){
    
    grid <- expand.grid(get.vertex.attribute(net,att),get.vertex.attribute(net,att), stringsAsFactors = F)
    grid$match <- 0
    grid$match[grid$Var1==grid$Var2] <- 1
    
    categories <- unique(grid$Var1[is.na(grid$Var1)==F])
    
    for(cat in categories){
      grid[,paste(cat,"sender",sep="_")] <- 0
      grid[,paste(cat,"receiver",sep="_")] <- 0
      grid[grid$Var1==cat & is.na(grid$Var1)==F,paste(cat,"sender",sep="_")] <- 1
      grid[grid$Var2==cat & is.na(grid$Var2)==F,paste(cat,"receiver",sep="_")] <- 1
    }
    
    res <- list()
    aa <- matrix(grid$match,nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
    rownames(aa) <- as.character(net%v%"vertex.names")
    colnames(aa) <- as.character(net%v%"vertex.names")
    res[["match"]] <- aa
    
    
    names <- c(paste(categories,"sender",sep="_"),paste(categories,"receiver",sep="_"))
    
    for(hop in names){
      aa <- matrix(grid[,hop],nrow = length(get.vertex.attribute(net,att)), ncol = length(get.vertex.attribute(net,att)))
      rownames(aa) <- as.character(net%v%"vertex.names")
      colnames(aa) <- as.character(net%v%"vertex.names")
      res[[hop]] <- aa
      
    }
    return(res)
  }
}

get_cstat <- function(net, formula){
  
  # Parameter names from the formula :
  parameters <- names(summary(formula))
  
  # Array to stock the results :
  change_stats <- array(data = 0, dim = c(length(parameters),nrow(as.sociomatrix(net)),nrow(as.sociomatrix(net))))
  
  # Turn the network into a matrix (so we don t have to do it a bazilion times in the loop) + some
  # objects that are worth pre-setting (for saving computation time) :
  net_as_mat <- as.sociomatrix(net)
  gg <- vector(mode="numeric")
  vect_net <- summary(formula)
  vect_newnet <- vector(mode="numeric")
  
  # Matrix that tells if the tie should be toggled from o to 1 or 1 to 0 :
  mat_toggle <- net_as_mat
  mat_toggle[net_as_mat==0] <- (-1)
  
  # Toggle on/off each tie and compute the corresponding change statistics :
  print("Computing change stats...")
  for(i in 1:nrow(net_as_mat)){
    print(paste("node",i, "out of",nrow(net_as_mat)))
    
    for(j in 1:nrow(net_as_mat)){
      # Switch the tie from 0 to 1 or from 1 to 0 (doing minus 1 then squaring it):
      newnet <- net_as_mat
      gg[[1]] <- net_as_mat[i,j]
      newnet[i,j] <- (gg-1)*(gg-1)
      
      # UPDATE : if the network is undirected, then we NEED to toggle both corresponding edges at once
      # (if not, the "network" function below will attribute an edge to both of them).
      if(is.directed(net)==F){
        newnet[j,i] <- (gg-1)*(gg-1)
      }
      
      # Back to a network object :
      newnet <- network(newnet,directed = is.directed(net))
      
      # Put all the other attributes that net had back in this new network :
      for(xx in list.vertex.attributes(net)){
        aa <- get.vertex.attribute(net,xx)
        set.vertex.attribute(newnet,xx,aa)
      }
      
      # Now get the change statistics, that is, the sufficient statistics when tie is equal to 1 minus that when it is equal to 0.
      # For this, we do net - newnet, and multiply by either 1 (the edge exists, so newnet is case 0) or -1 (the edges does not exist,
      # so newnet is case 1).
      
      vect_newnet[1:length(vect_net)] <- summary(reformulate(as.character(formula[3]), response = "newnet"))
      
      # We use array attribution (the change statistics are in the [i,j] cell of each matrix along the 1st dimension of the array.
      # Each element of that 1st dimension of the array corresponds to one parameter of the ERGM).
      
      change_stats[,i,j] <- (vect_net - vect_newnet)*mat_toggle[i,j]
      
      # NB : by default, this method sets any change statistic for any tie on the diagonal to 0
      # (because the sufficient statistics that we get from the call to summaries do not take the diagonal into
      # account, so toggling the tie makes no difference whatsoever)
      
    }}
  
  final <- list()
  final[["formula"]] <- formula
  final[["parameter_names"]] <- parameters
  final[["cstat"]] <- change_stats
  
  return(final)
  
}

get_cstat_handleNA <- function(net, formula){
  
  # First, create a fictive network with no NA on any attributes.
  # This is just to be able to run "summary(formula)" the first time and get the parameter names.
  dumb_net <- net
  
  for(xx in list.vertex.attributes(net)){
    aa <- get.vertex.attribute(net,xx)
    aa[is.na(aa)] <- aa[is.na(aa)==F][1] # replace all NAs by one non NA value of the attribute
    set.vertex.attribute(dumb_net,xx,aa)
  }
  
  # Parameter names from the formula :
  dumb_formula <-  reformulate(as.character(formula[3]), response = "dumb_net")
  parameters <- names(summary(dumb_formula))
  
  # Identify the parameters that involve an attribute, and which attribute :
  att_names <- list.vertex.attributes(net)
  parameter_att <- vector(mode = "character", length = length(parameters))
  for(xx in 1:length(parameters)){
    att <- intersect(unlist(strsplit(parameters[[xx]],".",fixed = T)),att_names)
    if(length(att)==0){parameter_att[[xx]] <- "ref"} else {parameter_att[[xx]] <- att}
  }
  
  # Identify the formula part corresponding to each parameter :
  parameter_term <- vector(mode = "character", length = length(parameters))
  terms <- as.character(formula)[[3]]
  terms <- unlist(strsplit(terms,"+",fixed=T))
  for(xx in 1:length(terms)){
    split_term <- unlist(strsplit(unlist(strsplit(unlist(strsplit(terms[[xx]],"\"")),"(",fixed=T))," "))
    split_param <- strsplit(parameters,".",fixed=T)
    
    match <- unlist(lapply(split_param,function(y){
      y[1]%in%split_term & (y[2]%in%split_term | is.na(y[2]))
    }))
    
    parameter_term[parameters%in%parameters[match]] <- terms[[xx]]
  } # mais quel enfer putain...
  
  # Ca marche !
  # cbind(parameters,parameter_att, parameter_term)
  
  ###
  
  # Now, extract the attribute names :
  att_names <- list.vertex.attributes(net)
  
  # Find the ones with some missing values :
  missing_att <- vector()
  w <- 0
  for(att in att_names){
    if(NA%in%get.vertex.attribute(net,att)){
      w <- w+1
      missing_att[[w]] <- att
    }
  }
  
  # Only keep those that are part of the formula, of course :
  missing_att <- missing_att[missing_att%in%unique(parameter_att)]
  
  # Create the list where I ll stock the ouputs for each call to get_cstat.
  big_stuff <- list()
  
  # Now for each of those...
  for(gg in 1:length(missing_att)){
    
    # Create a network with off_NA.
    subnet <- off_NA(net,missing_att[gg])
    
    # Get a new formula with the relevant terms only.
    subterms <- parameter_term[parameter_att==missing_att[gg]]
    subterms <- unique(subterms)
    subformula <- paste("subnet", "~", paste(subterms,collapse="+"), collapse = "")
    subformula <- as.formula(subformula)
    
    # Tell the user about how f... long it s going to take.
    print(paste("Running get_cstat on subnetwork",gg,"out of",length(missing_att)))
    
    # Go for get_cstat on that subnetwork.
    rondoudou <- get_cstat(net = subnet, formula = subformula)
    
    # We are going to change the format output a little bit : make it a data frame
    # using the nodes' id, with one variable per parameter.
    
    grosdoudou <- data.frame(row.names = 1:(length(subnet%v%"vertex.names")*length(subnet%v%"vertex.names")))
    grosdoudou$i <- rep(subnet%v%"vertex.names", each = length(subnet%v%"vertex.names"))
    grosdoudou$j <- rep(subnet%v%"vertex.names", times = length(subnet%v%"vertex.names"))
    grosdoudou$edge_id <- paste(grosdoudou$i,grosdoudou$j,sep="_")
    
    # Get the array info into there :
    for(crapaud in 1:length(rondoudou$parameter_names)){
      varname <- rondoudou$parameter_names[crapaud]
      grosdoudou[,varname] <- as.vector(t(rondoudou$cstat[crapaud,,]))
    }
    
    # Diagonal off :
    off <- which(grosdoudou$i == grosdoudou$j)
    grosdoudou <- grosdoudou[-off,]
    
    # This data frame is stocked in the list, then removed (not to overload the RAM) :
    big_stuff[[gg]] <- grosdoudou
    rm(grosdoudou,rondoudou)
  }
  
  # Also get the "normal" network with only terms that do not involve any missing attribute.
  refterms <- parameter_term[(parameter_att%in%missing_att)==F]
  refterms <- unique(refterms)
  refformula <- paste("net", "~", paste(refterms,collapse="+"), collapse = "")
  refformula <- as.formula(refformula)
  
  # And the last call to get_cstat.
  print("One last time, kids !")
  rondoudou <- get_cstat(net = net, formula = refformula)
  
  # Again, should be a data frame.
  
  grosdoudou <- data.frame(row.names = 1:(length(net%v%"vertex.names")*length(net%v%"vertex.names")))
  grosdoudou$i <- rep(net%v%"vertex.names", each = length(net%v%"vertex.names"))
  grosdoudou$j <- rep(net%v%"vertex.names", times = length(net%v%"vertex.names"))
  grosdoudou$edge_id <- paste(grosdoudou$i,grosdoudou$j,sep="_")
  
  # Get the array info into there :
  for(crapaud in 1:length(rondoudou$parameter_names)){
    varname <- rondoudou$parameter_names[crapaud]
    grosdoudou[,varname] <- as.vector(t(rondoudou$cstat[crapaud,,]))
  }
  
  # Diagonal off :
  off <- which(grosdoudou$i == grosdoudou$j)
  grosdoudou <- grosdoudou[-off,]
  
  # This data frame is stocked in the list, then removed (not to overload the RAM) :
  big_stuff[[1+length(big_stuff)]] <- grosdoudou
  rm(grosdoudou,rondoudou)
  
  # And we return the list of data frames (the user is big enough to merge those by himself).
  return(big_stuff)
  
}

########## 3. NA management ##########

# What is the situation ? For now, the composition (and ordering) of ind data frames and of
# network objects is the exact same. For e.g., we have :

table(pm3$id==gf_pm[[3]]%v%"vertex.names")

# In those, there are only kids who where present for sure at that wave (so for e.g. kids that were on the school
# lists but that in fact arrived later in the year have been taken off those lists)(see the "creation networks" and
# "ind data bases" scripts).
# We still have students marked as "missing" :
pm2$missing

# But these are the ones that did not answer the questionnaires (either on that wave, or because they never did, especially 
# in PM).
# In other words, for this wave at least they made no nominations, but they could be nominated by others !

# Easy to see : in the nomination matrices some rows are only NAs, but never the columns.
as.sociomatrix(gf_pm[[3]])["2P10",]
as.sociomatrix(gf_pm[[3]])[,"2P10"]

# To make things clear, we have to differentiate between 3 types of missing kids altogether :
  # missing type 0 : the kid was not present at that wave. In this case, he should ALREADY be out of the data
    # (done in earlier data cleaning on scripts that directly treat raw data)
  # missing type 1 : the kid did not answer the questionnaire at one particular time.
  # missing type 2 : the kid never answered any questionnaire (mostly PM)

# For missing type 1, it seems reasonable to use missing value imputation.
# For missing type 2, unclear, but probably safer to exclude them (although we do have information about how other
# kids nominated them, but we have no covariates whatsoever so it is pretty problematic).

###

# So what to do with that ? Well, there are 2 approaches.

## PLAN A

# In one case we just want to remove these missings altogether. This gives us clean, unambiguous networks, though
# we lose quite a bit of data. This is perhaps the best option for treatments related to ERGMs, because it is not
# easy to impute missing values from a cross-sectional model.

## PLAN B

# However, in longitudinal treatments, most softwares (includind siena) can handle missing observations pretty well.
# There is more info about each kid so transwave imputation can be performed and be relatively reliable.
# Still, I think we should only impute missing of type 1 : missing type 2 is more problematic and will be excluded
# still.

# It would be a mess to have the 2 types of data coexist (with suffixes on all object names to differentiate,
# making them really long) (plus I d need to code everything that comes later on this script TWICE, one per rule
# for missing).

# So instead, the WHOLE SCRIPT can be defined to run in "plan A" or "plan B" mode (logical condition is
# at the very beginning, even before the setup).
# No matter what condition, everything else will roll - simply, if "remove_all_nas" is selected, it will run
# the 2 next sections that will remove individuals that did not answer the questionnaires from any subsquent treatment.
# Otherwise, still some individuals will be removed, but only those that NEVER answered any questionnaire (mostly PM).

# So depending on what I need to do on a particular treatment script, I can run this setup script in either mode,
# which will remove the needed individuals (all missings, or only missings type 2).

# Note that I do not implement missing values yet, meaning when I leave some NA people, the NAs just stay as such in
# the matrices.
# For now at least, implementation should be done on the treatment scripts (maybe later I will centralize here,
# I ll see).

## PLAN C (update august 2020)

# Plan C is to use imputed networks (see the scripts on network imputation). However, this is only possible for
# gf as of now, whereas the missing value imputation is also necessary on the other networks.

# In order to get coherent network sizes, the imputed version goes along with plan B, i.e. keeping time-specific NAs
# in all networks beyond gf.
# Note that in gf, only these time-specific NAs have been imputed - not the others.

##### a. remove_all_nas - Take off anyone who did not answer the questionnaires #####

# Only under "remove_all_nas".

if(NA_management=="remove_all_nas"){
  # For now, I ll be simple : exclude all individuals who did not answer the questionnaires (i.e the whole row is NA).
  # Exceptional NA edges are set to 0.
  # We also want to take them off the ind data frames, so as to keep both types of data aligned.
  
  miss_rvl1 <- which(rvl1$missing==1 & is.na(rvl1$missing)==F)
  miss_rvl2 <- which(rvl2$missing==1 & is.na(rvl2$missing)==F)
  miss_rvl3 <- which(rvl3$missing==1 & is.na(rvl3$missing)==F)
  miss_rvl4 <- which(rvl4$missing==1 & is.na(rvl4$missing)==F)
  miss_rvl6 <- which(rvl6$missing==1 & is.na(rvl6$missing)==F)
  
  if(length(miss_rvl1)>0){rvl1 <- rvl1[-miss_rvl1,]}
  if(length(miss_rvl2)>0){rvl2 <- rvl2[-miss_rvl2,]}
  if(length(miss_rvl3)>0){rvl3 <- rvl3[-miss_rvl3,]}
  if(length(miss_rvl4)>0){rvl4 <- rvl4[-miss_rvl4,]}
  if(length(miss_rvl6)>0){rvl6 <- rvl6[-miss_rvl6,]}
  
  gf_rvl[[1]] <- network(as.sociomatrix(gf_rvl[[1]])[rvl1$id,rvl1$id])
  gf_rvl[[2]] <- network(as.sociomatrix(gf_rvl[[2]])[rvl2$id,rvl2$id])
  gf_rvl[[3]] <- network(as.sociomatrix(gf_rvl[[3]])[rvl3$id,rvl3$id])
  gf_rvl[[4]] <- network(as.sociomatrix(gf_rvl[[4]])[rvl4$id,rvl4$id])
  gf_rvl[[5]] <- network(as.sociomatrix(gf_rvl[[5]])[rvl6$id,rvl6$id])
  
  f_rvl[[1]] <- network(as.sociomatrix(f_rvl[[1]])[rvl1$id,rvl1$id])
  f_rvl[[2]] <- network(as.sociomatrix(f_rvl[[2]])[rvl2$id,rvl2$id])
  f_rvl[[3]] <- network(as.sociomatrix(f_rvl[[3]])[rvl3$id,rvl3$id])
  f_rvl[[4]] <- network(as.sociomatrix(f_rvl[[4]])[rvl4$id,rvl4$id])
  f_rvl[[5]] <- network(as.sociomatrix(f_rvl[[5]])[rvl6$id,rvl6$id])
  
  allf_rvl[[1]] <- network(as.sociomatrix(allf_rvl[[1]])[rvl1$id,rvl1$id])
  allf_rvl[[2]] <- network(as.sociomatrix(allf_rvl[[2]])[rvl2$id,rvl2$id])
  allf_rvl[[3]] <- network(as.sociomatrix(allf_rvl[[3]])[rvl3$id,rvl3$id])
  allf_rvl[[4]] <- network(as.sociomatrix(allf_rvl[[4]])[rvl4$id,rvl4$id])
  allf_rvl[[5]] <- network(as.sociomatrix(allf_rvl[[5]])[rvl6$id,rvl6$id])
  
  dl_rvl[[1]] <- network(as.sociomatrix(dl_rvl[[1]])[rvl1$id,rvl1$id])
  dl_rvl[[2]] <- network(as.sociomatrix(dl_rvl[[2]])[rvl2$id,rvl2$id])
  dl_rvl[[3]] <- network(as.sociomatrix(dl_rvl[[3]])[rvl3$id,rvl3$id])
  dl_rvl[[4]] <- network(as.sociomatrix(dl_rvl[[4]])[rvl4$id,rvl4$id])
  dl_rvl[[5]] <- network(as.sociomatrix(dl_rvl[[5]])[rvl6$id,rvl6$id])
  
  top5_rvl[[1]] <- network(as.sociomatrix(top5_rvl[[1]])[rvl1$id,rvl1$id])
  top5_rvl[[2]] <- network(as.sociomatrix(top5_rvl[[2]])[rvl2$id,rvl2$id])
  top5_rvl[[3]] <- network(as.sociomatrix(top5_rvl[[3]])[rvl3$id,rvl3$id])
  top5_rvl[[4]] <- network(as.sociomatrix(top5_rvl[[4]])[rvl4$id,rvl4$id])
  top5_rvl[[5]] <- network(as.sociomatrix(top5_rvl[[5]])[rvl6$id,rvl6$id])
  
  bully_rvl[[1]] <- network(as.sociomatrix(bully_rvl[[1]])[rvl1$id,rvl1$id])
  bully_rvl[[2]] <- network(as.sociomatrix(bully_rvl[[2]])[rvl2$id,rvl2$id])
  bully_rvl[[3]] <- network(as.sociomatrix(bully_rvl[[3]])[rvl3$id,rvl3$id])
  bully_rvl[[4]] <- network(as.sociomatrix(bully_rvl[[4]])[rvl4$id,rvl4$id])
  bully_rvl[[5]] <- network(as.sociomatrix(bully_rvl[[5]])[rvl6$id,rvl6$id])
  
  cool_rvl[[1]] <- network(as.sociomatrix(cool_rvl[[1]])[rvl1$id,rvl1$id])
  cool_rvl[[2]] <- network(as.sociomatrix(cool_rvl[[2]])[rvl2$id,rvl2$id])
  cool_rvl[[3]] <- network(as.sociomatrix(cool_rvl[[3]])[rvl3$id,rvl3$id])
  cool_rvl[[4]] <- network(as.sociomatrix(cool_rvl[[4]])[rvl4$id,rvl4$id])
  cool_rvl[[5]] <- network(as.sociomatrix(cool_rvl[[5]])[rvl6$id,rvl6$id])
  
  known_rvl[[2]] <- network(as.sociomatrix(known_rvl[[2]])[rvl2$id,rvl2$id])
  known_rvl[[3]] <- network(as.sociomatrix(known_rvl[[3]])[rvl3$id,rvl3$id])
  known_rvl[[4]] <- network(as.sociomatrix(known_rvl[[4]])[rvl4$id,rvl4$id])
  known_rvl[[5]] <- network(as.sociomatrix(known_rvl[[5]])[rvl6$id,rvl6$id])
  
  fear_rvl[[2]] <- network(as.sociomatrix(fear_rvl[[2]])[rvl2$id,rvl2$id])
  fear_rvl[[3]] <- network(as.sociomatrix(fear_rvl[[3]])[rvl3$id,rvl3$id])
  fear_rvl[[4]] <- network(as.sociomatrix(fear_rvl[[4]])[rvl4$id,rvl4$id])
  fear_rvl[[5]] <- network(as.sociomatrix(fear_rvl[[5]])[rvl6$id,rvl6$id])
  
  liked_rvl[[2]] <- network(as.sociomatrix(liked_rvl[[2]])[rvl2$id,rvl2$id])
  liked_rvl[[3]] <- network(as.sociomatrix(liked_rvl[[3]])[rvl3$id,rvl3$id])
  liked_rvl[[4]] <- network(as.sociomatrix(liked_rvl[[4]])[rvl4$id,rvl4$id])
  liked_rvl[[5]] <- network(as.sociomatrix(liked_rvl[[5]])[rvl6$id,rvl6$id])
  
  ##
  
  miss_jbs1 <- which(jbs1$missing==1 & is.na(jbs1$missing)==F)
  miss_jbs2 <- which(jbs2$missing==1 & is.na(jbs2$missing)==F)
  miss_jbs3 <- which(jbs3$missing==1 & is.na(jbs3$missing)==F)
  miss_jbs4 <- which(jbs4$missing==1 & is.na(jbs4$missing)==F)
  miss_jbs6 <- which(jbs6$missing==1 & is.na(jbs6$missing)==F)
  
  if(length(miss_jbs1)>0){jbs1 <- jbs1[-miss_jbs1,]}
  if(length(miss_jbs2)>0){jbs2 <- jbs2[-miss_jbs2,]}
  if(length(miss_jbs3)>0){jbs3 <- jbs3[-miss_jbs3,]}
  if(length(miss_jbs4)>0){jbs4 <- jbs4[-miss_jbs4,]}
  if(length(miss_jbs6)>0){jbs6 <- jbs6[-miss_jbs6,]}
  
  gf_jbs[[1]] <- network(as.sociomatrix(gf_jbs[[1]])[jbs1$id,jbs1$id])
  gf_jbs[[2]] <- network(as.sociomatrix(gf_jbs[[2]])[jbs2$id,jbs2$id])
  gf_jbs[[3]] <- network(as.sociomatrix(gf_jbs[[3]])[jbs3$id,jbs3$id])
  gf_jbs[[4]] <- network(as.sociomatrix(gf_jbs[[4]])[jbs4$id,jbs4$id])
  gf_jbs[[5]] <- network(as.sociomatrix(gf_jbs[[5]])[jbs6$id,jbs6$id])
  
  f_jbs[[1]] <- network(as.sociomatrix(f_jbs[[1]])[jbs1$id,jbs1$id])
  f_jbs[[2]] <- network(as.sociomatrix(f_jbs[[2]])[jbs2$id,jbs2$id])
  f_jbs[[3]] <- network(as.sociomatrix(f_jbs[[3]])[jbs3$id,jbs3$id])
  f_jbs[[4]] <- network(as.sociomatrix(f_jbs[[4]])[jbs4$id,jbs4$id])
  f_jbs[[5]] <- network(as.sociomatrix(f_jbs[[5]])[jbs6$id,jbs6$id])
  
  allf_jbs[[1]] <- network(as.sociomatrix(allf_jbs[[1]])[jbs1$id,jbs1$id])
  allf_jbs[[2]] <- network(as.sociomatrix(allf_jbs[[2]])[jbs2$id,jbs2$id])
  allf_jbs[[3]] <- network(as.sociomatrix(allf_jbs[[3]])[jbs3$id,jbs3$id])
  allf_jbs[[4]] <- network(as.sociomatrix(allf_jbs[[4]])[jbs4$id,jbs4$id])
  allf_jbs[[5]] <- network(as.sociomatrix(allf_jbs[[5]])[jbs6$id,jbs6$id])
  
  dl_jbs[[1]] <- network(as.sociomatrix(dl_jbs[[1]])[jbs1$id,jbs1$id])
  dl_jbs[[2]] <- network(as.sociomatrix(dl_jbs[[2]])[jbs2$id,jbs2$id])
  dl_jbs[[3]] <- network(as.sociomatrix(dl_jbs[[3]])[jbs3$id,jbs3$id])
  dl_jbs[[4]] <- network(as.sociomatrix(dl_jbs[[4]])[jbs4$id,jbs4$id])
  dl_jbs[[5]] <- network(as.sociomatrix(dl_jbs[[5]])[jbs6$id,jbs6$id])
  
  top5_jbs[[1]] <- network(as.sociomatrix(top5_jbs[[1]])[jbs1$id,jbs1$id])
  top5_jbs[[2]] <- network(as.sociomatrix(top5_jbs[[2]])[jbs2$id,jbs2$id])
  top5_jbs[[3]] <- network(as.sociomatrix(top5_jbs[[3]])[jbs3$id,jbs3$id])
  top5_jbs[[4]] <- network(as.sociomatrix(top5_jbs[[4]])[jbs4$id,jbs4$id])
  top5_jbs[[5]] <- network(as.sociomatrix(top5_jbs[[5]])[jbs6$id,jbs6$id])
  
  bully_jbs[[1]] <- network(as.sociomatrix(bully_jbs[[1]])[jbs1$id,jbs1$id])
  bully_jbs[[2]] <- network(as.sociomatrix(bully_jbs[[2]])[jbs2$id,jbs2$id])
  bully_jbs[[3]] <- network(as.sociomatrix(bully_jbs[[3]])[jbs3$id,jbs3$id])
  bully_jbs[[4]] <- network(as.sociomatrix(bully_jbs[[4]])[jbs4$id,jbs4$id])
  bully_jbs[[5]] <- network(as.sociomatrix(bully_jbs[[5]])[jbs6$id,jbs6$id])
  
  cool_jbs[[1]] <- network(as.sociomatrix(cool_jbs[[1]])[jbs1$id,jbs1$id])
  cool_jbs[[2]] <- network(as.sociomatrix(cool_jbs[[2]])[jbs2$id,jbs2$id])
  cool_jbs[[3]] <- network(as.sociomatrix(cool_jbs[[3]])[jbs3$id,jbs3$id])
  cool_jbs[[4]] <- network(as.sociomatrix(cool_jbs[[4]])[jbs4$id,jbs4$id])
  cool_jbs[[5]] <- network(as.sociomatrix(cool_jbs[[5]])[jbs6$id,jbs6$id])
  
  known_jbs[[2]] <- network(as.sociomatrix(known_jbs[[2]])[jbs2$id,jbs2$id])
  known_jbs[[3]] <- network(as.sociomatrix(known_jbs[[3]])[jbs3$id,jbs3$id])
  known_jbs[[4]] <- network(as.sociomatrix(known_jbs[[4]])[jbs4$id,jbs4$id])
  known_jbs[[5]] <- network(as.sociomatrix(known_jbs[[5]])[jbs6$id,jbs6$id])
  
  fear_jbs[[2]] <- network(as.sociomatrix(fear_jbs[[2]])[jbs2$id,jbs2$id])
  fear_jbs[[3]] <- network(as.sociomatrix(fear_jbs[[3]])[jbs3$id,jbs3$id])
  fear_jbs[[4]] <- network(as.sociomatrix(fear_jbs[[4]])[jbs4$id,jbs4$id])
  fear_jbs[[5]] <- network(as.sociomatrix(fear_jbs[[5]])[jbs6$id,jbs6$id])
  
  liked_jbs[[2]] <- network(as.sociomatrix(liked_jbs[[2]])[jbs2$id,jbs2$id])
  liked_jbs[[3]] <- network(as.sociomatrix(liked_jbs[[3]])[jbs3$id,jbs3$id])
  liked_jbs[[4]] <- network(as.sociomatrix(liked_jbs[[4]])[jbs4$id,jbs4$id])
  liked_jbs[[5]] <- network(as.sociomatrix(liked_jbs[[5]])[jbs6$id,jbs6$id])
  
  ##
  
  miss_pm1 <- which(pm1$missing==1 & is.na(pm1$missing)==F)
  miss_pm2 <- which(pm2$missing==1 & is.na(pm2$missing)==F)
  miss_pm3 <- which(pm3$missing==1 & is.na(pm3$missing)==F)
  miss_pm4 <- which(pm4$missing==1 & is.na(pm4$missing)==F)
  miss_pm6 <- which(pm6$missing==1 & is.na(pm6$missing)==F)
  
  if(length(miss_pm1)>0){pm1 <- pm1[-miss_pm1,]}
  if(length(miss_pm2)>0){pm2 <- pm2[-miss_pm2,]}
  if(length(miss_pm3)>0){pm3 <- pm3[-miss_pm3,]}
  if(length(miss_pm4)>0){pm4 <- pm4[-miss_pm4,]}
  if(length(miss_pm6)>0){pm6 <- pm6[-miss_pm6,]}
  
  gf_pm[[1]] <- network(as.sociomatrix(gf_pm[[1]])[pm1$id,pm1$id])
  gf_pm[[2]] <- network(as.sociomatrix(gf_pm[[2]])[pm2$id,pm2$id])
  gf_pm[[3]] <- network(as.sociomatrix(gf_pm[[3]])[pm3$id,pm3$id])
  gf_pm[[4]] <- network(as.sociomatrix(gf_pm[[4]])[pm4$id,pm4$id])
  gf_pm[[5]] <- network(as.sociomatrix(gf_pm[[5]])[pm6$id,pm6$id])
  
  f_pm[[1]] <- network(as.sociomatrix(f_pm[[1]])[pm1$id,pm1$id])
  f_pm[[2]] <- network(as.sociomatrix(f_pm[[2]])[pm2$id,pm2$id])
  f_pm[[3]] <- network(as.sociomatrix(f_pm[[3]])[pm3$id,pm3$id])
  f_pm[[4]] <- network(as.sociomatrix(f_pm[[4]])[pm4$id,pm4$id])
  f_pm[[5]] <- network(as.sociomatrix(f_pm[[5]])[pm6$id,pm6$id])
  
  allf_pm[[1]] <- network(as.sociomatrix(allf_pm[[1]])[pm1$id,pm1$id])
  allf_pm[[2]] <- network(as.sociomatrix(allf_pm[[2]])[pm2$id,pm2$id])
  allf_pm[[3]] <- network(as.sociomatrix(allf_pm[[3]])[pm3$id,pm3$id])
  allf_pm[[4]] <- network(as.sociomatrix(allf_pm[[4]])[pm4$id,pm4$id])
  allf_pm[[5]] <- network(as.sociomatrix(allf_pm[[5]])[pm6$id,pm6$id])
  
  dl_pm[[1]] <- network(as.sociomatrix(dl_pm[[1]])[pm1$id,pm1$id])
  dl_pm[[2]] <- network(as.sociomatrix(dl_pm[[2]])[pm2$id,pm2$id])
  dl_pm[[3]] <- network(as.sociomatrix(dl_pm[[3]])[pm3$id,pm3$id])
  dl_pm[[4]] <- network(as.sociomatrix(dl_pm[[4]])[pm4$id,pm4$id])
  dl_pm[[5]] <- network(as.sociomatrix(dl_pm[[5]])[pm6$id,pm6$id])
  
  top5_pm[[1]] <- network(as.sociomatrix(top5_pm[[1]])[pm1$id,pm1$id])
  top5_pm[[2]] <- network(as.sociomatrix(top5_pm[[2]])[pm2$id,pm2$id])
  top5_pm[[3]] <- network(as.sociomatrix(top5_pm[[3]])[pm3$id,pm3$id])
  top5_pm[[4]] <- network(as.sociomatrix(top5_pm[[4]])[pm4$id,pm4$id])
  top5_pm[[5]] <- network(as.sociomatrix(top5_pm[[5]])[pm6$id,pm6$id])
  
  bully_pm[[1]] <- network(as.sociomatrix(bully_pm[[1]])[pm1$id,pm1$id])
  bully_pm[[2]] <- network(as.sociomatrix(bully_pm[[2]])[pm2$id,pm2$id])
  bully_pm[[3]] <- network(as.sociomatrix(bully_pm[[3]])[pm3$id,pm3$id])
  bully_pm[[4]] <- network(as.sociomatrix(bully_pm[[4]])[pm4$id,pm4$id])
  bully_pm[[5]] <- network(as.sociomatrix(bully_pm[[5]])[pm6$id,pm6$id])
  
  cool_pm[[1]] <- network(as.sociomatrix(cool_pm[[1]])[pm1$id,pm1$id])
  cool_pm[[2]] <- network(as.sociomatrix(cool_pm[[2]])[pm2$id,pm2$id])
  cool_pm[[3]] <- network(as.sociomatrix(cool_pm[[3]])[pm3$id,pm3$id])
  cool_pm[[4]] <- network(as.sociomatrix(cool_pm[[4]])[pm4$id,pm4$id])
  cool_pm[[5]] <- network(as.sociomatrix(cool_pm[[5]])[pm6$id,pm6$id])
  
  known_pm[[2]] <- network(as.sociomatrix(known_pm[[2]])[pm2$id,pm2$id])
  known_pm[[3]] <- network(as.sociomatrix(known_pm[[3]])[pm3$id,pm3$id])
  known_pm[[4]] <- network(as.sociomatrix(known_pm[[4]])[pm4$id,pm4$id])
  known_pm[[5]] <- network(as.sociomatrix(known_pm[[5]])[pm6$id,pm6$id])
  
  fear_pm[[2]] <- network(as.sociomatrix(fear_pm[[2]])[pm2$id,pm2$id])
  fear_pm[[3]] <- network(as.sociomatrix(fear_pm[[3]])[pm3$id,pm3$id])
  fear_pm[[4]] <- network(as.sociomatrix(fear_pm[[4]])[pm4$id,pm4$id])
  fear_pm[[5]] <- network(as.sociomatrix(fear_pm[[5]])[pm6$id,pm6$id])
  
  liked_pm[[2]] <- network(as.sociomatrix(liked_pm[[2]])[pm2$id,pm2$id])
  liked_pm[[3]] <- network(as.sociomatrix(liked_pm[[3]])[pm3$id,pm3$id])
  liked_pm[[4]] <- network(as.sociomatrix(liked_pm[[4]])[pm4$id,pm4$id])
  liked_pm[[5]] <- network(as.sociomatrix(liked_pm[[5]])[pm6$id,pm6$id])
  
  ##
  
  miss_dmz1 <- which(dmz1$missing==1 & is.na(dmz1$missing)==F)
  miss_dmz2 <- which(dmz2$missing==1 & is.na(dmz2$missing)==F)
  miss_dmz3 <- which(dmz3$missing==1 & is.na(dmz3$missing)==F)
  miss_dmz4 <- which(dmz4$missing==1 & is.na(dmz4$missing)==F)
  miss_dmz6 <- which(dmz6$missing==1 & is.na(dmz6$missing)==F)
  
  if(length(miss_dmz1)>0){dmz1 <- dmz1[-miss_dmz1,]}
  if(length(miss_dmz2)>0){dmz2 <- dmz2[-miss_dmz2,]}
  if(length(miss_dmz3)>0){dmz3 <- dmz3[-miss_dmz3,]}
  if(length(miss_dmz4)>0){dmz4 <- dmz4[-miss_dmz4,]}
  if(length(miss_dmz6)>0){dmz6 <- dmz6[-miss_dmz6,]}
  
  gf_dmz[[1]] <- network(as.sociomatrix(gf_dmz[[1]])[dmz1$id,dmz1$id])
  gf_dmz[[2]] <- network(as.sociomatrix(gf_dmz[[2]])[dmz2$id,dmz2$id])
  gf_dmz[[3]] <- network(as.sociomatrix(gf_dmz[[3]])[dmz3$id,dmz3$id])
  gf_dmz[[4]] <- network(as.sociomatrix(gf_dmz[[4]])[dmz4$id,dmz4$id])
  gf_dmz[[5]] <- network(as.sociomatrix(gf_dmz[[5]])[dmz6$id,dmz6$id])
  
  f_dmz[[1]] <- network(as.sociomatrix(f_dmz[[1]])[dmz1$id,dmz1$id])
  f_dmz[[2]] <- network(as.sociomatrix(f_dmz[[2]])[dmz2$id,dmz2$id])
  f_dmz[[3]] <- network(as.sociomatrix(f_dmz[[3]])[dmz3$id,dmz3$id])
  f_dmz[[4]] <- network(as.sociomatrix(f_dmz[[4]])[dmz4$id,dmz4$id])
  f_dmz[[5]] <- network(as.sociomatrix(f_dmz[[5]])[dmz6$id,dmz6$id])
  
  allf_dmz[[1]] <- network(as.sociomatrix(allf_dmz[[1]])[dmz1$id,dmz1$id])
  allf_dmz[[2]] <- network(as.sociomatrix(allf_dmz[[2]])[dmz2$id,dmz2$id])
  allf_dmz[[3]] <- network(as.sociomatrix(allf_dmz[[3]])[dmz3$id,dmz3$id])
  allf_dmz[[4]] <- network(as.sociomatrix(allf_dmz[[4]])[dmz4$id,dmz4$id])
  allf_dmz[[5]] <- network(as.sociomatrix(allf_dmz[[5]])[dmz6$id,dmz6$id])
  
  dl_dmz[[1]] <- network(as.sociomatrix(dl_dmz[[1]])[dmz1$id,dmz1$id])
  dl_dmz[[2]] <- network(as.sociomatrix(dl_dmz[[2]])[dmz2$id,dmz2$id])
  dl_dmz[[3]] <- network(as.sociomatrix(dl_dmz[[3]])[dmz3$id,dmz3$id])
  dl_dmz[[4]] <- network(as.sociomatrix(dl_dmz[[4]])[dmz4$id,dmz4$id])
  dl_dmz[[5]] <- network(as.sociomatrix(dl_dmz[[5]])[dmz6$id,dmz6$id])
  
  top5_dmz[[1]] <- network(as.sociomatrix(top5_dmz[[1]])[dmz1$id,dmz1$id])
  top5_dmz[[2]] <- network(as.sociomatrix(top5_dmz[[2]])[dmz2$id,dmz2$id])
  top5_dmz[[3]] <- network(as.sociomatrix(top5_dmz[[3]])[dmz3$id,dmz3$id])
  top5_dmz[[4]] <- network(as.sociomatrix(top5_dmz[[4]])[dmz4$id,dmz4$id])
  top5_dmz[[5]] <- network(as.sociomatrix(top5_dmz[[5]])[dmz6$id,dmz6$id])
  
  bully_dmz[[1]] <- network(as.sociomatrix(bully_dmz[[1]])[dmz1$id,dmz1$id])
  bully_dmz[[2]] <- network(as.sociomatrix(bully_dmz[[2]])[dmz2$id,dmz2$id])
  bully_dmz[[3]] <- network(as.sociomatrix(bully_dmz[[3]])[dmz3$id,dmz3$id])
  bully_dmz[[4]] <- network(as.sociomatrix(bully_dmz[[4]])[dmz4$id,dmz4$id])
  bully_dmz[[5]] <- network(as.sociomatrix(bully_dmz[[5]])[dmz6$id,dmz6$id])
  
  cool_dmz[[1]] <- network(as.sociomatrix(cool_dmz[[1]])[dmz1$id,dmz1$id])
  cool_dmz[[2]] <- network(as.sociomatrix(cool_dmz[[2]])[dmz2$id,dmz2$id])
  cool_dmz[[3]] <- network(as.sociomatrix(cool_dmz[[3]])[dmz3$id,dmz3$id])
  cool_dmz[[4]] <- network(as.sociomatrix(cool_dmz[[4]])[dmz4$id,dmz4$id])
  cool_dmz[[5]] <- network(as.sociomatrix(cool_dmz[[5]])[dmz6$id,dmz6$id])
  
  known_dmz[[2]] <- network(as.sociomatrix(known_dmz[[2]])[dmz2$id,dmz2$id])
  known_dmz[[3]] <- network(as.sociomatrix(known_dmz[[3]])[dmz3$id,dmz3$id])
  known_dmz[[4]] <- network(as.sociomatrix(known_dmz[[4]])[dmz4$id,dmz4$id])
  known_dmz[[5]] <- network(as.sociomatrix(known_dmz[[5]])[dmz6$id,dmz6$id])
  
  fear_dmz[[2]] <- network(as.sociomatrix(fear_dmz[[2]])[dmz2$id,dmz2$id])
  fear_dmz[[3]] <- network(as.sociomatrix(fear_dmz[[3]])[dmz3$id,dmz3$id])
  fear_dmz[[4]] <- network(as.sociomatrix(fear_dmz[[4]])[dmz4$id,dmz4$id])
  fear_dmz[[5]] <- network(as.sociomatrix(fear_dmz[[5]])[dmz6$id,dmz6$id])
  
  liked_dmz[[2]] <- network(as.sociomatrix(liked_dmz[[2]])[dmz2$id,dmz2$id])
  liked_dmz[[3]] <- network(as.sociomatrix(liked_dmz[[3]])[dmz3$id,dmz3$id])
  liked_dmz[[4]] <- network(as.sociomatrix(liked_dmz[[4]])[dmz4$id,dmz4$id])
  liked_dmz[[5]] <- network(as.sociomatrix(liked_dmz[[5]])[dmz6$id,dmz6$id])
  
  rm(miss_rvl1,miss_rvl2,miss_rvl3,
     miss_jbs1,miss_jbs2,miss_jbs3,
     miss_pm1,miss_pm2,miss_pm3,
     miss_dmz1,miss_dmz2,miss_dmz3) 
}

##### b. remove_all_nas - Set the remaining NAs to 0 #####

# Only under "remove_all_nas".

if(NA_management=="remove_all_nas"){
  # NB : it seems as if now, ergm does handle NA for edges.
  # So I ll just leave it like that for now ?? (see in Zotero, methodo, missing network data, the statnet tutorial on this)
  
  # UPDATE : ergm handles it, but not bterm because terms that takes previous networks as independant variables (e.g.
  # autoregression) get fucked by the NAs.
  
  # So NAs are 0.
  
  lapply(gf_rvl,function(x){length(which(is.na(as.sociomatrix(x))))})
  lapply(gf_jbs,function(x){length(which(is.na(as.sociomatrix(x))))})
  lapply(gf_pm,function(x){length(which(is.na(as.sociomatrix(x))))})
  lapply(gf_dmz,function(x){length(which(is.na(as.sociomatrix(x))))}) # A few cells, but not so much, so it s okay.
  
  gf_rvl <- lapply(gf_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  gf_jbs <- lapply(gf_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  gf_pm <- lapply(gf_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  gf_dmz <- lapply(gf_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  f_rvl <- lapply(f_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  f_jbs <- lapply(f_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  f_pm <- lapply(f_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  f_dmz <- lapply(f_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  allf_rvl <- lapply(allf_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  allf_jbs <- lapply(allf_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  allf_pm <- lapply(allf_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  allf_dmz <- lapply(allf_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  dl_rvl <- lapply(dl_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  dl_jbs <- lapply(dl_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  dl_pm <- lapply(dl_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  dl_dmz <- lapply(dl_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  top5_rvl <- lapply(top5_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  top5_jbs <- lapply(top5_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  top5_pm <- lapply(top5_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  top5_dmz <- lapply(top5_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  bully_rvl <- lapply(bully_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  bully_jbs <- lapply(bully_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  bully_pm <- lapply(bully_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  bully_dmz <- lapply(bully_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  cool_rvl <- lapply(cool_rvl,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  cool_jbs <- lapply(cool_jbs,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  cool_pm <- lapply(cool_pm,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  cool_dmz <- lapply(cool_dmz,function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  known_rvl[2:length(known_rvl)] <- lapply(known_rvl[2:length(known_rvl)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  known_jbs[2:length(known_jbs)] <- lapply(known_jbs[2:length(known_jbs)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  known_pm[2:length(known_pm)] <- lapply(known_pm[2:length(known_pm)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  known_dmz[2:length(known_dmz)] <- lapply(known_dmz[2:length(known_dmz)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  fear_rvl[2:length(fear_rvl)] <- lapply(fear_rvl[2:length(fear_rvl)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  fear_jbs[2:length(fear_jbs)] <- lapply(fear_jbs[2:length(fear_jbs)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  fear_pm[2:length(fear_pm)] <- lapply(fear_pm[2:length(fear_pm)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  fear_dmz[2:length(fear_dmz)] <- lapply(fear_dmz[2:length(fear_dmz)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  ##
  
  liked_rvl[2:length(liked_rvl)] <- lapply(liked_rvl[2:length(liked_rvl)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  liked_jbs[2:length(liked_jbs)] <- lapply(liked_jbs[2:length(liked_jbs)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  liked_pm[2:length(liked_pm)] <- lapply(liked_pm[2:length(liked_pm)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
  
  liked_dmz[2:length(liked_dmz)] <- lapply(liked_dmz[2:length(liked_dmz)],function(x){
    a <- which(is.na(as.sociomatrix(x)),arr.ind=T)
    if(nrow(a)>0){
      for(i in 1:nrow(a)){
        x[a[i,1],a[i,2]] <- 0
      }
    }
    return(x)
  })
}

##### c. leave_time_specific_nas - take off anyone who NEVER answered ANY questionnaire #####

# It s basically like a, but with a more narrow definition of the "missing" vectors.

if(NA_management=="leave_time_specific_nas" | NA_management=="load_imputed_networks"){

  miss_rvl1 <- rvl1$id[rvl1$missing==1 & is.na(rvl1$missing)==F]
  miss_rvl2 <- rvl2$id[rvl2$missing==1 & is.na(rvl2$missing)==F]
  miss_rvl3 <- rvl3$id[rvl3$missing==1 & is.na(rvl3$missing)==F]
  miss_rvl4 <- rvl4$id[rvl4$missing==1 & is.na(rvl4$missing)==F]
  miss_rvl6 <- rvl6$id[rvl6$missing==1 & is.na(rvl6$missing)==F]
  
  here_rvl1 <- rvl1$id[rvl1$missing==0 | is.na(rvl1$missing)==T]
  here_rvl2 <- rvl2$id[rvl2$missing==0 | is.na(rvl2$missing)==T]
  here_rvl3 <- rvl3$id[rvl3$missing==0 | is.na(rvl3$missing)==T]
  here_rvl4 <- rvl4$id[rvl4$missing==0 | is.na(rvl4$missing)==T]
  here_rvl6 <- rvl6$id[rvl6$missing==0 | is.na(rvl6$missing)==T]
  
  here_any <- union(here_rvl1,union(here_rvl2,union(here_rvl3,union(here_rvl4,here_rvl6))))
  here_any
  
  # The only "true" missings are those that were never here, ie never part of "here_any".
  miss_rvl1 <- miss_rvl1[(miss_rvl1%in%here_any)==F]
  miss_rvl2 <- miss_rvl2[(miss_rvl2%in%here_any)==F]
  miss_rvl3 <- miss_rvl3[(miss_rvl3%in%here_any)==F]
  miss_rvl4 <- miss_rvl4[(miss_rvl4%in%here_any)==F]
  miss_rvl6 <- miss_rvl6[(miss_rvl6%in%here_any)==F]
  
  miss_rvl1
  miss_rvl2
  miss_rvl3
  miss_rvl4
  miss_rvl6 # not much left (if any)
  
  # Now we dont want the id but the row number in ind data :
  miss_rvl1 <- which(rvl1$id%in%miss_rvl1)
  miss_rvl2 <- which(rvl2$id%in%miss_rvl2)
  miss_rvl3 <- which(rvl3$id%in%miss_rvl3)
  miss_rvl4 <- which(rvl4$id%in%miss_rvl4)
  miss_rvl6 <- which(rvl6$id%in%miss_rvl6)
    
  if(length(miss_rvl1)>0){rvl1 <- rvl1[-miss_rvl1,]}
  if(length(miss_rvl2)>0){rvl2 <- rvl2[-miss_rvl2,]}
  if(length(miss_rvl3)>0){rvl3 <- rvl3[-miss_rvl3,]}
  if(length(miss_rvl4)>0){rvl4 <- rvl4[-miss_rvl4,]}
  if(length(miss_rvl6)>0){rvl6 <- rvl6[-miss_rvl6,]}
  
  gf_rvl[[1]] <- network(as.sociomatrix(gf_rvl[[1]])[rvl1$id,rvl1$id])
  gf_rvl[[2]] <- network(as.sociomatrix(gf_rvl[[2]])[rvl2$id,rvl2$id])
  gf_rvl[[3]] <- network(as.sociomatrix(gf_rvl[[3]])[rvl3$id,rvl3$id])
  gf_rvl[[4]] <- network(as.sociomatrix(gf_rvl[[4]])[rvl4$id,rvl4$id])
  gf_rvl[[5]] <- network(as.sociomatrix(gf_rvl[[5]])[rvl6$id,rvl6$id])
  
  f_rvl[[1]] <- network(as.sociomatrix(f_rvl[[1]])[rvl1$id,rvl1$id])
  f_rvl[[2]] <- network(as.sociomatrix(f_rvl[[2]])[rvl2$id,rvl2$id])
  f_rvl[[3]] <- network(as.sociomatrix(f_rvl[[3]])[rvl3$id,rvl3$id])
  f_rvl[[4]] <- network(as.sociomatrix(f_rvl[[4]])[rvl4$id,rvl4$id])
  f_rvl[[5]] <- network(as.sociomatrix(f_rvl[[5]])[rvl6$id,rvl6$id])
  
  allf_rvl[[1]] <- network(as.sociomatrix(allf_rvl[[1]])[rvl1$id,rvl1$id])
  allf_rvl[[2]] <- network(as.sociomatrix(allf_rvl[[2]])[rvl2$id,rvl2$id])
  allf_rvl[[3]] <- network(as.sociomatrix(allf_rvl[[3]])[rvl3$id,rvl3$id])
  allf_rvl[[4]] <- network(as.sociomatrix(allf_rvl[[4]])[rvl4$id,rvl4$id])
  allf_rvl[[5]] <- network(as.sociomatrix(allf_rvl[[5]])[rvl6$id,rvl6$id])
  
  dl_rvl[[1]] <- network(as.sociomatrix(dl_rvl[[1]])[rvl1$id,rvl1$id])
  dl_rvl[[2]] <- network(as.sociomatrix(dl_rvl[[2]])[rvl2$id,rvl2$id])
  dl_rvl[[3]] <- network(as.sociomatrix(dl_rvl[[3]])[rvl3$id,rvl3$id])
  dl_rvl[[4]] <- network(as.sociomatrix(dl_rvl[[4]])[rvl4$id,rvl4$id])
  dl_rvl[[5]] <- network(as.sociomatrix(dl_rvl[[5]])[rvl6$id,rvl6$id])
  
  top5_rvl[[1]] <- network(as.sociomatrix(top5_rvl[[1]])[rvl1$id,rvl1$id])
  top5_rvl[[2]] <- network(as.sociomatrix(top5_rvl[[2]])[rvl2$id,rvl2$id])
  top5_rvl[[3]] <- network(as.sociomatrix(top5_rvl[[3]])[rvl3$id,rvl3$id])
  top5_rvl[[4]] <- network(as.sociomatrix(top5_rvl[[4]])[rvl4$id,rvl4$id])
  top5_rvl[[5]] <- network(as.sociomatrix(top5_rvl[[5]])[rvl6$id,rvl6$id])
  
  bully_rvl[[1]] <- network(as.sociomatrix(bully_rvl[[1]])[rvl1$id,rvl1$id])
  bully_rvl[[2]] <- network(as.sociomatrix(bully_rvl[[2]])[rvl2$id,rvl2$id])
  bully_rvl[[3]] <- network(as.sociomatrix(bully_rvl[[3]])[rvl3$id,rvl3$id])
  bully_rvl[[4]] <- network(as.sociomatrix(bully_rvl[[4]])[rvl4$id,rvl4$id])
  bully_rvl[[5]] <- network(as.sociomatrix(bully_rvl[[5]])[rvl6$id,rvl6$id])
  
  cool_rvl[[1]] <- network(as.sociomatrix(cool_rvl[[1]])[rvl1$id,rvl1$id])
  cool_rvl[[2]] <- network(as.sociomatrix(cool_rvl[[2]])[rvl2$id,rvl2$id])
  cool_rvl[[3]] <- network(as.sociomatrix(cool_rvl[[3]])[rvl3$id,rvl3$id])
  cool_rvl[[4]] <- network(as.sociomatrix(cool_rvl[[4]])[rvl4$id,rvl4$id])
  cool_rvl[[5]] <- network(as.sociomatrix(cool_rvl[[5]])[rvl6$id,rvl6$id])
  
  known_rvl[[2]] <- network(as.sociomatrix(known_rvl[[2]])[rvl2$id,rvl2$id])
  known_rvl[[3]] <- network(as.sociomatrix(known_rvl[[3]])[rvl3$id,rvl3$id])
  known_rvl[[4]] <- network(as.sociomatrix(known_rvl[[4]])[rvl4$id,rvl4$id])
  known_rvl[[5]] <- network(as.sociomatrix(known_rvl[[5]])[rvl6$id,rvl6$id])
  
  fear_rvl[[2]] <- network(as.sociomatrix(fear_rvl[[2]])[rvl2$id,rvl2$id])
  fear_rvl[[3]] <- network(as.sociomatrix(fear_rvl[[3]])[rvl3$id,rvl3$id])
  fear_rvl[[4]] <- network(as.sociomatrix(fear_rvl[[4]])[rvl4$id,rvl4$id])
  fear_rvl[[5]] <- network(as.sociomatrix(fear_rvl[[5]])[rvl6$id,rvl6$id])
  
  liked_rvl[[2]] <- network(as.sociomatrix(liked_rvl[[2]])[rvl2$id,rvl2$id])
  liked_rvl[[3]] <- network(as.sociomatrix(liked_rvl[[3]])[rvl3$id,rvl3$id])
  liked_rvl[[4]] <- network(as.sociomatrix(liked_rvl[[4]])[rvl4$id,rvl4$id])
  liked_rvl[[5]] <- network(as.sociomatrix(liked_rvl[[5]])[rvl6$id,rvl6$id])
  
  ##
  
  miss_jbs1 <- jbs1$id[jbs1$missing==1 & is.na(jbs1$missing)==F]
  miss_jbs2 <- jbs2$id[jbs2$missing==1 & is.na(jbs2$missing)==F]
  miss_jbs3 <- jbs3$id[jbs3$missing==1 & is.na(jbs3$missing)==F]
  miss_jbs4 <- jbs4$id[jbs4$missing==1 & is.na(jbs4$missing)==F]
  miss_jbs6 <- jbs6$id[jbs6$missing==1 & is.na(jbs6$missing)==F]
  
  here_jbs1 <- jbs1$id[jbs1$missing==0 | is.na(jbs1$missing)==T]
  here_jbs2 <- jbs2$id[jbs2$missing==0 | is.na(jbs2$missing)==T]
  here_jbs3 <- jbs3$id[jbs3$missing==0 | is.na(jbs3$missing)==T]
  here_jbs4 <- jbs4$id[jbs4$missing==0 | is.na(jbs4$missing)==T]
  here_jbs6 <- jbs6$id[jbs6$missing==0 | is.na(jbs6$missing)==T]
  
  here_any <- union(here_jbs1,union(here_jbs2,union(here_jbs3,union(here_jbs4,here_jbs6))))
  here_any
  
  # The only "true" missings are those that were never here, ie never part of "here_any".
  miss_jbs1 <- miss_jbs1[(miss_jbs1%in%here_any)==F]
  miss_jbs2 <- miss_jbs2[(miss_jbs2%in%here_any)==F]
  miss_jbs3 <- miss_jbs3[(miss_jbs3%in%here_any)==F]
  miss_jbs4 <- miss_jbs4[(miss_jbs4%in%here_any)==F]
  miss_jbs6 <- miss_jbs6[(miss_jbs6%in%here_any)==F]
  
  miss_jbs1
  miss_jbs2
  miss_jbs3
  miss_jbs4
  miss_jbs6 # not much left (if any)
  
  # Now we dont want the id but the row number in ind data :
  miss_jbs1 <- which(jbs1$id%in%miss_jbs1)
  miss_jbs2 <- which(jbs2$id%in%miss_jbs2)
  miss_jbs3 <- which(jbs3$id%in%miss_jbs3)
  miss_jbs4 <- which(jbs4$id%in%miss_jbs4)
  miss_jbs6 <- which(jbs6$id%in%miss_jbs6)
  
  if(length(miss_jbs1)>0){jbs1 <- jbs1[-miss_jbs1,]}
  if(length(miss_jbs2)>0){jbs2 <- jbs2[-miss_jbs2,]}
  if(length(miss_jbs3)>0){jbs3 <- jbs3[-miss_jbs3,]}
  if(length(miss_jbs4)>0){jbs4 <- jbs4[-miss_jbs4,]}
  if(length(miss_jbs6)>0){jbs6 <- jbs6[-miss_jbs6,]}
  
  gf_jbs[[1]] <- network(as.sociomatrix(gf_jbs[[1]])[jbs1$id,jbs1$id])
  gf_jbs[[2]] <- network(as.sociomatrix(gf_jbs[[2]])[jbs2$id,jbs2$id])
  gf_jbs[[3]] <- network(as.sociomatrix(gf_jbs[[3]])[jbs3$id,jbs3$id])
  gf_jbs[[4]] <- network(as.sociomatrix(gf_jbs[[4]])[jbs4$id,jbs4$id])
  gf_jbs[[5]] <- network(as.sociomatrix(gf_jbs[[5]])[jbs6$id,jbs6$id])
  
  f_jbs[[1]] <- network(as.sociomatrix(f_jbs[[1]])[jbs1$id,jbs1$id])
  f_jbs[[2]] <- network(as.sociomatrix(f_jbs[[2]])[jbs2$id,jbs2$id])
  f_jbs[[3]] <- network(as.sociomatrix(f_jbs[[3]])[jbs3$id,jbs3$id])
  f_jbs[[4]] <- network(as.sociomatrix(f_jbs[[4]])[jbs4$id,jbs4$id])
  f_jbs[[5]] <- network(as.sociomatrix(f_jbs[[5]])[jbs6$id,jbs6$id])
  
  allf_jbs[[1]] <- network(as.sociomatrix(allf_jbs[[1]])[jbs1$id,jbs1$id])
  allf_jbs[[2]] <- network(as.sociomatrix(allf_jbs[[2]])[jbs2$id,jbs2$id])
  allf_jbs[[3]] <- network(as.sociomatrix(allf_jbs[[3]])[jbs3$id,jbs3$id])
  allf_jbs[[4]] <- network(as.sociomatrix(allf_jbs[[4]])[jbs4$id,jbs4$id])
  allf_jbs[[5]] <- network(as.sociomatrix(allf_jbs[[5]])[jbs6$id,jbs6$id])
  
  dl_jbs[[1]] <- network(as.sociomatrix(dl_jbs[[1]])[jbs1$id,jbs1$id])
  dl_jbs[[2]] <- network(as.sociomatrix(dl_jbs[[2]])[jbs2$id,jbs2$id])
  dl_jbs[[3]] <- network(as.sociomatrix(dl_jbs[[3]])[jbs3$id,jbs3$id])
  dl_jbs[[4]] <- network(as.sociomatrix(dl_jbs[[4]])[jbs4$id,jbs4$id])
  dl_jbs[[5]] <- network(as.sociomatrix(dl_jbs[[5]])[jbs6$id,jbs6$id])
  
  top5_jbs[[1]] <- network(as.sociomatrix(top5_jbs[[1]])[jbs1$id,jbs1$id])
  top5_jbs[[2]] <- network(as.sociomatrix(top5_jbs[[2]])[jbs2$id,jbs2$id])
  top5_jbs[[3]] <- network(as.sociomatrix(top5_jbs[[3]])[jbs3$id,jbs3$id])
  top5_jbs[[4]] <- network(as.sociomatrix(top5_jbs[[4]])[jbs4$id,jbs4$id])
  top5_jbs[[5]] <- network(as.sociomatrix(top5_jbs[[5]])[jbs6$id,jbs6$id])
  
  bully_jbs[[1]] <- network(as.sociomatrix(bully_jbs[[1]])[jbs1$id,jbs1$id])
  bully_jbs[[2]] <- network(as.sociomatrix(bully_jbs[[2]])[jbs2$id,jbs2$id])
  bully_jbs[[3]] <- network(as.sociomatrix(bully_jbs[[3]])[jbs3$id,jbs3$id])
  bully_jbs[[4]] <- network(as.sociomatrix(bully_jbs[[4]])[jbs4$id,jbs4$id])
  bully_jbs[[5]] <- network(as.sociomatrix(bully_jbs[[5]])[jbs6$id,jbs6$id])
  
  cool_jbs[[1]] <- network(as.sociomatrix(cool_jbs[[1]])[jbs1$id,jbs1$id])
  cool_jbs[[2]] <- network(as.sociomatrix(cool_jbs[[2]])[jbs2$id,jbs2$id])
  cool_jbs[[3]] <- network(as.sociomatrix(cool_jbs[[3]])[jbs3$id,jbs3$id])
  cool_jbs[[4]] <- network(as.sociomatrix(cool_jbs[[4]])[jbs4$id,jbs4$id])
  cool_jbs[[5]] <- network(as.sociomatrix(cool_jbs[[5]])[jbs6$id,jbs6$id])
  
  known_jbs[[2]] <- network(as.sociomatrix(known_jbs[[2]])[jbs2$id,jbs2$id])
  known_jbs[[3]] <- network(as.sociomatrix(known_jbs[[3]])[jbs3$id,jbs3$id])
  known_jbs[[4]] <- network(as.sociomatrix(known_jbs[[4]])[jbs4$id,jbs4$id])
  known_jbs[[5]] <- network(as.sociomatrix(known_jbs[[5]])[jbs6$id,jbs6$id])
  
  fear_jbs[[2]] <- network(as.sociomatrix(fear_jbs[[2]])[jbs2$id,jbs2$id])
  fear_jbs[[3]] <- network(as.sociomatrix(fear_jbs[[3]])[jbs3$id,jbs3$id])
  fear_jbs[[4]] <- network(as.sociomatrix(fear_jbs[[4]])[jbs4$id,jbs4$id])
  fear_jbs[[5]] <- network(as.sociomatrix(fear_jbs[[5]])[jbs6$id,jbs6$id])
  
  liked_jbs[[2]] <- network(as.sociomatrix(liked_jbs[[2]])[jbs2$id,jbs2$id])
  liked_jbs[[3]] <- network(as.sociomatrix(liked_jbs[[3]])[jbs3$id,jbs3$id])
  liked_jbs[[4]] <- network(as.sociomatrix(liked_jbs[[4]])[jbs4$id,jbs4$id])
  liked_jbs[[5]] <- network(as.sociomatrix(liked_jbs[[5]])[jbs6$id,jbs6$id])
  
  ##
  
  miss_pm1 <- pm1$id[pm1$missing==1 & is.na(pm1$missing)==F]
  miss_pm2 <- pm2$id[pm2$missing==1 & is.na(pm2$missing)==F]
  miss_pm3 <- pm3$id[pm3$missing==1 & is.na(pm3$missing)==F]
  miss_pm4 <- pm4$id[pm4$missing==1 & is.na(pm4$missing)==F]
  miss_pm6 <- pm6$id[pm6$missing==1 & is.na(pm6$missing)==F]
  
  here_pm1 <- pm1$id[pm1$missing==0 | is.na(pm1$missing)==T]
  here_pm2 <- pm2$id[pm2$missing==0 | is.na(pm2$missing)==T]
  here_pm3 <- pm3$id[pm3$missing==0 | is.na(pm3$missing)==T]
  here_pm4 <- pm4$id[pm4$missing==0 | is.na(pm4$missing)==T]
  here_pm6 <- pm6$id[pm6$missing==0 | is.na(pm6$missing)==T]
  
  here_any <- union(here_pm1,union(here_pm2,union(here_pm3,union(here_pm4,here_pm6))))
  here_any
  
  # The only "true" missings are those that were never here, ie never part of "here_any".
  miss_pm1 <- miss_pm1[(miss_pm1%in%here_any)==F]
  miss_pm2 <- miss_pm2[(miss_pm2%in%here_any)==F]
  miss_pm3 <- miss_pm3[(miss_pm3%in%here_any)==F]
  miss_pm4 <- miss_pm4[(miss_pm4%in%here_any)==F]
  miss_pm6 <- miss_pm6[(miss_pm6%in%here_any)==F]
  
  miss_pm1
  miss_pm2
  miss_pm3
  miss_pm4
  miss_pm6 # not much left (if any)
  
  # Now we dont want the id but the row number in ind data :
  miss_pm1 <- which(pm1$id%in%miss_pm1)
  miss_pm2 <- which(pm2$id%in%miss_pm2)
  miss_pm3 <- which(pm3$id%in%miss_pm3)
  miss_pm4 <- which(pm4$id%in%miss_pm4)
  miss_pm6 <- which(pm6$id%in%miss_pm6)
  
  if(length(miss_pm1)>0){pm1 <- pm1[-miss_pm1,]}
  if(length(miss_pm2)>0){pm2 <- pm2[-miss_pm2,]}
  if(length(miss_pm3)>0){pm3 <- pm3[-miss_pm3,]}
  if(length(miss_pm4)>0){pm4 <- pm4[-miss_pm4,]}
  if(length(miss_pm6)>0){pm6 <- pm6[-miss_pm6,]}
  
  gf_pm[[1]] <- network(as.sociomatrix(gf_pm[[1]])[pm1$id,pm1$id])
  gf_pm[[2]] <- network(as.sociomatrix(gf_pm[[2]])[pm2$id,pm2$id])
  gf_pm[[3]] <- network(as.sociomatrix(gf_pm[[3]])[pm3$id,pm3$id])
  gf_pm[[4]] <- network(as.sociomatrix(gf_pm[[4]])[pm4$id,pm4$id])
  gf_pm[[5]] <- network(as.sociomatrix(gf_pm[[5]])[pm6$id,pm6$id])
  
  f_pm[[1]] <- network(as.sociomatrix(f_pm[[1]])[pm1$id,pm1$id])
  f_pm[[2]] <- network(as.sociomatrix(f_pm[[2]])[pm2$id,pm2$id])
  f_pm[[3]] <- network(as.sociomatrix(f_pm[[3]])[pm3$id,pm3$id])
  f_pm[[4]] <- network(as.sociomatrix(f_pm[[4]])[pm4$id,pm4$id])
  f_pm[[5]] <- network(as.sociomatrix(f_pm[[5]])[pm6$id,pm6$id])
  
  allf_pm[[1]] <- network(as.sociomatrix(allf_pm[[1]])[pm1$id,pm1$id])
  allf_pm[[2]] <- network(as.sociomatrix(allf_pm[[2]])[pm2$id,pm2$id])
  allf_pm[[3]] <- network(as.sociomatrix(allf_pm[[3]])[pm3$id,pm3$id])
  allf_pm[[4]] <- network(as.sociomatrix(allf_pm[[4]])[pm4$id,pm4$id])
  allf_pm[[5]] <- network(as.sociomatrix(allf_pm[[5]])[pm6$id,pm6$id])
  
  dl_pm[[1]] <- network(as.sociomatrix(dl_pm[[1]])[pm1$id,pm1$id])
  dl_pm[[2]] <- network(as.sociomatrix(dl_pm[[2]])[pm2$id,pm2$id])
  dl_pm[[3]] <- network(as.sociomatrix(dl_pm[[3]])[pm3$id,pm3$id])
  dl_pm[[4]] <- network(as.sociomatrix(dl_pm[[4]])[pm4$id,pm4$id])
  dl_pm[[5]] <- network(as.sociomatrix(dl_pm[[5]])[pm6$id,pm6$id])
  
  top5_pm[[1]] <- network(as.sociomatrix(top5_pm[[1]])[pm1$id,pm1$id])
  top5_pm[[2]] <- network(as.sociomatrix(top5_pm[[2]])[pm2$id,pm2$id])
  top5_pm[[3]] <- network(as.sociomatrix(top5_pm[[3]])[pm3$id,pm3$id])
  top5_pm[[4]] <- network(as.sociomatrix(top5_pm[[4]])[pm4$id,pm4$id])
  top5_pm[[5]] <- network(as.sociomatrix(top5_pm[[5]])[pm6$id,pm6$id])
  
  bully_pm[[1]] <- network(as.sociomatrix(bully_pm[[1]])[pm1$id,pm1$id])
  bully_pm[[2]] <- network(as.sociomatrix(bully_pm[[2]])[pm2$id,pm2$id])
  bully_pm[[3]] <- network(as.sociomatrix(bully_pm[[3]])[pm3$id,pm3$id])
  bully_pm[[4]] <- network(as.sociomatrix(bully_pm[[4]])[pm4$id,pm4$id])
  bully_pm[[5]] <- network(as.sociomatrix(bully_pm[[5]])[pm6$id,pm6$id])
  
  cool_pm[[1]] <- network(as.sociomatrix(cool_pm[[1]])[pm1$id,pm1$id])
  cool_pm[[2]] <- network(as.sociomatrix(cool_pm[[2]])[pm2$id,pm2$id])
  cool_pm[[3]] <- network(as.sociomatrix(cool_pm[[3]])[pm3$id,pm3$id])
  cool_pm[[4]] <- network(as.sociomatrix(cool_pm[[4]])[pm4$id,pm4$id])
  cool_pm[[5]] <- network(as.sociomatrix(cool_pm[[5]])[pm6$id,pm6$id])
  
  known_pm[[2]] <- network(as.sociomatrix(known_pm[[2]])[pm2$id,pm2$id])
  known_pm[[3]] <- network(as.sociomatrix(known_pm[[3]])[pm3$id,pm3$id])
  known_pm[[4]] <- network(as.sociomatrix(known_pm[[4]])[pm4$id,pm4$id])
  known_pm[[5]] <- network(as.sociomatrix(known_pm[[5]])[pm6$id,pm6$id])
  
  fear_pm[[2]] <- network(as.sociomatrix(fear_pm[[2]])[pm2$id,pm2$id])
  fear_pm[[3]] <- network(as.sociomatrix(fear_pm[[3]])[pm3$id,pm3$id])
  fear_pm[[4]] <- network(as.sociomatrix(fear_pm[[4]])[pm4$id,pm4$id])
  fear_pm[[5]] <- network(as.sociomatrix(fear_pm[[5]])[pm6$id,pm6$id])
  
  liked_pm[[2]] <- network(as.sociomatrix(liked_pm[[2]])[pm2$id,pm2$id])
  liked_pm[[3]] <- network(as.sociomatrix(liked_pm[[3]])[pm3$id,pm3$id])
  liked_pm[[4]] <- network(as.sociomatrix(liked_pm[[4]])[pm4$id,pm4$id])
  liked_pm[[5]] <- network(as.sociomatrix(liked_pm[[5]])[pm6$id,pm6$id])
  
  ##
  
  miss_dmz1 <- dmz1$id[dmz1$missing==1 & is.na(dmz1$missing)==F]
  miss_dmz2 <- dmz2$id[dmz2$missing==1 & is.na(dmz2$missing)==F]
  miss_dmz3 <- dmz3$id[dmz3$missing==1 & is.na(dmz3$missing)==F]
  miss_dmz4 <- dmz4$id[dmz4$missing==1 & is.na(dmz4$missing)==F]
  miss_dmz6 <- dmz6$id[dmz6$missing==1 & is.na(dmz6$missing)==F]
  
  here_dmz1 <- dmz1$id[dmz1$missing==0 | is.na(dmz1$missing)==T]
  here_dmz2 <- dmz2$id[dmz2$missing==0 | is.na(dmz2$missing)==T]
  here_dmz3 <- dmz3$id[dmz3$missing==0 | is.na(dmz3$missing)==T]
  here_dmz4 <- dmz4$id[dmz4$missing==0 | is.na(dmz4$missing)==T]
  here_dmz6 <- dmz6$id[dmz6$missing==0 | is.na(dmz6$missing)==T]
  
  here_any <- union(here_dmz1,union(here_dmz2,union(here_dmz3,union(here_dmz4,here_dmz6))))
  here_any
  
  # The only "true" missings are those that were never here, ie never part of "here_any".
  miss_dmz1 <- miss_dmz1[(miss_dmz1%in%here_any)==F]
  miss_dmz2 <- miss_dmz2[(miss_dmz2%in%here_any)==F]
  miss_dmz3 <- miss_dmz3[(miss_dmz3%in%here_any)==F]
  miss_dmz4 <- miss_dmz4[(miss_dmz4%in%here_any)==F]
  miss_dmz6 <- miss_dmz6[(miss_dmz6%in%here_any)==F]
  
  miss_dmz1
  miss_dmz2
  miss_dmz3
  miss_dmz4
  miss_dmz6 # not much left (if any)
  
  # Now we dont want the id but the row number in ind data :
  miss_dmz1 <- which(dmz1$id%in%miss_dmz1)
  miss_dmz2 <- which(dmz2$id%in%miss_dmz2)
  miss_dmz3 <- which(dmz3$id%in%miss_dmz3)
  miss_dmz4 <- which(dmz4$id%in%miss_dmz4)
  miss_dmz6 <- which(dmz6$id%in%miss_dmz6)
  
  if(length(miss_dmz1)>0){dmz1 <- dmz1[-miss_dmz1,]}
  if(length(miss_dmz2)>0){dmz2 <- dmz2[-miss_dmz2,]}
  if(length(miss_dmz3)>0){dmz3 <- dmz3[-miss_dmz3,]}
  if(length(miss_dmz4)>0){dmz4 <- dmz4[-miss_dmz4,]}
  if(length(miss_dmz6)>0){dmz6 <- dmz6[-miss_dmz6,]}
  
  gf_dmz[[1]] <- network(as.sociomatrix(gf_dmz[[1]])[dmz1$id,dmz1$id])
  gf_dmz[[2]] <- network(as.sociomatrix(gf_dmz[[2]])[dmz2$id,dmz2$id])
  gf_dmz[[3]] <- network(as.sociomatrix(gf_dmz[[3]])[dmz3$id,dmz3$id])
  gf_dmz[[4]] <- network(as.sociomatrix(gf_dmz[[4]])[dmz4$id,dmz4$id])
  gf_dmz[[5]] <- network(as.sociomatrix(gf_dmz[[5]])[dmz6$id,dmz6$id])
  
  f_dmz[[1]] <- network(as.sociomatrix(f_dmz[[1]])[dmz1$id,dmz1$id])
  f_dmz[[2]] <- network(as.sociomatrix(f_dmz[[2]])[dmz2$id,dmz2$id])
  f_dmz[[3]] <- network(as.sociomatrix(f_dmz[[3]])[dmz3$id,dmz3$id])
  f_dmz[[4]] <- network(as.sociomatrix(f_dmz[[4]])[dmz4$id,dmz4$id])
  f_dmz[[5]] <- network(as.sociomatrix(f_dmz[[5]])[dmz6$id,dmz6$id])
  
  allf_dmz[[1]] <- network(as.sociomatrix(allf_dmz[[1]])[dmz1$id,dmz1$id])
  allf_dmz[[2]] <- network(as.sociomatrix(allf_dmz[[2]])[dmz2$id,dmz2$id])
  allf_dmz[[3]] <- network(as.sociomatrix(allf_dmz[[3]])[dmz3$id,dmz3$id])
  allf_dmz[[4]] <- network(as.sociomatrix(allf_dmz[[4]])[dmz4$id,dmz4$id])
  allf_dmz[[5]] <- network(as.sociomatrix(allf_dmz[[5]])[dmz6$id,dmz6$id])
  
  dl_dmz[[1]] <- network(as.sociomatrix(dl_dmz[[1]])[dmz1$id,dmz1$id])
  dl_dmz[[2]] <- network(as.sociomatrix(dl_dmz[[2]])[dmz2$id,dmz2$id])
  dl_dmz[[3]] <- network(as.sociomatrix(dl_dmz[[3]])[dmz3$id,dmz3$id])
  dl_dmz[[4]] <- network(as.sociomatrix(dl_dmz[[4]])[dmz4$id,dmz4$id])
  dl_dmz[[5]] <- network(as.sociomatrix(dl_dmz[[5]])[dmz6$id,dmz6$id])
  
  top5_dmz[[1]] <- network(as.sociomatrix(top5_dmz[[1]])[dmz1$id,dmz1$id])
  top5_dmz[[2]] <- network(as.sociomatrix(top5_dmz[[2]])[dmz2$id,dmz2$id])
  top5_dmz[[3]] <- network(as.sociomatrix(top5_dmz[[3]])[dmz3$id,dmz3$id])
  top5_dmz[[4]] <- network(as.sociomatrix(top5_dmz[[4]])[dmz4$id,dmz4$id])
  top5_dmz[[5]] <- network(as.sociomatrix(top5_dmz[[5]])[dmz6$id,dmz6$id])
  
  bully_dmz[[1]] <- network(as.sociomatrix(bully_dmz[[1]])[dmz1$id,dmz1$id])
  bully_dmz[[2]] <- network(as.sociomatrix(bully_dmz[[2]])[dmz2$id,dmz2$id])
  bully_dmz[[3]] <- network(as.sociomatrix(bully_dmz[[3]])[dmz3$id,dmz3$id])
  bully_dmz[[4]] <- network(as.sociomatrix(bully_dmz[[4]])[dmz4$id,dmz4$id])
  bully_dmz[[5]] <- network(as.sociomatrix(bully_dmz[[5]])[dmz6$id,dmz6$id])
  
  cool_dmz[[1]] <- network(as.sociomatrix(cool_dmz[[1]])[dmz1$id,dmz1$id])
  cool_dmz[[2]] <- network(as.sociomatrix(cool_dmz[[2]])[dmz2$id,dmz2$id])
  cool_dmz[[3]] <- network(as.sociomatrix(cool_dmz[[3]])[dmz3$id,dmz3$id])
  cool_dmz[[4]] <- network(as.sociomatrix(cool_dmz[[4]])[dmz4$id,dmz4$id])
  cool_dmz[[5]] <- network(as.sociomatrix(cool_dmz[[5]])[dmz6$id,dmz6$id])
  
  known_dmz[[2]] <- network(as.sociomatrix(known_dmz[[2]])[dmz2$id,dmz2$id])
  known_dmz[[3]] <- network(as.sociomatrix(known_dmz[[3]])[dmz3$id,dmz3$id])
  known_dmz[[4]] <- network(as.sociomatrix(known_dmz[[4]])[dmz4$id,dmz4$id])
  known_dmz[[5]] <- network(as.sociomatrix(known_dmz[[5]])[dmz6$id,dmz6$id])
  
  fear_dmz[[2]] <- network(as.sociomatrix(fear_dmz[[2]])[dmz2$id,dmz2$id])
  fear_dmz[[3]] <- network(as.sociomatrix(fear_dmz[[3]])[dmz3$id,dmz3$id])
  fear_dmz[[4]] <- network(as.sociomatrix(fear_dmz[[4]])[dmz4$id,dmz4$id])
  fear_dmz[[5]] <- network(as.sociomatrix(fear_dmz[[5]])[dmz6$id,dmz6$id])
  
  liked_dmz[[2]] <- network(as.sociomatrix(liked_dmz[[2]])[dmz2$id,dmz2$id])
  liked_dmz[[3]] <- network(as.sociomatrix(liked_dmz[[3]])[dmz3$id,dmz3$id])
  liked_dmz[[4]] <- network(as.sociomatrix(liked_dmz[[4]])[dmz4$id,dmz4$id])
  liked_dmz[[5]] <- network(as.sociomatrix(liked_dmz[[5]])[dmz6$id,dmz6$id])
  
}

##### d. imputed networks #####

if(NA_management=="load_imputed_networks"){
  
  # For now, we only use a single imputed network, the last one (by convention).
  
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w1_rvl.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w2_rvl.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w3_rvl.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w4_rvl.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w6_rvl.RData") 
  
  new_net_rvl_w1 <- network(wave1imp[[length(wave1imp)]])
  new_net_rvl_w2 <- network(wave2imp[[length(wave2imp)]])
  new_net_rvl_w3 <- network(wave3imp[[length(wave3imp)]])
  new_net_rvl_w4 <- network(wave4imp[[length(wave4imp)]])
  new_net_rvl_w6 <- network(wave6imp[[length(wave6imp)]])
  
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w1_jbs.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w2_jbs.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w3_jbs.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w4_jbs.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w6_jbs.RData") 
  
  new_net_jbs_w1 <- network(wave1imp[[length(wave1imp)]])
  new_net_jbs_w2 <- network(wave2imp[[length(wave2imp)]])
  new_net_jbs_w3 <- network(wave3imp[[length(wave3imp)]])
  new_net_jbs_w4 <- network(wave4imp[[length(wave4imp)]])
  new_net_jbs_w6 <- network(wave6imp[[length(wave6imp)]])
  
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w1_pm.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w2_pm.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w3_pm.RData") 
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w4_pm.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w6_pm.RData")
  
  new_net_pm_w1 <- network(wave1imp[[length(wave1imp)]])
  new_net_pm_w2 <- network(wave2imp[[length(wave2imp)]])
  new_net_pm_w3 <- network(wave3imp[[length(wave3imp)]])
  new_net_pm_w4 <- network(wave4imp[[length(wave4imp)]])
  new_net_pm_w6 <- network(wave6imp[[length(wave6imp)]])
  
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w1_dmz.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w2_dmz.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w3_dmz.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w4_dmz.RData")
  load("./data/7 - saved models/fixed data treatment/imputed_networks_gf_w6_dmz.RData")

  new_net_dmz_w1 <- network(wave1imp[[length(wave1imp)]])
  new_net_dmz_w2 <- network(wave2imp[[length(wave2imp)]])
  new_net_dmz_w3 <- network(wave3imp[[length(wave3imp)]])
  new_net_dmz_w4 <- network(wave4imp[[length(wave4imp)]])
  new_net_dmz_w6 <- network(wave6imp[[length(wave6imp)]])
  
  ### A few manual checks (just in case)
  # Should have the exact same nodes as existing networks:
  setdiff(new_net_rvl_w1%v%"vertex.names",gf_rvl[[1]]%v%"vertex.names")
  setdiff(new_net_rvl_w2%v%"vertex.names",gf_rvl[[2]]%v%"vertex.names")
  setdiff(new_net_rvl_w3%v%"vertex.names",gf_rvl[[3]]%v%"vertex.names")
  setdiff(new_net_rvl_w4%v%"vertex.names",gf_rvl[[4]]%v%"vertex.names")
  setdiff(new_net_rvl_w6%v%"vertex.names",gf_rvl[[5]]%v%"vertex.names")
  
  # And the non-NA parts should be identical:
  table(as.sociomatrix(new_net_rvl_w1)[is.na(as.sociomatrix(gf_rvl[[1]]))==F]==as.sociomatrix(gf_rvl[[1]])[is.na(as.sociomatrix(gf_rvl[[1]]))==F])
  table(as.sociomatrix(new_net_rvl_w2)[is.na(as.sociomatrix(gf_rvl[[2]]))==F]==as.sociomatrix(gf_rvl[[2]])[is.na(as.sociomatrix(gf_rvl[[2]]))==F])
  table(as.sociomatrix(new_net_rvl_w3)[is.na(as.sociomatrix(gf_rvl[[3]]))==F]==as.sociomatrix(gf_rvl[[3]])[is.na(as.sociomatrix(gf_rvl[[3]]))==F])
  table(as.sociomatrix(new_net_rvl_w4)[is.na(as.sociomatrix(gf_rvl[[4]]))==F]==as.sociomatrix(gf_rvl[[4]])[is.na(as.sociomatrix(gf_rvl[[4]]))==F])
  table(as.sociomatrix(new_net_rvl_w6)[is.na(as.sociomatrix(gf_rvl[[5]]))==F]==as.sociomatrix(gf_rvl[[5]])[is.na(as.sociomatrix(gf_rvl[[5]]))==F])
  
  # Cool. Now transfer attributes and we have our new objects !
  # NB: actually, for now there is simply no attributes in gf networks, so no need.
  
  ###
  
  gf_rvl <- list(new_net_rvl_w1,new_net_rvl_w2,new_net_rvl_w3,new_net_rvl_w4,new_net_rvl_w6)
  gf_jbs <- list(new_net_jbs_w1,new_net_jbs_w2,new_net_jbs_w3,new_net_jbs_w4,new_net_jbs_w6)
  gf_pm <- list(new_net_pm_w1,new_net_pm_w2,new_net_pm_w3,new_net_pm_w4,new_net_pm_w6)
  gf_dmz <- list(new_net_dmz_w1,new_net_dmz_w2,new_net_dmz_w3,new_net_dmz_w4,new_net_dmz_w6)
  
  rm(new_net_rvl_w1,new_net_rvl_w2,new_net_rvl_w3,new_net_rvl_w4,new_net_rvl_w6)
  rm(new_net_jbs_w1,new_net_jbs_w2,new_net_jbs_w3,new_net_jbs_w4,new_net_jbs_w6)
  rm(new_net_pm_w1,new_net_pm_w2,new_net_pm_w3,new_net_pm_w4,new_net_pm_w6)
  rm(new_net_dmz_w1,new_net_dmz_w2,new_net_dmz_w3,new_net_dmz_w4,new_net_dmz_w6)
  
  rm(wave1imp,wave2imp,wave3imp,wave4imp,wave6imp)
  
}

##### e. lockdown (special w6) #####

# The case of lockdown networks is a bit special because, by construction, it is at the intersection of w4
# and w6 samples.
# w4 people that left in w6 might have been nominated (though probably more rarely), but cannot emit nomination.
# people arrived in w6 should not have nominated, nor have been nominated.

# There are also types 1 and 2 NA, as usual : type 1 being those that were here in w4 and w6, but missed the w6
# questionnaires, and type 2 those that never answered any questionnaire (mostly rvl).

# So, 2 steps : first, narrow down the networks to only the w4/6 intersection.
# Then handle type 1 NAs based on either "remove_all_nas" or "leave_time_specific_nas", just like other
# networks.
# Type 2 NAs are removed either way.

###

# NB : for now, the ids of the matrices are those of both w4 and w6.

# So, the intersection ids :
w46_rvl <- intersect(gf_rvl[[4]]%v%"vertex.names",gf_rvl[[5]]%v%"vertex.names")
w46_jbs <- intersect(gf_jbs[[4]]%v%"vertex.names",gf_jbs[[5]]%v%"vertex.names")
w46_pm <- intersect(gf_pm[[4]]%v%"vertex.names",gf_pm[[5]]%v%"vertex.names")
w46_dmz <- intersect(gf_dmz[[4]]%v%"vertex.names",gf_dmz[[5]]%v%"vertex.names")

# NB : since gf networks have already been corrected from their NAs (with either method),
# using their ids means we also handle type 1 and 2 missings at the same time, so we re happy.

# And resize :
lock_pers_rvl <- network(as.sociomatrix(lock_pers_rvl)[w46_rvl,w46_rvl])
lock_vocal_rvl <- network(as.sociomatrix(lock_vocal_rvl)[w46_rvl,w46_rvl])
lock_write_rvl <- network(as.sociomatrix(lock_write_rvl)[w46_rvl,w46_rvl])
lock_work_rvl <- network(as.sociomatrix(lock_work_rvl)[w46_rvl,w46_rvl])

lock_pers_jbs <- network(as.sociomatrix(lock_pers_jbs)[w46_jbs,w46_jbs])
lock_vocal_jbs <- network(as.sociomatrix(lock_vocal_jbs)[w46_jbs,w46_jbs])
lock_write_jbs <- network(as.sociomatrix(lock_write_jbs)[w46_jbs,w46_jbs])
lock_work_jbs <- network(as.sociomatrix(lock_work_jbs)[w46_jbs,w46_jbs])

lock_pers_pm <- network(as.sociomatrix(lock_pers_pm)[w46_pm,w46_pm])
lock_vocal_pm <- network(as.sociomatrix(lock_vocal_pm)[w46_pm,w46_pm])
lock_write_pm <- network(as.sociomatrix(lock_write_pm)[w46_pm,w46_pm])
lock_work_pm <- network(as.sociomatrix(lock_work_pm)[w46_pm,w46_pm])

lock_pers_dmz <- network(as.sociomatrix(lock_pers_dmz)[w46_dmz,w46_dmz])
lock_vocal_dmz <- network(as.sociomatrix(lock_vocal_dmz)[w46_dmz,w46_dmz])
lock_write_dmz <- network(as.sociomatrix(lock_write_dmz)[w46_dmz,w46_dmz])
lock_work_dmz <- network(as.sociomatrix(lock_work_dmz)[w46_dmz,w46_dmz])
  
##### f. individual attributes #####

# New from may 2021.

# The following function takes a data frame, find all variables with name "imputed" and "imputed_type" in it
# (the variable names used for imputed attributes), extract the corresponding variable names and proceed to transfer
# the values following either the "soft" or "hard" rule.
my_impute_function <- function(mydata, soft_or_hard){
  
  # NB : "hard" means we should impute both soft and hard types.
  if(soft_or_hard=="hard"){soft_or_hard <- c("soft","hard")}
  
  var_names <- names(mydata)[sapply(strsplit(names(mydata),"_"),"[",1)=="imputed"]
  var_names <- var_names[sapply(strsplit(var_names,"_"),"[",2)!="type"]
  original_vars <- sapply(strsplit(var_names,"mputed_"),"[",2)
  
  if(length(original_vars)>0){
    for(var in original_vars){
      aa <- paste("imputed",var,sep="_")
      bb <- paste("imputed","type",var,sep="_")
      mydata[mydata[,bb]%in%soft_or_hard & is.na(mydata[,bb])==F,var] <- mydata[mydata[,bb]%in%soft_or_hard & is.na(mydata[,bb])==F,aa]
    }
    
  }
  
  return(mydata)
}

# If "NA_impute" is empty, nothing happens.
# If it s "hard" or "soft", we start the procedure.

if(NA_impute%in%c("soft","hard")){
  
  # NB: if there is nothing to impute in a data frame, the function does nothing.
  
  occupation <- my_impute_function(occupation,NA_impute)
  fxd <- my_impute_function(fxd,NA_impute)
  
  rvl1 <- my_impute_function(rvl1,NA_impute)
  rvl2 <- my_impute_function(rvl2,NA_impute)
  rvl3 <- my_impute_function(rvl3,NA_impute)
  rvl4 <- my_impute_function(rvl4,NA_impute)
  rvl6 <- my_impute_function(rvl6,NA_impute) 
  
  jbs1 <- my_impute_function(jbs1,NA_impute)
  jbs2 <- my_impute_function(jbs2,NA_impute)
  jbs3 <- my_impute_function(jbs3,NA_impute)
  jbs4 <- my_impute_function(jbs4,NA_impute)
  jbs6 <- my_impute_function(jbs6,NA_impute) 
  
  pm1 <- my_impute_function(pm1,NA_impute)
  pm2 <- my_impute_function(pm2,NA_impute)
  pm3 <- my_impute_function(pm3,NA_impute)
  pm4 <- my_impute_function(pm4,NA_impute)
  pm6 <- my_impute_function(pm6,NA_impute) 
  
  dmz1 <- my_impute_function(dmz1,NA_impute)
  dmz2 <- my_impute_function(dmz2,NA_impute)
  dmz3 <- my_impute_function(dmz3,NA_impute)
  dmz4 <- my_impute_function(dmz4,NA_impute)
  dmz6 <- my_impute_function(dmz6,NA_impute) 
  
}

# A quick check just in case : we know for sure that "factbrev" at elast has imputed values in occupation.
# If it s not there, it means the script "impute missing data - ind" has not been used to create the imputed values in the data
# frames (ie the script creating ind data frames has been re-run, creating new data frames without imputed values, and then the imputation
# script has not been run).
# In this case, calling NA_impute soft or hard should return an error, because it means there is nothing to impute.

if(NA_impute%in%c("soft","hard")){
  
  if("imputed_factbrev"%in%names(occupation)==F){
    print("ERROR. IMPUTED VALUES ARE NOT AVALAIBLE FOR INDIVIDUAL ATTRIBUTES.THE SCRIPT *impute missing data - ind* HAS PROBABLY NOT BEEN RUN")
    break
  }
  
}

########## 4. Out-of-school networks (optional) ########## 

if(process_out_net==TRUE){
  
  ## Load Out-of-school networks (not done in the general setup)
  
  load("./data/5 - final data/net/myhome_rvl.RData")
  load("./data/5 - final data/net/hishome_rvl.RData")
  load("./data/5 - final data/net/outside_rvl.RData")
  load("./data/5 - final data/net/parents_rvl.RData")
  
  load("./data/5 - final data/net/myhome_jbs.RData")
  load("./data/5 - final data/net/hishome_jbs.RData")
  load("./data/5 - final data/net/outside_jbs.RData")
  load("./data/5 - final data/net/parents_jbs.RData")
  
  load("./data/5 - final data/net/myhome_pm.RData")
  load("./data/5 - final data/net/hishome_pm.RData")
  load("./data/5 - final data/net/outside_pm.RData")
  load("./data/5 - final data/net/parents_pm.RData")
  
  load("./data/5 - final data/net/myhome_dmz.RData")
  load("./data/5 - final data/net/hishome_dmz.RData")
  load("./data/5 - final data/net/outside_dmz.RData")
  load("./data/5 - final data/net/parents_dmz.RData")
  
  # Remove NAs (same rule than for the rest) :
  myhome_rvl[[1]] <- network(as.sociomatrix(myhome_rvl[[1]])[rvl1$id,rvl1$id])
  myhome_rvl[[2]] <- network(as.sociomatrix(myhome_rvl[[2]])[rvl3$id,rvl3$id])
  hishome_rvl[[1]] <- network(as.sociomatrix(hishome_rvl[[1]])[rvl1$id,rvl1$id])
  hishome_rvl[[2]] <- network(as.sociomatrix(hishome_rvl[[2]])[rvl3$id,rvl3$id])
  outside_rvl[[1]] <- network(as.sociomatrix(outside_rvl[[1]])[rvl1$id,rvl1$id])
  outside_rvl[[2]] <- network(as.sociomatrix(outside_rvl[[2]])[rvl3$id,rvl3$id])
  parents_rvl[[1]] <- network(as.sociomatrix(parents_rvl[[1]])[rvl1$id,rvl1$id])
  parents_rvl[[2]] <- network(as.sociomatrix(parents_rvl[[2]])[rvl3$id,rvl3$id])
  
  myhome_jbs[[1]] <- network(as.sociomatrix(myhome_jbs[[1]])[jbs1$id,jbs1$id])
  myhome_jbs[[2]] <- network(as.sociomatrix(myhome_jbs[[2]])[jbs3$id,jbs3$id])
  hishome_jbs[[1]] <- network(as.sociomatrix(hishome_jbs[[1]])[jbs1$id,jbs1$id])
  hishome_jbs[[2]] <- network(as.sociomatrix(hishome_jbs[[2]])[jbs3$id,jbs3$id])
  outside_jbs[[1]] <- network(as.sociomatrix(outside_jbs[[1]])[jbs1$id,jbs1$id])
  outside_jbs[[2]] <- network(as.sociomatrix(outside_jbs[[2]])[jbs3$id,jbs3$id])
  parents_jbs[[1]] <- network(as.sociomatrix(parents_jbs[[1]])[jbs1$id,jbs1$id])
  parents_jbs[[2]] <- network(as.sociomatrix(parents_jbs[[2]])[jbs3$id,jbs3$id])
  
  myhome_pm[[1]] <- network(as.sociomatrix(myhome_pm[[1]])[pm1$id,pm1$id])
  myhome_pm[[2]] <- network(as.sociomatrix(myhome_pm[[2]])[pm3$id,pm3$id])
  hishome_pm[[1]] <- network(as.sociomatrix(hishome_pm[[1]])[pm1$id,pm1$id])
  hishome_pm[[2]] <- network(as.sociomatrix(hishome_pm[[2]])[pm3$id,pm3$id])
  outside_pm[[1]] <- network(as.sociomatrix(outside_pm[[1]])[pm1$id,pm1$id])
  outside_pm[[2]] <- network(as.sociomatrix(outside_pm[[2]])[pm3$id,pm3$id])
  parents_pm[[1]] <- network(as.sociomatrix(parents_pm[[1]])[pm1$id,pm1$id])
  parents_pm[[2]] <- network(as.sociomatrix(parents_pm[[2]])[pm3$id,pm3$id])
  
  myhome_dmz[[1]] <- network(as.sociomatrix(myhome_dmz[[1]])[dmz1$id,dmz1$id])
  myhome_dmz[[2]] <- network(as.sociomatrix(myhome_dmz[[2]])[dmz3$id,dmz3$id])
  hishome_dmz[[1]] <- network(as.sociomatrix(hishome_dmz[[1]])[dmz1$id,dmz1$id])
  hishome_dmz[[2]] <- network(as.sociomatrix(hishome_dmz[[2]])[dmz3$id,dmz3$id])
  outside_dmz[[1]] <- network(as.sociomatrix(outside_dmz[[1]])[dmz1$id,dmz1$id])
  outside_dmz[[2]] <- network(as.sociomatrix(outside_dmz[[2]])[dmz3$id,dmz3$id])
  parents_dmz[[1]] <- network(as.sociomatrix(parents_dmz[[1]])[dmz1$id,dmz1$id])
  parents_dmz[[2]] <- network(as.sociomatrix(parents_dmz[[2]])[dmz3$id,dmz3$id])
  
}

# Now we don t want the miss objects anymore.
rm(miss_rvl1,miss_rvl2,miss_rvl3,
   miss_jbs1,miss_jbs2,miss_jbs3,
   miss_pm1,miss_pm2,miss_pm3,
   miss_dmz1,miss_dmz2,miss_dmz3) 

########## 5. Place individual attributes in the network objects ########## 
##### a. "Untemporal" information - transfer from "occupation" #####

# These are variables that are fixed attributes of students (e.g. gender) and thus can be passed
# similarly to all waves and networks.

names(occupation)

# For having more convenient coefficients to interpret, isei2 is standardized (in occupation).
# We save the non-standardized version though.

occupation$isei2_nonstd <- occupation$isei2

occupation$isei2 <- (occupation$isei2 - mean(occupation$isei2,na.rm=T))/sd(occupation$isei2,na.rm=T)
occupation$isei2 <- round(occupation$isei2,digits = 3)

a <- c("factbrev","pcs_chef","isei2","pcs_chef_ag") # name in occupation
b <- c("factbrev","pcs_chef","isei2","pcs_chef_ag") # name we want in the network objects

# I took off suli and isei for now, they haven t been coded using w4 info yet. Just add them to a and b if needed.

# Use imputed values or not ? (usually, they re not precise at all -.-).


for(i in 1:length(gf_rvl)){
  gf_rvl[[i]] <- attribute_in_network(net = gf_rvl[[i]], data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
  
  gf_jbs[[i]] <- attribute_in_network(net = gf_jbs[[i]], data = occupation, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  gf_pm[[i]] <- attribute_in_network(net = gf_pm[[i]], data = occupation, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  gf_dmz[[i]] <- attribute_in_network(net = gf_dmz[[i]], data = occupation, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  f_rvl[[i]] <- attribute_in_network(net = f_rvl[[i]], data = occupation, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)
  
  f_jbs[[i]] <- attribute_in_network(net = f_jbs[[i]], data = occupation, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  f_pm[[i]] <- attribute_in_network(net = f_pm[[i]], data = occupation, id_data = c("id"), 
                                    var_names = a,
                                    new_names = b)
  
  f_dmz[[i]] <- attribute_in_network(net = f_dmz[[i]], data = occupation, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  allf_rvl[[i]] <- attribute_in_network(net = allf_rvl[[i]], data = occupation, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  allf_jbs[[i]] <- attribute_in_network(net = allf_jbs[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  allf_pm[[i]] <- attribute_in_network(net = allf_pm[[i]], data = occupation, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  allf_dmz[[i]] <- attribute_in_network(net = allf_dmz[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  dl_rvl[[i]] <- attribute_in_network(net = dl_rvl[[i]], data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
  
  dl_jbs[[i]] <- attribute_in_network(net = dl_jbs[[i]], data = occupation, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  dl_pm[[i]] <- attribute_in_network(net = dl_pm[[i]], data = occupation, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  dl_dmz[[i]] <- attribute_in_network(net = dl_dmz[[i]], data = occupation, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  top5_rvl[[i]] <- attribute_in_network(net = top5_rvl[[i]], data = occupation, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  top5_jbs[[i]] <- attribute_in_network(net = top5_jbs[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  top5_pm[[i]] <- attribute_in_network(net = top5_pm[[i]], data = occupation, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  top5_dmz[[i]] <- attribute_in_network(net = top5_dmz[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  bully_rvl[[i]] <- attribute_in_network(net = bully_rvl[[i]], data = occupation, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  bully_jbs[[i]] <- attribute_in_network(net = bully_jbs[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  bully_pm[[i]] <- attribute_in_network(net = bully_pm[[i]], data = occupation, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  bully_dmz[[i]] <- attribute_in_network(net = bully_dmz[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  cool_rvl[[i]] <- attribute_in_network(net = cool_rvl[[i]], data = occupation, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  cool_jbs[[i]] <- attribute_in_network(net = cool_jbs[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  cool_pm[[i]] <- attribute_in_network(net = cool_pm[[i]], data = occupation, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  cool_dmz[[i]] <- attribute_in_network(net = cool_dmz[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
}

for(i in 2:length(gf_rvl)){
  known_rvl[[i]] <- attribute_in_network(net = known_rvl[[i]], data = occupation, id_data = c("id"),
                                         var_names = a,
                                         new_names = b)
  
  known_jbs[[i]] <- attribute_in_network(net = known_jbs[[i]], data = occupation, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  known_pm[[i]] <- attribute_in_network(net = known_pm[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  known_dmz[[i]] <- attribute_in_network(net = known_dmz[[i]], data = occupation, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  fear_rvl[[i]] <- attribute_in_network(net = fear_rvl[[i]], data = occupation, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  fear_jbs[[i]] <- attribute_in_network(net = fear_jbs[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  fear_pm[[i]] <- attribute_in_network(net = fear_pm[[i]], data = occupation, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  fear_dmz[[i]] <- attribute_in_network(net = fear_dmz[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  liked_rvl[[i]] <- attribute_in_network(net = liked_rvl[[i]], data = occupation, id_data = c("id"),
                                         var_names = a,
                                         new_names = b)
  
  liked_jbs[[i]] <- attribute_in_network(net = liked_jbs[[i]], data = occupation, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  liked_pm[[i]] <- attribute_in_network(net = liked_pm[[i]], data = occupation, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  liked_dmz[[i]] <- attribute_in_network(net = liked_dmz[[i]], data = occupation, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
}

# Lockdown on its own :
lock_pers_rvl <- attribute_in_network(net = lock_pers_rvl, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_rvl <- attribute_in_network(net = lock_vocal_rvl, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_write_rvl <- attribute_in_network(net = lock_write_rvl, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_work_rvl <- attribute_in_network(net = lock_work_rvl, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

lock_pers_jbs <- attribute_in_network(net = lock_pers_jbs, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_jbs <- attribute_in_network(net = lock_vocal_jbs, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_jbs <- attribute_in_network(net = lock_write_jbs, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_jbs <- attribute_in_network(net = lock_work_jbs, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

lock_pers_pm <- attribute_in_network(net = lock_pers_pm, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_pm <- attribute_in_network(net = lock_vocal_pm, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_pm <- attribute_in_network(net = lock_write_pm, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_pm <- attribute_in_network(net = lock_work_pm, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

lock_pers_dmz <- attribute_in_network(net = lock_pers_dmz, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_dmz <- attribute_in_network(net = lock_vocal_dmz, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_dmz <- attribute_in_network(net = lock_write_dmz, data = occupation, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_dmz <- attribute_in_network(net = lock_work_dmz, data = occupation, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

##### b. "Untemporal" information - transfer from "fxd" #####

# NB : I use "primary1", i.e the last primary school the kid was in (when he declared several).
# NB2 : "single" means the kid was in a primary school with no one else.

# Also transfer OIB (appart from rvl it s all NAs).

# We also include students first name.
names(fxd)

a <- c("sex","primary1","name","ethn","ethn_ag","ethn5","oib")
b <- c("sex","primary1","name","ethn","ethn_ag","ethn5","oib")

for(i in 1:length(gf_rvl)){
  gf_rvl[[i]] <- attribute_in_network(net = gf_rvl[[i]], data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
  
  gf_jbs[[i]] <- attribute_in_network(net = gf_jbs[[i]], data = fxd, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  gf_pm[[i]] <- attribute_in_network(net = gf_pm[[i]], data = fxd, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  gf_dmz[[i]] <- attribute_in_network(net = gf_dmz[[i]], data = fxd, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  #
  
  f_rvl[[i]] <- attribute_in_network(net = f_rvl[[i]], data = fxd, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)
  
  f_jbs[[i]] <- attribute_in_network(net = f_jbs[[i]], data = fxd, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  f_pm[[i]] <- attribute_in_network(net = f_pm[[i]], data = fxd, id_data = c("id"), 
                                    var_names = a,
                                    new_names = b)
  
  f_dmz[[i]] <- attribute_in_network(net = f_dmz[[i]], data = fxd, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  #
  
  allf_rvl[[i]] <- attribute_in_network(net = allf_rvl[[i]], data = fxd, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  allf_jbs[[i]] <- attribute_in_network(net = allf_jbs[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  allf_pm[[i]] <- attribute_in_network(net = allf_pm[[i]], data = fxd, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  allf_dmz[[i]] <- attribute_in_network(net = allf_dmz[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  #
  
  dl_rvl[[i]] <- attribute_in_network(net = dl_rvl[[i]], data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
  
  dl_jbs[[i]] <- attribute_in_network(net = dl_jbs[[i]], data = fxd, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  dl_pm[[i]] <- attribute_in_network(net = dl_pm[[i]], data = fxd, id_data = c("id"), 
                                     var_names = a,
                                     new_names = b)
  
  dl_dmz[[i]] <- attribute_in_network(net = dl_dmz[[i]], data = fxd, id_data = c("id"), 
                                      var_names = a,
                                      new_names = b)
  
  #
  
  top5_rvl[[i]] <- attribute_in_network(net = top5_rvl[[i]], data = fxd, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  top5_jbs[[i]] <- attribute_in_network(net = top5_jbs[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  top5_pm[[i]] <- attribute_in_network(net = top5_pm[[i]], data = fxd, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  top5_dmz[[i]] <- attribute_in_network(net = top5_dmz[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  #
  
  bully_rvl[[i]] <- attribute_in_network(net = bully_rvl[[i]], data = fxd, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  bully_jbs[[i]] <- attribute_in_network(net = bully_jbs[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  bully_pm[[i]] <- attribute_in_network(net = bully_pm[[i]], data = fxd, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  bully_dmz[[i]] <- attribute_in_network(net = bully_dmz[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  #
  
  cool_rvl[[i]] <- attribute_in_network(net = cool_rvl[[i]], data = fxd, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  cool_jbs[[i]] <- attribute_in_network(net = cool_jbs[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  cool_pm[[i]] <- attribute_in_network(net = cool_pm[[i]], data = fxd, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  cool_dmz[[i]] <- attribute_in_network(net = cool_dmz[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
}

for(i in 2:length(gf_rvl)){
  known_rvl[[i]] <- attribute_in_network(net = known_rvl[[i]], data = fxd, id_data = c("id"),
                                         var_names = a,
                                         new_names = b)
  
  known_jbs[[i]] <- attribute_in_network(net = known_jbs[[i]], data = fxd, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  known_pm[[i]] <- attribute_in_network(net = known_pm[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  known_dmz[[i]] <- attribute_in_network(net = known_dmz[[i]], data = fxd, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  #
  
  fear_rvl[[i]] <- attribute_in_network(net = fear_rvl[[i]], data = fxd, id_data = c("id"),
                                        var_names = a,
                                        new_names = b)
  
  fear_jbs[[i]] <- attribute_in_network(net = fear_jbs[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  fear_pm[[i]] <- attribute_in_network(net = fear_pm[[i]], data = fxd, id_data = c("id"), 
                                       var_names = a,
                                       new_names = b)
  
  fear_dmz[[i]] <- attribute_in_network(net = fear_dmz[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  #
  
  liked_rvl[[i]] <- attribute_in_network(net = liked_rvl[[i]], data = fxd, id_data = c("id"),
                                         var_names = a,
                                         new_names = b)
  
  liked_jbs[[i]] <- attribute_in_network(net = liked_jbs[[i]], data = fxd, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
  
  liked_pm[[i]] <- attribute_in_network(net = liked_pm[[i]], data = fxd, id_data = c("id"), 
                                        var_names = a,
                                        new_names = b)
  
  liked_dmz[[i]] <- attribute_in_network(net = liked_dmz[[i]], data = fxd, id_data = c("id"), 
                                         var_names = a,
                                         new_names = b)
}

# Lockdown on its own :
lock_pers_rvl <- attribute_in_network(net = lock_pers_rvl, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_rvl <- attribute_in_network(net = lock_vocal_rvl, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_rvl <- attribute_in_network(net = lock_write_rvl, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_rvl <- attribute_in_network(net = lock_work_rvl, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

lock_pers_jbs <- attribute_in_network(net = lock_pers_jbs, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_jbs <- attribute_in_network(net = lock_vocal_jbs, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_jbs <- attribute_in_network(net = lock_write_jbs, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_jbs <- attribute_in_network(net = lock_work_jbs, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

lock_pers_pm <- attribute_in_network(net = lock_pers_pm, data = fxd, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)
lock_vocal_pm <- attribute_in_network(net = lock_vocal_pm, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_write_pm <- attribute_in_network(net = lock_write_pm, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_work_pm <- attribute_in_network(net = lock_work_pm, data = fxd, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

lock_pers_dmz <- attribute_in_network(net = lock_pers_dmz, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)
lock_vocal_dmz <- attribute_in_network(net = lock_vocal_dmz, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_write_dmz <- attribute_in_network(net = lock_write_dmz, data = fxd, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)
lock_work_dmz <- attribute_in_network(net = lock_work_dmz, data = fxd, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

##### c. "Temporal" information - some corrections (grades) #####

# For grades, the variable names are all over the place : it can be "g_mean_T1" or "T_2" depending on
# the wave, and it s "g_ACP_Tx" for dmz.
# So in each data frames, I will create a simple "grade" variable that gives the appropriate info (i.e the right trimester
# and variable names).

# Each variable will be standardized so I can compare them across schools.

# In the data frames with no grade info, I will still create an empty NA variable that will be passed to the networks (this
# way I won t have to change too many things when I finally get the info).

# Finally, for demotz w1, I lack the T2 info for some kids so I will use T1 instead.

rvl1$grade <- NA
rvl2$grade <- NA
rvl3$grade <- NA
rvl4$grade <- NA
rvl6$grade <- NA

jbs1$grade <- NA
jbs2$grade <- NA
jbs3$grade <- NA
jbs4$grade <- NA
jbs6$grade <- NA

pm1$grade <- NA
pm2$grade <- NA
pm3$grade <- NA
pm4$grade <- NA
pm6$grade <- NA

dmz1$grade <- NA
dmz2$grade <- NA
dmz3$grade <- NA
dmz4$grade <- NA
dmz6$grade <- NA

# Transfer and standardize :
rvl1$grade <- (rvl1$g_mean_T2 - mean(rvl1$g_mean_T2,na.rm=T))/sd(rvl1$g_mean_T2,na.rm=T)
jbs1$grade <- (jbs1$g_mean_T2 - mean(jbs1$g_mean_T2,na.rm=T))/sd(jbs1$g_mean_T2,na.rm=T)
pm1$grade <- (pm1$g_mean_T2 - mean(pm1$g_mean_T2,na.rm=T))/sd(pm1$g_mean_T2,na.rm=T)
dmz1$grade <- (dmz1$g_acp_T2 - mean(dmz1$g_acp_T2,na.rm=T))/sd(dmz1$g_acp_T2,na.rm=T)

rvl2$grade <- (rvl2$g_mean_T1 - mean(rvl2$g_mean_T1,na.rm=T))/sd(rvl2$g_mean_T1,na.rm=T)
#jbs2$grade <- (jbs2$g_mean_T1 - mean(jbs2$g_mean_T1,na.rm=T))/sd(jbs2$g_mean_T1,na.rm=T)
pm2$grade <- (pm2$g_mean_T1 - mean(pm2$g_mean_T1,na.rm=T))/sd(pm2$g_mean_T1,na.rm=T)
dmz2$grade <- (dmz2$g_acp_T1 - mean(dmz2$g_acp_T1,na.rm=T))/sd(dmz2$g_acp_T1,na.rm=T)

rvl3$grade <- (rvl3$g_mean_T2 - mean(rvl3$g_mean_T2,na.rm=T))/sd(rvl3$g_mean_T2,na.rm=T)
#jbs3$grade <- (jbs3$g_mean_T2 - mean(jbs3$g_mean_T2,na.rm=T))/sd(jbs3$g_mean_T2,na.rm=T)
pm3$grade <- (pm3$g_mean_T2 - mean(pm3$g_mean_T2,na.rm=T))/sd(pm3$g_mean_T2,na.rm=T)
dmz3$grade <- (dmz3$g_acp_T2 - mean(dmz3$g_acp_T2,na.rm=T))/sd(dmz3$g_acp_T2,na.rm=T)

rvl4$grade <- (rvl4$g_mean_T1 - mean(rvl4$g_mean_T1,na.rm=T))/sd(rvl4$g_mean_T1,na.rm=T)
#jbs4$grade <- (jbs4$g_mean_T1 - mean(jbs4$g_mean_T1,na.rm=T))/sd(jbs4$g_mean_T1,na.rm=T)
pm4$grade <- (pm4$g_mean_T1 - mean(pm4$g_mean_T1,na.rm=T))/sd(pm4$g_mean_T1,na.rm=T)
dmz4$grade <- (dmz4$g_acm_T1 - mean(dmz4$g_acm_T1,na.rm=T))/sd(dmz4$g_acm_T1,na.rm=T)

rvl6$grade <- (rvl6$g_mean_T1 - mean(rvl6$g_mean_T1,na.rm=T))/sd(rvl4$g_mean_T1,na.rm=T)
#jbs6$grade <- (jbs6$g_mean_T1 - mean(jbs6$g_mean_T1,na.rm=T))/sd(jbs6$g_mean_T1,na.rm=T)
pm6$grade <- (pm6$g_mean_T1 - mean(pm6$g_mean_T1,na.rm=T))/sd(pm6$g_mean_T1,na.rm=T)
dmz6$grade <- (dmz6$g_acm_T1 - mean(dmz6$g_acm_T1,na.rm=T))/sd(dmz6$g_acm_T1,na.rm=T)

# For dmz only we correct T2 with T1 in w1 :
dmz1$g_acp_T1[is.na(dmz1$g_acp_T2)] # Well, no case where it would help anyway...

##### d. "Temporal" information - transfer from ind data frames #####

rvl6$discipline <- NA
jbs6$discipline <- NA
pm6$discipline <- NA
dmz6$discipline <- NA # so the loop find those, even if they re empty.

a <- c("classroom","grade","discipline","factatt")
b <- c("classroom","grade","discipline","factatt")

gf_rvl[[1]] <- attribute_in_network(net = gf_rvl[[1]], data = rvl1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_rvl[[2]] <- attribute_in_network(net = gf_rvl[[2]], data = rvl2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_rvl[[3]] <- attribute_in_network(net = gf_rvl[[3]], data = rvl3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_rvl[[4]] <- attribute_in_network(net = gf_rvl[[4]], data = rvl4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_rvl[[5]] <- attribute_in_network(net = gf_rvl[[5]], data = rvl6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

#

gf_jbs[[1]] <- attribute_in_network(net = gf_jbs[[1]], data = jbs1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_jbs[[2]] <- attribute_in_network(net = gf_jbs[[2]], data = jbs2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_jbs[[3]] <- attribute_in_network(net = gf_jbs[[3]], data = jbs3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_jbs[[4]] <- attribute_in_network(net = gf_jbs[[4]], data = jbs4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_jbs[[5]] <- attribute_in_network(net = gf_jbs[[5]], data = jbs6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

#

gf_pm[[1]] <- attribute_in_network(net = gf_pm[[1]], data = pm1, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

gf_pm[[2]] <- attribute_in_network(net = gf_pm[[2]], data = pm2, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

gf_pm[[3]] <- attribute_in_network(net = gf_pm[[3]], data = pm3, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

gf_pm[[4]] <- attribute_in_network(net = gf_pm[[4]], data = pm4, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

gf_pm[[5]] <- attribute_in_network(net = gf_pm[[5]], data = pm6, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

#

gf_dmz[[1]] <- attribute_in_network(net = gf_dmz[[1]], data = dmz1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_dmz[[2]] <- attribute_in_network(net = gf_dmz[[2]], data = dmz2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_dmz[[3]] <- attribute_in_network(net = gf_dmz[[3]], data = dmz3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_dmz[[4]] <- attribute_in_network(net = gf_dmz[[4]], data = dmz4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

gf_dmz[[5]] <- attribute_in_network(net = gf_dmz[[5]], data = dmz6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

###

f_rvl[[1]] <- attribute_in_network(net = f_rvl[[1]], data = rvl1, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_rvl[[2]] <- attribute_in_network(net = f_rvl[[2]], data = rvl2, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_rvl[[3]] <- attribute_in_network(net = f_rvl[[3]], data = rvl3, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_rvl[[4]] <- attribute_in_network(net = f_rvl[[4]], data = rvl4, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_rvl[[5]] <- attribute_in_network(net = f_rvl[[5]], data = rvl6, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

#

f_jbs[[1]] <- attribute_in_network(net = f_jbs[[1]], data = jbs1, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_jbs[[2]] <- attribute_in_network(net = f_jbs[[2]], data = jbs2, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_jbs[[3]] <- attribute_in_network(net = f_jbs[[3]], data = jbs3, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_jbs[[4]] <- attribute_in_network(net = f_jbs[[4]], data = jbs4, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_jbs[[5]] <- attribute_in_network(net = f_jbs[[5]], data = jbs6, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

#

f_pm[[1]] <- attribute_in_network(net = f_pm[[1]], data = pm1, id_data = c("id"),
                                  var_names = a,
                                  new_names = b)

f_pm[[2]] <- attribute_in_network(net = f_pm[[2]], data = pm2, id_data = c("id"),
                                  var_names = a,
                                  new_names = b)

f_pm[[3]] <- attribute_in_network(net = f_pm[[3]], data = pm3, id_data = c("id"),
                                  var_names = a,
                                  new_names = b)

f_pm[[4]] <- attribute_in_network(net = f_pm[[4]], data = pm4, id_data = c("id"),
                                  var_names = a,
                                  new_names = b)

f_pm[[5]] <- attribute_in_network(net = f_pm[[5]], data = pm6, id_data = c("id"),
                                  var_names = a,
                                  new_names = b)

#

f_dmz[[1]] <- attribute_in_network(net = f_dmz[[1]], data = dmz1, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_dmz[[2]] <- attribute_in_network(net = f_dmz[[2]], data = dmz2, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_dmz[[3]] <- attribute_in_network(net = f_dmz[[3]], data = dmz3, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_dmz[[4]] <- attribute_in_network(net = f_dmz[[4]], data = dmz4, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

f_dmz[[5]] <- attribute_in_network(net = f_dmz[[5]], data = dmz6, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

###

allf_rvl[[1]] <- attribute_in_network(net = allf_rvl[[1]], data = rvl1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_rvl[[2]] <- attribute_in_network(net = allf_rvl[[2]], data = rvl2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_rvl[[3]] <- attribute_in_network(net = allf_rvl[[3]], data = rvl3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_rvl[[4]] <- attribute_in_network(net = allf_rvl[[4]], data = rvl4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_rvl[[5]] <- attribute_in_network(net = allf_rvl[[5]], data = rvl6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

allf_jbs[[1]] <- attribute_in_network(net = allf_jbs[[1]], data = jbs1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_jbs[[2]] <- attribute_in_network(net = allf_jbs[[2]], data = jbs2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_jbs[[3]] <- attribute_in_network(net = allf_jbs[[3]], data = jbs3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_jbs[[4]] <- attribute_in_network(net = allf_jbs[[4]], data = jbs4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_jbs[[5]] <- attribute_in_network(net = allf_jbs[[5]], data = jbs6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)


#

allf_pm[[1]] <- attribute_in_network(net = allf_pm[[1]], data = pm1, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

allf_pm[[2]] <- attribute_in_network(net = allf_pm[[2]], data = pm2, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

allf_pm[[3]] <- attribute_in_network(net = allf_pm[[3]], data = pm3, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

allf_pm[[4]] <- attribute_in_network(net = allf_pm[[4]], data = pm4, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

allf_pm[[5]] <- attribute_in_network(net = allf_pm[[5]], data = pm6, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

#

allf_dmz[[1]] <- attribute_in_network(net = allf_dmz[[1]], data = dmz1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_dmz[[2]] <- attribute_in_network(net = allf_dmz[[2]], data = dmz2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_dmz[[3]] <- attribute_in_network(net = allf_dmz[[3]], data = dmz3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_dmz[[4]] <- attribute_in_network(net = allf_dmz[[4]], data = dmz4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

allf_dmz[[5]] <- attribute_in_network(net = allf_dmz[[5]], data = dmz6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

###

dl_rvl[[1]] <- attribute_in_network(net = dl_rvl[[1]], data = rvl1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_rvl[[2]] <- attribute_in_network(net = dl_rvl[[2]], data = rvl2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_rvl[[3]] <- attribute_in_network(net = dl_rvl[[3]], data = rvl3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_rvl[[4]] <- attribute_in_network(net = dl_rvl[[4]], data = rvl4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_rvl[[5]] <- attribute_in_network(net = dl_rvl[[5]], data = rvl6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

#

dl_jbs[[1]] <- attribute_in_network(net = dl_jbs[[1]], data = jbs1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_jbs[[2]] <- attribute_in_network(net = dl_jbs[[2]], data = jbs2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_jbs[[3]] <- attribute_in_network(net = dl_jbs[[3]], data = jbs3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_jbs[[4]] <- attribute_in_network(net = dl_jbs[[4]], data = jbs4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_jbs[[5]] <- attribute_in_network(net = dl_jbs[[5]], data = jbs6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

#

dl_pm[[1]] <- attribute_in_network(net = dl_pm[[1]], data = pm1, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

dl_pm[[2]] <- attribute_in_network(net = dl_pm[[2]], data = pm2, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

dl_pm[[3]] <- attribute_in_network(net = dl_pm[[3]], data = pm3, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

dl_pm[[4]] <- attribute_in_network(net = dl_pm[[4]], data = pm4, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

dl_pm[[5]] <- attribute_in_network(net = dl_pm[[5]], data = pm6, id_data = c("id"),
                                   var_names = a,
                                   new_names = b)

#

dl_dmz[[1]] <- attribute_in_network(net = dl_dmz[[1]], data = dmz1, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_dmz[[2]] <- attribute_in_network(net = dl_dmz[[2]], data = dmz2, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_dmz[[3]] <- attribute_in_network(net = dl_dmz[[3]], data = dmz3, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_dmz[[4]] <- attribute_in_network(net = dl_dmz[[4]], data = dmz4, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

dl_dmz[[5]] <- attribute_in_network(net = dl_dmz[[5]], data = dmz6, id_data = c("id"),
                                    var_names = a,
                                    new_names = b)

###

top5_rvl[[1]] <- attribute_in_network(net = top5_rvl[[1]], data = rvl1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_rvl[[2]] <- attribute_in_network(net = top5_rvl[[2]], data = rvl2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_rvl[[3]] <- attribute_in_network(net = top5_rvl[[3]], data = rvl3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_rvl[[4]] <- attribute_in_network(net = top5_rvl[[4]], data = rvl4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_rvl[[5]] <- attribute_in_network(net = top5_rvl[[5]], data = rvl6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

top5_jbs[[1]] <- attribute_in_network(net = top5_jbs[[1]], data = jbs1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_jbs[[2]] <- attribute_in_network(net = top5_jbs[[2]], data = jbs2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_jbs[[3]] <- attribute_in_network(net = top5_jbs[[3]], data = jbs3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_jbs[[4]] <- attribute_in_network(net = top5_jbs[[4]], data = jbs4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_jbs[[5]] <- attribute_in_network(net = top5_jbs[[5]], data = jbs6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

top5_pm[[1]] <- attribute_in_network(net = top5_pm[[1]], data = pm1, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

top5_pm[[2]] <- attribute_in_network(net = top5_pm[[2]], data = pm2, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

top5_pm[[3]] <- attribute_in_network(net = top5_pm[[3]], data = pm3, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

top5_pm[[4]] <- attribute_in_network(net = top5_pm[[4]], data = pm4, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

top5_pm[[5]] <- attribute_in_network(net = top5_pm[[5]], data = pm6, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

#

top5_dmz[[1]] <- attribute_in_network(net = top5_dmz[[1]], data = dmz1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_dmz[[2]] <- attribute_in_network(net = top5_dmz[[2]], data = dmz2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_dmz[[3]] <- attribute_in_network(net = top5_dmz[[3]], data = dmz3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_dmz[[4]] <- attribute_in_network(net = top5_dmz[[4]], data = dmz4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

top5_dmz[[5]] <- attribute_in_network(net = top5_dmz[[5]], data = dmz6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

###

bully_rvl[[1]] <- attribute_in_network(net = bully_rvl[[1]], data = rvl1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_rvl[[2]] <- attribute_in_network(net = bully_rvl[[2]], data = rvl2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_rvl[[3]] <- attribute_in_network(net = bully_rvl[[3]], data = rvl3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_rvl[[4]] <- attribute_in_network(net = bully_rvl[[4]], data = rvl4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_rvl[[5]] <- attribute_in_network(net = bully_rvl[[5]], data = rvl6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

bully_jbs[[1]] <- attribute_in_network(net = bully_jbs[[1]], data = jbs1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_jbs[[2]] <- attribute_in_network(net = bully_jbs[[2]], data = jbs2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_jbs[[3]] <- attribute_in_network(net = bully_jbs[[3]], data = jbs3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_jbs[[4]] <- attribute_in_network(net = bully_jbs[[4]], data = jbs4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_jbs[[5]] <- attribute_in_network(net = bully_jbs[[5]], data = jbs6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

bully_pm[[1]] <- attribute_in_network(net = bully_pm[[1]], data = pm1, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

bully_pm[[2]] <- attribute_in_network(net = bully_pm[[2]], data = pm2, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

bully_pm[[3]] <- attribute_in_network(net = bully_pm[[3]], data = pm3, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

bully_pm[[4]] <- attribute_in_network(net = bully_pm[[4]], data = pm4, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

bully_pm[[5]] <- attribute_in_network(net = bully_pm[[5]], data = pm6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

bully_dmz[[1]] <- attribute_in_network(net = bully_dmz[[1]], data = dmz1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_dmz[[2]] <- attribute_in_network(net = bully_dmz[[2]], data = dmz2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_dmz[[3]] <- attribute_in_network(net = bully_dmz[[3]], data = dmz3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_dmz[[4]] <- attribute_in_network(net = bully_dmz[[4]], data = dmz4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

bully_dmz[[5]] <- attribute_in_network(net = bully_dmz[[5]], data = dmz6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

###

cool_rvl[[1]] <- attribute_in_network(net = cool_rvl[[1]], data = rvl1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_rvl[[2]] <- attribute_in_network(net = cool_rvl[[2]], data = rvl2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_rvl[[3]] <- attribute_in_network(net = cool_rvl[[3]], data = rvl3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_rvl[[4]] <- attribute_in_network(net = cool_rvl[[4]], data = rvl4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_rvl[[5]] <- attribute_in_network(net = cool_rvl[[5]], data = rvl6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

cool_jbs[[1]] <- attribute_in_network(net = cool_jbs[[1]], data = jbs1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_jbs[[2]] <- attribute_in_network(net = cool_jbs[[2]], data = jbs2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_jbs[[3]] <- attribute_in_network(net = cool_jbs[[3]], data = jbs3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_jbs[[4]] <- attribute_in_network(net = cool_jbs[[4]], data = jbs4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_jbs[[5]] <- attribute_in_network(net = cool_jbs[[5]], data = jbs6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

cool_pm[[1]] <- attribute_in_network(net = cool_pm[[1]], data = pm1, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

cool_pm[[2]] <- attribute_in_network(net = cool_pm[[2]], data = pm2, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

cool_pm[[3]] <- attribute_in_network(net = cool_pm[[3]], data = pm3, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

cool_pm[[4]] <- attribute_in_network(net = cool_pm[[4]], data = pm4, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

cool_pm[[5]] <- attribute_in_network(net = cool_pm[[5]], data = pm6, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

#

cool_dmz[[1]] <- attribute_in_network(net = cool_dmz[[1]], data = dmz1, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_dmz[[2]] <- attribute_in_network(net = cool_dmz[[2]], data = dmz2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_dmz[[3]] <- attribute_in_network(net = cool_dmz[[3]], data = dmz3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_dmz[[4]] <- attribute_in_network(net = cool_dmz[[4]], data = dmz4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

cool_dmz[[5]] <- attribute_in_network(net = cool_dmz[[5]], data = dmz6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

###

known_rvl[[2]] <- attribute_in_network(net = known_rvl[[2]], data = rvl2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_rvl[[3]] <- attribute_in_network(net = known_rvl[[3]], data = rvl3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_rvl[[4]] <- attribute_in_network(net = known_rvl[[4]], data = rvl4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_rvl[[5]] <- attribute_in_network(net = known_rvl[[5]], data = rvl6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

known_jbs[[2]] <- attribute_in_network(net = known_jbs[[2]], data = jbs2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_jbs[[3]] <- attribute_in_network(net = known_jbs[[3]], data = jbs3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_jbs[[4]] <- attribute_in_network(net = known_jbs[[4]], data = jbs4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_jbs[[5]] <- attribute_in_network(net = known_jbs[[5]], data = jbs6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

known_pm[[2]] <- attribute_in_network(net = known_pm[[2]], data = pm2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

known_pm[[3]] <- attribute_in_network(net = known_pm[[3]], data = pm3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

known_pm[[4]] <- attribute_in_network(net = known_pm[[4]], data = pm4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

known_pm[[5]] <- attribute_in_network(net = known_pm[[5]], data = pm6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

known_dmz[[2]] <- attribute_in_network(net = known_dmz[[2]], data = dmz2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_dmz[[3]] <- attribute_in_network(net = known_dmz[[3]], data = dmz3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_dmz[[4]] <- attribute_in_network(net = known_dmz[[4]], data = dmz4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

known_dmz[[5]] <- attribute_in_network(net = known_dmz[[5]], data = dmz6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

###

fear_rvl[[2]] <- attribute_in_network(net = fear_rvl[[2]], data = rvl2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_rvl[[3]] <- attribute_in_network(net = fear_rvl[[3]], data = rvl3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_rvl[[4]] <- attribute_in_network(net = fear_rvl[[4]], data = rvl4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_rvl[[5]] <- attribute_in_network(net = fear_rvl[[5]], data = rvl6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

fear_jbs[[2]] <- attribute_in_network(net = fear_jbs[[2]], data = jbs2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_jbs[[3]] <- attribute_in_network(net = fear_jbs[[3]], data = jbs3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_jbs[[4]] <- attribute_in_network(net = fear_jbs[[4]], data = jbs4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_jbs[[5]] <- attribute_in_network(net = fear_jbs[[5]], data = jbs6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

fear_pm[[2]] <- attribute_in_network(net = fear_pm[[2]], data = pm2, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

fear_pm[[3]] <- attribute_in_network(net = fear_pm[[3]], data = pm3, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

fear_pm[[4]] <- attribute_in_network(net = fear_pm[[4]], data = pm4, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

fear_pm[[5]] <- attribute_in_network(net = fear_pm[[5]], data = pm6, id_data = c("id"),
                                     var_names = a,
                                     new_names = b)

#

fear_dmz[[2]] <- attribute_in_network(net = fear_dmz[[2]], data = dmz2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_dmz[[3]] <- attribute_in_network(net = fear_dmz[[3]], data = dmz3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_dmz[[4]] <- attribute_in_network(net = fear_dmz[[4]], data = dmz4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

fear_dmz[[5]] <- attribute_in_network(net = fear_dmz[[5]], data = dmz6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

###

liked_rvl[[2]] <- attribute_in_network(net = liked_rvl[[2]], data = rvl2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_rvl[[3]] <- attribute_in_network(net = liked_rvl[[3]], data = rvl3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_rvl[[4]] <- attribute_in_network(net = liked_rvl[[4]], data = rvl4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_rvl[[5]] <- attribute_in_network(net = liked_rvl[[5]], data = rvl6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

liked_jbs[[2]] <- attribute_in_network(net = liked_jbs[[2]], data = jbs2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_jbs[[3]] <- attribute_in_network(net = liked_jbs[[3]], data = jbs3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_jbs[[4]] <- attribute_in_network(net = liked_jbs[[4]], data = jbs4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_jbs[[5]] <- attribute_in_network(net = liked_jbs[[5]], data = jbs6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

#

liked_pm[[2]] <- attribute_in_network(net = liked_pm[[2]], data = pm2, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

liked_pm[[3]] <- attribute_in_network(net = liked_pm[[3]], data = pm3, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

liked_pm[[4]] <- attribute_in_network(net = liked_pm[[4]], data = pm4, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

liked_pm[[5]] <- attribute_in_network(net = liked_pm[[5]], data = pm6, id_data = c("id"),
                                      var_names = a,
                                      new_names = b)

#

liked_dmz[[2]] <- attribute_in_network(net = liked_dmz[[2]], data = dmz2, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_dmz[[3]] <- attribute_in_network(net = liked_dmz[[3]], data = dmz3, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_dmz[[4]] <- attribute_in_network(net = liked_dmz[[4]], data = dmz4, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

liked_dmz[[5]] <- attribute_in_network(net = liked_dmz[[5]], data = dmz6, id_data = c("id"),
                                       var_names = a,
                                       new_names = b)

##### e. "Temporal" information - transfer from culture #####

names(culture)

a1 <- c("art_lvl_w1","music_lvl_w1","book_lvl_w1","cult_lvl_w1","bin_music_w1","bin_sport_w1","bin_book_w1")
a3 <- c("art_lvl_w3","music_lvl_w3","book_lvl_w3","cult_lvl_w3","bin_music_w3","bin_sport_w3","bin_book_w3")

a2 <- c("m_fact1_w2","m_fact2_w2","m_fact3_w2","y_fact1_w2","y_fact2_w2","ym_fact1_w2","ym_fact2_w2","tree_w2","forest_w2")
a4 <- c("m_fact1_w4","m_fact2_w4","m_fact3_w4","y_fact1_w4","y_fact2_w4","ym_fact1_w4","ym_fact2_w4","tree_w4","forest_w4")

b1 <- c("art_lvl","music_lvl","book_lvl","cult_lvl","bin_music","bin_sport","bin_book")
b3 <- c("art_lvl","music_lvl","book_lvl","cult_lvl","bin_music","bin_sport","bin_book")

b2 <- c("m_fact1","m_fact2","m_fact3","y_fact1","y_fact2","ym_fact1","ym_fact2","tree","forest")
b4 <- c("m_fact1","m_fact2","m_fact3","y_fact1","y_fact2","ym_fact1","ym_fact2","tree","forest")


# NB : the cultural variables from w1 are also placed in the w2 networks, and those of w3 in the w4 networks.
# Similarly, w2 variables are put in the w1 networks, and w4 in w3 (yeah, it goes backward, but otherwise I wouldn t get
# any info in w1...)
# Finally, in wave 6, the infos are the same as in wave 4... Should be changed one day, when I code the cultural category
# of wave 6 (but for now I don t have them).

gf_rvl[[1]] <- attribute_in_network(net = gf_rvl[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_rvl[[2]] <- attribute_in_network(net = gf_rvl[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_rvl[[3]] <- attribute_in_network(net = gf_rvl[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_rvl[[4]] <- attribute_in_network(net = gf_rvl[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_rvl[[5]] <- attribute_in_network(net = gf_rvl[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

#

gf_jbs[[1]] <- attribute_in_network(net = gf_jbs[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_jbs[[2]] <- attribute_in_network(net = gf_jbs[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_jbs[[3]] <- attribute_in_network(net = gf_jbs[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_jbs[[4]] <- attribute_in_network(net = gf_jbs[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_jbs[[5]] <- attribute_in_network(net = gf_jbs[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

#

gf_pm[[1]] <- attribute_in_network(net = gf_pm[[1]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

gf_pm[[2]] <- attribute_in_network(net = gf_pm[[2]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

gf_pm[[3]] <- attribute_in_network(net = gf_pm[[3]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

gf_pm[[4]] <- attribute_in_network(net = gf_pm[[4]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

gf_pm[[5]] <- attribute_in_network(net = gf_pm[[5]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

#

gf_dmz[[1]] <- attribute_in_network(net = gf_dmz[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_dmz[[2]] <- attribute_in_network(net = gf_dmz[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

gf_dmz[[3]] <- attribute_in_network(net = gf_dmz[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_dmz[[4]] <- attribute_in_network(net = gf_dmz[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

gf_dmz[[5]] <- attribute_in_network(net = gf_dmz[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

###

f_rvl[[1]] <- attribute_in_network(net = f_rvl[[1]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_rvl[[2]] <- attribute_in_network(net = f_rvl[[2]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_rvl[[3]] <- attribute_in_network(net = f_rvl[[3]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_rvl[[4]] <- attribute_in_network(net = f_rvl[[4]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_rvl[[5]] <- attribute_in_network(net = f_rvl[[5]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

#

f_jbs[[1]] <- attribute_in_network(net = f_jbs[[1]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_jbs[[2]] <- attribute_in_network(net = f_jbs[[2]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_jbs[[3]] <- attribute_in_network(net = f_jbs[[3]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_jbs[[4]] <- attribute_in_network(net = f_jbs[[4]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_jbs[[5]] <- attribute_in_network(net = f_jbs[[5]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

#

f_pm[[1]] <- attribute_in_network(net = f_pm[[1]], data = culture, id_data = c("id"),
                                  var_names = c(a1,a2),
                                  new_names = c(b1,b2))

f_pm[[2]] <- attribute_in_network(net = f_pm[[2]], data = culture, id_data = c("id"),
                                  var_names = c(a1,a2),
                                  new_names = c(b1,b2))

f_pm[[3]] <- attribute_in_network(net = f_pm[[3]], data = culture, id_data = c("id"),
                                  var_names = c(a3,a4),
                                  new_names = c(b3,b4))

f_pm[[4]] <- attribute_in_network(net = f_pm[[4]], data = culture, id_data = c("id"),
                                  var_names = c(a3,a4),
                                  new_names = c(b3,b4))

f_pm[[5]] <- attribute_in_network(net = f_pm[[5]], data = culture, id_data = c("id"),
                                  var_names = c(a3,a4),
                                  new_names = c(b3,b4))

#

f_dmz[[1]] <- attribute_in_network(net = f_dmz[[1]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_dmz[[2]] <- attribute_in_network(net = f_dmz[[2]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

f_dmz[[3]] <- attribute_in_network(net = f_dmz[[3]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_dmz[[4]] <- attribute_in_network(net = f_dmz[[4]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

f_dmz[[5]] <- attribute_in_network(net = f_dmz[[5]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

###

allf_rvl[[1]] <- attribute_in_network(net = allf_rvl[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_rvl[[2]] <- attribute_in_network(net = allf_rvl[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_rvl[[3]] <- attribute_in_network(net = allf_rvl[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_rvl[[4]] <- attribute_in_network(net = allf_rvl[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_rvl[[5]] <- attribute_in_network(net = allf_rvl[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

allf_jbs[[1]] <- attribute_in_network(net = allf_jbs[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_jbs[[2]] <- attribute_in_network(net = allf_jbs[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_jbs[[3]] <- attribute_in_network(net = allf_jbs[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_jbs[[4]] <- attribute_in_network(net = allf_jbs[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_jbs[[5]] <- attribute_in_network(net = allf_jbs[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))


#

allf_pm[[1]] <- attribute_in_network(net = allf_pm[[1]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

allf_pm[[2]] <- attribute_in_network(net = allf_pm[[2]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

allf_pm[[3]] <- attribute_in_network(net = allf_pm[[3]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

allf_pm[[4]] <- attribute_in_network(net = allf_pm[[4]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

allf_pm[[5]] <- attribute_in_network(net = allf_pm[[5]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

#

allf_dmz[[1]] <- attribute_in_network(net = allf_dmz[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_dmz[[2]] <- attribute_in_network(net = allf_dmz[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

allf_dmz[[3]] <- attribute_in_network(net = allf_dmz[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_dmz[[4]] <- attribute_in_network(net = allf_dmz[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

allf_dmz[[5]] <- attribute_in_network(net = allf_dmz[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

###

dl_rvl[[1]] <- attribute_in_network(net = dl_rvl[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_rvl[[2]] <- attribute_in_network(net = dl_rvl[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_rvl[[3]] <- attribute_in_network(net = dl_rvl[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_rvl[[4]] <- attribute_in_network(net = dl_rvl[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_rvl[[5]] <- attribute_in_network(net = dl_rvl[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

#

dl_jbs[[1]] <- attribute_in_network(net = dl_jbs[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_jbs[[2]] <- attribute_in_network(net = dl_jbs[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_jbs[[3]] <- attribute_in_network(net = dl_jbs[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_jbs[[4]] <- attribute_in_network(net = dl_jbs[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_jbs[[5]] <- attribute_in_network(net = dl_jbs[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

#

dl_pm[[1]] <- attribute_in_network(net = dl_pm[[1]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

dl_pm[[2]] <- attribute_in_network(net = dl_pm[[2]], data = culture, id_data = c("id"),
                                   var_names = c(a1,a2),
                                   new_names = c(b1,b2))

dl_pm[[3]] <- attribute_in_network(net = dl_pm[[3]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

dl_pm[[4]] <- attribute_in_network(net = dl_pm[[4]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

dl_pm[[5]] <- attribute_in_network(net = dl_pm[[5]], data = culture, id_data = c("id"),
                                   var_names = c(a3,a4),
                                   new_names = c(b3,b4))

#

dl_dmz[[1]] <- attribute_in_network(net = dl_dmz[[1]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_dmz[[2]] <- attribute_in_network(net = dl_dmz[[2]], data = culture, id_data = c("id"),
                                    var_names = c(a1,a2),
                                    new_names = c(b1,b2))

dl_dmz[[3]] <- attribute_in_network(net = dl_dmz[[3]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_dmz[[4]] <- attribute_in_network(net = dl_dmz[[4]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

dl_dmz[[5]] <- attribute_in_network(net = dl_dmz[[5]], data = culture, id_data = c("id"),
                                    var_names = c(a3,a4),
                                    new_names = c(b3,b4))

###

top5_rvl[[1]] <- attribute_in_network(net = top5_rvl[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_rvl[[2]] <- attribute_in_network(net = top5_rvl[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_rvl[[3]] <- attribute_in_network(net = top5_rvl[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

top5_rvl[[4]] <- attribute_in_network(net = top5_rvl[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

top5_rvl[[5]] <- attribute_in_network(net = top5_rvl[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

top5_jbs[[1]] <- attribute_in_network(net = top5_jbs[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_jbs[[2]] <- attribute_in_network(net = top5_jbs[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_jbs[[3]] <- attribute_in_network(net = top5_jbs[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

top5_jbs[[4]] <- attribute_in_network(net = top5_jbs[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))


top5_jbs[[5]] <- attribute_in_network(net = top5_jbs[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

top5_pm[[1]] <- attribute_in_network(net = top5_pm[[1]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

top5_pm[[2]] <- attribute_in_network(net = top5_pm[[2]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

top5_pm[[3]] <- attribute_in_network(net = top5_pm[[3]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

top5_pm[[4]] <- attribute_in_network(net = top5_pm[[4]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

top5_pm[[5]] <- attribute_in_network(net = top5_pm[[5]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

#

top5_dmz[[1]] <- attribute_in_network(net = top5_dmz[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_dmz[[2]] <- attribute_in_network(net = top5_dmz[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

top5_dmz[[3]] <- attribute_in_network(net = top5_dmz[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

top5_dmz[[4]] <- attribute_in_network(net = top5_dmz[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

top5_dmz[[5]] <- attribute_in_network(net = top5_dmz[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

###

bully_rvl[[1]] <- attribute_in_network(net = bully_rvl[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_rvl[[2]] <- attribute_in_network(net = bully_rvl[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_rvl[[3]] <- attribute_in_network(net = bully_rvl[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_rvl[[4]] <- attribute_in_network(net = bully_rvl[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_rvl[[5]] <- attribute_in_network(net = bully_rvl[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

bully_jbs[[1]] <- attribute_in_network(net = bully_jbs[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_jbs[[2]] <- attribute_in_network(net = bully_jbs[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_jbs[[3]] <- attribute_in_network(net = bully_jbs[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_jbs[[4]] <- attribute_in_network(net = bully_jbs[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_jbs[[5]] <- attribute_in_network(net = bully_jbs[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

bully_pm[[1]] <- attribute_in_network(net = bully_pm[[1]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

bully_pm[[2]] <- attribute_in_network(net = bully_pm[[2]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

bully_pm[[3]] <- attribute_in_network(net = bully_pm[[3]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

bully_pm[[4]] <- attribute_in_network(net = bully_pm[[4]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

bully_pm[[5]] <- attribute_in_network(net = bully_pm[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

bully_dmz[[1]] <- attribute_in_network(net = bully_dmz[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_dmz[[2]] <- attribute_in_network(net = bully_dmz[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

bully_dmz[[3]] <- attribute_in_network(net = bully_dmz[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_dmz[[4]] <- attribute_in_network(net = bully_dmz[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

bully_dmz[[5]] <- attribute_in_network(net = bully_dmz[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

###

cool_rvl[[1]] <- attribute_in_network(net = cool_rvl[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_rvl[[2]] <- attribute_in_network(net = cool_rvl[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_rvl[[3]] <- attribute_in_network(net = cool_rvl[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_rvl[[4]] <- attribute_in_network(net = cool_rvl[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_rvl[[5]] <- attribute_in_network(net = cool_rvl[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

cool_jbs[[1]] <- attribute_in_network(net = cool_jbs[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_jbs[[2]] <- attribute_in_network(net = cool_jbs[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_jbs[[3]] <- attribute_in_network(net = cool_jbs[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_jbs[[4]] <- attribute_in_network(net = cool_jbs[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_jbs[[5]] <- attribute_in_network(net = cool_jbs[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

cool_pm[[1]] <- attribute_in_network(net = cool_pm[[1]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

cool_pm[[2]] <- attribute_in_network(net = cool_pm[[2]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

cool_pm[[3]] <- attribute_in_network(net = cool_pm[[3]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

cool_pm[[4]] <- attribute_in_network(net = cool_pm[[4]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

cool_pm[[5]] <- attribute_in_network(net = cool_pm[[5]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

#

cool_dmz[[1]] <- attribute_in_network(net = cool_dmz[[1]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_dmz[[2]] <- attribute_in_network(net = cool_dmz[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

cool_dmz[[3]] <- attribute_in_network(net = cool_dmz[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_dmz[[4]] <- attribute_in_network(net = cool_dmz[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

cool_dmz[[5]] <- attribute_in_network(net = cool_dmz[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

###

known_rvl[[2]] <- attribute_in_network(net = known_rvl[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

known_rvl[[3]] <- attribute_in_network(net = known_rvl[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_rvl[[4]] <- attribute_in_network(net = known_rvl[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_rvl[[5]] <- attribute_in_network(net = known_rvl[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

known_jbs[[2]] <- attribute_in_network(net = known_jbs[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

known_jbs[[3]] <- attribute_in_network(net = known_jbs[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_jbs[[4]] <- attribute_in_network(net = known_jbs[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_jbs[[5]] <- attribute_in_network(net = known_jbs[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

known_pm[[2]] <- attribute_in_network(net = known_pm[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

known_pm[[3]] <- attribute_in_network(net = known_pm[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

known_pm[[4]] <- attribute_in_network(net = known_pm[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

known_pm[[5]] <- attribute_in_network(net = known_pm[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

known_dmz[[2]] <- attribute_in_network(net = known_dmz[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

known_dmz[[3]] <- attribute_in_network(net = known_dmz[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_dmz[[4]] <- attribute_in_network(net = known_dmz[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

known_dmz[[5]] <- attribute_in_network(net = known_dmz[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

###

fear_rvl[[2]] <- attribute_in_network(net = fear_rvl[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

fear_rvl[[3]] <- attribute_in_network(net = fear_rvl[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_rvl[[4]] <- attribute_in_network(net = fear_rvl[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_rvl[[5]] <- attribute_in_network(net = fear_rvl[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

fear_jbs[[2]] <- attribute_in_network(net = fear_jbs[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

fear_jbs[[3]] <- attribute_in_network(net = fear_jbs[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_jbs[[4]] <- attribute_in_network(net = fear_jbs[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_jbs[[5]] <- attribute_in_network(net = fear_jbs[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

fear_pm[[2]] <- attribute_in_network(net = fear_pm[[2]], data = culture, id_data = c("id"),
                                     var_names = c(a1,a2),
                                     new_names = c(b1,b2))

fear_pm[[3]] <- attribute_in_network(net = fear_pm[[3]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

fear_pm[[4]] <- attribute_in_network(net = fear_pm[[4]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

fear_pm[[5]] <- attribute_in_network(net = fear_pm[[5]], data = culture, id_data = c("id"),
                                     var_names = c(a3,a4),
                                     new_names = c(b3,b4))

#

fear_dmz[[2]] <- attribute_in_network(net = fear_dmz[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

fear_dmz[[3]] <- attribute_in_network(net = fear_dmz[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_dmz[[4]] <- attribute_in_network(net = fear_dmz[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

fear_dmz[[5]] <- attribute_in_network(net = fear_dmz[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

###

liked_rvl[[2]] <- attribute_in_network(net = liked_rvl[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

liked_rvl[[3]] <- attribute_in_network(net = liked_rvl[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_rvl[[4]] <- attribute_in_network(net = liked_rvl[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_rvl[[5]] <- attribute_in_network(net = liked_rvl[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

liked_jbs[[2]] <- attribute_in_network(net = liked_jbs[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

liked_jbs[[3]] <- attribute_in_network(net = liked_jbs[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_jbs[[4]] <- attribute_in_network(net = liked_jbs[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_jbs[[5]] <- attribute_in_network(net = liked_jbs[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

#

liked_pm[[2]] <- attribute_in_network(net = liked_pm[[2]], data = culture, id_data = c("id"),
                                      var_names = c(a1,a2),
                                      new_names = c(b1,b2))

liked_pm[[3]] <- attribute_in_network(net = liked_pm[[3]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

liked_pm[[4]] <- attribute_in_network(net = liked_pm[[4]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

liked_pm[[5]] <- attribute_in_network(net = liked_pm[[5]], data = culture, id_data = c("id"),
                                      var_names = c(a3,a4),
                                      new_names = c(b3,b4))

#

liked_dmz[[2]] <- attribute_in_network(net = liked_dmz[[2]], data = culture, id_data = c("id"),
                                       var_names = c(a1,a2),
                                       new_names = c(b1,b2))

liked_dmz[[3]] <- attribute_in_network(net = liked_dmz[[3]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_dmz[[4]] <- attribute_in_network(net = liked_dmz[[4]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

liked_dmz[[5]] <- attribute_in_network(net = liked_dmz[[5]], data = culture, id_data = c("id"),
                                       var_names = c(a3,a4),
                                       new_names = c(b3,b4))

##### f. "Temporal" information - transfer from psycho #####

# UPDATE august 2022 : I now have another version of the factor scores, called "xx1" (instead of x1),
# which uses a factorial space built with data from w2, w4 and w6 (instead of just w2 and w4).
# Same for myneruo and myextrav, replaced by "new_myneuro" and "new_myextrav".
# Results are very close for w2 and w4. For w6, of course, only xx is available, not x.

# Choose the option here for using either x or xx. xx is more up to date, but all the treatments I did
# up to august 2022 (including the PhD) were based on x, so it may be useful to load the x variables instead
# for the reproductibility of these treatments.

use_psycho <- "new" # either "old" (x) or "new" (xx and attribution in w6).

names(psycho)

if(use_psycho=="old"){
  # Some students were here only for a given wave, in which case we want to use the info on that wave
  # for all waves.
  # Meaning when a wave 2 variable is in NA, we attribute the w4 value, and vice versa.
  
  psycho$x1_2[is.na(psycho$x1_2)==T & is.na(psycho$x1_4)==F] <- psycho$x1_4[is.na(psycho$x1_2)==T & is.na(psycho$x1_4)==F]
  psycho$x2_2[is.na(psycho$x2_2)==T & is.na(psycho$x2_4)==F] <- psycho$x2_4[is.na(psycho$x2_2)==T & is.na(psycho$x2_4)==F]
  psycho$x3_2[is.na(psycho$x3_2)==T & is.na(psycho$x3_4)==F] <- psycho$x3_4[is.na(psycho$x3_2)==T & is.na(psycho$x3_4)==F]
  
  psycho$x1_4[is.na(psycho$x1_2)==F & is.na(psycho$x1_4)==T] <- psycho$x1_2[is.na(psycho$x1_2)==F & is.na(psycho$x1_4)==T]
  psycho$x2_4[is.na(psycho$x2_2)==F & is.na(psycho$x2_4)==T] <- psycho$x2_2[is.na(psycho$x2_2)==F & is.na(psycho$x2_4)==T]
  psycho$x3_4[is.na(psycho$x3_2)==F & is.na(psycho$x3_4)==T] <- psycho$x3_2[is.na(psycho$x3_2)==F & is.na(psycho$x3_4)==T]
  
  psycho$myextrav_2[is.na(psycho$myextrav_2)==T & is.na(psycho$myextrav_4)==F] <- psycho$myextrav_4[is.na(psycho$myextrav_2)==T & is.na(psycho$myextrav_4)==F]
  psycho$myneuro_2[is.na(psycho$myneuro_2)==T & is.na(psycho$myneuro_4)==F] <- psycho$myneuro_4[is.na(psycho$myneuro_2)==T & is.na(psycho$myneuro_4)==F]
  
  psycho$myextrav_4[is.na(psycho$myextrav_2)==F & is.na(psycho$myextrav_4)==T] <- psycho$myextrav_2[is.na(psycho$myextrav_2)==F & is.na(psycho$myextrav_4)==T]
  psycho$myneuro_4[is.na(psycho$myneuro_2)==F & is.na(psycho$myneuro_4)==T] <- psycho$myneuro_2[is.na(psycho$myneuro_2)==F & is.na(psycho$myneuro_4)==T]
  
  a2 <-c("x1_2","x2_2","x3_2","extrav_w2","agreeab_w2","consc_w2","neuro_w2","open_w2","myextrav_2","myneuro_2")
  a4 <-c("x1_4","x2_4","x3_4","extrav_w4","agreeab_w4","consc_w4","neuro_w4","open_w4","myextrav_4","myneuro_4")
  
  b2 <-c("psy_fact1","psy_fact2","psy_fact3","extrav","agreeab","consc","neuro","open","myextrav","myneuro")
  b4 <-c("psy_fact1","psy_fact2","psy_fact3","extrav","agreeab","consc","neuro","open","myextrav","myneuro")
}

if(use_psycho=="new"){
  # Some students were here only for a given wave, in which case we want to use the info on that wave
  # for all waves.
  # Meaning when a wave 2 variable is in NA, we attribute the w4 value, and vice versa.
  # As for w6, we use data from wave 4.
  
  psycho$xx1_2[is.na(psycho$xx1_2)==T & is.na(psycho$xx1_4)==F] <- psycho$xx1_4[is.na(psycho$xx1_2)==T & is.na(psycho$xx1_4)==F]
  psycho$xx2_2[is.na(psycho$xx2_2)==T & is.na(psycho$xx2_4)==F] <- psycho$xx2_4[is.na(psycho$xx2_2)==T & is.na(psycho$xx2_4)==F]
  psycho$xx3_2[is.na(psycho$xx3_2)==T & is.na(psycho$xx3_4)==F] <- psycho$xx3_4[is.na(psycho$xx3_2)==T & is.na(psycho$xx3_4)==F]
  
  psycho$xx1_4[is.na(psycho$xx1_2)==F & is.na(psycho$xx1_4)==T] <- psycho$xx1_2[is.na(psycho$xx1_2)==F & is.na(psycho$xx1_4)==T]
  psycho$xx2_4[is.na(psycho$xx2_2)==F & is.na(psycho$xx2_4)==T] <- psycho$xx2_2[is.na(psycho$xx2_2)==F & is.na(psycho$xx2_4)==T]
  psycho$xx3_4[is.na(psycho$xx3_2)==F & is.na(psycho$xx3_4)==T] <- psycho$xx3_2[is.na(psycho$xx3_2)==F & is.na(psycho$xx3_4)==T]
  
  psycho$xx1_6[is.na(psycho$xx1_6)==T & is.na(psycho$xx1_4)==F] <- psycho$xx1_4[is.na(psycho$xx1_6)==T & is.na(psycho$xx1_4)==F]
  psycho$xx2_6[is.na(psycho$xx2_6)==T & is.na(psycho$xx2_4)==F] <- psycho$xx2_4[is.na(psycho$xx2_6)==T & is.na(psycho$xx2_4)==F]
  psycho$xx3_6[is.na(psycho$xx3_6)==T & is.na(psycho$xx3_4)==F] <- psycho$xx3_4[is.na(psycho$xx3_6)==T & is.na(psycho$xx3_4)==F]
  
  psycho$new_myextrav_2[is.na(psycho$new_myextrav_2)==T & is.na(psycho$new_myextrav_4)==F] <- psycho$new_myextrav_4[is.na(psycho$new_myextrav_2)==T & is.na(psycho$new_myextrav_4)==F]
  psycho$new_myneuro_2[is.na(psycho$new_myneuro_2)==T & is.na(psycho$new_myneuro_4)==F] <- psycho$new_myneuro_4[is.na(psycho$new_myneuro_2)==T & is.na(psycho$new_myneuro_4)==F]
  
  psycho$new_myextrav_4[is.na(psycho$new_myextrav_2)==F & is.na(psycho$new_myextrav_4)==T] <- psycho$new_myextrav_2[is.na(psycho$new_myextrav_2)==F & is.na(psycho$new_myextrav_4)==T]
  psycho$new_myneuro_4[is.na(psycho$new_myneuro_2)==F & is.na(psycho$new_myneuro_4)==T] <- psycho$new_myneuro_2[is.na(psycho$new_myneuro_2)==F & is.na(psycho$new_myneuro_4)==T]
  
  psycho$new_myextrav_6[is.na(psycho$new_myextrav_6)==T & is.na(psycho$new_myextrav_4)==F] <- psycho$new_myextrav_4[is.na(psycho$new_myextrav_6)==T & is.na(psycho$new_myextrav_4)==F]
  psycho$new_myneuro_6[is.na(psycho$new_myneuro_6)==T & is.na(psycho$new_myneuro_4)==F] <- psycho$new_myneuro_4[is.na(psycho$new_myneuro_6)==T & is.na(psycho$new_myneuro_4)==F]
  
  a2 <-c("xx1_2","xx2_2","xx3_2","extrav_w2","agreeab_w2","consc_w2","neuro_w2","open_w2","new_myextrav_2","new_myneuro_2")
  a4 <-c("xx1_4","xx2_4","xx3_4","extrav_w4","agreeab_w4","consc_w4","neuro_w4","open_w4","new_myextrav_4","new_myneuro_4")
  a6 <-c("xx1_6","xx2_6","xx3_6","extrav_w6","agreeab_w6","consc_w6","neuro_w6","open_w6","new_myextrav_6","new_myneuro_6")
  
  b2 <-c("psy_fact1","psy_fact2","psy_fact3","extrav","agreeab","consc","neuro","open","new_myextrav","new_myneuro")
  b4 <-c("psy_fact1","psy_fact2","psy_fact3","extrav","agreeab","consc","neuro","open","new_myextrav","new_myneuro")
  b6 <-c("psy_fact1","psy_fact2","psy_fact3","extrav","agreeab","consc","neuro","open","new_myextrav","new_myneuro")
}

# NB : the variables from w2 are also placed in the w1 networks, and those of w4 in the w3 networks.
# (yeah, it goes backward, but otherwise I wouldn t get any info in w1...)

gf_rvl[[1]] <- attribute_in_network(net = gf_rvl[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_rvl[[2]] <- attribute_in_network(net = gf_rvl[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_rvl[[3]] <- attribute_in_network(net = gf_rvl[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

gf_rvl[[4]] <- attribute_in_network(net = gf_rvl[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  gf_rvl[[5]] <- attribute_in_network(net = gf_rvl[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

#

gf_jbs[[1]] <- attribute_in_network(net = gf_jbs[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_jbs[[2]] <- attribute_in_network(net = gf_jbs[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_jbs[[3]] <- attribute_in_network(net = gf_jbs[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

gf_jbs[[4]] <- attribute_in_network(net = gf_jbs[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  gf_jbs[[5]] <- attribute_in_network(net = gf_jbs[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

#

gf_pm[[1]] <- attribute_in_network(net = gf_pm[[1]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

gf_pm[[2]] <- attribute_in_network(net = gf_pm[[2]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

gf_pm[[3]] <- attribute_in_network(net = gf_pm[[3]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

gf_pm[[4]] <- attribute_in_network(net = gf_pm[[4]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

if(use_psycho=="new"){
  gf_pm[[5]] <- attribute_in_network(net = gf_pm[[5]], data = psycho, id_data = c("id"),
                                     var_names = a6,
                                     new_names = b6)
}

#

gf_dmz[[1]] <- attribute_in_network(net = gf_dmz[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_dmz[[2]] <- attribute_in_network(net = gf_dmz[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

gf_dmz[[3]] <- attribute_in_network(net = gf_dmz[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

gf_dmz[[4]] <- attribute_in_network(net = gf_dmz[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  gf_dmz[[5]] <- attribute_in_network(net = gf_dmz[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

###

f_rvl[[1]] <- attribute_in_network(net = f_rvl[[1]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_rvl[[2]] <- attribute_in_network(net = f_rvl[[2]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_rvl[[3]] <- attribute_in_network(net = f_rvl[[3]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

f_rvl[[4]] <- attribute_in_network(net = f_rvl[[4]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

if(use_psycho=="new"){
  f_rvl[[5]] <- attribute_in_network(net = f_rvl[[5]], data = psycho, id_data = c("id"),
                                     var_names = a6,
                                     new_names = b6)
}

#

f_jbs[[1]] <- attribute_in_network(net = f_jbs[[1]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_jbs[[2]] <- attribute_in_network(net = f_jbs[[2]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_jbs[[3]] <- attribute_in_network(net = f_jbs[[3]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

f_jbs[[4]] <- attribute_in_network(net = f_jbs[[4]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

if(use_psycho=="new"){
  f_jbs[[5]] <- attribute_in_network(net = f_jbs[[5]], data = psycho, id_data = c("id"),
                                     var_names = a6,
                                     new_names = b6)
}

#

f_pm[[1]] <- attribute_in_network(net = f_pm[[1]], data = psycho, id_data = c("id"),
                                  var_names = a2,
                                  new_names = b2)

f_pm[[2]] <- attribute_in_network(net = f_pm[[2]], data = psycho, id_data = c("id"),
                                  var_names = a2,
                                  new_names = b2)

f_pm[[3]] <- attribute_in_network(net = f_pm[[3]], data = psycho, id_data = c("id"),
                                  var_names = a4,
                                  new_names = b4)

f_pm[[4]] <- attribute_in_network(net = f_pm[[4]], data = psycho, id_data = c("id"),
                                  var_names = a4,
                                  new_names = b4)

if(use_psycho=="new"){
  f_pm[[5]] <- attribute_in_network(net = f_pm[[5]], data = psycho, id_data = c("id"),
                                    var_names = a6,
                                    new_names = b6)
}

#

f_dmz[[1]] <- attribute_in_network(net = f_dmz[[1]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_dmz[[2]] <- attribute_in_network(net = f_dmz[[2]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

f_dmz[[3]] <- attribute_in_network(net = f_dmz[[3]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

f_dmz[[4]] <- attribute_in_network(net = f_dmz[[4]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

if(use_psycho=="new"){
  f_dmz[[5]] <- attribute_in_network(net = f_dmz[[5]], data = psycho, id_data = c("id"),
                                     var_names = a6,
                                     new_names = b6)
}

###

allf_rvl[[1]] <- attribute_in_network(net = allf_rvl[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_rvl[[2]] <- attribute_in_network(net = allf_rvl[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_rvl[[3]] <- attribute_in_network(net = allf_rvl[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

allf_rvl[[4]] <- attribute_in_network(net = allf_rvl[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  allf_rvl[[5]] <- attribute_in_network(net = allf_rvl[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

allf_jbs[[1]] <- attribute_in_network(net = allf_jbs[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_jbs[[2]] <- attribute_in_network(net = allf_jbs[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_jbs[[3]] <- attribute_in_network(net = allf_jbs[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

allf_jbs[[4]] <- attribute_in_network(net = allf_jbs[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  allf_jbs[[5]] <- attribute_in_network(net = allf_jbs[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

allf_pm[[1]] <- attribute_in_network(net = allf_pm[[1]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

allf_pm[[2]] <- attribute_in_network(net = allf_pm[[2]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

allf_pm[[3]] <- attribute_in_network(net = allf_pm[[3]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

allf_pm[[4]] <- attribute_in_network(net = allf_pm[[4]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

if(use_psycho=="new"){
  allf_pm[[5]] <- attribute_in_network(net = allf_pm[[5]], data = psycho, id_data = c("id"),
                                       var_names = a6,
                                       new_names = b6)
}

#

allf_dmz[[1]] <- attribute_in_network(net = allf_dmz[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_dmz[[2]] <- attribute_in_network(net = allf_dmz[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

allf_dmz[[3]] <- attribute_in_network(net = allf_dmz[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

allf_dmz[[4]] <- attribute_in_network(net = allf_dmz[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  allf_dmz[[5]] <- attribute_in_network(net = allf_dmz[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

###

dl_rvl[[1]] <- attribute_in_network(net = dl_rvl[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_rvl[[2]] <- attribute_in_network(net = dl_rvl[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_rvl[[3]] <- attribute_in_network(net = dl_rvl[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

dl_rvl[[4]] <- attribute_in_network(net = dl_rvl[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  dl_rvl[[5]] <- attribute_in_network(net = dl_rvl[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

#

dl_jbs[[1]] <- attribute_in_network(net = dl_jbs[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_jbs[[2]] <- attribute_in_network(net = dl_jbs[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_jbs[[3]] <- attribute_in_network(net = dl_jbs[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

dl_jbs[[4]] <- attribute_in_network(net = dl_jbs[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  dl_jbs[[5]] <- attribute_in_network(net = dl_jbs[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

#

dl_pm[[1]] <- attribute_in_network(net = dl_pm[[1]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

dl_pm[[2]] <- attribute_in_network(net = dl_pm[[2]], data = psycho, id_data = c("id"),
                                   var_names = a2,
                                   new_names = b2)

dl_pm[[3]] <- attribute_in_network(net = dl_pm[[3]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

dl_pm[[4]] <- attribute_in_network(net = dl_pm[[4]], data = psycho, id_data = c("id"),
                                   var_names = a4,
                                   new_names = b4)

if(use_psycho=="new"){
  dl_pm[[5]] <- attribute_in_network(net = dl_pm[[5]], data = psycho, id_data = c("id"),
                                     var_names = a6,
                                     new_names = b6)
  
}

#

dl_dmz[[1]] <- attribute_in_network(net = dl_dmz[[1]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_dmz[[2]] <- attribute_in_network(net = dl_dmz[[2]], data = psycho, id_data = c("id"),
                                    var_names = a2,
                                    new_names = b2)

dl_dmz[[3]] <- attribute_in_network(net = dl_dmz[[3]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

dl_dmz[[4]] <- attribute_in_network(net = dl_dmz[[4]], data = psycho, id_data = c("id"),
                                    var_names = a4,
                                    new_names = b4)

if(use_psycho=="new"){
  dl_dmz[[5]] <- attribute_in_network(net = dl_dmz[[5]], data = psycho, id_data = c("id"),
                                      var_names = a6,
                                      new_names = b6)
}

###

top5_rvl[[1]] <- attribute_in_network(net = top5_rvl[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_rvl[[2]] <- attribute_in_network(net = top5_rvl[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_rvl[[3]] <- attribute_in_network(net = top5_rvl[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

top5_rvl[[4]] <- attribute_in_network(net = top5_rvl[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  top5_rvl[[5]] <- attribute_in_network(net = top5_rvl[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

top5_jbs[[1]] <- attribute_in_network(net = top5_jbs[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_jbs[[2]] <- attribute_in_network(net = top5_jbs[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_jbs[[3]] <- attribute_in_network(net = top5_jbs[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

top5_jbs[[4]] <- attribute_in_network(net = top5_jbs[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  top5_jbs[[5]] <- attribute_in_network(net = top5_jbs[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

top5_pm[[1]] <- attribute_in_network(net = top5_pm[[1]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

top5_pm[[2]] <- attribute_in_network(net = top5_pm[[2]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

top5_pm[[3]] <- attribute_in_network(net = top5_pm[[3]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

top5_pm[[4]] <- attribute_in_network(net = top5_pm[[4]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

if(use_psycho=="new"){
  top5_pm[[5]] <- attribute_in_network(net = top5_pm[[5]], data = psycho, id_data = c("id"),
                                       var_names = a6,
                                       new_names = b6)
}

#

top5_dmz[[1]] <- attribute_in_network(net = top5_dmz[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_dmz[[2]] <- attribute_in_network(net = top5_dmz[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

top5_dmz[[3]] <- attribute_in_network(net = top5_dmz[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

top5_dmz[[4]] <- attribute_in_network(net = top5_dmz[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  top5_dmz[[5]] <- attribute_in_network(net = top5_dmz[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
  
}

###

bully_rvl[[1]] <- attribute_in_network(net = bully_rvl[[1]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_rvl[[2]] <- attribute_in_network(net = bully_rvl[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_rvl[[3]] <- attribute_in_network(net = bully_rvl[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

bully_rvl[[4]] <- attribute_in_network(net = bully_rvl[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  bully_rvl[[5]] <- attribute_in_network(net = bully_rvl[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

bully_jbs[[1]] <- attribute_in_network(net = bully_jbs[[1]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_jbs[[2]] <- attribute_in_network(net = bully_jbs[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_jbs[[3]] <- attribute_in_network(net = bully_jbs[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

bully_jbs[[4]] <- attribute_in_network(net = bully_jbs[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  bully_jbs[[5]] <- attribute_in_network(net = bully_jbs[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

bully_pm[[1]] <- attribute_in_network(net = bully_pm[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

bully_pm[[2]] <- attribute_in_network(net = bully_pm[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

bully_pm[[3]] <- attribute_in_network(net = bully_pm[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

bully_pm[[4]] <- attribute_in_network(net = bully_pm[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  bully_pm[[5]] <- attribute_in_network(net = bully_pm[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

bully_dmz[[1]] <- attribute_in_network(net = bully_dmz[[1]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_dmz[[2]] <- attribute_in_network(net = bully_dmz[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

bully_dmz[[3]] <- attribute_in_network(net = bully_dmz[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

bully_dmz[[4]] <- attribute_in_network(net = bully_dmz[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  bully_dmz[[5]] <- attribute_in_network(net = bully_dmz[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

###

cool_rvl[[1]] <- attribute_in_network(net = cool_rvl[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_rvl[[2]] <- attribute_in_network(net = cool_rvl[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_rvl[[3]] <- attribute_in_network(net = cool_rvl[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

cool_rvl[[4]] <- attribute_in_network(net = cool_rvl[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  cool_rvl[[5]] <- attribute_in_network(net = cool_rvl[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

cool_jbs[[1]] <- attribute_in_network(net = cool_jbs[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_jbs[[2]] <- attribute_in_network(net = cool_jbs[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_jbs[[3]] <- attribute_in_network(net = cool_jbs[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

cool_jbs[[4]] <- attribute_in_network(net = cool_jbs[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  cool_jbs[[5]] <- attribute_in_network(net = cool_jbs[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

cool_pm[[1]] <- attribute_in_network(net = cool_pm[[1]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

cool_pm[[2]] <- attribute_in_network(net = cool_pm[[2]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

cool_pm[[3]] <- attribute_in_network(net = cool_pm[[3]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

cool_pm[[4]] <- attribute_in_network(net = cool_pm[[4]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

if(use_psycho=="new"){
  cool_pm[[5]] <- attribute_in_network(net = cool_pm[[5]], data = psycho, id_data = c("id"),
                                       var_names = a6,
                                       new_names = b6)
}

#

cool_dmz[[1]] <- attribute_in_network(net = cool_dmz[[1]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_dmz[[2]] <- attribute_in_network(net = cool_dmz[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

cool_dmz[[3]] <- attribute_in_network(net = cool_dmz[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

cool_dmz[[4]] <- attribute_in_network(net = cool_dmz[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  cool_dmz[[5]] <- attribute_in_network(net = cool_dmz[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

###

known_rvl[[2]] <- attribute_in_network(net = known_rvl[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

known_rvl[[3]] <- attribute_in_network(net = known_rvl[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

known_rvl[[4]] <- attribute_in_network(net = known_rvl[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  known_rvl[[5]] <- attribute_in_network(net = known_rvl[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

known_jbs[[2]] <- attribute_in_network(net = known_jbs[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

known_jbs[[3]] <- attribute_in_network(net = known_jbs[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

known_jbs[[4]] <- attribute_in_network(net = known_jbs[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  known_jbs[[5]] <- attribute_in_network(net = known_jbs[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

known_pm[[2]] <- attribute_in_network(net = known_pm[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

known_pm[[3]] <- attribute_in_network(net = known_pm[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

known_pm[[4]] <- attribute_in_network(net = known_pm[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  known_pm[[5]] <- attribute_in_network(net = known_pm[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

known_dmz[[2]] <- attribute_in_network(net = known_dmz[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

known_dmz[[3]] <- attribute_in_network(net = known_dmz[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

known_dmz[[4]] <- attribute_in_network(net = known_dmz[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  known_dmz[[5]] <- attribute_in_network(net = known_dmz[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

###

fear_rvl[[2]] <- attribute_in_network(net = fear_rvl[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

fear_rvl[[3]] <- attribute_in_network(net = fear_rvl[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

fear_rvl[[4]] <- attribute_in_network(net = fear_rvl[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  fear_rvl[[5]] <- attribute_in_network(net = fear_rvl[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

fear_jbs[[2]] <- attribute_in_network(net = fear_jbs[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

fear_jbs[[3]] <- attribute_in_network(net = fear_jbs[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

fear_jbs[[4]] <- attribute_in_network(net = fear_jbs[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  fear_jbs[[5]] <- attribute_in_network(net = fear_jbs[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

fear_pm[[2]] <- attribute_in_network(net = fear_pm[[2]], data = psycho, id_data = c("id"),
                                     var_names = a2,
                                     new_names = b2)

fear_pm[[3]] <- attribute_in_network(net = fear_pm[[3]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

fear_pm[[4]] <- attribute_in_network(net = fear_pm[[4]], data = psycho, id_data = c("id"),
                                     var_names = a4,
                                     new_names = b4)

if(use_psycho=="new"){
  fear_pm[[5]] <- attribute_in_network(net = fear_pm[[5]], data = psycho, id_data = c("id"),
                                       var_names = a6,
                                       new_names = b6)
}

#

fear_dmz[[2]] <- attribute_in_network(net = fear_dmz[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

fear_dmz[[3]] <- attribute_in_network(net = fear_dmz[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

fear_dmz[[4]] <- attribute_in_network(net = fear_dmz[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  fear_dmz[[5]] <- attribute_in_network(net = fear_dmz[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

###

liked_rvl[[2]] <- attribute_in_network(net = liked_rvl[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

liked_rvl[[3]] <- attribute_in_network(net = liked_rvl[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

liked_rvl[[4]] <- attribute_in_network(net = liked_rvl[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  liked_rvl[[5]] <- attribute_in_network(net = liked_rvl[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

liked_jbs[[2]] <- attribute_in_network(net = liked_jbs[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

liked_jbs[[3]] <- attribute_in_network(net = liked_jbs[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

liked_jbs[[4]] <- attribute_in_network(net = liked_jbs[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  liked_jbs[[5]] <- attribute_in_network(net = liked_jbs[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

#

liked_pm[[2]] <- attribute_in_network(net = liked_pm[[2]], data = psycho, id_data = c("id"),
                                      var_names = a2,
                                      new_names = b2)

liked_pm[[3]] <- attribute_in_network(net = liked_pm[[3]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

liked_pm[[4]] <- attribute_in_network(net = liked_pm[[4]], data = psycho, id_data = c("id"),
                                      var_names = a4,
                                      new_names = b4)

if(use_psycho=="new"){
  liked_pm[[5]] <- attribute_in_network(net = liked_pm[[5]], data = psycho, id_data = c("id"),
                                        var_names = a6,
                                        new_names = b6)
}

#

liked_dmz[[2]] <- attribute_in_network(net = liked_dmz[[2]], data = psycho, id_data = c("id"),
                                       var_names = a2,
                                       new_names = b2)

liked_dmz[[3]] <- attribute_in_network(net = liked_dmz[[3]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

liked_dmz[[4]] <- attribute_in_network(net = liked_dmz[[4]], data = psycho, id_data = c("id"),
                                       var_names = a4,
                                       new_names = b4)

if(use_psycho=="new"){
  liked_dmz[[5]] <- attribute_in_network(net = liked_dmz[[5]], data = psycho, id_data = c("id"),
                                         var_names = a6,
                                         new_names = b6)
}

########## 6. Matrices for homophily and activity terms (unactive now - use ctrl+c to activate again) ##########

# btergm needs matrices if we are to interact those with timecov().
# I use my custom function "mat_att".

### Factbrev matrices

# mat_diff_factbrev_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"factbrev")$diff})
# mat_diff_factbrev_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"factbrev")$diff})
# mat_diff_factbrev_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"factbrev")$diff})
# mat_diff_factbrev_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"factbrev")$diff})
# 
# mat_sender_factbrev_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"factbrev")$sender})
# mat_sender_factbrev_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"factbrev")$sender})
# mat_sender_factbrev_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"factbrev")$sender})
# mat_sender_factbrev_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"factbrev")$sender})
# 
# mat_receiver_factbrev_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"factbrev")$receiver})
# mat_receiver_factbrev_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"factbrev")$receiver})
# mat_receiver_factbrev_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"factbrev")$receiver})
# mat_receiver_factbrev_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"factbrev")$receiver})
# 
# ### Sex matrices (for factbrev networks)
# 
# mat_match_sex_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"sex")$match})
# mat_match_sex_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"sex")$match})
# mat_match_sex_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"sex")$match})
# mat_match_sex_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"sex")$match})
# 
# # NB : for sender and receiver, we just need the boy matrices ; girls will be the ref category.
# mat_Mreceiver_sex_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"sex")$M_receiver})
# mat_Mreceiver_sex_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"sex")$M_receiver})
# mat_Mreceiver_sex_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"sex")$M_receiver})
# mat_Mreceiver_sex_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"sex")$M_receiver})
# 
# mat_Msender_sex_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"sex")$M_sender})
# mat_Msender_sex_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"sex")$M_sender})
# mat_Msender_sex_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"sex")$M_sender})
# mat_Msender_sex_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"sex")$M_sender})
# 
# ### Classroom matrices (for factbrev networks)
# 
# mat_match_classroom_rvl <- lapply(gf_rvl_cfb,function(x){mat_att(x,"classroom")$match})
# mat_match_classroom_jbs <- lapply(gf_jbs_cfb,function(x){mat_att(x,"classroom")$match})
# mat_match_classroom_pm <- lapply(gf_pm_cfb,function(x){mat_att(x,"classroom")$match})
# mat_match_classroom_dmz <- lapply(gf_dmz_cfb,function(x){mat_att(x,"classroom")$match})

########## 7. Ind data frames for all schools at once ##########

# NB : for now some parts don t work because isei and sios are not there yet.
# SO for now it s all in # in the loops.

# First i ll make sure there is a "school" variable everywhere.
rvl1$school <- "rvl"
jbs1$school <- "jbs"
pm1$school <- "pm"
dmz1$school <- "dmz"

rvl2$school <- "rvl"
jbs2$school <- "jbs"
pm2$school <- "pm"
dmz2$school <- "dmz"

rvl3$school <- "rvl"
jbs3$school <- "jbs"
pm3$school <- "pm"
dmz3$school <- "dmz"

rvl4$school <- "rvl"
jbs4$school <- "jbs"
pm4$school <- "pm"
dmz4$school <- "dmz"

rvl6$school <- "rvl"
jbs6$school <- "jbs"
pm6$school <- "pm"
dmz6$school <- "dmz"

##

names1 <- intersect(intersect(intersect(names(rvl1),names(jbs1)),names(pm1)),names(dmz1))
names2 <- intersect(intersect(intersect(names(rvl2),names(jbs2)),names(pm2)),names(dmz2))
names3 <- intersect(intersect(intersect(names(rvl3),names(jbs3)),names(pm3)),names(dmz3))
names4 <- intersect(intersect(intersect(names(rvl4),names(jbs4)),names(pm4)),names(dmz4))
names6 <- intersect(intersect(intersect(names(rvl6),names(jbs6)),names(pm6)),names(dmz6))

ind1 <- rbind(rvl1[,names1],jbs1[,names1],pm1[,names1],dmz1[,names1])
ind2 <- rbind(rvl2[,names2],jbs2[,names2],pm2[,names2],dmz2[,names2])
ind3 <- rbind(rvl3[,names3],jbs3[,names3],pm3[,names3],dmz3[,names3])
ind4 <- rbind(rvl4[,names4],jbs4[,names4],pm4[,names4],dmz4[,names4])
ind6 <- rbind(rvl6[,names6],jbs6[,names6],pm6[,names6],dmz6[,names6])

###

# Put the occupation variables in those.
# NB : for quant scores we also want the squared version (for non-linear effects).

ind1$factbrev <- NA
ind1$isei <- NA
ind1$siops <- NA
ind1$factbrev_sq <- NA
ind1$isei_sq <- NA
ind1$siops_sq <- NA
ind1$suli <- NA

ind2$factbrev <- NA
ind2$isei <- NA
ind2$siops <- NA
ind2$factbrev_sq <- NA
ind2$isei_sq <- NA
ind2$siops_sq <- NA
ind2$suli <- NA

ind3$factbrev <- NA
ind3$isei <- NA
ind3$siops <- NA
ind3$factbrev_sq <- NA
ind3$isei_sq <- NA
ind3$siops_sq <- NA
ind3$suli <- NA

ind4$factbrev <- NA
ind4$isei <- NA
ind4$siops <- NA
ind4$factbrev_sq <- NA
ind4$isei_sq <- NA
ind4$siops_sq <- NA
ind4$suli <- NA

ind6$factbrev <- NA
ind6$isei <- NA
ind6$siops <- NA
ind6$factbrev_sq <- NA
ind6$isei_sq <- NA
ind6$siops_sq <- NA
ind6$suli <- NA

for(i in 1:nrow(ind1)){
  id <- ind1$id[i]
  ind1$factbrev[i] <- occupation$factbrev[occupation$id==id]
  ind1$isei[i] <- occupation$isei2[occupation$id==id]
  # ind1$siops[i] <- occupation$siops1[occupation$id==id]
  # ind1$suli[i] <- occupation$suli[occupation$id==id]
}

for(i in 1:nrow(ind2)){
  id <- ind2$id[i]
  ind2$factbrev[i] <- occupation$factbrev[occupation$id==id]
  ind2$isei[i] <- occupation$isei2[occupation$id==id]
  # ind2$siops[i] <- occupation$siops1[occupation$id==id]
  # ind2$suli[i] <- occupation$suli[occupation$id==id]
}

for(i in 1:nrow(ind3)){
  id <- ind3$id[i]
  ind3$factbrev[i] <- occupation$factbrev[occupation$id==id]
  ind3$isei[i] <- occupation$isei2[occupation$id==id]
  # ind3$siops[i] <- occupation$siops1[occupation$id==id]
  # ind3$suli[i] <- occupation$suli[occupation$id==id]
}

for(i in 1:nrow(ind4)){
  id <- ind4$id[i]
  ind4$factbrev[i] <- occupation$factbrev[occupation$id==id]
  ind4$isei[i] <- occupation$isei2[occupation$id==id]
  # ind4$siops[i] <- occupation$siops1[occupation$id==id]
  # ind4$suli[i] <- occupation$suli[occupation$id==id]
}

for(i in 1:nrow(ind6)){
  id <- ind6$id[i]
  ind6$factbrev[i] <- occupation$factbrev[occupation$id==id]
  ind6$isei[i] <- occupation$isei2[occupation$id==id]
  # ind6$siops[i] <- occupation$siops1[occupation$id==id]
  # ind6$suli[i] <- occupation$suli[occupation$id==id]
}

ind1$factbrev_sq <- ind1$factbrev * ind1$factbrev
ind1$isei_sq <- ind1$isei * ind1$isei
ind1$siops_sq <- ind1$siops * ind1$siops

ind2$factbrev_sq <- ind2$factbrev * ind2$factbrev
ind2$isei_sq <- ind2$isei * ind2$isei
ind2$siops_sq <- ind2$siops * ind2$siops

ind3$factbrev_sq <- ind3$factbrev * ind3$factbrev
ind3$isei_sq <- ind3$isei * ind3$isei
ind3$siops_sq <- ind3$siops * ind3$siops

ind4$factbrev_sq <- ind4$factbrev * ind4$factbrev
ind4$isei_sq <- ind4$isei * ind4$isei
ind4$siops_sq <- ind4$siops * ind4$siops

ind6$factbrev_sq <- ind6$factbrev * ind6$factbrev
ind6$isei_sq <- ind6$isei * ind6$isei
ind6$siops_sq <- ind6$siops * ind6$siops

###

# And also put occupation in the per-school data frames :

rvl1$factbrev <- NA
rvl1$isei <- NA
rvl1$siops <- NA
rvl1$suli <- NA

for(i in 1:nrow(rvl1)){
  id <- rvl1$id[i]
  rvl1$factbrev[i] <- occupation$factbrev[occupation$id==id]
  rvl1$isei[i] <- occupation$isei2[occupation$id==id]
  # rvl1$siops[i] <- occupation$siops1[occupation$id==id]
  # rvl1$suli[i] <- occupation$suli[occupation$id==id]
}

rvl1$factbrev_sq <- rvl1$factbrev * rvl1$factbrev
rvl1$isei_sq <- rvl1$isei * rvl1$isei
rvl1$siops_sq <- rvl1$siops * rvl1$siops

rvl2$factbrev <- NA
rvl2$isei <- NA
rvl2$siops <- NA
rvl2$suli <- NA

for(i in 1:nrow(rvl2)){
  id <- rvl2$id[i]
  rvl2$factbrev[i] <- occupation$factbrev[occupation$id==id]
  rvl2$isei[i] <- occupation$isei2[occupation$id==id]
  # rvl2$siops[i] <- occupation$siops1[occupation$id==id]
  # rvl2$suli[i] <- occupation$suli[occupation$id==id]
}

rvl2$factbrev_sq <- rvl2$factbrev * rvl2$factbrev
rvl2$isei_sq <- rvl2$isei * rvl2$isei
rvl2$siops_sq <- rvl2$siops * rvl2$siops

rvl3$factbrev <- NA
rvl3$isei <- NA
rvl3$siops <- NA
rvl3$suli <- NA

for(i in 1:nrow(rvl3)){
  id <- rvl3$id[i]
  rvl3$factbrev[i] <- occupation$factbrev[occupation$id==id]
  rvl3$isei[i] <- occupation$isei2[occupation$id==id]
  # rvl3$siops[i] <- occupation$siops1[occupation$id==id]
  # rvl3$suli[i] <- occupation$suli[occupation$id==id]
}

rvl3$factbrev_sq <- rvl3$factbrev * rvl3$factbrev
rvl3$isei_sq <- rvl3$isei * rvl3$isei
rvl3$siops_sq <- rvl3$siops * rvl3$siops

rvl4$factbrev <- NA
rvl4$isei <- NA
rvl4$siops <- NA
rvl4$suli <- NA

for(i in 1:nrow(rvl4)){
  id <- rvl4$id[i]
  rvl4$factbrev[i] <- occupation$factbrev[occupation$id==id]
  rvl4$isei[i] <- occupation$isei2[occupation$id==id]
  # rvl4$siops[i] <- occupation$siops1[occupation$id==id]
  # rvl4$suli[i] <- occupation$suli[occupation$id==id]
}

rvl4$factbrev_sq <- rvl4$factbrev * rvl4$factbrev
rvl4$isei_sq <- rvl4$isei * rvl4$isei
rvl4$siops_sq <- rvl4$siops * rvl4$siops

rvl6$factbrev <- NA
rvl6$isei <- NA
rvl6$siops <- NA
rvl6$suli <- NA

for(i in 1:nrow(rvl6)){
  id <- rvl6$id[i]
  rvl6$factbrev[i] <- occupation$factbrev[occupation$id==id]
  rvl6$isei[i] <- occupation$isei2[occupation$id==id]
  # rvl6$siops[i] <- occupation$siops1[occupation$id==id]
  # rvl6$suli[i] <- occupation$suli[occupation$id==id]
}

rvl6$factbrev_sq <- rvl6$factbrev * rvl6$factbrev
rvl6$isei_sq <- rvl6$isei * rvl6$isei
rvl6$siops_sq <- rvl6$siops * rvl6$siops

##


jbs1$factbrev <- NA
jbs1$isei <- NA
jbs1$siops <- NA
jbs1$suli <- NA

for(i in 1:nrow(jbs1)){
  id <- jbs1$id[i]
  jbs1$factbrev[i] <- occupation$factbrev[occupation$id==id]
  jbs1$isei[i] <- occupation$isei2[occupation$id==id]
  # jbs1$siops[i] <- occupation$siops1[occupation$id==id]
  # jbs1$suli[i] <- occupation$suli[occupation$id==id]
}

jbs1$factbrev_sq <- jbs1$factbrev * jbs1$factbrev
jbs1$isei_sq <- jbs1$isei * jbs1$isei
jbs1$siops_sq <- jbs1$siops * jbs1$siops

jbs2$factbrev <- NA
jbs2$isei <- NA
jbs2$siops <- NA
jbs2$suli <- NA

for(i in 1:nrow(jbs2)){
  id <- jbs2$id[i]
  jbs2$factbrev[i] <- occupation$factbrev[occupation$id==id]
  jbs2$isei[i] <- occupation$isei2[occupation$id==id]
  # jbs2$siops[i] <- occupation$siops1[occupation$id==id]
  # jbs2$suli[i] <- occupation$suli[occupation$id==id]
}

jbs2$factbrev_sq <- jbs2$factbrev * jbs2$factbrev
jbs2$isei_sq <- jbs2$isei * jbs2$isei
jbs2$siops_sq <- jbs2$siops * jbs2$siops

jbs3$factbrev <- NA
jbs3$isei <- NA
jbs3$siops <- NA
jbs3$suli <- NA

for(i in 1:nrow(jbs3)){
  id <- jbs3$id[i]
  jbs3$factbrev[i] <- occupation$factbrev[occupation$id==id]
  jbs3$isei[i] <- occupation$isei2[occupation$id==id]
  # jbs3$siops[i] <- occupation$siops1[occupation$id==id]
  # jbs3$suli[i] <- occupation$suli[occupation$id==id]
}

jbs3$factbrev_sq <- jbs3$factbrev * jbs3$factbrev
jbs3$isei_sq <- jbs3$isei * jbs3$isei
jbs3$siops_sq <- jbs3$siops * jbs3$siops

jbs4$factbrev <- NA
jbs4$isei <- NA
jbs4$siops <- NA
jbs4$suli <- NA

for(i in 1:nrow(jbs4)){
  id <- jbs4$id[i]
  jbs4$factbrev[i] <- occupation$factbrev[occupation$id==id]
  jbs4$isei[i] <- occupation$isei2[occupation$id==id]
  # jbs4$siops[i] <- occupation$siops1[occupation$id==id]
  # jbs4$suli[i] <- occupation$suli[occupation$id==id]
}

jbs4$factbrev_sq <- jbs4$factbrev * jbs4$factbrev
jbs4$isei_sq <- jbs4$isei * jbs4$isei
jbs4$siops_sq <- jbs4$siops * jbs4$siops

jbs6$factbrev <- NA
jbs6$isei <- NA
jbs6$siops <- NA
jbs6$suli <- NA

for(i in 1:nrow(jbs6)){
  id <- jbs6$id[i]
  jbs6$factbrev[i] <- occupation$factbrev[occupation$id==id]
  jbs6$isei[i] <- occupation$isei2[occupation$id==id]
  # jbs6$siops[i] <- occupation$siops1[occupation$id==id]
  # jbs6$suli[i] <- occupation$suli[occupation$id==id]
}

jbs6$factbrev_sq <- jbs6$factbrev * jbs6$factbrev
jbs6$isei_sq <- jbs6$isei * jbs6$isei
jbs6$siops_sq <- jbs6$siops * jbs6$siops

##


pm1$factbrev <- NA
pm1$isei <- NA
pm1$siops <- NA
pm1$suli <- NA

for(i in 1:nrow(pm1)){
  id <- pm1$id[i]
  pm1$factbrev[i] <- occupation$factbrev[occupation$id==id]
  pm1$isei[i] <- occupation$isei2[occupation$id==id]
  # pm1$siops[i] <- occupation$siops1[occupation$id==id]
  # pm1$suli[i] <- occupation$suli[occupation$id==id]
}

pm1$factbrev_sq <- pm1$factbrev * pm1$factbrev
pm1$isei_sq <- pm1$isei * pm1$isei
pm1$siops_sq <- pm1$siops * pm1$siops

pm2$factbrev <- NA
pm2$isei <- NA
pm2$siops <- NA
pm2$suli <- NA

for(i in 1:nrow(pm2)){
  id <- pm2$id[i]
  pm2$factbrev[i] <- occupation$factbrev[occupation$id==id]
  pm2$isei[i] <- occupation$isei2[occupation$id==id]
  # pm2$siops[i] <- occupation$siops1[occupation$id==id]
  # pm2$suli[i] <- occupation$suli[occupation$id==id]
}

pm2$factbrev_sq <- pm2$factbrev * pm2$factbrev
pm2$isei_sq <- pm2$isei * pm2$isei
pm2$siops_sq <- pm2$siops * pm2$siops

pm3$factbrev <- NA
pm3$isei <- NA
pm3$siops <- NA
pm3$suli <- NA

for(i in 1:nrow(pm3)){
  id <- pm3$id[i]
  pm3$factbrev[i] <- occupation$factbrev[occupation$id==id]
  pm3$isei[i] <- occupation$isei2[occupation$id==id]
  # pm3$siops[i] <- occupation$siops1[occupation$id==id]
  # pm3$suli[i] <- occupation$suli[occupation$id==id]
}

pm3$factbrev_sq <- pm3$factbrev * pm3$factbrev
pm3$isei_sq <- pm3$isei * pm3$isei
pm3$siops_sq <- pm3$siops * pm3$siops

pm4$factbrev <- NA
pm4$isei <- NA
pm4$siops <- NA
pm4$suli <- NA

for(i in 1:nrow(pm4)){
  id <- pm4$id[i]
  pm4$factbrev[i] <- occupation$factbrev[occupation$id==id]
  pm4$isei[i] <- occupation$isei2[occupation$id==id]
  # pm4$siops[i] <- occupation$siops1[occupation$id==id]
  # pm4$suli[i] <- occupation$suli[occupation$id==id]
}

pm4$factbrev_sq <- pm4$factbrev * pm4$factbrev
pm4$isei_sq <- pm4$isei * pm4$isei
pm4$siops_sq <- pm4$siops * pm4$siops

pm6$factbrev <- NA
pm6$isei <- NA
pm6$siops <- NA
pm6$suli <- NA

for(i in 1:nrow(pm6)){
  id <- pm6$id[i]
  pm6$factbrev[i] <- occupation$factbrev[occupation$id==id]
  pm6$isei[i] <- occupation$isei2[occupation$id==id]
  # pm6$siops[i] <- occupation$siops1[occupation$id==id]
  # pm6$suli[i] <- occupation$suli[occupation$id==id]
}

pm6$factbrev_sq <- pm6$factbrev * pm6$factbrev
pm6$isei_sq <- pm6$isei * pm6$isei
pm6$siops_sq <- pm6$siops * pm6$siops

##


dmz1$factbrev <- NA
dmz1$isei <- NA
dmz1$siops <- NA
dmz1$suli <- NA

for(i in 1:nrow(dmz1)){
  id <- dmz1$id[i]
  dmz1$factbrev[i] <- occupation$factbrev[occupation$id==id]
  dmz1$isei[i] <- occupation$isei2[occupation$id==id]
  # dmz1$siops[i] <- occupation$siops1[occupation$id==id]
  # dmz1$suli[i] <- occupation$suli[occupation$id==id]
}

dmz1$factbrev_sq <- dmz1$factbrev * dmz1$factbrev
dmz1$isei_sq <- dmz1$isei * dmz1$isei
dmz1$siops_sq <- dmz1$siops * dmz1$siops

dmz2$factbrev <- NA
dmz2$isei <- NA
dmz2$siops <- NA
dmz2$suli <- NA

for(i in 1:nrow(dmz2)){
  id <- dmz2$id[i]
  dmz2$factbrev[i] <- occupation$factbrev[occupation$id==id]
  dmz2$isei[i] <- occupation$isei2[occupation$id==id]
  # dmz2$siops[i] <- occupation$siops1[occupation$id==id]
  # dmz2$suli[i] <- occupation$suli[occupation$id==id]
}

dmz2$factbrev_sq <- dmz2$factbrev * dmz2$factbrev
dmz2$isei_sq <- dmz2$isei * dmz2$isei
dmz2$siops_sq <- dmz2$siops * dmz2$siops

dmz3$factbrev <- NA
dmz3$isei <- NA
dmz3$siops <- NA
dmz3$suli <- NA

for(i in 1:nrow(dmz3)){
  id <- dmz3$id[i]
  dmz3$factbrev[i] <- occupation$factbrev[occupation$id==id]
  dmz3$isei[i] <- occupation$isei2[occupation$id==id]
  # dmz3$siops[i] <- occupation$siops1[occupation$id==id]
  # dmz3$suli[i] <- occupation$suli[occupation$id==id]
}

dmz3$factbrev_sq <- dmz3$factbrev * dmz3$factbrev
dmz3$isei_sq <- dmz3$isei * dmz3$isei
dmz3$siops_sq <- dmz3$siops * dmz3$siops

dmz4$factbrev <- NA
dmz4$isei <- NA
dmz4$siops <- NA
dmz4$suli <- NA

for(i in 1:nrow(dmz4)){
  id <- dmz4$id[i]
  dmz4$factbrev[i] <- occupation$factbrev[occupation$id==id]
  dmz4$isei[i] <- occupation$isei2[occupation$id==id]
  # dmz4$siops[i] <- occupation$siops1[occupation$id==id]
  # dmz4$suli[i] <- occupation$suli[occupation$id==id]
}

dmz4$factbrev_sq <- dmz4$factbrev * dmz4$factbrev
dmz4$isei_sq <- dmz4$isei * dmz4$isei
dmz4$siops_sq <- dmz4$siops * dmz4$siops

dmz6$factbrev <- NA
dmz6$isei <- NA
dmz6$siops <- NA
dmz6$suli <- NA

for(i in 1:nrow(dmz6)){
  id <- dmz6$id[i]
  dmz6$factbrev[i] <- occupation$factbrev[occupation$id==id]
  dmz6$isei[i] <- occupation$isei2[occupation$id==id]
  # dmz6$siops[i] <- occupation$siops1[occupation$id==id]
  # dmz6$suli[i] <- occupation$suli[occupation$id==id]
}

dmz6$factbrev_sq <- dmz6$factbrev * dmz6$factbrev
dmz6$isei_sq <- dmz6$isei * dmz6$isei
dmz6$siops_sq <- dmz6$siops * dmz6$siops

########## 8. Tie change statistics in the networks ##########

# Use the function "get_cstat" to have a data frame with the change stat of each tie-variable, on a number
# of possible model parameters.

# Note that this is long to run, so I exported the result and just reload them here (original code in comments below).
# Also, I get the cstats for terms I know will end up in many models (e.g. absdiff factbrev or mutual) but not for
# ones that need to be adjusted to specific models (e.g. gwesp with specific decays).

# This does not demand to run estimations, so I could include terms that are generally not converging in ERGMs (triangle
# in particular). Having the change stat will still allow me to see how triangle descriptively correlates with certain
# homophilies, for e.g.

###

# Make the environment a bit lighter :
# rm(culture,psycho,rvl1,rvl2,rvl3,rvl4,jbs1,jbs2,jbs3,jbs4,pm1,pm2,pm3,pm4,dmz1,dmz2,dmz3,dmz4)
# rm(ind1,ind2,ind3,ind4)
# rm(allf_rvl,allf_jbs,allf_pm,allf_jbs)
# rm(allf_rvl_cfb,allf_jbs_cfb,allf_pm_cfb,allf_jbs_cfb)
# rm(dl_rvl,dl_jbs,dl_pm,dl_jbs)
# rm(dl_rvl_cfb,dl_jbs_cfb,dl_pm_cfb,dl_jbs_cfb)
# rm(top5_rvl,top5_jbs,top5_pm,top5_jbs)
# rm(top5_rvl_cfb,top5_jbs_cfb,top5_pm_cfb,top5_jbs_cfb)
# rm(cool_rvl,cool_jbs,cool_pm,cool_jbs)
# rm(cool_rvl_cfb,cool_jbs_cfb,cool_pm_cfb,cool_jbs_cfb)
# rm(mat_diff_factbrev_rvl,mat_diff_factbrev_jbs,mat_diff_factbrev_pm,mat_diff_factbrev_dmz)
# rm(mat_match_classroom_rvl,mat_match_classroom_jbs,mat_match_classroom_pm,mat_match_classroom_dmz)
# rm(mat_match_sex_rvl,mat_match_sex_jbs,mat_match_sex_pm,mat_match_sex_dmz)
# 
# # Hell on earth :
# 
# myformula <- net ~ edges +
#   absdiff("factbrev") + nodeicov("factbrev") + nodeocov("factbrev") +
#   absdiff("grade") + nodeicov("grade") + nodeocov("grade") +
#   nodematch("ethn_ag",diff=T) + nodeifactor("ethn_ag") + nodeofactor("ethn_ag") +
#   nodematch("sex") + nodeifactor("sex") + nodeofactor("sex") +
#   nodematch("primary1") +
#   nodematch("classroom") + nodeifactor("classroom") + nodeofactor("classroom") +
#   absdiff("art_lvl") + nodeicov("art_lvl") + nodeocov("art_lvl") +
#   absdiff("music_lvl") + nodeicov("music_lvl") + nodeocov("music_lvl") +
#   absdiff("book_lvl") + nodeicov("book_lvl") + nodeocov("book_lvl") +
#   absdiff("cult_lvl") + nodeicov("cult_lvl") + nodeocov("cult_lvl") +
#   absdiff("m_fact1") + nodeicov("m_fact1") + nodeocov("m_fact1") +
#   absdiff("m_fact2") + nodeicov("m_fact2") + nodeocov("m_fact2") +
#   absdiff("y_fact1") + nodeicov("y_fact1") + nodeocov("y_fact1") +
#   absdiff("y_fact2") + nodeicov("y_fact2") + nodeocov("y_fact2") +
#   absdiff("ym_fact1") + nodeicov("ym_fact1") + nodeocov("ym_fact1") +
#   absdiff("ym_fact2") + nodeicov("ym_fact2") + nodeocov("ym_fact2") +
#   absdiff("tree") + nodeicov("tree") + nodeocov("tree") +
#   absdiff("forest") + nodeicov("forest") + nodeocov("forest") +
#   absdiff("psy_fact1") + nodeicov("psy_fact1") + nodeocov("psy_fact1") +
#   absdiff("psy_fact2") + nodeicov("psy_fact2") + nodeocov("psy_fact2") +
#   absdiff("psy_fact3") + nodeicov("psy_fact3") + nodeocov("psy_fact3") +
#   absdiff("extrav") + nodeicov("extrav") + nodeocov("extrav") +
#   absdiff("agreeab") + nodeicov("agreeab") + nodeocov("agreeab") +
#   absdiff("consc") + nodeicov("consc") + nodeocov("consc") +
#   absdiff("neuro") + nodeicov("neuro") + nodeocov("neuro") +
#   absdiff("open") + nodeicov("open") + nodeocov("open") +
#   mutual + triangle
# 
# myform <- update.formula(myformula,gf_rvl[[1]]~.)
# cstat_rvl1 <- get_cstat_handleNA(net = gf_rvl[[1]], formula = myform)
# save(cstat_rvl1, file = "./cstat_rvl1.RData")
# rm(cstat_rvl1)
# 
# myform <- update.formula(myformula,gf_jbs[[1]]~.)
# cstat_jbs1 <- get_cstat_handleNA(net = gf_jbs[[1]], formula = myform)
# save(cstat_jbs1, file = "./cstat_jbs1.RData")
# rm(cstat_jbs1)
# 
# myform <- update.formula(myformula,gf_pm[[1]]~.)
# cstat_pm1 <- get_cstat_handleNA(net = gf_pm[[1]], formula = myform)
# save(cstat_pm1, file = "./cstat_pm1.RData")
# rm(cstat_pm1)
# 
# myform <- update.formula(myformula,gf_dmz[[1]]~.)
# cstat_dmz1 <- get_cstat_handleNA(net = gf_dmz[[1]], formula = myform)
# save(cstat_dmz1, file = "./cstat_dmz1.RData")
# rm(cstat_dmz1)

###

# Now merge those lists to a single data frame per school (with NAs when a tie does not
# have a certain change stat).

### rvl 

# d <- cstat_rvl1[[1]]
# for(i in 2:length(cstat_rvl1)){
#   d <- merge(d,cstat_rvl1[[i]],all.x=T,all.y = T)
# }
# 
# a <- gf_rvl[[1]]%v%"vertex.names"
# b <- paste(rep(a,each=length(a)),rep(a,times=length(a)),sep="_")
# b <- b[rep(a,each=length(a))!=rep(a,times=length(a))] # take the diagonal off
# length(b)
# 
# # Use this to properly order d.
# rownames(d) <- d$edge_id
# d <- d[b,]
# cstat_rvl1 <- d
# 
# ### jbs 
# 
# d <- cstat_jbs1[[1]]
# for(i in 2:length(cstat_jbs1)){
#   d <- merge(d,cstat_jbs1[[i]],all.x=T,all.y = T)
# }
# 
# a <- gf_jbs[[1]]%v%"vertex.names"
# b <- paste(rep(a,each=length(a)),rep(a,times=length(a)),sep="_")
# b <- b[rep(a,each=length(a))!=rep(a,times=length(a))] # take the diagonal off
# length(b)
# 
# # Use this to properly order d.
# rownames(d) <- d$edge_id
# d <- d[b,]
# cstat_jbs1 <- d
# 
# ### pm 
# 
# d <- cstat_pm1[[1]]
# for(i in 2:length(cstat_pm1)){
#   d <- merge(d,cstat_pm1[[i]],all.x=T,all.y = T)
# }
# 
# a <- gf_pm[[1]]%v%"vertex.names"
# b <- paste(rep(a,each=length(a)),rep(a,times=length(a)),sep="_")
# b <- b[rep(a,each=length(a))!=rep(a,times=length(a))] # take the diagonal off
# length(b)
# 
# # Use this to properly order d.
# rownames(d) <- d$edge_id
# d <- d[b,]
# cstat_pm1 <- d
# 
# ### dmz 
# 
# d <- cstat_dmz1[[1]]
# for(i in 2:length(cstat_dmz1)){
#   d <- merge(d,cstat_dmz1[[i]],all.x=T,all.y = T)
# }
# 
# a <- gf_dmz[[1]]%v%"vertex.names"
# b <- paste(rep(a,each=length(a)),rep(a,times=length(a)),sep="_")
# b <- b[rep(a,each=length(a))!=rep(a,times=length(a))] # take the diagonal off
# length(b)
# 
# # Use this to properly order d.
# rownames(d) <- d$edge_id
# d <- d[b,]
# cstat_dmz1 <- d

### 

# # And export those.
# write.csv(cstat_rvl1,file = "./data/7 - saved models/cstat_rvl1.csv",fileEncoding = "UTF-8")
# write.csv(cstat_jbs1,file = "./data/7 - saved models/cstat_jbs1.csv",fileEncoding = "UTF-8")
# write.csv(cstat_pm1,file = "./data/7 - saved models/cstat_pm1.csv",fileEncoding = "UTF-8")
# write.csv(cstat_dmz1,file = "./data/7 - saved models/cstat_dmz1.csv",fileEncoding = "UTF-8")

########## 9. Layouts for plots ##########

# The code for creating plots and saving the layout is in comments.
# Run it only once, export the layouts, and then load them for all firther analyses. This way
# the netyork plots from gplot will always place the nodes in the same spot.ß

# layout_gf_rvl <- list()
# layout_gf_rvl[[1]] <- gplot(gf_rvl[[1]])
# layout_gf_rvl[[2]] <- gplot(gf_rvl[[2]])
# layout_gf_rvl[[3]] <- gplot(gf_rvl[[3]])
# layout_gf_rvl[[4]] <- gplot(gf_rvl[[4]])
# layout_gf_rvl[[6]] <- gplot(gf_rvl[[5]])
# 
# gplot(gf_rvl[[1]],coord = layout_gf_rvl[[1]],displaylabels = T, label = tolower(gf_rvl[[1]]%v%"name"), 
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_rvl[[2]],coord = layout_gf_rvl[[2]],displaylabels = T, label = tolower(gf_rvl[[2]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_rvl[[3]],coord = layout_gf_rvl[[3]],displaylabels = T, label = tolower(gf_rvl[[3]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_rvl[[4]],coord = layout_gf_rvl[[4]],displaylabels = T, label = tolower(gf_rvl[[4]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# 
# save(layout_gf_rvl,file = "./data/5 - final data/layout/layout_gf_rvl.RData")

# layout_gf_jbs <- list()
# layout_gf_jbs[[1]] <- gplot(gf_jbs[[1]])
# layout_gf_jbs[[2]] <- gplot(gf_jbs[[2]])
# layout_gf_jbs[[3]] <- gplot(gf_jbs[[3]])
# layout_gf_jbs[[4]] <- gplot(gf_jbs[[4]])
# 
# gplot(gf_jbs[[1]],coord = layout_gf_jbs[[1]],displaylabels = T, label = tolower(gf_jbs[[1]]%v%"name"), 
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_jbs[[2]],coord = layout_gf_jbs[[2]],displaylabels = T, label = tolower(gf_jbs[[2]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_jbs[[3]],coord = layout_gf_jbs[[3]],displaylabels = T, label = tolower(gf_jbs[[3]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_jbs[[4]],coord = layout_gf_jbs[[4]],displaylabels = T, label = tolower(gf_jbs[[4]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# 
# save(layout_gf_jbs,file = "./data/5 - final data/layout/layout_gf_jbs.RData")

# layout_gf_pm <- list()
# layout_gf_pm[[1]] <- gplot(gf_pm[[1]])
# layout_gf_pm[[2]] <- gplot(gf_pm[[2]])
# layout_gf_pm[[3]] <- gplot(gf_pm[[3]])
# layout_gf_pm[[4]] <- gplot(gf_pm[[4]])
# 
# gplot(gf_pm[[1]],coord = layout_gf_pm[[1]],displaylabels = T, label = tolower(gf_pm[[1]]%v%"name"), 
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_pm[[2]],coord = layout_gf_pm[[2]],displaylabels = T, label = tolower(gf_pm[[2]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_pm[[3]],coord = layout_gf_pm[[3]],displaylabels = T, label = tolower(gf_pm[[3]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_pm[[4]],coord = layout_gf_pm[[4]],displaylabels = T, label = tolower(gf_pm[[4]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# 
# save(layout_gf_pm,file = "./data/5 - final data/layout/layout_gf_pm.RData")

# layout_gf_dmz <- list()
# layout_gf_dmz[[1]] <- gplot(gf_dmz[[1]])
# layout_gf_dmz[[2]] <- gplot(gf_dmz[[2]])
# layout_gf_dmz[[3]] <- gplot(gf_dmz[[3]])
# layout_gf_dmz[[4]] <- gplot(gf_dmz[[4]])
# 
# gplot(gf_dmz[[1]],coord = layout_gf_dmz[[1]],displaylabels = T, label = tolower(gf_dmz[[1]]%v%"name"), 
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_dmz[[2]],coord = layout_gf_dmz[[2]],displaylabels = T, label = tolower(gf_dmz[[2]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_dmz[[3]],coord = layout_gf_dmz[[3]],displaylabels = T, label = tolower(gf_dmz[[3]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# gplot(gf_dmz[[4]],coord = layout_gf_dmz[[4]],displaylabels = T, label = tolower(gf_dmz[[4]]%v%"name"),
#       label.pos = 3, edge.col="grey50",usearrows = F,vertex.col="blue") 
# 
# save(layout_gf_dmz,file = "./data/5 - final data/layout/layout_gf_dmz.RData")

load( "./data/5 - final data/layout/layout_gf_rvl.RData")
load( "./data/5 - final data/layout/layout_gf_jbs.RData")
load( "./data/5 - final data/layout/layout_gf_pm.RData")
load( "./data/5 - final data/layout/layout_gf_dmz.RData")

###

# Layout for top5 (with rm NA!!).
# NB : was based on the top5 plots tried out for rvl in the script "chapter 7".
# Haven t done the other schools yet.

# save(layout_top5_rvl,file = "./data/5 - final data/layout/layout_top5_rvl.RData")

load(file = "./data/5 - final data/layout/layout_top5_rvl.RData")


########## 10. geographic data ##########

# Pairwise walking distances between all students' home have been computed
# from OSM data on a separate script. Here they are simply loaded.
# Can be obtained from the author on demand.

add_rvl <- read.csv("./data/5 - final data/special/addresses_rvl.csv",fileEncoding="UTF-8",as.is=T)
add_jbs <- read.csv("./data/5 - final data/special/addresses_jbs.csv",fileEncoding="UTF-8",as.is=T)
add_pm <- read.csv("./data/5 - final data/special/addresses_pm.csv",fileEncoding="UTF-8",as.is=T)
add_dmz <- read.csv("./data/5 - final data/special/addresses_dmz.csv",fileEncoding="UTF-8",as.is=T)

load("./data/5 - final data/special/walk_rvl.RData")
load("./data/5 - final data/special/walk_jbs.RData")
load("./data/5 - final data/special/walk_pm.RData")
load("./data/5 - final data/special/walk_dmz.RData")

# re-size matrices per wave (using the valid_ids defined by gf,
# depending on whether or not imputed data has been used).

walk_rg_rvl <- list() # rg for "regular"
walk_rg_jbs <- list()
walk_rg_pm <- list()
walk_rg_dmz <- list()

walk_sq_rvl <- list()
walk_sq_jbs <- list()
walk_sq_pm <- list()
walk_sq_dmz <- list()

for(i in 1:length(gf_rvl)){
  walk_rg_rvl[[i]] <- walk_rvl[["regular"]][as.character(gf_rvl[[i]]%v%"vertex.names"),as.character(gf_rvl[[i]]%v%"vertex.names")]
  walk_rg_jbs[[i]] <- walk_jbs[["regular"]][as.character(gf_jbs[[i]]%v%"vertex.names"),as.character(gf_jbs[[i]]%v%"vertex.names")]
  walk_rg_pm[[i]] <- walk_pm[["regular"]][as.character(gf_pm[[i]]%v%"vertex.names"),as.character(gf_pm[[i]]%v%"vertex.names")]
  walk_rg_dmz[[i]] <- walk_dmz[["regular"]][as.character(gf_dmz[[i]]%v%"vertex.names"),as.character(gf_dmz[[i]]%v%"vertex.names")]
  
  walk_sq_rvl[[i]] <- walk_rvl[["squared"]][as.character(gf_rvl[[i]]%v%"vertex.names"),as.character(gf_rvl[[i]]%v%"vertex.names")]
  walk_sq_jbs[[i]] <- walk_jbs[["squared"]][as.character(gf_jbs[[i]]%v%"vertex.names"),as.character(gf_jbs[[i]]%v%"vertex.names")]
  walk_sq_pm[[i]] <- walk_pm[["squared"]][as.character(gf_pm[[i]]%v%"vertex.names"),as.character(gf_pm[[i]]%v%"vertex.names")]
  walk_sq_dmz[[i]] <- walk_dmz[["squared"]][as.character(gf_dmz[[i]]%v%"vertex.names"),as.character(gf_dmz[[i]]%v%"vertex.names")]
}

valid_rvl <- unique(c(as.character(gf_rvl[[1]]%v%"vertex.names"),as.character(gf_rvl[[3]]%v%"vertex.names"),
                      as.character(gf_rvl[[2]]%v%"vertex.names"),as.character(gf_rvl[[4]]%v%"vertex.names"),
                      as.character(gf_rvl[[5]]%v%"vertex.names")))
valid_jbs <- unique(c(as.character(gf_jbs[[1]]%v%"vertex.names"),as.character(gf_jbs[[3]]%v%"vertex.names"),
                      as.character(gf_jbs[[2]]%v%"vertex.names"),as.character(gf_jbs[[4]]%v%"vertex.names"),
                      as.character(gf_jbs[[5]]%v%"vertex.names")))
valid_pm <- unique(c(as.character(gf_pm[[1]]%v%"vertex.names"),as.character(gf_pm[[3]]%v%"vertex.names"),
                      as.character(gf_pm[[2]]%v%"vertex.names"),as.character(gf_pm[[4]]%v%"vertex.names"),
                      as.character(gf_pm[[5]]%v%"vertex.names")))
valid_dmz <- unique(c(as.character(gf_dmz[[1]]%v%"vertex.names"),as.character(gf_dmz[[3]]%v%"vertex.names"),
                      as.character(gf_dmz[[2]]%v%"vertex.names"),as.character(gf_dmz[[4]]%v%"vertex.names"),
                      as.character(gf_dmz[[5]]%v%"vertex.names")))


valid_rvl <- rownames(walk_rvl[["regular"]])[rownames(walk_rvl[["regular"]])%in%valid_rvl]
valid_jbs <- rownames(walk_jbs[["regular"]])[rownames(walk_jbs[["regular"]])%in%valid_jbs]
valid_pm <- rownames(walk_pm[["regular"]])[rownames(walk_pm[["regular"]])%in%valid_pm]
valid_dmz <- rownames(walk_dmz[["regular"]])[rownames(walk_dmz[["regular"]])%in%valid_dmz]

walk_rvl[["regular"]] <- walk_rvl[["regular"]][valid_rvl,valid_rvl]
walk_jbs[["regular"]] <- walk_jbs[["regular"]][valid_jbs,valid_jbs]
walk_pm[["regular"]] <- walk_pm[["regular"]][valid_pm,valid_pm]
walk_dmz[["regular"]] <- walk_dmz[["regular"]][valid_dmz,valid_dmz]

walk_rvl[["squared"]] <- walk_rvl[["squared"]][valid_rvl,valid_rvl]
walk_jbs[["squared"]] <- walk_jbs[["squared"]][valid_jbs,valid_jbs]
walk_pm[["squared"]] <- walk_pm[["squared"]][valid_pm,valid_pm]
walk_dmz[["squared"]] <- walk_dmz[["squared"]][valid_dmz,valid_dmz]

########## 11. Remove useless data ##########
if(keep_dat_gf==F){rm(gf_rvl,gf_jbs,gf_pm,gf_dmz)}
if(keep_dat_f==F){rm(f_rvl,f_jbs,f_pm,f_dmz,allf_rvl,allf_jbs,allf_pm,allf_dmz)}
if(keep_dat_dl==F){rm(dl_rvl,dl_jbs,dl_pm,dl_dmz)}
if(keep_dat_top5==F){rm(top5_rvl,top5_jbs,top5_pm,top5_dmz)}
if(keep_dat_popularities==F){rm(bully_rvl,bully_jbs,bully_pm,bully_dmz,
                                cool_rvl,cool_jbs,cool_pm,cool_dmz,
                                fear_rvl,fear_jbs,fear_pm,fear_dmz,
                                liked_rvl,liked_jbs,liked_pm,liked_dmz,
                                known_rvl,known_jbs,known_pm,known_dmz)}
if(keep_dat_culture==F){rm(culture)}
if(keep_dat_psycho==F){rm(psycho)}
if(keep_dat_ind_all==F){rm(ind1,ind2,ind3,ind4)}
if(keep_dat_ind_per_school==F){rm(rvl1,rvl2,rvl3,rvl4,rvl6,
                                  jbs1,jbs2,jbs3,jbs4,jbs6,
                                  pm1,pm2,pm3,pm4,pm6,
                                  dmz1,dmz2,dmz3,dmz4,dmz6)}
if(keep_osm_addresses==F){rm(add_rvl,add_jbs,add_pm,add_dmz)}
if(keep_lockdown==F){rm(lock_pers_rvl,lock_pers_jbs,lock_pers_pm,lock_pers_dmz,
                        lock_vocal_rvl,lock_vocal_jbs,lock_vocal_pm,lock_vocal_dmz,
                        lock_write_rvl,lock_write_jbs,lock_write_pm,lock_write_dmz,
                        lock_work_rvl,lock_work_jbs,lock_work_pm,lock_work_dmz)}

time2 <- Sys.time()
time2-time1
beep()