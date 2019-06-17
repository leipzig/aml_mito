# local dir (Marijana)
# setwd("K:\\Dropbox\\CHOP\\MitoSeq\\phylotree\\")
# rCRS.df <- read.csv("rCRS.db.v16.csv", header = T, stringsAsFactors = F, check.names = F)
# rCRS.df[, 1] <- NULL
# nonregAncestor <- read.csv("non-regularAncestor.csv", header = T, stringsAsFactors = F, check.names = F)
# 
# List of functions
# 
# load.phylo()
# node.exists()
# node.variants()
# node.difference()
# node.lineage()
# node.ancestor()
#

# load.phylo
# Description: Import Phylotree v16 if not in workspace already
# Usage: load.phylo()
# Arguments:
#   - None
# Value:
#   - Phylotree v16 (rCRS.df) and a table of non-regular Ancestors (nonregAncestor) are imported into the workspace

load.phylo <- function()
{
  options(warn = -1)
  
  if(!exists("rCRS.df") | !exists("nonregAncestor"))
  {
    rCRSfile = "rawData/rCRS.db.v16.csv.gz" 
    nonregAncestorFile = "rawData/non-regularAncestor.csv" 
    if(!file.exists(rCRSfile) | !file.exists(nonregAncestorFile))
    {
      return("Phylotree and/or ancestor lookup file(s) not present in specified directory")
    }
    if(file.exists(rCRSfile) & file.exists(nonregAncestorFile))
    {
      rCRS.df <- read.csv(rCRSfile, header = T, stringsAsFactors = F, check.names = F)
      rCRS.df[, 1] <- NULL
      nonregAncestor <- read.csv(nonregAncestorFile, header = T, stringsAsFactors = F, check.names = F)
      assign('rCRS.df',  rCRS.df, envir = parent.frame())
      assign('nonregAncestor', nonregAncestor, envir = parent.frame())
    }
  }
}

# node.exists
# Description: Check whether node exists in Phylotree v16
# Usage: node.exists(nodeName)
# Arguments:
#   nodeName: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
# Value:
#   - logical(1) depending on whether node exists

node.exists <- function(nodeName = NULL)
{
  if(is.null(nodeName))
  {
    return("Error: A haplogroup name is required, e.g. node.exists('L3')")
  }
  
  ifelse(!(nodeName %in% colnames(rCRS.df)), return(F), return (T))
}

# node.variants
# Description: Retrieve list of variants from a node in Phylotree v16
# Usage: node.variants(nodeName, reference, descend)
# Arguments:
#   nodeName: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
# reference: character() - The reference sequence used (rCRS or RSRS)
# Value:
#   A list() of variant objects, each itself a list() consisting of:
#   - integer() The locus position
#   - character() The reference allele
#   - character() The variant allele

node.variants <- function(nodeName = NULL, reference = "rCRS")
{
  options(warn = -1)
  
  if(is.null(nodeName))
  {
    return("Error: A haplogroup name is required, e.g. node.variants('L3')")
  }
  
  if(node.exists(nodeName) == F)
  {
    return("Error: Haplogroup name not present in Phylotree v16")
  }
  
  # Return empty dataframe if nodeName H2a2a1 (rCRS reference)
  if(nodeName == "H2a2a1")
  {
    node <- data.frame(POS    = numeric(),
                       rCRS   = character(), 
                       H2a2a1 = character(), 
                      stringsAsFactors = FALSE) 
    return(node)
  }
  
  # Select variant
  node <- rCRS.df[which(is.na(rCRS.df[, nodeName]) == F), c("POS", "rCRS", nodeName)]
  row.names(node) <- NULL
    
  # Translate variant with regard to rCRS reference
  for(i in 1:nrow(node))
  {
    # Transition
    if(is.na(as.numeric(node[i, nodeName])) == F)
    {
      if(node[i, "rCRS"] == "A") { node[i, nodeName] <- "G"; next}
      if(node[i, "rCRS"] == "G") { node[i, nodeName] <- "A"; next}
      if(node[i, "rCRS"] == "C") { node[i, nodeName] <- "T"; next}
      if(node[i, "rCRS"] == "T") { node[i, nodeName] <- "C"; next}
    }
    
    # Deletion
    if(regexpr("d", node[i, nodeName])[[1]] > 0) 
    { 
      node[i, nodeName] <- "-"
      next
    }
    
    # Insertion
    if(regexpr("\\.", node[i, nodeName])[[1]] > 0)
    {
      # Known insertions (denoted with .1)
      if(substr(node[i, nodeName], regexpr("\\.", node[i, nodeName]) + 1, regexpr("\\.", node[i, nodeName]) + 1) == "1")
      {
        node[i, nodeName] <- paste(node[i, "rCRS"], substr(node[i, nodeName], regexpr("\\.", node[i, nodeName]) + 2, nchar(node[i, nodeName])), sep = "")
        next
      }
      # polynucleotide stretches of unknown length are noted with XC
      if(substr(node[i, nodeName], regexpr("\\.", node[i, nodeName]) + 1, regexpr("\\.", node[i, nodeName]) + 1) == "X")
      {
        node[i, nodeName] <- paste(node[i, "rCRS"], "XC", sep = "")
        next
      }  
    }
      
    # Transversion
    if(substr(node[i, nodeName], nchar(node[i, nodeName]), nchar(node[i, nodeName])) %in% c("A", "C", "G", "T"))
    {
      node[i, nodeName] <- substr(node[i, nodeName], nchar(node[i, nodeName]), nchar(node[i, nodeName]))
      next
    }
  } 
  return(node)
}

# node.difference
# Description: Retrieve the difference in variants in Phylotree v16
# Important: 
# If two node names are supplied, then their difference in variants is returned
# If one node name is supplied, then the difference between the node and it's direct ancestor is returned
# Usage: node.difference(nodeName1, nodeName2 (optional))
# Arguments:
#   nodeName1: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
#   nodeName2: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
# Value:
#   A list() of variant objects, each itself a list() consisting of:
#   - integer() The locus position
#   - character() The reference allele
#   - character() The variant allele

node.difference <- function(nodeName1 = NULL, nodeName2 = NULL)
{
  options(warn = -1)
  
  if(is.null(nodeName1))
  {
    return("Error: One or two haplogroup names are required, e.g. node.difference('L3') or node.difference('L3', 'D4b1b')")
  }
  
  # Two nodes (difference between two nodes)
  if(!is.null(nodeName2))
  {
    if(node.exists(nodeName1) == F | node.exists(nodeName2) == F)
    {
      return("Error: Haplogroup name(s) not present in Phylotree v16")
    }
    
    # get variants in node 1
    node1 = node.variants(nodeName1)
    colnames(node1)[3] <- "Variant"
    if(nrow(node1) != 0)
    {
      node1$Haplogroup   <- nodeName1
    }
    
    # get variants in node 2
    node2 = node.variants(nodeName2)
    colnames(node2)[3] <- "Variant"
    if(nrow(node1) != 0)
    {
      node2$Haplogroup <- nodeName2
    }
  }
  
  # One node (difference node and ancestor)
  if(is.null(nodeName2))
  {
    if(node.exists(nodeName1) == F)
    {
      return("Error: Haplogroup name(s) not present in Phylotree v16")
    }
    
    if(nodeName1 == "mt-MRCA")
    {
      return("Error: mt-MRCA has no ancestor")
    }
    
    # get variants in node 1
    node1 = node.variants(nodeName1)
    colnames(node1)[3] <- "Variant"
    if(nrow(node1) != 0)
    {
      node1$Haplogroup <- nodeName1
    }
    
    # get variants in ancestor
    node2 = node.variants(node.ancestor(nodeName1))
    colnames(node2)[3] <- "Variant"
    node2$Haplogroup <- node.ancestor(nodeName1)
  }
  
  # combine variants from both nodes (non-duplicated)
  nodeDiff <- rbind(node1, node2)[!duplicated(rbind(node1, node2)[1:3]), ]
  
  # return variants that exist in one and not in other
  node1Var     <- do.call("paste", node1[1:3])
  node2Var     <- do.call("paste", node2[1:3])
  nodeDiff     <- rbind(nodeDiff[!node1Var %in% node2Var, ], nodeDiff[!node2Var %in% node1Var, ])
  return(nodeDiff)
}

# node.lineage
# Description: Retrieve entire ancestral lineage of a node from Phylotree
# Usage: node.lineage(nodeName)
# Arguments:
#   nodeName: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
# Value:
#   A list() containing all ancestors

node.lineage <- function(nodeName = NULL)
{
  options(warn = -1)
  
  if(is.null(nodeName))
  {
    return("Error: A haplogroup name is required, e.g. node.lineage('L3')")
  }
  
  if(node.exists(nodeName) == F)
  {
    return("Error: Haplogroup name not present in Phylotree v16")
  }
  
  # initialize ancestral tree (list)
  ancestralLineage <- nodeName
  
  # Obtain next ancestor until the root of the tree (mt-MRCA)
  while(nodeName != "mt-MRCA")
  {
    ancestralLineage  <- c(node.ancestor(nodeName), ancestralLineage)
    nodeName <- node.ancestor(nodeName)
  }  
  return(ancestralLineage)
}

# node.ancestor
# Description: Retrieve the parent of a node from Phylotree
# Usage: node.ancestor(nodeName, reference)
# Arguments:
#   nodeName: character() - Any node of Phylotree (e.g. H, H2a, H2a21a, etc...)
#   reference: character() - The reference sequence used (rCRS or RSRS)
# Value:
#   A character() representing the node's ancestor

node.ancestor <- function(nodeName = NULL, reference = "rCRS")
{
  options(warn = -1)
  
  if(is.null(nodeName))
  {
    return("Error: A haplogroup name is required, e.g. node.ancestor('L3')")
  }
  
  if(node.exists(nodeName) == F)
  {
    return("Error: Haplogroup name not present in Phylotree v16")
  }
  
  if(nodeName == "mt-MRCA")
  {
    return("Error: mt-MRCA has no ancestor")
  }
  
  # Initialize parent
  ancestorName = NA
  
  # Deal with dashes
  if(gregexpr("-", nodeName)[[1]] > 0)
  {
    nodeName <- ifelse(gregexpr("[0-9]", substr(nodeName, regexpr("-", nodeName)[[1]] - 1, regexpr("-", nodeName)[[1]] - 1))[[1]] > 0,
           paste(substr(nodeName, 1, regexpr("-", nodeName)[[1]] - 1), "a", sep = ""), 
           paste(substr(nodeName, 1, regexpr("-", nodeName)[[1]] - 1), "1", sep = ""))
  }
  
  # Deal with node that have a non-regularry named ancestor (look-up)
  for(i in 1:nrow(nonregAncestor))
  {
    if(nonregAncestor$Node[i] == nodeName)
    {
      ancestorName <- nonregAncestor$Parent[i]
      next
    }
  }
  
  # Regular node-names
  if(is.na(ancestorName) == T)
  {
    seq <- c()
    for (i in 1:nchar(nodeName))
    {
      if(gregexpr("[0-9]", substr(nodeName, i, i))[[1]] > 0)
      {
        seq <- c(seq, "N")
      }
      if(gregexpr("[a-zA-Z]", substr(nodeName, i, i))[[1]] > 0)
      {
        seq <- c(seq, "C")
      }
      if(gregexpr("'", substr(nodeName, i, i))[[1]] > 0)
      {
        seq <- c(seq, "A")
      }
    }
    # Position of Numerics
    num.start <- gregexpr("[0-9]", nodeName)
    if(length(num.start[[1]]) >= 2)
    {
      # Remove Double-Digit Position
      for(i in length(num.start[[1]]):2)
      {
        if(num.start[[1]][i] == num.start[[1]][i - 1] + 1)
        {
          num.start[[1]] <- num.start[[1]][-i]
        }
      }
    }
    
    # Position of Characters
    char.start <- gregexpr("[a-zA-Z]", nodeName)
    if(char.start[[1]] >= 2)
    {
      for(i in length(char.start[[1]]):2)
      {
        if(char.start[[1]][i] == char.start[[1]][i - 1] + 1)
        {
          char.start[[1]] <- char.start[[1]][-i]
        }
      }
    }
    
    # Combine Numeric and Character Positions
    char.num <- as.data.frame(cbind(char.start[[1]], num.start[[1]]))
    
    # Remove last element if vectors are of unequal length 
    if(nrow(char.num) >= 2)
    {  
      if(char.num[nrow(char.num), 2] <= char.num[nrow(char.num) - 1, 2])
      {
        char.num[nrow(char.num), 2] <- NA
      }
    }
    colnames(char.num)[1] <- "LETTER"
    colnames(char.num)[2] <- "NUMBER"
    
    # Initialize the tree
    tree <- c()
    
    # Node name without apestrophe
    if(!("A" %in% seq))
    {
      for(i in 1:nrow(char.num))
      {
        # Fill in the tree (until last branch)
        if(i != nrow(char.num))
        {
          tree <- c(tree, substr(nodeName, char.num[i, 1], char.num[i, 2] - 1))
          tree <- c(tree, substr(nodeName, char.num[i, 2], char.num[i + 1, 1] - 1))
        }
        # Last branch
        if(i == nrow(char.num))
        {
          # Last element is letter
          if(is.na(char.num[i, 2]) == T)
          {
            tree <- c(tree, substr(nodeName, char.num[i, 1], nchar(nodeName)))
          }
          # Last element is number
          if(is.na(char.num[i, 2]) == F)
          {
            tree <- c(tree, substr(nodeName, char.num[i, 1], char.num[i, 2] - 1))
            tree <- c(tree, substr(nodeName, char.num[i, 2], nchar(nodeName)))
          }
        }
      }
    }
    ancestorName <- paste(tree[-length(tree)], collapse = '')
    
    # Node contains apestrophe 
    if("A" %in% seq)
    {
      # M73'79 should become M
      if(seq[grep("A", seq)[1] - 1] != seq[grep("A", seq)[1] - 2])
      {
        ancestorName <- substr(nodeName, 1, grep("A", seq)[1] - 2)
      }
      if(seq[grep("A", seq)[1] - 1] == seq[grep("A", seq)[1] - 2])
      {
        ancestorName <- substr(nodeName, 1, grep("A", seq)[1] - 3)
      }
    }
  }
  return(ancestorName)
}


#######
####### 
####### 

# refSEQ <- as.data.frame(strsplit(readChar('rCRS.txt', file.info('rCRS.txt')$size), ""))
# colnames(refSEQ)[1] <- "rCRS"
# refSEQ$rCRS <- as.character(refSEQ$rCRS)
# 
# rCRS.df <- as.data.frame(matrix(seq(16569)), nrow = 16569)
# colnames(rCRS.df)[1] <- "POS"
# rCRS.df <- cbind(rCRS.df, refSEQ)
# 
# rCRS <- read.csv("rCRS-build16.csv", F, fill = T, stringsAsFactors = F)
# rCRS[rCRS == ""] <- NA
# 
# # FILL OUT MATRIX
# for(i in 1:nrow(rCRS))
# {
#   rCRS.df[, i + 2] <- NA
#   colnames(rCRS.df)[i + 2] <- rCRS[i, 1]
#   for(j in 2:ncol(rCRS))
#   {
#     # TRANSITION
#     if(is.na(as.numeric(rCRS[i, j])) == F)
#     {
#       rCRS.df[rCRS[i, j], i + 2] <- rCRS[i, j]
#       next
#     }
#     # DELETIONS: SNPs that contain the letter d, e.g. 531d, or 527-563d
#     if(regexpr("d", rCRS[i, j])[[1]] > 0 & is.na(rCRS[i, j]) == F)
#     {
#       # Region deleted
#       if(regexpr("-", rCRS[i, j])[[1]] > 0)
#       {
#         for(k in as.numeric(unlist(strsplit(substr(rCRS[i, j], 1, nchar(rCRS[i, j])), "-"))[[1]]):as.numeric(unlist(strsplit(substr(rCRS[i, j], 1, nchar(rCRS[i, j]) - 1), "-"))[[2]]))
#         {
#           rCRS.df[k, i + 2] <- paste(k, "d", sep = "")
#         }
#         next
#       }
#       # Single Variant Deletion
#       if(regexpr("-", rCRS[i, j])[[1]] == -1)
#       {
#         if(is.na(as.numeric(substr(rCRS[i, j], 1, 1))) == F)
#         {
#           rCRS.df[as.numeric(substr(rCRS[i, j], 1, nchar(rCRS[i, j]) - 1)), i + 2] <- rCRS[i, j]
#           next
#         }
#         # There is one exception, C123d (starts with character)
#         if(is.na(as.numeric(substr(rCRS[i, j], 1, 1))) == F)
#         {
#           rCRS.df[as.numeric(substr(rCRS[i, j], 2, nchar(rCRS[i, j]) - 1)), i + 2] <- rCRS[i, j]
#           next
#         }
#       }
#     }
#     # INSERTIONS: SNPs that contain a dot, e.g. 123.1XT
#     if(regexpr("\\.", rCRS[i, j])[[1]] > 0 & is.na(regexpr("\\.", rCRS[i, j])[[1]]) == F)
#     {
#       rCRS.df[as.numeric(substr(rCRS[i, j], 1, regexpr("\\.", rCRS[i, j])[[1]] - 1)), i + 2] <- rCRS[i, j]
#       next
#     }
#     # TRANSVERSIONS: SNPs that end in the letters, A, C, G, or T
#     if(substr(rCRS[i, j], nchar(rCRS[i, j]), nchar(rCRS[i, j])) %in% c("A", "C", "G", "T") & is.na(rCRS[i, j]) == F)
#     {
#       rCRS.df[as.numeric(substr(rCRS[i, j], 1, nchar(rCRS[i, j]) - 1)), i + 2] <- rCRS[i, j]
#       next
#     }
#     # CONTROL: Empty fields
#     if(is.na(rCRS[i, j]) == T)
#     {
#       next
#     }
#   }
# }
# write.csv(rCRS.df, "rCRS.db.v16.csv")

#######
####### 
####### 
