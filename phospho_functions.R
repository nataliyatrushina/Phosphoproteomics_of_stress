# Phosphoproteomics custom functions

# Important: sometimes operations work on matrices only then:
# dat_proc_out$Uniprot_accession_matrix <- as.matrix(dat_proc_out[,"Uniprot_accession"])

# Example usage:
# Get protein sequence for a given Uniprot ID
# get_protein_sequence("O70593")
# Get list of phosphomarks for a given peptide
# get_list_phosphomarks("TSSRES(+79.97)LNQVDLDC(+57.02)AVATFDPPS(+79.97)")
# Get list of phosphosites for a given peptide and start position
# get_list_phosphosites("TSSRES(+79.97)LNQVDLDC(+57.02)AVATFDPPS(+79.97)", 10)
# Get list of residues for a given protein sequence and phosphomarks
# get_list_residues("TSSRESLNQVDLDCAVATFDPPSDMESEAEDAPWSSDSLSR", "6, 20")
# Paste phosphosites for given residues and positions
# paste_phosphosites("S, T", "10, 15")

get_protein_sequence <- function(ac) {
  if(ac != ""){
    #ac = "A0A0G2JSW2"
    print(ac)
    requestURL <- paste("https://www.ebi.ac.uk/proteins/api/proteins/", as.character(ac), sep = "")
    r <- GET(requestURL, accept("application/xml"))
    
    try.test<-try(stop_for_status(r), silent = TRUE) 
    if (class(try.test) == "try-error"){
      return("")
    }
    else {
      xml <- read_xml(content(r,as = "raw"))
      return((as_list(xml))$entry$sequence[[1]])
    }
  } else {
    return("")
  }
} #get_protein_sequence("O70593")

#counter for skipping (+79.97) marks in phosphopeptide
counter <- list(-1,-9,-17,-25,-33) #because maximum allowed PTMs in PEAKS is 5

get_list_phosphomarks <- function(phosphopept) {
  return(toString(unlist(lapply(seq_along(as.list(unlist(gregexpr('\\(', phosphopept)))),function(i)
    unlist(as.list(unlist(gregexpr('\\(', phosphopept)))[i])+unlist(counter[i])))))
}

get_list_phosphosites <- function(phosphopept, start) {
  if (start != "") {
    return(toString(unlist(lapply(seq_along(as.list(unlist(gregexpr('\\(', phosphopept)))),function(i)
      unlist(as.list(unlist(gregexpr('\\(', phosphopept)))[i])+unlist(counter[i]))) + start - 1))
  }
  else {
    return("")
  }
}

get_list_residues <- function(pept,marks) {
  return(letter(pept, c(as.numeric(unlist(as.list(strsplit(marks, ", ")[[1]]))))))
}

paste_phosphosites <- function(res,site) {
  result <- unlist(lapply(seq_along(unlist(as.list(unlist(strsplit(res,", "))))),function(i)
    paste(unlist(as.list(unlist(strsplit(res,", ")))[i]),unlist(as.list(unlist(strsplit(site,", ")))[i]), sep = "")))
  return(unlist(result))
}
