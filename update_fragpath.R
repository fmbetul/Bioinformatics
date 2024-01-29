# Update the fragments file path

new_path = ".../fragments/GSE228120_atacfragments.tsv.gz"

frags <- Fragments(SeuFile)  # get list of fragment objects
Fragments(SeuFile) <- NULL  # remove fragment information from assay

frags[[1]] <- UpdatePath(frags[[1]], new.path = new_path) # update path

Fragments(SeuFile) <- frags # assign updated list back to the object

# check the result
SeuFile@assays[["ATAC"]]@fragments[[1]]@path