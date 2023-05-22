
.geneSetSignatures=fromJSON(read_json("geneSetSignatures.json",simplifyVector=T))

usethis::use_data(.geneSetSignatures,internal=T,overwrite=T)