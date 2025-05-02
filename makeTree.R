# retrieve sequences, align, infer tree
require(refseqR)
require(msa)
require(seqinr)
require(phangorn)
makeTree <- function(seq_list) {
  seqs <- refseq_AAseq(seq_list)
  aln <- msa(seqs)
  dist <- dist.ml(as.phyDat(msaConvert(aln, type = "seqinr::alignment") , type="AA"))
  tree <- bionj(dist)
  plot(tree)
  return(list(tree=tree, aln=aln, dist=dist, seqs=seqs))
}