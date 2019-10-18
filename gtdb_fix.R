require("castor")
treefile="bac120_r89.tree"
tree=castor::read_tree(file=treefile,interpret_quotes=TRUE)
cat(sprintf("Tree contains %d tips\n",length(tree$tip.label)))
tree$node.label=NULL

fixed_treefile="bac120_r89_fixed.tre"
castor::write_tree(tree,file=fixed_treefile)
