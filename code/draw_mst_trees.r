#' Code for drawing igraph trees

library(igraph)
source("code/read_abundance.R")

#' Read in the data
y = rowSums(read_abundance())
y.red = rowSums(read.csv("no_error_abundance.csv", row.names=1))
d = as.matrix(round(read.csv("sero_dist15.csv", row.names=1) * 284))

#' Function to scale the abundances by so that the node size isn't
#' too giant
scale_abund <- function(x) {
  x^(1/3)/10
}
set.seed(2)
d.g = graph_from_adjacency_matrix(d, mode="undirected", weighted=TRUE)
d.mst = mst(d.g)
d.mst = set_vertex_attr(d.mst, "size", value=scale_abund(y))
d.mst = set_edge_attr(d.mst, "weight", value=edge_attr(d.mst, "weight"))

# label the big ones
V(d.mst)$label = substring(V(d.mst)$name, pmax(nchar(V(d.mst)$name) - 4, 0))
V(d.mst)$label[V(d.mst)$size < 1.5] = NA

layout = layout_nicely(d.mst)
pdf("figures/tree_full.pdf", width=20, height=20)
plot(d.mst, vertex.label.family="sans", vertex.label.cex=0.3, vertex.label.color = "black", layout=layout)
dev.off()

wch = names(y.red)

indices_to_remove <- function(tree, vertices_to_remain) {
  vertex_names = V(tree)$name
  dont_want = setdiff(vertex_names, vertices_to_remain)
  match(dont_want, vertex_names)
}

highlight_removed <- function(tree, vertices_to_remain, col = "orange", highlight="red") {
  ind_dont_want = indices_to_remove(tree, vertices_to_remain)

  # right, now highlight these
  col = rep(col, length(V(tree)))
  vsize = V(tree)$size
  vsize[ind_dont_want] = pmax(1, vsize[ind_dont_want])
  col[ind_dont_want] = "red"
  tree.high = set_vertex_attr(tree, "size", value=vsize)
  tree.high = set_vertex_attr(tree.high, "color", value=col)
  tree.high
}

reduce_tree <- function(tree, vertices_to_remain) {

  ind_dont_want = indices_to_remove(tree, vertices_to_remain)

  # the basic idea is we iterate removing one or more nodes at a time
  # until they're all gone.
  # If at any time we have leaves, we remove these first as they're easy.
  # If we have no leaves, we remove the node with smallest degree by replacing
  # the subgraph formed by the complete graph on it's neighbours, weighted by
  # shortest path length, with it's MST.
  to_remove = V(tree)$name[ind_dont_want]
  tree.red = tree
  while(length(to_remove) > 0) {
    degree = unlist(lapply(to_remove, function(x, g) { length(incident(g, x)) }, tree.red))
    remove_now = which(degree == 1)
    if (length(remove_now) > 0) { # we have leaves, so just remove them
      tree.red = delete_vertices(tree.red, to_remove[remove_now])
      to_remove = to_remove[-remove_now]
    } else {
      # no leaves, so remove the one of lowest degree
      remove_now = which.min(degree)

      # construct a complete graph from remove_now's neighbours
      v = neighbors(tree.red, to_remove[remove_now])
      sub = make_full_graph(length(v))
      # add vertex names
      V(sub)$name = names(v)
      V(sub)$label = names(v)
      # now add edge weights. This is a bit trickier...
      dis_sub = distances(tree.red, v = v, to = v)
      # fill in the weights
      E(sub)$weight = NA
      for (i in 1:(length(v)-1)) {
        for (j in (i+1):length(v))
          E(sub, P = names(v)[c(i,j)])$weight = dis_sub[i,j]
      }
      # ok, now we want to replace this graph with it's MST
      sub.mst = mst(sub)
      V(tree)$name
      # now add the edges from that MST subtree to our tree and remove the original vertex
      tree.red = add_edges(tree.red, t(ends(sub.mst, E(sub.mst))), attr=list(weight = E(sub.mst)$weight))
      tree.red = delete_vertices(tree.red, to_remove[remove_now])
      to_remove = to_remove[-remove_now]
    }
  }
  tree.red
}

d.mst.high <- highlight_removed(d.mst, names(y.red))
pdf("figures/tree_full_highlight.pdf", width=20, height=20)
plot(d.mst.high, vertex.label.family="sans", vertex.label.cex=0.3, vertex.label.color = "black", layout=layout)
dev.off()

d.mst.red = reduce_tree(d.mst, names(y.red))
pdf("figures/tree_simplified.pdf", width=20, height=20)
layout.red = layout[match(V(d.mst.red)$name, V(d.mst)$name),]
scale.y.red = scale_abund(y.red)
V(d.mst.red)[names(y.red)]$size <- scale.y.red
plot(d.mst.red, vertex.label.family="sans", vertex.label.cex=0.3, vertex.label.color = "black", layout=layout.red)
dev.off()

#' read in Patricks dataset
y.patrick = rowSums(read.csv("temp/patrick_cdhit_abundance.csv", row.names=1))

d.mst.patrick.high <- highlight_removed(d.mst, names(y.patrick))
pdf("figures/tree_full_highlight_patrick.pdf", width=20, height=20)
plot(d.mst.patrick.high, vertex.label.family="sans", vertex.label.cex=0.3, vertex.label.color = "black", layout=layout)
dev.off()

d.mst.patrick.red = reduce_tree(d.mst, names(y.patrick))
pdf("figures/tree_simplified_patrick.pdf", width=20, height=20)
layout.red.patrick = layout[match(V(d.mst.patrick.red)$name, V(d.mst)$name),]
scale.y.patrick.red = scale_abund(y.patrick)
V(d.mst.patrick.red)[names(y.patrick)]$size <- scale.y.patrick.red
plot(d.mst.patrick.red, vertex.label.family="sans", vertex.label.cex=0.3, vertex.label.color = "black", layout=layout.red.patrick)
dev.off()
