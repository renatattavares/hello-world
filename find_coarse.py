a=0
    b=0
    for nodes, edges in zip(self.nodes_neighbors.values(), self.edges_neighbors.values()):
        self.all_nodes_neighbors = rng.unite(self.all_nodes_neighbors,nodes)
        self.all_edges_neighbors = rng.unite(self.all_edges_neighbors,edges)
        a = a + 1

    for faces, volumes in zip(self.faces_neighbors.values(), self.volumes_neighbors.values()):
        self.all_faces_neighbors = rng.unite(self.all_faces_neighbors,nodes)
        self.all_volumes_neighbors = rng.unite(self.all_volumes_neighbors,edges)
        b = b + 1
