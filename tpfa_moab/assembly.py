class Assembly(): # Montagem da matriz de coeficientes do sistema.
    def __init__(num_elements, M):
        print("Assembly")
        start = time.time()
        coef = lil_matrix((num_elements, num_elements), dtype=np.float_)

        for i in range(num_elements):
            adjacencies = M.volumes.bridge_adjacencies(i, 2, 3)
            length = np.shape(adjacencies)
            for j in range(length[1]):
                id = np.array([adjacencies[0][j]],  dtype= np.int)
                coef[i,id] = equiv_perm(M.permeability[i], M.permeability[id])/centroid_dist(M.volumes.center[i], M.volumes.center[id])
        coef[i,i] = (-1)*coef[i].sum()
        end = time.time()
        print("This step lasted {0}s".format(end-start))
        return coef

    def equiv_perm(k1, k2):
        return (2*k1*k2)/(k1 + k2)

    def centroid_dist(c1, c2):
        return ((c1-c2)**2).sum()
