load('intertwiner.sage')


def partial_trace(M, n=2):
    """
    Takes the partial trace of `M` over CC^n
    e.g. replaces each `n` by `n` block with its trace
    """

    if n.divides(M.ncols()) and n.divides(M.nrows()):
        p_trace = []
        M.subdivide(range(n, M.nrows()-n+1,n), range(n, M.ncols()-n+1,n))

        for i in range(0, M.nrows()-n+1,n):
            row = []
            for j in range(0,M.ncols()-n+1,n):
                row.append(M.subdivision(i/n,j/n).trace())
            p_trace.append(row)

        return matrix(p_trace)
    else:
        raise Exception

def pre_partial_trace(V1, V2, normalization_constant):
    """
    Computes an intermediate matrix used in the derivation of the density matrix.

    Normalization constant is generally either 1/2 or 1/3.
    """
    return normalization_constant*matrix.identity(4).tensor_product(V2)*V1*V1.conjugate_transpose()*matrix.identity(4).tensor_product(V2.conjugate_transpose())

def density_matrix():
    V2 = intertwiner(1/2, 3/2, ZZ(1))
    V1 = intertwiner(ZZ(1), 3/2, 1/2)

    M1 = partial_trace(pre_partial_trace(V1, V2, 1/3), 3)
    M2 = partial_trace(pre_partial_trace(V2, V1, 1/2), 2)

    return [M1, M2]


