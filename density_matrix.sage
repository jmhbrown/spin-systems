load('intertwiner.sage')
load('utilities.sage')
import logging

logging.basicConfig(level=logging.INFO)

def pre_partial_trace(V1, V2, normalization_constant, k=2):
    """
    Computes an intermediate matrix used in the derivation of the density matrix.

    Normalization constant is generally either 1/2 or 1/3.
    """
    tmp = normalization_constant*V1*V1.conjugate_transpose()
    logging.info("Building pre-partial trace version of the density matrix.")

    for i in range(2,k+1):
        logging.info("Starting iteration %d / %d of pre-partial trace matrix build." % (i-1, k-1))
        logging.info("Rank: %d" % tmp.rank())
        idty = matrix.identity(4**(i-1), sparse=True)
        if i % 2 == 0:
            tmp = idty.tensor_product(V2)*tmp*idty.tensor_product(V2.conjugate_transpose())
        if i % 2 == 1:
            tmp = idty.tensor_product(V1)*tmp*idty.tensor_product(V1.conjugate_transpose())

    return tmp

def density_matrix(k=2):
    """
    Builds the density matrix associated with the (type B) dimerized
    ground states. See Nachtergaele '96 for more information.

    INPUT::

        - `k` --    Integer. Optional, default: 2. Number of lattice sites.
                    Note that the density matrix grows as 4^k.

    OUTPUT::

        Array<Matrix> Density matrices associated with the two (type B) dimerized ground states.
        The two associated states differ only by a shift.
    """
    V2 = intertwiner(1/2, 3/2, Integer(1))
    V1 = intertwiner(Integer(1), 3/2, 1/2)

    if k % 2 == 0:
        M1 = partial_trace(pre_partial_trace(V1, V2, 1/3, k), 3)
        M2 = partial_trace(pre_partial_trace(V2, V1, 1/2, k), 2)
    if k % 2 == 1:
        M1 = partial_trace(pre_partial_trace(V1, V2, 1/2, k), 2)
        M2 = partial_trace(pre_partial_trace(V2, V1, 1/3, k), 3)

    return [M1, M2]