# Based on B. Nachtergaele's code.
load('utilities.sage')

import logging
FORMAT = '%(asctime)s %(levelname)s %(funcName)s : %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

def intertwiner(j, j1, j2):
    """
    Computes the SU(2) intertwining isometry

        V : D^(j) \\to D^(j1)\otimes D^(j2)

    TODO: Still needs to check whether D^(k) is contained in D^(j1)\otimes D^(j2), and should return warning and V=0 if not.

    INPUT::

        - `j` --    (Half) Integer. Spin of the original representation. This program assumes that |j1 - j2| \leq j \leq j1 + j2.
        - `j1` --   (Half) Integer. Spin of one of the target representations.
        - `j2` --   (Half) Integer. Spin of one of the target representations.

    OUTPUT::

        The SU(2) intertwining isometry, V. Each column of V contains the coefficients of a vector |j,m>, m=j,j-1,...,-j expressed in the tensor product basis |j1, m1> \otimes |j2, m2>. (where m1 = j1,j1-1,...,-j1 and m2 = j2,j2-1,...,-j2).
    """

    V = []
    logging.info("Building intertwiner V^{%s; %s, %s}" % (j, j1, j2))
    for spins in cartesian_product_iterator([spin_range(j1), spin_range(j2)]):
        #   If j, j1, j2 don't satisfy the triangle relations:
        #   | j1 - j2 | \leq j \leq j1 + j2
        # then we know that the Clebsch Gordan coefficients will be zero
        # and can skip actually calculating them.
        if abs(j1 - j2) > j or j > j1 + j2:
            V.append([0]*Integer(2*j+1))
        else:
            V.append(
                map(lambda j_tmp:
                    safe_simplify(clebsch_gordan(j1, j2, j, spins[0], spins[1], j_tmp)), spin_range(j)
                    )
                )

    #XXX:   Sage has some extra methods for dealing with dense symbolic matrices.
    #       not sure if any of them are worth the performance hit though.
    return matrix(V, sparse=True)
