# Based on B. Nachtergaele's code.


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
        The potentially non-zero elements of each column, e.g., for each m, are computed solving a linear recursion (as described, e.g., in the paper by Straub.
    """

    tensor_dim = Integer((2*j1+1)*(2*j2+1))
    original_dim = Integer(2*j+1)
    V = []
    for spins in cartesian_product_iterator([spin_range(j1), spin_range(j2)]):
        V.append(
                map(lambda j_tmp:
                    safe_simplify(clebsch_gordan(j1, j2, j, spins[0], spins[1], j_tmp)), spin_range(j)
                    )
                )

    return matrix(V)


def safe_simplify(a):

    try:
        a_tmp = a.canonicalize_radical()
    except AttributeError:
        pass

    return a_tmp.simplify_full()

def spin_range(j):
    """
    Computes the array [-j, -j+1, -j+2,\\ldots, j-2, j-1, j]

    INPUT::

        - `j` --    (Half) Integer.

    OUTPUT::

        An `array` with entries [-j, -j+1, \\ldots, j-1, j]

    """

    return [i-j for i in range(0, Integer(2*j+1))]
