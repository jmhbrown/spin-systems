load('intertwiner.sage')
load('utilities.sage')
import logging
FORMAT = '%(asctime)s %(levelname)s %(funcName)s : %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)

from sympy import SparseMatrix as sympySM
from sympy import Matrix as m
from sympy import diag, floor
from sympy.physics.quantum import TensorProduct as tensor


def projection(vectors):
    """
    Creates the orthogonal projection onto the space spanned by `vectors`.

    INPUT::

        - `vectors` --  Required a list of vectors.

    OUTPUT::

        A matrix. The orthogonal projection onto the space spanned by
        the provided vectors. If V has the provided vectors as columns,
        then this computes V (V* V)^{-1} V*
    """


    vec_matrix = sympySM(map(
        lambda vec: list(vec), vectors
    )).T

    logging.info("Constructing projection from %d by %d matrix of vectors" % vec_matrix.shape)

    tmp = (vec_matrix.C.T*vec_matrix).inv()

    return vec_matrix*tmp*vec_matrix.C.T

def hamiltonian(proj, **kwargs):
    """
    Constructs the hamiltonian.

    INPUT::

    - `proj` --     Sage Matrix. Required. Projection onto the ground state space.
    - `length` --   Integer. Optional, default: 4. The number of sites in the spin chain.

    OUTPUT::

        Matrix representation of Hamlitonian.
    """

    d = 4   # 2*spin + 1 - dimension of the state space at each site
    interaction_length = Integer(Integer(proj.ncols()).log(d))
    length = parse_kwargs(kwargs, 'length', d)

    # initialize a nested array.
    h_array = [ e[:] for e in [[0]*(d**length)]*(d**length) ]

    logging.info("Computing entries for %d by %d hamiltonian!" % (len(h_array), len(h_array[0])))
    for site in range(0, length-interaction_length+1):
        logging.info("Computing contribution from site %d" % site)
        for i,j in cartesian_product_iterator([range(0,d**length), range(0,d**length)]):
            # This is identity(d**site) \otimes proj \otimes identity(d**(site+interaction_length))
            i0 = i % d**site
            j0 = j % d**site
            i1 = i % d**(site + interaction_length)
            j1 = j % d**(site + interaction_length)

            if (i-i1 == j-j1)*(i0 == j0):
                try:
                    h_array[i][j] += proj[(i1-i0) % proj.nrows()][(j1-j0) % proj.ncols()]
                except IndexError:
                    logging.error("Index error at h_array[%d][%d], proj[%d][%d]" %(i,j, (i1-i0) % proj.nrows(), (j1-j0) % proj.ncols()))
                    raise

    logging.info("Converting hamiltonian from array to matrix!")

    hamiltonian = sympySM(h_array)

    return hamiltonian
