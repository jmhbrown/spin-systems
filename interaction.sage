load('intertwiner.sage')
load('utilities.sage')
import logging
logging.basicConfig(level=logging.DEBUG)

from sympy import SparseMatrix as sympySM


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

def interaction(proj, **kwargs):
    """
    Constructions the interaction.

    INPUT::

        - `proj` --     Matrix. Required. Projection onto the ground state space. Should be a 4^n by 4^n matrix,
                        with `n` being the number of sites affected by the projection.

        - `length` --   Integer. Optional, default: 4. The number of sites in the spin chain.
        - `site` --     Integer between 0 and `length`-1. The first site at which the interaction acts.

    OUTPUT::

        Matrix which acts non-trivially on `site`, `site`+1, ... , `site`+ `n`
    """
