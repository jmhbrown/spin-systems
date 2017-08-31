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
        the provided vectors.
    """


    vec_matrix = sympySM(map(
        lambda vec: list(vec), vectors
    )).T

    logging.info("Constructing projection from %d by %d matrix of vectors" % vec_matrix.shape)

    tmp = (vec_matrix.conjugate()*vec_matrix).inv()

    return vec_matrix*tmp*vec_matrix.conjugate()
