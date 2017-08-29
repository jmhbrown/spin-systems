load('intertwiner.sage')
load('utilities.sage')
import logging
import sympy


def projection(vectors):
    """
    Creates the orthogonal projection onto the space spanned by vectors.

    INPUT::

        - `vectors` --  Required a list of vectors.

    OUTPUT::

        A matrix. The orthogonal projection onto the space spanned by
        the provided vectors.
    """


    vec_matrix = sympy.SparseMatrix(map(
        lambda col:
            list(col),

        vectors
    )).T

    logging.info("Constructing projection from %d by %d matrix of vectors" % vec_matrix.shape)

    tmp = (vec_matrix.conjugate()*vec_matrix).inv()

    return vec_matrix*tmp*vec_matrix.conjugate()
