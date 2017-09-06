load('intertwiner.sage')
load('utilities.sage')
import logging
logging.basicConfig(level=logging.DEBUG)

from sympy import SparseMatrix as sympySM
from sympy import Matrix as m
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

def interaction(proj, **kwargs):
    """
    Constructs the interaction.

    INPUT::

        - `proj` --     Sympy Matrix. Required. Projection onto the ground state space. Should be a 4^n by 4^n matrix,
                        with `n` being the number of sites affected by the projection.

        - `length` --   Integer. Optional, default: 4. The number of sites in the spin chain.
        - `site` --     Integer between 0 and `length`-1. The first site at which the interaction acts.

    OUTPUT::

        Matrix which acts non-trivially on `site`, `site`+1, ... , `site`+ `n`
    """

    length = parse_kwargs(kwargs, 'length', 4)
    site = parse_kwargs(kwargs, 'site', 0)

    interaction_length = Integer(Integer(proj.cols).log(4))

    if length - site - interaction_length < 0:
        raise Exception("Interaction can't start at site %d!" % site)
    else:
        return tensor(tensor(m.eye(4**site), m.eye(proj.cols) - m(proj)), m.eye(4**(length - site - interaction_length)))

def hamiltonian(proj, **kwargs):
    """
    Constructs the hamiltonian.

    INPUT::

    - `proj` --     Sympy Matrix. Required. Projection onto the ground state space.
    - `length` --   Integer. Optional, default: 4. The number of sites in the spin chain.

    OUTPUT::

        Matrix representation of Hamlitonian.
    """

    interaction_length = Integer(Integer(proj.cols).log(4))

    length = parse_kwargs(kwargs, 'length', 4)

    hamiltonian = m.zeros(4**length)

    for site in range(0, length-interaction_length+1):
        logging.debug("Constructing interaction on sites %d to %d (of %d)" % ( site, site+interaction_length, length))
        hamiltonian += interaction(proj, length=length, site=site)

    return hamiltonian




