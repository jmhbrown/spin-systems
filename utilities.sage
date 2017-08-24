from sage.symbolic.expression_conversions import algebraic

import logging
logging.basicConfig(level=logging.INFO)


def partial_trace(M, n=2):
    """
    Takes the partial trace of `M` over CC^n
    e.g. replaces each `n` by `n` block with its trace
    """

    if n.divides(M.ncols()) and n.divides(M.nrows()):
        logging.info("Taking partial trace over CC^%d" % n)
        p_trace = []
        M.subdivide(range(n, M.nrows()-n+1,n), range(n, M.ncols()-n+1,n))

        for i in range(0, M.nrows()-n+1,n):
            row = []
            for j in range(0,M.ncols()-n+1,n):
                row.append(M.subdivision(i/n,j/n).trace())
                logging.debug("Tracing over the subdivision starting at (%d, %d)." % ( i, j) )
            p_trace.append(row)

        logging.info("Done computing partial trace. Converting result to matrix form.")
        mat = matrix(p_trace, sparse=True)
        logging.info("Post-trace rank: %d" % mat.rank())
        return mat
    else:
        raise Exception("Can't take partial trace over CC^%d of %d by %d matrix!" % (n, M.nrows(), M.ncols()))

def tensor_exponential(M, n):
    """
    Computes M tensored with itself n times.

    INPUT::

        - `M` --    Matrix.
        - `n` --    Integer. Tensor exponent.

    OUTPUT::

        A new matrix, M \otimes M \otimes \ldots \otimes M.

    """
    tmp = M

    for i in range(1,n):
        tmp = tmp.tensor_product(M)

    return tmp

def symbolic_eigenvectors(matrix):
    """
    Finds the left eigenvectors for a matrix with Symbolic entries.

    Works with either sparse or dense matrices.

    INPUT::

        - `matrix` --   Matrix. Must be over a symbolic ring.

    OUTPUT::

        A nested array
            [
                (first eigenvalue, [vector, vector,..., vector]),
                (second eigenvalue, [vector, vector,..., vector]),
                ...
                (last eigenvalue, [vector, vector,..., vector])
            ]
    """

    if matrix.base_ring() != SR:
        raise Exception("Matrix must be over a Symbolic Ring! Base ring is: %s" % matrix.base_ring())

    logging.info("Finding eigenvalues/vectors for %d by %d matrix." % (matrix.nrows(), matrix.ncols()))
    ev_array = matrix._maxima_(maxima).eigenvectors()

    # Fix the formating:
    #   1. convert elements to symbolic expressions
    #   2. convert funky maximal-arrays into actual vectors
    vectors_array = []
    for vec_list in ev_array[1]:
        vectors_array.append( map(
            lambda vec: vector( map(
                lambda elem: SR(elem),
                vec
            )),
            vec_list
        ))

    return zip( ev_array[0][0], vectors_array )



def safe_simplify(a):
    """
    A wrapper for .canonicalize_radical() which swallows AttributeError.
    """

    try:
        a_tmp = a.canonicalize_radical()
    except AttributeError:
        return a

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

def to_field_element(number, field=QQ):
    """
    Converts a (symbolic) element or expression to a field element.

    INPUT::

        - `number` --   A symbolic expression.
        - `field`  --   Optional. Field over which to define element.
                        Note that we automatically adjoin `number` to
                        the provided field if it isn't already there.

    OUTPUT::

        A version of `number` defined over a field.

    """

    if number.__class__ is sage.symbolic.expression.Expression:
        # Check if the provided field is big enough.
        poly = number.minpoly().change_ring(field)
        if poly.degree() > 1 and poly.is_irreducible():
            logging.info("%s is not an element in %s. Building new field." % ( number, field))
            # K = field.extension(number.minpoly(), "x%d" % (len(field.variable_names())+1))
            K = poly.splitting_field("x%d" % (len(field.variable_names())+1))

            logging.info("Polynomial now factors as %s" % poly.change_ring(K).factor())
            logging.info("Converting %s to an algebraic element over new field %s" % (number, K))
            #TODO - this doesn't work. How to convert constant symbolic expression into an element in its own splitting field?
            return K(number)
        else:
            logging.info("%s is already an element in %s. Converting it to an algebraic element." % (number, field))
            return algebraic(number, field)
    else:
        logging.info("%s is already a(n) %s" %(number, number.__class__))
        return number

def base_field(M):
    """
    Builds an extension field over QQ which contains every entry in the provided matrix.
    """

    polys = []
    # First build an array of non-linear minimal polynomials
    for row in M:
        for elem in row:
            p = elem.minpoly()
            # if the minimal polynomial in linear, then this element is rational.
            if p.degree() > 1:
                if polys.count(p.monic()) == 0:
                    polys.append(p.monic())
    logging.info("Non-linear minimal polynomials: %s" % polys)

    # Then build the field
    alph = list(build_alphabet(len(polys), 'x'))
    alph.reverse()
    K = QQ
    for p in polys:
        if p.change_ring(K).is_irreducible():
            K = K.extension(p, alph.pop())
        else:
            logging.info("The polynomial %s factors as %s over %s" % ( p, p.change_ring(K).factor(), K))

    return K
