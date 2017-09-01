from sage.symbolic.expression_conversions import algebraic

import logging
from sympy import SparseMatrix as sympySM
logging.basicConfig(level=logging.DEBUG)

def parse_kwargs(kwargs, option, default):
    """
    Parses keyword arguments. Used by other methods.

    INPUT::

        - `kwargs` --   A dictionary of the form {'option1': 'value1', 'option2': 'value2'}
        - `option` --   The relevant option's keyword.
        - `default` --  The default value. Used if `option` is not a key in `kwargs`

    OUPUT::

        Returns `kwargs`[`option`] if it exists, otherwise returns default value.
    """

    return kwargs[option] if kwargs.has_key(option) else default

def partial_trace(mat, n=2):
    """
    Takes the partial trace of `mat` over CC^n
    e.g. replaces each `n` by `n` block with its trace
    """

    if n.divides(mat.ncols()) and n.divides(mat.nrows()):
        logging.info("Taking partial trace over CC^%d" % n)
        p_trace = []
        mat.subdivide(range(n, mat.nrows()-n+1,n), range(n, mat.ncols()-n+1,n))

        for i in range(0, mat.nrows()-n+1,n):
            row = []
            for j in range(0,mat.ncols()-n+1,n):
                row.append(mat.subdivision(i/n,j/n).trace())
                logging.debug("Tracing over the subdivision starting at (%d, %d)." % ( i, j) )
            p_trace.append(row)

        logging.info("Done computing partial trace. Converting result to matrix form.")
        mat = matrix(p_trace, sparse=True)
        logging.info("Post-trace rank: %d" % mat.rank())
        return mat
    else:
        raise Exception("Can't take partial trace over CC^%d of %d by %d matrix!" % (n, mat.nrows(), mat.ncols()))


def convert_to_sympy(mat):
    """
    Converts a sage matrix into a sympy sparse matrix.

    INPUT::

        - `mat` --  Matrix.

    OUTPUT::

        A sympy SparseMatrix with the same entries as mat.
    """

    if isinstance(mat, sympy.MatrixBase):
        logging.debug("Provided matrix already in sympy.")
        return mat
    else:
        logging.debug("Converting matrix to sympy.SparseMatrix")
        return sympySM(map( lambda r: list(r), mat.rows()))

def convert_to_sage(mat):
    """
    Converts a sympy matrix into a sparse sage matrix over a symbolic ring.

    INPUT::

        - `mat` --  Matrix.

    OUPUT::

        A Sparse sage matrix.
    """

    if sage.matrix.matrix.is_Matrix(mat):
        logging.debug("Proved matrix already in sage.")
        return mat
    else:
        logging.debug("Converting matrix to Sage representation.")
        return matrix(SR, map(lambda r: list(mat.row(r)), range(0, mat.rows)), sparse=True)


def support_basis(mat, **kwargs):
    """
    Returns a basis for the support of `mat`, i.e. a list of all
    eigenvectors corresponding to non-zero eigenvalues.

    INPUT::

    - `mat` --  Either a Matrix or a list of them. Required. Sage and Sympy matrices
                are both fine. (A mixture is also fine.)

    OUTPUT::

        A list of vectors spanning the support(s) of (the matrices in) `mat`.
    """

    if isinstance(mat, list):
        evecs = []
        for m in mat:
            evecs += eigenspaces(m)
    else:
        evecs = eigenspaces(mat)

    return flatten(map(lambda space: space[2], evecs), max_level=1)


def eigenspaces(mat, **kwargs):
    """
    Computes the eigenvales and vectors of the provided matrix.

    For each distinct eigenvalue returns a list of the form (e, V, n),
    where e is the eigenvalue, V is an array of vectors, and n is the
    algebraic multiplicity.


    INPUT::

    - `mat` --  Sage or Sympy Matrix. Required.
    - `include_zero` -- Boolean. Optional, default: False. Eigenvectors
                        corresponding to the nullspace will only be included
                        in the output if this is set to True.

    OUTPUT::

        Returns a list with entries of the form (e, n, V), where `e` is the eigenvalue,
        `n` the algebraic multiplicity, and `V` is a list of corresponding eigenvectors.

    """

    include_zero = parse_kwargs(kwargs, 'include_zero', False)
    symat = convert_to_sympy(mat)

    logging.info("Computing eigenvectors & values.")

    if include_zero:
        return symat.eigenvects()
    else:
        return filter(lambda x: x[0] != 0, symat.eigenvects())


def tensor_exponential(mat, n):
    """
    Computes mat tensored with itself n times.

    INPUT::

        - `mat` --    Matrix.
        - `n` --    Integer. Tensor exponent.

    OUTPUT::

        A new matrix, mat \otimes mat \otimes \ldots \otimes mat.

    """
    tmp = mat

    for i in range(1,n):
        tmp = tmp.tensor_product(mat)

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

def base_field(mat):
    """
    Builds an extension field over QQ which contains every entry in the provided matrix.
    """

    polys = []
    # First build an array of non-linear minimal polynomials
    for row in mat:
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
