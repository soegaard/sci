#lang scribble/manual
@(require (for-label racket/base ffi/vector ffi/unsafe flomat))

@; raco scribble +m --dest html --redirect-main http://docs.racket-lang.org manual-flomat.scrbl && open html/manual-flomat.html
@(require scribble/example scribble-math )
@(require racket/format)

@; Used to reference other manuals.
@(define reference.scrbl '(lib "scribblings/reference/reference.scrbl"))
@(define math.scrbl      '(lib "math/scribblings/math.scrbl"))

@; Convenient shortcuts for writing TeX.
@(use-mathjax)
@(define mxn @${m\times n})
@(define mxm @${m\times m})
@(define nxn @${n\times n})
@(define mx1 @${m\times 1})
@(define nx1 @${n\times 1})
@(define mxk @${m\times k})
@(define nxk @${n\times k})
@(define 1xm @${1\times m})
@(define 1xn @${1\times n})
@(define A   @${A})
@(define B   @${B})
@(define B.  @~a{@${B}.})
@(define C   @${C})
@(define L   @${L})
@(define Q   @${Q})
@(define R   @${R})
@(define S   @${S})
@(define T   @${T})
@(define U   @${U})
@(define V   @${V})
@(define X   @${X})
@(define b   @${b})
@(define i   @${i})
@(define j   @${j})
@(define k   @${k})
@(define l   @${l})
@(define m   @${m})
@(define n   @${n})
@(define x   @${x})
@(define y   @${y})

@(define ith (list @${i} "'th"))
@(define jth (list @${j} "'th"))
@(define kth (list @${k} "'th"))
@(define mth (list @${m} "'th"))
@(define nth (list @${n} "'th"))

@; Note: Without `with-html5` MathJax isn't loaded.
@title[#:tag "flomat" #:style (with-html5 manual-doc-style)]{Flomat: Floating Point Matrices}

@defmodule[flomat]

This manual documents the matrix library @racketmodname[flomat].

@author[@author+email["Jens Axel Søgaard" "jensaxel@soegaard.net"]]

@local-table-of-contents[#:style 'immediate-only]

@section{Introduction}

A matrix is a rectangular arrangements of numbers in rows and columns.
This library provides functions to construct and compute with
matrices whose elements are IEEE double precision floating point numbers.
These numbers are referred to as @tech[#:doc reference.scrbl]{flonums}
in the Racket manual, but the most common name for these numbers are 
simply @emph{doubles}.

Restricting the scope of the library to dense matrices with floating numbers
allow the implementation to use routines implemented in Fortran and C.
The low-level routines consists of calls to functions in BLAS and LAPACK.
BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra PACKage)
are industry standard libraries and are available on all major platforms.

If you are in need of matrix algebra over more general numbers then
look at the functional matrix library in @secref["matrices" #:doc math.scrbl].

This library can be used in a functional manner, but imperative operations
are available. There are at least two reasons to use the imperative approach:
1) text books often describe matrix algorithms using imperative operations,
and, 2) imperative operations can reduce the amount of memory needed
during a computation.

The available operations can be divided into rough categories:
@itemlist[
  @item{Level 1:  High   level - do what I mean}
  @item{Level 2:  Medium level - do what I mean this way}
  @item{Level 3:  Low    level - do it using this underlying C-function}]

To use the library one can get by with level 1 operations, but if you understand
the underlying representation, you can improve your algorithms using 
level 2 operations. For those that really want to squeeze out the last bit of
performance we have made level 3 operations available as well. The Quick Tutorial
only describes level 1 and level 2 operations.


@section{Quick Tutorial}

@(define quick-eval (let ([e (make-base-eval)]) (e '(require flomat)) e))

This section shows how to do simple matrix computations.
The beginning part of the tutorial describes working with matrices simply as arrays
of numbers.
The end shows how to do linear algebra.

@subsection{Basic Properties}

An @racket[flomat] consists conceptually of a two-dimensional
array of floating point numbers. An @mxn (@m by @n) matrix is
divided into @m rows and @n columns. The rows are numbered @${0, \ldots, m-1}
and the columns are numbered @${0, \ldots, n-1}.

@(define (barrow elm) @~a{{\Rule 20pt 0.5pt 0pt} \mspace 2pt @elm \mspace 2pt {\Rule 20pt 0.5pt 0pt}})
@(define (barcol elm) @~a{\begin{matrix} {\Rule 0.5pt 20pt 0pt} \\ @elm \\ {\Rule 0.5pt 20pt 0pt} \end{matrix}})
@$${A = 
    \begin{bmatrix}
     a_{0,0}   & a_{0,1} & \cdots  &  a_{0,n-1}   \\ 
     a_{1,0}   & a_{1,1} & \cdots  & a_{1,n-1}    \\
     \vdots    & \vdots  & \cdots  & \vdots       \\
     a_{m-1,0} & a_{m-1,1} & \cdots & a_{m-1,n-1} \\
    \end{bmatrix}
    =
    \begin{bmatrix}  @barrow{a_{0,*}} \\ @barrow{a_{1,*}} \\ \cdots \\ @barrow{a_{{m-1},*}} \end{bmatrix}
    =
    \begin{bmatrix}  @barcol{A_0} &  @barcol{A_1}   \cdots   @barcol{A_{n-1}} \end{bmatrix}}

The basic properties of an @racket[flomat] can be examined using these functions:

@bold[@racket[(shape A)]]
  return a list of with the number of rows and columns @linebreak[]
@bold[@racket[(size A)]]
  the number of elements in the matrix                 @linebreak[]
@bold[@racket[(nrows A)]]
  the number of rows                                   @linebreak[]
@bold[@racket[(ncols A)]]
  the number of columns                                 

@examples[#:label #f #:eval quick-eval
          (define A (flomat: [[1 2 3]
                                [4 5 5]]))
          (shape A)   ; dimensions 
          (size  A)
          (nrows A)
          (ncols A)]

@subsection{Basic Indexing}

Since a matrix is divided into rows and columns we can refer to an
element in the matrix by row and column numbers. The element on the
@ith row and @jth column is referred to as the element with index @${(i,j)}.

The indices are zero-based so a matrix with @m rows and @n columns
has row-indices @${0, 1, \ldots, m-1} and column-indices @${0, 1, \ldots n-1}.

@$${A = 
    \begin{bmatrix}
     a_{0,0}   & a_{0,1} & \cdots  &  a_{0,n-1}   \\ 
     a_{1,0}   & a_{1,1} & \cdots  & a_{1,n-1}    \\
     \vdots    & \vdots  & \cdots  & \vdots       \\
     a_{m-1,0} & a_{m-1,1} & \cdots & a_{m-1,n-1} \\
    \end{bmatrix}}

@bold[@racket[(ref A i j)]]
   the element in @A with index @${(i,j)}  @linebreak[]
@bold[@racket[(row A i)]]
   the @ith row of @A                      @linebreak[]
@bold[@racket[(col A j)]]
   the @jth column of @A

Notice that row and column vectors are simply matrices with
a single row and a single column respectively

@examples[#:label #f
          #:eval quick-eval
          (define A (flomat: [[1 2 3]
                                [4 5 5]]))
          (ref A 0 1)
          (row A 0)
          (col A 1)]


@subsection{Matrix Creation}

There are several ways of creating matrices.

Use @racket[matrix] to create an @racket[flomat] from existing Racket data.
It can convert vector-of-vector-of and list-of-lists representation of matrices
into the @racket[flomat] representation. A vector of numbers or a list
of numbers will be converted into a column vector (a matrix with only one column).

Any non-floating point numbers will be converted to floating point.
The function @racket[matrix] also accepts @racket[f64vectors] as input.

@bold[@racket[(matrix obj)]]
  create a matrix with values from @racket[obj]

@examples[#:label #f #:eval quick-eval
          (matrix '[[1/2 1/3] [4 5]])
          (matrix #[#[1 2 3] #[4 5 6]])
          (matrix (list 1 2 3))
          (matrix (vector 1 2 3))
          (matrix (f64vector 1 2 3))]

After conversion the created @racket[flomat] will contain a pointer to
a newly allocated piece of memory containing the floating point numbers.
If you happen to work with data in the form of @racket[f64vector]s, then
you can avoid the allocation, if you use @racket[matrix!] instead.
If the same @racket[f64vector] is used to create two matrices with @racket[matrix!]
they will share the same backing array - so setting an element one matrix
will affect the other.

@bold[@racket[(matrix! obj)]]
  create a matrix with values from @racket[obj] avoid allocation of
  backing array if possible

@examples[#:label #f #:eval quick-eval
          (define v (f64vector 1 2 3))
          (define A (matrix! v))
          (define B (matrix! v))
          (list A B)
          (mset! A 0 0 42)
          (list A B)]
For comparision the same example with @racket[matrix]:
@examples[#:label #f #:eval quick-eval
          (define v (f64vector 1 2 3))
          (define A (matrix v))
          (define B (matrix v))
          (list A B)
          (mset! A 0 0 42)
          (list A B)]

In order to create a matrix of specific size with all zeros or all ones,
use the functions @racket[zeros] and @racket[ones]. Use @racket[eye]
to make matrix with ones on a diagonal.


@bold[@racket[(zeros n)]]
  create a square @nxn matrix with all zeros   @linebreak[]
@bold[@racket[(zeros m n)]]
  create a @mxn matrix with all zeros          @linebreak[] 
@bold[@racket[(ones n)]]
  create a square @nxn matrix with all ones    @linebreak[]
@bold[@racket[(ones m n)]]
  create a @mxn matrix with all ones           @linebreak[]
@bold[@racket[(eye m n k)]]
  create a @mxn matrix with ones on the @kth diagonal 

The arguments @n and @k are optional for @racket[eye]
and defaults to @m and @${0} respectively.

@examples[#:label #f #:eval quick-eval
          (zeros 2)
          (zeros 2 3)
          (ones 2)
          (ones 2 3)
          (list (eye 3) (eye 3 4) (eye 3 3 1) (eye 3 3 -1))]

To create ranges of values use @racket[arange] or @racket[colarange] which both work like
@racket[(matrix (range start stop step))], but avoids build an intermediary list.
The functions @racket[arange] and @racket[colarange] produce row and column vectors respectively.
The vector created has length @racket[(ceiling (/ (- stop start) step))].

@bold[@racket[(arange start stop step)]]
  create a row vector with values from start to stop (exclusively),
  here step is the gap between values     @linebreak[]
@bold[@racket[(arange start stop)]]
  like @racket[(arange start stop 1.0)]   @linebreak[]
@bold[@racket[(arange stop)]]
  like @racket[(arange 0.0 stop 1.0)]     @linebreak[]

@bold[@racket[(colarange start stop step)]] @linebreak[]
@bold[@racket[(colarange start stop)]]      @linebreak[]
@bold[@racket[(colarange start)]]           @linebreak[]
  like @racket[arange] but produces a column vector.


@examples[#:label #f #:eval quick-eval
          (arange 5 10 2)
          (arange 5 10)
          (arange 5)]
@examples[#:label #f #:eval quick-eval
          (colarange 5 10 2)
          (colarange 5 10)
          (colarange 5)]

Sometimes it is possible to keep the elements of matrix, but change its shape.

@bold[@racket[(reshape A m n)]]  @linebreak[]
  return a matrix with shape @mxn using the elements of @A,  @linebreak[]
@bold[@racket[(reshape! A m n)]] @linebreak[]
  return a matrix with shape @mxn using the elements of @A, share the backing area with @A

@examples[#:label #f #:eval quick-eval
          (arange 9)
          (reshape (arange 9) 3 3)
          (transpose (reshape (arange 9) 3 3))]

As an alternative to @racket[arange] consider using @racket[linspace], which
allow you to provide an exact endpoint.

@bold[@racket[(linspace start stop num)]] @linebreak[]
  return a column vector with @racket[num] numbers evenly spaced from
  @racket[start] to @racket[stop]   @linebreak[]
@bold[@racket[(linspace start stop num #f)]] @linebreak[]
   like @racket[(linspace start stop num)] but omit the last number

@examples[#:label #f #:eval quick-eval
          (linspace 2 4 6)
          (linspace 2 4 6 #f)]


@subsection{Elementwise Operations}

Elementwise operations (also called @emph{pointwise} operations) work on each element.
The operations are named with a beginning point.
Besides the elementwise versions of the standard arithmetic operations, 
the standard numerical functions also have elementwise counterparts.
Binary operators work both on matrices (of the same side)
and on a number and matrix.

Formally for a function @${f} of one or two arguments, the corresponding pointwise
function @${.f} satisfy: 

@$${.f( \begin{bmatrix} a_{ij} \end{bmatrix})
     = \begin{bmatrix} f(a_{ij}) \end{bmatrix}
    \textrm{ or  }  
    .f( \begin{bmatrix} a_{ij} \end{bmatrix},
       \begin{bmatrix} b_{ij} \end{bmatrix})
     = \begin{bmatrix} f(a_{ij},b_{ij}) \end{bmatrix}}


@bold[@racket[(.+ A B)]]                            @linebreak[]
@bold[@racket[(.- A B)]] and @bold[@racket[(.- A)]] @linebreak[]
@bold[@racket[(.* A B)]] @linebreak[]
@bold[@racket[(./ A B)]] and @bold[@racket[(./ A)]] @linebreak[]
  Elementwise version of the arithmetical operations. 
  The operations returns the result as a new matrix.

Note that @racket[.*] is elementwise multiplication. Use @racket[times]
to multiply two matrices in the linear algebra sense.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.- A)
          (./ A)
          (.+ A B)
          (.- A B)
          (.* A B)
          (./ A B)]
One of the arguments can be a number:
@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.+ A 1)
          (.- A 2)
          (.* 3 B)
          (./ A 4)]

The elementwise versions of the standard numerical functions are:

@bold[@racket[(.sin A)]]
@bold[@racket[(.cos A)]]
@bold[@racket[(.tan A)]]
@bold[@racket[(.exp A)]]
@bold[@racket[(.log A)]]
@bold[@racket[(.sqr A)]]
@bold[@racket[(.sqrt A)]]
@bold[@racket[(.expt A B)]]

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.sqr A)
          (.expt A 2)
          (.expt 2 A)]


The elementwise operations above all, allocate a new matrix.
If instead you want to modify the elements of an existing matrix,
the following functions are for you.

@bold[@racket[(.-! A)]]
@bold[@racket[(./! A)]]
@bold[@racket[(.sin! A)]]
@bold[@racket[(.cos! A)]]
@bold[@racket[(.tan! A)]]
@bold[@racket[(.exp! A)]]
@bold[@racket[(.log! A)]]
@bold[@racket[(.sqr! A)]]
@bold[@racket[(.sqrt! A)]]
@bold[@racket[(.expt! A B)]]
@bold[@racket[(.+! A B)]]      
@bold[@racket[(.-! A B)]]
@bold[@racket[(.*! A B)]]
@bold[@racket[(./! A B)]] @linebreak[]


For binary operations, the result is stored in the first argument.
@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.-! A)
          (.-! A)
          (.expt! A 2)
          (.expt! A 2)]

Also, if you want to store the result of an elementwise in another 
matrix @C, you can do as follows for the unary operations:


@bold[@racket[(.sin!  A C)]]
@bold[@racket[(.cos!  A C)]]
@bold[@racket[(.tan!  A C)]]
@bold[@racket[(.exp!  A C)]]
@bold[@racket[(.log!  A C)]]
@bold[@racket[(.sqr!  A C)]]
@bold[@racket[(.sqrt! A C)]]

And for the binary operations:

@bold[@racket[(.+! A B C)]]      
@bold[@racket[(.-! A B C)]]
@bold[@racket[(.*! A B C)]]
@bold[@racket[(./! A B C)]] @linebreak[]
@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.sqr! B A)
          A]

Finally, for @racket[.-!] and @racket[./!] which are both unary and binary operations
at once, use @racket[#f] as @B to get the unary version.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.-! B #f A)
          A]


@subsection{Indexing, Submatrices and Iterating}

From the section on Basic Indexing we know that the element on row @i in column @j,
has index @${(i,j)} and can be extracted with the function @racket[ref].

The @ith row and the @jth column can be extraced with @racket[row] and @racket[col]
respectively.

To get a submatrix use @racket[sub] and @racket[sub!].

@bold[@racket[(sub A i j m n)]] @linebreak[]
  Make a copy of the submatrix of @A with upper left corner in @${(i,j)} and with size mxn.

@bold[@racket[(sub! A i j m n)]] @linebreak[]
  Same as @racket[sub], but the elements are copied - the underlying 
  array of flonums are shared.


@examples[#:label #f #:eval quick-eval
          (define A (transpose (reshape (arange 25) 5 5)))
          A
          (sub A 0 0 3 2)
          (sub A 1 1 2 2)]

The function @racket[sub!] can be used to mutate part of a larger submatrix.

Let's say we have a matrix, in which we want to zero out all elements except
those on the edges. We can use @racket[sub!] to get a submatrix of the inner part,
then use @racket[zeros!] to clear the elements.

@examples[#:label #f #:eval quick-eval
          (define A (transpose (reshape (arange 10 35) 5 5)))
          A
          (define B (sub! A 1 1 3 3))
          B
          (zeros! B)
          A]

To iterate over a row or a column use @racket[in-flrow] and @racket[in-flcolumn].

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((11 22) (33 44))))
          (for/list ([   x  (in-flrow A 0)]) x)
          (for/list ([(i x) (in-flrow A 0)]) (list x i))
          (for/list ([   x  (in-flcolumn A 0)]) x)
          (for/list ([(i x) (in-flcolumn A 0)]) (list x i))]

@subsection{Basic Linear Algebra}

The basic linear algebra are @racket[plus], @racket[minus] and @racket[times],
which compute the sum, difference and product of a series of matrices.

@bold[@racket[(plus  A ...)]]
@bold[@racket[(minus A ...)]]
@bold[@racket[(times A ...)]] @linebreak[]
Computes the sum, difference and product of a series of matrices and/or numbers.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((2 0) (0 2))))
          (define B (matrix '((1 2) (3 4))))
          (define C (column 4 5))
          (plus A B)
          (plus A 10)
          (plus A 10 B)
          (minus A)
          (minus A B)
          (times A B)
          (times A 2 B)
          (times A C)]

As usual, there are variants that mutate the first given matrix
instead of allocating a new backing array of flonums.

@bold[@racket[(plus!  A B ...)]]
@bold[@racket[(minus! A B ...)]]  @linebreak[]
Like @racket[plus] and @racket[minus]  but stores
the result in @A, which must be a matrix.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((2 0) (0 2))))
          (define B (matrix '((0 2) (2 0))))
          (plus! A B)
          A]

@bold[@racket[(power A n)]] @linebreak[]
Computes the @nth power of a matrix @A, where @n is a natural number.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 1) (0 1))))
          (list (power A 0) (power A 1) (power A 2) (power A 3))]

@subsection{Matrix and Vector Products}

The inner product (also known as the dot product) of two column vectors
can be computed by @racket[dot].

@bold[@racket[(dot v w)]] @linebreak[]
Computes the inner product of two column vectors (i.e. matrices with only one column).

@examples[#:label #f #:eval quick-eval
          (define v (column -1 1))
          (define w (matrix '((2) (2))))
          (dot v w)]

The outer product of a column vector @A with @m rows and an row @B with @n columns
is an @mxn matrix @${O} with elements @${o_{i,j} = a_i\cdot b_j}.

@bold[@racket[(outer A B)]] @linebreak[]
Computes the outer product of the first column of @A and the first row of @${B}.

@examples[#:label #f #:eval quick-eval
          (define A (column 2 3))
          (define B (transpose (column 5 7)))
          (outer A B)]


The Kronecker product between two matrices @racket[A] and @racket[B] replaces
each element @racket[a] of @racket[A] with a copy of @racket[B] scaled with @racket[A].
The Kronecker product is a generalization of the outer product.

@bold[@racket[(kron A B)]] @linebreak[]
Computes the Kronecker product of the matrices @A and @${B}.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((1 1) (1 1))))
          (kron A B)]


@subsection{Matrix Decompositions}

@bold[@racket[(cholesky A)]] @bold[@racket[(qr A)]] @bold[@racket[(svd A)]] @linebreak[]
Computes the Cholesky, QR and SVD decompositions respectively.


The Singular Value Decomposition (SVD) returns three matrices: 
a unitary matrix @U, a column vector of singular values @S and 
a unitary matrix @${V^T} (@V transposed). The function @racket[diag] constructs
a diagonal matrix from the singular values.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define-values (U S VT) (svd A))
          (define Σ (diag S))
          (list U Σ VT S)
          (times U Σ VT)]

The QR Decomposition of @A consists of two matrices: an orthogonal matrix @Q
and an upper triangular matrix @R such that @${A=QR}.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define-values (Q R) (qr A))
          (list Q R)
          (times Q R)]

If the matrix @A is symmetric and positive-definite, then
the Cholesky decomposition can be computed. 
It comes in two forms

@$${A = L L^T  \textrm{ or } A = U^T U,}

where @L and @U are lower and upper triangular matrices.

Note: @racket[cholesky] does not check that the input matrix @A
is symmetric and positive definite.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (2 4))))
          (define L (cholesky A))
          (list L (transpose L))
          (times L (transpose L))
          (define U (cholesky A 'upper))
          (list (transpose U) U)          
          (times (transpose U) U)]


@subsection{Matrix Eigenvalues and Eigenvectors}

Eigenvalues and eigenvectors of a square matrix can be computed with @racket[eig] 
or, if only the eigenvalues are needed, with @racket[eigvals].
Note that even if all elements of a matrix are real, the eigenvalues in some
cases are complex. Therefore the eigenvalues are returned as a standard
Racket vector.


@bold[@racket[(eig A)]] @linebreak[]
Compute eigenvalues and right eigenvectors.

@bold[@racket[(eigvals A)]] @linebreak[]
Compute eigenvalues.

@examples[#:label #f #:eval quick-eval
          (eig     (diag '(1 2)))
          (eigvals (diag '(1 2)))
          (eig     (matrix '((1 -1) (1 1))))]
          


@subsection{Norms and Invariants}

The standard 2-norm @${|\cdot|_2} can be computed by @racket[norm].
For a column vector the norm is sometimes referred to as the length.

@bold[@racket[(norm A)]] @linebreak[]
Compute the square root om the sum of the square of all elements.

@examples[#:label #f #:eval quick-eval
          (norm    (matrix '((1 1))))
          (norm    (matrix '((1 -1) (-1 1))))]


@bold[@racket[(det A)]] @linebreak[]
Computes the determinant of a square matrix @${A}.

@examples[#:label #f #:eval quick-eval
          (det  (matrix '((1 2) (0 4))))
          (det  (matrix '((1 1) (2 2))))]

@bold[@racket[(trace A)]] @linebreak[]
Computes the trace, the sum along a diagonal, of a matrix.

@examples[#:label #f #:eval quick-eval
          (trace (matrix '((1 2) (0 4))))]


@bold[@racket[(rank A)]] @linebreak[]
Computes the rank of a square matrix.
The rank is the dimension of the column space,
which is equal to the dimension of the row space,
which is equal to the number  of non-zero singular values
in an SVD decomposition.

@examples[#:label #f #:eval quick-eval
          (rank  (matrix '((1 2) (0 4))))
          (rank  (matrix '((1 1) (2 2))))]


@subsection{Solving Equations and Inverting Matrices}

Solving linear equations are more or less the raison d'etre for matrices.
The main workhorse is @racket[mldivide], which can solve for @X
in the equation:

    @$${AX = B,}

where @A is a an @mxm matrix, and both @X and @B are @${mxn}.

Note that @A needs to be of full rank for the equation
to have a solution. The solver doesn't check that the input
matrix has full rank, it just runs it computation as usual.
To check that the output from @racket[solve] is indeed a solution,
you can evaluate @racket[(times A X)] and compare with @${B}.
The name @racket[mldivide] is short for "Matrix Left divide".
Although @racket[mldivide] doesn't find @X by 
multiplying @B with @${A^{-1}} on the left, 
it is a fitting analogy.

@bold[@racket[(mldivide A B)]] @linebreak[]
Solve the equation @${AX = B} using @${LU}-decomposition with
partial pivoting. The matrix @A must be square and of full rank, the number 
of rows in @A must be the same as the number columns in @${B}.

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((1) (0))))
          (define X (mldivide A B))
          (list X (times A X))]

@bold[@racket[(mrdivide B A)]] @linebreak[]
Solve the equation @${XA = B}.
The name @racket[mrdivide] is short for "Matrix Right divide".

@examples[#:label #f #:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((2 4) (6 8))))
          (define X (mrdivide B A))
          (list X (times X A))]


@bold[@racket[(inv A)]] @linebreak[]
Find the multiplicative inverse of a square matrix @${A}.

@examples[#:label #f #:eval quick-eval
          (define A    (matrix '((1 2) (3 4))))
          (define Ainv (inv A))
          (list Ainv (times A Ainv))]
An inverse of @A can be used to solve @${AX=B}, but
using @racket[mldivide] directly is normally better. However, let's
try to solve the equation from the previous example.

@examples[#:label #f #:eval quick-eval
          (define B    (matrix '((1) (0))))
          (define X    (times Ainv B))
          (list X (times A X))]


@bold[@racket[(pinv A)]] @linebreak[]
Find the Moore-Penrose pseudo-inverse @${A^+}of the matrix @${A}.
The matrix @A does not need to be square.
The pseudo inverse of an @mxn matrix is of size @${n\times m}.


@examples[#:label #f #:eval quick-eval
          (define A  (matrix '((1 2) (3 4))))
          (define A+ (pinv A))
          (list A+ (times A+ A A+) (times A A+ A))]
@examples[#:label #f #:eval quick-eval
          (define B  (matrix '((1 2 3) (4 5 6))))
          (define B+ (pinv B))
          (list B+ (times B+ B B+) (times B B+ B))]


@subsection{Least Square Problems}

Let @A be an @mxn matrix and let @b be an @nx1 column vector.
The equation @${Ax=b} may (depending on @A) may not have
an unique solution - or a solution at all.

As an alternative, one can look for the vector @x that minimizes:
    @$${|Ax-b|_2,}
where @${|\cdot|_2} is the Euclidean 2-norm.

The function @racket[lstsq] return the minimum norm solution @x
of the above the problem.

If @racket[lstsq] is given an @nxk matrix @B, then the
problem will be solved for each column @b of @${B}.



@bold[@racket[(lstsq A B)]] @linebreak[]
Find minimum norm solution to the least squares problem: @racket["minimize |Ax-b|"] ,
for each column @b of a larger matrix @${B}.


As an example, let's look at estimating @${b_0} and @${b_1} in the model:
    @$${y=b_0\cdot x+b_1}
given a data set consisting of corresponding @${x}- and @${y}-values.
The calculation reveals that the relation between @x and @y is @${y=2x+1}.

The matrix @X is called the @emph{design matrix} of the problem.
See @hyperlink["https://en.wikipedia.org/wiki/Design_matrix"]{Design Matrix} at Wikipedia.
In this case the design matrix has two columns: the first has the @${x}-values, the
second contains just ones.

@examples[#:label #f #:eval quick-eval
                (define xs (column 0 1 2 3))
                (define ys (column 1 3 5 7))
                (define X  (augment xs (flomat-ones (nrows xs) 1)))
                X
                (define B  (lstsq X ys))
                B]

@subsection{Matrix Functions}
@bold[@racket[(expm A)]] @linebreak[]
Compute the matrix exponential @${\exp(A)}.

@examples[#:label #f #:eval quick-eval
          (list (exp 1) (exp 2))
          (expm (matrix '((1 0) (0 2))))
          (expm (matrix '((1 2) (3 4))))]



@section{Installation}

The package @racket[flomat] is installed either in the terminal:

    @racket[raco pkg install flomat]

or using the Package Manager in DrRacket.

The package relies on the shared libraries CBLAS and LAPACK.
Depending on your OS, you might need to install these yourself.

On macOS both CBLAS and LAPACK is part of the Accelerate Framework
which is distributed by Apple. This means no extra installation is
needed.

On Linux you need copies of CBLAS and LAPACK. Since BLAS and LAPACK
exists in multiple versions, so a little care is needed. First
on most systems @racket[libblas] is used for the Fortran version,
and @racket[libcblas], so get the latter. However on Debian it turns
out @racket[libblas] is exporting the names used by CBLAS, so
(either?) ought to be fine.

On Windows: A tester is needed. Install CBLAS and LAPACK and let
me know if it works. Otherwise make an Issue at Github and we
will add the proper paths.


@section{Reference}


@subsection{Representation}

@(define (wikipedia name . preflow)
   (define url (string-append "https://en.wikipedia.org/wiki/" name))
   @margin-note{@hyperlink[url (list* @bold{Wikipedia: } " " preflow)]})

@wikipedia["Row-_and_column-major_order"]{Column Major Order}
An @mxn matrix is consist conceptually of @mxn floating points
arranged in @m rows and @n columns. Concretely the floating point
numbers are stored in an one-dimentional array in @emph{column major order}.
This means that the entries in each column are stored together in the array.

Given the address of an entry the
@deftech["leading dimension" #:key "leading dimension"] @racket[ld]
is the amount to add to get the address of the next entry in the same row.

For matrices with no gaps between columns in the array, the leading dimension
and the number of rows is the same @racket[ld=m].

For matrices with gaps between columns in the array, the leading dimension
might be larger than the number of rows @racket[ld>m].

Allowing gaps in the array allows submatrices of a larger matrix
to share the underlying array.

As an example, let's look at an @${2\times 3} matrix with leading dimension 5.

@$${A = 
    \begin{bmatrix}
     a_{00}   & a_{01}  &  a_{02}   \\ 
     a_{10}   & a_{11}  &  a_{12}   
     \end{bmatrix}}

The underlying array is:

@$${[\underbrace{\overbrace{a_{00},a_{10}}^{\text{first  column}},?,?}_{\text{ld}=5},
     \underbrace{\overbrace{a_{01},a_{11}}^{\text{second column}},?,?}_{\text{ld}=5},
     \underbrace{\overbrace{a_{02},a_{12}}^{\text{third  column}},?,?}_{\text{ld}=5}]}

Notice that the length of the underlying array is @${m\cdot\text{ld}=2\cdot 5=15}.

The main takeaway is that:

@itemlist[#:style 'ordered
  @item{A matrix has an underlying array.}
  @item{The entries in a column is stored together.}
  @item{The underlying array can be shared between matrices.}
  @item{There can be gaps between columns in the array.}]

The exact details of leading dimensions is mostly relevant if you need to
call BLAS or LAPACK functions directly.


@defthing[_flomat ctype?]
Pointer to an array of flonums.

@(defstruct flomat ([m natural?] [n natural?] [a _flomat] [lda natural?]) #:omit-constructor)
Strucure representing an @mxn matrix with leading dimension @racket[lda]
with an underlying array of floating points stored in @racket[a].

@defform*[[(define-param (m n)       A)
           (define-param (m n a)     A)
           (define-param (m n a lda) A)]]
Equivalent to @(linebreak)
    @racket[(match-define (flomat m n _ _)   A)] @(linebreak)
    @racket[(match-define (flomat m n a _)   A)] @(linebreak) 
    @racket[(match-define (flomat m n a lda) A)] @(linebreak)
respectively.

@defform[(index lda i j)]
Expands to @racket[(+ i (* j lda))] which is the array index
of the entry on row @i, column @${j}.

@defproc[(ptr-elm [a _flomat] [lda natural?] [i natural?] [j natural?]) _flomat]
Computes the address of the entry on row @i, column @${j}.

@defproc[(ptr-row [a _flomat] [i natural?]) _flomat]
Computes the adress of the beginning of row @${i}.

@defproc[(ptr-col [a _flomat] [lda natural?] [j natural?]) _flomat]
Computes the adress of the beginning of column @${j}.
Note that the leading dimenstion @racket[lda] is needed.


@defproc[(alloc-flomat [m natural?] [n natural?]) _flomat]
Allocate a floating point array with @${mn} elements and
return a tagged pointer to the array.

@defproc[(alloc-same-size-matrix [A flomat?]) _flomat]
Like @racket[alloc-flomat] but use the  dimensions @mxn
of @racket[A] to determine the size.

@subsection{Copying}

@defproc[(copy-flomat [A flomat?]) flomat?]
Return a newly allocated @racket[flomat] struct with a newly 
allocated backing array of flonums. If racket[A] has dimensions
@mxn then the newly allocated array will have length @${m\cdot n}.
The backing array of the copy can therefore be shorter than the original
backing array.


@defproc[(unsafe-vector-copy! [s natural?] [a _flomat] [lda natural?] [b natural?]) void?]
Copies @racket[s] elements from @racket[A] into @racket[B].
The elements copied has indices: @${0}, @${\text{lda}}, @${2\text{lda}} @${\ldots} .
No error checking is done.

Note that @racket[unsafe-vector-copy!] can be used to copy a column or a row
depending on the leading dimension used.

@defproc[(unsafe-matrix-copy! [m natural?] [n natural?] [a _flomat] [lda natural?] [b _flomat] [ldb natural?]) void?]
Copies the @mxn matrix @racket[A] into @racket[B].

If you need to copy @racket[A] into an index other than @${(0,0)} use
@racket[(ptr-elm b ldb i j)] to find the addres of the submatrix of @racket[B]
which has upper left corner in @${(i,j)}.

In the same manner you can use @racket[(ptr-elm a lda i j)] to find 
start of an submatrix in @racket[A].

@subsection{Simple Constructors}

@defproc[(make-flomat [m natural?] [n natural?] [x natural? 0.0]) flomat?]
Returns a flomat of dimension @mxn with an backing array of size @${mn}.
All entires are initialized to contain @racket[x].
@examples[#:label #f #:eval quick-eval (make-flomat 2 3 4)]

@defproc[(list->flomat [xss list-of-list-of-number]) flomat]
Given a matrix represented as list of rows (where a row is a list of numbers),
return a new matrix with the same entries.
Note: @racket[matrix] is usually simpler to use
@examples[#:label #f #:eval quick-eval (list->flomat '((1 2) (3 4)))]

@defproc[(vectors->flomat [xss vector-of-vector-of-number]) flomat]
Given a matrix represented as vector of rows (where a row is a vector of numbers),
return a new matrix with the same entries.
Note: @racket[matrix] is usually simpler to use
@examples[#:label #f #:eval quick-eval (vectors->flomat '#(#(1 2) #(3 4)))]

@defproc[(flomat->vectors [A flomat]) vector?]
Given a flomat @racket[A] return a matrix with the same entries represented
as vector of rows (where a row is a vector of numbers).
@examples[#:label #f #:eval quick-eval (flomat->vectors (matrix '[[1 2] [3 4]]))]

@defproc[(vector->flomat [m natural?] [n natural?] [v vector?]) flomat]
Given a vector @racket[v] of length @${mn} representing a matrix with
entries in row major order, return a matrix with the same dimensions
and entries represented as a @racket[flomat].
@examples[#:label #f #:eval quick-eval (vector->flomat 2 3 (vector 1 2 3 4 5 6))]

@defproc[(flomat->vector [A flomat]) vector?]
Return a vector of all entries in @racket[A] in row-major order.
@examples[#:label #f #:eval quick-eval (flomat->vector (matrix '[[1 2] [3 4]]))]

@defproc[(flomat/dim [m natural?] [n natural?] [xs list-of-numbers]) flomat?]
Construct a @racket[flomat?] matrix with entries from @racket[xs].
The numbers in @racket[xs] are expected to be in row major order.
@examples[#:label #f #:eval quick-eval (flomat/dim 2 3   1 2 3 4 5 6)]






@subsection{Printing}

@defproc[(flomat-print [A flomat?] [port port?] [mode boolean?]) void?]
Prints a flomat @racket[A] to the port @racket[port].
If the mode @racket[#t] means @emph{write} and
@racket[#f] means @emph{display}.

If the size of the matrix @racket[A] is less than the
value of the parameter @racket[current-max-flomat-print-size] the
entries are printed, otherwise an ellisis "..." is printed.

Currently there the output of @emph{write} and @emph{display} mode is the same.

@examples[#:label #f
          #:eval quick-eval
          (define A (matrix '[[1 2]]))
          (flomat-print A (current-output-port) #f)
          (display A)
          (flomat-print A (current-output-port) #t)
          (write A)]

@defparam[current-max-flomat-print-size n natural?]
Parameter that controls printing whether the entries of the matrix
are printed. See @racket[flomat-print].








