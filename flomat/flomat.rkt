#lang racket
(provide (all-defined-out))

;;; TODO
;;;        * Implement block. Handle mix between numbers and matrices.
;;;        * Improve matrix-expt! (avoid allocation)
;;;        * schur decomposition
;;;        * sqrtm See Higham paper.
;;;        * logm

;;; NOTES
;;;        * Contracts will be added before release
;;;        * See tests at bottom for examples.

;;; FEEDBACK
;;;        * Where is CBLAS and LAPACK on your platform
;;;          (Windows and Linux)
;;;        * What are the libraries named?
;;;        * Do all tests evaluate to #t on your platform?

;;;        * Mail: jensaxel@soegaard.net

;;;
;;; PLATFORMS TESTED       
;;;        * OS X Catalina (Working)

;;;
;;; IDEAS
;;;      Potential Improvements
;;;        * DONE Unsafe operations
;;;        * DONE add lda to the flomat structure
;;;        * DONE support shared submatrix without allocation
;;;        * DONE Improve equal?
;;;        * Use dgeequ before dgetrf (in matrix-lu!)
;;;        * Use an extra call with lwork=-1 in matrix-inverse!
;;;        * support different storage schemes 
;;;          http://www.netlib.org/lapack/lug/node121.html

;;; Useful routines to consider:

;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlazro.html
;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlaset.html
;;;   Constructs diagonal matrices. Use for flomat-identity
;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlaswp.html
;;;   Row interchanges
;;; * http://www.math.utah.edu/software/lapack/lapack-d/drscl.html
;;;   Scale by 1/a with correct rounding


(require ffi/vector
         ffi/unsafe
         ffi/unsafe/define
         racket/flonum
         (for-syntax 
          racket/format
          racket/string
          ffi/unsafe
          racket/syntax))


;;;
;;; LIBRARIES
;;;

; CBLAS and LAPACK are used.
; The first two are C-based whereas LAPACK is Fortran based.

; Note: Use trailing _ in names exported by LAPACK (in order to work both on macOS and Linux).


;; Find placement of libraries.

(define-values (cblas-lib lapack-lib)
  (case (system-type)
    ; MACOS
    [(macosx)
     (define veclib-lib 
       ; OS X: Contains CBLAS both CATLAS. CATLAS is not used here.
       ; https://developer.apple.com/library/mac/#documentation/Accelerate/
       ;         Reference/BLAS_Ref/Reference/reference.html
       (ffi-lib "/System/Library/Frameworks/vecLib.framework/Versions/Current/vecLib"))     
     (define cblas-lib veclib-lib)
     (define lapack-lib
       (ffi-lib 
        (string-append
         "/System/Library/Frameworks/Accelerate.framework/"
         "Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK")))
     (values cblas-lib lapack-lib)]
    ; UNIX
    [(unix)
     ; Note: The library names are different on Debian, Ubuntu and Arch.
     (define uname (string-downcase (system-type 'machine)))
     (define dist  (cond [(regexp-match "arch"           uname) 'arch]
                         [(regexp-match "debian"         uname) 'debian]
                         [(regexp-match "ubuntu"         uname) 'ubuntu]
                         [(regexp-match #px"fc\\d\\d"    uname) 'fedora]
                         [(regexp-match #px"raspberrypi" uname) 'rp400]  ; raspberry pi
                         [(regexp-match #px"rp400"       uname) 'rp400]  ; raspberry pi
                         [else                                  'other]))
     ; The lib order is important here. 
     ; Since cblas depends on gfortran, gfortran needs to come first.
     (define gfortran-lib (case dist
                            [(fedora) (ffi-lib "libgfortran" '("5" #f))]
                            [(rp400)  (ffi-lib "libgfortran" '("5" #f))]
                            [(ubuntu) (ffi-lib "libgfortran" '("5" "3" #f))]
                            [(debian) (ffi-lib "libgfortran" '("5" "3" #f))]
                            [else     (ffi-lib "libgfortran" '("3" #f))]))

     (define quadmath-lib (case dist
                            [(rp400) #f]
                            [else    (ffi-lib "libquadmath" '("0" #f))]))

     (define cblas-lib    (case dist
                            [(debian) (ffi-lib "libblas"  '("3" #f))]
                            [(arch)   (ffi-lib "libcblas" '("3" #f))]
                            [(ubuntu) (ffi-lib "libblas"  '("3" #f))] 
                            [(fedora) (ffi-lib "libcblas" '("3" #f))]
                            [(rp400)  (ffi-lib "libblas"  '("3" #f))]
                            [(other)  (ffi-lib "libblas"  '("3" #f))]))
                            


     (define lapack-lib   (ffi-lib "liblapack"   '("3" #f)))
     
     (values cblas-lib lapack-lib)]
    [(windows) ; Windows 10
     (define (use-openblas)
       (define cblas-lib (ffi-lib "libopenblas.dll"))
       (define lapack-lib #f)
       (values cblas-lib lapack-lib))

     ; If RACKET_SCI_USE_OPENBLAS is set, we don't look for the standard names.
     (case (getenv "RACKET_SCI_USE_OPENBLAS")       
       [(#f) (with-handlers ([exn:fail:filesystem?
                              ; the standard names weren't found, try openblas
                              (λ (x) (use-openblas))])
               (define cblas-lib  (ffi-lib "libblas.dll"))   ; place them in PATH
               (define lapack-lib (ffi-lib "liblapack.dll"))
               (values cblas-lib lapack-lib))]
       [else (use-openblas)])]))

;;; Load libraries

(define-ffi-definer define-cblas  cblas-lib)
(define-ffi-definer define-lapack lapack-lib)

;;;
;;; REFERENCES
;;; 

; LAPACK Naming scheme:
; http://www.netlib.org/lapack/lug/node24.html

; On macOS the header files are here:

; /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/System/Library/
; Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
; Versions/A/Headers/clapack.h


;;;
;;; CONFIGURATION
;;;

(define epsilon 1e-13) 
; If two flomats have the same size and
; the differences between two entries are
; smaller than epsilon, they are considered
; equal? . Furthermore if all entries are
; smaller than epsilon flomat-zero? 
; returns true.

(define current-max-flomat-print-size (make-parameter 100))
; For matrices with smaller size, all
; entries are printed. For larger matrices
; only the dimension is printed.


;;;
;;; REPRESENTATION
;;;

; BLAS/LAPACK represents matrices as one-dimensional arrays
; of numbers (S=single, D=double, X=complex or Z=double complex).
; This library uses arrays of doubles.

(define _flomat (_cpointer 'flomat))

; The array is wrapped in a struct, which besides
; a pointer to the array, holds the number of
; rows and columns. Future extension could be to
; allow different types of numbers, or perhaps
; choose specialized operations for triangular matrices.

(define (flomat-print A port mode)
  (define print (if mode write display))
  (print 
   (if (< (flomat-size A) (current-max-flomat-print-size))
       ; constructor style printing:
       (list 'flomat: ; (flomat-m A) (flomat-n A) 
             (flomat->lists A))
       ; omit actual elements
       (list 'flomat (flomat-m A) (flomat-n A) 
             "..."))
   port))



(define (flomat= A B [eps #f])
  ; TODO: Use (< (norm1 (.- A B)) eps)
  (define-param (m n a lda) A)
  (define-param (r c b ldb) B)
  (and (= m r) (= n c)
       (for*/and ([j (in-range n)]
                  [i (in-range m)])
         (define aij (unsafe-ref a lda i j))
         (define bij (unsafe-ref b ldb i j))
         (if eps
             (fl<= (flabs (fl- aij bij)) eps)
             (fl= aij bij)))))

; m = rows, n = cols, a = mxn array of doubles
; lda = leading dimension of a (see below)
(struct flomat (m n a lda)
  #:methods gen:custom-write 
  [(define write-proc flomat-print)]
  #:methods gen:equal+hash
  [(define equal-proc 
     (λ (A B rec)
       (and (= (flomat-m A) (flomat-m B))
            (= (flomat-n A) (flomat-n B))
            (or (equal? (flomat-a A) (flomat-a B))
                (flomat= A B epsilon)))))
   (define hash-proc  
     ; TODO: Avoid allocation in hash-proc.
     (λ (A rec) 
       (define-param (m n) A)
       (rec (cons m (cons n (flomat->vector A))))))
   (define hash2-proc 
     (λ (A rec) 
       (define-param (m n) A)
       (rec (cons n (cons m (flomat->vector A))))))]) 

; convenient destructuring
(define-syntax (define-param stx)
  (syntax-case stx ()
    [(_ (m n) A)
     #'(begin
         (define A1 A) 
         (define m (flomat-m A1))
         (define n (flomat-n A1)))]
    [(_ (m n a) A)
    #'(begin
        (define A1 A) 
        (define m (flomat-m A1))
        (define n (flomat-n A1))
        (define a (flomat-a A1)))]
    [(_ (m n a lda) A)
    #'(begin
        (define A1 A) 
        (define m   (flomat-m A1))
        (define n   (flomat-n A1))
        (define a   (flomat-a A1))
        (define lda (flomat-lda A1)))]
    [_
     (syntax/loc stx (error "Wrong number of arguments"))]))

;;;
;;; MEMORY LAYOUT
;;;

; The entries are layed out in column major order.
; This means that the entries in a column are
; contigious. LAPACK needs this order.

;  a[0]   a[0   +lda]  a[0   + 2*lda] ... a[0+(n-1)*lda]
;  a[1]   a[1   +lda]
;  a[2]
;  ...    ...
; a[m-1]  a[m-1 +lda]  a[m01 + 2*lda] ... a[m-1+(n-1)*lda] 

; For most matrices lda=m.

; For a submatrix it is possible that lda is larger than m.
; See http://stackoverflow.com/q/5009513/23567 
; Example:
;   If   ma=10, na=12, a=<some adress>, lda=10,
;   then mb=7,  nb=2,  b=a+3+4*lda, ldb=10 (=lda)
; represent a 7x2 submatrix whose upper, lefter
; corner in A is (3,4) (indices are 0-based).

; The array index of the (i,j)th entry is:
(define-syntax-rule (index lda i j)
    (+ i (* j lda)))

(define (ptr-elm a lda i j)
  ; address of (i,j)th element
  (ptr-add a (index lda i j) _double))

(define (shared-submatrix! A i j r s)
  ; return rxs matrix with upper left corner (i,j)
  ; entries are shared with A
  ; TODO: consider garbage collection
  (define-param (m n a lda) A)
  (flomat r s (ptr-elm a lda i j) lda))

(define (flsubmatrix A m n i j)
  ; TODO: argument order not consistent with shared-submatrix!
  ; return a the mxn submatrix of with upper 
  ; left corner in (i,j)
  (copy-flomat (shared-submatrix! A i j m n)))


(define (ptr-row a i)
  ; pointer to beginning of row a
  (ptr-add a i _double))

(define (ptr-col a lda j)
  ; address of column j
  (ptr-add a (* j lda) _double))


;;;
;;; CHECKS
;;;

(define (check-flomat who A)
  (unless (flomat? A)
    (raise-type-error who "expected flomat" A)))

(define (check-same-dimensions A B who)
  (unless (flomat-same-dimensions? A B)
    (raise-argument-error who "expected two matrices of the same size" A B)))

(define (check-all-matrices-same-size who AS)
  (set! AS (filter flomat? AS))
  (when (not( empty? AS))
    (unless (and (apply = (map flomat-m AS))
                 (apply = (map flomat-n AS)))
      (raise-argument-error 
       who 
       "All input matrices are expected to have the same dimensions."
       AS))))

(define (check-product-dimensions who A B [C #f] [transA #f] [transB #f])
  (define-values (ma na) (flomat-dimensions A))
  (define-values (mb nb) (flomat-dimensions B))
  (when transA (set!-values (ma na) (values na ma)))
  (when transB (set!-values (mb nb) (values nb mb)))
  (unless (if (not C)
              (= na mb)
              (and (= na mb)
                   (= ma (flomat-m C))
                   (= nb (flomat-n C))))
    (raise-argument-error 
     who 
     (if C
         "expected three matrices with compatible dimensions"
         "expected two matrices with compatible dimensions")
     (list (map (λ (A) (flomat-dimensions A #t)) (list A B C)) (list A B C)))))

(define (check-matrix-vector-product-dimensions who A X Y transpose-A)
  ; ma x na * mx x nx = ma x nx
  (define-param (ma na) A)
  (define-param (mx nx) X)
  (define-param (my ny) Y)
  (when transpose-A (set!-values (ma na) (values na ma)))
  (unless (= ny nx 1)
    (raise-argument-error 
     who "expected two column vectors, got: "
     (list (list mx nx) (list my ny))))  
  (unless (= na mx)
    (raise-argument-error 
     who "expected same number of columns in matrix as there are columns in X"
     (list (list ma na) (list mx nx) (list 'transpose-A transpose-A)))))

(define (check-legal-column who j A)
  (unless (< j (flomat-n A))
    (raise-argument-error 
     who "column index too large" j))
  (unless (<= 0 j)
    (raise-argument-error 
     who "column index must be non-negative")))

(define (check-legal-row who i A)
  (unless (< i (flomat-m A))
    (raise-argument-error 
     who "row index too large" i))
  (unless (<= 0 i)
    (raise-argument-error 
     who "row index must be non-negative")))

(define (check-square who A)
  (define-param (m n) A)
  (unless (= m n)
    (raise-argument-error 
     who "square matrix expected" A)))

(define (check-vector who v)
  (unless (vector? v) (raise-argument-error who "vector expected" v)))

(define (check-integer who x)
  (unless (integer? x) (raise-argument-error who "integer expected" x)))

(define (check-positive-integer who x)
  (unless (and (integer? x) (>= x 0))
    (raise-argument-error who "positive integer epected" x)))



;;;
;;; SIZE and DIMENSION
;;;

(define (flomat-size A)
  (check-flomat 'flomat-size A )
  (define-param (m n) A)
  (* m n))

(define (flomat-dimensions A [as-list? #f])
  (check-flomat 'flomat-dimensions A)
  (define-param (m n) A)
  (if as-list?
      (list m n)
      (values m n)))

(define (flomat-same-dimensions? A B)
  (define-param (ma na) A)
  (define-param (mb nb) B)
  (and (= ma mb) (= na nb)))

(define (flomat-row-vector? A)
  (= 1 (flomat-m A)))

(define (flomat-column-vector? A)
  (= 1 (flomat-n A)))

;;;
;;; ALLOCATIONS and CONSTRUCTORS
;;; 

(define (alloc-flomat m n)
  (if (or (= m 0) (= n 0))
      #f ; ~ NULL
      (cast (malloc (* m n) _double 'atomic)
            _pointer _flomat)))

; Note: Even though we use ptr-add we do not need to use 'atomic-interior
;       since the result of ptr-add contains both the base pointer and an offset.

(define (alloc-same-size-matrix A)
  (define-param (m n) A)
  (alloc-flomat m n))


(define-syntax (define-cblas* stx)
  (syntax-case stx ()
    [(def xname _x (c ...) body ...)
     (let ()       
       (define ((xname->name ctx xname) c)
         (datum->syntax 
          ctx
          (string->symbol
           (string-replace (~a xname) "x" (~a c) #:all? #f))))
       (define (c->_c c)
         (unless (symbol? c)
           (error (format "expected symbol, got: ~a" c)))
         (case c
           [(c) #'_double] ; TODO missing from ffi?
           [(z) #'_double] ; TODO
           [(d) #'_double]
           [(s) #'_float]
           [else (error "expected one of c, z, d, s")]))
       (with-syntax ([(name ...) 
                      (map (xname->name stx (syntax->datum #'xname))
                           (syntax->datum #'(c ...)))]
                     [(_c ...) (map c->_c (syntax->datum #'(c ...)))])
         #'(begin
             (define-cblas name 
               (let ([_x _c]) body ...))
             ...)))]))

;; Note: BLAS expects complex numbers to be passed as two doubles,
;        so _complex aren't per se missing from the FFI.
;   https://stackoverflow.com/questions/34103818/swift-blas-cblas-cgemv-complex-numbers
;  For now this library sticks with doubles.

(define-cblas* cblas_xcopy _x (s d c z)  ; dcopy
  ; copy n elements from vector X to vector Y
  (_fun (n : _int)  
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _void))

#;(define-cblas cblas_dcopy
  ; copy n elements from vector X to vector Y
  (_fun (n : _int)  
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _void))

(define (unsafe-vector-copy! s a lda b)
  ; copy s elements from A into B
  ; element 0, lda, 2*lda, ... is copied
  (cblas_dcopy s a lda b 1))

(define (unsafe-matrix-copy! m n a lda b ldb)
  ; Todo: The for loop currently copies each column.
  ;       If the rows are longer than columns then it would be
  ;       faster to copy each row.
  
  ; copy the mxn matrix A into B
  ; copy has upper left corner in (i,j)
  ; Note: use (ptr-elm b ldb i j) to
  ;       copy into a submatrix of b.
  (for ([j (in-range n)])
    (unsafe-vector-copy! 
     m (ptr-elm a lda 0 j) 1 
     (ptr-add b (* j ldb) _double))))

(define (copy-flomat A)
  (define-param (m n a lda) A)
  (define size (* m n))
  (define b (cast (malloc size _double 'atomic)
                  _pointer _flomat))
  (define ldb m)
  (cond 
    [(= lda m) ; elements in a are contigious
     (unsafe-vector-copy! size a 1 b)]
    [else ; copy each column separately
     (unsafe-matrix-copy! m n a lda b ldb)])
  (flomat m n b ldb))


;; (define (make-flomat m n [x 0.0])
;;   (define a (alloc-flomat m n))
;;   (define x* (real->double-flonum x))
;;   (if (= x 0.0)
;;       (memset a 0 (* m n) _double)
;;       (for ([i (* m n)]) (ptr-set! a _double i x*)))
;;   (flomat m n a m))

(define (make-flomat m n [x 0.0])
  (define a  (alloc-flomat m n))
  (define x* (cast (malloc 1 _double 'atomic) _pointer _flomat))
  (ptr-set! x* _double (real->double-flonum x))
  (if (= x 0.0)
      (memset a 0  (* m n) _double)
      (cblas_dcopy (* m n) x* 0 a 1))
  (flomat m n a m))


(define (flomat-zeros m n)
  (make-flomat m n 0.0))

(define (flomat-ones m n)
  (make-flomat m n 1.0))


(define (list->flomat xss)
  (define m (length xss))
  (define n (apply max (map length xss)))
  (for*/flomat m n
                 ([xs (in-list xss)]
                  [x  (in-list xs)])
                 x))

(define (vectors->flomat xss)
  (define m (vector-length xss))
  (define n (vector-length (vector-ref xss 0)))
  (for*/flomat m n
                 ([xs (in-vector xss)]
                  [x  (in-vector xs)])
                 x))

(define (flomat-identity m)
  (define A (make-flomat m m 0.0))
  (for ([i (in-range m)])
    (flomat-set! A i i 1.0))
  A)

(define (flomat-column A j)
  ; copy column j
  (check-legal-column 'flomat-column j A)  
  (define-param (m n) A)
  (copy-flomat (shared-submatrix! A 0 j m 1)))

(define (flomat-row A i)
  ; copy row i
  (define-param (m n) A)
  (check-legal-row 'flomat-row i A)
  (copy-flomat (shared-submatrix! A i 0 1 n)))

;;;
;;; CONVERSIONS MATRIX <-> VECTOR
;;;

(define (flomat->vector A)
  ; the result vector uses row-major order
  (define-param (m n a lda) A)
  (for*/vector #:length (* m n)
    ([i (in-range 0 m)]
     [j (in-range 0 n)])
    (unsafe-ref a lda i j)))

(define (flomat->vectors A)
  ; the result is a vector of rows
  (define-param (m n a lda) A)
  (for/vector #:length m
    ([i (in-range 0 m)])
    (for/vector #:length n
        ([j (in-range 0 n)])
      (ptr-ref a _double (+ i (* j lda)))
      #;(ptr-ref (ptr-elm a lda i j) _double))))

(define (vector->flomat m n v)
  (unless (= (* m n) (vector-length v))
    (raise-argument-error
     'vector->flomat
     "expected m*n to be the same as the length of the vector"))
  (define a (alloc-flomat m n))
  (define k 0)
  (for* ([j (in-range n)]
         [i (in-range m)])
    (ptr-set! a _double* k ; (index m i j) 
              (vector-ref v (+ (* i n) j)))
    (set! k (+ k 1)))  
  (flomat m n a m))

; (: matrix/dim : Integer Integer Number * -> (Matrix Number))
; construct a mxn flomat with elements from the values xs
; the length of xs must be m*n
(define (flomat/dim m n . xs)
  (vector->flomat m n (list->vector xs)))


;;;
;;; COMPREHENSIONS
;;;

; (for/flomat m n (clause ...) . defs+exprs)
; (for/matrix (i in m) (j in n) (clauses ...) . body)
;    Return an  m x n  flomat with elements from the last expr.
;    The first n values produced becomes the first row.
;    The next n values becomes the second row and so on.
;    The bindings in clauses run in parallel.
(define-syntax (for/flomat stx)
  (syntax-case stx (in)
    [(_for/matrix (i in m-expr)  (j in n-expr) #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define a (alloc-flomat m n))
         (define idx 0)
         (for* ([j (in-range n)]
                [i (in-range m)]                
                clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx 1)))
         (flomat m n a m)))]
    ; elements in column 0 are generated first, then column 1, ...
    [(_ m-expr n-expr #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n  n-expr])
         (define a (alloc-flomat m n))
         (define size (* m n))
         (for ([idx (in-range size)] clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x))
         (flomat m n a m)))]
    [(_for/matrix (i in m-expr) (j in n-expr) (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define size (* m n))
         (define a (alloc-flomat m n))
         (define idx 0)
         (for* ([i (in-range m)]
                [j (in-range n)]                                
                clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flomat m n a m)))]
    ; elements in row 0 are generated first, then row 1, ...
    [(_ m-expr n-expr (clause ...) . defs+exprs)
     (syntax/loc stx
       (let* ([m m-expr] [n n-expr])
         (define a (alloc-flomat m n))
         (define idx 0)
         (define size (* m n))
         (for ([k (in-range size)] clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flomat m n a m)))]))

; (for*/flomat m n (clause ...) . defs+exprs)
;    Return an  m x n  flomat with elements from the last expr.
;    The first n values produced becomes the first row.
;    The next n values becomes the second row and so on.
;    The bindings in clauses run nested.
; (for*/flomat m n #:column (clause ...) . defs+exprs)
;    Return an  m x n  flomat with elements from the last expr.
;    The first m values produced becomes the first column.
;    The next m values becomes the second column and so on.
;    The bindings in clauses run nested.

(define-syntax (for*/flomat stx)
  (syntax-case stx ()
    [(_ m-expr n-expr #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let* ([m  m-expr] [n  n-expr])
         (define a (alloc-flomat m n))
         (define idx 0)
         (define size (* m n))
         (for* (clause ... #:break (= idx size))
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx 1)))
         (flomat m n a m)))]
    [(_ m-expr n-expr (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define a (alloc-flomat m n))
         (define idx 0)
         (define size (* m n))
         (for* (clause ... #:final (= idx (- size 1)))
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flomat m n a m)))]))

(define-syntax (for/flomat-sum stx)
  (syntax-case stx ()
    [(_ (for:-clause ...) . defs+exprs)
     (syntax/loc stx
       (let ()
         (define sum #f)
         (for (for:-clause ...)
           (define a (let () . defs+exprs))
           (set! sum (if sum (flomat+ sum a) a)))
         sum))]))

;;;
;;; BINARY MATRIX OPERATIONS
;;; 

;;; MATRIX SUM AND DIFFERENCE

(define-cblas* cblas_xaxpy _x (s d #;c #;z)
  ; Y := αX+Y  ; X and Y are vectors
  ; If incX=3 then every 3rd element of X is used.
  (_fun (n : _int) (alpha : _x)
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _void))

#;(define-cblas cblas_daxpy
  ; Y := αX+Y  ; X and Y are vectors
  ; If incX=3 then every 3rd element of X is used.
  (_fun (n : _int) (alpha : _double) 
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _void))

(define (unsafe-vector-clear n a [lda 1])
  (cblas_daxpy n -1.0 a lda a lda))

; TODO: Allow adding row to different matrix!

(define (flomat-add-scaled-row! A i1 s i2)
  ; scale row i2 and add to row i1
  (check-legal-row 'matrix-add-scaled-row! i1 A)
  (check-legal-row 'matrix-add-scaled-row! i2 A)
  (define-param (m n a lda) A)
  (define rowi1 (ptr-row a i1))
  (define rowi2 (ptr-row a i2))
  (define s* (real->double-flonum s))
  (cblas_daxpy n s* rowi2 lda rowi1 lda)
  A)

(define (flomat-add-scaled-row A i1 s i2)
  (define B (copy-flomat A))
  (flomat-add-scaled-row! B i1 s i2)
  B)

(define (flomat-add-scaled-column! A j1 s j2)    
  (check-legal-row 'flomat-add-scaled-column! j1 A)
  (check-legal-row 'flomat-add-scaled-column! j2 A)
  (define-param (m n a lda) A)
  (define colj1 (ptr-col a lda j1))
  (define colj2 (ptr-col a lda j2))
  (define s* (real->double-flonum s))
  (cblas_daxpy m s* colj1 1 colj2 1)
  A)

(define (flomat-add-scaled-column A i1 s i2)
  (define B (copy-flomat A))
  (flomat-add-scaled-column! B i1 s i2)
  B)

(define (constant*flomat+flomat! alpha A B)
  ; B := αA+B  
  (define-param (m n a lda) A)
  (define-param (r s b ldb) B)
  (for ([j (in-range n)])
    (cblas_daxpy m alpha 
                 (ptr-col a lda j) 1 
                 (ptr-col b ldb j) 1))
  B)

(define (constant*flomat+flomat alpha A B)
  ; αA+B
  (define αA+B (copy-flomat B))
  (constant*flomat+flomat! alpha A αA+B)
  αA+B)

(define (flomat+! A B)
  ; B := A + B
  (check-same-dimensions A B 'flomat+!)
  (constant*flomat+flomat! 1.0 A B))

(define (flomat+ A B)
  ; A + B
  (check-same-dimensions A B 'flomat+)
  (constant*flomat+flomat 1.0 A B))

(define (flomat-! A B)
  ; A := A - B
  (check-same-dimensions A B 'flomat-!)
  (constant*flomat+flomat! -1.0 B A))

(define (flomat- A [B #f])
  (cond
    [B
     (check-same-dimensions A B 'flomat-)
     (constant*flomat+flomat -1.0 B A)]
    [else
     (flomat-scale -1.0 A)]))

;;; Matrix x Matrix Multiplication

(define _CBLAS_ORDER _int)
(define CblasRowMajor 101)
(define CblasColMajor 102)

(define _CBLAS_TRANSPOSE _int)
(define CblasNoTrans   111)
(define CblasTrans     112)
(define CblasConjTrans 113)

(define-cblas* cblas_xgemm _x (s d z c)
  ; C := α(A*B)+βC 
  ; 1. Multiplies A and B.
  ; 2. Scales result with alpha
  ; 3. Scales C with beta.
  ; 4. Stores sum in in C.
  (_fun (order : _CBLAS_ORDER) 
        (transa : _CBLAS_TRANSPOSE) ; transpose A?
        (transb : _CBLAS_TRANSPOSE) ; transpose B?
        (m : _int) ; rows in A and C
        (n : _int) ; cols in B and C
        (k : _int) ; cols in A = rows in B
        (alpha : _x) ; scaling factor for A and B
        (A : _flomat) 
        (lda : _int) ; size of first dim of A
        (B : _flomat) 
        (ldb : _int) ; size of first dim of B
        (beta : _double) ; scaling for C
        (C : _flomat) 
        (ldc : _int) ; size of first dim of C
        -> _void))

(define (constant*matrix*matrix+constant*matrix! alpha A B beta C transA transB)
  ; C := α(A*B)+βC, maybe transpose A and/or B first
  ; Note: the check fails when the matrices are transposed.
  ; todo: pass transA and transB to the checker.
  (check-product-dimensions 'constant*matrix*matrix+constant*matrix!
                            A B C transA transB)
  (define-param (m n a lda) A)
  (define-param (r s b ldb) B)
  (define-param (x y c ldc) C)
  (define alpha* (real->double-flonum alpha))
  (define beta*  (real->double-flonum beta))
  (cblas_dgemm CblasColMajor 
               (if transA CblasTrans CblasNoTrans)
               (if transB CblasTrans CblasNoTrans)
               (if transA n m) ; rows in A
               (if transB r s) ; cols in B
               (if transA m n) ; cols in A
               alpha* 
               a lda  b ldb  beta*  c ldc)
  C)

(define (flomat*! A B C 
                    [alpha 1.0] [beta 1.0] 
                    [transpose-A #f] [transpose-B #f])
  ; C := α(A*B)+βC, maybe transpose A and/or B first 
  (constant*matrix*matrix+constant*matrix! 
   alpha A B beta C transpose-A transpose-B))

(define (flomat* A B [C #f]
                   [alpha 1.0] [beta 1.0] 
                   [transpose-A #f] [transpose-B #f])
  ; C := α(A*B)+βC, maybe transpose A and/or B first
  (define-values (ma na) (flomat-dimensions A))
  (define-values (mb nb) (flomat-dimensions B))
  (when transpose-A (set!-values (ma na) (values na ma)))
  (when transpose-B (set!-values (mb nb) (values nb mb)))  
  (define C1 (or C (make-flomat ma nb)))
  (flomat*! A B C1 alpha beta transpose-A transpose-B))

;;; Matrix Power

(define (flomat-expt a n)
  (check-flomat 'flomat-expt a)
  (check-square 'matrix-expt a)
  (cond
    [(= n 0)  (flomat-identity (flomat-m a))]
    [(= n 1)  (copy-flomat a)]
    [(= n 2)  (flomat* a a)]
    [(even? n) (let ([a^n/2 (flomat-expt a (quotient n 2))])
                 (flomat* a^n/2 a^n/2))]
    [else     (flomat* a (flomat-expt a (sub1 n)))]))

;;; Matrix x Vector Multiplication

; NOTE: Functions accepting column vectors automatically
;       convert (standard) vectors into mx1 matrices.

(define-cblas* cblas_xgemv _x (s d c z) ; Double GEneral Matrix Vector multiplication
  ; Y := α(AX) +(βY) 
  (_fun (order : _CBLAS_ORDER) 
        (transa : _CBLAS_TRANSPOSE) ; transpose A?
        (m : _int) ; rows in A 
        (n : _int) ; cols in A 
        (alpha : _x) ; scaling factor for A 
        (A : _flomat) 
        (lda : _int) 
        (X : _flomat) ; vector
        (ldx : _int) 
        (beta : _x) ; scaling for Y
        (Y : _flomat) ; vector
        (ldy : _int) 
        -> _void))

(define (constant*matrix*vector+constant*vector! alpha A X beta Y transA)
  ; unsafe:  Y := α(AX) +(βY), maybe transpose A first 
  (define-param (m n a lda) A)
  (cblas_dgemv CblasColMajor 
               (if transA CblasTrans CblasNoTrans)
               m n
               (real->double-flonum alpha)
               a lda
               (flomat-a X) 1 
               (real->double-flonum beta)
               (flomat-a Y) 1)
  Y)


(define (flomat*vector! A X Y [alpha 1.0] [beta 1.0] 
                          [transpose-A #f])
  (define X1 (result-flcolumn X))
  (define Y1 (result-flcolumn Y))  
  (check-matrix-vector-product-dimensions 
   'constant*matrix*vector+constant*vector! A X1 Y1 transpose-A)
  ; Y := α(AX) +(βY), maybe transpose A first 
  (constant*matrix*vector+constant*vector! 
   alpha A X1 beta Y1 transpose-A))

(define (flomat*vector A X [Y #f] [alpha 1.0] [beta 1.0] 
                         [transpose-A #f] )
  ; Y := α(AX) +(βY), maybe transpose A first
  (define Y1 (or Y (make-flomat (if transpose-A (flomat-n A) (flomat-m A)) 1 0.0)))
  (flomat*vector! A X Y1 alpha 1.0 transpose-A))

;;;
;;; ELEMENT WISE OPERATIONS
;;;

;;; Ref

(define (unsafe-ref a lda i j)
  (ptr-ref (ptr-elm a lda i j) _double))

(define (flomat-ref A i j)
  (define-param (m n a lda) A)
  (unless (< -1 i m)
    (raise-arguments-error 
     'matrix-ref (format "expected row index between 0 and ~a, got ~a" m i)))
  (unless (< -1 j n)
    (error 'matrix-ref 
           (format "expected column index between 0 and ~a, got ~a" n j)))
  (unsafe-ref a lda i j))

;;; Set!

(define (unsafe-set! a lda i j x)
  (ptr-set! (ptr-elm a lda i j) _double x))

(define (flomat-set! A i j x)
  (check-legal-row    'flomat-set! i A)
  (check-legal-column 'flomat-set! j A)
  (define-param (m n a lda) A)
  (define x* (real->double-flonum x))
  (unsafe-set! a lda i j x*)
  A)

;;; Scaling

(define-cblas* cblas_xscal _x (s d c z)
  ; X := αX  vector
  (_fun (n : _int) (alpha : _x) 
        (X : _flomat) (incX : _int)        
        -> _void))

(define (constant*matrix! s A)
  ; A := s*A
  (define-param (m n a lda) A)
  (define s* (real->double-flonum s))
  (cond 
    [(= lda m)
     (cblas_dscal (* m n) s* a 1)]
    [else 
     (for ([j (in-range n)])
       (cblas_dscal m s* (ptr-col a lda j) 1))])
  A)

(define (flomat-scale! s A)
  ; A := s*A
  (constant*matrix! s A))

(define (flomat-scale s A)
  ; s*A
  (define sA (copy-flomat A))
  (flomat-scale! s sA))

(define (shared-column-flomat A j)
  (check-legal-column 'shared-column-flomat j A)
  (define-param (m n) A)
  (shared-submatrix! A 0 j m 1))

(define (shared-row-flomat A i)
  (check-legal-row 'shared-row-flomat i A)
  (shared-submatrix! A i 0 1 (flomat-n A)))

(define (flomat-scale-column! A j s)
  ; col_j := s * col_j
  (constant*matrix! s (shared-column-flomat A j))
  A)

(define (flomat-scale-column A j s)
  (define B (copy-flomat A))
  (flomat-scale-column! B j s)
  B)

(define (flomat-scale-row! A i s)
  ; row_i := s * rwo_i
  (check-legal-row 'flomat-scale-row! i A)
  (define-values (m n) (flomat-dimensions A))
  (constant*matrix! s (shared-row-flomat A i))
  A)

(define (flomat-scale-row A i s)
  (define B (copy-flomat A))
  (flomat-scale-row! B i s)
  B)

;;; Swapping

(define-cblas* cblas_xswap _x (s d c z)
  ; Swaps elements in the vectors x and y
  (_fun (n : _int) ; length of vector
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _void))

(define (flomat-swap-rows! A i1 i2)
  (check-legal-row 'flomat-swap-rows! i1 A)
  (check-legal-row 'flomat-swap-rows! i2 A)
  (unless (= i1 i2)
    (define-param (m n a lda) A)
    (define rowi1 (ptr-row a i1))
    (define rowi2 (ptr-row a i2))
    (cblas_dswap n rowi1 lda rowi2 lda))
  A)

(define (flomat-swap-rows A i1 i2)
  (define B (copy-flomat A))
  (flomat-swap-rows! B i1 i2)
  B)

(define (flomat-swap-columns! A j1 j2)
  (check-legal-row 'flomat-swap-columns! j1 A)
  (check-legal-row 'flomat-swap-columns! j2 A)
  (unless (= j1 j2)
    (define-param (m n a lda) A)
    (define colj1 (ptr-col a lda j1))
    (define colj2 (ptr-col a lda j2))
    (cblas_dswap m colj1 1 colj2 1))
  A)

(define (flomat-swap-columns A j1 j2)
  (define B (copy-flomat A))
  (flomat-swap-columns! B j1 j2)
  B)

(define (flomat-flip-left-to-right! A)
  (check-flomat 'flomat-flip-left-to-right! A)
  (define-param (m n a lda) A)
  (unless (= n 1)
    (for ([j (in-range (quotient n 2))])
      (flomat-swap-columns! A j (- n j 1))))
  A)

(define (flomat-flip-left-to-right A)
  (check-flomat 'flomat-flip-left-to-right A)
  (define B (copy-flomat A))
  (flomat-flip-left-to-right! B)
  B)

(define (flomat-flip-up-to-down! A)
  (check-flomat 'flomat-flip-up-to-down! A)
  (define-param (m n a lda) A)
  (unless (= m 1)
    (for ([i (in-range (quotient m 2))])
      (flomat-swap-rows! A i (- m i 1))))
  A)

(define (flomat-flip-up-to-down A)
  (check-flomat 'flomat-flip-up-to-down A)
  (define B (copy-flomat A))
  (flomat-flip-up-to-down! B)
  B)

(define (flomat-rotate90 A)
  ; 1 2 3          3 6 9
  ; 4 5 6 becomes  2 5 8
  ; 7 8 9          1 4 7
  (check-flomat 'flomat-rotate90 A)
  (check-square   'flomat-rotate90 A)
  (define-param (m n a lda) A)  
  (for*/flomat n m
                 ([i (in-range (- n 1) -1 -1)]
                  [j (in-range (- m 1) -1 -1)])
                 (flomat-ref A (- n j 1) i)))




;;; Max Absolute Value

(define-cblas* cblas_ixamax _x (s d c z)
  ; Returns the index of the element with the largest 
  ; absolute value in a vector.
  (_fun (n : _int) (X : _flomat) (incX : _int)
        -> _int))


(define (flomat-max-abs-index A)
  (define-param (m n a lda) A)
  (cond
    [(= m lda) 
     (define idx (cblas_idamax (* m n) a lda))
     (values (remainder idx m) (quotient idx m))]
    [(= n 1)
     (define idx (cblas_idamax m a 1))
     (values (- idx 1) 0)]
    [else
     (define idx (make-vector n))
     (for ([j (in-range n)])
       (define i (cblas_idamax m (ptr-col a lda j) 1))
       (vector-set! idx j (cons (cons i j) (unsafe-ref a lda i j))))
     (define ij (car (vector-argmax cdr idx)))
     (values (car ij) (cdr ij))]))

(define (flomat-max-abs-value A)
  (define-values (i j) (flomat-max-abs-index A))
  (flomat-ref A i j))

(define (flomat-zero? A [eps epsilon])
  ; set eps=#f to use normal equal?
  (define val (flomat-max-abs-value A))
  (if eps (< (abs val) eps) (zero? val)))

(define (flomat-ones? A [eps #f])
  ; is A a square matrix with ones on the main diaginal and zero elsewhere?
  (define-param (m n a lda) A)  
  (and (= m n)
       (for*/and ([j (in-range n)]
                  [i (in-range m)])
         (define aij (unsafe-ref a lda i j))
         (if (= i j)
             (if eps
                 (fl<= (flabs (fl- aij 1.0)) eps)
                 (fl= aij 1.0))
             (if eps
                 (fl<= (flabs aij) eps)
                 (fl= aij 0.0))))))



;;;
;;; BLOCK LEVEL OPERATIONS
;;;



(define (calculate-row-height Xs)
  ; Given a list of matrix-or-integer, calculate the row height.
  ; #f means no common height
  (define Ms (for/list ([X Xs] #:when (flomat? X)) X))
  (cond [(empty? Ms) 1] ; no matrices in Xs => no height => default is 1
        [else        (define hs (map flomat-m Ms))
                     (if (apply = hs) (first hs) #f)]))

  


(define (flomat-augment C . Cs)
  ; 1. Check that all have same number of rows.
  (define-param (mc nc c ldc) C)
  (define rows (map flomat-m (cons C Cs)))
  (unless (andmap (λ (r) (= mc r)) rows)
    (raise-arguments-error 
     'flomat-augment
     "all arguments must have same number of rows"))
  ; 2. Find size for result matrix and allocate
  (define m mc)
  (define n (apply + (map flomat-n (cons C Cs))))
  (define a (alloc-flomat m n))
  (define lda m)
  ; 3. Fill in blocks
  (define j 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define aj (ptr-col a lda j))
    (unsafe-matrix-copy! mb nb b ldb aj lda)
    (set! j (+ j nb)))
  (flomat m n a lda))

(define (flomat-stack C . Cs)
  ; 1. Check that all have same number of columns
  (define-param (mc nc c ldc) C)
  (define cols (map flomat-n (cons C Cs)))
  (unless (andmap (λ (x) (= x nc)) cols)
    (raise-arguments-error 
     'flomat-stack
     "all arguments must have same number of columns"))
  ; 2. Find size for result matrix and allocate
  (define rows (map flomat-m (cons C Cs)))
  (define m (apply + rows))
  (define n nc)
  (define a (alloc-flomat m n))
  (define lda m)  
  ; 3. Fill in blocks
  (define i 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define ai (ptr-row a i))
    (unsafe-matrix-copy! mb nb b ldb ai lda)
    (set! i (+ i mb)))
  (flomat m n a lda))

(define (flomat-block-diagonal C . Cs)
  (define rows (map flomat-m (cons C Cs)))
  (define cols (map flomat-n (cons C Cs)))
  ; 2. Find size for result matrix and allocate
  (define m (apply + rows))
  (define n (apply + cols))
  (define a (alloc-flomat m n))
  (define lda m)
  (unsafe-vector-clear (* m n) a)
  ; 3. Fill in blocks
  (define i 0)
  (define j 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define aij (ptr-elm a lda i j))
    (unsafe-matrix-copy! mb nb b ldb aij lda)
    (set! i (+ i mb))
    (set! j (+ j nb)))
  (flomat m n a lda))

(define (flomat-repeat A m [n m])
  ; Make a matrix with mxn blocks, each block is A.
  (define row (apply flomat-augment (for/list ([i m]) A)))
  (apply flomat-stack (for/list ([j n]) row)))
  
  

;;;
;;; NORMS
;;;

(define-cblas* cblas_xnrm2 _x (s d)
  ; L2-norm = (sqrt (sum (sqr X_i))), vector
  (_fun (n : _int) (X : _flomat) (incX : _int)
        -> _x))

#;(define (flomat-norm A)
  ; (sqrt (sum (sqr A_ij)))
  (define-param (m n a lda) A)
  (cond
    [(= lda m)
     (cblas_dnrm2 (* m n) a 1)]
    [(= n 1)
     (cblas_dnrm2 m a 1)]
    [else
     (sqrt
      (for/sum ([j (in-range n)])
        (expt (cblas_dnrm2 m (ptr-col a lda j) 1) 2)))]))

(define-lapack dlange_
  (_fun (norm : (_ptr i _byte)) ; char M, 1, I, F
        (m    : (_ptr i _int))
        (n    : (_ptr i _int))
        (a    : _flomat) 
        (lda  : (_ptr i _int))
        (work : _flomat)        ; used only if norm is #\M
        -> _double))



(define (flomat-norm A [norm-type 'frob])
  (define-param (m n a lda) A)
  (define norm (char->integer 
                (match norm-type
                  [1      #\1]
                  ['inf   #\I]
                  ['frob  #\F]
                  ['max   #\M]
                  [_ (error)])))
  (define lwork (if (equal? norm-type 'inf) (max 1 m) 1))
  (define W (make-flomat lwork 1))
  (define w (flomat-a W))
  (dlange_ norm m n a lda w))


(define (flomat-norm1 A) ; maximum column sum (absolute values)
  (flomat-norm A 1))

(define (flomat-norm-inf A) ; maximum row sum (absolute values)
  (flomat-norm A 'inf))

(define (flomat-norm-frob A) ; Frobenius (sqrt of sum of squares)
  (flomat-norm A 'frob))

(define (flomat-norm-max A) ; not real norm
  (flomat-norm A 'max))


;;; 
;;; UNARY MATRIX OPERATIONS
;;;

(define (flomat-transpose A)
  ; TODO: Measure: Is it faster to use 
  ;       a loop with unsafe-vector-copy ?
  (define-param (m n a lda) A)
  (define AT (make-flomat n m))
  (define at (flomat-a AT))
  (for* ([j (in-range n)]
         [i (in-range m)])
    (unsafe-set! at n j i (unsafe-ref a lda i j)))
  AT)

;;;
;;; Eigenvalues and Eigenvectors
;;;

(define-lapack dgeev_ 
  ; http://www.netlib.org/lapack/lapack-3.1.1/html/dgeev.f.html
  ; DGEEV computes for an N-by-N real nonsymmetric matrix A, the
  ; eigenvalues and, optionally, the left and/or right eigenvectors.
  ;
  ; The right eigenvector v(j) of A satisfies
  ;                   A * v(j) = lambda(j) * v(j)
  ; where lambda(j) is its eigenvalue.
  ; The left eigenvector u(j) of A satisfies
  ;                u(j)**H * A = lambda(j) * u(j)**H
  ; where u(j)**H denotes the conjugate transpose of u(j).
  ;
  ; The computed eigenvectors are normalized to have Euclidean norm
  ; equal to 1 and largest component real.

  (_fun (jobvl : (_ptr i _byte)) ; char 'N' or 'V'
        (jobvr : (_ptr i _byte)) ; char 'N' or 'V'
        (n     : (_ptr i _int))  ; order of a
        (a     : _flomat)      ; io: the matrix
        (lda   : (_ptr i _int))  
        (wr    : _flomat)      ; out: real part of eigenvalues
        (wi    : _flomat)      ; out: imag part of eigenvalues
        (vl    : _flomat)      ; left eigenvectors
        (ldvl  : (_ptr i _int))
        (vr    : _flomat)      ; right eigenvectors
        (ldvr  : (_ptr i _int))
        (work  : _flomat)      ; dim max(1,lwork)
        (lwork : (_ptr i _int))  ; dim >= 4n 
        (info  : (_ptr o _int))
        -> _void
        -> info))

(define (flomat-eigenvalues-and-vectors!
         A #:left [left? #f] #:right [right? #f] #:overwrite [overwrite? #f])
  (define A0 A)
  (set! A (if overwrite? A (copy-flomat A)))
  (define jobvl (char->integer (if left?  #\V #\N)))
  (define jobvr (char->integer (if right? #\V #\N)))
  (define-param (m n a lda) A)
  (define WR (make-flomat n 1))
  (define WI (make-flomat n 1))
  (define wr (flomat-a WR))
  (define wi (flomat-a WI))
  (define VL (if left?  (make-flomat n n) (make-flomat 1 1)))
  (define VR (if right? (make-flomat n n) (make-flomat 1 1)))
  (define-param (mvl nvl vl ldvl) VL)
  (define-param (mvr nvr vr ldvr) VR)
  (define lwork (max 1 (* 4 n)))
  (define WORK  (make-flomat lwork 1))
  (define work  (flomat-a WORK))
  (define info (dgeev_ jobvl jobvr n a lda wr wi vl ldvl vr ldvr work lwork))
  ; (when (> info 0) (displayln "Warning: no convergence"))
  (values A0 WR WI VL VR info))

(define (real+imaginary->vector X Y)
  ; Convert two vectors of same length with real and imaginary
  ; part to a Racket vector of imaginary numbers.
  (define who 'real+imaginary->vector)
  (set! X (result-flcolumn X))
  (set! Y (result-flcolumn Y))
  (define-param (mx nx x ldx) X)
  (define-param (my ny y ldy) Y)
  (unless (and (= mx my) (= nx ny 1))
    (raise-argument-error 
       who  "The two inputs must be vectors of the same length" (list X Y)))
  (for/vector #:length mx
      ([xi (in-col X)]
       [yi (in-col Y)])
    (if (zero? yi)
        xi
        (make-rectangular xi yi))))





;;;
;;; MATRIX DECOMPOSITIONS
;;;

;;; Pivots 

(struct pivots (ps)) ; ps is a u32vector
;   ps[i]=j  <=>  row i and row j-1 is swapped
; Note: Fortran counts from 1 !

(define (unsafe-pivot-ref ps i)
  ; Fortran indices are 1-based.
  (- (u32vector-ref ps i) 1))

(define (pivots-ref Ps i)
  (unsafe-pivot-ref (pivots-ps Ps) i))

(define (pivots-length Ps)
  (u32vector-length (pivots-ps Ps)))

(define (pivots->flomat Ps)
  ; return the permuation matrix
  (define ps (pivots-ps Ps))
  (define k (u32vector-length ps))  
  (define A (make-flomat k k 0.0))
  (define-param (m n a lda) A)
  ; introduce ones on diagonal
  (for ([i (in-range m)])
    (unsafe-set! a lda i i 1.0))
  ; perform row permutations
  (for ([i (in-range (- m 1) -1 -1)])
    (define i* (unsafe-pivot-ref ps i))
    (unless (= i i*)
      (flomat-swap-rows! A i i*)))
  A)

(define (pivots-sign Ps)
  ; return the sign of the corresponding permuation
  (define ps (pivots-ps Ps))
  (define n (u32vector-length ps))
  (for/product ([i (in-range n)])
    (define i* (unsafe-pivot-ref ps i))
    (if (= i i*) 1 -1)))

;;;
;;; PLU Factorization
;;;

; A = P L U
; where P is a permutation matrix,
;       L is lower triangular
;       U is upper triangular.
; Note: U is the result of Gauss elimation.

(define-lapack dgetrf_ 
  ; http://www.netlib.org/lapack/double/dgetrf.f
  ; DGETRF computes an LU factorization of a general M-by-N matrix A
  ; using partial pivoting with row interchanges.
  ; The factorization has the form
  ;     A = P * L * U
  ; where P is a permutation matrix, L is lower triangular with unit
  ; diagonal elements (lower trapezoidal if m > n), and U is upper
  ; triangular (upper trapezoidal if m < n).
  
  ; Algorithm: Gaussian elimination with partial pivoting
  (_fun (m : (_ptr i _int))
        (n : (_ptr i _int))
        (a : _flomat)
        (lda : (_ptr i _int))
        (ipiv : (_u32vector o (ptr-ref m _int)))
        (info : (_ptr o _int)) 
        -> _void
        -> (values (pivots ipiv) info)))

(define (flomat-lu! A)
  (define-param (m n a lda) A)
  (dgetrf_ m n a lda))

(define (flomat-plu A)
  (define B (copy-flomat A))
  (define-values (ps info) (flomat-lu! B))
  (define P (pivots->flomat ps))
  (define L (flomat-extract-lower/ones B))
  (define U (flomat-extract-upper B))
  ; TODO: What to do with info?
  (values P L U))

(define (flomat-extract-upper A)
  ; extract the upper matrix, 
  ; including the diagonal
  ; discard below diagonal
  (define-param (m n) A)
  (define k (min m n))
  (define U (make-flomat k k))
  ; TODO: use unsafe-ref or unsafe-vector-copy
  (for* ([j (in-range 0 k)]
         [i (in-range (min (+ 1 j) k))])
    (flomat-set! U i j (flomat-ref A i j)))
  U)

(define (flomat-extract-lower/ones A)
  ; extract the lower matrix, 
  ; and insert ones on diagonal
  (define L (copy-flomat A))
  (define-param (m n) A)
  ; TODO: use unsafe-ref or unsafe-vector-copy
  (for* ([j (in-range n)]
         [i (in-range 0 j)])
    (flomat-set! L i j 0))
  (for* ([j (in-range (min m n))])
    (flomat-set! L j j 1.0))
  L)

(define (flomat-extract-lower A)
  ; extract the lower matrix including the diagonal
  (define L (copy-flomat A))
  (define-param (m n) A)
  ; TODO: use unsafe-ref or unsafe-vector-copy
  (for* ([j (in-range 1 n)]
         [i (in-range 0 j)])
    (flomat-set! L i j 0))
  L)

;;; SVD - Singular Value Decomposition

(define-lapack dgesvd_ 
  ; compute SVD 
  ; A = U * SIGMA * transpose(V)
  ; SIGMA is an mxn matrix, 
  ; Algorith: QR used
  (_fun (jobu  : (_ptr i _byte)) ; char: a, s, o or n
        (jobvt : (_ptr i _byte)) ; char
        (m     : (_ptr i _int))  ; rows in A
        (n     : (_ptr i _int))  ; cols in A
        (a     :  _flomat)       ; io
        (lda   : (_ptr i _int))
        (s     :  _flomat)       ; min(m,n) x 1
        (u     :  _flomat)       ; mxm if jobu = a
        (ldu   : (_ptr i _int)) 
        (vt    :  _flomat)       ; nxn if jobvt = a
        (ldvt  : (_ptr i _int))  ;         
        (work  :  _flomat)       ; dim max(1,lwork)
        (lwork : (_ptr i _int))  ; 
        (info  : (_ptr o _int))
        -> _void
        -> info))

(define-lapack dgesdd_ 
  ; compute SVD 
  ; A = U * SIGMA * transpose(V)
  ; SIGMA is an mxm matrix, 
  ; Algorithm: Divide and conquer with QR used for small
  ; This is the recommended algorithm, but uses
  ; more work space.
  (_fun (jobu  : (_ptr i _byte)) ; char: a, s, o or n
        (jobvt : (_ptr i _byte))
        (m     : (_ptr i _int))  ; rows in A
        (n     : (_ptr i _int))  ; cols in A
        (a     :  _flomat)       ; io
        (lda   : (_ptr i _int))
        (s     :  _flomat)       ; min(m,n) x 1
        (u     :  _flomat)       ; mxm if jobu = a
        (ldu   : (_ptr i _int)) 
        (vt    :  _flomat)       ; nxn if jobvt = a
        (ldvt  : (_ptr i _int))  ;         
        (work  :  _flomat)       ; dim max(1,lwork)
        (lwork : (_ptr i _int))   
        (info  : (_ptr o _int))
        -> _void
        -> info))

(define (flomat-svd! A)
  ; TODO: Use lwork=-1 to get size of work
  (define-param (m n a lda) A)
  (define superb (- (min m n) 1))
  (define U  (make-flomat m m))
  (define S  (make-flomat (min m n) 1))
  (define VT (make-flomat n n))
  (define u  (flomat-a U))
  (define s  (flomat-a S))
  (define vt (flomat-a VT))
  (define lwork (* 10 (max m n))) ; conservative estimate
  (define W  (make-flomat lwork lwork))
  (define w  (flomat-a W))
  (define ca (char->integer #\A))
  (define info (dgesvd_ ca ca m n a lda s u m vt n w lwork))
  ; ? TODO: Best way to return error ?
  ; (displayln (list 'info-from-svd info))
  ; (when (> info 0) (displayln "Warning: no convergence"))
  ; S is column vector of singular values
  ; Turn S into SIGMA (mxn) by placing the values of S on the diagonal.
  (values U S VT)) 

(define (flomat-svd A)
  (flomat-svd! (copy-flomat A)))

(define (flomat-singular-values A)
  (define B (copy-flomat A))
  (define-values (S V D) (flomat-svd! B))
  (flomat->vector V))

(define (flomat-diagonal-from-singular-values m n S [reciproc? #f])
  (define A (flomat-zeros m n))
  (for ([i (in-range m)]
        [s (in-col S)])
    (flomat-set! A i i (if reciproc? (/ 1 s) s)))
  A)

(define (flomat-pseudo-inverse A)
  (define-param (m n) A)
  (define-values (U S VT) (flomat-svd A))
  (define Σ  (flomat-diagonal-from-singular-values m n S #t))
  (define A+ (flomat* (flomat-transpose VT) (flomat* (transpose Σ) (flomat-transpose U))))
  A+)


;;; QR Factorization
; dgeqrfp returns positive entries on the diagonal
; for some reason this is missing on macOS, so now dgeqrf is used instead
#;(define-lapack dgeqrfp_ 
  ; Compute A = Q*R  
  ; Use dorgqr to generate matrix from output
  (_fun (m : (_ptr i _int)) ; rows in A
        (n : (_ptr i _int)) ; cols in A
        (a : _flomat) ; io
        (lda : (_ptr i _int))
        (tau : _flomat) ; min(m,n)x1        
        (work : _flomat) ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int)) ; >=max(1,n) best with >=n * blocksize
        (info : (_ptr o _int))  ; 
        -> _void
        -> info))

(define-lapack dgeqrf_
  ; Compute A = Q*R  
  ; Use dorgqr to generate matrix from output
  (_fun (m     : (_ptr i _int))  ; rows in A
        (n     : (_ptr i _int))  ; cols in A
        (a     :  _flomat)       ; io
        (lda   : (_ptr i _int))
        (tau   :  _flomat)       ; min(m,n)x1        
        (work  :  _flomat)       ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int))  ; >=max(1,n) best with >=n * blocksize
        (info  : (_ptr o _int))   
        -> _void
        -> info))


(define-lapack dorgqr_
  ; generate matrix from output of dgeqrf
  (_fun (m : (_ptr i _int)) ; rows in Q
        (n : (_ptr i _int)) ; cols in Q m>=n>=0
        (k : (_ptr i _int)) ; number of reflectors
        (a : _flomat) ; io
        (lda : (_ptr i _int))
        (tau : _flomat)       ; min(m,n)x1        
        (work : _flomat)      ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int)) ; >=max(1,n) best with >=n * blocksize
        (info : (_ptr o _int))  ; 
        -> _void
        -> info))


(define (flomat-qr B)
  (define A (copy-flomat B))
  (define-param (m n a lda) A)
  (define k (min m n))
  (define tau (make-flomat k 1))
  (define atau (flomat-a tau))
  ; first call dgeqrf_ to get a working size
  (define work0 (make-flomat 1 1))
  (define info0 (dgeqrf_ m n a lda atau (flomat-a work0) -1))
  (define lwork (inexact->exact (flomat-ref work0 0 0))) ; 64 is a typical value
  ; now make the real call
  (define work  (make-flomat lwork 1))  
  (define awork (flomat-a work))  
  (define info (dgeqrf_ m n a lda atau awork lwork))
  (define R (flomat-extract-upper A))
  (define info1 (dorgqr_ m n k a lda atau awork lwork))  
  ; ? TODO: what to do with info  
  (values A R))

; old version used dgeqrfp
#;(define (flomat-qr B)
  (define A (copy-flomat B))
  (define-param (m n a lda) A)
  (define k (min m n))
  (define tau (make-flomat k k))
  (define atau (flomat-a tau))
  (define lwork (* 64 n)) ; 64 block size guess
  (define work (make-flomat lwork 1))
  (define awork (flomat-a work))
  ; TODO: Use lwork=-1 to get optimal lwork size
  (define info (dgeqrf_ m n a lda atau awork lwork))
  (define R (flomat-extract-upper A))
  (define info1 (dorgqr_ m n k a lda atau awork lwork))  
  ; ? TODO: what to do with info
  (values A R))

;;; 
;;; INVERSE
;;;

(define-lapack dgetri_
  ; http://www.netlib.org/lapack/double/dgetri.f
  ; DGETRI computes the inverse of a matrix using the LU factorization
  ; computed by DGETRF.
  ; This method inverts U and then computes inv(A) by solving the system
  ;   inv(A)*L = inv(U) for inv(A).
  (_fun (n     : (_ptr i _int))
        (a     :  _flomat)
        (lda   : (_ptr i _int))
        (ipiv  :  _u32vector)
        (work  : (_or-null _flomat)) ; output
        (lwork : (_ptr i _int))
        (info  : (_ptr o _int))
        -> _void
        -> (values info work)))

(define (flomat-inverse! A)
  (check-square 'flomat-inverse! A)
  ; TODO: this works, but call dgetri with lwork=-1
  ;       to get optimal size of workspace in first
  ;       entry of the work array.
  (define-param (m n a lda) A)
  (define work (copy-flomat A))
  (define-values (ipiv info) (flomat-lu! A))
  (dgetri_ m a lda (pivots-ps ipiv) (flomat-a work) (* m m))
  A)

(define (flomat-inverse A)
  (flomat-inverse! (copy-flomat A)))


;;;
;;; Cholesky Factorization
;;;

; A = U^T U or A = L L^T
; where L is lower triangular
;       U is upper triangular.
;  and  A is a square, real symmetric, positive definite matrix A.

(define-lapack dpotrf_ 
  ; http://www.netlib.org/lapack/double/dpotrf.f
  ;  DPOTRF computes the Cholesky factorization of a real symmetric
  ;  positive definite matrix A.
  ;
  ;  The factorization has the form
  ;     A = U^T * U     if UPLO = 'U', or
  ;     A = L   * L^T   if UPLO = 'L',
  ;  where U is an upper triangular matrix and L is lower triangular.
  (_fun (uplo : (_ptr i _byte))  ; 'U' or 'L'
        (n    : (_ptr i _int))
        (a    : _flomat)
        (lda  : (_ptr i _int))
        (info : (_ptr o _int)) 
        -> _void
        -> info))

(define ascii-U (char->integer #\U))
(define ascii-L (char->integer #\L))

(define (flomat-cholesky! A [upper? #f])
  (check-square flomat-cholesky! A)
  (define-param (m n a lda) A)
  (define uplo (if upper? ascii-U ascii-L))
  (define info (dpotrf_ uplo m a lda))
  info)

(define (flomat-cholesky A [upper? #f])
  (define B    (copy-flomat A))
  (define info (flomat-cholesky! B upper?))
  (if upper?
      (flomat-extract-upper B)
      (flomat-extract-lower B)))

;;;
;;; INVARIANTS
;;; 

(define (flomat-trace A)
  (check-square 'matrix-trace A)
  (for/sum ([i (in-range (flomat-m A))])
    (flomat-ref A i i)))

(define (flomat-determinant-from-plu LU pivots)
  ; compute determinant using output from PLU
  ; factorization
  (* (pivots-sign pivots)
     (for/product ([i (in-range (flomat-m LU))])
       (flomat-ref LU i i))))

(define (flomat-determinant A)
  (check-square 'matrix-determinant A)
  (define LU (copy-flomat A))
  (define-values (pivots info) (flomat-lu! LU))
  (flomat-determinant-from-plu LU pivots))

(define (count-columns-without-pivot pivots)
  ; TODO: Does this strategy work?
  (define ps (pivots-ps pivots))
  (define m (u32vector-length ps))
  (define with 
    (for/sum ([i (in-range m)])
      (define i* (- (u32vector-ref ps i) 1))
      (if (= i i*) 0 1)))
  (- m with))

(define (flomat-rank A)
  ; See answer: http://scicomp.stackexchange.com/questions/1861/understanding-how-numpy-does-svd
  ; rank = dimension of column space = dimension of row space  
  ;      = number of non-zero singular values
  (define-values (U Σ VT) (flomat-svd A))
  ; ? TODO: check info from -svd...
  (for/sum ([i (in-range (flomat-m Σ))]
            ; TODO: Which value for epsilon is correct?
            #:unless (< (abs (flomat-ref Σ i 0)) epsilon))
    1))

(define (flomat-nullity A)
  ; nullity = dimension of null space
  (define-param (m n) A)
  (- n (flomat-rank A)))

;;;
;;; VECTOR OPERATIONS
;;;

; Column vectors are represented as mx1 matrices.
; All operations working on column vectors accept
; standard vectors as input. Outputs are always
; in the form of a mx1 matrix.

(define (vector->flcolumn v)
  (define m (vector-length v))
  (vector->flomat m 1 v))

(define (vector->flrow v)
  (define n (vector-length v))
  (vector->flomat 1 n v))

(define (list->flcolumn xs)
  (vector->flcolumn (list->vector xs)))

(define (result-flcolumn c)
  ; convert output to mx1 matrix
  (if (vector? c)
      (vector->flcolumn c)
      (if (list? c)
          (list->flcolumn c)
          c)))

(define (flcolumn . xs)
  ; TODO: skip intermediary vector
  (vector->flcolumn 
   (list->vector xs)))

(define (flcolumn-size v)
  (if (vector? v)
      (vector-length v)
      (flomat-m v)))

;;; Dot Product

(define-cblas* cblas_xdot _x (s d)
  ; dot product, vectors
  (_fun (n : _int) 
        (X : _flomat) (incX : _int)
        (Y : _flomat) (incY : _int)
        -> _x))

(define (unsafe-vector-product n x y)
  (cblas_ddot n x 1 y 1))

(define (flcolumn-dot X Y)
  (set! X (result-flcolumn X))
  (set! Y (result-flcolumn Y))
  (define-param (m _ x ldx) X)
  (define-param (s __ y ldy) Y)
  (unless (= m s) 
    (error 
     'column-dot 
     "expected two mx1 matrices with same number of rows, got ~a and ~a"
     X Y))
  (unsafe-vector-product m x y))

(define fldot flcolumn-dot)

(define (flcolumn-norm v)
  (define-param (m _ a lda) (result-flcolumn v))
  (cblas_dnrm2 m a 1))

(define (flcolumn-unit m i)
  ; return i'th unit vector 
  (define U (make-flomat m 1 0.0))
  (flomat-set! U i 0 1.0)
  U)

(define (flscale-column s A)
  (define s* (real->double-flonum s))
  (cond
    [(vector? A)
     (define m (vector-length A))
     (vector->flcolumn 
      (for/vector #:length m 
        ([i (in-range m)])
        (* s* (vector-ref A i))))]
    [else
     (flomat-scale s A)]))

(define (flcolumn+ v w)
  (define m (flcolumn-size v))
  (define n (flcolumn-size w))
  (unless (= m n)
    (error 
     'flcolumn+ 
     "expected two column vectors of the same length, got ~a and ~a" v w))
  (cond 
    [(and (vector? v) (vector? w))
     (vector->flcolumn
      (for/vector #:length (+ m n) 
        ([i (in-range 0 m)]
         [x (in-vector v)]
         [y (in-vector w)])
        (+ x y)))]
    [else          
     (flomat+ (result-flcolumn v) (result-flcolumn w))]))

(define (flcolumn-projection v w)  
  ; Return the projection og vector v on vector w.
  (let ([w.w (fldot w w)])
    (if (zero? w.w)
        (error 'flcolumn-projection "projection on the zero vector not defined")
        (flscale-column (/ (fldot v w) w.w) w))))

(define (flcolumn-projection-on-unit v w)
  ; Return the projection of vector v on a unit vector w.  
  (flscale-column (flcolumn-dot v w) w))

(define (flcolumn-normalize w)
  ; Return unit vector with same direction as v.
  ; If v is the zero vector, the zero vector is returned.
  (define norm (flcolumn-norm w))
  (cond [(zero? norm) w]
        [else (flscale-column (/ norm) w)]))

(define (flzero-column-vector? v [eps #f])
  (define val (flomat-max-abs-value (result-flcolumn v)))
  (if eps (< (abs val) eps) (zero? val)))
  
; (flprojection-on-orthogonal-basis v bs)
;     Project the vector v on the orthogonal basis vectors in bs.
;     The basis bs must be either the column vectors of a matrix
;     or a sequence of column-vectors.
(define (flprojection-on-orthogonal-basis v bs)
  (if (empty? bs)
      (error 'flprojection-on-orthogonal-basis 
             "received empty list of basis vectors")
      (for/flomat-sum([b (in-list bs)])
                       (flcolumn-projection v (result-flcolumn b)))))

;     Project the vector v on the orthonormal basis vectors in bs.
;     The basis bs must be either the column vectors of a matrix
;     or a sequence of column-vectors.
(define (flprojection-on-orthonormal-basis v bs)
  (for/flomat-sum 
   ([b bs]) 
   (flomat-scale (flcolumn-dot v b) b)))
  

; (flgram-schmidt-orthogonal ws)
;     Given a list ws of flcolumn vectors, produce 
;     an orthogonal basis for the span of the
;     vectors in ws.
(define (flgram-schmidt-orthogonal ws1)
  (define ws (map result-flcolumn ws1))
  (cond 
    [(null? ws)       '()]
    [(null? (cdr ws)) (list (car ws))]
    [else 
     (define (loop vs ws)
       (cond [(null? ws) vs]
             [else
              (define w (car ws))
              (let ([w-proj (flprojection-on-orthogonal-basis w vs)])
                ; Note: We project onto vs (not on the original ws)
                ;       in order to get numerical stability.
                (let ([w-minus-proj (flomat- w w-proj)])
                  (if (flzero-column-vector? w-minus-proj)
                      (loop vs (cdr ws)) ; w in span{vs} => omit it
                      (loop (cons (flomat- w w-proj) vs) (cdr ws)))))]))
     (reverse (loop (list (car ws)) (cdr ws)))]))

; (flgram-schmidt-orthonormal ws)
;     Given a list ws of flcolumn vectors, produce 
;     an orthonormal basis for the span of the
;     vectors in ws.
(define (flgram-schmidt-orthonormal ws)
  (map flcolumn-normalize
       (flgram-schmidt-orthogonal ws)))

; (flprojection-on-subspace v ws)
;  Returns the projection of v on span{w_i}, w_i in ws.
(define (flprojection-on-subspace v ws)
  (flprojection-on-orthogonal-basis
   v (flgram-schmidt-orthogonal ws)))

;;;
;;; SUMS
;;;


(define unsafe-sum
  ;; (unsafe-sum n a lda)
  ;  Beginning from addrees a, sum every lda element in double array.
  ;
  ; There is no `sum` operation in BLAS or LAPCK, so we compute the dot product
  ; between the vector and a vector consisting only of ones.
  ; Also we don't need to allocate a vector with only ones, we simply use ld=0.  
  (let ()
    (define one (cast (f64vector->cpointer (f64vector 1.0)) _pointer _flomat))
    (λ (n a lda)
      (cblas_ddot n a lda one 0))))

(define (flomat-column-sum A j)
  ; sum of all entries in column j
  (define who 'flomat-column-sum)
  (check-flomat who A)
  (check-legal-column who j A)  
  (define-param (m n a lda) A)  
  (define aj (ptr-elm a lda 0 j)) ; address of j'th column
  (unsafe-sum m aj 1))

(define (flomat-row-sum A i)
  ; sum of all entries in row i
  (define who 'flomat-row-sum)
  (check-flomat who A)
  (check-legal-row who i A)  
  (define-param (m n a lda) A)  
  (define ai (ptr-elm a lda i 0))
  (unsafe-sum n ai lda))

(define (flomat-column-sums A)
  ; row vector of all column sums
  (check-flomat 'flomat-column-sums A)
  (define n (ncols A))
  (for/flomat 1 n ([j (in-range n)])
              (flomat-column-sum A j)))

(define (flomat-row-sums A)
  ; column vector of all row sums
  (check-flomat 'flomat-row-sums A)
  (define m (nrows A))
  (for/flomat m 1 ([i (in-range m)])
    (flomat-row-sum A i)))

;;;
;;; EQUATION SOLVING
;;;

(define-lapack dgesv_ ; Double, GEneral, Solve ...
  ; Compute solution to AX=B, where
  ; A is nxn and X and B are n x nrhs
  (_fun (n : (_ptr i _int))
        (nrhs : (_ptr i _int))
        (a : _flomat) ; io
        (lda : (_ptr i _int))
        (ipiv : (_u32vector o (ptr-ref n _int)))
        (b : _flomat) ; io
        (ldb : (_ptr i _int))
        (info : (_ptr o _int))
        -> _void
        -> info))

(define (flomat-solve A b)
  ; A matrix, b flcolumn
  ; called mldivide aka  matrix left divide in Matlab  
  (define-param (m n a lda) (copy-flomat A))
  (define bout (copy-flomat (result-flcolumn b)))
  (define info (dgesv_ n 1 a lda (flomat-a bout) m))
  ; ? TODO Handle info
  bout)

(define (flomat-solve-many! A B)  
  ; A matrix, b flcolumn  
  ; A and B are overwritten
  (define-param (m n a lda) A)
  (define-param (_ nrhs b ldb) B)
  (define info (dgesv_ n nrhs a lda b ldb))
  ; (displayln (list 'flomat-solve-many! info))
  ; ? TODO: handle info
  (values B))

(define (flomat-solve-many A bs-or-B)  
  ; A matrix, b flcolumn 
  (define-param (m n) A)
  (define B (if (list? bs-or-B)    
                (apply flomat-augment 
                       (map result-flcolumn bs-or-B))
                (copy-flomat bs-or-B)))
  (flomat-solve-many! (copy-flomat A) B))

(define (flomat->columns A)
  (define-param (m n) A)  
  (for/list ([j (in-range n)])
    (flomat-column A j)))

(define flomat-left-divide flomat-solve-many)


;;;
;;; LINEAR LEAST SQUARES PROBLEM
;;;

(define-lapack dgelsd_ ; Double, GEneral, Least Square 
  ; minimize 2-norm(| b - A*x |)  where b and x are columns in B and X
  ; See http://www.netlib.org/lapack/lapack-3.1.1/html/dgelsd.f.html
  ; B is overwritten with the result.
  (_fun (m     : (_ptr i _int))     ; m>=0 rows of A
        (n     : (_ptr i _int))     ; n>=0 cols of A
        (nrhs  : (_ptr i _int))     ; number of rhs (number of cols in B and X)
        (a     :  _flomat) ; io   ; is overwritten
        (lda   : (_ptr i _int))     ; lda >= max(1,m)
        (b     :  _flomat) ; io   ; mxnrhs 
        (ldb   : (_ptr i _int))     ; ldb >= max(1,max(m,n))
        (s     :  _flomat)        ; min(m,n) x 1  singular values of A in decreasing order
        ;                           ; the condition number of A is 2-norm of S(1)/S(min(m,n)
        (rcond : (_ptr i _double))  ; singular values S(i)<=rcond*S(1) is treated as 0.
        ;                           ; if rcond<0 then machine precision is used instead
        (rank  : (_ptr o _int))     ; effective rank of A (number of non-zero singular values 
        ;                           ; greater than rcond*s(1)
        (work  :  _flomat)        ; dim max(1,lwork) x1
        (lwork : (_ptr i _int))     ; dimension of work - use work query to find proper size
        (iwork : _pointer)          ; array of int   dim max(1,liwork)x1
        ;                           ;  LIWORK >= 3 * MINMN * NLVL + 11 * MINMN,
        ;                           ;  where MINMN = MIN( M,N ).
        (info : (_ptr o _int))      ; =0 succes exit, >0 svd failed to converge
        -> _void
        -> (values info rank)))

(define (lstsq A B)
  (define-param (mb nb aB ldB ) B) ; nb = nrhs
  ; Least squares solution to |Ax-b| minimizing |x| for each column in B.
  ; The common case is m>=n and rank(A)=n in which we get a solution to an overdetermined system.
  ; If m<n and rank(A)=m there are an infinite number of solutions.
  ; This routine will then find the solution x with minimal norm.

  ; We assume that A has full rank (that is rank min(m,n)).
  ;   A is mxn
  ;   b is mx1  B is mxk
  ;   x is nx1  X is nxk
  ; If the solution x is taller than b (that is n>m) then we need
  ; insert some zeros below B.

  ; In the text the case n>m one sees the driver dgelsd expects B
  ; to have ldb>=max(1,m,n) normally ldb>=rows in B = m ,
  ; so if n<m then B must be enlarged with extra rows.
  
  ; 1. Copy A, the driver will/can destroy A
  (define A0 (copy-flomat A))
  (define-param (m   n   a  lda)  A0)
  ; 2. Enlarge B into X
  (define X  (cond [(>= nb m) (copy-flomat B)]
                   [else      (make-flomat m nb)]))
  (define-param (_ nrhs b ldb) X)
  (when (< nb m)
    ; The result vectors are stored in B (or the copy of B).
    ; If the result vectors are B was enlarged, so we need to copy B into X
    (unsafe-matrix-copy! mb nb aB ldB b ldb))
  ; (displayln (list ldb (max m n) (= ldb (max m n))))
  ; 3. Allocate array for the singular values
  (define S (make-flomat (min m n) 1))
  (define s (flomat-a S))
  ; 4. The condition number -1.0 indicates machine precision.
  ;    Note: Numpy uses machine precision times mxn.
  (define rcond -1.0)
  ; 5. The size of the working are is determined by calling the driver dgelsd with lwork=-1.
  (define work0  (make-flomat 1 1)) ; first entry will contain optimal length of work
  (define awork0 (flomat-a work0))  
  (define lwork0 -1)  
  (define iwork0 (malloc 1 _int 'atomic))
  ; (displayln  (list m n nrhs a lda b ldb s rcond awork0 lwork0 iwork0))
  (define-values (__ ___) ; 
    (dgelsd_ m n nrhs a lda b ldb s rcond awork0 lwork0 iwork0))  
  (define lwork (inexact->exact (flomat-ref work0 0 0)))
  ; (displayln (list 'lwork lwork))
  ; 6. Prepare the actual call 
  ;(define nlvl  (max 0 (+ 1 (ceiling (log2 (/ (min m n) (+ SMLSIZ 1)))))))
  (define nlvl 32) ; 
  (define liwork (max 1 (+ (* 3 (min m n) nlvl) (* 11 (min m n)))))
  (define iwork  (malloc liwork _int 'atomic))    ; integer array 
  (define WORK   (make-flomat (max 1 lwork) 1))
  (define work   (flomat-a WORK))
  ;(displayln (list m n nrhs a lda b ldb s rcond work lwork iwork))
  (define-values (info rank) 
    (dgelsd_ m n nrhs a lda b ldb s rcond work lwork iwork))
  (shared-submatrix! X 0 0 n nrhs))



(define (linear-fit/plain xs ys)
  (set! xs (result-flcolumn xs))
  (set! ys (result-flcolumn ys))
  ; (check-same-length 'linear-fit xs ys)
  ; Here xs and ys are column vectors.

  ; Find a and b such that y=ax+b is a good fit to the input data.
  ; Minimize S = sum( (y - (ax+b))² )
  
  ; In vector form, the models is:
  ;  Y = X β + ε ; where β₁= a and β₂=b
  ;  X is an mx2 matrix where each is on the form  [x_i 1] (X is the design matrix)
  ; Multiply with X^T 
  ;  X^T Y = X^T X β + X^T ε
  ; Now left divide with X^T X which is a square matrix

  (define-param (mx nx) xs)
  (define-param (my ny) ys)  
  (define X      (flomat-augment xs (flomat-ones mx 1)))
  (define XT     (flomat-transpose X))
  (define XTX    (flomat* XT X))
  (define XTy    (flomat* XT ys))
  (define β      (flomat-left-divide XTX XTy))
  ; (define resids (flomat- ys (flomat* X β)))
  β)

(define (linear-fit/qr xs ys)
  (set! xs (result-flcolumn xs))
  (set! ys (result-flcolumn ys))
  ; (check-same-length 'linear-fit xs ys)
  ; Here xs and ys are column vectors.

  ; Find a and b such that y=ax+b is a good fit to the input data.
  ; Minimize S = sum( (y - (ax+b))² )
  
  ; In vector form, the models is:
  ;  Y = X β + ε ; where β₁= a and β₂=b
  ;  X is an mx2 matrix where each is on the form  [x_i 1] (X is the design matrix)

  ; With X = QR:
  ;   Q^T Y = Q^T X β + Q^T ε ~ Q^T X β = Q^T Q R β = R β
  ; Now solve for β.
  
  (define-param (mx nx) xs)
  (define-param (my ny) ys)  
  (define X      (flomat-augment xs (flomat-ones mx 1)))
  ; 1. Compute X = QR
  (define-values (Q R) (flomat-qr X))
  ; 2. Compute the vector Q^T Y
  (define QTY (flomat* Q ys #f 1.0 1.0 #t))
  ; 3. Solve upper triangular Rβ = Q^T Y for β.
  (define β (flomat-left-divide R QTY))
  ; (define resids (flomat- ys (flomat* X β)))
  β)

(define (linear-fit-residuals xs ys β)
  (set! xs (result-flcolumn xs))
  (set! ys (result-flcolumn ys))
  (define-param (mx nx) xs)
  (define X      (flomat-augment xs (flomat-ones mx 1)))
  (flomat- ys (flomat* X β)))

(define (values->list f)
  (call-with-values f list))

(define (linear-fit xs ys [method 'qr])
  (match method
    ['qr    (linear-fit/qr    xs ys)]
    ['plain (linear-fit/plain xs ys)]
    [_ (error linear-fit (~a "unknown method, given: " method))]))
  


;;;
;;; SEQUENCES
;;; 

(define (in-row/proc A r)
  (define-param (m n a lda) A)
  (make-do-sequence
   (λ ()
     (define (pos->elm j) (unsafe-ref a lda r j))
     (define next-pos add1)
     (define initial-pos 0)
     (define (continue? j) (< j n))
     (values pos->elm initial-pos continue? #f #f))))

; (in-row M i]
;     Returns a sequence of all elements of row i,
;     that is xi0, xi1, xi2, ...
(define-sequence-syntax in-row
  (λ () #'in-row/proc)
  (λ (stx)
    (syntax-case stx ()
      [[(x) (_ M-expr r-expr)]
       #'((x)
          (:do-in
           ([(M r m n a lda)
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 r-expr rd cd a lda))])
           (begin 
             (unless (flomat? M) 
               (raise-type-error 'in-row "expected flomat, got ~a" M))             
             (unless (and (integer? r) (and (<= 0 r ) (< r m))) 
               (raise-type-error 'in-row "expected row number" r)))
           ([j 0])
           (< j n)
           ([(x) (unsafe-ref a lda r j)])
           #true
           #true
           [(+ j 1)]))]
      [[(i x) (_ M-expr r-expr)]
       #'((i x)
          (:do-in
           ([(M r m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 r-expr rd cd a lda))])
           (begin 
             (unless (flomat? M) 
               (raise-type-error 'in-row "expected flomat, got ~a" M))             
             (unless (and (integer? r) (and (<= 0 r ) (< r m))) 
               (raise-type-error 'in-row "expected row number" r)))
           ([j 0])
           (< j n)
           ([(x) (unsafe-ref a lda r j)]
            [(i) j])
           #true
           #true
           [(+ j 1)]))]
      [[_ clause] (raise-syntax-error 
                   'in-row "expected (in-row <flomat> <row>)" #'clause #'clause)])))

; (in-flcol M j]
;     Returns a sequence of all elements of column j,
;     that is x0j, x1j, x2j, ...

(define (in-col/proc A s)
  (define-param (m n a lda) A)
  (make-do-sequence
   (λ ()
     (define (pos->elm i) (unsafe-ref a lda i s))
     (define next-pos add1)
     (define initial-pos 0)
     (define (continue? i) (< i m))
     (values pos->elm next-pos initial-pos #f #f))))

(define-sequence-syntax in-col
  (λ () #'in-col/proc)
  (λ (stx)
    (syntax-case stx ()
      ; M-expr evaluates to column
      [[(x) (_ M-expr)]
       #'((x)
          (:do-in
           ([(M n m a) 
             (let ([M1 (result-flcolumn M-expr)])
               (define-param (rd cd a) M1)
               (values M1 rd cd a))])
           (unless (flomat? M) 
             (raise-type-error 'in-column "expected matrix, got ~a" M))
           ([j 0])
           (< j n)
           ([(x) (ptr-ref a _double j)])
           #true
           #true
           [(+ j 1)]))]
      ; M-expr evaluates to matrix, s-expr to column index
      [[(x) (_ M-expr s-expr)]
       #'((x)
          (:do-in
           ([(M s m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 s-expr rd cd a lda))])
           (begin 
             (unless (flomat? M) 
               (raise-type-error 'in-col "expected matrix, got ~a" M))
             (unless (integer? s) 
               (raise-type-error 'in-col "expected column number, got ~a" s))
             (unless (and (integer? s) (and (<= 0 s ) (< s n))) 
               (raise-type-error 'in-col "expected column number, got ~a" s)))
           ([j 0])
           (< j m)
           ([(x) (unsafe-ref a lda j s)])
           #true
           #true
           [(+ j 1)]))]
      [[(i x) (_ M-expr s-expr)]
       #'((x)
          (:do-in
           ([(M s m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 s-expr rd cd a lda))])
           (begin 
             (unless (flomat? M) 
               (raise-type-error 'in-column "expected matrix, got ~a" M))
             (unless (integer? s) 
               (raise-type-error 'in-column "expected col number, got ~a" s))
             (unless (and (integer? s) (and (<= 0 s ) (< s n))) 
               (raise-type-error 'in-column "expected col number, got ~a" s)))
           ([j 0])
           (< j m)
           ([(x) (unsafe-ref a lda j s)]
            [(i) j])
           #true
           #true
           [(+ j 1)]))]
      [[_ clause] (raise-syntax-error 
                   'in-col "expected (in-col <flomat> <column>)" #'clause #'clause)])))


;;;
;;; Special Matrices
;;;

; Note: See Matlab gallery for ideas.

(define (flomat-vandermonde xs n)
  ; One row for each element of xs.
  ; Each row consist of the first 0..n-1 powers of x.
  (define m (length xs))
  (define αs (list->vector xs))
  (define α^j (make-vector m 1.0))
  (for*/flomat m n #:column
                 ([j (in-range 0 n)]
                  [i (in-range 0 m)])
                 (define αi^j (vector-ref α^j i))
                 (define αi   (vector-ref αs i ))
                 (vector-set! α^j i (* αi^j αi))
                 αi^j))


;;;
;;; SYNTAX
;;;

(require
 (for-syntax racket/base
             syntax/parse))

(define-syntax (flomat: stx)
  (syntax-parse stx 
    [(_ [[x0 xs0 ...] [x xs ...] ...])
     (syntax/loc stx (vectors->flomat (vector (vector x0 xs0 ...) (vector x xs ...) ...)))]
    [(_ [xs ... (~and [] r) ys ...])
     (raise-syntax-error 'flomat: "given empty row" stx #'r)]
    [(_ (~and [] c))
     (raise-syntax-error 'flomat: "given empty matrix" stx #'c)]
    [(_ x)
     (raise-syntax-error 'flomat: "expected two-dimensional data" stx)]))

(define-syntax (flrow-matrix stx)
  (syntax-parse stx
    [(_ [x xs ...]) (syntax/loc stx (flomat: [[x xs ...]]))]
    [(_ (~and [] r))
     (raise-syntax-error 'flrow-matrix "given empty row" stx #'r)]))

(define-syntax (flcol-matrix stx)
  (syntax-parse stx 
    [(_ [x xs ...])      (syntax/loc stx (flomat: [[x] [xs] ...]))]
    [(_ (~and [] c))
     (raise-syntax-error 'flrow-matrix "given empty column" stx #'c)]))


; TODO:
#;(provide 
   ; DONE matrix*
   ; DONE matrix-expt
   ; DONE matrix-ref
   ; DONE matrix-scale
   ; DONE matrix-row-vector?
   ; DONE matrix-column-vector?
   ; DONE matrix/dim     ; construct
   ; DONE matrix-augment ; horizontally
   ; DONE matrix-stack   ; vertically
   ; DONE matrix-block-diagonal
   ; norms
   ; DONE matrix-norm
   ; operators
   ; DONE matrix-transpose
   ; NO-COMPLEX matrix-conjugate
   ; NO-COMPLEX matrix-hermitian
   ; DONE matrix-inverse
   ; row and column
   ; DONE matrix-scale-row
   ; DONE matrix-scale-column
   ; DONE matrix-swap-rows
   ; DONE matrix-swap-columns
   ; DONE matrix-add-scaled-row
   ; DONE ADDED matrix-add-scaled-column
   ; reduction
   ; (DONE) matrix-gauss-eliminate          ; ? use upper in LU ?
   ; DONE matrix-gauss-jordan-eliminate     ; ? LU ? Use matrix-gauss-eliminate
   ; (DONE) matrix-row-echelon-form         ; ? LU ?
   ; DONE matrix-reduced-row-echelon-form    ; ? LU ? Use matrix-gauss-eliminate
   
   ; invariant
   ; DONE matrix-rank  (uses SVD!)
   ; DONE matrix-nullity
   ; DONE matrix-determinant
   ; DONE matrix-trace
   ; spaces
   ;matrix-column+null-space
   ; solvers
   ; DONE matrix-solve
   ; DONE matrix-solve-many
   ; spaces
   matrix-column-space  ; use SVD somehow
   ; column vectors
   ; DONE column        ; construct
   ; DONE unit-column
   ; DONE result-column ; convert to lazy
   ; DONE column-dimension
   ; DONE column-dot
   ; DONE column-norm
   ; DONE column-projection
   ; DONE column-normalize 
   ; DONE scale-column
   ; DONE column+
   ; projection
   ; DONE projection-on-orthogonal-basis
   ; DONE projection-on-orthonormal-basis
   ; DONE projection-on-subspace
   ; DONE gram-schmidt-orthogonal
   ; DONE gram-schmidt-orthonormal
   ; factorization
   ; DONE matrix-lu (renamed to matrix-plu)
   ; DONE matrix-qr
   ; comprehensions
   ; DONE for/matrix:
   ; DONE for*/matrix:
   ; DONE for/matrix-sum:
   ; sequences
   ; DONE in-row
   ; DONE in-column
   ; special matrices
   ; DONE vandermonde-matrix
   )


(define (flomat->lists A)
  (map vector->list
       (vector->list 
        (flomat->vectors A))))

(define (lists->flomat xss)
  (vectors->flomat
   (list->vector (map list->vector xss))))


(define (flomat-map! A f)
  (define-param (m n a lda) A)
  (for* ([i (in-range m)]
         [j (in-range n)])
    (define aij (unsafe-ref a lda i j))
    (define x   (f aij))
    (define x*  (real->double-flonum x))
    (unsafe-set! a lda i j x*))
  A)


(define (flomat-map A f)
  (flomat-map! (copy-flomat A) f))

(define (flomat-make-diagonal v [k 0])
  ; create a square matrix A with diagonal elements from v
  ; if k is given, place the elements on the kth digonal,
  ; k=0 is the main diagonal
  ; k>0 above the main diagonal
  ; k<0 below the main diagonal
  (check-vector  'flomat-make-diagonal v)
  (check-integer 'flomat-make-diagonal k)
  (define s (+ (vector-length v) (abs k)))
  (define A (make-flomat s s))
  (define-param (m n a lda) A)
  (cond
    [(= k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda i i x*))]
    [(> k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda i (+ i k) x*))]
    [(< k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda (- i k) i x*))])
  A)

(define (flomat-eye m0 [n0 m0] [k 0])
  (check-integer 'flomat-eye k)
  ; create an mxn matrix with 1 on the k'th diagonal
  (define A (make-flomat m0 n0 0.0))
  (define-param (m n a lda) A)
  (cond
    [(= k 0) (for ([i (in-range m)])
               (unsafe-set! a lda i i 1.0))]
    [(> k 0) (for ([i (in-range (- m k))])
               (unsafe-set! a lda i (+ i k) 1.0))]
    [(< k 0) (for ([i (in-range (+ m k))])
               (unsafe-set! a lda (- i k) i 1.0))])
  A)

(define (flomat-diagonal A [k 0])
  ; extract the k'th diagonal of the matrix A
  (check-flomat 'flomat-diagonal A)
  (check-integer  'flomat-diagonal k)
  (define-param (m n a lda) A)
  (define s (min m n))
  (cond
    [(= k 0) (for/vector ([i (in-range s)])
               (unsafe-ref a lda i i))]
    [(> k 0) (for/vector ([i (in-range (- s k))])
               (unsafe-ref a lda i (+ i k)))]
    [(< k 0) (for/vector ([i (in-range (+ s k))])
               (unsafe-ref a lda i (- i k)))]))

(define (flomat-lower-triangle A [k 0])
  ; return triangle with elements on or below the k'th diagonal
  (check-flomat 'flomat-lower-triangle A)
  (check-integer  'flomat-lower-triangle k)
  (define B (copy-flomat A))
  (define-param (m n a lda) A)
  (cond
    [(= k 0)  (for* ([i (in-range m)]
                     [j (in-range (+ i 1) n)])
                (flomat-set! B i j 0.0))]
    [(> k 0)  (for* ([i (in-range m)]
                     [j (in-range (+ i 1 k) n)])
                (flomat-set! B i j 0.0))]
    [(< k 0)  (for* ([i (in-range m)]
                     [j (in-range (max 0 (+ i 1 k)) n)])
                (flomat-set! B i j 0.0))])
  B)

(define (flomat-circulant-matrix v)
  ; A circulant matrix is a matrix in which each row
  ; is the previous row shifted one to the right.
  (define n (vector-length v))
  (for*/flomat n n 
                ([i (in-range n)]
                 [j (in-range n)])
     (vector-ref v (remainder (+ i j) n))))

  
(define (flomat-outer-product A B)
  ; accept standard vectors as input
  (define A1 (if (vector? A) (vector->flcolumn A) A))
  (define B1 (if (vector? B) (vector->flrow    B) B))
  ; compute outer product between first column of A and first row of B
  (define-values (am an) (flomat-dimensions A1))
  (define-values (bm bn) (flomat-dimensions B1))
  (for*/flomat am bn ([a (in-col A1 0)]
                        [b (in-row    B1 0)])
                 (* a b)))


; Since the matrix entries of a column are stored contigious,
; we can use cblas_ixamax with incX=1 to find pivots.
(define (flomat-find-partial-pivot A i j)
  ; Find the index k of the element a_kj with k>=i
  ; that has the largest absolute value.
  ; I.e. a partial pivot in row j.
  (define-param (m n a lda) A)
  (define ptr (ptr-elm a lda i j)) ; address of the (i,j)th element.
  (define idx (cblas_idamax (- m i) ptr 1))
  (+ i idx))


(define (flomat-gauss-elim! A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (define A-original A)
  
  (let loop ([A A] [i 0] [j 0] [without-pivot '()])
    (define-param (m n a lda) A)

    (define (eliminate-row! i pivot l)
      ; eliminate row l using row i which has pivot 
      (define x (unsafe-ref a lda l 0))                   ; first element in row l
      (flomat-add-scaled-row! A l (* -1 (/ x pivot)) i) ; scale and subtract
      (unsafe-set! a lda l 0 0.0))                        ; exact 0.0
    
    (define (eliminate-rows-below! i pivot) (for ([l (in-range (+ i 1) m)]) (eliminate-row! i pivot l)))
    (define (eliminate-rows-above! i pivot) (for ([l (in-range 0 i)])       (eliminate-row! i pivot l)))

    (define (A-without-first-column) (shared-submatrix! A 0 1 m (- n 1)))
    
    (cond
      [(= n 0) (values A-original (reverse without-pivot))]
      ;; None of the rest of the columns can have pivots
      [(= i m) (values A-original (append (reverse without-pivot) (range j (+ j n))))]
      [else    (define p  (case pivoting
                            [(partial) (flomat-find-partial-pivot A i 0)]
                            #;[(first)   (flomat-find-first-pivot   A 0 0)]
                            [else (error 'flomat-gauss-elim! "unknown pivoting type")]))
               (define pivot (flomat-ref A p 0))
               (cond                 
                 [(<= pivot epsilon) ;; no pivot
                  (loop (A-without-first-column) i (+ j 1) (cons j without-pivot))]
                 [else               ;; pivot found
                  (flomat-swap-rows! A i p)
                  (eliminate-rows-below! i pivot)
                  (when jordan?        (eliminate-rows-above! i pivot))
                  (when unitize-pivot? (flomat-scale-row! A i (/ 1. pivot)))
                  (loop (A-without-first-column) (+ i 1) (+ j 1) without-pivot)])])))

(define (flomat-gauss-elim A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (flomat-gauss-elim! (copy-flomat A) jordan? unitize-pivot? pivoting))

(define (matrix-row-echelon! A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (flomat-gauss-elim! A jordan? unitize-pivot? pivoting)
  A)

(define (matrix-row-echelon A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (matrix-row-echelon! (copy-flomat A) jordan? unitize-pivot? pivoting))

;;;
;;; RANDOM NUMBERS
;;;

; A seed is an integer array of dimension 4.
; On entry the elements must be between 0 and 4095,
; and the last element must be odd.
; The routines update the seed automatically.

(define _iseed (_cpointer 'flomat-iseed))

(define (generate-iseed)
  (define i* (cast (malloc 4 _int 'atomic)_pointer _iseed))
  (define r0 (random 4096))
  (define r1 (random 4096))
  (define r2 (random 4096))
  (define r3 (random 4096))
  (let loop ()
    (unless (odd? r3)
      (set! r3 (random 4096))
      (loop)))
  (ptr-set! i* _int 0 r0)
  (ptr-set! i* _int 1 r1)
  (ptr-set! i* _int 2 r2)
  (ptr-set! i* _int 3 r3)
  i*)

(define the-iseed (generate-iseed))

; DLARNV returns a vector of n random real numbers from a uniform or normal distribution.
(define-lapack dlarnv_
  (_fun (idist : (_ptr i _int)) ; i   1: uniform (0,1), 2: uniform (-1,1), 3: normal(0,1)
        (iseed : _iseed)        ; io
        (n    : (_ptr i _int))
        (x    : _flomat)
        -> _void))

; DLARND returns a random real number from a uniform or normal distribution.
; Not available on macOS (using the default lapack)
#;(define-lapack dlarnd_
  (_fun (idist : (_ptr i _int)) ; i   1: uniform (0,1), 2: uniform (-1,1), 3: normal(0,1)
        (iseed : _iseed)        ; io
        -> _double))


; uniform numbers: (0,1)
(define rand
  ; 1: uniform (0,1)
  (case-lambda
    [()
     (define a  (alloc-flomat 1 1))
     (dlarnv_ 1 the-iseed 1 a)
     (ptr-ref a _double 0)]
    [(n)
     (check-positive-integer 'rand n)
     (define a  (alloc-flomat n n))
     (dlarnv_ 1 the-iseed (* n n) a)
     (flomat n n a n)]
    [(m n)
     (check-positive-integer 'rand m)
     (check-positive-integer 'rand n)
     (define a  (alloc-flomat m n))
     (dlarnv_ 1 the-iseed (* m n) a)
     (flomat m n a m)]))


; standard normal distribution (0,1)
(define randn
  ; 1: normal distribution (0,1)
  (case-lambda
    [()
     (define a  (alloc-flomat 1 1))
     (dlarnv_ 3 the-iseed 1 a)
     (ptr-ref a _double 0)]
    [(n)
     (check-positive-integer 'randn n)
     (define a  (alloc-flomat n n))
     (dlarnv_ 3 the-iseed (* n n) a)
     (flomat n n a n)]
    [(m n)
     (check-positive-integer 'randn m)
     (check-positive-integer 'randn n)
     (define a  (alloc-flomat m n))
     (dlarnv_ 3 the-iseed (* m n) a)
     (flomat m n a m)]))


;;;
;;; HIGH LEVEL
;;;

(define (mat? A)  (flomat? A))
(define (row? A)  (and (flomat? A) (flomat-row-vector? A)))
(define (col? A)  (and (flomat? A) (flomat-column-vector? A)))

(define (row A i)   (flomat-row    A i))
(define (col A i)   (flomat-column A i))
(define (ref A i j) (flomat-ref    A i j))

(define (shape A)   (list (flomat-m A) (flomat-n A)))
(define (size  A)   (flomat-size A))
(define (nrows A)   (flomat-m A))
(define (ncols A)   (flomat-n A))

(define augment        flomat-augment)
(define stack          flomat-stack)
(define repeat         flomat-repeat)
(define block-diagonal flomat-block-diagonal)

(define (f64vector->flomat v [transpose? #f])
  (define m (f64vector-length v))
  (define n 1)
  (define lda 1)
  (define a (cast (f64vector->cpointer v) _pointer _flomat))
  (if transpose?
      (flomat n m a lda)
      (flomat m n a lda)))

(define (column . xs)
  (matrix (map list xs)))

(define (matrix x)
  (cond
    [(vector? x)    (cond
                      ; vector of vector
                      [(vector? (vector-ref x 0)) (vectors->flomat x)]
                      ; a single vector represents a column vector
                      [else                       (vector->flcolumn x)])]
    [(list? x)      (cond
                      ; list of lists
                      [(list? (first x))          (list->flomat x)]
                      [else                       (apply flcolumn x)])]
    [(f64vector? x) (copy-flomat (matrix! x))]
    [else (error)]))

(define (matrix! x)
  (cond
    [(f64vector? x) (f64vector->flomat x)]
    [else
     (error 'matrix! "expected an f64vector as input")]))

(provide f64vector)

(define (mset! A i j x) (flomat-set! A i j x))

(define (zeros m [n m])    (flomat-zeros m n))
(define (ones  m [n m])    (flomat-ones  m n))
(define (make m n [x 0.0]) (make-flomat m n [x 0.0]))
(define (constant! A x)
  (define x* (real->double-flonum x))
  (define-param (m n a lda) A)
  (for* ([i (in-range m)] [j (in-range n)])
    (unsafe-set! a lda i j x*)))
(define (zeros! A) (constant! A 0.0))
(define (ones!  A) (constant! A 1.0))


(define (do-arange start stop step transpose?)
  (define len (inexact->exact (ceiling (/ (- stop start) step))))
  (define v (make-f64vector len))
  (for ([x (in-range (* 1.0 start) (* 1.0 stop) (* 1. step))]
        [k (in-naturals)])
    (f64vector-set! v k x))
  (f64vector->flomat v transpose?))

(define (colarange start [stop #f] [step 1.0])
  (cond [stop (do-arange start stop step #f)]
        [else (colarange 0. start step)]))

(define (arange start [stop #f] [step 1.0])
  (cond [stop (do-arange start stop step #t)]
        [else (arange 0. start step)]))


(define (reshape A m n)
  (reshape! (copy-flomat A) m n))
  
(define (reshape! A m n)
  (define-param (M N a lda) A)
  (unless (<= (* m n) (* M N))
    (error 'reshape!
           "the size of the new shape is larger than the original matrix"))
  (cond
    [(= lda 1) ; A is a row vector => no gaps
     (flomat m n a m)]
    [(= lda M) ; no gap between each column
     (flomat m n a m)]
    ; [if ... M N m n ... has a nice relationship then ...]
    [else
     (error 'reshape! "not able to reshape")]))

(define (transpose A) (flomat-transpose A))

(define (linspace start stop [num 50] [endpoint #t])
  (define step (/ (- stop start) (* 1.0 (- num 1))))
  (define len (if endpoint num (- num 1)))
  (define v (make-f64vector len))
  (for/list ([k (in-range len)])
    (define x (+ start (* k step)))
    (f64vector-set! v k x))
  (f64vector->flomat v))

(define (sub! A i j r s) (shared-submatrix! A i j r s))
  ; return rxs matrix with upper left corner (i,j)
  ; entries are shared with A
(define (sub A i j m n) (flsubmatrix A m n i j))
  ; return a the mxn submatrix of with upper  left corner in (i,j)

(define (col! A j)
  (check-legal-column 'col j A)
  (define m (flomat-m A)) ; nrows
  (sub! A 0 j m 1))

(define (row! A i)
  (check-legal-row 'row i A)
  (define n (flomat-n A)) ; ncols
  (sub! A i 0 1 n))


;;;
;;; Pointwise Operations
;;;

(define-syntax (define-pointwise-unary stx)
  (syntax-parse stx
    [(_define-pointwise f:id)
     (with-syntax ([.f! (format-id #'f ".~a!" (syntax-e #'f))]
                   [.f  (format-id #'f ".~a"  (syntax-e #'f))])
       (syntax/loc stx
         (begin
           (define (.f! A [C #f])         
             (when (flomat? C) (check-same-dimensions A C '.f))
             (unless C (set! C A))
             (define-param (m n a lda) A)
             (define-param (_ __ c ldc) C)
             (for* ([i (in-range m)]
                    [j (in-range n)])
               (define aij (unsafe-ref a lda i j))
               (define x (f aij))
               (unsafe-set! c ldc i j x))
             C)
           (define (.f A)
             (.f! (copy-flomat A))))))]))

(define-syntax (define-pointwise-unaries stx)
  (syntax-parse stx
    [(_ f:id ...)
     (syntax/loc stx
       (begin
         (define-pointwise-unary f) ...))]))

(define-pointwise-unaries sin cos tan sqr sqrt log exp)

(define-syntax (define-pointwise-binary stx)
  (syntax-parse stx
    [(_define-pointwise f:id)
     (with-syntax ([.f! (format-id #'f ".~a!" (syntax-e #'f))]
                   [.f  (format-id #'f ".~a"  (syntax-e #'f))])
       (syntax/loc stx
         (begin
           (define (.f! A B [C #f])
             (cond
               [(and (flomat? A) (flomat? B))
                (check-same-dimensions A B '.f!)
                (when (flomat? C) (check-same-dimensions A C '.f!))
                (unless C (set! C A))
                (define-param (m n  a lda) A)
                (define-param (M N  b ldb) B)             
                (define-param (_ __ c ldc) C)             
                (for* ([i (in-range m)]
                       [j (in-range n)])
                  (define aij (unsafe-ref a lda i j))
                  (define bij (unsafe-ref b ldb i j))
                  (define x (f aij bij))
                  (unsafe-set! c ldc i j x))
                C]
               [(number? A) ; now C needs to be an flomat
                (unless (flomat? C)
                  (error (error '.! "if A is a constant, C must be a flomat")))                
                (cond [(flomat? B) (define-param (m n  b ldb) B)
                                     (.f! (make-flomat m n A) B C)]
                      [else (error '.! "wrong input types")])]
               [(number? B)
                (cond [(flomat? A)
                       (when (flomat? C) (check-same-dimensions A C '.f!))
                       (unless C (set! C A))
                       (define-param (m n  a lda) A)
                       (define-param (_ __ c ldc) C)             
                       (for* ([i (in-range m)]
                              [j (in-range n)])
                         (define aij (unsafe-ref a lda i j))
                         (define x (f aij B))
                         (unsafe-set! c ldc i j x))
                       C]
                      [else (error '.! "wrong input types")])]))
           (define (.f A B)
             (cond
               [(flomat? A) (.f! A B (copy-flomat A))]
               [(flomat? B) (.f! A B (copy-flomat B))]
               [else          (flomat 1 1 (f A B))])))))]))  

(define-syntax (define-pointwise-binaries stx)
  (syntax-parse stx
    [(_ f:id ...)
     (syntax/loc stx
       (begin
         (define-pointwise-binary f) ...))]))

(define-pointwise-binaries + * expt)

(define-syntax (define-pointwise-unary/binary stx)
  (syntax-parse stx
    [(_define-pointwise f:id)
     (with-syntax ([.f! (format-id #'f ".~a!" (syntax-e #'f))]
                   [.f  (format-id #'f ".~a"  (syntax-e #'f))])
       (syntax/loc stx
         (begin
           (define (.f! A [B #f] [C #f])
             (cond
               [B ; binary
                (cond
                  [(and (flomat? A) (flomat? B))
                   (check-same-dimensions A B '.f!)
                   (when (flomat? C) (check-same-dimensions A C '.f!))
                   (unless C (set! C A))
                   (define-param (m n  a lda) A)
                   (define-param (M N  b ldb) B)             
                   (define-param (_ __ c ldc) C)             
                   (for* ([i (in-range m)]
                          [j (in-range n)])
                     (define aij (unsafe-ref a lda i j))
                     (define bij (unsafe-ref b ldb i j))
                     (define x (f aij bij))
                     (unsafe-set! c ldc i j x))
                   C]
                  [(number? A) ; now C needs to be an flomat
                   (unless (flomat? C)
                     (error (error '.! "if A is a constant, C must be a flomat")))                
                   (cond [(flomat? B) (define-param (m n  b ldb) B)
                                        (.f! (make-flomat m n A) B C)]
                         [else (error '.! "wrong input types")])]
                  [(number? B)
                   (cond [(flomat? A)
                          (when (flomat? C) (check-same-dimensions A C '.f!))
                          (unless C (set! C A))
                          (define-param (m n  a lda) A)
                          (define-param (_ __ c ldc) C)             
                          (for* ([i (in-range m)]
                                 [j (in-range n)])
                            (define aij (unsafe-ref a lda i j))
                            (define x (f aij B))
                            (unsafe-set! c ldc i j x))
                          C]
                         [else (error '.! "wrong input types")])])]
               ; B is #f
               [else  ; unary with result to A or C
                (unless C (set! C A))
                (check-same-dimensions A C '.f!)
                (define-param (m n  a lda) A)
                (define-param (_ __ c ldc) C)
                (for* ([i (in-range m)]
                       [j (in-range n)])
                  (define aij (unsafe-ref a lda i j))
                  (define x (f aij))
                  (unsafe-set! c ldc i j x))
                C]))           
           (define (.f A [B #f])
             (if B
                 (.f! A B (copy-flomat A))
                 (.f! (copy-flomat A) B))))))]))
    
(define-syntax (define-pointwise-unary/binaries stx)
  (syntax-parse stx
    [(_ f:id ...)
     (syntax/loc stx
       (begin
         (define-pointwise-unary/binary f) ...))]))

(define-pointwise-unary/binaries - /)

(define (plus! A . BS)
  (check-flomat 'plus! A)
  (check-all-matrices-same-size 'plus! (cons A BS))
  (let loop ([BS BS])
    (cond
      [(empty? BS)            A]
      [(flomat? (first BS)) (flomat+! (first BS) A) (loop (rest BS))]
      [(number?   (first BS)) (.+! A (first BS))       (loop (rest BS))]))
  A)

(define (plus A . BS)
  (check-all-matrices-same-size 'plus (cons A BS))
  (let loop ([A A] [BS BS])
    (cond
      [(empty? BS)            A]
      [(flomat? (first BS)) (if (number? A)
                                  (loop (.+ A        (first BS)) (rest BS))
                                  (loop (flomat+ A (first BS)) (rest BS)))]
      [(number?   (first BS))     (loop (.+        A (first BS)) (rest BS))])))

(define (minus! A . BS)
  (check-flomat 'minus! A)
  (check-all-matrices-same-size 'minus! (cons A BS))
  (cond
    [(empty? BS) (.-! A)]
    [else        (let loop ([BS BS])
                   (cond
                     [(empty? BS)            A]
                     [(flomat? (first BS)) (flomat-! (first BS) A) (loop (rest BS))]
                     [(number?   (first BS)) (.-! A (first BS))        (loop (rest BS))]))
                 A]))

(define (minus A . BS)
  (check-all-matrices-same-size 'minus (cons A BS))
  (cond
    [(empty? BS) (.- A)]
    [else        (let loop ([A A] [BS BS])
                   (cond
                     [(empty? BS)            A]
                     [(flomat? (first BS)) (if (number? A)
                                                 (loop (.- A        (first BS)) (rest BS))
                                                 (loop (flomat- A (first BS)) (rest BS)))]
                     [(number?   (first BS))     (loop (.-        A (first BS)) (rest BS))]))]))

(define (times! A . BS)
  (check-flomat 'times! A)
  ; todo: check sizes
  (let loop ([BS BS])
    (cond
      [(empty? BS)            A]
      [(flomat? (first BS)) (flomat*! A (first BS) A) (loop (rest BS))]
      [(number?   (first BS)) (.*! A (first BS))          (loop (rest BS))]))
  A)

(define (times A . BS)
  ; todo: check sizes
  (let loop ([A A] [BS BS])
    (cond
      [(empty? BS)            A]
      [(flomat? (first BS)) (if (number? A)
                                  (loop (.* A        (first BS)) (rest BS))
                                  (loop (flomat* A (first BS)) (rest BS)))]
      [(number?   (first BS)) (if (number? A)
                                  (loop (*         A (first BS)) (rest BS))
                                  (loop (.*        A (first BS)) (rest BS)))])))

(define × times)

(define (dot A B)
  (flcolumn-dot A B))

(define (outer A B)
  ; compute outer product between first column of A and first row of B
  (flomat-outer-product A B))

(define (power A n)
  ; n natural
  (flomat-expt A n))

(define (kron A B)
  (define-param (m n   a lda) A)
  (define-param (mb nb b ldb) B)
  (define C (make-flomat (* m mb) (* n nb)))
  (define-param (mc nc c ldc) C)
  
  (for* ([i (in-range m)]
         [j (in-range n)])
    (define k (* i mb))
    (define l (* j nb))
    (define c_kl (ptr-elm c ldc k l))
    (unsafe-matrix-copy! mb nb b lda c_kl ldc)
    (define aij (unsafe-ref a lda i j))
    (flomat-scale! aij (sub! C k l mb nb)))
  C)


(define (diag X [m? #f] [n? #f] [reciproc? #f])
  (set! X (result-flcolumn X))
  ; todo: use optional arguments to denote diagonal like matlab?
  (define-param (m n) X)
  (flomat-diagonal-from-singular-values (or m? m) (or n? m) X reciproc?))


(define (cholesky A [triangle 'lower])
  ; triangle is 'lower or 'upper
  (define lower? (member triangle '(lower low l)))
  (flomat-cholesky A (not lower?)))

(define (svd A)                  (flomat-svd A))
(define (qr A)                   (flomat-qr A))

(define (eig A)
  (define who 'eig)
  (check-square who A)
  (define-values (A0 WR WI VL VR info)
    (flomat-eigenvalues-and-vectors! A #:right #t #:overwrite #f))
  (values (real+imaginary->vector WR WI) VR))

(define (eigvals A)
  (define who 'eigvals)
  (check-square who A)
  (define-values (A0 WR WI VL VR info)
    (flomat-eigenvalues-and-vectors! A #:right #f #:overwrite #f))
  (real+imaginary->vector WR WI))

(define (norm A [type 'frob])
  (flomat-norm A type))

(define (det A)
  (flomat-determinant A))

(define (trace A)
  (flomat-trace A))

(define (rank A)
  ; rank = dimension of column space = dimension of row space  
  ;      = number of non-zero singular values
  (flomat-rank A))

; todo: condition number

(define (mldivide A B)
  ; todo: also handle non-square A
  (flomat-solve-many A B))

(define (mrdivide B A)
  ; B/A = (A'\B')'
  (transpose (mldivide (transpose A) (transpose B))))

(define (inv A)
  (flomat-inverse A))

(define (pinv A)
  (flomat-pseudo-inverse A))

(define (eye m [n m] [k 0])
  (flomat-eye m n k))

(define rowsum  flomat-row-sum)
(define rowsums flomat-row-sums)
(define colsum  flomat-column-sum)
(define colsums flomat-column-sums)
(define (summ A) (rowsum (colsums A) 0))

;;;
;;; TEST
;;;


(module+ test
  (require rackunit)
  
  (define (flcheck-equal? a b)
    (< (abs (- b a)) 0.00001))
  

  (with-check-info
   (['test-case "flomat/dim"])
   (check-equal? (flomat->vector (flomat/dim 2 2 1 2 3 4))
                 #(1. 2. 3. 4.))
   (check-equal? (flomat->vector (vector->flomat 2 2 #(1 2 3 4)))
                 #(1. 2. 3. 4.))
   (check-equal? (flomat->vector (flomat/dim 2 2 1 2 3 4))
                 #(1. 2. 3. 4.))
   (check-equal? (flomat->vectors (flomat/dim 2 2 1 2 3 4))
                 #(#[1. 2.] #[3. 4.]))

   (let ()
     (define A  (flomat/dim 2 2 1 2 3 4))
     (define B  (flomat/dim 2 2 5 6 7 8))
     (define AB (flomat/dim 2 2 19 22 43 50))
     (check-equal? (flomat->vectors (flomat* A B))
                   (flomat->vectors AB)))

   (let ()
     (define C  (flomat/dim 2 2 1 2 3 4))
     (define D  (flomat/dim 2 3 5 6 7 8 9 10))
     (define CD (flomat/dim 2 3 21 24 27 47 54 61))
     (check-equal? (flomat->vectors (flomat* C D))
                   (flomat->vectors CD)))

   (check-equal? (flomat->vectors 
                  (flomat* (flomat/dim 2 3 0 0 1 0 0 0)
                             (flomat/dim 3 2 1 2 3 4 5 6)))
                 (flomat->vectors 
                  (flomat/dim 2 2 5 6 0 0))))


  (with-check-info
   (['test-group "matrix-constructors.rkt"])
   (with-check-info
    (['test-case 'flomat-identity])
          
    (check-equal? (flomat->lists (flomat-identity 1)) '[[1.]])
    (check-equal? (flomat->lists (flomat-identity 2)) '[[1. 0.] [0. 1.]])
    (check-equal? (flomat->lists (flomat-identity 3)) '[[1. 0. 0.] [0. 1. 0.] [0. 0. 1.]]) 
    (check-equal? (flomat->lists (flomat-identity 1)) '[[1.]])
    (check-equal? (flomat->lists (flomat-identity 2)) '[[1. 0.] [0. 1.]])
    (check-equal? (flomat->lists (flomat-identity 3)) '[[1. 0. 0.] [0. 1. 0.] [0. 0. 1.]]))
   (with-check-info
    (['test-case 'const-matrix])          
    (check-equal? (flomat->lists (make-flomat 2 3 0.)) '((0. 0. 0.) (0. 0. 0.))))
   (with-check-info
    (['test-case 'matrix->list])
    (check-equal? (flomat->lists (lists->flomat '((1. 2.) (3. 4.)))) '((1. 2.) (3. 4.))))
   (with-check-info
    (['test-case 'matrix->vectors])
    (check-equal? (flomat->vectors (vectors->flomat '#(#(1. 2.) #(3. 4.)))) '#(#(1. 2.) #(3. 4.))))
   (with-check-info
    (['test-case 'matrix-row])
    (check-equal? (flomat-row (flomat-identity 3) 0) (list->flomat '[[1 0 0]]))
    (check-equal? (flomat-row (flomat-identity 3) 1) (list->flomat '[[0 1 0]]))
    (check-equal? (flomat-row (flomat-identity 3) 2) (list->flomat '[[0 0 1]])))
   (with-check-info
    (['test-case 'matrix-col])
    (check-equal? (flomat-column (flomat-identity 3) 0) (list->flomat '[[1] [0] [0]]))
    (check-equal? (flomat-column (flomat-identity 3) 1) (list->flomat '[[0] [1] [0]]))
    (check-equal? (flomat-column (flomat-identity 3) 2) (list->flomat '[[0] [0] [1]])))
   (with-check-info
    (['test-case 'flsubmatrix])
    (check-equal? (flsubmatrix (flomat-identity 3) 1 2 0 0)
                  (list->flomat '[[1 0]]))
    (check-equal? (flsubmatrix (flomat-identity 3) 2 3 0 0)
                  (list->flomat '[[1 0 0] [0 1 0]]))))

  (with-check-info
    (['test-group "product"])
    (check-equal?
     (flomat*vector (list->flomat '[[1 2 3] [4 5 6]])
                      (list->flomat '[[10] [11] [12]]))
     (list->flomat '[[68] [167]]))
    (check-equal?
     (flomat*vector (list->flomat '[[1 4] [2 5] [3 6]])
                      (list->flomat '[[10] [11] [12]])
                      #f 1. 1. #t)
     (list->flomat '[[68] [167]])))

  (with-check-info
   (['test-group "flomat-pointwise.rkt"])
   (let ()
     (define A   (list->flomat '[[1 2] [3 4]]))
     (define ~A  (list->flomat '[[-1 -2] [-3 -4]]))
     (define B   (list->flomat '[[5 6] [7 8]]))
     (define A+B (list->flomat '[[6 8] [10 12]]))
     (define A-B (list->flomat '[[-4 -4] [-4 -4]]))         
     (with-check-info
      (['test-case 'flomat+])
      (check-equal? (flomat+ A B) A+B))
     (with-check-info
      (['test-case 'flomat-])
      (check-equal? (flomat- A B) A-B)
      (check-equal? (flomat- A)   ~A))))
  
  (with-check-info
   (['test-group "flomat-expt.rkt"])
   (define A (list->flomat '[[1 2] [3 4]]))
   (with-check-info
    (['test-case 'flomat-expt])
    (check-equal? (flomat-expt A 0) (flomat-identity 2))
    (check-equal? (flomat-expt A 1) A)
    (check-equal? (flomat-expt A 2) (list->flomat '[[7 10] [15 22]]))
    (check-equal? (flomat-expt A 3) (list->flomat '[[37 54] [81 118]]))
    (check-equal? (flomat-expt A 8) (list->flomat '[[165751 241570] [362355 528106]]))))

  (with-check-info
   (['test-group "flomat-operations.rkt"])
   (with-check-info
    (['test-case 'vandermonde-flomat])
    (check-equal? (flomat-vandermonde '(1 2 3) 5)
                  (list->flomat '[[1 1 1 1 1] [1 2 4 8 16] [1 3 9 27 81]])))
   (with-check-info
    (['test-case 'in-column])
    (check-equal? (for/list ([x (in-col (flomat/dim 2 2  1 2 3 4) 0)]) x)
                  '(1. 3.))
    (check-equal? (for/list ([x (in-col (flomat/dim 2 2  1 2 3 4) 1)]) x)
                  '(2. 4.))
    (check-equal? (for/list ([x (in-col (flcolumn 5 2 3))]) x)
                  '(5. 2. 3.)))
   (with-check-info
    (['test-case 'in-row])
    (check-equal? (for/list ([x (in-row (flomat/dim 2 2  1 2 3 4) 0)]) x)
                  '(1. 2.))
    (check-equal? (for/list ([x (in-row (flomat/dim 2 2  1 2 3 4) 1)]) x)
                  '(3. 4.)))
   (with-check-info
    (['test-case 'for/flomat:])
    (check-equal? (for/flomat 2 4 ([i (in-naturals)]) i)
                  (flomat/dim 2 4 
                                0 1 2 3
                                4 5 6 7))
    (check-equal? (for/flomat 2 4 #:column ([i (in-naturals)]) i)
                  (flomat/dim 2 4    
                                0 2 4 6
                                1 3 5 7))
    (check-equal? (for/flomat 3 3 ([i (in-range 10 100)]) i)
                  (flomat/dim 3 3 10 11 12 13 14 15 16 17 18)))
   (with-check-info
    (['test-case 'for*/flomat:])
    (check-equal? (for*/flomat 3 3 ([i (in-range 3)] [j (in-range 3)]) (+ (* i 10) j))
                  (flomat/dim 3 3 0 1 2 10 11 12 20 21 22)))    
   (with-check-info
    (['test-case 'flomat-block-diagonal])
    (check-equal? (flomat-block-diagonal (flomat/dim 2 2 1 2 3 4) (flomat/dim 1 3 5 6 7))
                  (list->flomat '[[1 2 0 0 0] [3 4 0 0 0] [0 0 5 6 7]])))
   (with-check-info
    (['test-case 'flomat-augment])
    (check-equal? (flomat-augment (flcolumn 1 2 3) (flcolumn 4 5 6) (flcolumn 7 8 9))
                  (flomat/dim 3 3  1 4 7  2 5 8  3 6 9)))
   (with-check-info
    (['test-case 'flomat-stack])
    (check-equal? (flomat-stack (flcolumn 1 2 3) (flcolumn 4 5 6) (flcolumn 7 8 9))
                  (flcolumn 1 2 3 4 5 6 7 8 9)))
   (with-check-info
    (['test-case 'column-dimension])
    (= (flcolumn-size #(1 2 3)) 3)
    (= (flcolumn-size (vector->flomat 1 2 #(1 2))) 1))
   (let ([flomat: vector->flomat])
     (with-check-info
      (['test-case 'column-dot])
      (= (flcolumn-dot (flcolumn 1 2)   (flcolumn 1 2)) 5)
      (= (flcolumn-dot (flcolumn 1 2)   (flcolumn 3 4)) 11)
      (= (flcolumn-dot (flcolumn 3 4)   (flcolumn 3 4)) 25)
      (= (flcolumn-dot (flcolumn 1 2 3) (flcolumn 4 5 6))
         (+ (* 1 4) (* 2 5) (* 3 6)))))
   (with-check-info
    (['test-case 'flomat-trace])
    (check-equal? (flomat-trace (vector->flomat 2 2 #(1 2 3 4))) 5.))
   (let ([flomat: vector->flomat])
     (with-check-info
      (['test-case 'column-norm])
      (= (flcolumn-norm (flcolumn 2 4)) (sqrt 20))))
   (with-check-info
    (['test-case 'column-projection])
    (check-equal? (flcolumn-projection #(1 2 3) #(4 5 6)) (flcolumn 128/77 160/77 192/77))
    (check-equal? (flcolumn-projection (flcolumn 1 2 3) (flcolumn 2 4 3))
                  (flomat-scale 19/29 (flcolumn 2 4 3))))
   (with-check-info
    (['test-case 'projection-on-orthogonal-basis])
    (check-equal? (flprojection-on-orthogonal-basis #(3 -2 2) (list #(-1 0 2) #( 2 5 1)))
                  (flcolumn -1/3 -1/3 1/3))
    (check-equal? (flprojection-on-orthogonal-basis 
                   (flcolumn 3 -2 2) (list #(-1 0 2) (flcolumn 2 5 1)))
                  (flcolumn -1/3 -1/3 1/3)))
   (with-check-info
    (['test-case 'projection-on-orthonormal-basis])
    (check-equal? (flprojection-on-orthonormal-basis 
                   #(1 2 3 4) 
                   (list (flomat-scale 1/2 (flcolumn  1  1  1 1))
                         (flomat-scale 1/2 (flcolumn -1  1 -1 1))
                         (flomat-scale 1/2 (flcolumn  1 -1 -1 1))))
                  (flcolumn 2 3 2 3)))
   (with-check-info
    (['test-case 'flgram-schmidt-orthogonal])
    (check-equal? (flgram-schmidt-orthogonal (list #(3 1) #(2 2)))
                  (list (flcolumn 3 1) (flcolumn -2/5 6/5))))
   (with-check-info
    (['test-case 'flvector-normalize])
    (check-equal? (flcolumn-normalize #(3 4)) 
                  (flcolumn 3/5 4/5)))
   (with-check-info
    (['test-case 'flgram-schmidt-orthonormal])
    (check-equal? (flgram-schmidt-orthonormal '(#(3 1) #(2 2)))
                  (list (flcolumn-normalize #(3 1))
                        (flcolumn-normalize #(-2/5 6/5)))))
    
   (with-check-info
    (['test-case 'projection-on-subspace])
    (check-equal? (flprojection-on-subspace #(1 2 3) '(#(2 4 3)))
                  (flomat-scale 19/29 (flcolumn 2 4 3))))
   (with-check-info
    (['test-case 'unit-vector])
    (check-equal? (flcolumn-unit 4 1) (flcolumn 0 1 0 0)))
    (with-check-info (['test-case 'flomat-qr])
      (let*-values ([(A) (flomat/dim 3 2  1 1 0 1 1 1)]
                    [(Q R) (flomat-qr A)])
        (check-true
         (flomat= (flomat* Q R)
                    A
                    epsilon))))
   (with-check-info
    (['test-case 'flomat-solve])
    (let* ([M (list->flomat '[[1 5] [2 3]])] 
           [b (list->flomat '[[5] [5]])])
      (check-equal? (flomat* M (flomat-solve M b)) b)))
   (with-check-info
    (['test-case 'flomat-inverse])
    (check-equal? (let ([M (list->flomat '[[1 2] [3 4]])]) (flomat* M (flomat-inverse M)))
                  (flomat-identity 2))
    (check-equal? (let ([M (list->flomat '[[1 2] [3 4]])]) (flomat* (flomat-inverse M) M))
                  (flomat-identity 2)))
   (with-check-info
    (['test-case 'flomat-determinant])
    (check-equal? (flomat-determinant (list->flomat '[[3]])) 3.)
    (check-equal? (flomat-determinant (list->flomat '[[1 2] [3 4]])) (- (* 1. 4.) (* 2. 3.)))
    (flcheck-equal? (flomat-determinant (list->flomat '[[1 2 3] [4  5 6] [7 8 9]])) 0.)
    (flcheck-equal? (flomat-determinant (list->flomat '[[1 2 3] [4 -5 6] [7 8 9]])) 120.)
    (flcheck-equal? (flomat-determinant 
                     (list->flomat '[[1 2 3 4] [-5 6 7 8] [9 10 -11 12] [13 14 15 16]])) 5280.))
   (with-check-info
    (['test-case 'flomat-scale])
    (check-equal? (flomat-scale 2 (list->flomat '[[1 2] [3 4]]))
                  (list->flomat '[[2 4] [6 8]])))
   (with-check-info
    (['test-case 'flomat-transpose])
    (check-equal? (flomat-transpose (list->flomat '[[1 2] [3 4]]))
                  (list->flomat '[[1 3] [2 4]])))
   ; TODO: Just use U from LU factorization
   #;(let ()
       (: gauss-eliminate : (flomat Number) Boolean Boolean -> (flomat Number))
       (define (gauss-eliminate M u? p?)
         (let-values ([(M wp) (flomat-gauss-eliminate M u? p?)])
           M))
       (with-check-info
        (['test-case 'flomat-gauss-eliminate])
        (check-equal? (let ([M (list->flomat '[[1 2] [3 4]])])
                        (gauss-eliminate M #f #f))
                      (list->flomat '[[1 2] [0 -2]]))
        (check-equal? (let ([M (list->flomatixix  '[[2 4] [3 4]])])
                        (gauss-eliminate M #t #f))
                      (list->flomatixixix '[[1 2] [0 1]]))
        (check-equal? (let ([M (list->flomatix  '[[2. 4.] [3. 4.]])])
                        (gauss-eliminate M #t #t))
                      (list->flomatix '[[1. 1.3333333333333333] [0. 1.]]))
        (check-equal? (let ([M (list->flomat  '[[1 4] [2 4]])])
                        (gauss-eliminate M #t #t))
                      (list->flomat '[[1 2] [0 1]]))
        (check-equal? (let ([M (list->flomat  '[[1 2] [2 4]])])
                        (gauss-eliminate M #f #t))
                      (list->flomat '[[2 4] [0 0]]))))
   (with-check-info
    (['test-case 'flomat-scale-row])
    (check-equal? (flomat-scale-row (flomat-identity 3) 0 2)
                  (lists->flomat '[[2 0 0] [0 1 0] [0 0 1]])))
   (with-check-info
    (['test-case 'flomat-swap-rows])
    (check-equal? (flomat-swap-rows (lists->flomat '[[1 2 3] [4 5 6] [7 8 9]]) 0 1)
                  (lists->flomat '[[4 5 6] [1 2 3] [7 8 9]])))
   (with-check-info
    (['test-case 'flomat-add-scaled-row])
    (check-equal? (flomat-add-scaled-row (lists->flomat '[[1 2 3] [4 5 6] [7 8 9]]) 0 2 1)
                  (lists->flomat '[[9 12 15] [4 5 6] [7 8 9]])))
   (let ()
     (define M (lists->flomat '[[1  1  0  3]
                                  [2  1 -1  1]
                                  [3 -1 -1  2]
                                  [-1  2  3 -1]]))
     (define-values (P L U) (flomat-plu M))
     (with-check-info
      (['test-case 'flomat-plu])
      (check-equal? (flomat* P (flomat* L U)) M)))
   (with-check-info
    (['test-case 'flomat-rank])
    (check-equal? (flomat-rank (list->flomat '[[0 0] [0 0]])) 0)
    (check-equal? (flomat-rank (list->flomat '[[1 0] [0 0]])) 1)
    (check-equal? (flomat-rank (list->flomat '[[1 0] [0 3]])) 2)
    (check-equal? (flomat-rank (list->flomat '[[1 2] [2 4]])) 1)
    (check-equal? (flomat-rank (list->flomat '[[1 2] [3 4]])) 2))
   (with-check-info
    (['test-case 'flomat-nullity])
    (check-equal? (flomat-nullity (list->flomat '[[0 0] [0 0]])) 2)
    (check-equal? (flomat-nullity (list->flomat '[[1 0] [0 0]])) 1)
    (check-equal? (flomat-nullity (list->flomat '[[1 0] [0 3]])) 0)
    (check-equal? (flomat-nullity (list->flomat '[[1 2] [2 4]])) 1)
    (check-equal? (flomat-nullity (list->flomat '[[1 2] [3 4]])) 0))
   ; Not implemented yet...
   #;(let ()
       (define-values (c1 n1) 
         (flomat-column+null-space (list->flomat '[[0 0] [0 0]])))
       (define-values (c2 n2) 
         (flomat-column+null-space (list->flomat '[[1 2] [2 4]])))
       (define-values (c3 n3) 
         (flomat-column+null-space (list->flomat '[[1 2] [2 5]])))
       (with-check-info
        (['test-case 'flomat-column+null-space])
        (check-equal? c1 '())
        (check-equal? n1 (list (list->flomat '[[0] [0]])
                               (list->flomat '[[0] [0]])))
        (check-equal? c2 (list (list->flomat '[[1] [2]])))
        ;(check-equal? n2 '([0 0]))
        (check-equal? c3 (list (list->flomat '[[1] [2]])
                               (list->flomat '[[2] [5]])))
        (check-equal? n3 '()))))
  
  (with-check-info
   (['test-group "matrix-multiply.rkt"])
   (with-check-info
    (['test-case 'flomat*])
    (let ()
      (define-values (A B AB) (values '[[1 2] [3 4]] '[[5 6] [7 8]] '[[19 22] [43 50]]))
      (check-equal? (flomat* (list->flomat A) (list->flomat B)) (list->flomat AB)))
    (let () 
      (define-values (A B AB) (values '[[1 2] [3 4]] '[[5 6 7] [8 9 10]] '[[21 24 27] [47 54 61]]))
      (check-equal? (flomat* (list->flomat A) (list->flomat B)) (list->flomat AB))))))

(define (build-flomat m n f)
  (for*/flomat m n 
                 ([i (in-range m)]
                  [j (in-range n)])
                 (f i j)))

