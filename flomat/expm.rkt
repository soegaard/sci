#lang racket/base
(provide expm)

;;;
;;; Matrix Exponential
;;;

; This module implements `expm` the so-called matrix exponential.
; The standard exponential function `exp` can be written as
; as an infinite sum (the Taylor series of exp).

;                x    x^2   x^3
;   exp(x) = 1 + -- + --- + --- + ...
;                1!    2!    3!

; The matrix exponential expm(A) of a square matrix A is
; defined as the limit of:

;                A    A^2   A^3
;   exp(A) = 1 + -- + --- + --- + ...
;                1!    2!    3!

; if the limit exist (otherwise expm(A) is undefined).

; In principle one could use this definition to calculate `expm`,
; but the convergence is slow. Instead we will use
; an approximation of `exp(x)` using a rational expression:
; a quotient of two polynomials.

;             p(x)
;  exp(x) ~ --------  where p and q are polynomials
;             q(x)

; Henri Padé studied how to find the best rational expression
; approximating a given function. Such an approximation is
; called an "Padé approximant" today.

; The precision of the approximation depends on the number
; of terms used in the numerator and denominator. More terms
; give higher precision.

; To keep things simple, we will look at the situation where
; we use the same number of terms in both numerator and denominator.
; This is called a diagonal Padé approximation (in tables they appear
; on the diagonal.

; The first few diagonal Padé approximations are:
;             1
;   exp(x) ~ ---
;             1

;             1 + 1/2 x
;   exp(x) ~ -----------
;             1 - 1/2 x

;             1 + 1/2 x + 1/12 x^2
;   exp(x) ~ -----------------------
;             1 - 1/2 x + 1/12 x^2

;             1 + 1/2 x + 1/10 x^2 + 1/120 x^3
;   exp(x) ~ ---------------------------------
;             1 - 1/2 x + 1/10 x^2 - 1/120 x^3

; The explicit formula for the i'th coefficient of the p'th diagonal approximant is:

;          (2p-i)!    p!
;  N_p  = ----------------
;          (2p)! i! (p-i)!

;          (2p-i)!    p!
;  D_p  = ---------------- *(-1)^i
;          (2p)! i! (p-i)!

; Given these formulas we can compute the coefficients we need,
; when we can experiment with different number of terms.

; This is fine - but what about our matrices?
; As a simple example, let's say we have chosen to use the approximation:

;             1 + 1/2 x
;   exp(x) ~ -----------
;             1 - 1/2 x

;  We can rewrite this as:
;     (1 - 1/2 x) exp(x) =  1 + 1/2 x

; Plugging in a matrix A, we get:
;     (1 - 1/2 A) exp(A) =  1 + 1/2 A

; This can be seen as a matrix equation in which the unknown is exp(A).
; The solution can be computed by `mldivide` ("left dividing with A").
;     exp(A) = mldvide( (1 - 1/2 A), 1 + 1/2 A)

; This basic idea is the back bone of the computation.

; There is a little observation that we can use to reduce
; amount of compuation in the denominator. The terms of the
; numerator and denominator are the same - except for different
; signs in the terms with odd degree. This suggests that we
; should calculate the even and odd terms of the numerator first:

;    numerator  = sum_of_even + sum_of_odd

; then the the denominator is a sumple sum:

;    denominator = sum_of_even - sum_of_odd

; The last piece of the puzzle concerns the domain the magnitude
; of the matrix entries. The approximation works best, if
; the entries |a_ij|<=0.5. (using norms, if |A|₁ <=0.5 ).
; See [1]

; Dividing the entries of A with 2 reduces the norm
;     exp(A) = exp(A/2)^2
; and if we do it s times, we get:
;     exp(A) = ... = exp(A/2^s)^(2s)

; That is we use the Padé approximation on exp(A/2^s),
; then we computer the 2s'th power.


; [1] The Pade Method for computing the Matrix Exponential
;     M. Arioli, B. Codenotti, C. Fassino
;     https://www.sciencedirect.com/science/article/pii/0024379594001901


(require racket/list)
(require "flomat.rkt")

;;;
;;; Coefficients
;;;

(define (fact n) (if (= n 0) 1 (* n (fact (- n 1)))))

(define (ncoef p i) ; for numerator
  (/ (* (fact (+ p p (- i))) (fact p))
     (* (fact (+ p p))       (fact i) (fact (- p i)))))


(define (dcoef p i) ; for denominator
  (* (ncoef p i)
     (if (odd? i) -1 1)))

; We can compute the coefficients we need.
(define n8 (let ([p 8]) (for/list ([j (+ p 1)]) (* 1.0 (ncoef p j)))))
(define d8 (let ([p 8]) (for/list ([j (+ p 1)]) (* 1.0 (dcoef p j)))))

(define (poly cs x)
  ; given (list c0 c1 c2 ... cn) compute c0 + c1*x + c2*x² + ... + cn x^n
  (cond
    [(empty? cs) 0.]
    [else       (+ (first cs)
                   (* x (poly (rest cs) x)))]))

; We can test the accuracy of the Pade approximation on a single number:

(define (exp-pade x)
  (/ (poly n8 x)
     (poly d8 x)))

; Now:  (exp 1.0) and (exp-pade 1.0) give the same result.
; Note: using n7,d7 give an error on the last decimal when using doubles.

; The matrix version of poly looks like this:

(define (mpoly cs A)
  (define n (nrows A))
  (define I (eye   n))
  (define Z (zeros n))
  
  (define (loop cs)
    (cond
      [(empty? cs) Z]
      [else        (plus (.* (first cs) I)
                         (times A (loop (rest cs))))]))
  (loop cs))

; We can directly expm as below, because we want to compute the
; the even and odd terms separately.

(define (expm-pade0 A)
  (mldivide (mpoly d8 A) (mpoly n8 A)))

(define (even-coefs cs)
  (if (empty? cs) '() (cons (first cs) (odd-coefs (rest cs)))))

(define (odd-coefs cs)
  (if (empty? cs) '() (even-coefs (rest cs))))

; Now we without scaling, we have:

(define (expm-pade/no-scaling A)
  (define AA   (times A A))
  (define even          (mpoly (even-coefs n8)  AA))
  (define odd  (times A (mpoly (odd-coefs  n8)  AA)))
  (define num  (plus  even odd))
  (define den  (minus even odd))
  (mldivide den num))

; With pre-scaling we get:

(define (expm A)
  (define s     (find-scaling-exponent A))
  (define 2^s   (expt 2 s))
  (define A/2^s (./ A 2^s))
  (repeated-squaring (expm-pade/no-scaling A/2^s) s))

(define (repeated-squaring A s)
  (if (= s 0) A (repeated-squaring (times A A) (- s 1))))

(define (find-scaling-exponent A)
  (define N (norm A 1))
  (define s (+ 1 (ceiling (log N 2))))
  s)


;; (exp-pade 1.0)
;; (exp 1.0)
;; (exp-pade 2.0)
;; (exp 2.0)

;; (expm-pade0 (matrix '([1])))
;; (expm-pade0 (matrix '([1  0] [0  2])))

;; (expm-pade/no-scaling (matrix '([1])))
;; (expm-pade/no-scaling (matrix '([1  0] [0  2])))

;; (expm-pade (matrix '([1])))
;; (expm-pade (matrix '([1  0] [0  2])))

;; (expm-pade (matrix '([1  2] [3  4])))


;;; Notes:
;; Found: The Scaling and Squaring Method for the Matrix Exponential Revisited,
;         Nicolas J. Higham
;  TODO: Read it and see if we can improve the algorithm.

