#lang info
;;; Info file for sci/flomat.

;; Name
;   The collection name can be different from the directory name,
;   Here they are the same.
(define collection "flomat-doc") 

;; Version
(define version "1.0")

;; Dependencies

(define deps       '("base"
                     "flomat"))

(define scribblings '(("manual-flomat.scrbl" () "Math and Science")))

(define build-deps '("math-doc"
                     "racket-doc"
                     "scribble-lib"
                     "scribble-math"
                     ("linux-shared-libraries"  #:platform "x86_64-linux-natipkg")))
