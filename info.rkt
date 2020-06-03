#lang info
;;; Info file for sci which contains multiple collections

;; Name
;   No name since `sci` consists of multiple collections.
(define collection 'multi) 

;; Version
(define version "1.0")

;; Dependencies
(define deps '("base"
               ))

;; 
(define implies '())

;; Package Description
(define pkg-desc "Scientific libraries: flomat (floating point matrices)")

(define pkg-authors '(soegaard))
(define build-deps '("rackunit-lib"
                     "scribble-lib"
                     "scribble-math"
                     "math-doc"
                     "racket-doc"
                     ("linux-shared-libraries"  #:platform "x86_64-linux-natipkg")))
