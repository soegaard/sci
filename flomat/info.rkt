#lang info
;;; Info file for sci/flomat.

;; Name
;   The collection name can be different from the directory name,
;   Here they are the same.
(define collection "flomat") 

;; Version
(define version "1.0")

;; Dependencies

(define deps       '("base"
                     ; The shared libraries
                     ("flomat-x86_64-linux-natipkg"  #:platform "x86_64-linux-natipkg")))

(define build-deps '("rackunit-lib"))
