# Improved Isogeny Finding

This repository contains the code that computes the isogeny of degree `d`
between two supersingular elliptic curves. It accompanies the paper "Improved
quantum algorithms for finding fixed-degree isogenies".

Finding isogenies between supersingular elliptic curves is a natural
algorithmic problem which is known to be equivalent to computing the curves'
endomorphism rings. Similarly, another natural algorithmic problem is to find
an isogeny of a specified degree `d`. This problem appears to be somewhat
different in nature, but it is also considered to be a hard problem in
isogeny-based cryptography.

The approach makes use of the equivalence of categories of isogenies between
supersingular elliptic curves and the left ideals of maximal orders in
quaternion algebra. Finding the isogeny of degree `d` is translated to the
problem of solving a quadratic form of four variables.

# Guessing Variables

We use Magma script `magma/quaternion.m` to generate example quadratic forms.
Some example quadratic forms are in `sample_quadratic_forms.sage`.

There are two algorithms implemented that can solve the quadratic form:

 * Coron, J.S.: Finding small roots of bivariate integer polynomial equations:
   A direct approach. In: CRYPTO. vol. 4622, pp. 379–394. Springer (2007)

 * Bauer, A., Joux, A.: Toward a rigorous variation of Coppersmith’s algorithm
   on three variables. In: Advances in Cryptology — EUROCRYPT 2007, Lec- ture
   Notes in Computer Science, vol. 4515, pp. 361–378. Springer Berlin
   Heidelberg (2007).

To find the roots of the forms you can either run:
```
sage sample_quadratic_forms.sage coron
```
or
```
sage sample_quadratic_forms.sage bauer-joux
```

For the bivariate forms, you can use `cornacchia.sage`.

`cornacchia.sage` contains the function `transform_to_cornacchia_and_solve`.
Just to check whether `D` in the translated form is always negative, the
function is called from `sample_quadratic_forms_sage`.

# The Order Embedding Problem

We use the Magma script `magma/order_embedding.m` together with the Sage script
`./generate_forms_order_embedding.sage` to automatically generate and test
solving forms, that we use for experiments.
