# STARK Proof System (w/ FRI)

This repository contains a rust port of [this great tutorial on STARKs from Alan Szepieniec](https://aszepieniec.github.io/stark-anatomy/). 

I worked on this project over the summer of 2024.


All tests pass for:
- field algebra
- polynomial algebra 
- merkle building
- FRI 

However, there is a bug somewhere within either the STARK prover or verifier that is causing proof verification to fail.


I will also note that this repository was not meant to be an optimized STARK prover. Rather, it was for my own understanding of the STARK protocol and FRI/IOPPs. 


    cargo test
