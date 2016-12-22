# prism
Allows SWI Prolog to drive an instance of the PRISM probabilistic programming system.

PRISM (PRogramming in Statistical Modelling) is a system invented by Taisuke Sato to do
probabilistic programming in a Prolog-like environment, using predicate tabling to implement
efficient inference and learning. It is built on top of Sato's B Prolog, which is a fast
Prolog but does not provide the libraries and comfortable development environment of
SWI Prolog. This pack allows SWI Prolog to manage and communicate with an instance of
PRISM, looking after its state, restoring it when it PRISM crashes, etc.
