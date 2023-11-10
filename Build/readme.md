A sample of results, can be generated using demo.ipynb
![image](https://github.com/StormieWormie/PlateDeformation/assets/46678214/4fa4860a-90eb-4da1-95ca-98253bc9d761)

- fdm_fast: fdm but with a lazy Neumann condition implementation
- fdm_accurate: fdm but with a better (second order) Neumann condition approximation
- fem_mixed: hatfunctions as basis functions, solve a couple system
- fem_exotic: cubic hermite basis functions, does something to the system directly but it sure isn't solving it (see, inconsistency with other methods)
