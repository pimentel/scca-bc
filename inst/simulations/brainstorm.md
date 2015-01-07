# Ideas for simulation framework

# Simulation creation
- Creating simulations should be easy
    - Provide a function for performing simulation
        - Function should return a list with the raw data (matrix), and a list
          with the row and column indicies for the solution.
          - Provide a wrapper function for solution i.e. `sccab_sol(rows,
            cols)` which returns an object of `sccab_sol``
    - Provide number of iterations
    - Provide a seed
    - Returns an unnamed list of simulations with raw data and `sccab_sol`
      objects

# Running the simulation
- Maybe this is best suited by just using a custom lapply for each method
- Want to parallelize
    - What will be better -- parallelize at the simulation level or at the
      iteration level?
    - Can delay for now -- pretty easy to implement after the fact

# Evaluating the simulation
- Have another function for computing some summary statistics
    - Given the truth provided above, compare the true `sccab_sol` to the
      computed `sccab_sol`
      - This requires having a comparison function that takes two `sccab_sol`
        objects, and returns a statistic
