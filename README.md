# RT_Battery_Control

This repository provides the code associated with the pdf that is at the root of the repository. 

I.The Matlab Control folder includes all the Matlab code associated with the 3 phases of the proposed framework described in the paper:
  1. Aggregation phase of distributed assets to determine optimal export quantities and times
  2. Bidding phase and associated real curve for the virtually aggregated assets (for example if the proposed bids are not validated, the load curve determined in phase 1 must be adapted)
  3. the real time control algorithm, that can run using Model Precitvie Control in Battery Controller that include Python capabilities (optimisation capability), or an approached version of MPC for Battery controllers that cannot run linear programming algorithms

II. The Smart Contract folder that provides:
 1. A general tutorial about Smart Contract for local replication
 2. the Python code to coordinate the interaction between the blockchain nodes
 3. the solidity code that includes all the Smart Contract code proposed in the paper.

