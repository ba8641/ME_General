# Insights on enhancing affordability and profit in a global non-cooperative coordinated pediatric vaccine market

This repository contains files related to the academic paper "Insights on enhancing affordability and profit in a global non-cooperative coordinated pediatric vaccine market".
The files are organized as follows:

- Figures.zip contains all Figures that were explored in our research, including all permutations of the four factors involved (number of market segments, number of coordinating entities, order of negotiation and discount policy). Those results are shown at three price points: LP refers to lowest prices, and shows results at the lowest prices that would still keep the problem feasible. HP refers to highest prices, which is the opposite. Figures not labeled with either LP or HP are shown at the mean price between both.
- ABP_base.mod, ABP_ME.mod and ABP_Distribution.mod are all AMPL model files that specifies the objective functions that were used in the research, as well as additional constraints or functions that were explored but did not make into the final version of this paper.
- ABP_FL_test.run is the main AMPL script that runs all experiments. Its variations ABP_FL_ME_test.run and ABP_FL_MEC_test.run are called in the code, respectively, when exploring no discounts or including discounts.
- ExpFileDetailsTestFL.txt and EntFileDetailsFL.txt are input files for the AMPL script that specify which and in what order the datafiles should be read to properly represent the configuration of markets and entities being simulated.
