# Motes-Rodrigo, Albery, et al. 2024 (Ecology Letters)

An analysis linking hurricane-induced changes in social structure with epidemic risk in free-ranging rhesus macaques (Macaca mulatta). It uses epidemiological simulations, applied to observed data from the Cayo Santiago population in Puerto Rico.

The code was initially written mostly by Camille Testard and Alba Motes-Rodrigo, and later rewritten, combined and smoothed by Greg Albery.

## Script structure:

- All the code can be found in the `R` folder. 
- Initial raw data can be found in the `Data` folder. 

- All code will run automatically using `000_Master Code.R`, but it will take a while to do so.

- The scripts are numbered in order.

- All the scripts beginning with `00` will import and clean the raw observational data and create a series of simulated networks using BISoN. The output of this code is large (~800MB) and therefore not stored here.

- The subsequent scripts will:
- 01. Simulate 1000 epidemics on each group-by-population combination, and run a population-level model examining the factors driving epidemic risk.
- 02. Summarise the outputs of these epidemics.
- 03. Run models that produce estimates of the factors driving epidemic risk, at the individual level.
- 04. Produce figures.
- 05. Write out results.

Enjoy! If you have any questions at all about the approach, feel free to email (Greg Albery)[gfalbery@gmail.com].