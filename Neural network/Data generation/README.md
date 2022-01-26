### DATA1_generator.ipynb
For generating DATA~1~, we chose four types of velocity models: *uniform* (constant velocity throughout), *2 parts* (spatial domain divided into two unequal parts with different velocities), *3 parts* (spatial domain divided into three unequal parts with different velocities), and *4 parts* (spatial domain divided into four unequal parts with different velocities). The values of velocity were random but within the physical constraints defined above. We then defined a number of source locations (randomly chosen between 1 to 10) and assigned each with some randomly chosen amplitude. This amplitude distribution was allowed to propagate for 1001 time steps using Finite Difference method (FDM), and we randomly chose 25 instances from the last 500 steps into our data. The next 10 time samples corresponding to the picked instances are the predictions we needed to make. Therefore, for each amplitude distribution and velocity distribution, we have 25 samples.

### DATA2_generator.ipynb
For generating DATA~2~ we used analytical Green’s function in propagating the impulse instead of the FDM. Impulse functions with different magnitudes were chosen at a number of randomly selected locations, and the corresponding Green’s function was computed for the 1001 time steps for each of them and superimposed. Again 25 randomly chosen time samples were taken, but this time it was from throughout the dataset. A constraint that this dataset had was that we only had to choose uniform velocity models because it is not easy to estimate the analytical solutions for Green’s function for inhomogeneous velocity models. We, therefore, used DATA~2~ to fine tune the model obtained from the model trained on DATA1, and for separate training as well.