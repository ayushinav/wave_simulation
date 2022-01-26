Presently, it's only a collection of models I tested on. Different scripts for different models will be added soon.

Temoral propagation was done using RNNs only. 

### Models
1) Using `Dense` layers: Only Dense layers were used in the architecture which were fed to RNNs.
2) Using `Conv1D` layers: 1D Convolutional layers were used to capture the features in the amplitudes and the velocity distributions, which were flattened and fed to RNNs.