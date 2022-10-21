import lib_functions
import numpy as np

# 2*bar(eta)>alpha and delta*beta > eta_h and alpha < delta/(math.sqrt(delta*beta/eta)-1)

#lib_functions.plot_graphs_gamma()
# lib_functions.plot_graphs_alpha(alpha_vector=np.linspace(.75, 1.75, 3))
#lib_functions.plot_graphs_alpha(alpha_vector=np.linspace(1.25, 1.75, 3))
# lib_functions.plot_graphs_eta()
# lib_functions.plot_graphs_beta()


lib_functions.plot_graphs_theta(alpha_local=0.50)
lib_functions.plot_graphs_theta(alpha_local=1.50)




