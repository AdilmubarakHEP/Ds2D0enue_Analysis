#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

# A simple example to use bayesian optimization for the hyperparameters of a FastBDT.
# The package used in this example is https://github.com/scikit-optimize
# and can be installed with
# pip3 install scikit-optimize

from basf2 import find_file
import basf2_mva
import basf2_mva_util
import skopt
from skopt import plots
import matplotlib.pyplot as plt


def f(x):
    """Returns the figure of merit for the optimization.
    The functions trains the classifier with the given hyperparameters on the training sample and
    calculates the AUC on the independent test sample.
    """
    g_options = general_options
    g_options.m_identifier = "test.xml"
    options = basf2_mva.FastBDTOptions()
    options.m_nTrees = int(x[0])
    options.m_nLevels = int(x[1])
    basf2_mva.teacher(g_options, options)
    m = basf2_mva_util.Method(g_options.m_identifier)
    p, t = m.apply_expert(test_data, general_options.m_treename)
    return -basf2_mva_util.calculate_auc_efficiency_vs_background_retention(p, t)


if __name__ == "__main__":

    train_data = '/home/belle2/amubarak/C02-MVA/Splits/Ds2D0e_train.root'
    test_data = '/home/belle2/amubarak/C02-MVA/Splits/Ds2D0e_test.root'

    training_data = basf2_mva.vector(train_data)
    testing_data = basf2_mva.vector(test_data)

    # Define the variables for training.
    variables = ['K_kaonID',
                 'pi_pionID',
                 'D0_decayAngle_0','D0_chiProb','D0_flightDistance','D0_useCMSFrame_p','D0_M']

    general_options = basf2_mva.GeneralOptions()
    general_options.m_datafiles = training_data
    general_options.m_treename = "Dstree"
    general_options.m_variables = basf2_mva.vector(*variables)
    general_options.m_target_variable = "D0_isSignal"

    # Start optimization
    res = skopt.gp_minimize(f,  # the function to minimize
                            [(10, 1000), (2, 6)],  # the bounds on each dimension of x
                            x0=[10, 2],  # initial guess
                            n_calls=20)  # number of evaluations of f

    # Give some results
    print(res)
    skopt.plots.plot_convergence(res)
    plt.savefig('/home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Optimization/convergence.png')
    plots.plot_evaluations(res)
    plt.savefig('/home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Optimization/evaluations.png')
    plots.plot_objective(res)
    plt.savefig('/home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Optimization/objective.png')

    # Store result of optimization
    skopt.dump(res, '/home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Optimization/opt-result.pkl')
